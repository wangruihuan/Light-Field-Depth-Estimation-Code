#include <opencv2/core/core.hpp>
#include <string>
#include <iostream>
#include <rslf.hpp>
#include <rslf_depth_computation.hpp>
#include <chrono>

/*
 * test_depth_computation.cpp
 * 
 * Compute 1D depth (on center line of an epi).
 */
 
int read_skysat_lasvegas_rectified(int inspected_row);
int read_mansion_image(int inspected_row);

int main(int argc, char* argv[])
{

    // Select row
    int inspected_row = 380;
    if (argc > 1)
    {
        std::cout << "Selected line: " << argv[1] << std::endl;
        inspected_row = std::stoi(argv[1]);
    }
    
    // Load all images in folder
    std::vector<cv::Mat> list_mat = rslf::read_imgs_from_folder("../data/skysat_lasvegas_rectified/rectified_equalized_resized_frames_step18/", "tif", CV_LOAD_IMAGE_UNCHANGED);
    //~ std::vector<cv::Mat> list_mat = rslf::read_imgs_from_folder("../data/mansion_image_resized/", "jpg", CV_LOAD_IMAGE_UNCHANGED);
    
    std::cout << list_mat.size() << " images read" << std::endl;
    
    cv::Mat epi = rslf::build_row_epi_from_imgs(list_mat, inspected_row);
    
    float d_min = -2.0;
    float d_max = 4.0;
    int dim_d = 120;
    
    std::cout << dim_d << " d values requested" << std::endl;
    
    rslf::Depth1DParameters<float>& parameters = rslf::Depth1DParameters<float>::get_default();
    
    //~ rslf::DepthComputer1D<cv::Vec3f> depth_computer_1d(epi, d_min, d_max, dim_d);
    rslf::Depth1DComputer<float> depth_computer_1d(epi, d_min, d_max, dim_d);
    depth_computer_1d.run();
    
    std::cout << "Plotting" << std::endl;
    
    cv::Mat coloured_epi = depth_computer_1d.get_coloured_epi();
    
    // Save imgs 
    // Base name
    // https://stackoverflow.com/questions/31255486/c-how-do-i-convert-a-stdchronotime-point-to-long-and-back
    auto now = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    auto ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now).time_since_epoch().count();
    
    std::string base_filename = std::to_string(ms);
    std::cout << "Base name: " << base_filename << std::endl;
    
    rslf::write_mat_to_imgfile
    (
        epi,
        "../output/", 
        base_filename + "_epi",
        "png"
    );
    rslf::write_mat_to_imgfile
    (
        coloured_epi,
        "../output/", 
        base_filename + "_epi_colored",
        "png"
    );
    
    rslf::plot_mat(epi, "EPI");
    rslf::plot_mat(coloured_epi, "EPI + depth");
    
    cv::waitKey();
    
    
    
    return 0;
}

