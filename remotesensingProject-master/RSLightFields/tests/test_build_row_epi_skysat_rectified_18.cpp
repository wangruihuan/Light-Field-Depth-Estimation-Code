#include <opencv2/core/core.hpp>
#include <string>
#include <iostream>
#include <rslf.hpp>
#include <chrono>

/*
 * test_build_row_epi_skysat_rectified_18.cpp
 * 
 * Given a folder containing the SkySat images, build an EPI corresponding
 * to a given row and plot it.
 * Provide a desired row in argument (otherwise default).
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
    
    cv::Mat epi = rslf::build_row_epi_from_imgs(list_mat, inspected_row);
    
    // Save imgs 
    // Base name
    // https://stackoverflow.com/questions/31255486/c-how-do-i-convert-a-stdchronotime-point-to-long-and-back
    auto now = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    auto ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now).time_since_epoch().count();
    
    std::string base_filename = std::to_string(ms);
    std::cout << "Base name: " << base_filename << std::endl;
    
    cv::Mat first_image_overlayed = rslf::draw_red_lines
    (
        list_mat[0],
        inspected_row,
        epi.rows
    );
    
    rslf::write_mat_to_imgfile
    (
        first_image_overlayed,
        "../output/", 
        base_filename + "_1st",
        "png"
    );
    rslf::write_mat_to_imgfile
    (
        epi,
        "../output/", 
        base_filename + "_epi",
        "png"
    );
    
    // Plot the first image overlayed with a red line corresponding
    // to the epipolar plane
    rslf::plot_mat(first_image_overlayed, "000.tif");
    rslf::plot_mat(epi, "EPI");
    
    cv::waitKey();
}

