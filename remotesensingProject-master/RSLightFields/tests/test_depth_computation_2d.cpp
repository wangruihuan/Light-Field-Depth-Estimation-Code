#include <opencv2/core/core.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <rslf.hpp>
#include <rslf_depth_computation.hpp>
#include <chrono>

/*
 * test_depth_computation_2d.cpp
 * 
 * Compute 1D depth (on center line of an epi).
 */
 
int read_skysat_lasvegas_rectified(int inspected_row);
int read_mansion_image(int inspected_row);

int main(int argc, char* argv[])
{
    
    std::cout << "Started depth computation..." << std::endl;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // Select row
    int inspected_row = 380;
    if (argc > 1)
    {
        std::cout << "Selected line: " << argv[1] << std::endl;
        inspected_row = std::stoi(argv[1]);
    }
    
    // Load all images in folder
    std::vector<cv::Mat> list_mat = rslf::read_imgs_from_folder("../data/skysat_lasvegas_rectified/rectified_equalized_resized_frames_step18/", "tif", CV_LOAD_IMAGE_UNCHANGED);
    //~ std::vector<cv::Mat> list_mat = rslf::read_imgs_from_folder("../data/sieleters_bande2_short/", "tif", CV_LOAD_IMAGE_UNCHANGED, true, true);
    //~ std::vector<cv::Mat> list_mat = rslf::read_imgs_from_folder("../data/mansion_image_resized/", "jpg", CV_LOAD_IMAGE_UNCHANGED);
    
    std::cout << list_mat.size() << " images read" << std::endl;
    
    std::vector<cv::Mat> epis = rslf::build_epis_from_imgs(list_mat);
    
    //~ rslf::plot_mat(epis[0], "EPI 0");
    //~ rslf::plot_mat(epis[500], "EPI 500");
    //~ cv::waitKey();
    
    float d_min = -1.0;
    float d_max = 4;
    //~ float d_min = -5.5;
    //~ float d_max = 5.5;
    //~ float d_max = 0.5;
    //~ int dim_d = 240;
    int dim_d = 120;
    
    std::cout << dim_d << " d values requested" << std::endl;
    
    rslf::Depth2DComputer_1ch depth_computer_2d(epis, d_min, d_max, dim_d);
    //~ rslf::Depth2DComputer<cv::Vec3f> depth_computer_2d(epis, d_min, d_max, dim_d);
    depth_computer_2d.run();
    
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    
    // Save imgs 
    // Base name
    // https://stackoverflow.com/questions/31255486/c-how-do-i-convert-a-stdchronotime-point-to-long-and-back
    auto now = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    auto ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now).time_since_epoch().count();
    
    std::string base_filename = std::to_string(ms);
    std::cout << "Base name: " << base_filename << std::endl;
    
    std::cout << "Plotting" << std::endl;
    
    for (int s=0; s<epis[0].rows; s++)
    {
        cv::Mat disparity_map = depth_computer_2d.get_disparity_map(s, cv::COLORMAP_JET);
        
        std::stringstream ss;
        ss << std::setw(3) << std::setfill('0') << s;
        std::string s_str = ss.str();
        
        rslf::write_mat_to_imgfile
        (
            disparity_map,
            "../output/", 
            base_filename + "_dmap_" + s_str,
            "png"
        );
    }
    
    cv::Mat coloured_epi = depth_computer_2d.get_coloured_epi(-1, cv::COLORMAP_JET);
    rslf::write_mat_to_imgfile
    (
        coloured_epi,
        "../output/", 
        base_filename + "_epi_colored",
        "png"
    );
    
    rslf::write_mat_to_yml
    (
        depth_computer_2d.m_edge_confidence_s_v_u[50],
        "../output/",
        "output_edge_confidence"
    );
    rslf::write_mat_to_yml
    (
        depth_computer_2d.m_disp_confidence_s_v_u[50],
        "../output/",
        "output_disp_confidence"
    );
    //~ rslf::write_mat_to_yml
    //~ (
        //~ depth_computer_2d.m_line_confidence_s_v_u[50],
        //~ "../output/",
        //~ "output_line_confidence"
    //~ );
    
    // Write csv
    std::ofstream out_e;
    out_e.open("../output/output_edge_confidence.csv");
    //~ std::ofstream out_l;
    //~ out_l.open("../output/output_line_confidence.csv");
    for (int i=0; i<depth_computer_2d.m_edge_confidence_s_v_u[50].rows; i++)
    {
        for (int j=0; j<depth_computer_2d.m_edge_confidence_s_v_u[50].cols; j++)
        {
            if (depth_computer_2d.m_edge_confidence_mask_s_v_u[50].at<uchar>(i, j))
            {
                out_e << depth_computer_2d.m_edge_confidence_s_v_u[50].at<float>(i, j) << std::endl;
                //~ out_l << depth_computer_2d.m_line_confidence_s_v_u[50].at<float>(i, j) << std::endl;
            }
        }
    }
    out_e.close();
    //~ out_l.close();
    
    
    //~ rslf::plot_mat(disparity_map, "Disparity map");
    rslf::plot_mat(coloured_epi, "EPI + depth");
    
    cv::waitKey();
    
    
    
    return 0;
}

