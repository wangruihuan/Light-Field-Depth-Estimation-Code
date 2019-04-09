#include <opencv2/core/core.hpp>
#include <string>
#include <iostream>
#include <rslf.hpp>
#include <chrono>

/*
 * test_build_row_epi_mansion_rectified_18.cpp
 * 
 * Given a folder containing the resized mansion images, build an EPI corresponding
 * to a given row and plot it.
 * Provide a desired row in argument (otherwise default).
 */
 

int main(int argc, char* argv[])
{
    
    // Select row
    int inspected_row = 380;
    if (argc > 1)
    {
        std::cout << "Selected line: " << argv[1] << std::endl;
        inspected_row = std::stoi(argv[1]);
    }
    
    // Read first image
    cv::Mat first_image = rslf::read_img_from_file("../data/mansion_image_resized/", "mansion_image_0000.jpg.resized", "jpg", CV_LOAD_IMAGE_UNCHANGED); 
    // Build epi
    cv::Mat epi = rslf::build_row_epi_from_path("../data/mansion_image_resized/", "jpg", inspected_row, CV_LOAD_IMAGE_UNCHANGED);
    
    // Save imgs 
    // Base name
    // https://stackoverflow.com/questions/31255486/c-how-do-i-convert-a-stdchronotime-point-to-long-and-back
    auto now = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    auto ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now).time_since_epoch().count();
    
    std::string base_filename = std::to_string(ms);
    std::cout << "Base name: " << base_filename << std::endl;
    
    cv::Mat first_image_overlayed = rslf::draw_red_lines
    (
        first_image,
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
    rslf::plot_mat(first_image_overlayed, "mansion_image_0000.jpg.resized");
    rslf::plot_mat(epi, "EPI");
    
    cv::waitKey();
    
}

