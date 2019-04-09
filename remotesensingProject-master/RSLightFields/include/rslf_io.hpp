#ifndef _RSLF_IO
#define _RSLF_IO


#include <rslf_types.hpp>
#include <opencv2/highgui/highgui.hpp>


/*!
 * \file
 * \brief Implement i/o functions (read, write matrices).
 */


namespace rslf 
{


/*
 * *****************************************************************
 * Reading/Writing matrices
 * *****************************************************************
 */

    
/**
 * Reads an image from the given path, name and extension.
 * 
 * @param path_to_folder Path to the folder containing the file. Should end with '/' (but the function will correct if not the case).
 * @param name_we Name of the file without path or extension.
 * @param extension The extension of the file. Should not start with '.' (but the algorithm will correct if not the case).
 * @param cv_read_mode The mode with which OpenCV will read the file.
 * @return The requested image in Mat format.
 */
Mat read_img_from_file
(
    std::string path_to_folder, 
    std::string name_we,
    std::string extension,
    int cv_read_mode = CV_LOAD_IMAGE_UNCHANGED,
    bool transpose = false,
    bool rotate_180 = false
);

/**
 * Reads all images with requested extension in requested folder.
 * 
 * @param path_to_folder Path to the folder containing the files. Should end with '/' (but the function will correct if not the case).
 * @param extension The extension of the files. Should not start with '.' (but the algorithm will correct if not the case).
 * @param cv_read_mode The mode with which OpenCV will read the files.
 * @return The requested images in a vector of Mat.
 */
std::vector<Mat> read_imgs_from_folder
(
    std::string path_to_folder,
    std::string extension,
    int cv_read_mode = CV_LOAD_IMAGE_UNCHANGED,
    bool transpose = false,
    bool rotate_180 = false
);

/**
 * Writes the given Mat to a .yml file in the given folder.
 * 
 * @param img The Mat to save.
 * @param path_to_folder Path to the folder where to save the file. Should end with '/' (but the function will correct if not the case).
 * @param name_we Name of the file without path or extension.
 * @param extension The extension of the file. Should not start with '.' (but the algorithm will correct if not the case).
 */
void write_mat_to_yml
(
    Mat img,
    std::string path_to_folder, 
    std::string name_we,
    std::string extension = "yml"
);

/**
 * Writes the given Mat to an image in the given folder.
 * 
 * @param img The Mat to save.
 * @param path_to_folder Path to the folder where to save the file. Should end with '/' (but the function will correct if not the case).
 * @param name_we Name of the file without path or extension.
 * @param extension The extension of the file. Should not start with '.' (but the algorithm will correct if not the case).
 * @param compression_params Compression parameters.
 */
void write_mat_to_imgfile
(
    Mat img,
    std::string path_to_folder, 
    std::string name_we,
    std::string extension = "png",
    std::vector<int> compression_params = std::vector<int>()
);

/**
 * Reads the given .yml file in the given folder to a Mat.
 * 
 * @param path_to_folder Path to the folder where to read the file.
 * @param name_we Name of the file without path or extension.
 * @param extension The extension of the file.
 * @return A Mat containing the image
 */
Mat read_mat_from_yml
(
    std::string path_to_folder, 
    std::string name_we,
    std::string extension = "yml"
);


/*
 * *****************************************************************
 * Building EPIs
 * *****************************************************************
 */

/**
 * \brief Assume we have a set of rectified images such that the epipolar
 * planes are on the lines. Build the EPI along the given row.
 * 
 * @param imgs The set of images from which to build the EPI.
 * @param row The row along which the EPI will be built.
 * @param transpose Indicates whether the algorithm should transpose the images (epipolar lines should be horizontal)
 * @return The requested EPI.
 */
Mat build_row_epi_from_imgs
(
    Vec<Mat> imgs,
    int row,
    bool transpose = false,
    bool rotate_180 = false
);

/**
 * \brief Assume we have a set of rectified images such that the epipolar
 * planes are on the lines. Build the EPIs along the given row.
 * 
 * @param imgs The set of images from which to build the EPIs.
 * @param transpose Indicates whether the algorithm should transpose the images (epipolar lines should be horizontal)
 * @return The requested EPIs.
 */
Vec<Mat> build_epis_from_imgs
(
    Vec<Mat> imgs,
    bool transpose = false,
    bool rotate_180 = false
);

/**
 * \brief Assume we have a set of rectified images such that the epipolar
 * planes are on the lines in a folder. Build the EPI along the given row.
 * Images are supposed to be ordered in alphanumerical order.
 * 
 * @param path_to_folder The path containing the images. Should end with '/' (but the function will correct if not the case).
 * @param extension The extension of the images to be read. Should not start with '.' (but the algorithm will correct if not the case).
 * @param row The row along which the EPI will be built.
 * @param cv_read_mode The mode with which OpenCV will read the files.
 * @param transpose Indicates whether the algorithm should transpose the images (epipolar lines should be horizontal).
 * @return The requested EPI.
 */
Mat build_row_epi_from_path
(
    std::string path_to_folder,
    std::string extension,
    int row,
    int cv_read_mode = CV_LOAD_IMAGE_UNCHANGED,
    bool transpose = false,
    bool rotate_180 = false
);

}

#endif
