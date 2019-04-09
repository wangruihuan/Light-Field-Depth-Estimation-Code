#include <iostream>

#include <rslf_io.hpp>
#include <experimental/filesystem>


#define _RSLF_IO_VERBOSE
//~ #define _RSLF_IO_DEBUG


rslf::Mat rslf::read_img_from_file
(
    std::string path_to_folder, 
    std::string name_we,
    std::string extension,
    int cv_read_mode,
    bool transpose,
    bool rotate_180
) 
{
    if (path_to_folder.back() != '/')
        path_to_folder = path_to_folder + "/";
    
    if (extension[0] == '.')
        extension = extension.substr(1, extension.size() - 1);
    
    Mat img = cv::imread(path_to_folder + name_we + "." + extension, cv_read_mode);
    if (img.empty())
        std::cout << ">>> rslf::read_img_from_file -- WARNING: Empty image: " <<
            path_to_folder + name_we + "." + extension +
            " -- cv_read_mode=" << std::to_string(cv_read_mode) << std::endl;
    else
#ifdef _RSLF_IO_DEBUG
        std::cout << img.at<float>(0, 0) << std::endl;
#endif
    
    if (transpose)
        cv::transpose(img, img);
    
    if (rotate_180)
        cv::rotate(img, img, cv::ROTATE_180);
    
    return img;
}

rslf::Vec<rslf::Mat> rslf::read_imgs_from_folder
(
    std::string path_to_folder,
    std::string extension,
    int cv_read_mode,
    bool transpose,
    bool rotate_180
) 
{
    if (path_to_folder.back() != '/')
        path_to_folder = path_to_folder + "/";
    
    if (extension[0] == '.')
        extension = extension.substr(1, extension.size() - 1);
    
    // Get all valid item names in directory
    Vec<std::string> list_files;
    for (auto & p : std::experimental::filesystem::directory_iterator(path_to_folder)) 
    {
        std::string file_name = p.path().string();
        std::string current_extension = file_name.substr(file_name.find_last_of(".") + 1);
        std::string tmp = file_name.substr(file_name.find_last_of("/") + 1);
        std::string current_name = tmp.substr(0, tmp.find_last_of("."));
        if (current_extension == extension) 
        {
            list_files.push_back(current_name);
        } 
    }
    
    std::sort(list_files.begin(), list_files.end());
    
    // Read items
    Vec<Mat> imgs;
    for (std::string current_name : list_files)
    {
#ifdef _RSLF_IO_VERBOSE
        std::cout << "Read " << current_name << "." << extension << std::endl;
#endif
        Mat img = read_img_from_file(path_to_folder, current_name, extension, cv_read_mode);
        
        if (transpose)
            cv::transpose(img, img);
        
        if (rotate_180)
            cv::rotate(img, img, cv::ROTATE_180);
        
        imgs.push_back(img);
    }
    
    return imgs;
}

void rslf::write_mat_to_yml
(
    Mat img,
    std::string path_to_folder, 
    std::string name_we,
    std::string extension
)
{
    if (path_to_folder.back() != '/')
        path_to_folder = path_to_folder + "/";
    
    if (extension[0] == '.')
        extension = extension.substr(1, extension.size() - 1);
    
#ifdef _RSLF_IO_VERBOSE
    std::cout << "Write " << path_to_folder + name_we + "." + extension << std::endl;
#endif
    cv::FileStorage storage(path_to_folder + name_we + "." + extension, cv::FileStorage::WRITE);
    storage << "img" << img;
    storage.release();  
}

void rslf::write_mat_to_imgfile
(
    Mat img,
    std::string path_to_folder, 
    std::string name_we,
    std::string extension,
    Vec<int> compression_params
)
{
#ifdef _RSLF_IO_VERBOSE
    std::cout << "Write " << path_to_folder + name_we + "." + extension << std::endl;
#endif
    cv::imwrite(path_to_folder + name_we + "." + extension, img);
}

rslf::Mat rslf::read_mat_from_yml
(
    std::string path_to_folder, 
    std::string name_we,
    std::string extension
)
{
    if (path_to_folder.back() != '/')
        path_to_folder = path_to_folder + "/";
    
    if (extension[0] == '.')
        extension = extension.substr(1, extension.size() - 1);
    
#ifdef _RSLF_IO_VERBOSE
    std::cout << "Read " << path_to_folder + name_we + "." + extension << std::endl;
#endif
    Mat img;
    cv::FileStorage storage(path_to_folder + name_we + "." + extension, cv::FileStorage::READ);
    storage.getFirstTopLevelNode() >> img;
    storage.release();
    return img;
}

rslf::Mat rslf::build_row_epi_from_imgs
(
    Vec<Mat> imgs,
    int row,
    bool transpose,
    bool rotate_180
)
{
    Mat img0 = imgs[0];
    Mat epi
    (
        imgs.size(), // rows
        img0.cols, // cols
        img0.type() // dtype
    );
#ifdef _RSLF_EPI_DEBUG
    std::cout << "Created Mat of size (" << imgs.size() << "x" << 
        img0.cols << "), dtype=" << rslf::type2str(img0.type()) << std::endl;
#endif

    // Builds rows by copy
#pragma omp parallel for
    for (int i = 0; i < imgs.size(); i++)
    {
        imgs[i].row(row).copyTo(epi.row(i));
    }
    
    if (transpose)
        cv::transpose(epi, epi);
    
    if (rotate_180)
        cv::rotate(epi, epi, cv::ROTATE_180);
    
    return epi;
}

rslf::Vec<rslf::Mat> rslf::build_epis_from_imgs
(
    Vec<Mat> imgs,
    bool transpose,
    bool rotate_180
)
{
    Vec<Mat> epis;//(imgs[0].rows, cv::Mat::zeros(imgs.size(), imgs[0].cols, imgs[0].type()));
    //~ #pragma omp parallel for
    for (int v=0; v<imgs[0].rows; v++)
    {   
        Mat epi(imgs.size(), imgs[0].cols, imgs[0].type());
        // Builds rows by copy
        for (int s=0; s<imgs.size(); s++)
        {
            imgs[s].row(v).copyTo(epi.row(s));
        }
        
        if (transpose)
            cv::transpose(epi, epi);
        
        if (rotate_180)
            cv::rotate(epi, epi, cv::ROTATE_180);
        
        epis.push_back(epi);
        
    }
    //~ for (int v=0; v<imgs[0].rows; v++)
    //~ {   
        //~ std::cout << epis[v].at<float>(0, 0) << std::endl;
    //~ }
    
    return epis;
}

rslf::Mat rslf::build_row_epi_from_path
(
    std::string path_to_folder,
    std::string extension,
    int row,
    int cv_read_mode,
    bool transpose,
    bool rotate_180
)
{
    if (path_to_folder.back() != '/')
        path_to_folder = path_to_folder + "/";
    
    if (extension[0] == '.')
        extension = extension.substr(1, extension.size() - 1);
    
    // Get all valid item names in directory
    Vec<std::string> list_files;
    for (auto & p : std::experimental::filesystem::directory_iterator(path_to_folder)) 
    {
        std::string file_name = p.path().string();
        std::string current_extension = file_name.substr(file_name.find_last_of(".") + 1);
        std::string tmp = file_name.substr(file_name.find_last_of("/") + 1);
        std::string current_name = tmp.substr(0, tmp.find_last_of("."));
        if (current_extension == extension) 
        {
            list_files.push_back(current_name);
        } 
    }
    
    std::sort(list_files.begin(), list_files.end());
    
    // Sample image
    Mat tmp = read_img_from_file
        (
            path_to_folder, 
            list_files[0], 
            extension, 
            cv_read_mode
        );
    
    // Read items
    Mat epi
    (
        list_files.size(), // rows
        tmp.cols, // cols
        tmp.type() // dtype
    );
    
    // For each file, read the corresponding row
#pragma omp parallel for
    for (int i=0; i<list_files.size(); i++)
    {
#ifdef _RSLF_EPI_DEBUG
        std::cout << "Read " << current_name << "." << extension << std::endl;
#endif
        Mat img = read_img_from_file(path_to_folder, list_files[i], extension, cv_read_mode);
        img.row(row).copyTo(epi.row(i));
    }
    
    if (transpose)
        cv::transpose(epi, epi);
    
    if (rotate_180)
        cv::rotate(epi, epi, cv::ROTATE_180);
    
    return epi;
}
