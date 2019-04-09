#include <iostream>

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <rslf_plot.hpp>


#define _RSLF_MAX_WINDOW_DIMENSION 800


void rslf::plot_mat
(
    Mat img, 
    std::string window_name
)
{
    // Copy and scale values
    Mat tmp = rslf::copy_and_scale_uchar(img);
    
    // Plot in window
    // Ratio
    float ratio = (tmp.cols + 0.0) / tmp.rows;
    int window_cols, window_rows;
    if (tmp.cols >= tmp.rows)
    {
        window_cols = std::min(_RSLF_MAX_WINDOW_DIMENSION, tmp.cols);
        window_rows = (int)std::ceil(window_cols / ratio);
    }
    else
    {
        window_rows = std::min(_RSLF_MAX_WINDOW_DIMENSION, tmp.rows);
        window_cols = (int)std::ceil(window_rows * ratio);
    }
    cv::namedWindow(window_name, CV_WINDOW_NORMAL);
    cv::resizeWindow(window_name, window_cols, window_rows);
    cv::imshow(window_name, tmp);
}

rslf::Mat rslf::copy_and_scale_uchar
(
    Mat img
)
{
    Mat res;
    img.copyTo(res);
    // If the dtype is not a multiple of uchar
    if (img.type() % 8 != 0) {
        // Copy and scale values
        // Find min and max values
        double min, max;
        cv::minMaxLoc(img, &min, &max);
        res -= min;
        // Copy and scale values to uchar
        res.convertTo(res, CV_8U, 255.0/(max-min));
    } 
    //~ else 
    //~ {
        //~ // Only copy values
        //~ img.copyTo(res);
    //~ }
    return res;
}


void rslf::ImageConverter_uchar::fit(const Mat& img, bool saturate)
{
    Mat reshaped_img = img.reshape(1, img.rows * img.cols);
    std::cout << "reshaped_img.size()=" << reshaped_img.size() << std::endl;
    Mat sorted_img;
    cv::sort(reshaped_img, sorted_img, CV_SORT_EVERY_COLUMN + CV_SORT_ASCENDING);
    
    if (saturate)
    {
        // min=0.02% quantile
        int idx_min = (int)std::floor(0.02 * reshaped_img.rows);
        min = sorted_img.at<float>(idx_min);
        // max=0.98% quantile
        int idx_max = (int)std::floor(0.98 * reshaped_img.rows);
        max = sorted_img.at<float>(idx_max);
    
    
        std::cout << "idx_min=" << idx_min << ", idx_max=" << idx_max << std::endl;
        std::cout << "min=" << min << ", max=" << max << std::endl;
    }
    else
    {
        // Find min and max values
        double true_min, true_max;
        cv::minMaxLoc(img, &true_min, &true_max);
        cv::Scalar mean, std;
        cv::meanStdDev(img, mean, std);
        min = std::max(0.0, true_min);
        min = true_min;
        max = std::min(mean[0] + 12 * std[0], true_max); 
    }
    m_initialized = true;
}

void rslf::ImageConverter_uchar::copy_and_scale(const Mat& src, Mat& dst)
{
    assert(m_initialized && "ImageConverter should be fitted before usage");
    
    // Copy and scale values to uchar
    float alpha = 255.0/(max-min);
    src.convertTo(dst, CV_8U, alpha, -alpha * min);
}


rslf::Mat rslf::draw_red_lines
(
    Mat img,
    int fill_row_red,
    int max_height,
    int fill_col_red,
    int max_width
)
{
    // Copy and scale values
    Mat res = rslf::copy_and_scale_uchar(img);
    
    if (fill_row_red < 0 && fill_col_red < 0)
        // Nothing to do
        return res;
    
    if (img.channels() == 1)
        // Convert to RGB
        cv::cvtColor(res, res, CV_GRAY2RGB);
    
    // Red color
    cv::Vec3b red;
    red.val[2] = 255;
    
    // Fill row?
    if (fill_row_red > -1) 
    {
        for (int j=0; j<res.cols; j++)
            res.at<cv::Vec3b>(fill_row_red, j) = red;
    }
    
    // Fill col?
    if (fill_col_red > -1) 
    {
        for (int i=0; i<res.rows; i++)
            res.at<cv::Vec3b>(i, fill_col_red) = red;
    }
    
    // Truncate the image along rows
    if (fill_row_red > -1 && max_height > 0)
    {
        int first_row, last_row;
        // Get the first and last rows to be copied
        // First row
        if (fill_row_red - max_height < 0)
            first_row = 0;
        else
            first_row = fill_row_red - max_height / 2;
        // Last row
        if (first_row + max_height < res.rows)
            last_row = first_row + max_height;
        else
            last_row = res.rows - 1;
        
        // Get a copy of the image between first and last row
        Mat tmp;
        cv::Range ranges[2];
        ranges[0] = cv::Range(first_row, last_row); 
        ranges[1] = cv::Range::all();
        res(ranges).copyTo(tmp);
        res = tmp;
    }
    
    // Truncate the image along cols
    if (fill_col_red > -1 && max_width > 0)
    {
        int first_col, last_col;
        // Get the first and last cols to be copied
        // First col
        if (fill_col_red - max_width < 0)
            first_col = 0;
        else
            first_col = fill_col_red - max_width / 2;
        // Last row
        if (first_col + max_width < res.cols)
            last_col = first_col + max_width;
        else
            last_col = res.cols - 1;
        
        // Get a copy of the image between first and last row
        Mat tmp;
        cv::Range ranges[2];
        ranges[0] = cv::Range::all();
        ranges[0] = cv::Range(first_col, last_col); 
        res(ranges).copyTo(tmp);
        res = tmp;
    }
    
    return res;
}
