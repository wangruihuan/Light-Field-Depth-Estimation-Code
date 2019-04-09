#include <iostream>

#include <rslf_depth_computation_core.hpp>


template<>
void rslf::_square_sum_channels_into<float>(const Mat& src, Mat& dst, Mat& buffer)
{
    cv::pow(src, 2, buffer);
    dst += buffer;
}

template<>
void rslf::_square_sum_channels_into<cv::Vec3f>(const Mat& src, Mat& dst, Mat& buffer)
{
    Vec<Mat> channels;
    cv::split(src, channels);
    for (int c=0; c<3; c++) 
    {
        cv::pow(channels[c], 2, buffer);
        dst += buffer;
    }
}

template<>
void rslf::_multiply_multi_channel<float>(const Mat& line_mat, const Mat& vec_mat, Mat& res_mat, Mat& buffer)
{
    cv::multiply(vec_mat, line_mat, res_mat);
}

template<>
void rslf::_multiply_multi_channel<cv::Vec3f>(const Mat& line_mat, const Mat& vec_mat, Mat& res_mat, Mat& buffer)
{
    Vec<Mat> channels(3, line_mat);
    cv::merge(channels, buffer);
    cv::multiply(vec_mat, buffer, res_mat);
}

template<>
void rslf::_divide_multi_channel<float>(const Mat& line_mat, const Mat& vec_mat, Mat& res_mat, Mat& buffer)
{
    cv::divide(vec_mat, line_mat, res_mat);
} 

template<>
void rslf::_divide_multi_channel<cv::Vec3f>(const Mat& line_mat, const Mat& vec_mat, Mat& res_mat, Mat& buffer)
{
    Vec<Mat> channels(3, line_mat);
    cv::merge(channels, buffer);
    cv::divide(vec_mat, buffer, res_mat);
}
