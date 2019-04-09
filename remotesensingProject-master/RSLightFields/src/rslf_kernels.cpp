#include <rslf_kernels.hpp>


template<>
float rslf::BandwidthKernel<float>::evaluate(float x) 
{
    // If none, return 0
    if (rslf::is_nan_type<float>(x))
        return 0;
    // Else return 1 - (x/h)^2 if (x/h)^2 < 1, else 0
    // In order to be coherent with 3channel thresholds: x3
    float tmp = 3 * std::pow(x / m_h_, 2);
    return (tmp > 1 ? 0 : 1 - tmp);
}

template<>
void rslf::BandwidthKernel<float>::evaluate_mat(const Mat& src, Mat& dst) 
{
    // (x/h)^2
    // In order to be coherent with 3channel thresholds: x3
    cv::multiply(src, src, dst, 3 * inv_m_h_sq);
    // 1 - (x/h)^2
    cv::subtract(1.0, dst, dst);
    // max( 1 - (x/h)^2, 0 ) ; this operation removes nan's
    cv::max(dst, 0.0, dst);
}

template<>
float rslf::BandwidthKernel<cv::Vec3f>::evaluate(cv::Vec3f x) 
{
    // If none, return 0
    if (rslf::is_nan_type<cv::Vec3f>(x))
        return 0;
    // Else return 1 - (x/h)^2 if (x/h)^2 < 1, else 0
    float tmp = std::pow(cv::norm(x) / m_h_, 2);
    return (tmp > 1 ? 0 : 1 - tmp);
}

template<>
void rslf::BandwidthKernel<cv::Vec3f>::evaluate_mat(const Mat& src, Mat& dst) 
{
    // (x/h)^2
    cv::multiply(src, src, dst, inv_m_h_sq);
    // sum over channels
    int rows = src.rows;
    int cols = src.cols;
    dst = dst.reshape(1, rows * cols);
    cv::reduce(dst, dst, 1, cv::REDUCE_SUM);
    dst = dst.reshape(1, rows);
    // 1 - (x/h)^2
    cv::subtract(1.0, dst, dst);
    // max( 1 - (x/h)^2, 0 ) ; this operation removes nan's
    cv::max(dst, 0.0, dst);
}
