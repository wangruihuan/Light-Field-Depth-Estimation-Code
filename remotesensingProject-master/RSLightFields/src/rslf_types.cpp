#include <limits>

#include <rslf_types.hpp>


/*
 * NaN
 */
template<>
float rslf::nan_type<float>() 
{
    return std::numeric_limits<float>::quiet_NaN();
}

template<>
cv::Vec3f rslf::nan_type<cv::Vec3f>() 
{
    return cv::Vec3f(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
}

template<>
bool rslf::is_nan_type<float>(float x) 
{
    return x != x;
}

template<>
bool rslf::is_nan_type<cv::Vec3f>(cv::Vec3f x) 
{
    return x[0] != x[0] || x[1] != x[1] ||x[2] != x[2];
}

/*
 * Zero
 */
template<>
cv::Scalar rslf::zero_scalar<float>() 
{
    return cv::Scalar(0.0);
}

template<>
cv::Scalar rslf::zero_scalar<cv::Vec3f>() 
{
    return cv::Scalar(0.0, 0.0, 0.0);
}

/*
 * OpenCV type explanation
 */
// Courtesy of Octopus https://stackoverflow.com/a/17820615
std::string rslf::type2str(int type) 
{
    std::string r;

    uchar depth = type & CV_MAT_DEPTH_MASK;
    uchar chans = 1 + (type >> CV_CN_SHIFT);

    switch ( depth ) 
    {
    case CV_8U:  r = "8U"; break;
    case CV_8S:  r = "8S"; break;
    case CV_16U: r = "16U"; break;
    case CV_16S: r = "16S"; break;
    case CV_32S: r = "32S"; break;
    case CV_32F: r = "32F"; break;
    case CV_64F: r = "64F"; break;
    default:     r = "User"; break;
    }

    r += "C";
    r += (chans+'0');

    return r;
}

/*
 * Template norm
 */
template<>
float rslf::norm<float>(float x) 
{
    // In order to be consistent with 3channel thresholds: times sqrt(3)
    return std::abs(x) * 1.73205080757;
}

template<>
float rslf::norm<cv::Vec3f>(cv::Vec3f x) 
{
    return cv::norm(x);
}

template<>
void rslf::sin_mat<float>(const Mat& src, Mat& dst)
{
    src.copyTo(dst);
    dst.forEach<float>
    (
        [](float &pixel, const int* position) -> void
        {
            pixel = std::sin(pixel);
        }
    );
}

template<>
void rslf::cos_mat<float>(const Mat& src, Mat& dst)
{
    src.copyTo(dst);
    dst.forEach<float>
    (
        [](float &pixel, const int* position) -> void
        {
            pixel = std::cos(pixel);
        }
    );
}

template<>
void rslf::sin_mat<cv::Vec3f>(const Mat& src, Mat& dst)
{
    src.copyTo(dst);
    dst.forEach<cv::Vec3f>
    (
        [](cv::Vec3f &pixel, const int* position) -> void
        {
            pixel[0] = std::sin(pixel[0]);
            pixel[1] = std::sin(pixel[1]);
            pixel[2] = std::sin(pixel[2]);
        }
    );
}

template<>
void rslf::cos_mat<cv::Vec3f>(const Mat& src, Mat& dst)
{
    src.copyTo(dst);
    dst.forEach<cv::Vec3f>
    (
        [](cv::Vec3f &pixel, const int* position) -> void
        {
            pixel[0] = std::cos(pixel[0]);
            pixel[1] = std::cos(pixel[1]);
            pixel[2] = std::cos(pixel[2]);
        }
    );
}

template<>
void rslf::fft_htranslate<float>(const Mat& src, float tr_amount, Mat& dst)
{
    //~ Mat padded;                              // expand input image to optimal size
    //~ int m = cv::getOptimalDFTSize( src.rows );
    //~ int n = cv::getOptimalDFTSize( src.cols ); // on the border add zero values
    //~ cv::copyMakeBorder(src, padded, 0, m - src.rows, 0, n - src.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));

    //~ Mat planes[] = {cv::Mat_<float>(padded), cv::Mat::zeros(padded.size(), CV_32F)};
    Mat planes[] = {cv::Mat_<float>(src), cv::Mat::zeros(src.size(), CV_32F)};
    Mat complexI;
    cv::merge(planes, 2, complexI);         // Add to the expanded another plane with zeros

    cv::dft(complexI, complexI);            // this way the result may fit in the source matrix
    cv::split(complexI, planes);            // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))

    Mat real;
    Mat imag;

    Mat ph = Mat(planes[0].rows, planes[0].cols, CV_32FC1, cv::Scalar(- 2 * 3.1415926535898 * tr_amount / planes[0].cols));
    for (int j=0; j<planes[0].cols; j++)
    {
        cv::multiply(ph.col(j), j, ph.col(j));
    }
    Mat cos_ph;
    Mat sin_ph;
    cos_mat<float>(ph, cos_ph);
    sin_mat<float>(ph, sin_ph);
    
    Mat tmp1, tmp2;
    cv::multiply(planes[0], cos_ph, tmp1);
    cv::multiply(planes[1], sin_ph, tmp2);
    cv::subtract(tmp1, tmp2, real);
    
    cv::multiply(planes[1], cos_ph, tmp1);
    cv::multiply(planes[0], sin_ph, tmp2);
    cv::add(tmp1, tmp2, imag);
    
    // Replace planes
    planes[0] = real;
    planes[1] = imag;
    cv::merge(planes, 2, complexI);

    //calculating the idft
    cv::dft(complexI, dst, cv::DFT_INVERSE|cv::DFT_REAL_OUTPUT);
}

template<>
void rslf::fft_htranslate<cv::Vec3f>(const Mat& src, float tr_amount, Mat& dst)
{
    // Separate the channels
    Vec<Mat> channels;
    cv::split(src, channels);
    
    Vec<Mat> res_channels(channels.size());
    fft_htranslate<float>(channels[0], tr_amount, res_channels[0]);
    fft_htranslate<float>(channels[1], tr_amount, res_channels[1]);
    fft_htranslate<float>(channels[2], tr_amount, res_channels[2]);
    
    cv::merge(res_channels, dst);
}
