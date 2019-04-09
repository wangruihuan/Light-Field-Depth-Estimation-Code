#ifndef _RSLF_TYPES
#define _RSLF_TYPES    

//~ #include <vector>
//~ #include <string>

#include <opencv2/core/core.hpp>


/*! 
 * \file
 * \brief Implement a toolbox of utilitary functions. 
 * */ 


namespace rslf
{

/*
 * Aliases
 */
 
/**
 * \brief RSLF Matrix class (cv::Mat)
 */
using Mat = cv::Mat;

/**
 * \brief RSLF Vector class (std::vector)
 */
template<typename T>
using Vec = std::vector<T>;

/**
 * \brief Get an explicit form of an OpenCV type.
 * 
 * @param type The type to disambiguate.
 * @return A string describing the type.
 */
std::string type2str(int type);

/**
 * \brief Implement NaN
 */
template<typename DataType>
DataType nan_type();

/**
 * \brief Check whether the element is NaN
 */
template<typename DataType>
bool is_nan_type(DataType x);

template<>
bool is_nan_type<float>(float x);

template<>
bool is_nan_type<cv::Vec3f>(cv::Vec3f x);

/**
 * \brief Impelment zero
 */
template<typename DataType>
cv::Scalar zero_scalar();

template<>
cv::Scalar zero_scalar<float>();

template<>
cv::Scalar zero_scalar<cv::Vec3f>();


/**
 * \brief "No array" static value
 */
static Mat no_array = Mat();

/**
 * \brief "No Vec<Mat>" static value
 */
static Vec<Mat> no_vec = Vec<Mat>();

/**
 * \brief Implement the euclidean norm
 */
template<typename DataType>
float norm(DataType x);

/**
 * \brief Elementwise sinus
 */
template<typename DataType>
void sin_mat(const Mat& src, Mat& dst);

/**
 * \brief Elementwise cosinus
 */
template<typename DataType>
void cos_mat(const Mat& src, Mat& dst);

/**
 * \brief Performs I(x) <- I(x-t)
 */
template<typename DataType>
void fft_htranslate(const Mat& src, float tr_amount, Mat& dst);

}

#endif
