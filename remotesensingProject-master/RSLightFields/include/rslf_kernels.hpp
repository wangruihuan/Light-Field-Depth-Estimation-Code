#ifndef _RSLF_KERNELS
#define _RSLF_KERNELS


#include <rslf_types.hpp>


namespace rslf
{

/*
 * *****************************************************************
 * KERNEL CLASSES
 * *****************************************************************
 */

/**
 * \brief Pure virtual class implementing a generic kernel instance.
 */
template<typename DataType>
class KernelClass
{
    public:
        /**
         * \brief Get the value of the kernel.
         */
        virtual float evaluate(DataType x) = 0;
        /**
         * \brief Get a matrix with kernel evaluated at each element.
         */
        virtual void evaluate_mat(const Mat& src, Mat& dst) = 0;
};

/**
 * \brief This kernel returns value:
 * 1 - norm(x/h)^2 if norm(x/h) < 1,
 * 0 else
 */
template<typename DataType>
class BandwidthKernel: public KernelClass<DataType>
{
    public:
        BandwidthKernel(float h): m_h_(h) { inv_m_h_sq = 1.0 / (m_h_ * m_h_); }
        float evaluate(DataType x);
        void evaluate_mat(const Mat& src, Mat& dst);

    private:
        float m_h_;
        float inv_m_h_sq;
};

}


#endif
