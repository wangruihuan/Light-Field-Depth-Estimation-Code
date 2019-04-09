#ifndef _RSLF_DEPTH_COMPUTATION
#define _RSLF_DEPTH_COMPUTATION

#include <rslf_plot.hpp>
#include <rslf_depth_computation_core.hpp>


/*!
 * \file 
 * \brief Implement high level classes for disparity computation.
 */


namespace rslf
{

/*
 * *****************************************************************
 * Depth1DComputer
 * *****************************************************************
 */

/**
 * \brief Template class with depth computation using 1d slices of the EPI.
 */
template<typename DataType>
class Depth1DComputer
{
public:
    Depth1DComputer
    (
        const Mat& epi, 
        float dmin,
        float dmax,
        int dim_d,
        int s_hat = -1, // default s_hat will be s_max / 2,
        float epi_scale_factor = -1,
        const Depth1DParameters<DataType>& parameters = Depth1DParameters<DataType>::get_default()
    );
    
    /**
     * \brief Runs the algorithm
     */
    void run();
    
    /**
     * \brief Gets an EPI with colors corresponding to the computed slopes.
     */
    Mat get_coloured_epi(int a_cv_colormap = cv::COLORMAP_JET);

private:
    Mat m_epi;

    int m_dim_d;
    Mat m_dmin_u;
    Mat m_dmax_u;
    
    Mat m_edge_confidence_u;
    Mat m_edge_confidence_mask_u;
    Mat m_disp_confidence_u;
    Mat m_rbar_u;
    Mat m_best_depth_u;
    
    /**
     * Line on which to compute the depth
     */
    int m_s_hat;
    
    const Depth1DParameters<DataType>& m_parameters;
};


/**
 * \brief Specialization of Depth1DComputer for 1-channel matrices
 */
using Depth1DComputer_1ch = Depth1DComputer<float>;

/**
 * \brief Specialization of Depth1DComputer for 3-channel matrices
 */
using Depth1DComputer_3ch = Depth1DComputer<cv::Vec3f>;


/*
 * *****************************************************************
 * Depth1DComputer_pile
 * *****************************************************************
 */

/**
 * Template class with depth computation using 1d slices of a pile of EPIs.
 */
template<typename DataType>
class Depth1DComputer_pile
{
public:
    Depth1DComputer_pile
    (
        const Vec<Mat>& epis, 
        float dmin,
        float dmax,
        int dim_d,
        int s_hat = -1, // default s_hat will be s_max / 2,
        float epi_scale_factor = -1,
        const Depth1DParameters<DataType>& parameters = Depth1DParameters<DataType>::get_default()
  );
    
    /**
     * \brief Runs the algorithm
     */
    void run();
    
    /**
     * \brief Gets an EPI with colors corresponding to the computed slopes.
     */
    Mat get_coloured_epi(int a_v = -1, int a_cv_colormap = cv::COLORMAP_JET);
    
    /**
     * \brief Gets a disparity map with colors corresponding to the computed disparity.
     */
    Mat get_disparity_map(int a_cv_colormap = cv::COLORMAP_JET);
    int get_s_hat() { return m_s_hat; }

private:
    Vec<Mat> m_epis;

    int m_dim_d;
    Mat m_dmin_v_u;
    Mat m_dmax_v_u;
    
    Mat m_edge_confidence_v_u;
    Mat m_edge_confidence_mask_v_u;
    Mat m_disp_confidence_v_u;
    Mat m_rbar_v_u;
    Mat m_best_depth_v_u;
    
    /**
     * Line on which to compute the depth
     */
    int m_s_hat;
    
    const Depth1DParameters<DataType>& m_parameters;
};


/**
 * \brief Specialization of Depth1DComputer_pile for 1-channel matrices
 */
using Depth1DComputer_pile_1ch = Depth1DComputer_pile<float>;

/**
 * \brief Specialization of Depth1DComputer_pile for 3-channel matrices
 */
using Depth1DComputer_pile_3ch = Depth1DComputer_pile<cv::Vec3f>;


/*
 * *****************************************************************
 * Depth2DComputer
 * *****************************************************************
 */

/**
 * Template class with depth computation using a pile of EPIs.
 */
template<typename DataType>
class Depth2DComputer
{
public:
    Depth2DComputer
    (
        const Vec<Mat>& epis, 
        float dmin,
        float dmax,
        int dim_d,
        float epi_scale_factor = -1,
        const Depth1DParameters<DataType>& parameters = Depth1DParameters<DataType>::get_default(),
        bool verbose = true
  );
    
    /**
     * \brief Runs the algorithm
     */
    void run();
    
    /**
     * \brief Gets an EPI with colors corresponding to the computed slopes.
     */
    Mat get_coloured_epi(int a_v = -1, int a_cv_colormap = cv::COLORMAP_JET);
    
    /**
     * \brief Gets a disparity map with colors corresponding to the computed disparity.
     */
    Mat get_disparity_map(int a_s = -1, int a_cv_colormap = cv::COLORMAP_JET);
    
    const Vec<Mat>& get_depths_s_v_u() { return m_best_depth_s_v_u; }
    
    const Vec<Mat>& get_valid_depths_mask_s_v_u();

    const Vec<Mat>& get_epis() { return m_epis; }
    
    Vec<Mat>& edit_dmin() { return m_dmin_s_v_u; }
    
    Vec<Mat>& edit_dmax() { return m_dmax_s_v_u; }
    
    void set_accept_all(bool b) { m_accept_all = b; }
    
    Vec<Mat> m_edge_confidence_s_v_u;
    Vec<Mat> m_edge_confidence_mask_s_v_u;
    Vec<Mat> m_disp_confidence_s_v_u;
#ifdef _USE_LINE_CONFIDENCE_SCORE
    Vec<Mat> m_line_confidence_s_v_u;
#endif
    Vec<Mat> m_rbar_s_v_u;
    Vec<Mat> m_best_depth_s_v_u;

private:
    Vec<Mat> m_epis;

    int m_dim_d;
    Vec<Mat> m_dmin_s_v_u;
    Vec<Mat> m_dmax_s_v_u;
    
    bool m_accept_all = false;
    
    const Depth1DParameters<DataType>& m_parameters;
    
    bool m_verbose;
};


/**
 * \brief Specialization of Depth1DComputer for 1-channel matrices
 */
using Depth2DComputer_1ch = Depth2DComputer<float>;

/**
 * \brief Specialization of Depth1DComputer for 3-channel matrices
 */
using Depth2DComputer_3ch = Depth2DComputer<cv::Vec3f>;



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



/*
 * *****************************************************************
 * IMPLEMENTATION
 * Depth1DComputer
 * *****************************************************************
 */

template<typename DataType>
Depth1DComputer<DataType>::Depth1DComputer
(
    const Mat& epi, 
    float dmin,
    float dmax,
    int dim_d,
    int s_hat,
    float epi_scale_factor,
    const Depth1DParameters<DataType>& parameters
) : 
    m_parameters(parameters)
{
    // If the input epi is a uchar, scale uchar to 1.0
    if (epi.depth() == CV_8U)
    {
        epi.convertTo(m_epi, CV_32F, 1.0/255.0);
    }
    else
    {
        // If provided scale factor is invalid, scale from max in all channels
        if (epi_scale_factor < 0)
        {
            Vec<Mat> channels;
            cv::split(epi, channels);
            for (int c=0; c<channels.size(); c++) 
            {
                double min, max;
                cv::minMaxLoc(epi, &min, &max);
                epi_scale_factor = std::max((float)max, epi_scale_factor);
            }
        }
        epi.convertTo(m_epi, CV_32F, 1.0/epi_scale_factor);
    }
    
    
    // Dimensions
    int dim_s = m_epi.rows;
    int dim_u = m_epi.cols;
    
    // d's
    m_dim_d = dim_d;
    m_dmin_u = Mat(1, dim_u, CV_32FC1, cv::Scalar(dmin));
    m_dmax_u = Mat(1, dim_u, CV_32FC1, cv::Scalar(dmax));
    
    // s_hat
    if (s_hat < 0 || s_hat > dim_s - 1) 
    {
        // Default: s_hat is the center horizontal line index of the epi
        m_s_hat = (int) std::floor((0.0 + dim_s) / 2);
    }
    else
    {
        m_s_hat = s_hat;
    }
    
    // Edge confidence
    m_edge_confidence_u = cv::Mat::zeros(1, m_epi.cols, CV_32FC1);
    
    // Disparity confidence
    m_disp_confidence_u = cv::Mat::zeros(1, m_epi.cols, CV_32FC1);
    
    // Scores and best scores & depths
    m_best_depth_u = cv::Mat::zeros(1, dim_u, CV_32FC1);
    
    // rbar
    m_rbar_u = cv::Mat::zeros(1, dim_u, m_epi.type());
}

template<typename DataType>
void Depth1DComputer<DataType>::run() 
{
    // Dimension
    int dim_s = m_epi.rows;
    int dim_u = m_epi.cols;
    
    std::cout << "Max num of threads: " << omp_get_max_threads() << std::endl;
    //~ omp_set_nested(1);
    
    // Empty mask -> no masked point
    Mat mask;
    
    // Buffer
    BufferDepth1D<DataType> buffer(
        dim_s,
        m_dim_d, 
        dim_u, 
        m_epi.type(),
        m_parameters
    );
    
    compute_1D_edge_confidence<DataType>(
        m_epi,
        m_s_hat,
        m_edge_confidence_u,
        m_edge_confidence_mask_u,
        m_parameters,
        buffer
    );
    
    compute_1D_depth_epi<DataType>(
        m_epi,
        m_dmin_u,
        m_dmax_u,
        m_dim_d,
        m_s_hat,
        m_edge_confidence_u,
        m_edge_confidence_mask_u,
        m_disp_confidence_u,
        m_best_depth_u,
        m_rbar_u,
        m_parameters,
        buffer,
        mask
    );
}

template<typename DataType>
Mat Depth1DComputer<DataType>::get_coloured_epi(int a_cv_colormap) {
    
    // Dimensions
    int dim_s = m_epi.rows;
    int dim_u = m_epi.cols;
    
    // Build a matrix of occlusions: each element is the max observed depth
    Mat occlusion_map(dim_s, dim_u, CV_32FC1, -std::numeric_limits<float>::infinity());
    
    // Build a correspondance depth->color: scale to uchar and map to 3-channel matrix
    Mat coloured_depth = rslf::copy_and_scale_uchar(m_best_depth_u);
    cv::applyColorMap(coloured_depth.clone(), coloured_depth, a_cv_colormap);
    
    // Construct an EPI with overlay
    Mat coloured_epi = cv::Mat::zeros(m_epi.rows, m_epi.cols, CV_8UC3);
    
    // For each column of the s_hat row, draw the line, taking overlays into account
    for (int u=0; u<dim_u; u++)
    {
        // Only paint if the confidence threshold was high enough
        if (m_edge_confidence_mask_u.at<uchar>(u))
        {
            float current_depth_value = m_best_depth_u.at<float>(u);
            for (int s=0; s<dim_s; s++)
            {
            
                int requested_index = u + (int)std::round(m_best_depth_u.at<float>(u) * (m_s_hat - s));
                if 
                (
                    requested_index > 0 && 
                    requested_index < dim_u && 
                    occlusion_map.at<float>(s, requested_index) < current_depth_value // only draw if the current depth is higher
                )
                {
                    coloured_epi.at<cv::Vec3b>(s, requested_index) = coloured_depth.at<cv::Vec3b>(u);
                    occlusion_map.at<float>(s, requested_index) = current_depth_value;
                }
            }
        }
    }
    
    return coloured_epi;
}

/*
 * *****************************************************************
 * IMPLEMENTATION
 * Depth1DComputer_pile
 * *****************************************************************
 */

template<typename DataType>
Depth1DComputer_pile<DataType>::Depth1DComputer_pile
(
    const Vec<Mat>& epis, 
    float dmin,
    float dmax,
    int dim_d,
    int s_hat,
    float epi_scale_factor,
    const Depth1DParameters<DataType>& parameters
) : 
    m_parameters(parameters)
{
    m_epis = Vec<Mat>(epis.size(), Mat(epis[0].rows, epis[0].cols, epis[0].type()));

    // Look for the max and min values across all epis
    // If provided scale factor is invalid, scale from max in all channels
    if (epis[0].depth() != CV_8U && epi_scale_factor < 0)
    {
#pragma omp parallel for
        for (int v=0; v<epis.size(); v++)
        {
            Mat epi = epis[v];
            Vec<Mat> channels;
            cv::split(epi, channels);
            for (int c=0; c<channels.size(); c++) 
            {
                double min, max;
                cv::minMaxLoc(epi, &min, &max);
#pragma omp critical
{
                epi_scale_factor = std::max((float)max, epi_scale_factor);
}
            }
        }
    }
    
    // If the input epi is a uchar, scale uchar to 1.0
#pragma omp parallel for
    for (int v=0; v<epis.size(); v++) 
    {
        Mat epi = epis[v];
        Mat epi2;
        if (epi.depth() == CV_8U)
        {
            epi.convertTo(epi2, CV_32F, 1.0/255.0);
        }
        else
        {
            epi.convertTo(epi2, CV_32F, 1.0/epi_scale_factor);
        }
        m_epis[v] = epi2;
    }
    
    // Dimensions
    int dim_s = m_epis[0].rows;
    int dim_u = m_epis[0].cols;
    int dim_v = m_epis.size();
    
    // d's
    m_dim_d = dim_d;
    m_dmin_v_u = Mat(dim_v, dim_u, CV_32FC1, cv::Scalar(dmin));
    m_dmax_v_u = Mat(dim_v, dim_u, CV_32FC1, cv::Scalar(dmax));
    
    // s_hat
    if (s_hat < 0 || s_hat > dim_s - 1) 
    {
        // Default: s_hat is the center horizontal line index of the epi
        m_s_hat = (int) std::floor((0.0 + dim_s) / 2);
    }
    else
    {
        m_s_hat = s_hat;
    }
    
    // Edge confidence
    m_edge_confidence_v_u = Mat(dim_v, dim_u, CV_32FC1);
    
    // Disparity confidence
    m_disp_confidence_v_u = Mat(dim_v, dim_u, CV_32FC1);
    
    // Scores and best scores & depths
    m_best_depth_v_u = cv::Mat::zeros(dim_v, dim_u, CV_32FC1);
    
    // rbar
    m_rbar_v_u = cv::Mat::zeros(dim_v, dim_u, m_epis[0].type());
}

template<typename DataType>
void Depth1DComputer_pile<DataType>::run() 
{
    // No mask
    Mat mask_v_u;

    // Dimension
    int dim_s = m_epis[0].rows;
    int dim_v = m_epis.size();
    int dim_u = m_epis[0].cols;

    int thr_max = omp_get_max_threads();
    std::cout << "Max num of threads: " << thr_max << std::endl;
    
    Vec<BufferDepth1D<DataType>*> buffers;
    for (int t=0; t<thr_max; t++)
        buffers.push_back(new BufferDepth1D<DataType>(
            dim_s,
            m_dim_d, 
            dim_u, 
            m_epis[0].type(),
            m_parameters
            )
        );

    compute_1D_edge_confidence_pile<DataType>(
        m_epis,
        m_s_hat,
        m_edge_confidence_v_u,
        m_edge_confidence_mask_v_u,
        m_parameters,
        buffers
    );

    compute_1D_depth_epi_pile<DataType>(
        m_epis,
        m_dmin_v_u,
        m_dmax_v_u,
        m_dim_d,
        m_s_hat,
        m_edge_confidence_v_u,
        m_edge_confidence_mask_v_u,
        m_disp_confidence_v_u,
        m_best_depth_v_u,
        m_rbar_v_u,
        m_parameters,
        buffers,
        mask_v_u
    );

    for (int t=0; t<thr_max; t++)
        delete buffers[t];
}

template<typename DataType>
Mat Depth1DComputer_pile<DataType>::get_coloured_epi(int a_v, int a_cv_colormap) {
    
    // Dimensions
    int dim_s = m_epis[0].rows;
    int dim_u = m_epis[0].cols;
    int dim_v = m_epis.size();
    
    if (a_v < 0)
        a_v = (int)std::floor(dim_v/2.0);
    
    Mat epi = m_epis[a_v];
    Mat best_depth_u = m_best_depth_v_u.row(a_v);
    Mat edge_confidence_u = m_edge_confidence_v_u.row(a_v);
    Mat edge_confidence_mask_u = m_edge_confidence_mask_v_u.row(a_v);
    
    // Build a matrix of occlusions: each element is the max observed depth
    Mat occlusion_map(dim_s, dim_u, CV_32FC1, -std::numeric_limits<float>::infinity());
    
    // Build a correspondance depth->color: scale to uchar and map to 3-channel matrix
    Mat coloured_depth = rslf::copy_and_scale_uchar(best_depth_u);
    cv::applyColorMap(coloured_depth.clone(), coloured_depth, a_cv_colormap);
    
    // Construct an EPI with overlay
    Mat coloured_epi = cv::Mat::zeros(epi.rows, epi.cols, CV_8UC3);
    
    // For each column of the s_hat row, draw the line, taking overlays into account
#pragma omp parallel for
    for (int u=0; u<dim_u; u++)
    {
        // Only paint if the confidence threshold was high enough
        if (edge_confidence_mask_u.at<uchar>(u))
        {
            float current_depth_value = best_depth_u.at<float>(u);
            for (int s=0; s<dim_s; s++)
            {
            
                int requested_index = u + (int)std::round(best_depth_u.at<float>(u) * (m_s_hat - s));
                if 
                (
                    requested_index > -1 && 
                    requested_index < dim_u && 
                    occlusion_map.at<float>(s, requested_index) < current_depth_value // only draw if the current depth is higher
                )
                {
                    coloured_epi.at<cv::Vec3b>(s, requested_index) = coloured_depth.at<cv::Vec3b>(u);
                    occlusion_map.at<float>(s, requested_index) = current_depth_value;
                }
            }
        }
    }
    
    return coloured_epi;
}

template<typename DataType>
Mat Depth1DComputer_pile<DataType>::get_disparity_map(int a_cv_colormap)
{
    // Dimensions
    int dim_s = m_epis[0].rows;
    int dim_u = m_epis[0].cols;
    int dim_v = m_epis.size();
    
    Mat disparity_map;
    
    disparity_map = rslf::copy_and_scale_uchar(m_best_depth_v_u);
    cv::applyColorMap(disparity_map, disparity_map, a_cv_colormap);
    
    // Threshold scores
    Mat disparity_map_with_scores = cv::Mat::zeros(dim_v, dim_u, disparity_map.type());
    
    cv::add(disparity_map, disparity_map_with_scores, disparity_map_with_scores, m_edge_confidence_mask_v_u);
    
    return disparity_map_with_scores;
}


/*
 * *****************************************************************
 * IMPLEMENTATION
 * Depth2DComputer
 * *****************************************************************
 */

template<typename DataType>
Depth2DComputer<DataType>::Depth2DComputer
(
    const Vec<Mat>& epis, 
    float dmin,
    float dmax,
    int dim_d,
    float epi_scale_factor,
    const Depth1DParameters<DataType>& parameters,
    bool verbose
) : 
    m_parameters(parameters),
    m_verbose(verbose)
{
    m_epis = Vec<Mat>(epis.size(), Mat(epis[0].rows, epis[0].cols, epis[0].type()));

    // Look for the max and min values across all epis
    // If provided scale factor is invalid, scale from max in all channels
    if (epis[0].depth() != CV_8U && epi_scale_factor < 0)
    {
#pragma omp parallel for
        for (int v=0; v<epis.size(); v++)
        {
            Mat epi = epis[v];
            Vec<Mat> channels;
            cv::split(epi, channels);
            for (int c=0; c<channels.size(); c++) 
            {
                double min, max;
                cv::minMaxLoc(epi, &min, &max);
#pragma omp critical
{
                epi_scale_factor = std::max((float)max, epi_scale_factor);
}
            }
        }
    }
    
    // If the input epi is a uchar, scale uchar to 1.0
#pragma omp parallel for
    for (int v=0; v<epis.size(); v++) 
    {
        Mat epi = epis[v];
        Mat epi2;
        if (epi.depth() == CV_8U)
        {
            epi.convertTo(epi2, CV_32F, 1.0/255.0);
        }
        else
        {
            epi.convertTo(epi2, CV_32F, 1.0/epi_scale_factor);
        }
        m_epis[v] = epi2;
    }
    
    // Dimensions
    int dim_s = m_epis[0].rows;
    int dim_u = m_epis[0].cols;
    int dim_v = m_epis.size();
    
    // d's
    m_dim_d = dim_d;
    for (int s=0; s<dim_s; s++)
    {
        m_dmin_s_v_u.push_back(Mat(dim_v, dim_u, CV_32FC1, cv::Scalar(dmin)));
        m_dmax_s_v_u.push_back(Mat(dim_v, dim_u, CV_32FC1, cv::Scalar(dmax)));
    }
    
    m_edge_confidence_s_v_u = Vec<Mat>(dim_s);
    m_disp_confidence_s_v_u = Vec<Mat>(dim_s);
#ifdef _USE_LINE_CONFIDENCE_SCORE
    m_line_confidence_s_v_u = Vec<Mat>(dim_s);
#endif
    m_best_depth_s_v_u = Vec<Mat>(dim_s);
    m_rbar_s_v_u = Vec<Mat>(dim_s);

    for (int s=0; s<dim_s; s++)
    {
        // Edge confidence
        m_edge_confidence_s_v_u[s] = Mat(dim_v, dim_u, CV_32FC1);
    
        // Disparity confidence
        m_disp_confidence_s_v_u[s] = Mat(dim_v, dim_u, CV_32FC1);

#ifdef _USE_LINE_CONFIDENCE_SCORE
        // Line confidence
        m_line_confidence_s_v_u[s] = Mat(dim_v, dim_u, CV_32FC1);
#endif
        
        // Scores and best scores & depths
        m_best_depth_s_v_u[s] = cv::Mat::zeros(dim_v, dim_u, CV_32FC1);
        
        // rbar
        m_rbar_s_v_u[s] = cv::Mat::zeros(dim_v, dim_u, m_epis[0].type());
    }
}

template<typename DataType>
void Depth2DComputer<DataType>::run() 
{
    // Dimension
    int dim_s = m_epis[0].rows;
    int dim_v = m_epis.size();
    int dim_u = m_epis[0].cols;

    int thr_max = omp_get_max_threads();
    
    if (m_verbose)
    {
        std::cout << "Max num of threads: " << thr_max << std::endl;
        std::cout << "Slope factor: " << m_parameters.par_slope_factor << std::endl;
    }
    
    Vec<BufferDepth1D<DataType>*> m_buffers_;
    for (int t=0; t<thr_max; t++)
        m_buffers_.push_back(new BufferDepth1D<DataType>(
            dim_s,
            m_dim_d, 
            dim_u, 
            m_epis[0].type(),
            m_parameters
            )
        );

    compute_2D_edge_confidence<DataType>(
        m_epis,
        m_edge_confidence_s_v_u,
        m_edge_confidence_mask_s_v_u,
        m_parameters,
        m_buffers_
    );

    compute_2D_depth_epi<DataType>(
        m_epis,
        m_dmin_s_v_u,
        m_dmax_s_v_u,
        m_dim_d,
        m_edge_confidence_s_v_u,
        m_edge_confidence_mask_s_v_u,
        m_disp_confidence_s_v_u,
#ifdef _USE_LINE_CONFIDENCE_SCORE
        m_line_confidence_s_v_u,
#else
        rslf::no_vec,
#endif
        m_best_depth_s_v_u,
        m_rbar_s_v_u,
        m_parameters,
        m_buffers_,
        m_verbose
    );

    for (int t=0; t<thr_max; t++)
        delete m_buffers_[t];
}

template<typename DataType>
Mat Depth2DComputer<DataType>::get_coloured_epi(int a_v, int a_cv_colormap) {
    
    // Dimensions
    int dim_s = m_epis[0].rows;
    int dim_u = m_epis[0].cols;
    int dim_v = m_epis.size();
    
    if (a_v < 0)
        a_v = (int)std::floor(dim_v/2.0);
    
    // Construct an EPI with overlay
    Mat best_depth_s_u = cv::Mat::zeros(dim_s, dim_u, CV_32FC1);
    
    // For each column of the s_hat row, get the corresponding line
#pragma omp parallel for
    for (int s=0; s<dim_s; s++)
    {
        m_best_depth_s_v_u[s].row(a_v).copyTo(best_depth_s_u.row(s));
    }
    
    // Build a correspondance depth->color: scale to uchar and map to 3-channel matrix
    Mat coloured_depth = rslf::copy_and_scale_uchar(best_depth_s_u);
    cv::applyColorMap(coloured_depth.clone(), coloured_depth, a_cv_colormap);
    
    Mat coloured_epi = cv::Mat::zeros(dim_s, dim_u, CV_8UC3);
    for (int s=0; s<dim_s; s++)
    {
        for (int u=0; u<dim_u; u++)
        {
            // Only paint if the confidence threshold was high enough
#ifdef _USE_DISP_CONFIDENCE_SCORE
            Mat disp_confidence_mask_u = m_disp_confidence_s_v_u[s].row(a_v) > m_parameters.par_disp_score_threshold;
            if (disp_confidence_mask_u.at<uchar>(u))
#elseif define(_USE_LINE_CONFIDENCE_SCORE)
            Mat line_confidence_mask_u = m_line_confidence_s_v_u[s].row(a_v) > m_parameters.par_line_score_threshold;
            if (line_confidence_mask_u.at<uchar>(u))
#else
            Mat edge_confidence_mask_u = m_edge_confidence_mask_s_v_u[s].row(a_v);
            if (edge_confidence_mask_u.at<uchar>(u))
#endif
            {
                coloured_epi.at<cv::Vec3b>(s, u) = coloured_depth.at<cv::Vec3b>(s, u);
            }
        }
    }
    
    return coloured_epi;
}

template<typename DataType>
Mat Depth2DComputer<DataType>::get_disparity_map(int a_s, int a_cv_colormap)
{
    // Dimensions
    int dim_s = m_epis[0].rows;
    int dim_u = m_epis[0].cols;
    int dim_v = m_epis.size();
    
    if (a_s < 0)
        a_s = (int)std::floor(dim_s/2.0);
    
    Mat disparity_map;
    Mat best_depth_v_u = m_best_depth_s_v_u[a_s];    
    
    disparity_map = rslf::copy_and_scale_uchar(best_depth_v_u);
    cv::applyColorMap(disparity_map, disparity_map, a_cv_colormap);
    
    // Threshold scores
    Mat disparity_map_with_scores = cv::Mat::zeros(dim_v, dim_u, disparity_map.type());
    
#ifdef _USE_DISP_CONFIDENCE_SCORE
    Mat disp_confidence_v_u = m_disp_confidence_s_v_u[a_s];   
    Mat disp_confidence_mask_v_u = disp_confidence_v_u > m_parameters.par_disp_score_threshold; 
    cv::add(disparity_map, disparity_map_with_scores, disparity_map_with_scores, disp_confidence_mask_v_u);
#elseif defined(_USE_LINE_CONFIDENCE_SCORE)
    Mat line_confidence_v_u = m_line_confidence_s_v_u[a_s];   
    Mat line_confidence_mask_v_u = line_confidence_v_u > m_parameters.par_line_score_threshold; 
    cv::add(disparity_map, disparity_map_with_scores, disparity_map_with_scores, line_confidence_mask_v_u);
#else
    Mat edge_confidence_mask_v_u = m_edge_confidence_mask_s_v_u[a_s];   
    cv::add(disparity_map, disparity_map_with_scores, disparity_map_with_scores, edge_confidence_mask_v_u);
#endif
    
    return disparity_map_with_scores;
}

template<typename DataType>
const Vec<Mat>& Depth2DComputer<DataType>::get_valid_depths_mask_s_v_u()
{
    Vec<Mat>* validity_maps = new Vec<Mat>();
    for (int s=0; s<m_edge_confidence_s_v_u.size(); s++)
    {
        if (!m_accept_all)
        {
#ifdef _USE_DISP_CONFIDENCE_SCORE
            validity_maps->push_back(m_disp_confidence_s_v_u[s] > m_parameters.par_disp_score_threshold);
#elseif defined(_USE_LINE_CONFIDENCE_SCORE)
            validity_maps->push_back(m_line_confidence_s_v_u[s] > m_parameters.par_line_score_threshold);
#else
            validity_maps->push_back(m_edge_confidence_s_v_u[s] > m_parameters.par_edge_score_threshold);
#endif
        }
        else
        {
            validity_maps->push_back(m_edge_confidence_s_v_u[s] > -1);
        }
    }
    return *validity_maps;
}

}


#endif
