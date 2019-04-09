#ifndef _RSLF_FINE_TO_COARSE
#define _RSLF_FINE_TO_COARSE


#include <rslf_fine_to_coarse_core.hpp>


#define _MIN_SPATIAL_DIM 10


/*!
 * \file 
 * \brief Implement high level classes for fine-to-coarse.
 */


namespace rslf
{

/*
 * *****************************************************************
 * FineToCoarse
 * *****************************************************************
 */

template<typename DataType>
class FineToCoarse
{
public:
    FineToCoarse
    (
        const Vec<Mat>& epis, 
        float d_min,
        float d_max,
        int dim_d,
        float epi_scale_factor = -1,
        const Depth1DParameters<DataType>& parameters = Depth1DParameters<DataType>::get_default(),
        int max_pyr_depth = -1,
        bool accept_all_last_scale = true
    );
    ~FineToCoarse();
    
    /**
     * \brief Runs the algorithm
     */
    void run();
    
    /**
     * \brief Get the resulting disparity maps at the finest scale and the corresponding validity mask.
     */
    void get_results(Vec<Mat>& out_map_s_v_u_, Vec<Mat>& out_validity_s_v_u_);
    
    /**
     * \brief Get disparity maps with colors corresponding to the computed disparities.
     */
    void get_coloured_depth_maps(Vec<Mat>& out_plot_depth_s_v_u_, int cv_colormap = cv::COLORMAP_JET, bool saturate = true);
    
    /**
     * \brief Get disparity maps with colors corresponding to the computed disparities, juxtaposed with the original image.
     */
    void get_coloured_depth_maps_and_imgs(Vec<Mat>& out_plot_depth_s_v_u_, int cv_colormap = cv::COLORMAP_JET, bool saturate = true);
    
    /**
     * \brief Get the EPI pyramid with colors corresponding to computed slopes.
     */
    void get_coloured_epi_pyr(Vec<Mat>& out_plot_epi_pyr_p_s_u_, int v = -1, int cv_colormap = cv::COLORMAP_JET, bool saturate = true);
    
    /**
     * \brief Get the disparity map pyramid with colors corresponding to computed disparities.
     */
    void get_coloured_depth_pyr(Vec<Mat>& out_plot_depth_pyr_p_v_u_, int s = -1, int cv_colormap = cv::COLORMAP_JET, bool saturate = true);
    
    //~ Mat get_coloured_epi(int v = -1, int cv_colormap = cv::COLORMAP_JET);
    //~ Mat get_disparity_map(int s = -1, int cv_colormap = cv::COLORMAP_JET);

private:
    Vec<Depth2DComputer<DataType>* > m_computers;
    const Depth1DParameters<DataType>& m_parameters;
    Vec<Depth1DParameters<DataType>* > m_parameter_instances;
    bool m_accept_all_last_scale;
};


/*
 * Aliases
 */

/**
 * \brief Specialization of FineToCoarse for 1-channel matrices
 */
using FineToCoarse_1ch = FineToCoarse<float>;

/**
 * \brief Specialization of FineToCoarse for 3-channel matrices
 */
using FineToCoarse_3ch = FineToCoarse<cv::Vec3f>;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


template<typename DataType>
FineToCoarse<DataType>::FineToCoarse
(
    const Vec<Mat>& epis, 
    float d_min,
    float d_max,
    int dim_d,
    float epi_scale_factor,
    const Depth1DParameters<DataType>& parameters,
    int max_pyr_depth,
    bool accept_all_last_scale
) :
    m_parameters(parameters)
{
    Vec<Mat> tmp_epis = epis;
    
    int start_dim_v = tmp_epis.size();
    int start_dim_u = tmp_epis[0].cols;
    
    int dim_v = tmp_epis.size();
    int dim_u = tmp_epis[0].cols;
    
    if (max_pyr_depth < 1)
        max_pyr_depth = std::numeric_limits<int>::max();
    
    int counter = 0;
    
    while (dim_v > _MIN_SPATIAL_DIM && dim_u > _MIN_SPATIAL_DIM && counter < max_pyr_depth)
    {
        counter++;
        
        std::cout << "Creating Depth2DComputer with sizes (" << dim_v << ", " << dim_u << ")" << std::endl;
        
        Depth1DParameters<DataType>* new_parameters = new Depth1DParameters<DataType>(m_parameters);
        
        // Compute scale factor
        new_parameters->par_slope_factor = (0.0 + dim_u) / start_dim_u;
        
        // Create a new Depth2DComputer
        Depth2DComputer<DataType>* computer = new Depth2DComputer<DataType>(tmp_epis, d_min, d_max, dim_d, epi_scale_factor, *new_parameters);
        
        // Downsample
        Vec<Mat> downsampled_epis;
        downsample_EPIs(tmp_epis, downsampled_epis);
        tmp_epis = downsampled_epis;
        
        dim_v = tmp_epis.size();
        dim_u = tmp_epis[0].cols;
        
        m_computers.push_back(computer);
        m_parameter_instances.push_back(new_parameters);
    }
    
    // The last level should accept all disparity measures
    if (accept_all_last_scale)
        m_computers.back()->set_accept_all(true);
}

template<typename DataType>
FineToCoarse<DataType>::~FineToCoarse()
{
    for (int p=0; p<m_computers.size(); p++)
    {
        delete m_computers[p];
        delete m_parameter_instances[p];
    }
}

template<typename DataType>
void FineToCoarse<DataType>::run()
{
    for (int p=0; p<m_computers.size(); p++)
    {
        m_computers[p]->run();
        if (p < m_computers.size() - 1)
        {
            std::cout << "Setting new depth bounds..." << std::endl;
            
            // Depth map
            Vec<Mat> depth_map_up = m_computers[p]->get_depths_s_v_u();
            Vec<Mat> depth_map_down = m_computers[p+1]->get_depths_s_v_u();
            // Validity mask
            Vec<Mat> mask_up = m_computers[p]->get_valid_depths_mask_s_v_u();
            Vec<Mat> mask_down = m_computers[p+1]->get_valid_depths_mask_s_v_u();
            
            // dmin map
            Vec<Mat> dmin_map = m_computers[p+1]->edit_dmin();
            // dmax map
            Vec<Mat> dmax_map = m_computers[p+1]->edit_dmax();
            
            int dim_s = depth_map_up.size();
            
            int dim_v_up = depth_map_up[0].rows;
            int dim_u_up = depth_map_up[0].cols;
            
            int dim_v_down = depth_map_down[0].rows;
            int dim_u_down = depth_map_down[0].cols;
            
            // Edit the min/max d's
#pragma omp parallel for
            for (int s=0; s<dim_s; s++)
            {
                for (int v=0; v<dim_v_down; v++)
                {
                    for (int u=0; u<dim_u_down; u++)
                    {
                        Vec<float> candidate_ds;
                        // 1st line
                        // Get the upscaled u, v
                        int v_up = std::min(2 * v, dim_v_up - 1);
                        int u_up = std::min(2 * u, dim_u_up - 1);
                        // Get the point at the left
                        int u_left = u_up;
                        float d_left = nan_type<float>();
                        while (u_left > 1)
                        {
                            u_left -= 1;
                            if (mask_up[s].at<uchar>(v_up, u_left) > 0)
                            {
                                d_left = depth_map_up[s].at<float>(v_up, u_left);
                                break;
                            }
                        }
                        int u_right = u_up;
                        float d_right = nan_type<float>();
                        while (u_right < dim_u_up-1)
                        {
                            u_right += 1;
                            if (mask_up[s].at<uchar>(v_up, u_right) > 0)
                            {
                                d_right = depth_map_up[s].at<float>(v_up, u_right);
                                break;
                            }
                        }
                        //~ if (!is_nan_type<float>(d_left) && !is_nan_type<float>(d_right))
                        //~ {
                            //~ dmin_map[s].at<float>(v, u) = std::min(d_left, d_right);
                            //~ dmax_map[s].at<float>(v, u) = std::max(d_left, d_right);//candidate_ds[mid_idx+1];
                        //~ }
                        if (!is_nan_type<float>(d_left) && !is_nan_type<float>(d_right))
                        {
                            candidate_ds.push_back(d_left);
                            candidate_ds.push_back(d_right);
                        }
                        
                        // 2nd line
                        if (v_up+1 < dim_v_up)
                        {
                            // Get the upscaled u, v
                            v_up += 1;
                            int u_up = std::min(2 * u, dim_u_up - 1);
                            // Get the point at the left
                            int u_left = u_up;
                            float d_left = nan_type<float>();
                            while (u_left > 1)
                            {
                                u_left -= 1;
                                if (mask_up[s].at<uchar>(v_up, u_left) > 0)
                                {
                                    d_left = depth_map_up[s].at<float>(v_up, u_left);
                                    break;
                                }
                            }
                            int u_right = u_up;
                            float d_right = nan_type<float>();
                            while (u_right < dim_u_up-1)
                            {
                                u_right += 1;
                                if (mask_up[s].at<uchar>(v_up, u_right) > 0)
                                {
                                    d_right = depth_map_up[s].at<float>(v_up, u_right);
                                    break;
                                }
                            }
                            if (!is_nan_type<float>(d_left) && !is_nan_type<float>(d_right))
                            {
                                candidate_ds.push_back(d_left);
                                candidate_ds.push_back(d_right);
                            }
                        }
                        if (candidate_ds.size() > 1)
                        {
                            std::sort(candidate_ds.begin(), candidate_ds.end());
                            int nb_candidates = candidate_ds.size();
                            //int mid_idx = (int)std::floor((nb_candidates-1.0)/2);
                            // Modify dmin and dmax
                            dmin_map[s].at<float>(v, u) = candidate_ds[0];//candidate_ds[mid_idx];
                            dmax_map[s].at<float>(v, u) = candidate_ds.back();//candidate_ds[mid_idx+1];
                        }
                    }
                }
            }
            
            //~ std::cout << dmax_map[50] << std::endl;
        }
    }
}

template<typename DataType>
void FineToCoarse<DataType>::get_results(Vec<Mat>& out_map_s_v_u_, Vec<Mat>& out_validity_s_v_u_)
{
    std::cout << "Getting results..." << std::endl;
    
    Vec<Vec<Mat >> disp_pyr_p_s_v_u_;
    Vec<Vec<Mat >> validity_indicators_p_s_v_u_;
    
    for (int p=0; p<m_computers.size(); p++)
    {
        disp_pyr_p_s_v_u_.push_back(m_computers[p]->get_depths_s_v_u());
        validity_indicators_p_s_v_u_.push_back(m_computers[p]->get_valid_depths_mask_s_v_u());
    }
    
    fuse_disp_maps
    (
        disp_pyr_p_s_v_u_, 
        validity_indicators_p_s_v_u_, 
        out_map_s_v_u_, 
        out_validity_s_v_u_
    );
}

template<typename DataType>
void FineToCoarse<DataType>::get_coloured_depth_maps(Vec<Mat>& out_plot_depth_s_v_u_, int cv_colormap, bool saturate)
{
    std::cout << "Plot depth results..." << std::endl;
    
    Vec<Mat> out_map_s_v_u_;
    Vec<Mat> out_validity_s_v_u_;
    
    get_results
    (
        out_map_s_v_u_, 
        out_validity_s_v_u_
    );
    
    int dim_s = out_map_s_v_u_.size();
    int dim_v = out_map_s_v_u_[0].rows;
    int dim_u = out_map_s_v_u_[0].cols;
    
    ImageConverter_uchar image_converter;
    image_converter.fit(out_map_s_v_u_[(int)std::round(dim_s/2.0)], saturate);

    for (int s=0; s<dim_s; s++)
    {
        
        Mat disparity_map;
        image_converter.copy_and_scale(out_map_s_v_u_[s], disparity_map);
        cv::applyColorMap(disparity_map, disparity_map, cv_colormap);
        
        int dim_v = disparity_map.rows;
        int dim_u = disparity_map.cols;
        
        // Threshold scores
        disparity_map.setTo(0.0, out_validity_s_v_u_[s] == 0);
        
        // Cut shadows
        if (m_parameters.par_cut_shadows)
        {
            // Get image norm view
            Mat im_norm = Mat(dim_v, dim_u, CV_32FC1);
            for (int v=0; v<dim_v; v++)
            {
                for (int u=0; u<dim_u; u++)
                {
                    Mat tmp = m_computers[0]->get_epis()[v];
                    im_norm.at<float>(v, u) = norm<DataType>(tmp.at<DataType>(s, u));
                }
            }
            disparity_map.setTo(0.0, im_norm < _SHADOW_NORMALIZED_LEVEL);
        }
        
        out_plot_depth_s_v_u_.push_back(disparity_map);
    }

}

template<typename DataType>
void FineToCoarse<DataType>::get_coloured_depth_maps_and_imgs(Vec<Mat>& out_plot_depth_s_v_u_, int cv_colormap, bool saturate)
{
    Vec<Mat> tmp;
    get_coloured_depth_maps(tmp, cv_colormap, saturate);
    Vec<Mat> epis = m_computers[0]->get_epis();
    
    int dim_s = tmp.size();
    int dim_v = tmp[0].rows;
    int dim_u = tmp[0].cols;
    
    out_plot_depth_s_v_u_.clear();
    
    ImageConverter_uchar converter;
    
    for (int s=0; s<dim_s; s++)
    {
        // Get the corresponding image
        Mat img = Mat(dim_v, dim_u, epis[0].type());
        for (int v=0; v<dim_v; v++)
        {
            epis[v].row(s).copyTo(img.row(v));
        }
        
        if (s==0)
            converter.fit(img, false);
        
        converter.copy_and_scale(img, img);
        
        if (img.channels() == 1)
            // Convert to RGB
            cv::cvtColor(img, img, CV_GRAY2RGB);
        
        std::cout << "img: " << img.size() << ", " << rslf::type2str(img.type()) << std::endl;
        std::cout << "tmp[s]: " << tmp[s].size() << ", " << rslf::type2str(tmp[s].type()) << std::endl;
        
        // Concatenate with the depth map
        if (dim_u > dim_v)
        {
            // Concatenate rows
            cv::vconcat(img, tmp[s], img);
        }
        else
        {
            // Concatenate cols
            cv::hconcat(img, tmp[s], img);
        }
        
        out_plot_depth_s_v_u_.push_back(img);
    }
}

template<typename DataType>
void FineToCoarse<DataType>::get_coloured_epi_pyr(Vec<Mat>& out_plot_epi_pyr_p_s_u_, int v, int cv_colormap, bool saturate) 
{
    
    int m_pyr_size_ = m_computers.size();
    int m_dim_v_orig_ = m_computers[0]->get_depths_s_v_u()[0].rows;
    
    if (v == -1)
        v = (int)std::round(m_dim_v_orig_/2.0);
    
    ImageConverter_uchar image_converter;
    
    out_plot_epi_pyr_p_s_u_.clear();
    for (int p=0; p<m_pyr_size_; p++)
    {
        int dim_s = m_computers[p]->get_depths_s_v_u().size();
        int dim_v = m_computers[p]->get_depths_s_v_u()[0].rows;
        int dim_u = m_computers[p]->get_depths_s_v_u()[0].cols;
        
        Mat tmp(dim_s, dim_u, CV_32FC1);
        
        int v_scaled = (int)std::round(1.0 * v * dim_v / m_dim_v_orig_);
        
        Vec<Mat> masks = m_computers[p]->get_valid_depths_mask_s_v_u();
        
        for (int s=0; s<dim_s; s++)
        {
            m_computers[p]->get_depths_s_v_u()[s].row(v_scaled).copyTo(tmp.row(s));
            tmp.row(s).setTo(0.0, masks[s].row(v_scaled) == 0);
        }
        
        if (p == 0)
            image_converter.fit(tmp, saturate);
        
        image_converter.copy_and_scale(tmp, tmp);
        cv::applyColorMap(tmp, tmp, cv_colormap);
        
        // Cut shadows
        if (m_parameters.par_cut_shadows)
        {
            // Get epi norm view
            Mat epi_norm = Mat(dim_s, dim_u, CV_32FC1);
            for (int s=0; s<dim_s; s++)
            {
                for (int u=0; u<dim_u; u++)
                {
                    Mat tmp = m_computers[p]->get_epis()[v_scaled];
                    epi_norm.at<float>(s, u) = norm<DataType>(tmp.at<DataType>(s, u));
                }
            }
            tmp.setTo(cv::Scalar(0.0), epi_norm < _SHADOW_NORMALIZED_LEVEL);
        }
        
        out_plot_epi_pyr_p_s_u_.push_back(tmp);
    }
    
}


template<typename DataType>
void FineToCoarse<DataType>::get_coloured_depth_pyr(Vec<Mat>& out_plot_depth_pyr_p_v_u_, int s, int cv_colormap, bool saturate) {
    
    int m_pyr_size_ = m_computers.size();
    int dim_s = m_computers[0]->get_depths_s_v_u().size();
    
    if (s == -1)
        s = (int)std::round(dim_s/2.0);
    
    ImageConverter_uchar image_converter;
    
    out_plot_depth_pyr_p_v_u_.clear();
    for (int p=0; p<m_pyr_size_; p++)
    {
        Mat tmp = m_computers[p]->get_depths_s_v_u()[s].clone();
        
        if (p == 0)
            image_converter.fit(tmp, saturate);
        
        image_converter.copy_and_scale(tmp, tmp);
        cv::applyColorMap(tmp, tmp, cv_colormap);
        
        Vec<Mat> masks = m_computers[p]->get_valid_depths_mask_s_v_u();
        
        tmp.setTo(0.0, masks[s] == 0);
        
        out_plot_depth_pyr_p_v_u_.push_back(tmp);
    }
}

}

#endif
