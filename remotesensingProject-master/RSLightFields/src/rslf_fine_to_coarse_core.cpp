#include <rslf_fine_to_coarse_core.hpp>


#define _GAUSSIAN_KSIZE 7
#define _FINAL_MEDIAN_FILTER_SIZE 3


/*
 * *****************************************************************
 * IMPLEMENTATION
 * downsample_EPIs
 * *****************************************************************
 */
void rslf::downsample_EPIs(const Vec<Mat>& in_epis, Vec<Mat>& out_epis)
{
    out_epis.clear();
    
    int dim_s = in_epis[0].rows;
    int dim_u = in_epis[0].cols;
    int dim_v = in_epis.size();
    
    int dtype = in_epis[0].type();
    
    // Get a copy of the output epis in the dimensions s, v, u
    Vec<Mat> images_s_v_u_;
    Mat tmp = Mat(dim_v, dim_u, dtype, cv::Scalar(0.0));
    Mat tmp2;
    for (int s=0; s<dim_s; s++)
    {
        Mat tmp3;
        for (int v=0; v<dim_v; v++)
        {
            in_epis[v].row(s).copyTo(tmp.row(v));
        }
        // Filter by a gaussian filter along the spatial dimensions
        // TODO: better border condition useful?
        cv::GaussianBlur(tmp, tmp2, cv::Size(_GAUSSIAN_KSIZE, _GAUSSIAN_KSIZE), 0, 0, cv::BORDER_REFLECT);
        //~ std::cout << "filter ok" << std::endl;
        
        // Subsample by a factor 2 (since we already filtered, nn should be ok)
        cv::resize(tmp2, tmp3, cv::Size(), 0.5, 0.5, cv::INTER_LINEAR);
        //~ std::cout << "resize ok" << std::endl;
            
        // Store the result
        images_s_v_u_.push_back(tmp3);
    }
    
    // Rearrange lines to out_epis
    for (int v=0; v<images_s_v_u_[0].rows; v++)
    {
        Mat tmp = Mat(dim_s, images_s_v_u_[0].cols, dtype, cv::Scalar(0.0));
        for (int s=0; s<dim_s; s++)
        {
            //~ std::cout << "sizes: " << images_s_v_u_[s].row(v).size() << ", " << tmp.row(s).size() << std::endl;
            images_s_v_u_[s].row(v).copyTo(tmp.row(s));
            //~ std::cout << "copy ok" << std::endl;
        }
        out_epis.push_back(tmp);
    }
}


/*
 * *****************************************************************
 * IMPLEMENTATION
 * downsample_EPIs
 * *****************************************************************
 */
void rslf::fuse_disp_maps
(
    const Vec<Vec<Mat >>& in_disp_pyr_p_s_v_u, 
    const Vec<Vec<Mat> >& in_validity_indicators_p_s_v_u, 
    Vec<Mat>& out_map_s_v_u, 
    Vec<Mat>& out_validity_s_v_u)
{
    int dim_u = in_disp_pyr_p_s_v_u[0][0].cols;
    int dim_v = in_disp_pyr_p_s_v_u[0][0].rows;
    int dim_s = in_disp_pyr_p_s_v_u[0].size();
    
    int m_pyr_size = in_disp_pyr_p_s_v_u.size();
    
    out_map_s_v_u.clear();
    out_validity_s_v_u.clear();
    
    for (int s=0; s<dim_s; s++)
    {
        Vec<Mat> disp_pyr_v_u_ = Vec<Mat>();
        Vec<Mat> mask_pyr_v_u_ = Vec<Mat>();
        for (int p=0; p<m_pyr_size; p++)
        {
            disp_pyr_v_u_.push_back(in_disp_pyr_p_s_v_u[p][s]);
            mask_pyr_v_u_.push_back(in_validity_indicators_p_s_v_u[p][s]);
        }
    
        Mat tmp_map_down = disp_pyr_v_u_[m_pyr_size-1];
        Mat tmp_mask_down = mask_pyr_v_u_[m_pyr_size-1];
        Mat tmp_map_up;
        Mat tmp_mask_up;
        
        for (int p = m_pyr_size-1; p > 0; p--)
        {
            
            Mat upscaled_disp;
            
            // Upscale
            // TODO: interpolation type?
            cv::resize(tmp_map_down, tmp_map_up, disp_pyr_v_u_[p-1].size(), 0, 0, cv::INTER_LINEAR);
            //~ cv::resize(tmp_map_down, tmp_map_up, disp_pyr_v_u_[p-1].size(), 0, 0, cv::INTER_NEAREST);
            cv::resize(tmp_mask_down, tmp_mask_up, mask_pyr_v_u_[p-1].size(), 0, 0, cv::INTER_NEAREST);
            
            // Fill unknown points 
            
            // TODO violent... in the article, suggested to put boundaries on downscaled disp calculation
            // did not yet implement this since the structure of the code would be quite modified by this
            
            // One could refine the way the unknown disp are affected, for instance by checking the bounds
            // on the horizontal lines as described in the article
            tmp_map_down = disp_pyr_v_u_[p-1].clone();
            
            Mat tmp_mask = mask_pyr_v_u_[p-1] == 0;
            
            tmp_map_down.setTo(0.0, tmp_mask);
            cv::add(tmp_map_down, tmp_map_up, tmp_map_down, tmp_mask);
            
            cv::bitwise_or(mask_pyr_v_u_[p-1], tmp_mask_up, tmp_mask_down);
        }
        
        // Apply median filter
        Mat tmp;
        cv::medianBlur(tmp_map_down, tmp, _FINAL_MEDIAN_FILTER_SIZE);
        
        out_map_s_v_u.push_back(tmp);
        out_validity_s_v_u.push_back(tmp_mask_down);
    }
}
