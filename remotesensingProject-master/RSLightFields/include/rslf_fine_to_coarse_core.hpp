#ifndef _RSLF_FINE_TO_COARSE_CORE
#define _RSLF_FINE_TO_COARSE_CORE


#include <limits>
#include <iostream>


#include <rslf_depth_computation.hpp>


/*!
 * \file
 * \brief Implement low-level fine-to-coarse functions.
 */


namespace rslf
{


/**
 * \brief Downsample given EPIs by a factor 2 in the spatial dimensions
 * 
 * @param in_epis Input EPIs
 * @param out_epis Output EPIs
 */
void downsample_EPIs
(
    const Vec<Mat>& in_epis, 
    Vec<Mat>& out_epis
);


/**
 * \brief Fuse a pyramid of depths into one depth map and apply a median filter on top of it
 * 
 * @param in_disp_pyr_p_s_v_u A pyramid of disparity maps (finest first)
 * @param in_validity_indicators_p_s_v_u A pyramid of validity indicators for the respective disparity maps
 * @param out_map_s_v_u The output fused map
 * @param out_validity_s_v_u The output disp validity mask
 */
void fuse_disp_maps
(
    const Vec<Vec<Mat >>& in_disp_pyr_p_s_v_u, 
    const Vec<Vec<Mat> >& in_validity_indicators_p_s_v_u, 
    Vec<Mat>& out_map_s_v_u, 
    Vec<Mat>& out_validity_s_v_u
);


}



#endif
