#include <rslf_io.hpp>
#include <rslf_fine_to_coarse.hpp>


/*!
 * \file
 * \brief Including this header will load all the functions and classes.
 */


/*!
 * \mainpage RSLightFields
 * 
 * \author Quentin CHAN-WAI-NAM (14chanwa)
 * 
 * 
 * \section Introduction
 * 
 * This program implements a disparity computation method from 3D light fields
 * (images taken along a linear path), originally presented by Kim et al. in
 * Scene reconstruction from high spatio-angular resolution light fields. ACM Trans. Graph, July 2017.
 * 
 * 
 * It was made in the framework of the course ["Remote sensing: from sensor to large-scale
 * geospatial data exploitation"](https://mvaisat.wp.imt.fr/) of the Master 2 MVA 
 * (ENS Paris-Saclay, Telecom ParisTech).
 * 
 * 
 * It is distributed under the GNU-GPLv3 license. The corresponding GitHub can be found here:
 * https://github.com/14chanwa/remotesensingProject
 * 
 * 
 * \section Dependencies
 * 
 * * OpenCV 3.x
 * * OpenMP
 * * C++11 with support for stdc++fs (experimental/filesystem)
 * 
 * 
 * \section Structure
 * 
 * The project is organized under the form of high-level classes (implemented in the files
 * `rslf_depth_computation.hpp` and `rslf_fine_to_coarse.hpp`, and low-level functions
 * that are defined in the other headers (mainly `rslf_depth_computation_core.hpp` and 
 * `rslf_fine_to_coarse_core.hpp`). The documentation for the latter can be found
 * in the documentation of the files in which they are defined. 
 * 
 */
