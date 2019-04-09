#include "mex.h"
#include "math.h"
#include <stdlib.h>

/* Input Arguments */
#define	X_size_fin			prhs[0]
#define	Y_size_fin          prhs[1]
#define	UV_diameter_fin     prhs[2]
// buffer
#define	Im_in_remap         prhs[3]

/* Output Arguments */
// buffer
#define	Output_image        prhs[4]

void blockavg
        (
        double * im_in_remap,
        double * output_image,
        unsigned short width,
        unsigned short height,
        unsigned short window_side
        )
    {
    int                 x,y                                             ;
    int                 i,j                                             ;
    int                 x_index,y_index                                 ;
    double              output_color                                    ;
    unsigned int        height_of_remap, width_of_remap, pixels_of_remap;
    int                 window_size                                     ;
    
    window_size = window_side*window_side               ;
    
    height_of_remap = height*window_side                ;
    width_of_remap  = width*window_side                 ;
    pixels_of_remap = height_of_remap*width_of_remap    ;
    
    for (x = 0; x < width; ++x)
        for (y = 0; y < height; ++y)
        {
            output_color = 0;

            for (i = 0; i < window_side; ++i)
                for (j = 0; j < window_side; ++j)
                {

                    x_index = (x * window_side) + i;
                    y_index = (y * window_side) + j;

                    output_color = output_color + im_in_remap[y_index + (x_index * height_of_remap)];
            
                }

            output_image[y + (x * height)] = output_color/window_size;
        
        }
    }

// /* Input Arguments */
// #define	X_size_fin			prhs[0]
// #define	Y_size_fin          prhs[1]
// #define	UV_diameter_fin     prhs[2]
// // buffer
// #define	Im_in_remap         prhs[3]

// /* Output Arguments */
// // buffer
// #define	Output_image        prhs[4]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])

	{ 
	double *x_size_fin_pt, *y_size_fin_pt, *uv_diameter_fin_pt, *im_in_remap_pt, *output_image_pt;

	/* Check for proper number of arguments */

	if (nrhs != 5)  
		mexErrMsgTxt("Five input arguments required."); 
	else if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments."); 
 
	
    
	/* Assign pointers to the various parameters */ 

	x_size_fin_pt           = (double *) mxGetPr(X_size_fin)        ;
	y_size_fin_pt           = (double *) mxGetPr(Y_size_fin)        ;
	uv_diameter_fin_pt      = (double *) mxGetPr(UV_diameter_fin)   ;
    im_in_remap_pt          = (double *) mxGetPr(Im_in_remap)       ;
    output_image_pt         = (double *) mxGetPr(Output_image)      ;

	/* Do the actual computations in a subroutine */
    blockavg
        (
        im_in_remap_pt,
        output_image_pt,
        *x_size_fin_pt,
        *y_size_fin_pt,
        *uv_diameter_fin_pt
        );
            
    
	return;
	}