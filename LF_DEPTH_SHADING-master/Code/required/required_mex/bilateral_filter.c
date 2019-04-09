#include "mex.h"
#include "math.h"
#include <stdlib.h>

/* Input Arguments */
#define	Image_InC			prhs[0]
#define	Image_OutC          prhs[1]
#define	CANVAS_HEIGHT       prhs[2]
#define	CANVAS_WIDTH	    prhs[3]
#define	RADIUS      	    prhs[4]
#define	SIGMA_SPATIAL	    prhs[5]
#define SIMGA_RGB           prhs[6]

// PADDING
int index_y(int y, int height)
	{
	if (0 <= y && y < height)
		return y;
    else if (y < 0)
        return 0;
    else
        return height-1;
	}
int index_x(int x, int width)
	{
	if (0 <= x && x < width)
		return x;
    else if (x < 0)
        return 0;
    else
        return width-1;
    }

// 3D GAUSSIAN
double gaussian_vector_filter
	(
	double x,
	double y,
	double r,
	double g,
	double b,
	double sigma_spatial,
    double sigma_rgb
	)
	{
	return exp(
		-(x*x+y*y)/(2*sigma_spatial*sigma_spatial)
		-(r*r+g*g+b*b)/(2*sigma_rgb*sigma_rgb));
	}

void bilateral_filter
	(
	double* in_array,
	double* out_array,
	unsigned short width,
	unsigned short height,
    unsigned short radius,
    double sigma_spatial,
    double sigma_rgb
	)
    {
    double coef;
    double coef_sum;
    double R_sum;
    double G_sum;
    double B_sum;
    double base_R;
    double base_G;
    double base_B;
    double current_R;
    double current_G;
    double current_B;
	unsigned short x;
	unsigned short y;
    unsigned short i;
	unsigned short j;
	for(x = 0; x < width; ++x)
        {
		for(y = 0; y < height; ++y)
            {
            base_R = in_array[(y+x*height)];
            base_G = in_array[(y+x*height)+1*width*height];
            base_B = in_array[(y+x*height)+2*width*height];
            
            coef_sum    = 0;
            R_sum       = 0;
            G_sum       = 0;
            B_sum       = 0;
            
            for(j = index_y(y-radius,height); j < index_y(y+radius,height); ++j)
                for(i = index_x(x-radius,width); i < index_x(x+radius,width); ++i)
                    {
                    current_R = in_array[(j+i*height)]                 - base_R;
                    current_G = in_array[(j+i*height)+1*width*height]  - base_G;
                    current_B = in_array[(j+i*height)+2*width*height]  - base_B;
                    
                    coef = gaussian_vector_filter
                            (
                            i-x,
                            j-y,
                            current_R,
                            current_G,
                            current_B,
                            sigma_spatial,
                            sigma_rgb
                            );
                    
                    coef_sum += coef;
                    R_sum    += coef*in_array[(j+i*height)];
                    G_sum    += coef*in_array[(j+i*height)+1*width*height];
                    B_sum    += coef*in_array[(j+i*height)+2*width*height];
                    }
            if (coef_sum > 0.01f)
                {
                out_array[(y+x*height)]                 = R_sum/coef_sum;
                out_array[(y+x*height)+1*width*height]  = G_sum/coef_sum;
                out_array[(y+x*height)+2*width*height]  = B_sum/coef_sum;
                }
            else
                {
                out_array[(y+x*height)]                 = in_array[(y+x*height)]                ;
                out_array[(y+x*height)+1*width*height]  = in_array[(y+x*height)+1*width*height] ;
                out_array[(y+x*height)+2*width*height]  = in_array[(y+x*height)+2*width*height] ;
                }
        }   
     
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])

	{ 
	double *cnvsout, *cnvsin, *dimens, *canvasheight, *canvaswidth, *radius, *sigma_spatial, *sigma_rgb;

	/* Check for proper number of arguments */

	if (nrhs != 7)  
		mexErrMsgTxt("Seven input arguments required."); 
	else if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments."); 
 
	
    
	/* Assign pointers to the various parameters */ 

	cnvsin          = (double *) mxGetPr(Image_InC);
	cnvsout         = (double *) mxGetPr(Image_OutC);
	canvasheight    = (double *) mxGetPr(CANVAS_HEIGHT);
	canvaswidth     = (double *) mxGetPr(CANVAS_WIDTH);
    radius          = (double *) mxGetPr(RADIUS);
	sigma_spatial   = (double *) mxGetPr(SIGMA_SPATIAL);
    sigma_rgb       = (double *) mxGetPr(SIMGA_RGB);

	/* Do the actual computations in a subroutine */
	bilateral_filter
		(
		cnvsin,
		cnvsout,
		*canvaswidth,
		*canvasheight,
        *radius,
        *sigma_spatial,
        *sigma_rgb
		);
	return;
	}