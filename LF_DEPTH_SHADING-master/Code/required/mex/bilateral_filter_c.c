#include "mex.h"
#include "math.h"
#include <stdlib.h>

/* Constant Time Bilateral Filter */
// ********************************* USAGE IN MATLAB
// double* output = bilteral_filter_c(double* input, double filter_radius, double sigma_color)
//
// ********************************* EXAMPLE
// IN  = im2double(imread('xx.png'));
// % blur radius is 5; color threshold is 0.5
// OUT = bilateral_filter_c(IN, 5, 0.5);

/* Output Arguments */
#define	Image_OutC          plhs[0]

/* Input Arguments */
#define	Image_InC			prhs[0]
#define	RADIUS      	    prhs[1]
#define	SIGMA_COLOR    	    prhs[2]

/* Intermediate Buffers */
// original image
int         CANVAS_WIDTH             ;
int         CANVAS_HEIGHT            ;
// padding
double *    Image_In_Padded          ;
double *    Image_Out_Padded         ;
int         PAD_WIDTH                ;
int         PAD_HEIGHT               ;
int         image_pad_size           ;

/* Padding Stuff */
//#define ROW_OFFSET          prhs[6]
//#define START_OFFSET        prhs[7]
//#define PAD_HEIGHT          prhs[8]
//#define PAD_WIDTH           prhs[9]

void generate_weights(
                     double *weights,
                     double sigma_color
                     )
    {
    int i;
    if (sigma_color == 0.0f)
        {
        weights[0] = 1;
        for (i = 1; i < 256; ++i)
            weights[i] = 0;
        }
    else
        {
        for (i = 0; i < 256; ++i)
            weights[i] = exp(-((double) i)/(2*sigma_color*sigma_color));
        }
    }
                        

void RefBNBilateralRect(
                       double *Src,
                       double *Dst,
                       int rows,
                       int cols,
                       int radius,
                       int row_offset,
                       double *rangeWeightsPtr
                       )
    {
    //iterators
    int   COL, ROW, w, h, m;
    int   sRowBytes, dRowBytes;
    //summation and normalization
    float summation, normalization, rangeWeight;
    // histogram formation
    int   bins          ; 
    int   elements      ;
    int*  Hists         ;
    int   Hist[256]     ;
    int   val           ;
    
    double* in          ;
    double* out         ;
    
    int*    past        ;
    
    sRowBytes = dRowBytes = row_offset  ;
    
    bins     = 256                      ;
    elements = (cols+2*radius)*bins     ;
    
    Hists = (int *) malloc(elements*sizeof(int));
    memset(Hists, 0, elements*sizeof(int));
    
   
    --cols; //shift by one for easier looping
    
    //this is where the magic happens
    
    /* 1. Use Huang's method to parse through the first row, and save the column histograms */
    memset(Hist, 0 , sizeof(Hist));
    
    for (w = -radius; w <= radius; ++w)
    {
        for (h = -radius; h <= radius; ++h)
        {
            val = Src[ h * sRowBytes + w ];
            //add to current histogram
            Hist[ val ]+=1;
            //add to column histogram
            Hists [ val ]+=1;
        }
        Hists += bins;
    }
    
    for (COL = 0; COL < cols; ++COL)
    {
        summation = normalization = 0;
        val = Src[COL];
        //storing
        for (m = 0; m <bins; ++m)
        {
            if(Hist[m])
            {
                rangeWeight = rangeWeightsPtr[(int)abs(val-m)];
                rangeWeight *= Hist[m];
                summation += m*rangeWeight;
                normalization += rangeWeight;
            }
        }
        //store into destination
        Dst[COL] = (summation/normalization);
        
        //add the next column, subtract the previous column
        for (h=-radius; h <= radius; ++h)
        {
            //subtract
            Hist[(int) Src[ h * sRowBytes - radius + COL] ] -= 1;
            //add
            val = Src[ h * sRowBytes + radius + 1 + COL];
            Hist[ val ] += 1;
            Hists[ val ] += 1;
        }
        Hists += bins;
    }
    
    //the last one - done without shifting
    
    //storing
    val = Src[cols];
    summation = normalization = 0;
    for (m = 0; m <bins; ++m)
    {
        if(Hist[m])
        {
            rangeWeight = rangeWeightsPtr[(int)abs(val-m)];
            rangeWeight *= Hist[m];
            summation += m*rangeWeight;
            normalization += rangeWeight;
        }
    }
    //store into destination
    Dst[cols] = (summation/normalization);
    
    
    out = Src-sRowBytes*radius;
    in = Src+sRowBytes*(radius+1);
    
    //increment
    Src+=sRowBytes;
    Dst+=dRowBytes;
    
    /* 2. Use new algorithm for the rest of the rows */
    for (ROW = 1; ROW< rows; ++ROW)
    {
        //reset
        memset(Hist,0,sizeof(Hist));
        Hists -= elements;
        //first block
        for (w=-radius; w<=radius; ++w)
        {
            //shift this column down
            Hists[(int) out[w]]-=1; //remove top
            Hists[(int) in[w]]+=1; //add bottom
            //adding the first-radius histograms up to radius together to form first
            for (m=0; m<256; ++m)
            {
                if (Hists[m])
                { //if nonzero value
                    Hist[m]+=Hists[m];
                }
            }
            Hists+=bins;
        }
        
        //rest of the blocks
        for (COL = 0; COL < cols; ++COL)
        {
            //reset
            summation = normalization = 0;
            //storing
            val = Src[COL];
            for (m = 0; m <bins; ++m)
            {
                if(Hist[m])
                {
                    rangeWeight = rangeWeightsPtr[(int)abs(val-m)];
                    rangeWeight *= Hist[m];
                    summation += m*rangeWeight;
                    normalization += rangeWeight;
                }
            }
            //store into destination
            Dst[COL] = (summation/normalization);
            
            //add the next column, subtract the previous column
            //shift next column down
            past = Hists - (2*radius+1)*bins;
            Hists[(int) out[COL+radius+1]]-=1;
            Hists[(int) in[COL+radius+1]]+=1;
            for (m=0; m<bins; ++m)
            {
                Hist[m]+=-past[m]+Hists[m];
            }
            Hists+=bins;
            past+=bins;
        }
        
        //last one of the row
        //reset
        summation = normalization = 0;
        //storing
        val = Src[cols];
        for (m = 0; m <bins; ++m)
        {
            if(Hist[m])
            {
                rangeWeight = rangeWeightsPtr[(int)abs(val-m)];
                rangeWeight *= Hist[m];
                summation += m*rangeWeight;
                normalization += rangeWeight;
            }
        }
        //store into destination
        Dst[cols] = (summation/normalization);
        //increment
        Src+=sRowBytes;
        Dst+=dRowBytes;
        out+=sRowBytes;
        in +=sRowBytes;
    }
}

void image2pad_symmetric
       (
       double * image_in,
       double * image_in_pad,
       int      rows,
       int      cols,
       int      radius,
       int      pad_offset
       )
    {
    int w, h;
    
    // left-top corner
    for (h = 0; h < radius; ++h)
        for (w = 0; w < radius; ++w)
            image_in_pad[w + pad_offset*h]  = image_in[radius - w - 1 + cols*(radius - h - 1)];  
    // left-bottom corner
    for (h = rows + radius; h < rows + 2*radius; ++h)
        for (w = 0; w < radius; ++w)
            image_in_pad[w + pad_offset*h]  = image_in[radius - w - 1 + cols*(2*rows + radius - h - 1)];
    // right-top corner
    for (h = 0; h < radius; ++h)
        for (w = cols + radius; w < cols + 2*radius; ++w)
            image_in_pad[w + pad_offset*h]  = image_in[2*cols + radius - w - 1 + cols*(radius - h - 1)];  
    // right-bottom corner
    for (h = rows + radius; h < rows + 2*radius; ++h)
        for (w = cols + radius; w < cols + 2*radius; ++w)
            image_in_pad[w + pad_offset*h]  = image_in[2*cols + radius - w - 1 + cols*(2*rows + radius - h - 1)];
    
    // left
    for (h = radius; h < rows+radius; ++h)
        for (w = 0; w < radius; ++ w)
            image_in_pad[w + pad_offset*h] = image_in[radius - w - 1 + cols*(h-radius)];
    // right
    for (h = radius; h < rows+radius; ++h)
        for (w = cols + radius; w < cols + 2*radius; ++ w)
            image_in_pad[w + pad_offset*h] = image_in[2*cols - w + radius - 1 + cols*(h-radius)];
    // top
    for (h = 0; h < radius; ++h)
        for (w = radius; w < cols+radius; ++ w)
            image_in_pad[w + pad_offset*h] = image_in[(w-radius) + cols*(radius - h - 1)];
    // bottom
    for (h = rows+radius; h < rows + 2*radius; ++h)
        for (w = radius; w < cols+radius; ++ w)
            image_in_pad[w + pad_offset*h] = image_in[(w-radius) + cols*(2*rows + radius - h - 1)];
    
    // body
    for (h = 0; h < rows; ++h)
        for (w = 0; w < cols; ++w)
            image_in_pad[w+radius + pad_offset*(h+radius)] = image_in[w + cols*h];
    
    }

void pad2image
       (
       double* im_out_pad,
       double* im_out,
       int     rows,
       int     cols,
       int     radius,
       int     pad_offset
       )
    {
    int w,h;
    
    for (h = 0; h < rows; ++h)
        for (w = 0; w < cols; ++w)
            im_out[w + cols*h] = im_out_pad[w+radius + pad_offset*(h+radius)];
    
    }

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])

	{ 
    
    /* Output Arguments */
    //#define       Image_OutC          plhs[0]

    /* Input Arguments */
    //#define       Image_InC			prhs[0]
    //#define       RADIUS      	    prhs[1]
    //#define       SIGMA_COLOR    	    prhs[2]

    /* Intermediate Buffers */
    //int           CANVAS_WIDTH            ;
    //int           CANVAS_HEIGHT           ;
    //double *      Image_In_Padded         ;

	double *cnvsout, *cnvsin, *canvasheight, *canvaswidth, *radius, *sigma_color;
    double *pad_height, *pad_width;
    double *row_offset;
    double weights[256];
    int    image_size;
    int    pad_start_offset;
    int    dim[3];

	/* Check for proper number of arguments */

	if (nlhs != 1 && nrhs != 3)  
		mexErrMsgTxt("double* output = bilteral_filter_c(double* input, double filter_radius, double sigma_color)"); 
	
	/* Assign pointers to the various parameters */ 
    
    // Dimensions
    CANVAS_HEIGHT   = (int)      mxGetM(Image_InC)              ;
    CANVAS_WIDTH    = (int)      mxGetN(Image_InC)/3            ;
    dim[0]          =            CANVAS_HEIGHT                  ;
    dim[1]          =            CANVAS_WIDTH                   ;
    dim[2]          =            3                              ;
    image_size      = (CANVAS_WIDTH)*(CANVAS_HEIGHT)            ;
    // Create Output Buffer
    Image_OutC      = mxCreateNumericArray(3,dim, mxDOUBLE_CLASS, mxREAL) ;
    // Extract pointers
	cnvsin          = (double *) mxGetPr(Image_InC)     ;
	cnvsout         = (double *) mxGetPr(Image_OutC)    ;
    radius          = (double *) mxGetPr(RADIUS)        ;
	sigma_color     = (double *) mxGetPr(SIGMA_COLOR)   ;
    // Pad Image_In
    PAD_WIDTH           = CANVAS_WIDTH  + 2*(*radius)                                  ;
    PAD_HEIGHT          = CANVAS_HEIGHT + 2*(*radius)                                  ;
    image_pad_size      = PAD_WIDTH * PAD_HEIGHT                                       ;
    pad_start_offset    = (*radius) + (*radius)*PAD_HEIGHT                             ;
    Image_In_Padded     = (double *) malloc(sizeof(double)*PAD_WIDTH*PAD_HEIGHT)       ;
    Image_Out_Padded    = (double *) malloc(sizeof(double)*PAD_WIDTH*PAD_HEIGHT)       ;
    
    // generate weights
    
    generate_weights
            (
            weights,
            *sigma_color
            );
    
    /* Do the actual computations in a subroutine */
    // CHANNEL 1
    image2pad_symmetric
       (
       cnvsin,
       Image_In_Padded,
       CANVAS_WIDTH,
       CANVAS_HEIGHT,
       *radius,
       PAD_HEIGHT
       );
	RefBNBilateralRect
		(
		Image_In_Padded     + pad_start_offset,
		Image_Out_Padded    + pad_start_offset,
		CANVAS_WIDTH,
		CANVAS_HEIGHT,
        *radius,
        PAD_HEIGHT,
        weights
		);
    pad2image
       (
       Image_Out_Padded,
       cnvsout,
       CANVAS_WIDTH,
       CANVAS_HEIGHT,
       *radius,
       PAD_HEIGHT
       );
   
    // CHANNEL 2
    image2pad_symmetric
       (
       cnvsin+image_size,
       Image_In_Padded,
       CANVAS_WIDTH,
       CANVAS_HEIGHT,
       *radius,
       PAD_HEIGHT
       );
    RefBNBilateralRect
		(
		Image_In_Padded     + pad_start_offset,
		Image_Out_Padded    + pad_start_offset,
		CANVAS_WIDTH,
		CANVAS_HEIGHT,
        *radius,
        PAD_HEIGHT,
        weights
        );
    pad2image
       (
       Image_Out_Padded,
       cnvsout+image_size,
       CANVAS_WIDTH,
       CANVAS_HEIGHT,
       *radius,
       PAD_HEIGHT
      );
    
    // CHANNEL 3
    image2pad_symmetric
       (
       cnvsin+2*image_size,
       Image_In_Padded,
       CANVAS_WIDTH,
       CANVAS_HEIGHT,
       *radius,
       PAD_HEIGHT
       );
    RefBNBilateralRect
		(
		Image_In_Padded     + pad_start_offset,
		Image_Out_Padded    + pad_start_offset,
		CANVAS_WIDTH,
		CANVAS_HEIGHT,
        *radius,
        PAD_HEIGHT,
        weights
        );
    pad2image
       (
       Image_Out_Padded,
       cnvsout+2*image_size,
       CANVAS_WIDTH,
       CANVAS_HEIGHT,
       *radius,
       PAD_HEIGHT
       );
    
    // demalloc
    free(Image_In_Padded);
    free(Image_Out_Padded);
	return;
	}