#include "mex.h"
#include "math.h"
#include <stdlib.h>

/* Input Arguments */
#define	Image_InC			prhs[0]
#define	Image_OutC          prhs[1]
#define	CANVAS_HEIGHT       prhs[2]
#define	CANVAS_WIDTH	    prhs[3]
#define	RADIUS      	    prhs[4]
#define	SIGMA_COLOR    	    prhs[5]
/* Padding Stuff */
#define ROW_OFFSET          prhs[6]
#define START_OFFSET        prhs[7]
#define PAD_HEIGHT          prhs[8]
#define PAD_WIDTH           prhs[9]

void generate_weights(
                     double *weights,
                     double sigma_color
                     )
    {
    int i;
    
    for (i = 0; i < 256; ++i)
        {
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


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])

	{ 
    
    /* Input Arguments */
    // #define	Image_InC			prhs[0]
    // #define	Image_OutC          prhs[1]
    // #define	CANVAS_HEIGHT       prhs[2]
    // #define	CANVAS_WIDTH	    prhs[3]
    // #define	RADIUS      	    prhs[4]
    // #define	SIGMA_COLOR    	    prhs[5]
    // #define  ROW_OFFSET          prhs[6]

	double *cnvsout, *cnvsin, *canvasheight, *canvaswidth, *radius, *sigma_color;
    double *pad_height, *pad_width;
    double *row_offset, *start_offset;
    double weights[256];
    int    image_size;
    int    start_offset_const;

	/* Check for proper number of arguments */

	if (nrhs != 10)  
		mexErrMsgTxt("(image_in,image_out,height,width,radius,sigma_color, row_offset, start_offset)"); 
	
	/* Assign pointers to the various parameters */ 

	cnvsin          = (double *) mxGetPr(Image_InC);
	cnvsout         = (double *) mxGetPr(Image_OutC);
	canvasheight    = (double *) mxGetPr(CANVAS_HEIGHT);
	canvaswidth     = (double *) mxGetPr(CANVAS_WIDTH);
    radius          = (double *) mxGetPr(RADIUS);
	sigma_color     = (double *) mxGetPr(SIGMA_COLOR);
    row_offset      = (double *) mxGetPr(ROW_OFFSET);
    start_offset    = (double *) mxGetPr(START_OFFSET);
    pad_height      = (double *) mxGetPr(PAD_HEIGHT);
    pad_width       = (double *) mxGetPr(PAD_WIDTH);
    
    image_size = (*pad_height)*(*pad_width);
    
    start_offset_const = *start_offset;
    
    generate_weights
            (
            weights,
            *sigma_color
            );

	/* Do the actual computations in a subroutine */
	RefBNBilateralRect
		(
		cnvsin + start_offset_const,
		cnvsout + start_offset_const,
		*canvaswidth,
		*canvasheight,
        *radius,
        *row_offset,
        weights
		);
    
    RefBNBilateralRect
		(
		cnvsin+image_size+ start_offset_const,
		cnvsout+image_size+ start_offset_const,
		*canvaswidth,
		*canvasheight,
        *radius,
        *row_offset,
        weights
        );
    RefBNBilateralRect
		(
		cnvsin+2*image_size+ start_offset_const,
		cnvsout+2*image_size+ start_offset_const,
		*canvaswidth,
		*canvasheight,
        *radius,
        *row_offset,
        weights
        );
	return;
	}