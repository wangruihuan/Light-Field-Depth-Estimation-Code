#ifndef _RSLF_PLOT
#define _RSLF_PLOT


#include <rslf_types.hpp>


/*!
 * \file
 * \brief Implement functions to scale and plot matrices.
 */


namespace rslf
{
    
/**
 * \brief Plot the provided cv::Mat in a cv::namedWindow.
 * 
 * @param img The matrix to be plotted.
 * @param window_name Name of the created window.
 */
void plot_mat
(
    Mat img, 
    std::string window_name = "Image"
);

/**
 * \brief Build a copy of the matrix scaled to uchar (0..255).
 * 
 * @param img The matrix to be scaled.
 * @return A copy of the matrix, scaled to uchar.
 */
Mat copy_and_scale_uchar
(
    Mat img
);

/**
 * \brief Image normalizer to be fitted on a given image, that will further
 * scale images given the preceding fitted values. Enables to maintain a
 * coherent color scale in plots for instance.
 */
class ImageConverter_uchar
{
public:
    /**
     * \brief Given an image, fit the internal normalization parameters.
     * @param img Input image
     * @param saturate Flag which indicate whether the normalizer should
     * cut the lower 2% and higher 2% of the range of values (in order to
     * be robust to outliers)
     */
    void fit(const Mat& img, bool saturate = true);
    /**
     * \brief Scales the input image to uchar and return a copy.
     * @param src Input image
     * @param dst Output image
     */
    void copy_and_scale(const Mat& src, Mat& dst);
private:
    double min;
    double max;
    bool m_initialized = false;
};

/**
 * \brief Scale values to plottable (uchar 0..255) and overlays the given 
 * matrix with red lines.
 * 
 * @param img The matrix to be plotted.
 * @param fill_row Fill a row in blank? -1 for no, else row index
 * @param max_height Truncate the image to this height.
 * @param fill_col Fill a col in blank? -1 for no, else col index
 * @param max_width Truncate the image to this width.
 */
Mat draw_red_lines
(
    Mat img,
    int fill_row_red = -1,
    int max_height = -1,
    int fill_col_red = -1,
    int max_width = -1
);
    
}


#endif
