#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <string>
#include <rslf_io.hpp>

/*
 * test_read_tiff.cpp
 * 
 * Reads example tiff files.
 * Provide argument 0 or 1 to read different files.
 */

std::string type2str(int type);


int main(int argc, char* argv[])
{
    
    std::cout << "remoteSensingProject Hello World" << std::endl;
    
    /*
     * Read image from file
     */
    cv::Mat img;
    
    // Select image to read
    int chosen_image = 0;
    if (argc > 1)
    {
        std::cout << argv[1] << std::endl;
        chosen_image = std::stoi(argv[1]);
    }
    
    switch(chosen_image)
    {
        case 0:
            img = rslf::read_img_from_file("../data/", "CCITT_1", "TIF", CV_LOAD_IMAGE_COLOR);
            // cv::imread("../data/CCITT_1.TIF", CV_LOAD_IMAGE_COLOR);
            break;
        default:
            img = rslf::read_img_from_file("../data/", "000", "tif", CV_LOAD_IMAGE_UNCHANGED);
            // cv::imread("../data/000.tif", CV_LOAD_IMAGE_UNCHANGED);
            break;
    }
    
    /*
     * Display image
     */
    if (img.empty())
    {
        std::cout<< "Image not loaded" << std::endl;
    }
    else
    {
        int height = img.rows;
        int width = img.cols;
        
        std::cout << "Image size " << height << "x" << width << " type " <<
            type2str(img.type()) << std::endl;
        std::cout << "Depth " << img.depth() << " channels " << img.channels() << std::endl;
        
       
        std::cout << "First elements:" << std::endl;
        for (int i=0; i<std::min(3,height); i++)
        {
            for (int j=0; j<std::min(3,width); j++)
            {
                switch(chosen_image)
                {
                    case 0:
                        std::cout << img.at<cv::Vec3b>(i, j) << "\t";
                        break;
                    case 1:
                        std::cout << img.at<float>(i, j) << "\t";
                        break;
                }
            }
            std::cout << std::endl;
        }
        
        if (chosen_image == 1)
        {
            double min, max;
            cv::minMaxLoc(img, &min, &max);
            cv::Mat tmp;
            img.convertTo(tmp,CV_8U,255.0/(max-min));
            img = tmp;
        }

        cv::namedWindow("Image", CV_WINDOW_AUTOSIZE);
        cv::imshow("Image", img);
        cv::waitKey();
    }   
    
    return 0;
}


std::string type2str(int type) {
    std::string r;

    uchar depth = type & CV_MAT_DEPTH_MASK;
    uchar chans = 1 + (type >> CV_CN_SHIFT);

    switch ( depth ) {
    case CV_8U:  r = "8U"; break;
    case CV_8S:  r = "8S"; break;
    case CV_16U: r = "16U"; break;
    case CV_16S: r = "16S"; break;
    case CV_32S: r = "32S"; break;
    case CV_32F: r = "32F"; break;
    case CV_64F: r = "64F"; break;
    default:     r = "User"; break;
    }

    r += "C";
    r += (chans+'0');

    return r;
}

/*
 * Other paths...
 */
//~ img = cv::imread("../data/skysat_lasvegas_rectified/rectified_equalized_resized_frames_step18/000.tif", CV_LOAD_IMAGE_UNCHANGED);
//~ img = cv::imread("../data/mansion_image/mansion_image_0000.jpg", CV_LOAD_IMAGE_UNCHANGED);
//~ img = cv::imread("../data/skysat_video_port_hedland/video_frames/s02_20150507T02055427Z.tif", CV_LOAD_IMAGE_COLOR);

