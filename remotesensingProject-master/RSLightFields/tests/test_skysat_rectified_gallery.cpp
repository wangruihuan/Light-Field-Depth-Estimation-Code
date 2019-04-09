#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <string>
#include <rslf.hpp>


/*
 * test_skysat_rectified_gallery.cpp
 * 
 * Reads skysat files from folder and display as a gallery.
 * Press keys left/right to move to previous/next image
 * Press any other key to exit
 * Provide argument 0 or 1 to read galleries with steps 18 or 1.
 */


#define PATH_TO_SKYSAT_RECTIFIED_STEP1 "../data/skysat_lasvegas_rectified/rectified_equalized_resized_frames_step1/"
#define PATH_TO_SKYSAT_RECTIFIED_STEP18 "../data/skysat_lasvegas_rectified/rectified_equalized_resized_frames_step18/"


int main(int argc, char* argv[])
{
    
    std::cout << "remoteSensingProject gallery" << std::endl;
    
    // Choose mode:
    // 0: step18 (default)
    // 1: step1
    int chosen_mode = 0;
    if (argc > 1)
    {
        std::cout << argv[1] << std::endl;
        chosen_mode = std::stoi(argv[1]);
    }
    
    std::string path_to_images;
    switch (chosen_mode)
    {
        case 0:
            path_to_images = PATH_TO_SKYSAT_RECTIFIED_STEP18;
            break;
        case 1:
            path_to_images = PATH_TO_SKYSAT_RECTIFIED_STEP1;
            break;
        default:
            throw std::string("Error invalid chosen_mode");
            break;
    }
    
    std::cout << "Reading path " << path_to_images << std::endl;
    
    /*
     * Read image from file
     */
    cv::Mat img;


    /*
     * Display image
     */
    int current_image_counter = 0;
    bool stop_flag = false;
    bool first_flag = true;
    
    rslf::ImageConverter_uchar converter;
    do 
    {
        // Get next image
        int n_zero = 3;
        std::string old_string = std::to_string(current_image_counter);
        std::string new_string = 
                    path_to_images + 
                    std::string(n_zero - old_string.length(), '0') + 
                    old_string + 
                    ".tif";;
        img = cv::imread(new_string, CV_LOAD_IMAGE_UNCHANGED);
        int height = img.rows;
        int width = img.cols;
        
        if (!img.empty())
        {
            if (first_flag)
                converter.fit(img);
        
            std::cout << "#" << current_image_counter << " Image size " << height << "x" << width << " type " <<
                rslf::type2str(img.type()) << std::endl;
            //~ std::cout << "Depth " << img.depth() << " channels " << img.channels() << std::endl;
           
            //~ std::cout << "First elements:" << std::endl;
            //~ for (int i=0; i<std::min(3,height); i++)
            //~ {
                //~ for (int j=0; j<std::min(3,width); j++)
                //~ {
                    //~ std::cout << img.at<float>(i, j) << "\t";
                //~ }
                //~ std::cout << std::endl;
            //~ }
            
            //~ double min, max;
            //~ cv::minMaxLoc(img, &min, &max);
            //~ cv::Mat tmp;
            //~ img.convertTo(tmp,CV_8U,255.0/(max-min));
            //~ img = tmp;
            converter.copy_and_scale(img, img);

            cv::namedWindow("Image", CV_WINDOW_AUTOSIZE);
            cv::imshow("Image", img);
            int key_pressed = cv::waitKey();
            //~ std::cout << "Read key " << key_pressed << std::endl;
            switch (key_pressed)
            {
                case 81: // Left arrow
                    if (current_image_counter > 0)
                        current_image_counter--;
                    break;
                case 83: // Right arrow
                    current_image_counter++;
                    break;
                default:
                    stop_flag = true;
                    break;
            }
            
        }
        else
        {
            if (current_image_counter > 0)
                current_image_counter--;
            else
            {
                std::cout<< "Image not loaded" << std::endl;
                stop_flag = true;
            }
        }
    } while (!stop_flag);
    
    return 0;
}

/*
 * Other paths...
 */
//~ img = cv::imread("../data/skysat_lasvegas_rectified/rectified_equalized_resized_frames_step18/000.tif", CV_LOAD_IMAGE_UNCHANGED);
//~ img = cv::imread("../data/mansion_image/mansion_image_0000.jpg", CV_LOAD_IMAGE_UNCHANGED);
//~ img = cv::imread("../data/skysat_video_port_hedland/video_frames/s02_20150507T02055427Z.tif", CV_LOAD_IMAGE_COLOR);

