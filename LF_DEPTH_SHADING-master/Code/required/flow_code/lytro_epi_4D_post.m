clear all;
addpath('flow_code');
addpath(genpath('flow_code/utils'));
path_name = 'buddha';
file_name = 'buddha';

%% cleanup

folder_out_name = ['analysis_dataset' path_name '/'];
if (exist(folder_out_name) == 0)
    mkdir(folder_out_name);
end
folder_out_name_2 = ['analysis_dataset' path_name '/refocus_images/'];
if (exist(folder_out_name_2) == 0)
    mkdir(folder_out_name_2);
end

file_name = ['input/hci_dataset/' path_name '/' file_name '.h5'];
info = hdf5info(file_name);

% get window size, etc.
x_size_fin  = double(info.GroupHierarchy.Attributes(1).Value);
y_size_fin  = double(info.GroupHierarchy.Attributes(2).Value);
window_side = double(info.GroupHierarchy.Attributes(3).Value);
v_baseline  = double(info.GroupHierarchy.Attributes(5).Value);
h_baseline  = double(info.GroupHierarchy.Attributes(6).Value);
focal_len   = double(info.GroupHierarchy.Attributes(7).Value);
cam_dist    = double(info.GroupHierarchy.Attributes(7).Value);


% their estimated depth
depth_address = info.GroupHierarchy.Datasets(1).Name;
gt_address    = info.GroupHierarchy.Datasets(2).Name;
lf_address    = info.GroupHierarchy.Datasets(3).Name;

depth_image   = hdf5read(file_name,depth_address)   ;
gt_iamage     = hdf5read(file_name,gt_address)      ;
lf_image      = hdf5read(file_name,lf_address)      ;

%% REMAP
stereo_diff = floor(window_side./2)  ;
window_size = window_side*window_side   ;

lf_image      = double(lf_image)./255;

%% pinhole image
im_pinhole = zeros(y_size_fin,x_size_fin,3);
im1 = zeros(y_size_fin,x_size_fin,3);
im2 = zeros(y_size_fin,x_size_fin,3);

for y = 1:y_size_fin
    for x = 1:x_size_fin
        im_pinhole(y,x,1) = lf_image(1,x,y,stereo_diff+1,stereo_diff+1);
        im_pinhole(y,x,2) = lf_image(2,x,y,stereo_diff+1,stereo_diff+1);
        im_pinhole(y,x,3) = lf_image(3,x,y,stereo_diff+1,stereo_diff+1);
        
        im1(y,x,1) = lf_image(1,x,y,stereo_diff,stereo_diff+1);
        im1(y,x,2) = lf_image(2,x,y,stereo_diff,stereo_diff+1);
        im1(y,x,3) = lf_image(3,x,y,stereo_diff,stereo_diff+1);
        
        im2(y,x,1) = lf_image(1,x,y,stereo_diff+2,stereo_diff+1);
        im2(y,x,2) = lf_image(2,x,y,stereo_diff+2,stereo_diff+1);
        im2(y,x,3) = lf_image(3,x,y,stereo_diff+2,stereo_diff+1);
    end
end
%% FOLDER

dataset = folder_out_name;

%% READ

open([dataset '/10-pinhole.mat']);
im_pinhole = ans.im_pinhole;
open([dataset '/11-other_variables.mat']);
depth_resolution = ans.depth_resolution;
curve_exponent = ans.curve_exponent;
alpha_min = ans.alpha_min;
alpha_max = ans.alpha_max;
x_size_fin = ans.x_size_fin;
y_size_fin = ans.y_size_fin;
open([dataset '/12-remap_and_pinhole.mat']);
output_image_buffer = ans.output_image_buffer;
clear ans

initial_shear       = im2double(imread([dataset '/1-shear_depth_estimate.png']));
initial_shear_conf  = im2double(imread([dataset '/2-shear_depth_confidence.png']));

initial_corr        = im2double(imread([dataset '/3-corre_depth_estimate.png']));
initial_corr_conf   = im2double(imread([dataset '/4-corre_depth_confidence.png']));

combined_estimate   = (initial_shear + initial_corr) ./ 2;

initial = zeros(size(combined_estimate,1),size(combined_estimate,2),2);
initial(:,:,1) = combined_estimate;


uv = estimate_flow_interface(im1, im2, 'classic+nl-fast', 0,initial);
figure; subplot(1,2,1);imshow(uint8(flowToColor(uv))); title('Middlebury color coding');
subplot(1,2,2); plotflow(uv);   title('Vector plot');