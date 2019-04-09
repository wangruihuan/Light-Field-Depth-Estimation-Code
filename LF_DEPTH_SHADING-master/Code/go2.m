% Input Data (MODIFY FOR YOUR INPUT)
input_file       = 'IMG_0254-1.png';
input_folder     = 'input/';

% housekeeping
addpath(genpath('required'));
[pathstr,name,ext] = fileparts([input_folder '/' input_file]); 
if (exist(['output/' name], 'dir') == 0)
    mkdir(['output/' name]);
end

%% CAMERA
% 1 : LYTRO1
% 2 : ILLUM
% 3 : SYNTHETIC (OURS)
camera   = 2; 

%% REGULARIZER PARAMETERS

reg_generate_all        = 1;
reg_lf_data             = 1;
reg_smoothing           = 1;
reg_shading             = 1;

reg_iteration_thres     = 0.01;


%% INTERNAL PARAMETERS

%%% LF sizes                        --------------
UV_radius           = 3;
UV_diameter         = (2*UV_radius+1);
UV_size             = UV_diameter^2;

%%% Shearing                        --------------
depth_resolution        = 256;
alpha_min               = 0.2;
alpha_max               = 2;

%%% Analysis                        --------------
% defocus analysis radius
defocus_radius          = 3;
% correspondence analysis radius
corresp_radius          = UV_radius;

%%% Regularize                      --------------
lambda_data             = 1;
lambda_flat             = 2;
lambda_smooth           = 1;
lambda_shading          = 0;
confidence_scale        = 0.6; 
iter_max                = 20;
convergence_ratio       = 0.00001;
err_epsilon             = 0.001;
smooth_kernel           = [0 -1 0;-1 4 -1;0 -1 0];
flat_kernel             = [0 0 0;-1 1 0;0 0 0];

%%% Reprojection
pixel_length            = 0.001;

%%% Shading
tau_r                   = 0.0005;
normals_num_clusters    = 200;
texture_num_clusters    = 200;
k_means_iter_max        = 500;
k_means_spatial_weight  = 0.0001;
max_connections         = 4;    
w_local_shading         = 1;
w_local_reflectance     = 1;
w_nonlocal_shading      = 0.05;
w_nonlocal_reflectance  = 0.01;

if (camera == 1)
    % LOAD CAMERA DATA
    load('required/camera_data/lytro_2');
    image_cords         = image_cords_m;
    x_size              = size(image_cords_m,2);
    y_size              = size(image_cords_m,1);
    % JPEG (RAW IMAGE)
    LF_x_size           = size(image_cords,2)*UV_diameter;
    LF_y_size           = size(image_cords,1)*UV_diameter;                                                               
    % read out file
    Lytro_RAW           = (imread([input_folder '/' input_file]));
    if (size(Lytro_RAW,3) == 3)
        Lytro_RAW_Demosaic  = im2double(Lytro_RAW);
    else
        Lytro_RAW_Demosaic  = im2double(demosaic(Lytro_RAW,'bggr'));
    end
end
if (camera == 2)
    %%%%%%%%%% PUT ILLUM CODE HERE
    load('required/camera_data/illum_1');
    image_cords = round(image_cords);
    x_size              = size(image_cords,2);
    y_size              = size(image_cords,1);
    % JPEG (RAW IMAGE)
    LF_x_size           = size(image_cords,2)*UV_diameter;
    LF_y_size           = size(image_cords,1)*UV_diameter;
    
    Lytro_RAW           = (imread([input_folder '/' input_file]));
    Lytro_RAW_Demosaic  = im2double(Lytro_RAW);
end
if (camera == 3)
    %%%%%%%%%% PUT SYNTHETIC (OUR) CODE HERE
    load('required/camera_data/synth');
    image_cords         = synth;
    x_size              = size(image_cords,2);
    y_size              = size(image_cords,1);
    % JPEG (RAW IMAGE)
    LF_x_size           = size(image_cords,2)*UV_diameter;
    LF_y_size           = size(image_cords,1)*UV_diameter;                                                               
    % read out file
    Lytro_RAW           = (imread([input_folder '/' input_file]));
    if (size(Lytro_RAW,3) == 3)
        Lytro_RAW_Demosaic  = im2double(Lytro_RAW);
    else
        Lytro_RAW_Demosaic  = im2double(demosaic(Lytro_RAW,'bggr'));
    end
end

% GATHER PARAMTERS
    LF_parameters       = struct('LF_x_size',LF_x_size,...
        'LF_y_size',LF_y_size,...
        'x_size',x_size,...
        'y_size',y_size,...
        'UV_radius',UV_radius,...
        'UV_diameter',UV_diameter,...
        'UV_size',UV_size,...
        'depth_resolution',depth_resolution,...
        'alpha_min',alpha_min,...
        'alpha_max',alpha_max,...
        'defocus_radius',defocus_radius,...
        'corresp_radius',corresp_radius,...
        'lambda_data',lambda_data,...
        'lambda_flat',lambda_flat,...
        'lambda_smooth',lambda_smooth,...
        'lambda_shading',lambda_shading,...
        'confidence_scale',confidence_scale,...
        'iter_max',iter_max,...
        'convergence_ratio',convergence_ratio,...
        'err_epsilon',err_epsilon,...
        'smooth_kernel',smooth_kernel,...
        'flat_kernel',flat_kernel,...
        'pixel_length',pixel_length,...
        'tau_r',tau_r,...
        'normals_num_clusters',normals_num_clusters,...
        'texture_num_clusters',texture_num_clusters,...
        'k_means_iter_max',k_means_iter_max,...
        'k_means_spatial_weight',k_means_spatial_weight,...
        'max_connections',max_connections,...
        'w_local_shading',w_local_shading,...
        'w_local_reflectance',w_local_reflectance,...
        'w_nonlocal_shading',w_nonlocal_shading,...
        'w_nonlocal_reflectance',w_nonlocal_reflectance...
        ) ;

%% GATHER NECESSARY DATA
fprintf('I. Remapping LF JPEG to our standard                  *******\n');
tic;
% RAW to Remap
LF_Remap            = RAW2REMAP(Lytro_RAW_Demosaic,image_cords,LF_parameters);
                    
IM_Pinhole          = REMAP2PINHOLE(LF_Remap,LF_parameters);
              
LF_Remap_YCbCr      = rgb2ycbcr(LF_Remap);
LF_Remap_LOG        = log(LF_Remap+0.00000001);

%% RESPONSES
fprintf('II. Local Estimation (Defocus and Correspondence)     *******\n');

[defocus_response, corresp_response] = compute_LFdepth(LF_Remap, IM_Pinhole, LF_parameters);

%% CONFIDENCES 

defocus_confi = ALM_CONFIDENCE(defocus_response,IM_Pinhole,LF_parameters,0);
corresp_confi = ALM_CONFIDENCE(corresp_response,IM_Pinhole,LF_parameters,0);

[defocus_confi,corresp_confi] = NORMALIZE_CONFIDENCE(defocus_confi,corresp_confi);     

combined_response = COMBINE_RESPONSES(defocus_response, corresp_response, defocus_confi, corresp_confi);

combined_depth = DEPTH_ESTIMATION(combined_response,0);

combined_confi = ALM_CONFIDENCE(combined_response,IM_Pinhole,LF_parameters,0);
combined_confi = (combined_confi - min(combined_confi(:)))/(max(combined_confi(:)) - min(combined_confi(:)));

% if (reg_generate_all)
%     % depth
%     imwrite(defocus_depth/256, ['output/' name '/0a_defocus_depth.png']);
%     imwrite(corresp_depth/256, ['output/' name '/0b_corresp_depth.png']);
%     % confidence
%     imwrite(defocus_confi, ['output/' name '/0c_defocus_confi.png']);
%     imwrite(corresp_confi, ['output/' name '/0d_corresp_confi.png']);
% end

%% INITIAL DEPTH REGULARIZATION
[initial_depth] = depth_regularization(IM_Pinhole, combined_depth, combined_confi, LF_parameters);

imwrite(initial_depth, ['F:\科研\光场相机\顶会顶刊资料汇总\Depth from Shading, Defocus, and Correspondence Using Light-Field Angular Coherence\LF_DEPTH_SHADING\Code\output\Initial_Depth\' num2str(i) '.jpg']);