camera_idx = 2;
 if (camera_idx == 1)
    load('required/camera_data/image_cords')                              ;
else
    load('required/camera_data/image_cords_2')                             ;
end
% if (camera_idx == 1)
%     Lytro_RAW_Demosaic  = im2double(demosaic(Lytro_RAW,'bggr'))               ;
% else
%     Lytro_RAW_Demosaic  = im2double(Lytro_RAW)               ;
% end
% UV_radius           = 3                                                   ;
% UV_diameter         = (2*UV_radius+1)                                     ;
% UV_size             = UV_diameter^2                                       ;
% 
% x_size              = size(image_cords, 2)                                ;
% y_size              = size(image_cords, 1)                                ; 
% LF_x_size           = size(image_cords,2)*UV_diameter                     ;
% LF_y_size           = size(image_cords,1)*UV_diameter                     ;
% % GATHER PARAMTERS
% LF_parameters       = struct('LF_x_size',LF_x_size,...
%                              'LF_y_size',LF_y_size,...
%                              'x_size',x_size,...
%                              'y_size',y_size,...
%                              'UV_radius',UV_radius,...
%                              'UV_diameter',UV_diameter,...
%                              'UV_size',UV_size,...
%                              'depth_resolution',depth_resolution,...
%                              'alpha_min',alpha_min,...
%                              'alpha_max',alpha_max,...
%                              'defocus_radius',defocus_radius,...
%                              'corresp_radius',corresp_radius,...
%                              'WS_PENALTY_W1',WS_PENALTY_W1,...
%                              'WS_PENALTY_W2',WS_PENALTY_W2,...
%                              'lambda_flat',lambda_flat,...
%                              'lambda_smooth',lambda_smooth,...
%                              'ROBUSTIFY_SMOOTHNESS',ROBUSTIFY_SMOOTHNESS,...
%                              'gradient_thres',gradient_thres,...
%                              'SOFTEN_EPSILON',SOFTEN_EPSILON,...
%                              'CONVERGE_FRACTION',CONVERGE_FRACTION...
%                              )               ;

% file_dir = '../Code_Old/input/';
% for f = 1
% file_name = ['13-01-11/IMG_' num2str(f, '%04d')];
% file_path = [file_dir file_name '.jpg'];
% Lytro_RAW           = imread(file_path);
%                          
% % denoise
% color_sigma             = 0.01*255                                  ;
% radius                  = 5                                         ;
% brightness_mult         = 1.5                                       ;
% im_in_p   = (max(min(round(Lytro_RAW_Demosaic*255),255),0))                    ;
% im_in_p_d = bilateral_filter_c(im_in_p,radius,color_sigma)          ;
% im_in_p   = im_in_p_d/255                                           ;
% % sharpening
% im_in_p   = imfilter(im_in_p,fspecial('unsharp'))                   ;
% % brightness
% Lytro_RAW_Demosaic   = min(im_in_p * brightness_mult,1)                        ;
% LF_Remap_o = RAW2REMAP(Lytro_RAW_Demosaic,image_cords,LF_parameters)    ;
% LF_Refocus = RemapAverage(LF_Remap_o, 1, LF_parameters);
% cur_out_folder_dir = ['analysis/' file_name '_2'];
% imwrite(LF_Refocus, [cur_out_folder_dir '/IM_original.png']);
% end

for f = 11
    file_path = ['analysis/test_sets/IMG_' num2str(f, '%04d') '_2/'];
    file_d_name = [file_path 'LF_Remap_no_spec.png'];
    file_s_name = [file_path 'LF_Remap_spec.png'];
    im_d = im2double(imread(file_d_name));
    im_s = im2double(imread(file_s_name));
    im_d = RemapAverage(im_d, 1, LF_parameters);
    im_d = BFilMap(im_d, ones(y_size, x_size), 5, 1, 1);
    im_s = RemapAverage(im_s, 1, LF_parameters);
    im_s = BFilMap(im_s, ones(y_size, x_size), 5, 1, 1);
    imwrite(im_d, [file_path 'LF_Refocus_no_spec.png']);
    imwrite(im_s, [file_path 'LF_Refocus_spec.png']);
end