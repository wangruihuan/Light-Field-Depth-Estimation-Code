function [depth_reg_image] = depth_regularization_shading(IM_Pinhole, depth, combined_confi, shading, lighting, LF_parameters)

%Compute regularized depth with shading constraint (with iterative optimization)

%Parameters
lambda_data = LF_parameters.lambda_data; 
lambda_flat  = LF_parameters.lambda_flat;
lambda_smooth = LF_parameters.lambda_smooth;
lambda_shading = 0.5;%LF_parameters.lambda_shading;
confidence_scale = 0.5;%LF_parameters.confidence_scale; 
smooth_kernel  = LF_parameters.smooth_kernel;
flat_kernel = LF_parameters.flat_kernel;
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;

%Initial depth
[A_d, b_d] = reg_data_term_combined(depth, combined_confi.^confidence_scale);
[A_s, b_s] = reg_smooth_term(IM_Pinhole, smooth_kernel);
[A_f1, b_f1] = reg_smooth_term(IM_Pinhole, flat_kernel);
[A_f2, b_f2] = reg_smooth_term(IM_Pinhole, flat_kernel');
%shading_confi = (1-im2bw(combined_confi,0.1)).*combined_confi;
shading_confi = (1-combined_confi);
%shading_confi = (shading_confi - min(shading_confi(:))) / (max(shading_confi(:)) - min(shading_confi(:)));
[A_sh, b_sh] = reg_shading_harmonics_gradient_term(shading, lighting, shading_confi.^confidence_scale, LF_parameters);
%[A_sh, b_sh] = reg_shading_gradient_term(shading, lighting, shading_confi.^confidence_scale, LF_parameters);
%[A_sh, b_sh] = reg_shading_term(shading, lighting, shading_confi.^confidence_scale, LF_parameters);
A = [lambda_data*A_d;lambda_smooth*A_s;lambda_flat*A_f1;lambda_flat*A_f2;lambda_shading*A_sh];
b = [lambda_data*b_d;lambda_smooth*b_s;lambda_flat*b_f1;lambda_flat*b_f2;lambda_shading*b_sh];

%depth_reg_image = linear_least_squares_iteration(A, b, LF_parameters);
depth_reg_image = lsqlin(A, b);

depth_reg_image = (depth_reg_image - min(depth_reg_image(:))) / (max(depth_reg_image(:)) - min(depth_reg_image(:)));

depth_reg_image = reshape(depth_reg_image, [y_size x_size]);

end