function [depth_reg_image] = depth_regularization(IM_Pinhole, combined_depth, combined_confi, LF_parameters)

%Compute regularized depth (with iterative optimization)

%Parameters
lambda_data = LF_parameters.lambda_data;
lambda_flat  = LF_parameters.lambda_flat;
lambda_smooth = LF_parameters.lambda_smooth;
confidence_scale = LF_parameters.confidence_scale; 
smooth_kernel  = LF_parameters.smooth_kernel;
flat_kernel = LF_parameters.flat_kernel;

%Initial depth
[A_d, b_d] = reg_data_term_combined(combined_depth, combined_confi.^confidence_scale);
[A_s, b_s] = reg_smooth_term(IM_Pinhole, smooth_kernel);
[A_f1, b_f1] = reg_smooth_term(IM_Pinhole, flat_kernel);
[A_f2, b_f2] = reg_smooth_term(IM_Pinhole, flat_kernel');
A = [lambda_data*A_d;lambda_smooth*A_s;lambda_flat*A_f1;lambda_flat*A_f2];
b = [lambda_data*b_d;lambda_smooth*b_s;lambda_flat*b_f1;lambda_flat*b_f2];

depth_reg_image = linear_least_squares_iteration(A, b, LF_parameters);
depth_reg_image = (depth_reg_image - min(depth_reg_image(:))) / (max(depth_reg_image(:)) - min(depth_reg_image(:)));

end