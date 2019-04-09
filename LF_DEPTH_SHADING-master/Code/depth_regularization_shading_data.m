function [depth_reg_image] = depth_regularization_shading_data(depth, combined_confi, shading, lighting, LF_parameters)

%Compute regularized depth with shading constraint (with iterative optimization)

%Parameters
lambda_data = LF_parameters.lambda_data; 
lambda_shading = 0.5;%LF_parameters.lambda_shading;
confidence_scale = LF_parameters.confidence_scale; 

%Initial depth
depth_confi_diag = spdiags(combined_confi(:),0,length(combined_confi(:)),length(combined_confi(:)));
A_d = depth_confi_diag;
b_d = depth_confi_diag * depth(:);
shading_confi = (1-combined_confi);
[A_sh, b_sh] = reg_shading_harmonics_gradient_term(shading, lighting, shading_confi.^confidence_scale, LF_parameters);
A = [lambda_data*A_d;lambda_shading*A_sh];
b = [lambda_data*b_d;lambda_shading*b_sh];

depth_reg_image = linear_least_squares_iteration(A, b, LF_parameters);

depth_reg_image = (depth_reg_image - min(depth_reg_image(:))) / (max(depth_reg_image(:)) - min(depth_reg_image(:)));

end