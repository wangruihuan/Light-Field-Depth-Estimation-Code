function [depth_shading] = depth_shading_gradient_constraint(depth, shading, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;

shading_gradient = gradient(shading);
shading_gradient = shading_gradient(:);

data_term = @(depth_output) ((depth(:)-depth_output));
shading_term = @(depth_output) ((DEPTH2NORMALSGRADIENT(depth_output, LF_parameters)-shading_gradient));
myCost = @(depth_output) data_term(depth_output) + shading_term(depth_output);

depth_shading = lsqnonlin(myCost, zeros(y_size*x_size,1));

end