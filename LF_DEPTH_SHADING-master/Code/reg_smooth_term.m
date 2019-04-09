function [A_s, b_s] = reg_smooth_term(IM_Pinhole, kernel)

%Note: for 3x3 kernels only

%Parameters
x_size = size(IM_Pinhole,2);
y_size = size(IM_Pinhole,1);
num_pixels = x_size * y_size;

%Calculate Sum of Color Derivative in X Direction
h           = kernel/sum(sum(abs(kernel)));
im_gradient = imfilter(IM_Pinhole,h,'symmetric');
im_gradient = sqrt((im_gradient(:,:,1).^2+im_gradient(:,:,2).^2+im_gradient(:,:,3).^2)./3);
im_gradient_filt = (1-min(1,im_gradient))+0.01;
kernel_color_sum = im_gradient_filt;

%Weight Matrix of Color Differences
weights = spdiags(kernel_color_sum(:),0,length(kernel_color_sum(:)),length(kernel_color_sum(:)));

kernel_mat_same = CONVKERNEL2MATRIX(kernel, y_size, x_size);

%Outputs
A_s = weights * kernel_mat_same;
b_s = zeros(num_pixels,1);

end