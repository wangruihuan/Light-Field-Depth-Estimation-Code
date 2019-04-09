function [kernel_mat_same] = CONVKERNEL2MATRIX(kernel, y_size, x_size)

kernel_mirror = rot90(kernel,2);

num_pixels = y_size * x_size;

%Convolution matrix to multiply with depth vector
kernel_mat = convmtx2(kernel, y_size, x_size);
[rr, cc] = meshgrid(1:y_size, 1:x_size);
rows = reshape((rr+(cc)*(y_size+2)+1)', 1, y_size*x_size);
kernel_mat_same = kernel_mat(rows,:);

%Make convolution matrix symmetrically-padded
kernel_mat_same(1,1) = kernel_mat_same(1,1) + kernel_mirror(1,1) + kernel_mirror(2,1) + kernel_mirror(1,2); %top left corner
kernel_mat_same(1,2) = kernel_mat_same(1,2) + kernel_mirror(3,1);
kernel_mat_same(1,1+y_size) = kernel_mat_same(1,1+y_size) + kernel_mirror(1,3);
for i = 2:(y_size-1) %left boundary
    kernel_mat_same(i,i) = kernel_mat_same(i,i) + kernel_mirror(2,1);
    kernel_mat_same(i,i-1) = kernel_mat_same(i,i-1) + kernel_mirror(1,1);
    kernel_mat_same(i,i+1) = kernel_mat_same(i,i+1) + kernel_mirror(3,1);
end
kernel_mat_same(y_size,y_size) = kernel_mat_same(y_size,y_size) + kernel_mirror(2,1) + kernel_mirror(3,1) + kernel_mirror(3,2); %bottom left corner
kernel_mat_same(y_size,y_size-1) = kernel_mat_same(y_size,y_size-1) + kernel_mirror(1,1);
kernel_mat_same(y_size,2*y_size) = kernel_mat_same(y_size,2*y_size) + kernel_mirror(3,3);
for i = y_size+1:y_size:((x_size-2)*y_size)+1 %top boundary
    kernel_mat_same(i,i) = kernel_mat_same(i,i) + kernel_mirror(1,2);
    kernel_mat_same(i,i-y_size) = kernel_mat_same(i,i-y_size) + kernel_mirror(1,1);
    kernel_mat_same(i,i+y_size) = kernel_mat_same(i,i+y_size) + kernel_mirror(1,3);
end
kernel_mat_same(((x_size-1)*y_size)+1,((x_size-1)*y_size)+1) = kernel_mat_same(((x_size-1)*y_size)+1,((x_size-1)*y_size)+1) + kernel_mirror(1,2) + kernel_mirror(1,3) + kernel_mirror(2,3); %top right
kernel_mat_same(((x_size-1)*y_size)+1,((x_size-2)*y_size)+1) = kernel_mat_same(((x_size-1)*y_size)+1,((x_size-2)*y_size)+1) + kernel_mirror(1,1);
kernel_mat_same(((x_size-1)*y_size)+1,((x_size-1)*y_size)+2) = kernel_mat_same(((x_size-1)*y_size)+1,((x_size-1)*y_size)+2) + kernel_mirror(3,3);
for i = 2*y_size:y_size:(x_size-1)*y_size %bottom boundary
    kernel_mat_same(i,i) = kernel_mat_same(i,i) + kernel_mirror(3,2);
    kernel_mat_same(i,i-y_size) = kernel_mat_same(i,i-y_size) + kernel_mirror(3,1);
    kernel_mat_same(i,i+y_size) = kernel_mat_same(i,i+y_size) + kernel_mirror(3,3);
end
for i = ((x_size-1)*y_size)+2:num_pixels-1 %right boundary
    kernel_mat_same(i,i) = kernel_mat_same(i,i) + kernel_mirror(2,3);
    kernel_mat_same(i,i-1) = kernel_mat_same(i,i-1) + kernel_mirror(1,3);
    kernel_mat_same(i,i+1) = kernel_mat_same(i,i+1) + kernel_mirror(3,3);
end
kernel_mat_same(num_pixels,num_pixels) = kernel_mat_same(num_pixels,num_pixels) + kernel_mirror(2,3) + kernel_mirror(3,2) + kernel_mirror(3,3); %bottom right corner
kernel_mat_same(num_pixels,num_pixels-y_size) = kernel_mat_same(num_pixels,num_pixels-y_size) + kernel_mirror(3,1);
kernel_mat_same(num_pixels,num_pixels-1) = kernel_mat_same(num_pixels,num_pixels-1) + kernel_mirror(1,3);

end