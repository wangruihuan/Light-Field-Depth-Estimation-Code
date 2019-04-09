function [depth_to_points_matrix] = DEPTH2POINTS_MATRIX(LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
pixel_length = LF_parameters.pixel_length;
num_pixels = y_size * x_size;

%Matrix to reproject pixels into scene
[x,y] = meshgrid(1:x_size,1:y_size);
points_x = (x - (x_size / 2)) .* pixel_length;
points_y = (y - (y_size / 2)) .* pixel_length;
points_z = ones(y_size,x_size);

%3n*n matrix (where n is the number of pixels in the depth map
depth_to_points_matrix = [spdiags(points_x(:),0,num_pixels,num_pixels);...
                          spdiags(points_y(:),0,num_pixels,num_pixels);...
                          spdiags(points_z(:),0,num_pixels,num_pixels)];

end