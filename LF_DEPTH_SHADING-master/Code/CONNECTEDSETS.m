function [connections_sets] = CONNECTEDSETS(LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
num_pixels = y_size * x_size;

connections_sets = spdiags(ones(num_pixels,8), [-y_size-1, -y_size, -y_size+1, -1, 1, y_size-1, y_size, y_size+1], num_pixels, num_pixels);

end