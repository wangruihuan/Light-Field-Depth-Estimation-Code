function [l_vec] = POINTLIGHTESTIMATION(normals, shading, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
num_pixels = y_size * x_size;

normals_reshape = reshape(normals, [num_pixels 3]);

l_vec = lsqlin(normals_reshape, shading(:));

l_vec = l_vec/norm(l_vec);

end