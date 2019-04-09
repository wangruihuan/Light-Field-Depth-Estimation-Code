function [lighting] = LIGHTINGESTIMATION(normals, shading, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
num_pixels = y_size * x_size;

normals_reshape = reshape(normals, [num_pixels 3]);

%Estimate lighting
A = zeros(num_pixels, 4);
A(:,1) = 1;
A(:,2) = normals_reshape(:,2);
A(:,3) = normals_reshape(:,3);
A(:,4) = normals_reshape(:,1);

lighting = lsqlin(A, shading(:));

end