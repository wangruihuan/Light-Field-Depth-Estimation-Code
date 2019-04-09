function [normals_gradient] = DEPTH2NORMALSGRADIENT(depth, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;

[points_x, points_y] = meshgrid(1:x_size,1:y_size);      
[Nx, Ny, Nz] = surfnorm(points_x, points_y, reshape(depth,[y_size x_size]));
normals = cat(3,Nx,Ny,Nz);

normals_gradient = sqrt(sum(normals,3));

normals_gradient = normals_gradient(:);

end