function [shading] = CALCULATE_SHADING(IM_Pinhole, depth, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
num_pixels = y_size * x_size;

pinhole_gray = rgb2gray(IM_Pinhole);
pinhole_log = log(pinhole_gray+0.000001);
pinhole_vec = pinhole_log(:);

[points_x, points_y] = meshgrid(1:x_size,1:y_size);      
[Nx, Ny, Nz] = surfnorm(points_x, points_y, depth);
depth_normals = cat(3, Nx, Ny, Nz);

texture_clusters = TEXTURECLUSTERS(IM_Pinhole, LF_parameters);
normal_clusters = NORMALCLUSTERS(depth_normals, LF_parameters);

f = @(shading) LOCALSHADING(shading, depth_normals, LF_parameters) + ...
               LOCALREFLECTANCE(shading, pinhole_vec, IM_Pinhole, LF_parameters) + ...
               NONLOCALREFLECTANCE(shading, texture_clusters, pinhole_vec, LF_parameters) + ...
               NONLOCALSHADING(shading, normal_clusters, LF_parameters);
           
J = spdiags(ones(num_pixels,8), [-y_size-1, -y_size, -y_size+1, -1, 1, y_size-1, y_size, y_size+1], num_pixels, num_pixels);
options = optimset('JacobPattern',J,'Algorithm','trust-region-reflective');
shading = lsqnonlin(f, zeros(num_pixels,1),[],[],options);
            
end