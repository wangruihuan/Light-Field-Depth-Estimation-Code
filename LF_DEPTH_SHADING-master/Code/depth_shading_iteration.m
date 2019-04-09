function [final_depth, final_shading] = depth_shading_iteration(IM_Pinhole, initial_depth, combined_confi, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;

iter_max = 10;

depth = initial_depth;
depth = medfilt2(depth, [11 11], 'symmetric');

pinhole_filt(:,:,1) = medfilt2(IM_Pinhole(:,:,1),'symmetric');
pinhole_filt(:,:,2) = medfilt2(IM_Pinhole(:,:,2),'symmetric');
pinhole_filt(:,:,3) = medfilt2(IM_Pinhole(:,:,3),'symmetric');

[points_x, points_y] = meshgrid(1:x_size,1:y_size);  
[Nx, Ny, Nz] = surfnorm(points_x, points_y, depth);
normals = cat(3, -Nx, Ny, Nz);    
normals(:,:,1) = imfilter(normals(:,:,1), fspecial('gaussian', [5 5]), 'symmetric');
normals(:,:,2) = imfilter(normals(:,:,2), fspecial('gaussian', [5 5]), 'symmetric');
normals(:,:,3) = imfilter(normals(:,:,3), fspecial('gaussian', [5 5]), 'symmetric');

normal_clusters_mat = NORMALCLUSTERS(normals, LF_parameters);
texture_clusters_mat = TEXTURECLUSTERS(pinhole_filt, LF_parameters);

shading = shading_regularization(pinhole_filt, normals, texture_clusters_mat, LF_parameters);

for iter = 1:iter_max
    fprintf('iteration: %i \n', iter);
    lighting = LIGHTINGESTIMATION(normals, shading, LF_parameters);
    depth = depth_regularization_shading(pinhole_filt, depth, combined_confi, shading, lighting, LF_parameters);
    depth = medfilt2(depth, [11 11], 'symmetric');
    [Nx, Ny, Nz] = surfnorm(points_x, points_y, depth);
    normals = cat(3, -Nx, Ny, Nz);
    normals(:,:,1) = imfilter(normals(:,:,1), fspecial('gaussian', [5 5]), 'symmetric');
    normals(:,:,2) = imfilter(normals(:,:,2), fspecial('gaussian', [5 5]), 'symmetric');
    normals(:,:,3) = imfilter(normals(:,:,3), fspecial('gaussian', [5 5]), 'symmetric');
    shading = shading_regularization(pinhole_filt, normals, texture_clusters_mat, LF_parameters);
    imwrite(depth, ['C:\Users\Pratul\Documents\Shading_Depth\0-NORMALS_REGULARIZATION\output\Iterate_Depth_Shading' '\depth' num2str(iter) '.jpg']);
    imwrite(shading, ['C:\Users\Pratul\Documents\Shading_Depth\0-NORMALS_REGULARIZATION\output\Iterate_Depth_Shading' '\shading' num2str(iter) '.jpg']);
end

final_depth = depth;
final_shading = shading;

end