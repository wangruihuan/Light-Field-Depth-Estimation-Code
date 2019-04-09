function [texture_clusters_mat] = TEXTURECLUSTERS(IM_Pinhole, LF_parameters)

%Parameters
x_size = LF_parameters.x_size;
y_size = LF_parameters.y_size;
texture_num_clusters = 3;%LF_parameters.texture_num_clusters;
nonlocal_radius = 5;
k_means_iter_max = LF_parameters.k_means_iter_max;
max_connections = LF_parameters.max_connections;
%k_means_spatial_weight = LF_parameters.k_means_spatial_weight;
num_pixels = x_size * y_size;

pinhole_chromaticity(:,:,1) = IM_Pinhole(:,:,1) ./ sum(IM_Pinhole,3);
pinhole_chromaticity(:,:,2) = IM_Pinhole(:,:,2) ./ sum(IM_Pinhole,3);
pinhole_chromaticity(:,:,3) = IM_Pinhole(:,:,3) ./ sum(IM_Pinhole,3);
pinhole_chromaticity_normalized(:,:,1) = pinhole_chromaticity(:,:,1) ./ sqrt(sum(pinhole_chromaticity.^2,3));
pinhole_chromaticity_normalized(:,:,2) = pinhole_chromaticity(:,:,2) ./ sqrt(sum(pinhole_chromaticity.^2,3));
pinhole_chromaticity_normalized(:,:,3) = pinhole_chromaticity(:,:,3) ./ sqrt(sum(pinhole_chromaticity.^2,3));

pinhole_chromaticity_reshape = reshape(pinhole_chromaticity_normalized, [num_pixels, 3]);

features = [pinhole_chromaticity_reshape];% k_means_spatial_weight*[points_x(:) points_y(:) depth(:)]];

texture_clusters = kmeans(features, texture_num_clusters, 'MaxIter', k_means_iter_max);

ii = zeros(num_pixels*(max_connections+1),1);
jj = zeros(num_pixels*(max_connections+1),1);
ss = zeros(num_pixels*(max_connections+1),1);

count = 1;

for i = 1:num_pixels
    similar_texture = find(texture_clusters==texture_clusters(i));
    similar_texture(similar_texture==i) = [];
    neighbors = zeros((nonlocal_radius*2)+1);
    for r = -nonlocal_radius:nonlocal_radius
        neighbors(:,r+nonlocal_radius+1) = (i-nonlocal_radius:i+nonlocal_radius) + (r * y_size);
    end
    similar_texture(~ismember(similar_texture, neighbors)) = [];
    if (isempty(similar_texture))
        similar_texture = i;
    end
    num_connections = min(max_connections, length(similar_texture));
    texture_indices = randi(length(similar_texture),1,num_connections);
    ii(count:count+num_connections) = i;
    jj(count:count+num_connections-1) = similar_texture(texture_indices);
    ss(count:count+num_connections-1) = -1;
    jj(count+num_connections) = i;
    ss(count+num_connections) = num_connections;
    count = count + num_connections + 1;
end

remove_indices = (ii<1);
ii(remove_indices) = [];
jj(remove_indices) = [];
ss(remove_indices) = [];

texture_clusters_mat = sparse(ii, jj, ss, num_pixels, num_pixels);

end