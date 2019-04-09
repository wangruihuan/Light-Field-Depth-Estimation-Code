function [normal_clusters_mat] = NORMALCLUSTERS(normals, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
normals_num_clusters = 3;%LF_parameters.texture_num_clusters;
nonlocal_radius = 5;
k_means_iter_max = LF_parameters.k_means_iter_max;
max_connections = LF_parameters.max_connections;
%k_means_spatial_weight = LF_parameters.k_means_spatial_weight;
num_pixels = y_size * x_size;

normals_reshape = reshape(normals, [num_pixels 3 1]);

features = [normals_reshape];% k_means_spatial_weight*[points_x(:) points_y(:) depth(:)]];

normal_clusters = kmeans(features, normals_num_clusters, 'MaxIter', k_means_iter_max);

ii = zeros(num_pixels*(max_connections+1),1);
jj = zeros(num_pixels*(max_connections+1),1);
ss = zeros(num_pixels*(max_connections+1),1);

count = 1;

for i = 1:num_pixels
    similar_normal = find(normal_clusters==normal_clusters(i));
    similar_normal(similar_normal==i) = [];
    neighbors = zeros((nonlocal_radius*2)+1);
    for r = -nonlocal_radius:nonlocal_radius
        neighbors(:,r+nonlocal_radius+1) = (i-nonlocal_radius:i+nonlocal_radius) + (r * y_size);
    end
    similar_normal(~ismember(similar_normal, neighbors)) = [];
    if (isempty(similar_normal))
        similar_normal = i;
    end
    num_connections = min(max_connections, length(similar_normal));
    normal_indices = randi(length(similar_normal),1,num_connections);
    ii(count:count+num_connections) = i;
    jj(count:count+num_connections-1) = similar_normal(normal_indices);
    ss(count:count+num_connections-1) = -1;
    jj(count+num_connections) = i;
    ss(count+num_connections) = num_connections;
    count = count + num_connections + 1;
end

remove_indices = (ii<1);
ii(remove_indices) = [];
jj(remove_indices) = [];
ss(remove_indices) = [];

normal_clusters_mat = sparse(ii, jj, ss, num_pixels, num_pixels);

end