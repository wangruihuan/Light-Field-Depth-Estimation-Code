function [nonlocal_shading_cost] = NONLOCALSHADING(shading, normal_clusters, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
w_nonlocal_shading = LF_parameters.w_nonlocal_shading;
num_pixels = y_size * x_size;

nonlocal_shading_cost = zeros(num_pixels,1);

for i = 1:num_pixels
    similar_normals = find(normal_clusters==normal_clusters(i));
    similar_normals(similar_normals==i) = [];
    for j = 1:length(similar_normals)
        nonlocal_shading_cost(i) = nonlocal_shading_cost(i) + (shading(i)-shading(similar_normals(j)));
    end
end

nonlocal_shading_cost = w_nonlocal_shading * nonlocal_shading_cost;

end