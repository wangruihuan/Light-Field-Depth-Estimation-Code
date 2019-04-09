function [nonlocal_reflectance_cost] = NONLOCALREFLECTANCE(shading, texture_clusters, pinhole_vec, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
w_nonlocal_reflectance = LF_parameters.w_nonlocal_reflectance;
num_pixels = y_size * x_size;

nonlocal_reflectance_cost = zeros(num_pixels,1);

for i = 1:num_pixels
    similar_texture = find(texture_clusters==texture_clusters(i));
    similar_texture(similar_texture==i) = [];
    for j = 1:length(similar_texture)
        nonlocal_reflectance_cost(i) = nonlocal_reflectance_cost(i) + ((pinhole_vec(i)-shading(i))-(pinhole_vec(similar_texture(j))-shading(similar_texture(j))));
    end
end

nonlocal_reflectance_cost = w_nonlocal_reflectance * nonlocal_reflectance_cost;

end