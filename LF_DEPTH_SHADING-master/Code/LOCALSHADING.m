function [local_shading_cost] = LOCALSHADING(shading, normals, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
tau_r = LF_parameters.tau_r;
w_local_shading = LF_parameters.w_local_shading;
num_pixels = y_size * x_size;

normals_reshape = reshape(normals, [num_pixels, 1, 3]);

local_shading_cost = zeros(num_pixels,1);
w_s = zeros(num_pixels,1);

for i = 1:num_pixels
    neighbors = [i-y_size-1, i-y_size, i-y_size+1, i-1, i+1, i+y_size-1, i+y_size, i+y_size+1];
    neighbors(neighbors==i) = [];
    neighbors(neighbors<1) = [];
    neighbors(neighbors>num_pixels) = [];
    for j = 1:length(neighbors)
        local_shading_cost(i) = local_shading_cost(i) + (shading(i)-shading(neighbors(j)));
        if ((1 - dot(normals_reshape(i,:), normals_reshape(neighbors(j),:))) < tau_r)
            w_s(i) = w_local_shading;
        end
    end
end

local_shading_cost = w_s .* local_shading_cost;

end