function [local_reflectance_cost] = LOCALREFLECTANCE(shading, pinhole_vec, IM_Pinhole, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
tau_r = LF_parameters.tau_r;
w_local_reflectance = LF_parameters.w_local_reflectance;
num_pixels = y_size * x_size;

pinhole_chromaticity(:,:,1) = IM_Pinhole(:,:,1) ./ sqrt(sum(abs(IM_Pinhole).^2,3));
pinhole_chromaticity(:,:,2) = IM_Pinhole(:,:,2) ./ sqrt(sum(abs(IM_Pinhole).^2,3));
pinhole_chromaticity(:,:,3) = IM_Pinhole(:,:,3) ./ sqrt(sum(abs(IM_Pinhole).^2,3));
pinhole_chromaticity_reshape = reshape(pinhole_chromaticity, [num_pixels, 1, 3]);

local_reflectance_cost = zeros(num_pixels,1);
w_r = zeros(num_pixels,1);

for i = 1:num_pixels
    neighbors = [i-y_size-1, i-y_size, i-y_size+1, i-1, i+1, i+y_size-1, i+y_size, i+y_size+1];
    neighbors(neighbors==i) = [];
    neighbors(neighbors<1) = [];
    neighbors(neighbors>num_pixels) = [];
    for j = 1:length(neighbors)
        local_reflectance_cost(i) = local_reflectance_cost(i) + ((pinhole_vec(i)-shading(i))-(pinhole_vec(neighbors(j))-shading(neighbors(j))));
        if ((1 - dot(pinhole_chromaticity_reshape(i,:), pinhole_chromaticity_reshape(neighbors(j),:))) < tau_r)
            w_r(i) = w_local_reflectance;
        end
    end
end

local_reflectance_cost = w_r .* local_reflectance_cost;

end