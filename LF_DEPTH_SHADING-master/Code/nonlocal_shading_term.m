function [A_nls, b_nls] = nonlocal_shading_term(normal_clusters_mat, w_norm_similarity, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
w_nonlocal_shading = LF_parameters.w_nonlocal_shading;
num_pixels = y_size * x_size;

A_nls = spdiags(w_nonlocal_shading * ones(num_pixels,1), 0, num_pixels, num_pixels) * normal_clusters_mat;
b_nls = sparse([],[],[],num_pixels,1);

end