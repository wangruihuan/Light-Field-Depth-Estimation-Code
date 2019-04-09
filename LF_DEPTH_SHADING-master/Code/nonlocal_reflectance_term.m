function [A_nlr, b_nlr] = nonlocal_reflectance_term(texture_clusters_mat, IM_Pinhole, w_chrom_similarity, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
w_nonlocal_reflectance = LF_parameters.w_nonlocal_reflectance;
num_pixels = y_size * x_size;

pinhole_log = log(IM_Pinhole + realmin);
pinhole_log_R = pinhole_log(:,:,1);
pinhole_log_G = pinhole_log(:,:,2);
pinhole_log_B = pinhole_log(:,:,3);

A_nlr = [spdiags(w_nonlocal_reflectance .* ones(num_pixels,1), 0, num_pixels, num_pixels) * texture_clusters_mat;
         spdiags(w_nonlocal_reflectance .* ones(num_pixels,1), 0, num_pixels, num_pixels) * texture_clusters_mat;
         spdiags(w_nonlocal_reflectance .* ones(num_pixels,1), 0, num_pixels, num_pixels) * texture_clusters_mat];
b_nlr = [spdiags(w_nonlocal_reflectance .* ones(num_pixels,1), 0, num_pixels, num_pixels) * texture_clusters_mat * sparse(pinhole_log_R(:));
         spdiags(w_nonlocal_reflectance .* ones(num_pixels,1), 0, num_pixels, num_pixels) * texture_clusters_mat * sparse(pinhole_log_G(:));
         spdiags(w_nonlocal_reflectance .* ones(num_pixels,1), 0, num_pixels, num_pixels) * texture_clusters_mat * sparse(pinhole_log_B(:))];

end