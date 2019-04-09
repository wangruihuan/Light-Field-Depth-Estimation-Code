function [shading] = shading_regularization(IM_Pinhole, normals, texture_clusters_mat, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;

normal_clusters_mat = NORMALCLUSTERS(normals, LF_parameters);
           
%Shading
[A_ls, b_ls, w_norm_similarity] = local_shading_term(normals, LF_parameters);
[A_lr, b_lr, w_chrom_similarity] = local_reflectance_term(IM_Pinhole, LF_parameters);
[A_nlr, b_nlr] = nonlocal_reflectance_term(texture_clusters_mat, IM_Pinhole, w_chrom_similarity, LF_parameters);
[A_nls, b_nls] = nonlocal_shading_term(normal_clusters_mat, w_norm_similarity, LF_parameters);
A = [A_ls;A_lr;A_nlr;A_nls];
b = [b_ls;b_lr;b_nlr;b_nls];
%shading_log = full(linear_least_squares_iteration(A, b, LF_parameters));
shading_log = full(lsqlin(A,b));

shading = exp(shading_log);
shading = (shading - min(shading(:))) / (max(shading(:))-min(shading(:)));

shading = reshape(shading, [y_size x_size]);
            
end