function [shading] = shading_regularization_LF(Angular_Coherence_Mat, angular_use_indices, LF_Remap, normals, normal_clusters_mat, texture_clusters_mat, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
UV_diameter = LF_parameters.UV_diameter;

num_angles = length(angular_use_indices);

pinholes = zeros(y_size,x_size,3,num_angles);

for i = 1:num_angles
    angular_index = angular_use_indices(i);
    angular_x = floor((angular_index-1)/UV_diameter)+1;
    angular_y = angular_index - (UV_diameter*(angular_x-1));
    pinholes(:,:,1,i) = medfilt2(LF_Remap(angular_y:UV_diameter:end,angular_x:UV_diameter:end,1),'symmetric');
    pinholes(:,:,2,i) = medfilt2(LF_Remap(angular_y:UV_diameter:end,angular_x:UV_diameter:end,2),'symmetric');
    pinholes(:,:,3,i) = medfilt2(LF_Remap(angular_y:UV_diameter:end,angular_x:UV_diameter:end,3),'symmetric');
end

%Shading
[A_ls1, b_ls1, w_norm_similarity] = local_shading_term(normals, LF_parameters);
[A_lr1, b_lr1, w_chrom_similarity] = local_reflectance_term(pinholes(:,:,:,1), LF_parameters);
[A_nlr1, b_nlr1] = nonlocal_reflectance_term(texture_clusters_mat, pinholes(:,:,:,1), w_chrom_similarity, LF_parameters);
[A_nls1, b_nls1] = nonlocal_shading_term(normal_clusters_mat, w_norm_similarity, LF_parameters);
A_ac = Angular_Coherence_Mat;
b_ac = sparse([], [], [], x_size*y_size, 1);
A_ls = [A_ls1 sparse([],[],[],size(A_ls1,1),(num_angles-1)*size(A_ls1,2))];
A_lr = [A_lr1 sparse([],[],[],size(A_lr1,1),(num_angles-1)*size(A_lr1,2))];
A_nlr = [A_nlr1 sparse([],[],[],size(A_nlr1,1),(num_angles-1)*size(A_nlr1,2))];
A_nls = [A_nls1 sparse([],[],[],size(A_nls1,1),(num_angles-1)*size(A_nls1,2))];
b_ls = b_ls1;
b_lr = b_lr1;
b_nlr = b_nlr1;
b_nls = b_nls1;
for i = 2:num_angles
    [A_lr_next, b_lr_next, ~] = local_reflectance_term(pinholes(:,:,:,i), LF_parameters);
    [A_nlr_next, b_nlr_next] = nonlocal_reflectance_term(texture_clusters_mat, pinholes(:,:,:,i), w_chrom_similarity, LF_parameters);
    A_lr = [A_lr;
            sparse([],[],[],size(A_lr1,1),(i-1)*size(A_lr1,2))  A_lr_next sparse([],[],[],size(A_lr1,1),(num_angles-i)*size(A_lr1,2))];
    A_nlr = [A_nlr;
            sparse([],[],[],size(A_nlr1,1),(i-1)*size(A_nlr1,2))  A_nlr_next sparse([],[],[],size(A_nlr1,1),(num_angles-i)*size(A_nlr1,2))];  
    A_ls = [A_ls;
            sparse([],[],[],size(A_ls1,1),(i-1)*size(A_ls1,2))  A_ls1 sparse([],[],[],size(A_ls1,1),(num_angles-i)*size(A_ls1,2))];
    A_nls = [A_nls;
            sparse([],[],[],size(A_nls1,1),(i-1)*size(A_nls1,2))  A_nls1 sparse([],[],[],size(A_nls1,1),(num_angles-i)*size(A_nls1,2))];
    b_ls = [b_ls;b_ls1];
    b_lr = [b_lr;b_lr_next];
    b_nlr = [b_nlr;b_nlr_next];
    b_nls = [b_nls;b_nls1];
end

A = [A_ls; A_lr; A_nlr; A_nls; A_ac];
b = [b_ls; b_lr; b_nlr; b_nls; b_ac];

shading_log = full(lsqlin(A,b));

shading_array = zeros(y_size,x_size,num_angles);
for i = 1:num_angles
    shading_temp = exp(shading_log(((i-1)*(x_size*y_size))+1:(i*x_size*y_size)));
    shading_vec = (shading_temp - min(shading_temp(:))) / (max(shading_temp(:))-min(shading_temp(:)));
    shading_array(:,:,i) = reshape(shading_vec, [y_size x_size]);
end

shading = shading_array(:,:,angular_use_indices==25);
            
end