function [A_lr, b_lr, w_chrom_similarity] = local_reflectance_term(IM_Pinhole, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
%tau_r = 0.295;
w_local_reflectance = LF_parameters.w_local_reflectance;
num_pixels = y_size * x_size;

pinhole_log = log(IM_Pinhole + 0.01);
%pinhole_ycbcr = rgb2ycbcr(IM_Pinhole);
%pinhole_lum_reshape = reshape(pinhole_ycbcr(:,:,1), [num_pixels 1]);
pinhole_log_R = pinhole_log(:,:,1);
pinhole_log_G = pinhole_log(:,:,2);
pinhole_log_B = pinhole_log(:,:,3);

pinhole_chromaticity(:,:,1) = IM_Pinhole(:,:,1);% ./ sum(IM_Pinhole,3);
pinhole_chromaticity(:,:,2) = IM_Pinhole(:,:,2);% ./ sum(IM_Pinhole,3);
pinhole_chromaticity(:,:,3) = IM_Pinhole(:,:,3);% ./ sum(IM_Pinhole,3);
%pinhole_chromaticity_normalized(:,:,1) = pinhole_chromaticity(:,:,1) ./ sqrt(sum(pinhole_chromaticity.^2,3));
%pinhole_chromaticity_normalized(:,:,2) = pinhole_chromaticity(:,:,2) ./ sqrt(sum(pinhole_chromaticity.^2,3));
%pinhole_chromaticity_normalized(:,:,3) = pinhole_chromaticity(:,:,3) ./ sqrt(sum(pinhole_chromaticity.^2,3));
pinhole_chromaticity_normalized = pinhole_chromaticity;

pinhole_chromaticity_reshape = reshape(pinhole_chromaticity_normalized, [num_pixels, 3]);

%w_chrom_similarity = zeros(num_pixels,8);
w_chrom_similarity = zeros(num_pixels,1);

for i = 1:num_pixels
    neighbors = [i-(2*y_size)-2, i-(2*y_size)-1, i-(2*y_size), i-(2*y_size)+1, ...
                 i-(2*y_size)+2, i-y_size-2, i-y_size-1, i-y_size, i-y_size+1, ...
                 i-y_size+2, i-2, i-1, i+1, i+2, i+y_size-2, i+y_size-1, i+y_size, ...
                 i+y_size+1, i+y_size+2, i+(2*y_size)-2, i+(2*y_size)-1, ...
                 i+(2*y_size), i+(2*y_size)+1, i+(2*y_size)+2];
    neighbors(neighbors<1) = [];
    neighbors(neighbors>num_pixels) = [];
    chrom_dist = zeros(length(neighbors),1);
    for j = 1:length(neighbors)
        chrom_dist(j) = norm(pinhole_chromaticity_reshape(i,:)-pinhole_chromaticity_reshape(neighbors(j),:));
    end
    %for j = 1:length(neighbors)
    %    w_chrom_similarity(i,j) = (1 - (chrom_dist(j)/(max(chrom_dist)+realmin)));
    %end
    w_chrom_similarity(i) = (1-(mean(chrom_dist)/(max(chrom_dist)+0.01)));
end

%for i = 1:8
%    w_chrom_similarity(:,i) = ((w_chrom_similarity(:,i))/(max(w_chrom_similarity(:,i))));
%end

w_chrom_similarity = ((w_chrom_similarity)/(max(w_chrom_similarity)));

% %Get rid of edge artifacts
% w_chrom_similarity(1:y_size,1) = 1;
% w_chrom_similarity(1:y_size:end,1) = 1;
% w_chrom_similarity(1:y_size,2) = 1;
% w_chrom_similarity(1:y_size,3) = 1;
% w_chrom_similarity(y_size:y_size:end,3) = 1;
% w_chrom_similarity(1:y_size:end,4) = 1;
% w_chrom_similarity(y_size:y_size:end,5) = 1;
% w_chrom_similarity(end-y_size+1:end,6) = 1;
% w_chrom_similarity(1:y_size:end,6) = 1;
% w_chrom_similarity(end-y_size+1:end,7) = 1;
% w_chrom_similarity(end-y_size+1:end,8) = 1;
% w_chrom_similarity(y_size:y_size:end,8) = 1;

w_r = w_local_reflectance * w_chrom_similarity;

% w_r_sparse = [CONVKERNEL2MATRIX([0 0 0; 0 1 0; 0 0 -1],y_size,x_size)*spdiags(w_r(:,1),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 0; 0 1 -1; 0 0 0],y_size,x_size)*spdiags(w_r(:,2),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 -1; 0 1 0; 0 0 0],y_size,x_size)*spdiags(w_r(:,3),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 0; 0 1 0; 0 -1 0],y_size,x_size)*spdiags(w_r(:,4),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 -1 0; 0 1 0; 0 0 0],y_size,x_size)*spdiags(w_r(:,5),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 0; 0 1 0; -1 0 0],y_size,x_size)*spdiags(w_r(:,6),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 0; -1 1 0; 0 0 0],y_size,x_size)*spdiags(w_r(:,7),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([-1 0 0; 0 1 0; 0 0 0],y_size,x_size)*spdiags(w_r(:,8),0,num_pixels,num_pixels)];

w_r_sparse = spdiags(w_r, 0, num_pixels, num_pixels);

A_lr = [w_r_sparse * CONVKERNEL2MATRIX([-1 -1 -1; -1 8 -1; -1 -1 -1], y_size, x_size);
        w_r_sparse * CONVKERNEL2MATRIX([-1 -1 -1; -1 8 -1; -1 -1 -1], y_size, x_size);
        w_r_sparse * CONVKERNEL2MATRIX([-1 -1 -1; -1 8 -1; -1 -1 -1], y_size, x_size)];
b_lr = [w_r_sparse * CONVKERNEL2MATRIX([-1 -1 -1; -1 8 -1; -1 -1 -1], y_size, x_size) * sparse(pinhole_log_R(:));
        w_r_sparse * CONVKERNEL2MATRIX([-1 -1 -1; -1 8 -1; -1 -1 -1], y_size, x_size) * sparse(pinhole_log_G(:));
        w_r_sparse * CONVKERNEL2MATRIX([-1 -1 -1; -1 8 -1; -1 -1 -1], y_size, x_size) * sparse(pinhole_log_B(:))];

end