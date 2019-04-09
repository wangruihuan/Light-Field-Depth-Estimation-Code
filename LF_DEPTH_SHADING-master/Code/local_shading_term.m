function [A_ls, b_ls, w_norm_similarity] = local_shading_term(normals, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
w_local_shading = LF_parameters.w_local_shading;
num_pixels = y_size * x_size;

normals_reshape = reshape(normals, [num_pixels, 1, 3]);

w_norm_similarity = zeros(num_pixels,1);

for i = 1:num_pixels
    neighbors = [i-(2*y_size)-2, i-(2*y_size)-1, i-(2*y_size), i-(2*y_size)+1, ...
                 i-(2*y_size)+2, i-y_size-2, i-y_size-1, i-y_size, i-y_size+1, ...
                 i-y_size+2, i-2, i-1, i+1, i+2, i+y_size-2, i+y_size-1, i+y_size, ...
                 i+y_size+1, i+y_size+2, i+(2*y_size)-2, i+(2*y_size)-1, ...
                 i+(2*y_size), i+(2*y_size)+1, i+(2*y_size)+2];
    normal_dist = zeros(length(neighbors),1);
    neighbors(neighbors<1) = [];
    neighbors(neighbors>num_pixels) = [];
    for j = 1:length(neighbors)
        normal_dist(j) = norm(normals_reshape(i,:)-normals_reshape(neighbors(j),:));
    end
    %for j = 1:length(neighbors)
    %    w_norm_similarity(i,j) = (1 - (norm_dist(j)/(max(norm_dist)+realmin)));
    %end
    w_norm_similarity(i) = (1-(mean(normal_dist)/(max(normal_dist)+0.01)));
end

w_norm_similarity = ((w_norm_similarity)/(max(w_norm_similarity(:))));

% for i = 1:8
%     w_norm_similarity(:,i) = ((w_norm_similarity(:,i))/(max(w_norm_similarity(:,i))));
% end

% %Get rid of edge artifacts
% w_norm_similarity(1:y_size,1) = 1;
% w_norm_similarity(1:y_size:end,1) = 1;
% w_norm_similarity(1:y_size,2) = 1;
% w_norm_similarity(1:y_size,3) = 1;
% w_norm_similarity(y_size:y_size:end,3) = 1;
% w_norm_similarity(1:y_size:end,4) = 1;
% w_norm_similarity(y_size:y_size:end,5) = 1;
% w_norm_similarity(end-y_size+1:end,6) = 1;
% w_norm_similarity(1:y_size:end,6) = 1;
% w_norm_similarity(end-y_size+1:end,7) = 1;
% w_norm_similarity(end-y_size+1:end,8) = 1;
% w_norm_similarity(y_size:y_size:end,8) = 1;

w_s = w_local_shading * w_norm_similarity;

% w_s_sparse = [CONVKERNEL2MATRIX([0 0 0; 0 1 0; 0 0 -1],y_size,x_size)*spdiags(w_s(:,1),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 0; 0 1 -1; 0 0 0],y_size,x_size)*spdiags(w_s(:,2),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 -1; 0 1 0; 0 0 0],y_size,x_size)*spdiags(w_s(:,3),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 0; 0 1 0; 0 -1 0],y_size,x_size)*spdiags(w_s(:,4),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 -1 0; 0 1 0; 0 0 0],y_size,x_size)*spdiags(w_s(:,5),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 0; 0 1 0; -1 0 0],y_size,x_size)*spdiags(w_s(:,6),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([0 0 0; -1 1 0; 0 0 0],y_size,x_size)*spdiags(w_s(:,7),0,num_pixels,num_pixels);
%               CONVKERNEL2MATRIX([-1 0 0; 0 1 0; 0 0 0],y_size,x_size)*spdiags(w_s(:,8),0,num_pixels,num_pixels)];

w_s_sparse = spdiags(w_s, 0, num_pixels, num_pixels);

A_ls = w_s_sparse * CONVKERNEL2MATRIX([-1 -1 -1; -1 8 -1; -1 -1 -1], y_size, x_size);
b_ls = sparse([],[],[],num_pixels,1);

end