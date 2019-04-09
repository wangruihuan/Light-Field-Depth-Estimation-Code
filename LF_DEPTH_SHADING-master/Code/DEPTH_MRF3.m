function depth_output  = DEPTH_MRF3(depth,confi,IM_Pinhole,LF_parameters)
%DEPTH_MRF
%           Takes defocus and correspondence and uses MRF
%           Input : defocus_depth,corresp_depth
%                   defocus_confi,corresp_confi
%           Output: depth_output (regularized)

%           EQUATION (9) (10) (11) in paper


%class = 1xN vector of initial estimates
%unary = depthlabels x N of potential terms.(robust norm diff from initial)
%pairwise = sparse NxN link cost, positive for links, confidence & edge
%labelcost = depthlabels x depthlabels , cost of adj. depth diff. robust
%expansion = 1
% tic;
confi(~isfinite(confi)) = 0;
% confi_ratio = mean(mean(defocus_confi)) / mean(mean(corresp_confi));
depth_buffer        = depth   ;
confi_buffer        = confi   ;
% depth_buffer = corresp_depth;
% confi_buffer = corresp_confi;

% confidence gaussian for UNARY
confi_gaus_sigma_mult     = 150                                           ;
% labelcost gaussian for LABELCOST
label_gaus_radius         = 255                                           ;
label_gaus_sigma          = 10                                            ;

y_size              = LF_parameters.y_size                                ;
x_size              = LF_parameters.x_size                                ;
t_size              = x_size*y_size                                       ;
depth_resolution    = LF_parameters.depth_resolution                      ;

%% Class: Initial Depth Labels
CLASS               = reshape(depth_buffer, 1, t_size)                   ;
%% Unary: Confidence
% pixel number
n                  = 1:t_size;
% confidence  
confi_val          = 1 - confi_buffer;
confi_val          = min(1, max(0.01, confi_val));

d = 1:depth_resolution;
buffer = depth_buffer(n)';
[D, BUFFER] = meshgrid(d, buffer);
confi_gaus_sigma = repmat(confi_gaus_sigma_mult*confi_val(:), [1 depth_resolution]);
UNARY = exp(-(D - BUFFER).^2./(2*confi_gaus_sigma.^2))'./(sqrt(2*pi)*confi_gaus_sigma');
UNARY = UNARY ./ repmat(max(UNARY), [depth_resolution 1]);
UNARY = 1 - UNARY;

%% PAIRWISE: Image pixel color smoothness enforcement
width = x_size;
height = y_size;

h = fspecial('sobel');

image = im2double(IM_Pinhole);
output(:,:,1) = conv2(image(:,:,1),h,'same');
output(:,:,2) = conv2(image(:,:,2),h,'same');
output(:,:,3) = conv2(image(:,:,3),h,'same');

vertical = abs(output(:,:,1)) + abs(output(:,:,2)) + abs(output(:,:,3));
output(:,:,1) = conv2(image(:,:,1),h','same');
output(:,:,2) = conv2(image(:,:,2),h','same');
output(:,:,3) = conv2(image(:,:,3),h','same');

horizontal = abs(output(:,:,1)) + abs(output(:,:,2)) + abs(output(:,:,3));
%the confidence is quite crappy now. maybe i shouldn't use it.
numentries = (width-1)*height + (height-1)*width;
numentries = numentries*2 ;% make it sym

ii = zeros(numentries,1);
jj = zeros(numentries,1);
ss = zeros(numentries,1);

con_bias = 0.0;
plane_multiply = 0;
edge_str = 0.1;
edge_padding = 0.1;

im_edge = edge(rgb2gray(IM_Pinhole), 'canny');
canny_weight = max(max(vertical+horizontal));

count = 1;
for col = 1:width
     for row = 1:height-1
%          total_edgestr = (vertical(row,col) + vertical(row+1,col) + edge_padding)*edge_str;
         total_edgestr = (vertical(row,col) + vertical(row+1,col) + ...
             canny_weight*abs(im_edge(row, col) - im_edge(row+1, col)) + edge_padding) * edge_str;
         local_con = 1.0;%confidence(row,col) + confidence(row+1,col);
         local_plane = 0.0;%planes(row,col) + planes(row+1,col);
         %the more confident it is, the more weight on the edge
         %if it is a strong edge in image, lower it
         local_weight = (local_con + con_bias + local_plane*plane_multiply)/total_edgestr;
%          local_weight = edge_mul * log((max(max(total_edgestr))+(2*edge_sigma^2+1)-total_edgestr)/(2*edge_sigma^2));
         ss(count) = local_weight;
         ii(count) = row + (col-1)*height;
         jj(count) = row + 1 + (col-1)*height;
         count = count + 1;
         ss(count) = local_weight;
         jj(count) = row + (col-1)*height;
         ii(count) = row + 1 + (col-1)*height;
         count = count + 1;

     end
end

for col = 1:width-1
     for row = 1:height

%          total_edgestr = (horizontal(row,col) + horizontal(row,col+1) + edge_padding)*edge_str;
         total_edgestr = (horizontal(row,col) + horizontal(row,col+1) + ...
            canny_weight*abs(im_edge(row, col) - im_edge(row, col+1)) + edge_padding) * edge_str;
         local_con = 1.0;%confidence(row,col) + confidence(row,col+1);
         local_plane = 0.0;%planes(row,col) + planes(row,col+1);
         %the more confident it is, the more weight on the edge
         %if it is a strong edge in image, lower it
         local_weight = (local_con + con_bias + local_plane*plane_multiply)/total_edgestr;
%          local_weight = edge_mul * log((max(max(total_edgestr))+(2*edge_sigma^2+1)-total_edgestr)/(2*edge_sigma^2));
         ss(count) = local_weight;
         ii(count) = row + (col-1)*height;
         jj(count) = row + (col)*height;
         count = count + 1;
         ss(count) = local_weight;
         jj(count) = row + (col-1)*height;
         ii(count) = row + (col)*height;
         count = count + 1;

     end
end

PAIRWISE = sparse(ii,jj,ss);

%% LABELCOST: Depth change smoothness enforcement
LABELCOST = eye(depth_resolution,depth_resolution)                        ;
h         = fspecial('gaussian',[label_gaus_radius 1],label_gaus_sigma)   ;
LABELCOST = imfilter(LABELCOST,h,'symmetric')                             ;
LABELCOST = LABELCOST ./ repmat(max(LABELCOST), [depth_resolution 1]);
LABELCOST = 1 - LABELCOST                                                 ;

[labels E_0 E_1]    =...
            GCMex(CLASS-1, single(UNARY), PAIRWISE, single(LABELCOST),1)    ;

depth_output        = reshape(labels,height,width)                        ;
% imshow([rgb2gray(IM_Pinhole) defocus_depth/255 corresp_depth/255 depth_output]);
% 
% fprintf('%d %d Completed in %.3f seconds\n', E_0, E_1, toc);
end

