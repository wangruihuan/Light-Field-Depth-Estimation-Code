% function [Ls, spec_n] = LsEstimate(LF_Remap_o, defocus_depth, corresp_depth, defocus_confi, ...
%     corresp_confi, LF_parameters, cur_out_folder_dir, noise_var)

fprintf('**** L_S Estimation\n')  

directoryPath = strcat(cur_out_folder_dir(1:end-1), '1/remap/');
Listing = dir(fullfile(directoryPath, '*.png'));
names = {Listing.name};
names = sort(names);

UV_radius = LF_parameters.UV_radius;
UV_diameter = LF_parameters.UV_diameter;
UV_size = LF_parameters.UV_size;
LF_y_size = LF_parameters.LF_y_size;
LF_x_size = LF_parameters.LF_x_size;
depth_resolution = LF_parameters.depth_resolution;
alpha_min = LF_parameters.alpha_min;
alpha_max = LF_parameters.alpha_max;
alpha_step       = (alpha_max-alpha_min)/depth_resolution                 ;

depth_est = defocus_depth.*(defocus_confi>=corresp_confi) + ...
    corresp_depth.*(corresp_confi>defocus_confi);
% depth_est = DEPTH_MRF2(defocus_depth,corresp_depth,...
%     defocus_confi,corresp_confi,IM_Pinhole_o,LF_parameters);

spec = repmat(zeros(size(depth_est)), [1 1 3]);
dist = zeros(size(depth_est));
bright = zeros(size(depth_est));
sigma_bright = 0.1;
sigma_dist = 30;

LF_Remap_no_spec = LF_Remap_o;
subtract = zeros(size(LF_Remap_o));
sub_num = zeros(size(LF_Remap_o, 1), size(LF_Remap_o, 2));

reverseStr = '';
for depth = 1:depth_resolution
    [I, J] = ind2sub(size(depth_est), find(depth_est == depth));
    im = im2double(imread(strcat(directoryPath, names{depth})));
    for idx = 1:size(I, 1)
        i = (I(idx)-1)*UV_diameter+1 : (I(idx)-1)*UV_diameter+UV_diameter;
        j = (J(idx)-1)*UV_diameter+1 : (J(idx)-1)*UV_diameter+UV_diameter;
        region = reshape(im(i, j, :), UV_size, 3);
        if (mean(var(region)) > 0.001)
            opts = statset('MaxIter', 10);
            [cen_idx, centroid, sumD] = kmeans(region, 2, 'start', ...
                [mean(region)-0.01; mean(region)+0.01], 'options', opts);
            spec(I(idx), J(idx), :) = reshape(centroid(2, :) - centroid(1, :), [1 1 3]);
            dist(I(idx), J(idx)) = exp(sumD(2)^2/(2*sigma_dist^2));
            bright(I(idx), J(idx)) = exp(-(sqrt(3)-norm(max(centroid))).^2/(2*sigma_bright^2));
            spec_n = squeeze((abs(spec(I(idx), J(idx), :)).*repmat(bright(I(idx), J(idx))...
                ./dist(I(idx), J(idx)), [1 1 3])));
            if (isfinite(norm(spec_n)) && norm(spec_n) > 0.01)
                region_sub = zeros(size(region));
%                 region(cen_idx == 2, :) = region(cen_idx == 2, :) - repmat(spec_n', [sum(cen_idx == 2) 1]);
%                 region = reshape(region, UV_diameter, UV_diameter, 3);
                region_sub(cen_idx == 2, :) = repmat(squeeze(abs(spec(I(idx), J(idx), :)))', [sum(cen_idx == 2) 1]);
                region_sub = reshape(region_sub, UV_diameter, UV_diameter, 3);
                v = -UV_radius : UV_radius;
                u = -UV_radius : UV_radius;
                alpha = alpha_min + (depth-1)*alpha_step;
                y_idx = (I(idx) - v*(1-1/alpha) - 1) * UV_diameter + v + UV_radius+1;
                x_idx = (J(idx) - u*(1-1/alpha) - 1) * UV_diameter + u + UV_radius+1;
%                 LF_Remap_no_spec(floor(y_idx), floor(x_idx), :) = region;
%                 LF_Remap_no_spec(ceil(y_idx), ceil(x_idx), :) = region;
%                 x_idx = max(1, min(LF_x_size, x_idx));
%                 y_idx = max(1, min(LF_y_size, y_idx));
%                 x1 = floor(x_idx); x2 = ceil(x_idx);
%                 y1 = floor(y_idx); y2 = ceil(y_idx);
%                 subtract(y1, x1, :) = subtract(y1, x1, :) + region_sub;
%                 subtract(y1, x2, :) = subtract(y1, x2, :) + region_sub;
%                 subtract(y2, x1, :) = subtract(y2, x1, :) + region_sub;
%                 subtract(y2, x2, :) = subtract(y2, x2, :) + region_sub;
%                 sub_num(y1, x1, :) = sub_num(y1, x1, :) + 1;
%                 sub_num(y1, x2, :) = sub_num(y1, x2, :) + 1;
%                 sub_num(y2, x1, :) = sub_num(y2, x1, :) + 1;
%                 sub_num(y2, x2, :) = sub_num(y2, x2, :) + 1;
                x_idx = max(1, min(LF_x_size, floor(x_idx)));
                y_idx = max(1, min(LF_y_size, floor(y_idx)));
                for m = 1:length(y_idx)
                    for n = 1:length(x_idx)
                        i = max(1, (y_idx(m)-14)) : min(LF_y_size, y_idx(m)+20);
                        j = max(1, (x_idx(n)-14)) : min(LF_x_size, x_idx(n)+20);
                        original_color  = LF_Remap_o(i, j, :);
                        sub = original_color - repmat(reshape(centroid(2, :), [1 1 3]), [size(original_color, 1) size(original_color, 2) 1]);
                        sub = (sum(sub > 0, 3) == 3) | (sqrt(sum(sub.^2, 3)) < 0.1);
%                         if (sum(original_color - centroid(2, :) > 0) == 3 || norm(original_color - centroid(2, :)) < 0.1)
%                             LF_Remap_no_spec(i, j, :) = reshape(centroid(1, :), 1, 1, 3);
                            subtract(i, j, :) = subtract(i, j, :) + repmat(sub, [1 1 3]) .* (original_color - ...
                                repmat(reshape(centroid(1, :), [1 1 3]), [size(original_color, 1) size(original_color, 2) 1]));
                            sub_num(i, j, :) = sub_num(i, j, :) + sub;
%                         end
                    end
                end
            end
        end
    end
    msg = sprintf('Processing: %d out of %d done!\n',depth,depth_resolution);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf(reverseStr);

subtract = subtract ./ repmat(sub_num, [1 1 3]);
subtract(~isfinite(subtract)) = 0;
LF_Remap_no_spec = LF_Remap_o - subtract;

LF_Refocus_no_spec_o = RemapAverage(LF_Remap_no_spec, LF_parameters);
map = sum(RemapAverage(subtract, LF_parameters), 3) > 0;
LF_Refocus_no_spec = BFilMap(LF_Refocus_no_spec_o, map, 5, 1, 1);

spec_n = (abs(spec).*repmat(bright./dist, [1 1 3]));
spec_n(~isfinite(spec_n)) = 0;
% spec_n = spec_n / max(max(max(spec_n)));
Ls = squeeze(mean(mean(spec_n .* (repmat(sqrt(sum(spec_n.^2, 3))>0.1, [1 1 3])))));
Ls = Ls' / norm(Ls);

% noise = (LF_Remap_o-im2double( (imread(file_path)) ));
% diffuse = im2double(imread('input\synthetic_new\diffuse.png'));
% RMSE = mean(mean(sqrt(sum((RemapAverage(diffuse+noise, LF_parameters)-LF_Refocus_no_spec).^2, 3))));

% imwrite(spec_n, [cur_out_folder_dir '/L_S_estimate_' num2str(noise_var) '.png'])
% end