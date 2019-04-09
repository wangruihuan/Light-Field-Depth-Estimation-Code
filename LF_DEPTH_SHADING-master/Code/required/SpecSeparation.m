function diffuse = SpecSeparation(defocus_depth, corresp_depth,...
        defocus_confi, corresp_confi, LF_parameters, cur_out_folder_dir)

fprintf('**** Specular Separation\n')  

directoryPath = strcat(cur_out_folder_dir(1:end-1), '1/remap/');
Listing = dir(fullfile(directoryPath, '*.png'));
names = {Listing.name};
names = sort(names);

UV_diameter = LF_parameters.UV_diameter;
depth_resolution = LF_parameters.depth_resolution;

diffuse = repmat(zeros(size(corresp_depth)), [1 1 3]);

depth_est = defocus_depth.*(defocus_confi>=corresp_confi) + ...
    corresp_depth.*(corresp_confi>defocus_confi);

reverseStr = '';
for depth = 1:depth_resolution
    [I J] = ind2sub(size(corresp_depth), find(depth_est == depth));
    im = im2double(imread(strcat(directoryPath, names{depth})));
    for idx = 1:size(I, 1)
        i = (I(idx)-1)*UV_diameter+1 : (I(idx)-1)*UV_diameter+UV_diameter;
        j = (J(idx)-1)*UV_diameter+1 : (J(idx)-1)*UV_diameter+UV_diameter;
        region = reshape(im(i, j, :), UV_diameter^2, 3);
        diffuse(I(idx), J(idx), :) = reshape(min(region), [1 1 3]);
    end
    msg = sprintf('Processing: %d out of %d done!\n',depth,depth_resolution);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf(reverseStr);

% imwrite(spec, [cur_out_folder_dir '/L_S_estimate.png'])
end