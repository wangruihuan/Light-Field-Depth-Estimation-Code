function corresp_response = CORRESP_ANALYSIS2(LF_Remap_alpha,IM_Pinhole,LF_parameters)
%CORRESP_ANALYSIS 
%           Takes the LF remapped to alpha, and outputs response
%           for each pixel. 
%           Input : IM_Refoc_alpha
%           Output: corresp_response

%           EQUATION (4) (5) in paper

UV_diameter       = LF_parameters.UV_diameter;
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;


%Angular patch absolute difference from pinhole image detection
pinhole_angular = imresize(IM_Pinhole, UV_diameter, 'nearest');
angular_diff_map = mean(abs(LF_Remap_alpha - pinhole_angular),3);

corresp_response = zeros(y_size, x_size);

BLOCKMEAN_mex(x_size, y_size, UV_diameter, angular_diff_map, corresp_response);

end

