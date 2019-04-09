function dim_response = DIMENSION_ANALYSIS(IM_Refoc_alpha,LF_parameters)
%CORRESP_ANALYSIS 
%           Takes the LF remapped to alpha, and outputs response
%           for each pixel. This is a rudimentry correspondence, minimum
%           variance of patches.
%           Input : IM_Refoc_alpha
%           Output: corresp_response

%           EQUATION (4) (5) in paper

corresp_radius    = LF_parameters.corresp_radius                          ;
x_size            = LF_parameters.x_size                                  ;
y_size            = LF_parameters.y_size                                  ;
UV_radius         = LF_parameters.UV_radius                               ;
UV_diameter       = LF_parameters.UV_diameter                             ;


for x = 1:x_size
    for y = 1:y_size
            RGB_Vals = IM_Refoc_alpha((y-1)*UV_diameter+1:(y-1)*UV_diameter+UV_diameter,...
                (x-1)*UV_diameter+1:(x-1)*UV_diameter+UV_diameter,:);
            
            diff_1 = RGB_Vals(:,:,1) - mean2(RGB_Vals(:,:,1));
            diff_2 = RGB_Vals(:,:,2) - mean2(RGB_Vals(:,:,2));
            diff_3 = RGB_Vals(:,:,3) - mean2(RGB_Vals(:,:,3));
            
            dim_response(y,x) = sqrt(sum(sum(diff_1.^2+diff_2.^2+diff_3.^2))/(3*UV_diameter*UV_diameter));
            
    end
end

end

