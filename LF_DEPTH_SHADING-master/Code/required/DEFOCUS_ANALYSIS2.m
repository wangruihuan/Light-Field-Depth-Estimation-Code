function defocus_response = DEFOCUS_ANALYSIS2(IM_Refoc_alpha,IM_Pinhole,LF_parameters)
%DEFOCUS_ANALYSIS 
%           Takes the image refocused to alpha, and outputs response
%           for each pixel. This is a rudimentry contrast-based peak detection.
%           Input : IM_Refoc_alpha
%           Output: defocus_response

%           EQUATION (2) (3) in paper


defocus_radius = LF_parameters.defocus_radius                             ;

%Spatial patch absolute difference detection
diff_map         = abs(IM_Refoc_alpha - IM_Pinhole);

h                = fspecial('average',[defocus_radius defocus_radius])   ;
diff_avg_map     = imfilter(diff_map,h,'symmetric')                      ;

diff_avg_map     = ((diff_avg_map(:,:,1).^2 ...
                    +diff_avg_map(:,:,2).^2 ...
                    +diff_avg_map(:,:,3).^2)/3).^(1/2)                   ;
                                
defocus_response = diff_avg_map                                          ;
end

