function [defocus_response corresp_response] ...
    = compute_LFdepth(LF_Remap, IM_Pinhole, LF_parameters)
% compute_LFdepth
%        : computes the depth estimation of a Lytro LF Image

% Input  : file_path    (file path to the .jpeg file)
% Output : depth_output (x,y) buffer that contains 0-255
%                       0 is closest, 255 is farthest

% NOTE   : this function supports our Lytro camera. Lytro cameras
%          have manufacturing inconsistencies for the micro-lens array.

% SYSTEM REQUIREMENTS:
% PC/MAC/LINUX
% MATLAB 2009B (tested)

% CONTACT:
% Michael W. Tao (mtao@eecs.berkeley.edu)

% TERMS OF USE : 
% Any scientific work that makes use of our code should appropriately
% mention this in the text and cite our ICCV 2013 paper. For commerical
% use, please contact us.

% PAPER TO CITE:
% Michael W. Tao, Sunil Hadap, Jitendra Malik, and Ravi Ramamoorthi. ?Depth
% from Combining Defocus and Correspondence Using Light-Field Cameras?. 
% In Proceedings of International Conference on Computer Vision (ICCV), 
% 2013.

% BIBTEX TO CITE:
% @inproceedings{Tao13,
% author={Tao, {M. W.} and Hadap, S. and Malik, J. and Ramamoorthi, R.},
% title={Depth from combining defocus and correspondence using light-field cameras},
% year={2013},
% booktitle={ICCV}
% }

% Copyright (c) 2013
% All rights reserved.
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution      
%     * Proper citation to the paper above
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


% LYTRO COPYRIGHT :
% We are not affiliated, associated, authorized, endorsed by, or in any way
% officially connected with Lytro, or any of its subsidiaries or its 
% affiliates. The official Lytro web site is available at www.lytro.com. 
% All Lytro hardware, software, etc. are registered trademarks of Lytro.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE DEFINITIONS AND EXPLANATIONS
% x'                : x' LF coordinate
% y'                : y' LF coordinate
% x                 : x image spatial coordinate
% y                 : y image spatial coordinate
% u                 : u angular coordinate
% v                 : v angular coordinate

% image_cords       : (x,y) for each image coordinate (x,y), gives you the
%                      LF_RAW data coordinate (x',y') that contains the
%                      center pixel coordinate that corresponds to (x,y)

% Lytro_RAW         : (x',y') sensor values (needs demosaic)
%                 x':  x RAW LF coordinate
%                 y':  y RAW LF coordinate
%             RANGE : [0,255]

% Lytro_RAW_Demosaic: (x',y',c) color values (demosaiced)
%                 x':  x RAW LF coordinate
%                 y':  y RAW LF coordinate
%                 c':  c color channel (1 2 3 -> R G B format)
%             RANGE : [0,1]

% LF_Remap          : (x',y',c) color values
%                 x':  x REMAPPED LF coordinate
%                 y':  y REMAPPED LF coordinate
%                 c':  c color channel (1 2 3 -> R G B format)
%             RANGE : [0,1]
%                 x'= u + UV_diameter * (x - 1)
%                 y'= v + UV_diameter * (y - 1)

% IM_pinhole        : (x,y,c) color values of center pinhole view
%                  x:  x image spatial coordinate
%                  y:  y image spatial coordinate
%                  c:  c color channel (1 2 3 -> R G B format)
%             RANGE : [0,1]

% defocus_response  : (x,y,d) defocus - higher the better
%                  x:  x image spatial coordinate
%                  y:  y image spatial coordinate
%                  d:  depth value [1,depth_resolution]
%             RANGE : [0,1]

% corresp_response  : (x,y,d) corresp - lower the better
%                  x:  x image spatial coordinate
%                  y:  y image spatial coordinate
%                  d:  depth value [1,depth_resolution]
%             RANGE : [0,1]

% defocus_depth     : (x,y) defocus confidence - higher the better
%                  x:  x image spatial coordinate
%                  y:  y image spatial coordinate
%             RANGE : [1,depth_resoltuion]

% corresp_depth     : (x,y) depth value
%                  x:  x image spatial coordinate
%                  y:  y image spatial coordinate
%             RANGE : [1,depth_resolution]

% defocus_confi     : (x,y,d) defocus confidence - higher the better
%                  x:  x image spatial coordinate
%                  y:  y image spatial coordinate
%             RANGE : [0,1]

% corresp_confi     : (x,y,d) corresp confidence - lower the better
%                  x:  x image spatial coordinate
%                  y:  y image spatial coordinate
%             RANGE : [0,1]
%% I. Gather Parameters

y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
LF_y_size = LF_parameters.LF_y_size;
LF_x_size = LF_parameters.LF_x_size;
depth_resolution = LF_parameters.depth_resolution;
alpha_max = LF_parameters.alpha_max;
alpha_min = LF_parameters.alpha_min;
UV_diameter = LF_parameters.UV_diameter;
UV_radius = LF_parameters.UV_radius;

%% II. Compute Defocus and Correspondence Responses                 
fprintf('*************************************************************\n');
fprintf('i. Computing Defocus and Correspondence Responses    *******\n');
tic                 
%%%%% Intermediate Buffers
defocus_response = zeros(y_size,x_size,depth_resolution)                  ;
corresp_response = zeros(y_size,x_size,depth_resolution)                  ;
LF_Remap_alpha   = zeros(LF_y_size,LF_x_size,3)                           ;
IM_Refoc_alpha   = zeros(y_size,x_size,3)                                 ;

% Analysis
alpha_step       = (alpha_max-alpha_min)/depth_resolution                 ;
alpha_num        = 1                                                      ;
reverseStr       = ''                                                     ;
for alpha = alpha_min:alpha_step:alpha_max-alpha_step
    
   %%%% i. refocus            : eqn. (1) -> (8)
   % refocus
   REMAP2REFOCUS_mex(x_size,y_size,...
        UV_diameter,UV_radius,LF_Remap,...
        LF_Remap_alpha,IM_Refoc_alpha,alpha)                              ;
   %%%% ii. analysis
   % depth from defocus       : eqn. (2) (3)
   defocus_response(:,:,alpha_num) =...
                      DEFOCUS_ANALYSIS2(IM_Refoc_alpha,IM_Pinhole,LF_parameters)      ; 
   % depth from correspondence: eqn. (4) (5)
   corresp_response(:,:,alpha_num) =...
                      CORRESP_ANALYSIS2(LF_Remap_alpha,IM_Pinhole,LF_parameters)      ;
                  
   %figure(1); imagesc(IM_Refoc_alpha);
               
   % housekeeping
   msg = sprintf('Processing: %d/%d done!\n',alpha_num,depth_resolution)  ;
   fprintf([reverseStr, msg])                                             ;
   reverseStr = repmat(sprintf('\b'), 1, length(msg))                     ;
   alpha_num = alpha_num + 1                                              ;
end
fprintf(reverseStr)                                                       ;

fprintf('                                Completed in %.3f seconds\n',toc);

end