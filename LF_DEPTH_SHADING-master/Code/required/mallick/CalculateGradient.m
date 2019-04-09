function G = CalculateGradient(img)
%   
%        G = CalculateGradient(img)
%
%        This function computes the gradient of an input image L
%
%        Input:
%           -img: a gray scale image
%
%        Output:
%           - G: the gradient of img
%
%     Copyright (C) 2011  Francesco Banterle
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

col = size(img,3);

if(col>1)
    error('The image needs to be a gray scale image!');
end

kernelX = [0,0,0;-1,0,1;0,0,0];
kernelY = [0,1,0;0,0,0;0,-1,0];
G = struct('fx',imfilter(img,kernelX,'same')/2,'fy',imfilter(img,kernelY,'same')/2);

end