function imgOut = imMAX(img1,img2)
%
%        imgOut = imMAX(img1,img2)
%
%        This function computes the maximum between two images
%
%        Input:
%           -img1: first image
%           -img2: second image
%
%        Output:
%           -imgOut: the maximum between img1 and img2
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

[r,c,col] = size(img1);

imgWork = zeros(r,c,2);

imgOut = zeros(r,c,col);

for i = 1:col
    imgWork(:,:,1) = img1(:,:,i);
    imgWork(:,:,2) = img2(:,:,i);
    imgOut(:,:,i) = max(imgWork,[],3);
end

end