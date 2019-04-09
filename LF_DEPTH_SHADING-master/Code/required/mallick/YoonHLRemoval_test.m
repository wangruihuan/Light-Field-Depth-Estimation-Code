function [imgOut,imgFree] = YoonHLRemoval_test(img, light, n)
%
%        [imgOut,imgFree] = YoonHLRemoval(img)
%
%        This function computes Yoon et al.'s method
%
%        Input:
%           -img: input image with highlights
%
%        Output:
%           -imgOut: 
%           -imgFree: the specular free image
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

check3Color(img);

% Normalize light color
normalize = light / sum(light) * 3;
normalize = normalize / min(normalize);
Normalize = repmat(reshape(normalize, 1, 1, 3), [size(img, 1), size(img, 2) 1]);
img = img ./ Normalize;

%Calculate boundary pixels
maskBoundary = TanColorBoundary(img, TanSpecularFree(img));

%Specular Free image
imgFree = YoonSpecularFree(img);
%iterative framework
[r,c,col] = size(img);
x = 1;
y = 0;
imgWork = img;
for i = 1:n
    x = mod((x+i),2);
    y = mod((y+i),2);
    k = 2:(r-1);
    j = 2:(c-1);
    
    %neighboor
    mask = repmat((maskBoundary(k,j) == 0), [1 1 3]); %excluding boundary pixels
    rd = sum(imgFree(k,j,:), 3) ./ sum(imgFree(k+y,j+x,:), 3);
    rds = sum(img(k,j,:), 3) ./ sum(img(k+y,j+x,:), 3);
    rd = RemoveSpecials(rd);
    rds = RemoveSpecials(rds);                
    
    rds_larger = (rds > rd);
    %update x1
    m = rds_larger .* (sum(img(k,j,:), 3) - rd.*sum(img(k+y,j+x,:), 3));
    imgWork(k,j,:) = img(k,j,:) - repmat(m/3, [1 1 3]);
    

    rds_smaller = (rds < rd);
    %update x2
    m = rds_smaller .* (sum(img(k+y,j+x,:), 3) - RemoveSpecials(sum(img(k,j,:), 3) ./ rd));
    imgWork(k+y,j+x,:) = img(k+y,j+x,:) - repmat(m/3, [1 1 3]);
  
    imgWork(k,j,:) = img(k,j,:).*(~mask) + imgWork(k,j,:).*mask; 
  
    img = imgWork;
end

imgOut = img;

imgOut  = imgOut  .* Normalize;
imgFree = imgFree .* Normalize;

end
