function [imgOut,imgFree] = YoonHLRemoval(img)
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

%Calculate boundary pixels
maskBoundary = TanColorBoundary(img, TanSpecularFree(img));

%Specular Free image
imgFree = YoonSpecularFree(img);
%iterative framework
[r,c,col]=size(img);
n= 20;
x = 1;
y = 0;
imgWork = img;
for i=1:n
    x = mod((x+i),2);
    y = mod((y+i),2);
    for k=2:(r-1)
        for j=2:(c-1)
            %neighboor
            if(maskBoundary(k,j)==0)%excluding boundary pixels
                rd = sum(imgFree(k,j,:))/sum(imgFree(k+y,j+x,:));
                rds = sum(img(k,j,:))/sum(img(k+y,j+x,:));
                rd = RemoveSpecials(rd);
                rds = RemoveSpecials(rds);                
                    if(rds>rd)
                        %update x1
                        m = sum(img(k,j,:))-rd*sum(img(k+y,j+x,:));
                        imgWork(k,j,:) = img(k,j,:)-m/3;
                    end

                    if(rds<rd)
                        %update x2
                        m = sum(img(k+y,j+x,:))-RemoveSpecials(sum(img(k,j,:))/rd);
                        imgWork(k+y,j+x,:) = img(k+y,j+x,:)-m/3;
                    end
            end
        end
    end
    tmpImg = img;
    img = imgWork;
    imgWork = tmpImg;
end

imgOut = img;

end
