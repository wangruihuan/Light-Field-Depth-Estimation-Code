function imgOut = MallickHLRemoval(img)
%
%        imgOut = MallickHLRemoval(img)
%
%        This function computes Mallick et al.'s method using the bilateral
%        filter for speeding computations up
%
%        Input:
%           -img: input image with highlights
%
%        Output:
%           -imgOut: the specular free image
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

%col = col';

% normalization
%col = col./sum(col);
%IlluminantEstimationTan(img);
col = [1,1,1]'/3;

%From RGB to SUV
[N,V,U] = OrthoBasis(col');
M = [N;V;U];
imgM = ConvertLinearSpace(img,M);

%From SUV to (theta,phi,rho)
rho = sqrt(imgM(:,:,2).^2+imgM(:,:,3).^2);
% t2 = RemoveSpecials(imgM(:,:,2)./rho);
% t3 = RemoveSpecials(imgM(:,:,3)./rho);
theta = real(atan2(imgM(:,:,2),imgM(:,:,3))+pi);%real(acos(t2+t3));
% rho2 = sqrt(imgM(:,:,1).^2+imgM(:,:,2).^2+imgM(:,:,3).^2);
phi = real(atan2(imgM(:,:,1),rho)+pi);%RemoveSpecials(acos(imgM(:,:,1)./rho2));

%Diffusion
phiNew = max(max(phi))*bilateralFilter(phi/max(max(phi)),(theta/max(max(theta))));

%Conversion back from (theta,phi,rho) to SUV to RGB
imgMNew = zeros(size(img));
imgMNew(:,:,1) = tan(phiNew-pi).*rho;
imgMNew(:,:,2) = imgM(:,:,2);
imgMNew(:,:,3) = imgM(:,:,3);

imgOut = ConvertLinearSpace(imgMNew,inv(M));
end