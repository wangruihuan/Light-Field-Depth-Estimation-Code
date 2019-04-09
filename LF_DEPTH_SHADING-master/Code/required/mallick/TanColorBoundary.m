function mask = TanColorBoundary(imgNor, imgSpec, thR, thG)
%
%       mask = TanColorBoundary(imgNor, imgSpec, thR, thG)
%
%        This function 
%
%        Input:
%           -imgNor: 
%           -imgSpec:
%           -thR:
%           -thG:
%
%        Output:
%           -mask:
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

if(~exist('thR'))
    thR = 0.05;
end

if(~exist('thG'))
    thG = 0.1;
end

[rN,gN] = TanChromaticity(imgNor);
[rS,gS] = TanChromaticity(imgSpec);

deltaR = rN - rS;
deltaG = gN - gS;

indx = find(deltaR>thR&deltaG>thG);

mask = zeros(size(rN));
mask(indx) = 1;

end