function img=RemoveSpecials(img)
%
%
%       img=RemoveSpecials(img)
%
%
%       This function removes specials: Inf and NaN
%
%       Input:
%           -img: an image which can contain float special values
%
%       Output:
%           -img: the image without float special values
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
indx=find(isnan(img)|isinf(img));

img(indx)=0;

end
