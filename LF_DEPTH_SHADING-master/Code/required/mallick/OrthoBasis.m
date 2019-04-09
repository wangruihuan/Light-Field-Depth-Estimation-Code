function [N,V,U] = OrthoBasis(N)
%
%     [N,V,U] = OrthoBasis(N)
%
%     Input:
%       -N:
%
%     Output:
%       -[N,V,U]: the 
%
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

if(abs(N(3))>1e-10)
    V = [0.0, N(3), -N(1)];
else
    if(abs(N(1))>1e-10)
        V = [N(3),0.0,-N(1)];
    else
        V = [N(2),0.0,0.0];
    end
end

V = V/norm(V,2);
N = N/norm(N,2);

U = cross(N,V);
U = U/norm(U,2);

end