function [rho] = getDensityROM(pos,jdate,romState,r,F_U,M_U,maxAtmAlt)
% getDensityROM - Compute density using reduced-order density model
%
% Syntax:  [rho] = getDensityROM(pos,jdate,romState,r,F_U,M_U,maxAtmAlt)
%
% Inputs:
%   pos         position vectors in J2000
%   jdate       Julian date
%   romState    reduced-order density state
%   r           number of reduced order modes [integer]
%   F_U         interpolant of gridded reduced-order modes
%   M_U         interpolant of gridded mean density
%   maxAtmAlt   maximum altitude of ROM density model
%
% Outputs:
%    rho        densities at positions giving by pos
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

% Number of position vectors
n = size(pos,2);

% Date and time
[yy, mm, dd, hh, mnmn, ss] = datevec(jdate-1721058.5); % Year, month, day, hour, minute, seconds in UTC
UThrs = hh + mnmn/60 + ss/3600; % Hour of day in UTC

% Convert ECI position to longitude, latitude and altitude
[lon,lat,alt]=gc2gd(pos',yy,mm,dd,hh,mnmn,ss,0,0,0);
lon(lon>180) = lon(lon>180) - 360;

% Local solar time
lst = UThrs+lon/15;
lst(lst>24) = lst(lst>24)-24;
lst(lst<0) = lst(lst<0)+24;

% Spatial modes
UhI = zeros(n,r);
for j = 1:r
    UhI(:,j) = F_U{j}(lst,lat,alt);
end
% Mean density
MI = M_U(lst,lat,alt);

% Density
rho = 10.^(sum(UhI'.*romState,1)+MI');
rho(alt>maxAtmAlt) = 0;

end

%------------- END OF CODE --------------
