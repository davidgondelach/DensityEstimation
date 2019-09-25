function [alt] = altitude(r)
%ALTITUDE Compute altitude w.r.t. oblate Earth
%   r = Cartesian position vector in km

% Flattening of Earth and Earth radius constants
f=1/298.257;
req=6378.14;

% Magnitude of the position
rmag=(r(:,1).^2+r(:,2).^2+r(:,3).^2).^(0.5);

% Altitude
delta=asin(r(:,3)./rmag);
alt=rmag-req*(1-f*sin(delta).^2-(f^2/2)*(sin(2*delta).^2).*(req./rmag-0.25));

end

