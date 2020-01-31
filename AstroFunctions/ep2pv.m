% This function converts equinoctial parameters to eci position
% and velocity vectors 

function [rr, vv] = ep2pv(EP, mu)

% input

%  mu     = gravitational constant (km**3/sec**2)
%  EP(1) = semilatus rectum of orbit (kilometers)
%  EP(2) = f equinoctial element
%  EP(3) = g equinoctial element
%  EP(4) = h equinoctial element
%  EP(5) = k equinoctial element
%  EP(6) = true longitude (radians)

% output

%  rr = eci position vector (kilometers)
%  vv = eci velocity vector (kilometers/second)

% Orbital Mechanics with MATLAB by David Eagle, 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unload equinoctial orbital elements

p = EP(1);
f = EP(2);
g = EP(3);
h = EP(4);
k = EP(5);
l = EP(6);

sqrtmup = sqrt(mu / p);

cosl = cos(l);

sinl = sin(l);

q = 1 + f * cosl + g * sinl;

r = p / q;

alphasqrd = h^2 - k^2;

ssqrd = 1 + h^2 + k^2;

% compute eci position vector

rr	= r / ssqrd * ...
            [   cosl + alphasqrd * cosl + 2 * h * k * sinl;
                sinl - alphasqrd * sinl + 2 * h * k * cosl;
                2 * (h * sinl - k * cosl);
            ];
        
% compute eci velocity vector

vv	= sqrtmup  / ssqrd * ...
            [   - sinl - alphasqrd * sinl + 2 * h * k * cosl - g ...
                	+ 2 * f * h * k - alphasqrd * g;
                cosl - alphasqrd * cosl - 2 * h * k * sinl + f ...
                    - 2 * g * h * k - alphasqrd * f;
                2 * (h * cosl + k * sinl + f * h ...
                    + g * k)
            ];
