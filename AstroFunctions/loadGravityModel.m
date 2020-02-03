function [ ] = loadGravityModel( gravmodeldegree )
%LOADGRAVITYMODEL - Load Earth gravity spherical harmonics constants and coefficients
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

% Initialize gravity model
global GM Re gravdegree C_gravmodel S_gravmodel sF_gravmod
[GM, Re, gravdegree, C_gravmodel, S_gravmodel, sF_gravmod]= initgravitysphericalharmonic(gravmodeldegree);

% Set Earth constants
global Earth_radius Earth_mass gravconst
gravconst    = 6.67259e-20; % [km^3/kg/s^2]
Earth_radius = Re/1e3; % [km]
Earth_mass   = GM*1e-9/gravconst; % [kg]

end

