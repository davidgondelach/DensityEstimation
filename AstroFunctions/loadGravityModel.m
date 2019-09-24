function [ ] = loadGravityModel( gravitymodel, gravmodeldegree )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% Initialize gravity model
global GM Re gravdegree C_gravmodel S_gravmodel sF_gravmod
[GM, Re, gravdegree, C_gravmodel, S_gravmodel, sF_gravmod]= initgravitysphericalharmonic(gravitymodel,gravmodeldegree);

% Set Earth constants
global Earth_radius Earth_mass gravconst
gravconst    = 6.67259e-20; % [km^3/kg/s^2]
Earth_radius = Re/1e3;%6378.1363; % [km]
Earth_mass   = GM*1e-9/gravconst;%5.9742e24; % [kg]

end

