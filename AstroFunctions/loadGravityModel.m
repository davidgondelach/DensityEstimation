function [ ] = loadGravityModel( gravmodeldegree )

% Initialize gravity model
global GM Re gravdegree C_gravmodel S_gravmodel sF_gravmod
[GM, Re, gravdegree, C_gravmodel, S_gravmodel, sF_gravmod]= initgravitysphericalharmonic(gravmodeldegree);

% Set Earth constants
global Earth_radius Earth_mass gravconst
gravconst    = 6.67259e-20; % [km^3/kg/s^2]
Earth_radius = Re/1e3; % [km]
Earth_mass   = GM*1e-9/gravconst; % [kg]

end

