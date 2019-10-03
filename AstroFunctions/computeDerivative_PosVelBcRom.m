function [ f ] = computeDerivative_PosVelBcRom(t,xp,AC,BC,SWinputs,r,noo,svs,F_U,M_U,maxAtmAlt,et0,jdate0)
% COMPUTEDERIVATIVE_POSVELBCROM - Compute derivatives of objects
% position, velocity and BC, and reduced-order state
%
% Syntax:  [ f ] = computeDerivative_PosVelBcRom(t,xp,AC,BC,Inp,rR,noo,svs,F_U,M_U,maxAtmAlt,et0,jdate0)
%
% Inputs:
%   t           current time: seconds since et0 [s]
%   xp          state vector: position and velocity (J2000) and BC of 
%               multiple objects and reduced order density state
%   AC          continuous-time state transition matrix
%   BC          continuous-time input matrix
%   SWinputs    Space weather inputs
%   r           number of reduced order modes [integer]
%   noo         number of objects [integer]
%   svs         state size per object [integer]
%   F_U         interpolant of gridded reduced-order modes
%   M_U         interpolant of gridded mean density
%   maxAtmAlt   maximum altitude of ROM density model
%   et0         initial ephemeris time (seconds since J2000 epoch)
%   jdate0      initial Julian date
%
% Outputs:
%    f          Time derivative of xp: dxp/dt

% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 03-Oct-2019

%------------- BEGIN CODE --------------

% Convert state to single column
x = reshape(xp,svs*noo+r,[]);

% Date and time
et = et0 + t; % Ephemeris time
jdate = jdate0 + t / 86400; % Julian date

% Space weather inputs for current time
SWinputs = interp1(SWinputs(1,:),SWinputs(2:end,:)',jdate)';

% m = state vector lenght = noo*svs + r
% n = number of states (i.e. number of sigma points in UKF)
[m,n]=size(x);

% State derivative f=dx/dt
f=zeros(m,n);

% Compute accelerations per object
for i = 1:noo
    
    % Position and velocity in ECEF frame
    xform = cspice_sxform('J2000', 'ITRF93', et ); % J2000 to ECEF transformation matrix
    x_ecef = xform*x(svs*(i-1)+1:svs*(i-1)+6,:); % State in ECEF
    rr_ecef = x_ecef(1:3,:); % Position in ECEF
    vv_ecef = x_ecef(4:6,:); % Velocity in ECEF
    mag_v_ecef = sqrt( sum( vv_ecef.^2, 1 )); % Magnitude of velocity in ECEF
    
    % Gravitational accelerations in ECEF [m/s^2]
    [accGrav_ecefx, accGrav_ecefy, accGrav_ecefz] = rungravitysphericalharmonic(rr_ecef'*1000);
    
    % Atmospheric densities [kg/m^3]
    % Position in J2000
    pos = x(svs*(i-1)+1:svs*(i-1)+3,:); 
    % Reduced order density state
    romState = x(end-r+1:end,:); 
    % Compute density using reduced-order density model [kg/m^3]
    rho = getDensityROM(pos,jdate,romState,r,F_U,M_U,maxAtmAlt);
    
    % Ballistic coefficients (BC) [m^2/(1000kg)]
    b_star = x(svs*i,:);
    
    % Accelerations in ECEF: Earth gravity + drag acceleration [km/s^2]
    aa_ecef(1,:) = accGrav_ecefx'/1000 - 1/2.*b_star.*rho.*mag_v_ecef.*vv_ecef(1,:);
    aa_ecef(2,:) = accGrav_ecefy'/1000 - 1/2.*b_star.*rho.*mag_v_ecef.*vv_ecef(2,:);
    aa_ecef(3,:) = accGrav_ecefz'/1000 - 1/2.*b_star.*rho.*mag_v_ecef.*vv_ecef(3,:);
    
    % Velocities in J2000 frame [km/s]
    f(svs*(i-1)+1,:) = x(svs*(i-1)+4,:);
    f(svs*(i-1)+2,:) = x(svs*(i-1)+5,:);
    f(svs*(i-1)+3,:) = x(svs*(i-1)+6,:);
    % Accelerations in J2000 frame [km/s^2]
    f(svs*(i-1)+4:svs*(i-1)+6,:) = xform(1:3,1:3)' * aa_ecef;
    % BC time derivative is zero
    f(svs*(i-1)+7,:) = 0;
  
end

% Time derivative of reduced-order density state
f(end-r+1:end,:) = AC * x(end-r+1:end,:) + BC * SWinputs;

% Convert state derivative to single column
f = reshape(f,[],1);

end

%------------- END OF CODE --------------
