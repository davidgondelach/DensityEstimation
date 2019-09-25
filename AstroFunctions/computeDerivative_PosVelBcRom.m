function [ f ] = computeDerivative_PosVelBcRom(t,xp,AC,BC,Inp,rR,noo,svs,F_U,M_U,maxAtmAlt,et0,jdate0)
% COMPUTEDERIVATIVE_POSVELBCROM - Compute derivatives of objects
% position, velocity and BC, and reduced-order state
%
% Syntax:  [ f ] = computeDerivative_PosVelBcRom(t,xp,AC,BC,Inp,rR,noo,svs,F_U,M_U,maxAtmAlt,et0,jdate0)
%
% Inputs:
%   t           current time w.r.t. et0 [s]
%   xp          current state of multiple objects and reduced order atmosphere state
%   AC          continuous state transition matrix
%   BC          continuous input matrix
%   Inp         inputs (F10.7, Kp, UT, doy)
%   rR          number of reduced order modes
%   noo         number of objects
%   svs         state size per object
%   F_U         interpolant of gridded reduced-order modes
%   M_U         interpolant of gridded mean density
%   maxAtmAlt   maximum altitude of ROM density model
%   et0         initial ephemeris time
%   jdate0      initial Julian date
%
% Outputs:
%    f          Time derivative of xp: dxp/dt

% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 24-Sep-2019

%------------- BEGIN CODE --------------

x = reshape(xp,svs*noo+rR,[]);

% Date and time
et = et0 + t;
jdate = jdate0 + t / 86400;
[yy, mm, dd, hh, mnmn, ss] = datevec(jdate-1721058.5);
UThrs = hh + mnmn/60 + ss/3600;

% Space weather
Inp = interp1(Inp(1,:),Inp(2:end,:)',jdate)';

[m,n]=size(x);
f=zeros(m,n);

% Compute densities
rho = zeros(noo,n);
for i = 1:noo
    [LON,LAT,ALT]=gc2gd(x(svs*(i-1)+1:svs*(i-1)+3,:)',yy,mm,dd,hh,mnmn,ss,0,0,0);
    LON(LON>180) = LON(LON>180) - 360;
    SLT = UThrs+LON/15;
    SLT(SLT>24) = SLT(SLT>24)-24;SLT(SLT<0) = SLT(SLT<0)+24;
    UhI = zeros(n,rR);
    for j = 1:rR
       UhI(:,j) = F_U{j}(SLT,LAT,ALT);
    end
    MI = M_U(SLT,LAT,ALT);
    rho(i,:) = 10.^(sum(UhI'.*x(end-rR+1:end,:),1)+MI');
    rho(i,ALT>maxAtmAlt) = 0;
end

b_star = zeros(noo,n);
for i = 1:noo
    
    % Compute Magnitudes
    r2=(x(svs*(i-1)+1,:).^2+x(svs*(i-1)+2,:).^2+x(svs*(i-1)+3,:).^2);
    r=r2.^(1/2);
    
    % Ballistic coefficient (BC)
    b_star(i,:) = x(svs*i,:);
    
    % Position and velocity in ECEF frame
    xform = cspice_sxform('J2000', 'ITRF93', et );
    x_ecef = xform*x(svs*(i-1)+1:svs*(i-1)+6,:);
    rr_ecef = x_ecef(1:3,:);
    vv_ecef = x_ecef(4:6,:);  
    mag_v_ecef = sqrt( sum( vv_ecef.^2, 1 ));
    
    % Gravitational accelerations
    [accGrav_ecefx, accGrav_ecefy, accGrav_ecefz] = rungravitysphericalharmonic(rr_ecef'*1000);
    
    % Accelerations in ECEF frame
    aa_ecef(1,:) = accGrav_ecefx'/1000 - 1/2.*b_star(i,:).*rho(i,:).*mag_v_ecef.*vv_ecef(1,:);
    aa_ecef(2,:) = accGrav_ecefy'/1000 - 1/2.*b_star(i,:).*rho(i,:).*mag_v_ecef.*vv_ecef(2,:);
    aa_ecef(3,:) = accGrav_ecefz'/1000 - 1/2.*b_star(i,:).*rho(i,:).*mag_v_ecef.*vv_ecef(3,:);
    
    % Velocity in J2000 frame
    f(svs*(i-1)+1,:) = x(svs*(i-1)+4,:);
    f(svs*(i-1)+2,:) = x(svs*(i-1)+5,:);
    f(svs*(i-1)+3,:) = x(svs*(i-1)+6,:);
    % Acceleration in J2000 frame
    f(svs*(i-1)+4:svs*(i-1)+6,:) = xform(1:3,1:3)' * aa_ecef;
    % BC
    f(svs*(i-1)+7,:) = 0;
  
end

% Rate of change of reduced-order modes
f(end-rR+1:end,:) = AC * x(end-rR+1:end,:) + BC * Inp;

f = reshape(f,[],1);

end

%------------- END OF CODE --------------
