function [ f ] = Propagation_FullGravDrag_ROM(t,x,et0,romStateTime,r,F_U,M_U,varargin)
% t         current time
% x         current state of multiple objects (incl BC)
% et0    initial ephemerides time
% 
% 
% % Space weather data
% persistent SWmatDaily SWmatMonthlyPred
% if isempty(SWmatDaily)
% %     SW = load('SW.txt');
%     [ SWmatDaily, SWmatMonthlyPred ] = inputSWnrlmsise( '/Users/davidgondelach/Dropbox (MIT)/Research-ROM/code/sqrt-ukf-rom-gps-2019/Data/SW-All.txt' );
% end

persistent Sun_mass Moon_mass
if isempty(Sun_mass)
    Sun_mass    = 1.9891e30; % [kg]
    Moon_mass   = 7.3477e22; % [kg]
end

% Date and time
et = et0 + t;
% jdatestr    = cspice_et2utc( et, 'J', 12 );
% jdate       = str2double(jdatestr(4:end)); % Cut trailing 'JD ' off from string
% % jdate = jdate0 + t / 86400;
% % date = datetime(jdate,'ConvertFrom','juliandate');
% % [yy, mm, dd, hh, mnmn, ss] = datevec(date);
% [yy, mm, dd, hh, mnmn, ss] = datevec(jdate-1721058.5);
% % doy = day(date,'dayofyear');
% UThrs = hh + mnmn/60 + ss/3600;
% % UTsec = hh*3600 + mnmn*60 + ss;
% % 
% % % Compute space weather
% % [ f107A, f107, ap ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jdate );
% % 
% 
% romState = interp1(romStateTime(:,end),romStateTime(:,1:r),et);
% 
% % Longitude, Latitude and Height of Satellite and Density
% [LON,LAT,ALT]=gc2gd(x(1:3)',yy,mm,dd,hh,mnmn,ss,0,0,0);
% LON(LON>180) = LON(LON>180) - 360;
% SLT = UThrs+LON/15;
% SLT(SLT>24) = SLT(SLT>24)-24; SLT(SLT<0) = SLT(SLT<0)+24;
% UhI = zeros(1,r);
% for j = 1:r
%    UhI(:,j) = F_U{j}(SLT,LAT,ALT);
% end
% MI = M_U(SLT,LAT,ALT);
% rho = 10.^(sum(UhI'.*romState,1)+MI');
[rho] = getDensityROM(x(1:3)',et,romStateTime,r,F_U,M_U);

% J2000 to ECEF
xform = cspice_sxform('J2000', 'ITRF93', et );

% Dynamical model
[m,n] = size(x);
f = zeros(m,n);

svs = 7;
b_star = zeros(1,n);
for i = 1:1
    
    b_star(i,:) = x(svs*i,:) * 1e-6; % [km^2/kg]
    
    x_ecef = xform * x(svs*(i-1)+1:svs*(i-1)+6,:);
    rr_ecef = x_ecef(1:3,:);
    vv_ecef = x_ecef(4:6,:);  
    mag_v_ecef = sqrt( sum( vv_ecef.^2, 1 ));
    
    % Drag acceleration
    dragAcc_ecef = [- 1/2.*b_star(i,:).*rho(i,:).*mag_v_ecef.*vv_ecef(1,:);
                    - 1/2.*b_star(i,:).*rho(i,:).*mag_v_ecef.*vv_ecef(2,:);
                    - 1/2.*b_star(i,:).*rho(i,:).*mag_v_ecef.*vv_ecef(3,:) ];
    
%     [accGrav_ecefx, accGrav_ecefy, accGrav_ecefz] = rungravitysphericalharmonic(rr_ecef'*1000);
    accGrav_ecef = acc_geopotentialFAST(rr_ecef);

    aa_ecef(1,:) = accGrav_ecef(1,:) + dragAcc_ecef(1);
    aa_ecef(2,:) = accGrav_ecef(2,:) + dragAcc_ecef(2);
    aa_ecef(3,:) = accGrav_ecef(3,:) + dragAcc_ecef(3);
%     
%     aa_ecef(1,:) = accGrav_ecefx'/1000 + dragAcc_ecef(1);
%     aa_ecef(2,:) = accGrav_ecefy'/1000 + dragAcc_ecef(2);
%     aa_ecef(3,:) = accGrav_ecefz'/1000 + dragAcc_ecef(3);

    aa_eci = zeros(3,n);
    
    if nargin > 7 && varargin{1}{1} == 1
        rr_sat = x(svs*(i-1)+1:svs*(i-1)+3,:);
        
        rr_sun = cspice_spkezr('Sun',et,'J2000','NONE', 'Earth');
        rr_sun = rr_sun(1:3,1);
        
        rr_moon = cspice_spkezr('Moon',et,'J2000','NONE', 'Earth');
        rr_moon = rr_moon(1:3,1);
        
        AoMSRP = b_star/2.2; % Approximate area-to-mass ratio (assuming BC=Cd*AoM and Cd=2.2)
        C_R = 1.2;
        [ a_SRP ] = acc_SRP( rr_sat, rr_sun, AoMSRP, C_R );
        aa_eci = aa_eci + a_SRP;
        
        [ a_Moon ] = acc_3rdbody( rr_sat, rr_moon, Moon_mass );
        aa_eci = aa_eci + a_Moon;

        [ a_Sun ] = acc_3rdbody( rr_sat, rr_sun, Sun_mass );
        aa_eci = aa_eci + a_Sun;
    end
    
    % State derivative
    f(svs*(i-1)+1,:)=x(svs*(i-1)+4,:);
    f(svs*(i-1)+2,:)=x(svs*(i-1)+5,:);
    f(svs*(i-1)+3,:)=x(svs*(i-1)+6,:);
    
    f(svs*(i-1)+4:svs*(i-1)+6,:) = aa_eci + xform(1:3,1:3)' * aa_ecef;
    
    f(svs*i,:) = 0; % Constant bstar
  
end

f = reshape(f,[],1);

end


