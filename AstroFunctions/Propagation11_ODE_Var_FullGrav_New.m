function [ f ] = Propagation11_ODE_Var_FullGrav_New(t,xp,AC,BC,Inp,rR,noo,svs,F_U,M_U,maxAtmAlt,et0,jdate0)
% t     current time
% xp    current state of multiple objects and reduced order atmosphere state
% AC    continuous state transition matrix
% BC    continuous input matrix
% Inp   inputs (F10.7, Kp, UT, doy)
% rR    number of reduced order modes
% noo   number of objects
% svs   state size per object
% F_U   interpolant of gridded reduced-order modes
% M_U   interpolant of gridded mean density
% yr    current year
% time  time grid
% fprintf('time step is %.4e\n',t)
x = reshape(xp,svs*noo+rR,[]);

et = et0 + t;
% [yy, mm, dd, hh, mnmn, ss] = datevec(datenum(yr,1,(Inp(4)+Inp(3)/24)));
% Date and time
jdate = jdate0 + t / 86400;
[yy, mm, dd, hh, mnmn, ss] = datevec(jdate-1721058.5);
% [doy] = dayofyear(yy,mm,dd);
UThrs = hh + mnmn/60 + ss/3600;
% UTsec = hh*3600 + mnmn*60 + ss;

% Space weather
Inp = interp1(Inp(1,:),Inp(2:end,:)',jdate)';

[m,n]=size(x);
f=zeros(m,n);

% % Conversion from kilometers to meters
% mu = 398600.4415;Re = 6378.1363;
% % Rotation Rate of Earth
% we = 7.2921158553*1e-5;
% % J2 Effect and Acceleration Coefficient
% J2 = 1.082626925638815*1e-3;

rho = zeros(noo,n);%SLT1 = zeros(s1,noo);
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

% %%% Compute Sun Moon %%%
% % Sun position in J2000 ref frame
% rr_sun = cspice_spkezr('Sun',et,'J2000','NONE', 'Earth');
% rr_sun = rr_sun(1:3,1);
% % Moon position in J2000 ref frame
% rr_moon = cspice_spkezr('Moon',et,'J2000','NONE', 'Earth');
% rr_moon = rr_moon(1:3,1);


b_star = zeros(noo,n);
for i = 1:noo
    
    % Compute Magnitudes
    r2=(x(svs*(i-1)+1,:).^2+x(svs*(i-1)+2,:).^2+x(svs*(i-1)+3,:).^2);
    r=r2.^(1/2);
%     r_3_2=r.^3;
    
    b_star(i,:) = x(svs*i,:);
    
%     % J2 Effect and Acceleration Coefficient
%     z_r=x(svs*(i-1)+3,:)./r;
%     z_r2=z_r.^2;
%     j2_coeff=3/2*J2*(mu./r2).*(Re^2./r2); % DG: NO MINUS SIGN HERE!
    
%     xform = cspice_sxform('J2000', 'ITRF93', et );
%     x_ecef = xform(1:3,1:3)*x(1:3,:);
    
    xform = cspice_sxform('J2000', 'ITRF93', et );
    x_ecef = xform*x(svs*(i-1)+1:svs*(i-1)+6,:);
    rr_ecef = x_ecef(1:3,:);
    vv_ecef = x_ecef(4:6,:);  
    mag_v_ecef = sqrt( sum( vv_ecef.^2, 1 ));
        
    [accGrav_ecefx, accGrav_ecefy, accGrav_ecefz] = rungravitysphericalharmonic(rr_ecef'*1000);
%     accGrav_ecef = zeros(3,n);
%     for j=1:n
%         accGrav_ecef(:,j) = acc_geopotentialFAST(rr_ecef(:,j));
%     end
    
    aa_ecef(1,:) = accGrav_ecefx'/1000 - 1/2.*b_star(i,:).*rho(i,:).*mag_v_ecef.*vv_ecef(1,:);
    aa_ecef(2,:) = accGrav_ecefy'/1000 - 1/2.*b_star(i,:).*rho(i,:).*mag_v_ecef.*vv_ecef(2,:);
    aa_ecef(3,:) = accGrav_ecefz'/1000 - 1/2.*b_star(i,:).*rho(i,:).*mag_v_ecef.*vv_ecef(3,:);
    
    
    % Density and Drag Terms
%     v_r1=x(svs*(i-1)+4,:)+we*x(svs*(i-1)+1,:);
%     v_r2=x(svs*(i-1)+5,:)-we*x(svs*(i-1)+2,:);
%     mag_v_r=(v_r1.^2+v_r2.^2+x(svs*(i-1)+6,:).^2).^(0.5); 


%     %%% SRP and lunisolar perturbations %%%
%     aa_eci = zeros(3,n);
%     for j=1:n
%         rr_sat = x(svs*(i-1)+1:svs*(i-1)+3,j);
%         % Solar radiation pressure
%         AoMSRP = b_star(i,j)/2.2 * 1e-9; % Approximate area-to-mass ratio (assuming BC=Cd*AoM and Cd=2.2)
%         C_R = 1.2;
%         [ aa_SRP ] = acc_SRP( rr_sat, rr_sun, AoMSRP, C_R );
%         aa_eci(:,j) = aa_SRP;
%         
%         % Moon gravitational acceleration
%         Moon_mass   = 7.3477e22; % [kg]
%         [ aa_Moon ] = acc_3rdbody( rr_sat, rr_moon, Moon_mass );
%         aa_eci(:,j) = aa_eci(:,j) + aa_Moon;
%         
%         % Sun gravitational acceleration
%         Sun_mass    = 1.9891e30; % [kg]
%         [ aa_Sun ] = acc_3rdbody( rr_sat, rr_sun, Sun_mass );
%         aa_eci(:,j) = aa_eci(:,j) + aa_Sun;
%     end

    
    % Function
    f(svs*(i-1)+1,:) = x(svs*(i-1)+4,:);
    f(svs*(i-1)+2,:) = x(svs*(i-1)+5,:);
    f(svs*(i-1)+3,:) = x(svs*(i-1)+6,:);
% 
%     f(svs*(i-1)+4,:)= accGrav(1,:) - 1/2.*b_star(i,:).*rho(i,:).*mag_v_r.*vv_ecef(1,:);
%     f(svs*(i-1)+5,:)= accGrav(2,:) - 1/2.*b_star(i,:).*rho(i,:).*mag_v_r.*vv_ecef(2,:);
%     f(svs*(i-1)+6,:)= accGrav(3,:) - 1/2.*b_star(i,:).*rho(i,:).*mag_v_r.*vv_ecef(3,:);
    f(svs*(i-1)+4:svs*(i-1)+6,:) = xform(1:3,1:3)' * aa_ecef;
%     f(svs*(i-1)+4:svs*(i-1)+6,:) = xform(1:3,1:3)' * aa_ecef + aa_eci;
    
    f(svs*(i-1)+7,:) = 0;
  
end

f(end-rR+1:end,:) = AC * x(end-rR+1:end,:) + BC * Inp;

f = reshape(f,[],1);

end


