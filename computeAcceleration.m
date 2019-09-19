function [ f ] = computeAcceleration(t,x,et0,settings,densityModelData)
% t         current time
% x         current state of multiple objects (incl BC)
% et0       initial ephemerides time


persistent Sun_mass Moon_mass
if isempty(Sun_mass)
    Sun_mass    = 1.9891e30; % [kg]
    Moon_mass   = 7.3477e22; % [kg]
end

% Date and time
et = et0 + t;

% Dynamical model
x = reshape(x,7,[]);
[m,n] = size(x);

% b_star = zeros(1,n);
aa_eci = zeros(3,n);
aa_ecef = zeros(3,n);

b_star = x(7,:) * 1e-6; % [km^2/kg]

% Position and velocity in ECEF ref frame
xform = cspice_sxform('J2000', 'ITRF93', et ); % J2000 to ECEF transformation matrix
x_ecef = xform * x(1:6,:);
rr_ecef = x_ecef(1:3,:);
vv_ecef = x_ecef(4:6,:);

if settings.thirdbody == 1
    
    % Sun position in J2000 ref frame
    rr_sun = cspice_spkezr('Sun',et,'J2000','NONE', 'Earth');
    rr_sun = rr_sun(1:3,1);
    
    % Moon position in J2000 ref frame
    rr_moon = cspice_spkezr('Moon',et,'J2000','NONE', 'Earth');
    rr_moon = rr_moon(1:3,1);
    
end

for i=1:n
    rr_sat = x(1:3,i);
    
    % Earth gravitational acceleration
    accGrav_ecef = acc_geopotentialFAST(rr_ecef(:,i));
    aa_ecef(:,i) = aa_ecef(:,i) + accGrav_ecef;
    
    % Drag acceleration
    if settings.drag ~= 0
        % Velocity w.r.t. Earth atmosphere
        mag_v_ecef = sqrt( sum( vv_ecef(:,i).^2, 1 ));
        
        % Atmospheric density
        if settings.drag == 1
            % Use reduced-order density model
            [rho] = getDensityROM(rr_sat',et,densityModelData.romStateTime,densityModelData.r,densityModelData.F_U,densityModelData.M_U);
        elseif settings.drag == 2
            % Use NRLMSISE-00 atmosphere model
            [rho] = getDensityNRLMSISE(rr_sat',et);
        elseif settings.drag == 3
            % Use JB2008 atmosphere model
            [rho] = getDensityJB2008(rr_sat',et);
        end
        
        % Drag acceleration
        accDrag_ecef = [- 1/2.*b_star(1,i).*rho.*mag_v_ecef.*vv_ecef(1,i);
                        - 1/2.*b_star(1,i).*rho.*mag_v_ecef.*vv_ecef(2,i);
                        - 1/2.*b_star(1,i).*rho.*mag_v_ecef.*vv_ecef(3,i) ];
        aa_ecef(:,i) = aa_ecef(:,i) + accDrag_ecef;
    end
    
    if settings.thirdbody == 1
        
        % Solar radiation pressure
        AoMSRP = b_star(1,i)/2.2; % Approximate area-to-mass ratio (assuming BC=Cd*AoM and Cd=2.2)
        C_R = 1.2;
        [ aa_SRP ] = acc_SRP( rr_sat, rr_sun, AoMSRP, C_R );
        aa_eci(:,i) = aa_eci(:,i) + aa_SRP;
        
        % Moon gravitational acceleration
        [ aa_Moon ] = acc_3rdbody( rr_sat, rr_moon, Moon_mass );
        aa_eci(:,i) = aa_eci(:,i) + aa_Moon;
        
        % Sun gravitational acceleration
        [ aa_Sun ] = acc_3rdbody( rr_sat, rr_sun, Sun_mass );
        aa_eci(:,i) = aa_eci(:,i) + aa_Sun;
    end
    
end

% State derivative
f        = zeros(m,n);
f(1:3,:) = x(4:6,:); % Velocity
f(4:6,:) = aa_eci + xform(1:3,1:3)'*aa_ecef; % Acceleration
f(7,:)   = 0; % Constant bstar

f = reshape(f,[],1);

end


