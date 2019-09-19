function [ semiMajorAxisChange, ballisticCoefficient ] = computeSemiMajorAxisChangeDueToDragAveraged( rr0, vv0, jdate0, jdatef, ballisticCoefficient )
%COMPUTESEMIMAJORAXISCHANGEDUETODRAGAVERAGED - This function computes the
% change in semi-major due to drag between two epochs for an object.
% Using the initial state and object parameters, the object is propagated
% between the epochs and the change in semi-major axis due to drag is
% computed. The mean semi-major axis due to drag is computed using
% polynomial fitting and the change in mean semi-major axis is returned.
% If the object re-enters during propagation then the propagation is
% repeating with a smaller value for the ballistic coefficient.
%
% Syntax:  [ semiMajorAxisChange ] = computeSemiMajorAxisChangeDueToDrag( satrec, jdate0, jdatef, ballisticCoefficient )
%
% Inputs:
%    rr0                        - Initial J2000 position (x,y,z) [km]
%    vv0                        - Initial J2000 velocity (Vx,Vy,Vz) [km]
%    jdate0                     - Initial Julian date
%    jdatef                     - Final Julian date
%    ballisticCoefficient       - Ballistic coefficient (current estimate)
%    reflectivity               - Reflectivity coefficient
%    areaToMassRatio            - Area to mass ratio [m^2/kg]
%    propagatorSettings         - Array containing propagator settings
%
% Outputs:
%    semiMajorAxisChange        - vector of coefficients of fitted polynomial
%    ballisticCoefficient       - value of applied ballistic coefficient estimate
%
% Other m-files required: sgp4, convertTEME2ECI, runAIDA, runPlanODyn

% Author: David Gondelach
% Astronautics, University of Southampton
% July 2015; Last revision: 31-July-2015

%% Compute semi-major axis change due to drag

% Perform propagation
noo = 1; svs = 7;
x0 = [rr0; vv0; ballisticCoefficient];
tt = [0 (jdatef-jdate0)*86400];
% [tout, Xp, objectReentered] = propagateStateJ2dragNRLMSISE(x0,tt,jdate0,noo,svs);
noo = 1; svs = 8;
x0 = [rr0; vv0; 0; ballisticCoefficient];
et0  = cspice_str2et(strcat([num2str(jed2date(jdate0),'%d %d %d %d %d %.10f') 'UTC']));
    
opts = odeset('Events', @(t,x) isdecayed(t,x,noo,svs), 'RelTol',1e-10,'AbsTol',1e-10);
% [tout, xout, ~, ~, objectReentered] = ode45(@(t,x) Propagation_J2drag_NRLMSISE_withSMAdrag(t,x,jdate0,noo,svs),tt,x0,opts);
[tout, xout, ~, ~, objectReentered] = ode113(@(t,x) Propagation_FullGravDrag_NRLMSISE_withSMAdrag(t,x,jdate0,et0,noo,svs),tt,x0,opts);
Xp = reshape(xout',noo*svs,[]);

if objectReentered == 1
    % The object re-entered, so the initial BC guess was too high.
    % Re-run the function with a BC guess that is twice as small.
    ballisticCoefficient = 0.9 * ballisticCoefficient;
    
    [ semiMajorAxisChange, ballisticCoefficient ] = computeSemiMajorAxisChangeDueToDragAveraged( rr0, vv0, jdate0, jdatef, ballisticCoefficient );
else
%     mu = 398600.4415;
    r = sqrt(sum(Xp(1:3,:).^2));
%     v = sqrt(sum(Xp(4:6,:).^2));
%     a = mu./(2.*(mu./r-v.^2./2));
    a = Xp(7,:)/1e6;
    
    [~,apogeeIndices] = findpeaks(r);
    nofFullOrbits = length(apogeeIndices)-1;
    
    if nofFullOrbits >= 2
        % Enough full orbits are available to compute the average
        % semi-major axis.
        
        % A 2nd-order polynomial fit is made. If only two data points
        % are available, then a linear fit is made.
        orderOfPolyFit = min(2,nofFullOrbits-1);
        
        t_mean = zeros(1,nofFullOrbits);
        a_mean = zeros(1,nofFullOrbits);
        for i=1:nofFullOrbits
            t_mean(i)   = mean(tout(apogeeIndices(i):apogeeIndices(i+1)));
            a_mean(i)   = mean(a(apogeeIndices(i):apogeeIndices(i+1)));
        end
        
        % 7 points are preferred to make the fit. Less points are used
        % when less than 9 are available.
        nofPointsForFit = min(9,length(a_mean));
        [semiMajorAxis_Averaged_Fit1,~,mu1]  = polyfit( t_mean(1:nofPointsForFit)', a_mean(1:nofPointsForFit)', orderOfPolyFit );
        [semiMajorAxis_Averaged_Fit2,~,mu2]  = polyfit( t_mean(end-nofPointsForFit+1:end)', a_mean(end-nofPointsForFit+1:end)', orderOfPolyFit );
        
%         figure;
%         hold on;
%         plot(tout, a,'c');
%         plot(t_mean, a_mean,'.b');
%         plot(tout, polyval(semiMajorAxis_Averaged_Fit1,tout),'b');
%         plot(tout, polyval(semiMajorAxis_Averaged_Fit2,tout),'b');
        
        semiMajorAxisChange = polyval(semiMajorAxis_Averaged_Fit2,tout(end),[],mu2)-polyval(semiMajorAxis_Averaged_Fit1,tout(1),[],mu1);
    else
        % Too few full orbits are available to compute the average
        % semi-major axis. Instead we use the osculating semi-major
        % axis.
        semiMajorAxisChange = a(end)-a(1);
    end
end

end

