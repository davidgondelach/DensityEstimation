function [gx, gy, gz]= rungravitysphericalharmonic( p )
%  GRAVITYSPHERICALHARMONIC Implement a spherical harmonic representation
%   of planetary gravity. 
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P ) implements the mathematical
%   representation of spherical harmonic planetary gravity based on
%   planetary gravitational potential. Using P, a M-by-3 array of
%   Planet-Centered Planet-Fixed coordinates, GX, GY and GZ, arrays of M 
%   gravity values in the x-axis, y-axis and z-axis of the Planet-Centered
%   Planet-Fixed coordinates are calculated for planet using 120th degree 
%   and order spherical coefficients for EGM2008 by default. 
%
%   Alternate formats for calling spherical harmonic gravity are:
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, DEGREE )   
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL )   
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL, DEGREE )   
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL, DEGREE, ACTION )   
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, 'Custom', DEGREE, {DATAFILE DFREADER}, ACTION )   
%
%   Inputs for spherical harmonic gravity are:
%   P        :a M-by-3 array of Planet-Centered Planet-Fixed coordinates in
%            meters where the z-axis is positive towards the North Pole. For
%            Earth this would be ECEF coordinates.
%
%   Output calculated for the spherical harmonic gravity includes:
%   GX     :an array of M gravity values in the x-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared.
%   GY     :an array of M gravity values in the y-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared. 
%   GZ     :an array of M gravity values in the z-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared. 
%
%   Limitations:                                                           
%
%   This function has the limitations of excluding the centrifugal effects
%   of planetary rotation, and the effects of a precessing reference frame.
%
%   Spherical harmonic gravity model is valid for radial positions greater
%   than the planet equatorial radius.  Using it near or at the planetary
%   surface can probably be done with negligible error.  The spherical
%   harmonic gravity model is not valid for radial positions less than
%   planetary surface. 

%   References:  
%   [1] Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
%       McGraw-Hill, New York, 1997.  
%   [2] NIMA TR8350.2: "Department of Defense World Geodetic System 1984,
%       Its Definition and Relationship with Local Geodetic Systems." 
%   [3] Konopliv, A. S., S. W. Asmar, E. Carranza, W. L. Sjogen, D. N.
%       Yuan., "Recent Gravity Models as a Result of the Lunar Prospector
%       Mission", Icarus, Vol. 150, no. 1, pp 1?18, 2001.                    
%   [4] Lemoine, F. G., D. E. Smith, D.D. Rowlands, M.T. Zuber, G. A.
%       Neumann, and D. S. Chinn, "An improved solution of the gravity
%       field of Mars (GMM-2B) from Mars Global Surveyor", J. Geophys. Res.,
%       Vol. 106, No. E10, pp 23359-23376, October 25, 2001.   
%   [5] Kenyon S., J. Factor, N. Pavlis, and S. Holmes, "Towards the Next
%       Earth Gravitational Model", Society of Exploration Geophysicists
%       77th Annual Meeting, San Antonio, Texas, September 23-28, 2007.
%   [6] Pavlis, N.K., S.A. Holmes, S.C. Kenyon, and J.K. Factor, "An Earth
%       Gravitational Model to Degree 2160: EGM2008", presented at the 2008
%       General Assembly of the European Geosciences Union, Vienna,
%       Austria, April 13-18, 2008. 
%   [7] Grueber, T., and A. Kohl, "Validation of the EGM2008 Gravity Field
%       with GPS-Leveling and Oceanographic Analyses", presented at the IAG
%       International Symposium on Gravity, Geoid & Earth Observation 2008,
%       Chania, Greece, June 23-27, 2008.
%   [8] F?rste, C., Flechtner, F., Schmidt, R., K?nig, R., Meyer, U.,
%       Stubenvoll, R., Rothacher, M., Barthelmes, F., Neumayer, H.,
%       Biancale, R., Bruinsma, S., Lemoine, J.M., Loyer, S., "A Mean
%       Global Gravity Field Model From the Combination of Satellite
%       Mission and Altimetry/Gravmetry Surface Data - EIGEN-GL04C",
%       Geophysical Research Abstracts, Vol. 8, 03462, 2006 
%       http://icgem.gfz-potsdam.de/ICGEM/

checkinputs( );

global GM Re C_gravmodel S_gravmodel gravdegree sF_gravmod;

C = C_gravmodel;
S = S_gravmodel;

% Compute geocentric radius
r = sqrt( sum( p.^2, 2 ));
                                                                           
% Check if geocentric radius is less than equatorial (reference) radius
if r < Re
    error('Radial position is less than equatorial radius of planetary model, %g.', Re);
end

% Compute geocentric latitude
phic = asin( p(:,3)./ r );

% Compute lambda                                                           
lambda = atan2( p(:,2), p(:,1) );

smlambda = zeros( size(p,1), gravdegree );
cmlambda = zeros( size(p,1), gravdegree );

slambda = sin(lambda);
clambda = cos(lambda);
smlambda(:,1) = 0;
cmlambda(:,1) = 1;
smlambda(:,2) = slambda;
cmlambda(:,2) = clambda;

for m=3:gravdegree+1
    smlambda(:,m) = 2.0.*clambda.*smlambda(:, m-1) - smlambda(:, m-2);
    cmlambda(:,m) = 2.0.*clambda.*cmlambda(:, m-1) - cmlambda(:, m-2);
end


% Compute normalized associated legendre polynomials
[P] = loc_gravLegendre( phic, gravdegree );

scaleFactor = sF_gravmod;

% Compute gravity in ECEF coordinates
[gx gy gz] = loc_gravityPCPF( p, gravdegree, P, C( 1:gravdegree+1, 1:gravdegree+1 ), ...
                                  S( 1:gravdegree+1, 1:gravdegree+1 ), smlambda, ...
                                  cmlambda, GM, Re, r,scaleFactor );
                              
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs( )
        if ~isnumeric( p )
            error(message('aero:gravitysphericalharmonic:notNumeric'));
        end
        
        if (size( p, 2) ~= 3)
            error(message('aero:gravitysphericalharmonic:wrongDimension'));
        end
   end

end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [P] = loc_gravLegendre( phi, maxdeg )
% loc_GRAVLEGENDRE internal function computing normalized associated 
% legendre polynomials, P, via recursion relations for spherical harmonic
% gravity 

P = zeros(maxdeg+3, maxdeg+3, length(phi));
cphi = cos(pi/2-phi);
sphi = sin(pi/2-phi);

% force numerically zero values to be exactly zero
cphi(abs(cphi)<=eps) = 0;
sphi(abs(sphi)<=eps) = 0;
 
% Seeds for recursion formula
P(1,1,:) = 1;            % n = 0, m = 0;
P(2,1,:) = sqrt(3)*cphi; % n = 1, m = 0;
P(2,2,:) = sqrt(3)*sphi; % n = 1, m = 1;

for n = 2:maxdeg+2
    k = n + 1;
    
%     sqrt2np1 = sqrt(2*n+1);
%     sqrt2nm1 = sqrt(2*n-1);
%     sqrt2nm3 = sqrt(2*n-3);
    
    for m = 0:n
        p = m + 1;
        % Compute normalized associated legendre polynomials, P, via recursion relations 
        % Scale Factor needed for normalization of dUdphi partial derivative
                
        if (n == m)           
            P(k,k,:) = sqrt(2*n+1)/sqrt(2*n)*sphi.*reshape(P(k-1,k-1,:),size(phi));
%             P(k,k) = sqrt2np1/sqrt(2*n)*sphi.*P(k-1,k-1);
        elseif (m == 0)
            P(k,p,:) = (sqrt(2*n+1)/n)*(sqrt(2*n-1)*cphi.*reshape(P(k-1,p,:),size(phi)) - (n-1)/sqrt(2*n-3)*reshape(P(k-2,p,:),size(phi)));
%             P(k,p) = (sqrt2np1/n)*(sqrt2nm1*cphi.*P(k-1,p) - (n-1)/sqrt(2*n-3)*P(k-2,p));
        else
            P(k,p,:) = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1)*cphi.*reshape(P(k-1,p,:),size(phi)) - sqrt(n+m-1)*sqrt(n-m-1)/sqrt(2*n-3)*reshape(P(k-2,p,:),size(phi)));
%             P(k,p) = sqrt2np1/(sqrt(n+m)*sqrt(n-m))*(sqrt2nm1*cphi.*P(k-1,p) - sqrt(n+m-1)*sqrt(n-m-1)/sqrt2nm3*P(k-2,p) );
        end
    end
end
end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [gx gy gz] = loc_gravityPCPF(p,maxdeg,P,C,S,smlambda,cmlambda,GM,Re,r,scaleFactor)
% loc_GRAVITYPCPF internal function computing gravity in planet-centered
% planet-fixed (PCEF) coordinates using PCPF position, desired
% degree/order, normalized associated legendre polynomials, normalized
% spherical harmonic coefficients, trigonometric functions of geocentric
% latitude and longitude, planetary constants, and radius to center of
% planet. Units are MKS.

rRatio   = Re./r;
rRatio_n = rRatio;

% initialize summation of gravity in radial coordinates
dUdrSumN      = 1;
dUdphiSumN    = 0;
dUdlambdaSumN = 0;

% summation of gravity in radial coordinates
for n = 2:maxdeg
    k = n+1;
    rRatio_n      = rRatio_n.*rRatio;
    dUdrSumM      = 0;
    dUdphiSumM    = 0;
    dUdlambdaSumM = 0;
    for m = 0:n
        j = m+1;
%         dUdrSumM      = dUdrSumM + P(k,j).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j)); 
        dUdrSumM      = dUdrSumM + reshape(P(k,j,:),size(r)).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j)); 
%         dUdphiSumM    = dUdphiSumM + ( (P(k,j+1).*scaleFactor(k,j)) - p(:,3)./(sqrt(p(:,1).^2 + p(:,2).^2)).*m.*P(k,j)).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j)); 
        dUdphiSumM    = dUdphiSumM + ( (reshape(P(k,j+1,:),size(r)).*scaleFactor(k,j,:)) - p(:,3)./(sqrt(p(:,1).^2 + p(:,2).^2)).*m.*reshape(P(k,j,:),size(r))).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j));
%         dUdlambdaSumM = dUdlambdaSumM + m*P(k,j,:).*(S(k,j).*cmlambda(:,j) - C(k,j).*smlambda(:,j)); 
        dUdlambdaSumM = dUdlambdaSumM + m*reshape(P(k,j,:), size(r)).*(S(k,j).*cmlambda(:,j) - C(k,j).*smlambda(:,j));
    end
    dUdrSumN      = dUdrSumN      + dUdrSumM.*rRatio_n.*k;
    dUdphiSumN    = dUdphiSumN    + dUdphiSumM.*rRatio_n;
    dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM.*rRatio_n;
end

% gravity in spherical coordinates
dUdr      = -GM./(r.*r).*dUdrSumN;
dUdphi    =  GM./r.*dUdphiSumN;
dUdlambda =  GM./r.*dUdlambdaSumN;

% gravity in ECEF coordinates
gx = ((1./r).*dUdr - (p(:,3)./(r.*r.*sqrt(p(:,1).^2 + p(:,2).^2))).*dUdphi).*p(:,1) ...
      - (dUdlambda./(p(:,1).^2 + p(:,2).^2)).*p(:,2); 
gy = ((1./r).*dUdr - (p(:,3)./(r.*r.*sqrt(p(:,1).^2 + p(:,2).^2))).*dUdphi).*p(:,2) ...
      + (dUdlambda./(p(:,1).^2 + p(:,2).^2)).*p(:,1); 
gz = (1./r).*dUdr.*p(:,3) + ((sqrt(p(:,1).^2 + p(:,2).^2))./(r.*r)).*dUdphi;

% special case for poles
atPole = abs(atan2(p(:,3),sqrt(p(:,1).^2 + p(:,2).^2)))==pi/2;
if any(atPole)
    gx(atPole) = 0;
    gy(atPole) = 0;
    gz(atPole) = (1./r(atPole)).*dUdr(atPole).*p((atPole),3);
end

end


