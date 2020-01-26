function [aa_grav_x, aa_grav_y, aa_grav_z] = computeEarthGravitationalAcceleration( rr_ecef )
%  COMPUTEEARTHGRAVITATIONALACCELERATION - Compute Earth gravitational 
%   acceleration using spherical harmonic expansion of Earth's geopotential
%   in Earth-Centered Earth-Fixed coordinates [m/s^2].
%
%   [aa_grav_x, aa_grav_y, aa_grav_z] = computeEarthGravitationalAcceleration( rr_ecef )
%
%   Inputs:
%   rr_ecef - position in Earth-Centered Earth-Fixed (ECEF) coordinates
%            in meters [Mx3].
%
%   Outputs:
%   aa_grav_x - Earth gravitational acceleration in the x-direction of the
%          Earth-Centered Earth-Fixed coordinates in meters per second
%          squared [Mx1].
%   aa_grav_y - Earth gravitational acceleration in the y-direction of the
%          Earth-Centered Earth-Fixed coordinates in meters per second
%          squared [Mx1].
%   aa_grav_z - Earth gravitational acceleration in the z-direction of the
%          Earth-Centered Earth-Fixed coordinates in meters per second
%          squared [Mx1].

%   References:
%   Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 2001.

% Check that number of columns is equal to 3.
if (size( rr_ecef, 2) ~= 3)
    error('Input matrix has wrong dimension. Matrix must have 3 columns.');
end

% Load geopotential constants and coefficients
global GM Re C_gravmodel S_gravmodel gravdegree sF_gravmod;

C = C_gravmodel;
S = S_gravmodel;

% Compute geocentric radius
r = sqrt( sum( rr_ecef.^2, 2 ));

% Check if geocentric radius is less than equatorial radius
if r < Re
    error('Radial position is less than Earth equatorial radius, %g.', Re);
end

% Compute geocentric latitude
phic = asin( rr_ecef(:,3)./ r );

% Compute lambda
lambda = atan2( rr_ecef(:,2), rr_ecef(:,1) );

smlambda = zeros( size(rr_ecef,1), gravdegree );
cmlambda = zeros( size(rr_ecef,1), gravdegree );

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
[P] = computeLegendrePolynomials( phic, gravdegree );

scaleFactor = sF_gravmod;

% Compute gravity in ECEF coordinates
[aa_grav_x, aa_grav_y, aa_grav_z] = computeGravity( rr_ecef, gravdegree, P, ...
    C( 1:gravdegree+1, 1:gravdegree+1 ), S( 1:gravdegree+1, 1:gravdegree+1 ), ...
    smlambda, cmlambda, GM, Re, r,scaleFactor );


%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function [P] = computeLegendrePolynomials( phi, maxdeg )
        % computeLegendrePolynomials - compute normalized associated legendre polynomials P
        
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
            
            for m_ = 0:n
                p = m_ + 1;
                % Compute normalized associated legendre polynomials, P, via recursion relations
                % Scale Factor needed for normalization of dUdphi partial derivative
                
                if (n == m_)
                    P(k,k,:) = sqrt(2*n+1)/sqrt(2*n)*sphi.*reshape(P(k-1,k-1,:),size(phi));
                elseif (m_ == 0)
                    P(k,p,:) = (sqrt(2*n+1)/n)*(sqrt(2*n-1)*cphi.*reshape(P(k-1,p,:),size(phi)) - (n-1)/sqrt(2*n-3)*reshape(P(k-2,p,:),size(phi)));
                else
                    P(k,p,:) = sqrt(2*n+1)/(sqrt(n+m_)*sqrt(n-m_))*(sqrt(2*n-1)*cphi.*reshape(P(k-1,p,:),size(phi)) - sqrt(n+m_-1)*sqrt(n-m_-1)/sqrt(2*n-3)*reshape(P(k-2,p,:),size(phi)));
                end
            end
        end
    end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function [aa_grav_x, aa_grav_y, aa_grav_z] = computeGravity(p,maxdeg,P,C,S,smlambda,cmlambda,GM,Re,r,scaleFactor)
        % computeGravity - compute Earth gravitational acceleration in Earth-centered Earth-fixed (ECEF) coordinates
        
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
            for m_ = 0:n
                j = m_+1;
                dUdrSumM      = dUdrSumM + reshape(P(k,j,:),size(r)).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j));
                dUdphiSumM    = dUdphiSumM + ( (reshape(P(k,j+1,:),size(r)).*scaleFactor(k,j,:)) - p(:,3)./(sqrt(p(:,1).^2 + p(:,2).^2)).*m_.*reshape(P(k,j,:),size(r))).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j));
                dUdlambdaSumM = dUdlambdaSumM + m_*reshape(P(k,j,:), size(r)).*(S(k,j).*cmlambda(:,j) - C(k,j).*smlambda(:,j));
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
        aa_grav_x = ((1./r).*dUdr - (p(:,3)./(r.*r.*sqrt(p(:,1).^2 + p(:,2).^2))).*dUdphi).*p(:,1) ...
            - (dUdlambda./(p(:,1).^2 + p(:,2).^2)).*p(:,2);
        aa_grav_y = ((1./r).*dUdr - (p(:,3)./(r.*r.*sqrt(p(:,1).^2 + p(:,2).^2))).*dUdphi).*p(:,2) ...
            + (dUdlambda./(p(:,1).^2 + p(:,2).^2)).*p(:,1);
        aa_grav_z = (1./r).*dUdr.*p(:,3) + ((sqrt(p(:,1).^2 + p(:,2).^2))./(r.*r)).*dUdphi;
        
        % special case for poles
        atPole = abs(atan2(p(:,3),sqrt(p(:,1).^2 + p(:,2).^2)))==pi/2;
        if any(atPole)
            aa_grav_x(atPole) = 0;
            aa_grav_y(atPole) = 0;
            aa_grav_z(atPole) = (1./r(atPole)).*dUdr(atPole).*p((atPole),3);
        end
        
    end

end
