function [GM, Re, maxdeg, C_grav, S_grav, sF_grav]= initgravitysphericalharmonic( maxdeg )
%  GRAVITYSPHERICALHARMONIC loads gravity model.
% 
%   [GX GY GZ] = INITGRAVITYSPHERICALHARMONIC( MAXDEG ) loads constants
%   and spherical coefficients for EGM2008 up to degree MAXDEG.
%
%   Alternate formats for calling spherical harmonic gravity are:
%   [GX GY GZ] = INITGRAVITYSPHERICALHARMONIC( DEGREE )
%
%   Inputs:
%   DEGREE   :a scalar value specifying the degree and order of the
%            harmonic gravity model. 

load('EGM2008.mat','GM', 'Re', 'degree', 'C', 'S'); % [GM, Re, degree, C, S]

if maxdeg > degree
    warning(['Degree must be less than or equal to maximum degree of EGM2008 model. ' ...
        'Setting maximum degree to the maximum degree of the EGM2008 model.']);
    maxdeg = degree;
end

C_grav = C(1:maxdeg+2,1:maxdeg+2);
S_grav = S(1:maxdeg+2,1:maxdeg+2);

sF_grav = loc_gravLegendre_scaleFactor(maxdeg);

    function [scaleFactor] = loc_gravLegendre_scaleFactor( maxdegree )
        % loc_GRAVLEGENDRE internal function computing normalized associated
        % legendre polynomials, P, via recursion relations for spherical harmonic
        % gravity
        
        scaleFactor = zeros(maxdegree+3, maxdegree+3);
        
        % Seeds for recursion formula
        scaleFactor(1,1) = 0;
        scaleFactor(2,1) = 1;
        scaleFactor(2,2) = 0;
        
        for n = 2:maxdegree+2
            k = n + 1;
            
            for m = 0:n
                p = m + 1;
                % Scale Factor needed for normalization of dUdphi partial derivative
                
                if (n == m)
                    scaleFactor(k,k) = 0;
                elseif (m == 0)
                    scaleFactor(k,p) = sqrt( (n+1)*(n)/2);
                else
                    scaleFactor(k,p) = sqrt( (n+m+1)*(n-m));
                end
            end
        end
    end

end
