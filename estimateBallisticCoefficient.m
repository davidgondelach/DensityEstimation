function [ estimate ] = estimateBallisticCoefficient( satrec1, satrec2, maxNofIterations, varargin )
%ESTIMATEBALLISTICCOEFFICIENT - This function computes an estimate for the ballistic coefficient based on two TLEs. 
% A first estimate is obtained by comparing the change in semi-major found
% from the TLE with the change in semi-major found by propagating the
% initial state (from SGP4). After the estimate is improved iteratively, by
% propagating again using the new estimates.
%
% Syntax:  [ estimate ] = estimateBallisticCoefficient( satrec1, satrec2, firstGuessBallisticCoefficient, reflectivity, areaToMassRatio, propagatorSettings, maxNofIterations )
%
% Inputs:
%    satrec1                    - Struct with TLE data, first TLE
%    satrec2                    - Struct with TLE data, second TLE
%    reflectivity               - Reflectivity [-]
%    areaToMassRatio            - Area to mass ratio [m^2/kg]
%    propagatorSettings         - Array containing propagator settings
%    maxNofIterations           - Number of allowed iterations, if value <=1 
%                                 then only one improvement is made.
% Optional inputs:
%    firstGuessBallisticCoefficient - First guess of ballistic coefficient,
%                                     if this 0, then a first guess is taken 
%                                     form the Bstar value of the first TLE
%    initialState               - Initial state for orbit propagation (ECI)
%
% Outputs:
%    estimate        - vector of estimates for BC, last elements is final
%                      estimate
%
% Other m-files required: computeSemiMajorAxisChangeBetweenTwoTLE, computeSemiMajorAxisChangeDueToDrag
%
% Reference: Saunders, A., Swinerd, G.G., Lewis, H.G., 2012. Deriving
% accurate satellite ballistic coefficients from two-line element data

% Author: David Gondelach
% Astronautics, University of Southampton
% July 2015; Last revision: 31-July-2015

%% TLE epochs
jdateTLE1 = satrec1.jdsatepoch;
jdateTLE2 = satrec2.jdsatepoch;

%% First guess for ballistic coefficient
if nargin<4
    % If no first guess is provided, then use Bstar as first guess.
    % BC = 2*Bstar/(Re*rho0) = 12.741621*Bstar (Re=6378.135, rho0=2.461E-5) [Ref: Vallado, p.106]
    ballisticCoefficientGuesses(1) = 12.741621*satrec1.bstar;
    
    % Check if bstar is negative
    if ballisticCoefficientGuesses(1) < 1.0e-05
        % If BC guess is negative, then use the assumed area-to-mass ratio
        % and a drag coefficient of 2.2.
%         ballisticCoefficientGuesses(1) = 2.2 * areaToMassRatio;
        ballisticCoefficientGuesses(1) = 0.001;
    end
else
    % Use provided first guess for BC
    ballisticCoefficientGuesses(1) = varargin{1};
end

%% Spacecraft initial state 
if nargin==5 && length(varargin)>1
    initialState = varargin{2};
    rr0_ECI      = initialState(1:3);
    vv0_ECI      = initialState(4:6);
else
    % If no initial state is provided, then use SGP4 state of TLE1
    [~, rr0_TEME, vv0_TEME] = sgp4 (satrec1, 0);
    [rr0_ECI, vv0_ECI] = convertTEMEtoJ2000(rr0_TEME', vv0_TEME', jdateTLE1);
end

%% Compute semi-major change according to TLE data
semiMajorAxisChangeTLE = computeSemiMajorAxisChangeBetweenTwoTLE( satrec1, satrec2 );

% Compute semi-major axis change due to drag and difference with TLE
% The value of the BC guess is changed inside 
% 'computeSemiMajorAxisChangeDueToDrag' if the value was so high that 
% the object re-entered during propagating.

[semiMajorAxisChangePropagation, ballisticCoefficientGuesses(1)] = ...
    computeSemiMajorAxisChangeDueToDragAveraged( rr0_ECI, vv0_ECI, jdateTLE1, jdateTLE2, ballisticCoefficientGuesses(1) );

diffSemiMajorAxisChange(1) = semiMajorAxisChangePropagation - semiMajorAxisChangeTLE;

ballisticCoefficientGuesses(2) = ballisticCoefficientGuesses(1) * semiMajorAxisChangeTLE / semiMajorAxisChangePropagation;

%% Estimate ballistic coefficient
i = 1;
while( abs(diffSemiMajorAxisChange(i)/semiMajorAxisChangeTLE) > 1E-5 && i < maxNofIterations )
    i = i + 1;
    
    if isinf(ballisticCoefficientGuesses(i))
        disp('error');
    end
    
    % Compute semi-major axis change due to drag and difference with TLE
    [semiMajorAxisChangePropagation, ballisticCoefficientGuesses(i)] = ...
        computeSemiMajorAxisChangeDueToDragAveraged( rr0_ECI, vv0_ECI, jdateTLE1, jdateTLE2, ballisticCoefficientGuesses(i) );
    
    diffSemiMajorAxisChange(i) = semiMajorAxisChangePropagation - semiMajorAxisChangeTLE;
    
    % Update ballistic coefficient estimate
    ballisticCoefficientGuesses(i+1) = ballisticCoefficientGuesses(i) ...
        - diffSemiMajorAxisChange(i) * ( ballisticCoefficientGuesses(i) - ballisticCoefficientGuesses(i-1) ) ...
                                     / ( diffSemiMajorAxisChange(i) - diffSemiMajorAxisChange(i-1) );
                                 
end

estimate = ballisticCoefficientGuesses;
end