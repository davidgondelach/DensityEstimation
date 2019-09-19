function [ semiMajorAxisChange ] = computeSemiMajorAxisChangeBetweenTwoTLE( satrec1, satrec2 )
%COMPUTESEMIMAJORAXISCHANGEBETWEENTWOTLE - This function computes the
%change in semi-major axis between two TLEs using the mean mean major
%available in the TLEs.
% 
% Syntax:  [ semiMajorAxisChange ] = computeSemiMajorAxisChangeBetweenTwoTLE( satrec1, satrec2 )
%
% Inputs:
%    satrec1                    - Struct with TLE data, first TLE
%    satrec2                    - Struct with TLE data, second TLE
%
% Outputs:
%    semiMajorAxisChange        - Change in semi-major axis [km]

% Author: David Gondelach
% Astronautics, University of Southampton
% July 2015; Last revision: 31-July-2015

global mu;

%% Compute semi-major change according to TLE data
% Semi-major axis at first TLE
meanMotion1 = satrec1.no/60; % rad/s
semiMajorAxisTLE1 = (mu/meanMotion1^2)^(1/3);

% Semi-major axis at second TLE
meanMotion2 = satrec2.no/60; % rad/s
semiMajorAxisTLE2 = (mu/meanMotion2^2)^(1/3);

% Change in semi-major axis between TLEs
semiMajorAxisChange = semiMajorAxisTLE2 - semiMajorAxisTLE1;

end

