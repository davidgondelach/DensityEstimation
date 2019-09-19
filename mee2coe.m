function coe = mee2coe(mee)

% convert modified equinoctial elements to classical orbit elements

% input

%  mee(1) = semiparameter (kilometers)
%  mee(2) = f equinoctial element
%  mee(3) = g equinoctial element
%  mee(4) = h equinoctial element
%  mee(5) = k equinoctial element
%  mee(6) = true longitude (radians)

% output

%  coe(1) = semimajor axis (kilometers)
%  coe(2) = eccentricity
%  coe(3) = inclination (radians)
%  coe(4) = argument of periapsis (radians)
%  coe(5) = right ascension of ascending node (radians)
%  coe(6) = true anomaly (radians)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unload modified equinoctial orbital elements

pmee = mee(1);
fmee = mee(2);
gmee = mee(3);
hmee = mee(4);
kmee = mee(5);
lmee = mee(6);

% compute classical orbital elements

tani2s = sqrt(hmee * hmee + kmee * kmee);

% orbital eccentricity

ecc = sqrt(fmee * fmee + gmee * gmee);

% semimajor axis

sma = pmee / (1.0 - ecc * ecc);

% orbital inclination

inc = 2.0 * atan(tani2s);

% right ascension of ascending node

raan = atan2(kmee, hmee);

% argument of periapsis

atopo = atan2(gmee, fmee);

argper = mod(atopo - raan, 2.0 * pi);

% true anomaly

tanom = mod(lmee - atopo, 2.0 * pi);

% load classical orbital element array

coe(1) = sma;
coe(2) = ecc;
coe(3) = inc;
coe(4) = argper;
coe(5) = raan;
coe(6) = tanom;