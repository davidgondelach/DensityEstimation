%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  Thermospheric Density Estimation Via Two-Line-Element Data Assimilation  %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 24-Sep-2019

%------------- BEGIN CODE --------------

clearvars;
clearvars -global;

%% SETTINGS

% Estimation window
yr      = 2002; % Year
mth     = 8;    % Month
dy      = 1;    % Day
nofDays = 30;   % Number of days

% Reduced-order model
DMDmodel = 'JB2008_1999_2010';  % Name of reduced-order density model
r  = 10;                        % Reduced order

% NORAD catalog IDs of objects used for estimation
selectedObjects = [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]; % 17 objects

%% SET PATHS
spicePath = fullfile('[SPICE TOOLKIT DIRECTORY]','mice');

addpath( 'AstroFunctions' );
addpath( 'Estimation' );
addpath( 'JB2008' );
addpath( 'ROMDensityModels' );
addpath( 'SpaceWeather' );
addpath( 'TLEdata' );
addpath( fullfile('..','..','DMDc_Density_example') );
addpath( fullfile(spicePath,'src','mice') );
addpath( fullfile(spicePath,'lib') );

%% LOAD KERNELS, GRAVITY MODEL, EARTH ORIENTATION PARAMETERS AND SGP4
% Load SPICE
kernelpath  = fullfile('Data','kernel.txt');
loadSPICE(kernelpath);

% Load gravity model
gravitymodel     = 'EGM2008';
gravmodeldegree  = 20;
loadGravityModel( gravitymodel, gravmodeldegree );

% Load Earth orientation parameters (needed for TEME to ECI transformation)
global EOPMat
EOPpath = fullfile('Data','EOP-All.txt');
[ EOPMat ] = inputEOP_Celestrak( EOPpath );

% Setup the SGP4 propagator.
loadSGP4();

%% PERFORM DENSITY ESTIMATION
plotFigures = true;

runDensityEstimationTLE(yr,mth,dy,nofDays,DMDmodel,r,selectedObjects,plotFigures);

%------------- END OF CODE --------------