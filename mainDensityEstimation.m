% clearvars -except DMDmodel selectedObjects;
clearvars;
clearvars -global;
% close all;

%% SETTINGS
SUPERCLOUD = false;

yr      = 2003;
mth     = 4;
dy      = 1;
nofDays = 22;

DMDmodel = 'JB2008_1999_2010';
r   = 10;

% selectedObjects = [16111;38710;35867;13770;15354;39452;11822;42874;07337;08744;12138;12388]; %year=2017; 42874,42880 are Planet sat with good BC
% selectedObjects = [00022;00060;00063;00165;02611;14483;22875;23853;26929]; %year= August 2002; 02611: decay @ 2002-11-30
% 
% selectedObjects = [02016;02153;02622;04221;06073;20774;23278;26405;27391;27392]; % Less accurate: 03553;01370;00614; High BC: 00750; No TLE in 2004: 01808;02129;04330;25233;
% selectedObjects = [02016;02153;02622;04221;20774;23278;07337;08744;12138;12388]; % Without GRACE: 27391;27392; CHAMP: 26405; Less accurate: 03553;01370;00614;06073; High BC: 00750;
% selectedObjects = [00614;27391;26405;20774;23278;07337;08744;12138;12388]; % InclGRACE, InclCHAMP;  BadUKF: 04221; Perigee>450: 02016;02153;02622; GRACE:27392; Less accurate: 03553;01370;06073; High BC: 00750;
% % selectedObjects = [00614;04221;26405;20774;23278;07337;08744;12138;12388;06073]; % NoGRACE, InclCHAMP; BadUKF: 04221; Perigee>450: 27391;02016;02153;02622; GRACE:27391;27392; Less accurate: 03553;01370; High BC: 00750;
% selectedObjects = [02016;02153;02622;04221;20774;23278;07337;08744;12138;12388;27391;26405]; % JB2008 on SuperCloud
% selectedObjects = [00614;04221;26405;20774;23278;07337;08744;12138;12388;06073]; % TIEGCM on SuperCloud
% % selectedObjects = [02016;02153;02622;04221;20774;23278;07337;08744;12138;12388]; % 2005
% 
% selectedObjects = [00614;04221;06073;07337;08744;12138;12388;20774;23278]; % hp <450km
% selectedObjects = [00614;04221;06073;07337;08744;12138;12388;20774;23278;26405]; % hp <450km incl CHAMP (official TIEGCM)
% % selectedObjects = [00614;02016;02153;02622;04221;06073;07337;08744;12138;12388;20774;23278]; % hp >450km
% selectedObjects = [00614;02016;02153;02622;04221;06073;07337;08744;12138;12388;20774;23278;27391;27392;26405]; % hp >450km incl CHAMP & GRACE (official JB&MSISE)
% selectedObjects = [00614;02622;04221;06073;07337;08744;12138;12388;20774;23278;27391;27392;26405]; % hp >450km incl CHAMP & GRACE, without 02016,02153
% selectedObjects = [00060;00063;00165;00614;02016;02153;02622;04221;06073;07337;08744;12138;12388;14483;20774;23278;26405;27391;27392]; % Full set 2007, without 00750
% % selectedObjects = [00060;00614;02622;04221;06073;07337;08744;12138;12388;14483;20774;23278;26405]; % Full set 2007, hp >450km without 00750
% selectedObjects = [63;165;614;2153;2622;4221;6073;8744;14483;20774;23278;26405;27392]; % Best spread of objects
% % selectedObjects = [614;4221;6073;8744;14483;20774;23278;26405]; % Best spread of objects, hp <450km
% 
% % selectedObjects = [60;63;165;614;750;2016;2153;2622;4221;6073;12138;14483;26405;27391;27392]; % More and better spread objects
selectedObjects = [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]; % 17 objects

%% SET PATHS

global resultsDirPath GPSdataPath
if SUPERCLOUD
    spicePath = fullfile('..','..','SPICE','mice');
    resultsDirPath = fullfile('Results');
    GPSdataPath = fullfile('GPSdata');
else
    spicePath = fullfile('/Users','davidgondelach','Documents','mice');
    resultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/';
    GPSdataPath = fullfile('/Users/davidgondelach/Documents/PostDoc','GPSdata');
end

addpath( 'AstroFunctions' );
addpath( 'nrlmsise_matlab' );
addpath( 'JacchiaBowmanAtmosphericDensityModel' );
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
gpsFitAndCompare = false;
runDensityEstimationTLE;
% runDensityEstimationTLEandDensityData;
% compareUKFestimateWithGPS;
% runDensityEstimationTLE(yr,mth,dy,nofDays,DMDmodel,r,selectedObjects,plotFigures,gpsFitAndCompare);
