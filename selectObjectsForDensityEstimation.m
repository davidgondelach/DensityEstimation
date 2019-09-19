
clearvars;
clearvars -global;

addpath( 'AstroFunctions' );
addpath( 'nrlmsise_matlab' );
spicePath = fullfile('/Users','davidgondelach','Documents','mice');
addpath( fullfile(spicePath,'src','mice') );
addpath( fullfile(spicePath,'lib') );

%% Settings
% Start date
yr = 2001;
mth= 10;
dy = 1;
hr=0; mn=0; sc=0;
% Number of days
nofDays = 60;
jd0 = juliandate(yr,mth,dy,0,0,0);
[yrf, mthf, dyf, ~, ~, ~] = datevec(jd0+nofDays-1721058.5);
% jdf = juliandate(yr,mth,dy+nofDays,0,0,0);

%% load standard kernels and reference frames
% Clear cspice memory
cspice_kclear;
% Load SPK, PCK, LSK kernels
kernelpath  = fullfile('Data','kernel.txt');
cspice_furnsh( kernelpath );

global gravconst GM
gravconst   = 6.67259e-20; % [km^3/kg/s^2]

gravitymodel     = 'EGM2008';
gravmodeldegree  = 20;
loadGravityModel( gravitymodel, gravmodeldegree );
GM_kms = GM*1e-9;

SWpath = fullfile('Data','SW-All.txt');
[ SWmatDaily, SWmatMonthlyPred ] = inputSWnrlmsise( SWpath );
[ SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM ] = inputSWtiegcm( SWpath );

% Load Earth orientation parameters (needed for TEME to ECI transformation)
global EOPMat
EOPpath = fullfile('Data','EOP-All.txt');
[ EOPMat ] = inputEOP_Celestrak( EOPpath );

%% Load TLEs
% Setup the SGP4 propagator.
global tumin mu radiusearthkm xke j2 j3 j4 j3oj2 opsmode whichconst
opsmode = 'i';
whichconst = 72;
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );
xpdotp   =  1440.0 / (2.0*pi); % Conversion factor between SGP4 and TLE mean motion units [rev/day]/[rad/min]

% Download TLEs
selectedObjects = [];
maxAlt = 450; % km
[objects] = downloadTLEsForEstimation(yr, mth, dy, yrf, mthf, dyf, maxAlt, selectedObjects);
jdate0TLEs = juliandate(yr,mth,dy,0,0,0);

% Filter TLEs
[objects] = filterTLEs(objects);

nop = length(objects);
objectIDlabels = cell(1, nop);
for i=1:nop
    objectIDlabels(i) = {num2str(objects(i).noradID)};
end

% Analyse TLE quality
[objectDataSorted, covMEEerrors] = checkSelfConsistencyTLEs(objects, mu, jdate0TLEs, objectIDlabels);
%%
maxMeanPosDiff = 5;
minNofTLEs = nofDays;
goodObjects = objectDataSorted(:,2)<maxMeanPosDiff & objectDataSorted(:,3)>=minNofTLEs;
selectedObjects = objectDataSorted(goodObjects,1);
% selectedObjects = selectedObjects(selectedObjects(:,3)>=minNofTLEs,1);
for i=1:length(selectedObjects)
    ID = selectedObjects(i);
    indexObject = find([objects.noradID]==ID);
    newObjects(i) = objects(indexObject);
end
objects = newObjects;
nop = length(objects);
%%
% Plot orbital elements
figure;
for i=1:nop
    subplot(2,3,1); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,[objects(i).satrecs.a],'.'); hold on;
    subplot(2,3,2); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,[objects(i).satrecs.ecco],'.'); hold on;
    subplot(2,3,3); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,rad2deg([objects(i).satrecs.inclo]),'.'); hold on;
    subplot(2,3,4); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,rad2deg([objects(i).satrecs.nodeo]),'.'); hold on;
    subplot(2,3,5); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,rad2deg([objects(i).satrecs.argpo]),'.'); hold on;
    subplot(2,3,6); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,rad2deg([objects(i).satrecs.bstar]),'.-'); hold on;
end
subplot(2,3,1); xlabel('Days since t_0'); ylabel('a [Earth radii]');
subplot(2,3,2); xlabel('Days since t_0'); ylabel('e [-]');
subplot(2,3,3); xlabel('Days since t_0'); ylabel('i [deg]'); legend(objectIDlabels,'Location','northeast');
subplot(2,3,4); xlabel('Days since t_0'); ylabel('\Omega [deg]');
subplot(2,3,5); xlabel('Days since t_0'); ylabel('\omega [deg]');
subplot(2,3,6); xlabel('Days since t_0'); ylabel('Bstar');

figure;
semilogy(objectDataSorted(:,3),objectDataSorted(:,2),'o'); hold on;
semilogy(median(objectDataSorted(:,3)),median(objectDataSorted(:,2)),'+','MarkerSize',10,'LineWidth',2);
