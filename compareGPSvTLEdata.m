% % Download GPS files
% for i=1:9
% websave(['gps_data_2018010',num2str(i),'.zip'],['http://ephemerides.planet-labs.com/gps_data_2018010',num2str(i),'.zip']);
% end
% for i=10:31
% websave(['gps_data_201801',num2str(i),'.zip'],['http://ephemerides.planet-labs.com/gps_data_201801',num2str(i),'.zip']);
% end

clearvars;

spicePath = fullfile('/Users','davidgondelach','Documents','mice');
addpath( 'AstroFunctions' );
addpath( fullfile(spicePath,'src','mice') );
addpath( fullfile(spicePath,'lib') );

% Constants
GM = 398600.4415;
% Setup the SGP4 propagator.
global tumin mu radiusearthkm xke j2 j3 j4 j3oj2 opsmode whichconst
opsmode = 'i';
whichconst = 72;
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );
xpdotp   =  1440.0 / (2.0*pi); % Conversion factor between SGP4 and TLE mean motion units [rev/day]/[rad/min]
% Clear cspice memory
cspice_kclear;
% Load SPK, PCK, LSK kernels
kernelpath  = fullfile('Data','kernel.txt');
cspice_furnsh( kernelpath );
% Load Earth orientation parameters (needed for TEME to ECI transformation)
global EOPMat
% EOPpath = fullfile(pwd,'AstroFunctions','finals-IAU80.eop');
EOPpath = fullfile(pwd,'Data','EOP-All.txt');
[ EOPMat ] = inputEOP_Celestrak( EOPpath );

%%
objectPlanetID = '0E0E'; objectNORADID = 41609;
% objectPlanetID = '0F4A'; objectNORADID = 42861;
% objectPlanetID = '0F06'; objectNORADID = 42995;

objectGPSfiles = {};
objectGPSpaths = {};
for i=1:31
    day = num2str(i,'%02d');
    GPSfiledirectory = fullfile('/Users','davidgondelach','Documents','PostDoc','GPSdata',['gps_data_201801',day]);
    cd(GPSfiledirectory);
    listing = dir;

    % Get GPS data
    objectIndices = startsWith({listing(:).name},objectPlanetID);
    objectGPSfiles(end+1:end+sum(objectIndices)) = {listing(objectIndices).name};
    objectGPSpaths(end+1:end+sum(objectIndices)) = {GPSfiledirectory};
end

cd(fullfile('/Users','davidgondelach','Dropbox (MIT)','Research-ROM','code','sqrt-ukf-rom-gps-2019'));

% Get TLE data
filename = fullfile('/Users','davidgondelach','Dropbox (MIT)','Research-ROM','code','sqrt-ukf-rom-gps-2019','TLEdata',[num2str(objectNORADID),'.tle']);
[objectTLEs] = getTLEs(filename);

%%
for i=1:length(objectGPSfiles)
    % Read GPS message and get state in J2000
    filename = objectGPSfiles{i};
    GPSfiledirectory = objectGPSpaths{i};
    [et,xx_GPS_J2000] = readGPSmessage(fullfile(GPSfiledirectory,filename));
    jdatestr    = cspice_et2utc( et, 'J', 12 );
    jdate       = str2double(jdatestr(4:end)); % Cut leading 'JD ' off from string
    
    % Find nearest-newer TLE and get state in J2000
    satrecIndex(i) = find([objectTLEs.satrecs.jdsatepoch] >= jdate, 1, 'first');
    diffEpochMinutes(i) = (jdate - objectTLEs.satrecs(satrecIndex(i)).jdsatepoch) * 24*60;
    [~, rr_TLE_TEME ,vv_TLE_TEME] = sgp4( objectTLEs.satrecs(satrecIndex(i)), diffEpochMinutes(i) );
    [rr_TLE_J2000, vv_TLE_J2000] = convertTEMEtoJ2000(rr_TLE_TEME', vv_TLE_TEME', jdate);
    coeTLE(i,:) = pv2po(rr_TLE_J2000, vv_TLE_J2000, GM);
    
    rr_Diff = rr_TLE_J2000 - xx_GPS_J2000(1:3);
    
    diffR(i) = norm(rr_Diff);
    jdatesGPS(i) = jdate;
    
    [cart2rtnMatrix] = computeCart2RTNMatrix(xx_GPS_J2000(1:3), xx_GPS_J2000(4:6));
    rr_Diff_RTN(i,:) = cart2rtnMatrix*rr_Diff;
end

[UsedTLEs] = find([objectTLEs.satrecs.jdsatepoch] >= min(jdatesGPS) & [objectTLEs.satrecs.jdsatepoch] <= max(jdatesGPS));

figure;
plot(diffEpochMinutes,diffR,'.');
xlabel('Minutes since TLE epoch'); ylabel('Position error [km]');

% figure;
% plot(jdatesGPS,diffR,'.');
% xlabel('Julian date'); ylabel('Position error [km]');
% 
figure;
plot(jdatesGPS,diffR);
xlabel('Julian date'); ylabel('Position error [km]');

figure;
plot(jdatesGPS,rr_Diff_RTN(:,1),jdatesGPS,rr_Diff_RTN(:,2),jdatesGPS,rr_Diff_RTN(:,3)); hold on;
plot([objectTLEs.satrecs(UsedTLEs).jdsatepoch],0,'b.');
xlabel('Julian date'); ylabel('Position error [km]');
legend('Radial','Transverse','Normal');

meanRadialErr = mean(rr_Diff_RTN(:,1))
stdRadialErr = std(rr_Diff_RTN(:,1))
meanTransvErr = mean(rr_Diff_RTN(:,2))
stdTransvErr = std(rr_Diff_RTN(:,2))
meanNormalErr = mean(rr_Diff_RTN(:,3))
stdNormalErr = std(rr_Diff_RTN(:,3))

figure;
% plot(jdatesGPS,coeTLE(:,1));
plot([objectTLEs.satrecs(UsedTLEs).jdsatepoch],[objectTLEs.satrecs(UsedTLEs).a]*radiusearthkm,'.');
xlabel('Julian date'); ylabel('a_{TLE} [km]');