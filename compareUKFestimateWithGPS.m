% load('ukf_rom_tle_MSISE2008_workspace_20170501_8d_13obj_190614130350.mat');

etf_ukf = et0 + time(end);

% noradID, planetID, BCestimate
testObjects(1,:) = {42874, '1056', 0.0154841705917428};
testObjects(2,:) = {42880, '1047', 0.0149184720940771};
testObjects(3,:) = {41609, '0E0E', 0.0285353336673982};
testObjects(4,:) = {42861, '0F4A', 0.0352047559783136};
testObjects(5,:) = {41970, '100B', 0.0403420104279107};

testID = 1;
noradID = cell2mat(testObjects(testID,1));
planetID = char(testObjects(testID,2));
BC0 = cell2mat(testObjects(testID,3));

% Load GPS data
GPSobsWindowDays = 10;
GPSobsWindowSeconds = GPSobsWindowDays*86400;
[xx_GPS_obs,et_gps_obs,jdate_gps_obs] = getGPSdataET(GPSdataPath,planetID,etf_ukf,GPSobsWindowSeconds);

% Initial state guess and observations
x0 = [xx_GPS_obs(:,1); BC0];
Obs = xx_GPS_obs(1:3,:);
et0_gps  = et_gps_obs(1);
time_obs = et_gps_obs-et0_gps;


%%
romState = X_est(end-r+1:end,end);

r_pred = r;
F_U_pred = F_U;
AC_pred = AC;
BC_pred = BC;
romState_pred = romState;
[romStateTime] = predictROMdensity(romState_pred,etf_ukf,et_gps_obs(end),AC_pred,BC_pred,DMDmodel,SWmatDaily,SWmatMonthlyPred,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);


%% PREDICT ORBIT USING DENSITY-UKF RESULT
objNum = find(selectedObjects==noradID);

% Extract estimated final state of object of interest
[rf_obj_ukf_rom,vf_obj_ukf_rom] = ep2pv(X_est(svs*(objNum-1)+1:svs*(objNum-1)+6,end),GM_kms);
BCf_obj_ukf_rom = X_est(svs*objNum,end)/1000;
xf_obj_ukf_rom = [rf_obj_ukf_rom;vf_obj_ukf_rom;BCf_obj_ukf_rom];
time_obs_2 = et_gps_obs-etf_ukf;

% Predict trajectory using ROM density model
% [x_obj_ukf_rom_pred_rom] = propagateStatePosVel_FullGravDrag_ROM(xf_obj_ukf_rom,[0,time_obs_2],etf_ukf,romStateTime,r_pred,F_U_pred,M_U,includeSRPMoonSunPert);
% residuals_obj_ukf_rom_pred_rom = xx_GPS_obs(1:3,:)-x_obj_ukf_rom_pred_rom(1:3,2:end);

settingsROM.drag = 1;
settingsROM.thirdbody = 1;
densityModelData.romStateTime=romStateTime; densityModelData.r=r_pred; densityModelData.F_U=F_U_pred; densityModelData.M_U=M_U;
% settingsJB2000=settings;settingsJB2000.drag=3;
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[~,x_obj_ukf_rom_pred_rom]=ode113(@(t,x) computeAcceleration(t,x,etf_ukf,settingsROM,densityModelData),[0,time_obs_2],xf_obj_ukf_rom,opts);
residuals_obj_ukf_rom_pred_rom = xx_GPS_obs(1:3,:)-x_obj_ukf_rom_pred_rom(2:end,1:3)';

% Predict trajectory using NRLMSISE-00 density model
% [x_obj_ukf_rom_pred_msise] = propagateStatePosVel_FullGravDrag_NRLMSISE(xf_obj_ukf_rom,[0,time_obs_2],etf_ukf,includeSRPMoonSunPert);
% residuals_obj_ukf_rom_pred_msise = xx_GPS_obs(1:3,:)-x_obj_ukf_rom_pred_msise(1:3,2:end);

settingsMSISE.drag = 2;
settingsMSISE.thirdbody = 1;
[~,x_obj_ukf_rom_pred_msise]=ode113(@(t,x) computeAcceleration(t,x,etf_ukf,settingsMSISE,[]),[0,time_obs_2],xf_obj_ukf_rom,opts);
residuals_obj_ukf_rom_pred_msise = xx_GPS_obs(1:3,:)-x_obj_ukf_rom_pred_msise(2:end,1:3)';

settingsJB2000.drag=3;
settingsJB2000.thirdbody = 1;
[~,x_obj_ukf_rom_pred_jb2008]=ode113(@(t,x) computeAcceleration(t,x,etf_ukf,settingsJB2000,[]),[0,time_obs_2],xf_obj_ukf_rom,opts);
residuals_obj_ukf_rom_pred_jb2008 = xx_GPS_obs(1:3,:)-x_obj_ukf_rom_pred_jb2008(2:end,1:3)';


%% ESTIMATE ORBIT WITH UKF AND PREDICT USING NRLMSISE
% Initial state guess of object
x0g_obj = x0g(svs*(objNum-1)+1:svs*(objNum-1)+6); % pos & vel
x0g_obj(end+1) = x0g(svs*objNum)/1000; % Change units BC

% Measurements and covariances for objects
Meas_obj = Meas(6*(objNum-1)+1:6*(objNum-1)+6,:);
P_obj = P(svs*(objNum-1)+1:svs*objNum,svs*(objNum-1)+1:svs*objNum);
P_obj(end,end) = P_obj(end,end)/1e6; % Change units BC cov
RM_obj = RM(6*(objNum-1)+1:6*objNum,6*(objNum-1)+1:6*objNum);
Q_obj = Q(svs*(objNum-1)+1:svs*objNum,svs*(objNum-1)+1:svs*objNum);

%% UKF MSISE
% % Set state propagation and measurement functions
% stateFnc_MEE_MSISE = @(xx,t0,tf) propagateStateMEE_FullGravDrag_NRLMSISE(xx,t0,tf,et0);
% Set state propagation and measurement functions
settingsMSISEukf.drag = 2;
settingsMSISEukf.thirdbody = 0;
stateFnc_MEE_MSISE = @(xx,t0,tf) propagateStateMEE(xx,t0,tf,et0,settingsMSISEukf,[]);
measurementFcn = @(xx) fullmee2mee(xx,1,svs);

% Run Unscented Kalman filter
[x_obj_ukf_msise,Pv_obj_ukf_msise] = UKF(x0g_obj,Meas_obj,time,stateFnc_MEE_MSISE,measurementFcn,P_obj,RM_obj,Q_obj);
[rf_obj_ukf_msise,vf_obj_ukf_msise] = ep2pv(x_obj_ukf_msise(1:6,end),GM_kms);
BC_obj_ukf_msise = x_obj_ukf_msise(end);
xf_obj_ukf_msise = [rf_obj_ukf_msise; vf_obj_ukf_msise; BC_obj_ukf_msise];

% [xf_obj_ukf_msise_pred_msise] = propagateStatePosVel_FullGravDrag_NRLMSISE(xf_obj_ukf_msise,[0,time_obs_2],etf_ukf,includeSRPMoonSunPert);
% residuals_obj_ukf_msise_pred_msise = xx_GPS_obs(1:3,:)-xf_obj_ukf_msise_pred_msise(1:3,2:end);
[~,xf_obj_ukf_msise_pred_msise]=ode113(@(t,x) computeAcceleration(t,x,etf_ukf,settingsMSISE,[]),[0,time_obs_2],xf_obj_ukf_msise,opts);
residuals_obj_ukf_msise_pred_msise = xx_GPS_obs(1:3,:)-xf_obj_ukf_msise_pred_msise(2:end,1:3)';

%% UKF JB2008
% Set state propagation and measurement functions
settingsJB2008ukf.drag = 3;
settingsJB2008ukf.thirdbody = 0;
stateFnc_MEE_JB2008 = @(xx,t0,tf) propagateStateMEE(xx,t0,tf,et0,settingsJB2008ukf,[]);
measurementFcn = @(xx) fullmee2mee(xx,1,svs);

% Run Unscented Kalman filter
[x_obj_ukf_jb2008,Pv_obj_ukf_jb2008] = UKF(x0g_obj,Meas_obj,time,stateFnc_MEE_JB2008,measurementFcn,P_obj,RM_obj,Q_obj);
[rf_obj_ukf_jb2008,vf_obj_ukf_jb2008] = ep2pv(x_obj_ukf_jb2008(1:6,end),GM_kms);
BC_obj_ukf_jb2008 = x_obj_ukf_jb2008(end);
xf_obj_ukf_jb2008 = [rf_obj_ukf_jb2008; vf_obj_ukf_jb2008; BC_obj_ukf_jb2008];

% [xf_obj_ukf_msise_pred_msise2] = propagateStatePosVel_FullGravDrag_NRLMSISE(xf_obj_ukf_jb2008,[0,time_obs_2],etf_ukf,includeSRPMoonSunPert);
% residuals_obj_ukf_msise_pred_msise2 = xx_GPS_obs(1:3,:)-xf_obj_ukf_msise_pred_msise2(1:3,2:end);
[~,x_obj_ukf_jb2008_pred_jb2008]=ode113(@(t,x) computeAcceleration(t,x,etf_ukf,settingsJB2000,[]),[0,time_obs_2],xf_obj_ukf_jb2008,opts);
residuals_obj_ukf_jb2008_pred_jb2008 = xx_GPS_obs(1:3,:)-x_obj_ukf_jb2008_pred_jb2008(2:end,1:3)';

%%

orbitPredPlot = figure; hold on;
plot(time_obs_2/3600,sqrt(sum(reshape(residuals_obj_ukf_rom_pred_rom,3,[]).^2,1)));
plot(time_obs_2/3600,sqrt(sum(reshape(residuals_obj_ukf_rom_pred_msise,3,[]).^2,1)));
plot(time_obs_2/3600,sqrt(sum(reshape(residuals_obj_ukf_rom_pred_jb2008,3,[]).^2,1)));
plot(time_obs_2/3600,sqrt(sum(reshape(residuals_obj_ukf_msise_pred_msise,3,[]).^2,1)));
plot(time_obs_2/3600,sqrt(sum(reshape(residuals_obj_ukf_jb2008_pred_jb2008,3,[]).^2,1)));
xlabel('Time [hrs]');
ylabel('Position error [km]');
title('Position error w.r.t. GPS data');
legend('Initial state from density-UKF, prediction using ROM','Initial state from density-UKF, prediction using NRLMSISE','Initial state from density-UKF, prediction using JB2008','Initial state from orbit fit with NRLMSISE, prediction using NRLMSISE','Initial state from orbit fit with JB2008, prediction using JB2008','Location','NorthWest');
set(gca,'FontSize',14);
savefig(orbitPredPlot,[filenameBase testCaseName nowTimeStr '_' planetID '_orbitPred_ROMvMSISE.fig']);

%%
rho_rom_pred = zeros(1,length(time_obs_2));
rho_msise_pred = zeros(1,length(time_obs_2));
rho_jb2008_pred = zeros(1,length(time_obs_2));
for j=1:length(rho_rom_pred)
    rho_rom_pred(j) = getDensityROM(x_obj_ukf_rom_pred_rom(j+1,1:3),time_obs_2(j)+et0_gps,romStateTime,r_pred,F_U_pred,M_U);
end
for j=1:length(rho_msise_pred)
    rho_msise_pred(j) = getDensityNRLMSISE(xf_obj_ukf_msise_pred_msise(j+1,1:3),time_obs_2(j)+et0_gps);
end
for j=1:length(rho_jb2008_pred)
    rho_jb2008_pred(j) = getDensityJB2008(x_obj_ukf_jb2008_pred_jb2008(j+1,1:3),time_obs_2(j)+et0_gps);
end
OrbitPredDensityPlot = figure;
plot(time_obs_2/3600,rho_rom_pred*1e-9,'LineWidth',1); hold on;
plot(time_obs_2/3600,rho_msise_pred*1e-9,'LineWidth',1);
plot(time_obs_2/3600,rho_jb2008_pred*1e-9,'LineWidth',1);
xlabel('Time [hr]'); ylabel('Density [kg/m^3]'); set(gca,'FontSize',14);
legend('ROM','NRLMSISE','JB2008');

%%

save([filenameBase 'workspace_' testCaseName nowTimeStr '_' planetID '_comparePredWithGPS'], ...
            'gravmodeldegree','planetID','GPSobsWindowDays','xx_GPS_obs','et_gps_obs','jdate_gps_obs', ...
            'time_obs_2','etf_ukf','et0_gps','settingsROM','settingsMSISE','settingsMSISEukf','settingsJB2000','settingsJB2008ukf', ...
            'x_obj_ukf_rom_pred_rom','xf_obj_ukf_jb2008','xf_obj_ukf_msise','x0g_obj', 'residuals_obj_ukf_rom_pred_rom', ...
            'xf_obj_ukf_msise_pred_msise','x_obj_ukf_jb2008_pred_jb2008', ...
            'residuals_obj_ukf_rom_pred_msise','residuals_obj_ukf_rom_pred_jb2008','residuals_obj_ukf_msise_pred_msise','residuals_obj_ukf_jb2008_pred_jb2008');