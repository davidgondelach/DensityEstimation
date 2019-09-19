%JSK

% clear;clc;close all

% TA = load('/Users/piyushmehta/Desktop/New_DMD_Model/DMDc_TGCM_1997_2008.mat');
% TI = load(['/Volumes/U0fMinn/TIE_GCM_Data_Full_State/MATFILES/' num2str(yr) '_TIEGCM.mat']);
clearvars;

addpath( 'AstroFunctions' );
addpath( 'DEPRICATED' );
% HEOreentryPath = fullfile('/Users','davidgondelach','Documents','PhD','HEO-reentry');
spicePath = fullfile('/Users','davidgondelach','Documents','mice');
addpath( fullfile(spicePath,'src','mice') );
addpath( fullfile(spicePath,'lib') );


%% load standard kernels and reference frames
% Clear cspice memory
cspice_kclear;
% Load SPK, PCK, LSK kernels
kernelpath  = fullfile('Data','kernel.txt');
cspice_furnsh( kernelpath );

%% load gravity model
global gravconst GM
gravconst   = 6.67259e-20; % [km^3/kg/s^2]
gravitymodel     = 'EGM2008';
gravmodeldegree  = 20;
loadGravityModel( gravitymodel, gravmodeldegree );
GM_kms = GM*1e-9;

% Constants
% mu = 398600.4415; Re = 6378.1363; J2 = 1.082626925638815*1e-3;

%% Load Earth orientation parameters (needed for TEME to ECI transformation)
global EOPMat
EOPpath = fullfile(pwd,'Data','EOP-All.txt');
[ EOPMat ] = inputEOP_Celestrak( EOPpath );

%% Load TLEs
% Setup the SGP4 propagator.
global tumin mu radiusearthkm xke j2 j3 j4 j3oj2 opsmode whichconst
opsmode = 'i';
whichconst = 72;
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );
xpdotp   =  1440.0 / (2.0*pi); % Conversion factor between SGP4 and TLE mean motion units [rev/day]/[rad/min]

%% Load atmosphere data
yr = 2005;
TA = load('DMDc_1997_2008_Var.mat');
TI = load([ num2str(yr) '_TIEGCM.mat']);

%% Generate test objects
noo = 8;
x0gRealCase = load('SimulatedTLE_x0g.mat');

rng(11);

%% Generate Orbits
%     clearvars -except TA TI yr
%     clc;close all;

nofDays = 20;
% Time Interval for measurements
% dt=300;tf=5*24*60*60;time=[0:dt:tf]';m=length(time);
dt=3600;tf=nofDays*24*60*60;time=[0:dt:tf]';m=length(time);

% Converting the dynamic and input matrices from discrete to continuous time
rt = 20;r = 10;
[PhiC,Uh,q] = C2D(TA,rt,r); % q is number of inputs

% yr=2005;
mth=7;dy=10;hr=0;mn=0;sc=0;
doy = day(datetime(yr,mth,dy),'dayofyear');

% jd = juliandate(yr,mth,dy,hr,mn,sc);
jd0 = juliandate(yr,mth,dy,0,0,0);
jdf = juliandate(yr,mth,dy+nofDays,0,0,0);

et0  = cspice_str2et(strcat([num2str(jed2date(jd0),'%d %d %d %d %d %.10f') 'UTC']));

dtInput = 3600;
[Inp1,Inp2,z0] = Comp_Inputs_Var(hr,mn,sc,doy,Uh,tf,TI,dt,q,PhiC,TA,r);

% addpath('sensor_sim_020917/')

% Setup of ROM Modal Interpolation
sltm=linspace(0,24,24);
latm=linspace(-90,90,20);
altm=linspace(100,450,16);
[SLTm,LATm,ALTm]=ndgrid(sltm,latm,altm);

F_U{r} = [];
for i = 1:r
    Uhr = reshape(Uh(:,i),24,20,16); % i-th left singular vector on grid
    F_U{i} = griddedInterpolant(SLTm,LATm,ALTm,Uhr,'linear','linear'); % Create interpolant of Uhr
end
Mr = reshape(TA.DenS_Mean,24,20,16);
M_U = griddedInterpolant(SLTm,LATm,ALTm,Mr,'linear','linear');

% True initial state
svs = 7;    %Size of state vector for each object [3xpos,3xvel,1xBC]
x0trueMEEBC = x0gRealCase.x0g(1:svs*noo);
x0_pv = x0trueMEEBC;
for k = 1:noo
    for j=1:size(x0trueMEEBC,2)
        [pos,vel] = ep2pv(x0trueMEEBC((k-1)*svs+1:(k-1)*svs+6,j),GM_kms);
        x0_pv((k-1)*svs+1:(k-1)*svs+3,j) = pos;
        x0_pv((k-1)*svs+4:(k-1)*svs+6,j) = vel;
    end
end

x0 = zeros(svs*noo+r,1);
x0(1:svs*noo) = x0_pv;
x0(end-r+1:end) = z0;

% True trajectories (in Cartesian coordinates)
AC = PhiC(1:r,1:r)/3600; BC = PhiC(1:r,r+1:end)/3600;
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'OutputFcn',@odewbar);
[~,temp]=ode113(@(t,x) Propagation11_ODE_Var_FullGrav(t,x,AC,BC,Inp2,r,noo,svs,F_U,M_U,yr,time,et0),time,x0,opts);
XP = temp';

% % Weed out the diverged trajectories
% wD = zeros(noo,1);
% k=1;
% for i = 1:noo
%     wD(i,1) = sqrt(XP(svs*(i-1)+1,end).^2+XP(svs*(i-1)+2,end).^2+XP(svs*(i-1)+3,end).^2);
%     if wD(i,1) < 1e10 && (wD(i,1)-6378) > 120
%         XP2(svs*(k-1)+1:svs*(k-1)+svs,:) = XP(svs*(i-1)+1:svs*(i-1)+svs,:);
% %         BC1(1,k) = b_star(1,i);
%         k=k+1;
%     end
% end
% XP2((k-1)*svs+1:(k-1)*svs+r,:) = XP(end-r+1:end,:); % TODO: check index
% clear XP b_star;
% XP = XP2; %b_star = BC1;
% noo = size(b_star,2);

% Longitude, Latitude and Height of Satellite
long = zeros(m,noo);lat = zeros(m,noo);height = zeros(m,noo);
for i = 1:noo
    [long(:,i),lat(:,i),height(:,i)]=gc2gd(XP(svs*(i-1)+1:svs*(i-1)+3,:)',yr,mth,dy,hr,mn,sc,dt,tf,1);
end

% True densities along trajectories
rho = zeros(size(long,1),noo);
slt = zeros(size(long,1),size(long,2));
for i = 1:noo
    long(long(:,i)>180,i) = long(long(:,i)>180,i) - 360;
    slt(:,i) = Inp2(3,:)'+long(:,i)/15;
    slt(slt(:,i)>24,i) = slt(slt(:,i)>24,i)-24;slt(slt(:,i)<0,i) = slt(slt(:,i)<0,i)+24;
    UhI = zeros(size(long,1),r);
    for j = 1:r
        UhI(:,j) = F_U{j}(slt(:,i),lat(:,i),height(:,i));
    end
    MI = M_U(slt(:,i),lat(:,i),height(:,i));
    rho(:,i) = 10.^(sum(UhI.*XP(end-r+1:end,:)',2)+MI);
end


%% Set uncertainties
% posErr = 1000/1000; %10/1000;
% velErr = 1./1000;

%% Generate measurements

% nop = 10;    %number of propagated orbits used for measurements
% Meas = zeros(nop*3,size(XP,2)); GPS_Error = 10/1000;
% Meas = zeros(nop,size(XP,2)); GPS_Error = 10/1000; SMA_Error = 0.0;
% RM = SMA_Error^2*eye(nop);
% RM = diag(repmat([0.02 2e-6 2e-6 7e-7 7e-7 1.3e-6]',nop,1));
% posErr = 500/1000; %10/1000;
% velErr = 0.01/1000;

% Meas = state2posNoise(XP,nop,svs,posErr);
% RM = posErr^2*eye(nop*3);
% Meas = state2posVelNoise(XP,nop,svs,posErr,velErr);
% RM = diag(repmat([posErr^2 posErr^2 posErr^2 velErr^2 velErr^2 velErr^2]',nop,1));
rng(5);
% Meas = state2meeNoise(XP,nop,svs,posErr,velErr);
mee = zeros(noo*6,size(XP,2));
Meas = zeros(noo*6,size(XP,2));
for k = 1:noo
    for j=1:size(XP,2)
        pos = XP(svs*(k-1)+1:svs*(k-1)+3,j);
        vel = XP(svs*(k-1)+4:svs*(k-1)+6,j);
        mee(6*(k-1)+1:6*k,j) = pv2ep(pos,vel,GM_kms)';
    end
end

rng(6);
% pErr = 0.014; fErr = 1.6e-5; gErr = 1.6e-5; hErr = 1.1e-5; kErr = 1.1e-5; LErr = 0.00011; % Test case 1
pErr = 0.045; fErr = 2.0e-5; gErr = 2.0e-5; hErr = 2.0e-5; kErr = 2.0e-5; LErr = 0.000125; % Test case 2: realistic std: 0.04522508	2.36E-05	1.90E-05	2.17E-05	1.92E-05	0.000127646
Meas(1:6:end,:) = mee(1:6:end,:) + pErr * randn(noo,size(XP,2));
Meas(2:6:end,:) = mee(2:6:end,:) + fErr * randn(noo,size(XP,2));
Meas(3:6:end,:) = mee(3:6:end,:) + gErr * randn(noo,size(XP,2));
Meas(4:6:end,:) = mee(4:6:end,:) + hErr * randn(noo,size(XP,2));
Meas(5:6:end,:) = mee(5:6:end,:) + kErr * randn(noo,size(XP,2));
Meas(6:6:end,:) = mee(6:6:end,:) + LErr * randn(noo,size(XP,2));

measErrors = Meas - mee;
measVar = var(reshape(measErrors,6,[]),0,2);
RM = diag(repmat(measVar,noo,1));
% RM = diag(var(measErrors,0,2));

for k = 1:noo
    for j=1:size(Meas,2)
        [posMeas,velMeas] = ep2pv(Meas((k-1)*6+1:(k-1)*6+6,j),GM_kms);
        posTrue = XP(svs*(k-1)+1:svs*(k-1)+3,j);
        velTrue = XP(svs*(k-1)+4:svs*(k-1)+6,j);
        posError(k,j) = sqrt(sum( (posMeas-posTrue) .^2,1));
        velError(k,j) = sqrt(sum( (velMeas-velTrue) .^2,1));
    end
end
posError2 = posError.^2;
posError2 = reshape(posError2,1,[]);
posStd = sqrt(mean(posError2));
velError2 = velError.^2;
velError2 = reshape(velError2,1,[]);
velStd = sqrt(mean(velError2));


%% Generate initial state guess

nop = noo;    %number of propagated orbits used for measurements
x0g = zeros(svs*nop+r,1);   % initial state guess

rng(2);
% x0state = state2meeNoise(XP(:,1),nop,svs,posErr,velErr);
x0state = Meas(:,1);

rng(3);
for i = 1:nop
    x0g(svs*(i-1)+1:svs*(i-1)+6,1) = x0state(6*(i-1)+1:6*(i-1)+6,1);
    x0g(svs*i) = XP(svs*i,1) - 0.20*(-1)^floor(randn)*XP(svs*i,1);
end

% GPS_Error = 10/1000;
% for i = 1:nop
%     x0f(svs*(i-1)+1:svs*i,1) = XP(svs*(i-1)+1:svs*i,1);
%     UC(3*(i-1)+1:3*i,:) = GPS_Error * randn(3,size(XP,2));
%     x0f(svs*(i-1)+1:svs*(i-1)+3,1) = XP(svs*(i-1)+1:svs*(i-1)+3,1) + UC(3*(i-1)+1:3*i,1);
% end

% % Compute the Initial Atmosphere State from MSIS and MSIS densities along
% % simulated orbits
% addpath('nrlmsise_matlab/')
% UT = hr*3600+mn*60+sc;
% [Rho_M,rho_M] = MSIS(yr,doy,UT,SLTm,LATm,ALTm,long,lat,height,slt,Inp2);
% z0_M = Uh'*(log10(Rho_M)-TA.DenS_Mean);


    % Compute the Initial Atmosphere State from MSIS
    UT = hr*3600+mn*60+sc;
    
    % Space weather data
    SWpath = fullfile('Data','SW-All.txt');
    [ SWmatDaily, SWmatMonthlyPred ] = inputSWnrlmsise( SWpath );
    [ f107A, f107, ap ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jd0 );
    sltx = reshape(SLTm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    latx = reshape(LATm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    altx = reshape(ALTm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    
    % Density
    Den = zeros(numel(sltx),1);
    for i = 1:numel(sltx)
        lon = 15*(sltx(i)-UT/3600);
        
        [output] = nrlmsise(altx(i),latx(i),lon,yr,doy,UT,sltx(i),...
            f107A,f107,ap);
        Den(i,1) = output.d(6);
    end
    
    % Initial guess for ROM
    z0_M = Uh'*(log10(Den)-TA.DenS_Mean); % NRLMSISE-00 initialization
    

% z0_I = Uh'*TI.DenS(:,1);
rng(4);
x0g(end-r+1:end,1) = z0_M;% + normrnd(0,0.25,[r,1])*(-1)^floor(randn).*z0_M;

clear TA TI



%% Optimization for estimating the process noise and Initial Covariance
close all
Pv = zeros(svs*nop+r,1); % state covariance
Qv = zeros(svs*nop+r,1); % process coveriance
for i = 1:nop
%     Pv(svs*(i-1)+1:svs*(i-1)+3) = 1e-4;
%     Pv(svs*(i-1)+4:svs*(i-1)+6) = 1e-5;
    Pv(svs*(i-1)+1) = measVar(1);
    Pv(svs*(i-1)+2) = measVar(2);
    Pv(svs*(i-1)+3) = measVar(3);
    Pv(svs*(i-1)+4) = measVar(4);
    Pv(svs*(i-1)+5) = measVar(5);
    Pv(svs*(i-1)+6) = measVar(6);
%     Pv(svs*(i-1)+7) = 7e0;
    Pv(svs*(i-1)+7) = (XP(svs*i,1) * 0.2)^2; % 20% BC error (std)
    % Assuming 10 m position and 0.01 m/s velocity (1-std) error after 1h
    % propagation (simulate 1000 true MEE and noisy MEE and compute cov)
    % cov(p,f,g,h,k,L) = 
    % [0.000720344406143554;6.89889751950605e-12;6.66859838015428e-12;1.52574451666460e-12;1.38725517061063e-12;3.30710442740959e-12]
%     Qv(svs*(i-1)+1:svs*(i-1)+3) = 1e-14;
%     Qv(svs*(i-1)+4:svs*(i-1)+6) = 1e-14;
%     Qv(svs*(i-1)+7) = 1e-14;
    Qv(svs*(i-1)+1) = 1.5e-8; %7.2e-4;
    Qv(svs*(i-1)+2) = 2e-14; %6.9e-12;
    Qv(svs*(i-1)+3) = 2e-14; %6.9e-12;
    Qv(svs*(i-1)+4) = 1e-14; %1.5e-12;
    Qv(svs*(i-1)+5) = 1e-14; %1.5e-12;
    Qv(svs*(i-1)+6) = 1e-12; %3.3e-12;
    Qv(svs*(i-1)+7) = 1e-12; % 1-std error is 1e-6 per 1 hour
%     [1.63233439347216e-08;2.41366395075001e-14;1.90396329360770e-14;9.08604764627490e-15;1.04736563945798e-14;1.12575443917672e-12]
end
Pv(end-r+1:end) = (5e0)*ones(r,1);
Qv(end-r+1:end) = 1e-14*ones(r,1);
Pv(end-r+1) = 2e1;
P = diag(Pv);
Q = diag(Qv);

X_est=x0g;
% % Set BC errors
% for i = 1:nop
%     X_est(svs*i,1) = X_est(svs*i,1) - 0.25*(-1)^floor(randn)*X_est(svs*i,1);
% end
% size_X = numel(x0f);Xp = zeros(size_X,2*size_X+1);

% Set state propagation and measurement functions
stateFnc = @(xx,t0,tf) propagateStateMEE_FullGravDrag(xx,t0,tf,AC,BC,Inp2,r,nop,svs,F_U,M_U,yr,time,et0);
% stateFnc = @(xx,t0,tf) propagateStateMEE(xx,t0,tf,AC,BC,Inp2,r,nop,svs,F_U,M_U,yr,time);
% measurementFcn = @(xx) state2measurement(xx,nop,svs);
% measurementFcn = @(xx) state2mee(xx,nop,svs);
% measurementFcn = @(xx) state2pos(xx,nop,svs);
measurementFcn = @(xx) fullmee2mee(xx,nop,svs);

% Run Unscented Kalman filter
[X_est,Pv] = UKF(X_est,Meas,time,stateFnc,measurementFcn,P,RM,Q);

%%
nowTimeStr = datestr(now,'yymmddHHMMSS');
filenameBase = '/Users/davidgondelach/Google Drive/PostDoc/sqrt-ukf-rom_';
testCaseName = 'simulatedTLEobs_8obj_';

save([filenameBase 'workspace_' testCaseName nowTimeStr]);
%%
Pv(:,1) = diag(P);

set(0,'defaultAxesFontSize',14)
set(0, 'DefaultLineLineWidth', 1);

X_est_pv = X_est;
X_est_coe = X_est;
for k = 1:nop
    for j=1:size(X_est,2)
        [pos,vel] = ep2pv(X_est((k-1)*svs+1:(k-1)*svs+6,j),GM_kms);
        X_est_pv((k-1)*svs+1:(k-1)*svs+3,j) = pos;
        X_est_pv((k-1)*svs+4:(k-1)*svs+6,j) = vel;
        X_est_coe((k-1)*svs+1:(k-1)*svs+6,j) = pv2po(pos,vel,GM_kms);
    end
end

long_est = zeros(m,nop);lat_est = zeros(m,nop);height_est = zeros(m,nop);
b_est = zeros(1,nop);
for i = 1:nop
    % Longitude, Latitude and Height of Satellite
    [long_est(:,i),lat_est(:,i),height_est(:,i)]=gc2gd(X_est_pv(svs*(i-1)+1:svs*(i-1)+3,:)',yr,mth,dy,hr,mn,sc,dt,tf,1);
    b_est(1,i) = X_est(svs*i,end);
end

rho_est = zeros(size(long_est,1),nop);
slt_est = zeros(size(long_est,1),size(long_est,2));
for i = 1:nop
    long_est(long_est(:,i)>180,i) = long_est(long_est(:,i)>180,i) - 360;
    slt_est(:,i) = Inp2(3,:)'+long_est(:,i)/15;
    slt_est(slt_est(:,i)>24,i) = slt_est(slt_est(:,i)>24,i)-24;slt_est(slt_est(:,i)<0,i) = slt_est(slt_est(:,i)<0,i)+24;
    UhI = zeros(size(long_est,1),r);
    for j = 1:r
        UhI(:,j) = F_U{j}(slt_est(:,i),lat_est(:,i),height_est(:,i));
    end
    MI_est = M_U(slt_est(:,i),lat_est(:,i),height_est(:,i));
    rho_est(:,i) = 10.^(sum(UhI.*X_est(end-r+1:end,:)',2)+MI_est);
    %rho(:,i) = 10.^(sum(UhI.*XP(end-r+1:end,:)',2)+MI);
end

densPlot = figure;
for i = 1:nop
    subplot(5,4,2*i-1)
    plot(time/3600,rho(:,i));hold on;plot(time/3600,rho_est(:,i))
    title(sprintf('Orbit %.0f',i));
    xlabel('Time, hrs');ylabel('Density [kg/m^3]');
    legend('True','Estimated','Location','SouthEast','Orientation','Horizontal');
    subplot(5,4,2*i)
    plot(time/3600,abs(rho(:,i)-rho_est(:,i))./rho(:,i).*100);
    title(sprintf('Orbit %.0f',i));
    xlabel('Time, hrs');ylabel('Density error [%]');
end
% savefig(densPlot,'/Users/davidgondelach/Google Drive/PostDoc/sqrt-ukf-rom_estVSmeasDens_10obj_5d_GPSmeeObs1h_500mObsErr.fig');
% savefig(densPlot,[filenameBase 'dens_' testCaseName nowTimeStr '.fig']);

ROMplot = figure;
for i = 1:r
    subplot(5,4,2*i-1)
    plot(time/3600,XP(end-r+i,:));hold on;plot(time/3600,X_est(end-r+i,:))
    xlabel('Time, hrs');ylabel('z');
    title(sprintf('Mode %.0f',i));
    legend('True','Estimated','Location','SouthEast','Orientation','Horizontal');
    subplot(5,4,2*i)
    plot(time/3600,abs(XP(end-r+i,:)-X_est(end-r+i,:))); hold on;
    plot(time/3600,3*Pv(end-r+i,:).^0.5,'--k');
    plot(time/3600,-3*Pv(end-r+i,:).^0.5,'--k');
    legend('Error','3\sigma','Location','SouthEast','Orientation','Horizontal');
    title(sprintf('Mode %.0f',i)); grid on;
    xlabel('Time, hrs');ylabel('Error');
end
% savefig(ROMplot,'/Users/davidgondelach/Google Drive/PostDoc/sqrt-ukf-rom_estVSmeasROMmodes_10obj_5d_GPSmeeObs1h_500mObsErr.fig');
% savefig(ROMplot,[filenameBase 'ROMmodes_' testCaseName nowTimeStr '.fig']);


BCplot = figure;
for i = 1:nop
    subplot(5,4,2*i-1)
    plot(time/3600,XP(svs*i,:)/1000);hold on;plot(time/3600,X_est(svs*i,:)/1000)
    xlabel('Time, hrs');ylabel('BC [m^2/kg]');
    title(sprintf('Orbit %.0f',i));
    legend('True','Estimated','Location','SouthEast','Orientation','Horizontal');
    subplot(5,4,2*i)
    plot(time/3600,abs(XP(svs*i,:)-X_est(svs*i,:))./XP(svs*i,:)*100); hold on;
    plot(time/3600,3*Pv(svs*i,:).^0.5./XP(svs*i,:)*100,'--k');
    legend('Error','3\sigma','Location','SouthEast','Orientation','Horizontal');
    xlabel('Time, hrs'); ylabel('Error BC [%]');
    title(sprintf('Orbit %.0f',i));
end
% savefig(BCplot,'/Users/davidgondelach/Google Drive/PostDoc/sqrt-ukf-rom_estVSmeasBC_10obj_5d_GPSmeeObs1h_500mObsErr.fig');
% savefig(BCplot,[filenameBase 'BC_' testCaseName nowTimeStr '.fig']);

posPlot = figure;
for i = 1:nop
    posErrors = sqrt(sum( (XP(svs*(i-1)+1:svs*(i-1)+3,:)-X_est_pv(svs*(i-1)+1:svs*(i-1)+3,:)) .^2,1));
    subplot(5,2,i)
    plot(time/3600,posErrors); hold on;
    xlabel('Time, hrs');ylabel('Position error [km]');
    title(sprintf('Orbit %.0f, mean= %.2f',i,mean(posErrors)));
end
% savefig(posPlot,'/Users/davidgondelach/Google Drive/PostDoc/sqrt-ukf-rom_posErr_10obj_5d_GPSmeeObs1h_500mObsErr.fig');
% savefig(posPlot,[filenameBase 'posErr_' testCaseName nowTimeStr '.fig']);

meePlot = figure;
for i = 1:nop
    meeDiff = Meas((i-1)*6+1:i*6,:) - X_est((i-1)*svs+1:(i-1)*svs+6,:);
    for j=1:6
        subplot(2,3,j); hold on;
        plot(time/3600,meeDiff(j,:));
    end
    meeDiffMean(i,:) = mean(meeDiff,2);
    meeDiffStd(i,:) = std(meeDiff,0,2);
end
subplot(2,3,1); xlabel('Time [hours]'); ylabel('p [km]');
subplot(2,3,2); xlabel('Time [hours]'); ylabel('f [-]');
subplot(2,3,3); xlabel('Time [hours]'); ylabel('g [-]');
subplot(2,3,4); xlabel('Time [hours]'); ylabel('h [-]');
subplot(2,3,5); xlabel('Time [hours]'); ylabel('k [-]');
subplot(2,3,6); xlabel('Time [hours]'); ylabel('L [rad]');
% savefig(meePlot,[filenameBase 'MEEestVmeas_' testCaseName nowTimeStr '.fig']);

figure;
subplot(2,3,1); plot(time/86400,X_est_coe(1:svs:end-r,:)); hold on;  xlabel('t'); ylabel('a [km]');
subplot(2,3,2); plot(time/86400,X_est_coe(2:svs:end-r,:)); hold on;  xlabel('t'); ylabel('e');
subplot(2,3,3); plot(time/86400,rad2deg(X_est_coe(3:svs:end-r,:))); hold on;  xlabel('t'); ylabel('i');
subplot(2,3,4); plot(time/86400,rad2deg(X_est_coe(4:svs:end-r,:))); hold on;  xlabel('t'); ylabel('\Omega');
subplot(2,3,5); plot(time/86400,rad2deg(X_est_coe(5:svs:end-r,:))); hold on; xlabel('t'); ylabel('\omega');
subplot(2,3,6); plot(time/86400,rad2deg(X_est_coe(6:svs:end-r,:))); hold on;  xlabel('t'); ylabel('M [deg]');

XP_coe = zeros(6*nop,size(XP,2));
for k = 1:nop
for j=1:size(XP,2)
[pos] = XP((k-1)*svs+1:(k-1)*svs+3,j);
[vel] = XP((k-1)*svs+4:(k-1)*svs+6,j);
XP_coe((k-1)*6+1:(k-1)*6+6,j) = pv2po(pos,vel,GM_kms);
end
end
figure;
subplot(2,3,1); plot(time/86400,XP_coe(1:6:end-r,:)); hold on;  xlabel('t'); ylabel('a [km]');
subplot(2,3,2); plot(time/86400,XP_coe(2:6:end-r,:)); hold on;  xlabel('t'); ylabel('e');
subplot(2,3,3); plot(time/86400,rad2deg(XP_coe(3:6:end-r,:))); hold on;  xlabel('t'); ylabel('i');
subplot(2,3,4); plot(time/86400,rad2deg(XP_coe(4:6:end-r,:))); hold on;  xlabel('t'); ylabel('\Omega');
subplot(2,3,5); plot(time/86400,rad2deg(XP_coe(5:6:end-r,:))); hold on; xlabel('t'); ylabel('\omega');
subplot(2,3,6); plot(time/86400,rad2deg(XP_coe(6:6:end-r,:))); hold on;  xlabel('t'); ylabel('M [deg]');

%%
addpath('/Users/davidgondelach/Documents/PhD/Doctorado/MATLAB Scripts/altmany-export_fig');

ROMplot = figure;
for i = 1:4
    subplot(2,4,2*i-1)
    plot(time/86400,XP(end-r+i,:));hold on;
    plot(time/86400,X_est(end-r+i,:))
    xlabel('Time [days]');ylabel('z [-]');xticks([0:5:20]); grid on;
    xlim([0 15])
    title(sprintf('Mode %.0f',i));
    legend('True','Estimated','Location','SouthEast','Orientation','Vertical');
    subplot(2,4,2*i)
    plot(time/86400,abs(XP(end-r+i,:)-X_est(end-r+i,:))); hold on;
    plot(time/86400,3*Pv(end-r+i,:).^0.5,'--k');
    plot(time/86400,-3*Pv(end-r+i,:).^0.5,'--k');
    legend('Error','3\sigma','Location','SouthEast','Orientation','Vertical');
    xlim([0 15])
    title(sprintf('Mode %.0f',i)); grid on;xticks([0:5:20]);
    xlabel('Time [days]');ylabel('Error [-]');
end
% export_fig('RomModesAndErrWithCov.pdf', '-transparent', '-nocrop');

%%
BCplot = figure;
for i = 1:nop
    plot(time/86400,abs(XP(svs*i,:)-X_est(svs*i,:))./XP(svs*i,:)*100); hold on; grid on;
end
xlabel('Time [days]'); ylabel('BC error [%]');
legend('Object 1','Object 2','Object 3','Object 4','Object 5',...
    'Object 6','Object 7','Object 8',...
    'Location','NorthEast','Orientation','Vertical');

%%
densPlot = figure;
for i = 1:nop
    plot(time/86400,abs(rho(:,i)-rho_est(:,i))./rho(:,i).*100); hold on; grid on;
end
xlabel('Time [days]'); ylabel('Density error [%]');
legend('Object 1','Object 2','Object 3','Object 4','Object 5',...
    'Object 6','Object 7','Object 8',...
    'Location','NorthEast','Orientation','Vertical');

%%
figure;
for i = 1:nop
    subplot(1,2,1);
    plot(time/86400,abs(rho(:,i)-rho_est(:,i))./rho(:,i).*100); hold on; grid on;
    subplot(1,2,2);
    plot(time/86400,abs(XP(svs*i,:)-X_est(svs*i,:))./XP(svs*i,:)*100); hold on; grid on;
end
subplot(1,2,1); xlabel('Time [days]'); ylabel('Density error [%]');
legend('Object 1','Object 2','Object 3','Object 4','Object 5',...
    'Object 6','Object 7','Object 8',...
    'Location','NorthEast','Orientation','Vertical');
subplot(1,2,2); xlabel('Time [days]'); ylabel('BC error [%]');
legend('Object 1','Object 2','Object 3','Object 4','Object 5',...
    'Object 6','Object 7','Object 8',...
    'Location','NorthEast','Orientation','Vertical');

%%
% sma_est = zeros(nop,size(X_est,2));
% orbE_est = zeros(nop,size(X_est,2));
% mu = 398600.4415; Re = 6378.1363; J2 = 1.082626925638815*1e-3;
% for k = 1:nop
%     coe_est = zeros(size(X_est,2),6);
%     for j=1:size(X_est,2)
%         coe_est(j,:) = pv2po(X_est(svs*(k-1)+1:svs*(k-1)+3,j),X_est(svs*(k-1)+4:svs*(k-1)+6,j),mu);
%         mee_est(j,:) = pv2ep(X_est(svs*(k-1)+1:svs*(k-1)+3,j),X_est(svs*(k-1)+4:svs*(k-1)+6,j),mu);
% %         sma_est(k,j) = coe_est(j,1);
% %         rr_est = norm(X_est(svs*(k-1)+1:svs*(k-1)+3,j));
% %         vv_est = norm(X_est(svs*(k-1)+4:svs*(k-1)+6,j));
%         %         orbE_est(k,j) = 0.5*vv_est^2 - mu/rr_est * ( 1 - 0.5*J2*Re^2/rr_est^2*(3*sin(coe(5)+coe(6))^2*sin(coe(3))^2-1) );
% %         orbE_est(k,j) = 0.5*vv_est^2 - mu/rr_est * ( 1 - 0.5*J2*Re^2/rr_est^2*(3*X_est(svs*(k-1)+3,j)^2/rr_est^2-1) );
%         %         orbE_est(k,j) = 0.5*vv_est^2 - mu/rr_est;
% %         avgSMA(k,j) = -mu/(2*orbE_est(k,j));
%     end
% %     figure(111);
% %     subplot(5,2,k);
%     %     plot(time/3600,sma_est(k,:),'-.b',time/3600,Meas(k,:),'-.m','Linewidth',2);grid on;
% %     plot(time/3600,Meas(k,:),time/3600,orbE_est(k,:));grid on;
%     %     plot(time/3600,orbE_est(k,:),'r');grid on;
%     %     xlabel('Time, hrs');ylabel('a [km]');
% %     xlabel('Time, hrs');ylabel('Orb. energy [km^2/s^2]');
%
% %         figure;
% %         for n=1:6
% %             subplot(3,2,n); plot(time/3600,coe_est(:,n));
% %         end
%         figure;
%         for n=1:6
%             subplot(3,2,n); plot(time/3600,mee_est(:,n));
%         end
% end
% legend('True','Estimated','Location','SouthEast','Orientation','Horizontal');

% for k = 1:nop
%     figure;
%     for n=1:3
%         subplot(1,3,n);
%         plot(time/3600,X_est(svs*(k-1)+n,:)-XP(svs*(k-1)+n,:)); hold on;
%         plot(time/3600,Meas(3*(k-1)+n,:)-XP(svs*(k-1)+n,:));
% %         plot(time/3600,XP(svs*(k-1)+n,:));
%     end
% end

%
% figure;set(gcf,'color','w');
% for i = 1:nop
%     subplot(5,2,i)
%     plot(time/3600,rho(:,i),'-.b',time/3600,rho_est(:,i),'-.m','Linewidth',2);grid on;
% %     xlim([0 24]);ylim([min(min(rho(1:round(24*3600/300),i)),min(rho_est(1:round(24*3600/300),i))) ...
% max(max(rho(1:round(24*3600/300),i)),max(rho_est(1:round(24*3600/300),i)))])
%     if i > nop-2;xlabel('Time, hrs');end;if i == 5;ylabel('Density (kg/m^3)');end
%     if i == 1;legend('True','Estimated','Orientation','horizontal');end
%     title(sprintf('Orbit %.0f',i));
%     set(gca,'fontsize',14)
% end

%
% figure(30);set(gcf,'color','w');
% for i = 1:r
%     subplot(5,2,i)
%
%     yyaxis left
%     r1 = plot(time/3600,XP(end-r+i,:),'b','Linewidth',2);hold on;
%     r2 = plot(time/3600,X_est(end-r+i,:),'m','Linewidth',2);hold on;
%     r3 = plot(time/3600,X_est(end-r+i,:)+3*Pv(end-r+i,:),'-.k','Linewidth',2);hold on;
%     plot(time/3600,X_est(end-r+i,:)-3*Pv(end-r+i,:),'-.k','Linewidth',2);grid on;
%     xlim([0 120])
%     if i==1
%         ylim([-2170 -2020])     %1 and 3 satellite
% %       ylim([-2120 -2070])   % 50 satellites
%       legend([r1 r2 r3],{'True','Estimated','+/- 3\sigma'},'Orientation','horizontal','Location','best','AutoUpdate','off')
%     elseif i == 2 || i == 3
%         ylim([-20 20])  %1 satellite
% %         ylim([-10 10])  %3 and 50 satellite
%     else
%       ylim([-10 10])
%     end
%     if i == 5;ylabel('z, (~)');end
%     if i > 8; xlabel('Time, hrs'); end
%     title(sprintf('Mode %.0f',i));
%     set(gca,'fontweight','bold','fontsize',20)
%
%     yyaxis right
%     plot(time/3600,XP(end-r+i,:)-X_est(end-r+i,:),'Linewidth',2);
%     if i == 6;ylabel('Absolute Difference, (~)');end
%     set(gca,'fontweight','bold','fontsize',20)
% end

% figure(301);set(gcf,'color','w');
% for i = 1:r
%     subplot(5,2,i)
%     plot(time/3600,XP(end-r+i,:)-X_est(end-r+i,:),'b','Linewidth',2);hold on;
%     plot(time/3600,3*Pv(end-r+i,:).^0.5,'-.k','Linewidth',2);hold on;
%     plot(time/3600,-3*Pv(end-r+i,:).^0.5,'-.k','Linewidth',2);grid on;
%     ylim([-5*3*Pv(end-r+i,end).^0.5 5*3*Pv(end-r+i,end).^0.5]);
%     if i == 1
%         legend('Error','+/- 3\sigma','Orientation','horizontal','Location','best','AutoUpdate','off')
%     end
%     if i < 9; set(gca,'Xticklabel',[]) ;end
%     ylabel('Error, (~)');
%     if i > 8; xlabel('Time, hrs'); end
%     title(sprintf('Mode %.0f',i));
%     set(gca,'fontweight','bold','fontsize',20)
% end
%
% BCi = zeros(nop,1);
% for i = 1:nop
%     BCi(i,1) = X_est(svs*i);
% end
% figure(31);set(gcf,'color','w');
% bar([b_star(1:nop)' BCi b_est']);grid on;
% legend('True','Initialized','Estimated','Orientation','horizontal','Location','best');
% xlim([0 nop+1]);
% xlabel('Number of Satellite');ylabel('Ballistic Coefficient, (~)')
% set(gca,'fontweight','bold','fontsize',20)

% figure(32);set(gcf,'color','w');
% subplot(2,1,1)
% plot(time/3600,Inp2(1,:),'Linewidth',2); grid on
% ylabel('F_{10.7} index, sfu')
% set(gca,'fontweight','bold','fontsize',20)
% subplot(2,1,2)
% plot(time/3600,Inp2(2,:),'Linewidth',2); grid on
% xlabel('Time, hrs'); ylabel('K_p index, (~)')
% set(gca,'fontweight','bold','fontsize',20)


% %% Propagating Orbits with MSIS
% x0M = XP(1:svs*nop,1);x0M = [x0M;z0_M];
% x0M_ROM = x0M;
% load('Kp2Ap.mat');
%
% opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'OutputFcn',@odewbar);
%
% [~,temp]=ode45(@(t,x) Propagation11_ODE(t,x,AC,BC,Inp2,r,nop,svs,F_U,yr,time),time,x0M_ROM,opts);
% XPM_ROM = temp';
%
% [~,temp]=ode45(@(t,x) Propagation11_ODE_M(t,x,AC,BC,Inp2,r,nop,svs,K2A,yr,time),time,x0M,opts);
% XPM = temp';

% Scaling the Ballistic Coefficients
% loni = long(1,1:nop)';lati = lat(1,1:nop)'; alti = height(1,1:nop)';
% loni(loni>180) = loni(loni>180) - 360;
% slti = Inp2(3,1)+loni/15;
% slti(slti>24) = slti(slti>24)-24;slti(slti<0) = slti(slti<0)+24;
% UhI = zeros(nop,r);
% for j = 1:r
%    UhI(:,j) = F_U{j}(slti,lati,alti);
% end
% rho_T_Ini(1,:) = exp(sum(UhI'.*z0,1));
% rho_M_Ini(1,:) = exp(sum(UhI'.*z0_M,1));
% BC_Scaled = b_star(1,1:nop) .* rho_T_Ini./rho_M_Ini;

% x0MS = XP(1:svs*nop,1);x0MS = [x0MS;z0_M];
% for i = 1:nop
%     x0MS(svs*i) = BC_Scaled(1,i);
% end
% x0MS_ROM = x0MS;
%
% [~,temp]=ode45(@(t,x) Propagation11_ODE(t,x,AC,BC,Inp2,r,nop,svs,F_U,yr,time),time,x0MS_ROM,opts);
% XPMS_ROM = temp';
%
% [~,temp]=ode45(@(t,x) Propagation11_ODE_M(t,x,AC,BC,Inp2,r,nop,svs,K2A,yr,time),time,x0MS,opts);
% XPMS = temp';
%
% Pos_Er_M = zeros(nop,m); Pos_Er_M_ROM = zeros(nop,m);
% Pos_Er_MS = zeros(nop,m); Pos_Er_MS_ROM = zeros(nop,m);
% for i = 1:nop
%     Pos_Er_M_ROM(i,:) = sqrt(sum((XP(svs*(i-1)+1:svs*(i-1)+3,:) - XPM_ROM(svs*(i-1)+1:svs*(i-1)+3,:)).^2,1));
%     Pos_Er_M(i,:) = sqrt(sum((XP(svs*(i-1)+1:svs*(i-1)+3,:) - XPM(svs*(i-1)+1:svs*(i-1)+3,:)).^2,1));
%     Pos_Er_MS_ROM(i,:) = sqrt(sum((XP(svs*(i-1)+1:svs*(i-1)+3,:) - XPMS_ROM(svs*(i-1)+1:svs*(i-1)+3,:)).^2,1));
%     Pos_Er_MS(i,:) = sqrt(sum((XP(svs*(i-1)+1:svs*(i-1)+3,:) - XPMS(svs*(i-1)+1:svs*(i-1)+3,:)).^2,1));
% end
%
% colors = distinguishable_colors(nop);
% figure(33);set(gcf,'color','w');
% for i = 1:nop
%     semilogy(time/3600,Pos_Er_M_ROM(i,:),'Color',colors(i,:),'Linewidth',2);hold on;
%     legendInfo{i} = ['Orbit ' num2str(i)];
% end
% for i = 1:nop
%     semilogy(time/3600,Pos_Er_M(i,:),'-.','Color',colors(i,:),'Linewidth',2);hold on;
% end
% grid on;
% legend(legendInfo,'Location','Best');
% xlabel('Time, hrs');ylabel('log_{10}(Relative Position Error, km)')
% set(gca,'fontweight','bold','fontsize',20)
%
% colors = distinguishable_colors(nop);
% figure(34);set(gcf,'color','w');
% for i = 1:nop
%     semilogy(time/3600,Pos_Er_MS_ROM(i,:)/1e3,'Color',colors(i,:),'Linewidth',2);hold on;
%     legendInfo{i} = ['Orbit ' num2str(i)];
% end
% for i = 1:nop
%     semilogy(time/3600,Pos_Er_MS(i,:)/1e3,'-.','Color',colors(i,:),'Linewidth',2);hold on;
% end
% grid on;
% legend(legendInfo,'Location','Best');
% xlabel('Time, hrs');ylabel('log_{10}(Relative Position Error, km)')
% set(gca,'fontweight','bold','fontsize',20)

%% Uncertainty on grid
% n_slt = length(sltm); n_lat = length(latm);
% [SLT,LAT]=ndgrid(sltm,latm);
%
% sltx = reshape(SLT,n_slt*n_lat,1);
% latx = reshape(LAT,n_slt*n_lat,1);
% H_SL = zeros(numel(sltx),r);
%
% len = 10;
% figure(8);set(gcf,'Color','w');
% for i = 1:nop
%
%     for ij = 1:numel(sltx)
%         for jk = 1:r
%             H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height(1,i));
%         end
%     end
%     Pyy = diag(H_SL(:,:) * diag(Pv(end-r+1:end,1)) * H_SL(:,:)');
%     Pyyr = reshape(Pyy,n_slt,n_lat);
%     Pyy1 = 100*Pyyr.^0.5;
%     max1 = max(max(Pyy1)); min1 = min(min(Pyy1));
%     subplot(nop,3,3*i-2)
%     contourf(sltm,latm,Pyy1',100,'LineStyle','none');
%     h = colorbar;caxis([min1 max1]);hold on;%ylabel(h,'1\sigma Error (%)','FontSize',20,'Fontweight','bold');hold on;
%     scatter(slt(1,i),lat(1,i),'og','Linewidth',20);hold on
%     scatter(slt(1:len,i),lat(1:len,i),'xr','Linewidth',2);
%     if i < 10; set(gca,'XTickLabel',[]);end
%     if i == 5; ylabel('Latitude'); end; if i == 10; xlabel('Local Time'); end
%     set(gca,'FontSize',24,'Fontweight','bold');
%
%     for ij = 1:numel(sltx)
%         for jk = 1:r
%             H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height(round(end-m/2),i));
%         end
%     end
%     Pyy = diag(H_SL(:,:) * diag(Pv(end-r+1:end,round(end-m/2))) * H_SL(:,:)');
%     Pyyr = reshape(Pyy,n_slt,n_lat);
%     Pyy1 = 100*Pyyr.^0.5;
%     max1 = max(max(Pyy1)); min1 = min(min(Pyy1));
%     subplot(nop,3,3*i-1)
%     contourf(sltm,latm,Pyy1',100,'LineStyle','none');
%     h = colorbar;caxis([min1 max1]);hold on;%ylabel(h,'1\sigma Error (%)','FontSize',20,'Fontweight','bold');hold on;
%     scatter(slt(round(end-m/2),i),lat(round(end-m/2),i),'og','Linewidth',20);hold on
%     scatter(slt(round(end-m/2)-len:round(end-m/2)+len,i),lat(round(end-m/2)-len:round(end-m/2)+len,i),'xr','Linewidth',2);
%     if i < 10; set(gca,'XTickLabel',[]);end
% %     ylabel('Latitude');xlabel('Local Time');
%     if i == 10; xlabel('Local Time'); end
%     set(gca,'FontSize',24,'Fontweight','bold');
%
%     for ij = 1:numel(sltx)
%         for jk = 1:r
%             H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height(end,i));
%         end
%     end
%     Pyy = diag(H_SL(:,:) * diag(Pv(end-r+1:end,end)) * H_SL(:,:)');
%     Pyyr = reshape(Pyy,n_slt,n_lat);
%     Pyy1 = 100*Pyyr.^0.5;
%     max1 = max(max(Pyy1)); min1 = min(min(Pyy1));
%     subplot(nop,3,3*i)
%     contourf(sltm,latm,Pyy1',100,'LineStyle','none');
%     h = colorbar;caxis([min1 max1]);hold on;
%     if i == 5; ylabel(h,'1\sigma Error (%)','FontSize',20,'Fontweight','bold'); end
%     scatter(slt(end,i),lat(end,i),'og','Linewidth',20);hold on
%     scatter(slt(end-len:end,i),lat(end-len:end,i),'xr','Linewidth',2);
%     if i < 10; set(gca,'XTickLabel',[]);end
% %     ylabel('Latitude');xlabel('Local Time');
%     if i == 10; xlabel('Local Time'); end
%     set(gca,'FontSize',24,'Fontweight','bold');
% end
%
