function runDensityEstimationTLE(yr,mth,dy,nofDays,DMDmodel,r,selectedObjects,varargin)
%runDensityEstimationTLE - Estimate thermospheric density using TLE data

% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 24-Sep-2019

%------------- BEGIN CODE --------------

if nargin > 7
    plotFigures = varargin{1};
else
    plotFigures = false;
end

try
    
    global mu  % Earth gravitational parameter according to SGP4 model [km^3 s^-2]
    global GM  % Earth gravitational parameter according to accurate gravity model [m^3 s^-2]
    GM_kms = GM*1e-9; % Earth gravitational parameter according to accurate gravity model [km^3 s^-2]
    
    %% Settings
    hr=0; mn=0; sc=0;
    
    jd0 = juliandate(datetime(yr,mth,dy,0,0,0));
    jdf = juliandate(datetime(yr,mth,dy+nofDays,0,0,0));
    
    % Time Interval for measurements
    dt = 3600;
    tf = (jdf-jd0)*24*60*60;
    time = [0:dt:tf]'; m=length(time);
    
    
    SWpath = fullfile('Data','SW-All.txt');
    [ SWmatDaily, SWmatMonthlyPred ] = inputSWnrlmsise( SWpath );
    [ SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM ] = inputSWtiegcm( SWpath );
    
    
    maxAlt = 10000; % km
    [yrf, mthf, dyf, ~, ~, ~] = datevec(jdf+30-1721058.5);
    
    downloadTLEs = true;
    if downloadTLEs
        username = "[USERNAME]"; % www.space-track.org username
        password = "[PASSWORD]"; % www.space-track.org password
        [objects] = downloadTLEsForEstimation(username, password, yr, mth, 1, yrf, mthf, dyf, maxAlt, selectedObjects);
    else
        [objects] = getTLEsForEstimation(yr, mth, 1, yrf, mthf, dyf, selectedObjects);
    end
    jdate0TLEs = juliandate(datetime(yr,mth,1,0,0,0));

    
    %% Load BC estimates
    medianBCestimates = ...
        [11822,	0.00907878223;      13770,	0.0119172584;   15354,	0.0104993335;
        16111,	0.00911779416;      24547,	0.0367465478;   25735,	0.0062344914;
        28822	0.01529856281;      35867,	0.0095876238;   38710,	0.0094608777;
        38998,	0.00797160419;      39452,	0.0052089979;   43020,	0.0190711952;
        24102,	0.028206654;        24487,	0.034716499;    27115,	0.030552165;
        36801,	0.069746786;        39029,	0.026673532;    39454,	0.007618455;
        42735,	0.019392673;        42911,	0.020974571;    43019,	0.010199617;
        41970,  0.04034201043;      42874,  0.015484171;    42880,  0.0149184721; %42735, 42911, 43019 not available at 1 Jan 2017
        01843,	0.0117362;          10363,	0.0134046;      11056,	0.0122993;
        11156,	0.0123989;          11216,	0.0115763;      11269,	0.0126358;
        24645,	0.0124980;          25018,	0.0141775;      25738,	0.0188309;
        07337,  0.01120;            08744,  0.01117;        12138,  0.01115;        12388,  0.01121; % Bowman (2004) - A Method For Computing Accurate Daily Atmospheric Density Values From Satellite Drag Data
    
        00022,  0.02338;            00060,  0.02266;        00063,  0.01486; % Yurasov 2005 - Density Corrections for the NRLMSIS-00
        00165,  0.05326;            00229,  0.05220;        02611,  0.02314;
        14483,  0.01130;            22875,  0.00824;        23853,  0.00846;
        25769,  0.00972;            26929,  0.01542;        26996,  0.01016;
    
        00614,  0.01463;        00750,  0.06846;        01370,  0.11858;        01808,  0.13997; % Emmert 2006 - Thermospheric density 2002?2004
        02016,  0.02879;        02129,  0.04307;        02153,  0.03329;        02622,  0.02240;
        03553,  0.13464;        04221,  0.02201;        04330,  0.02634;        06073,  0.00378;
        20774,  0.01168;        23278,  0.01168;        25233,  0.01734;        26405,  0.00477; %26405=0.00514 (Emmert 2006), 26405=0.00440 (Lu 2017 - Estimation of ballistic coefficients)
        27391,  0.00697;        27392,  0.00693];
    
    for i=1:length(selectedObjects)
        ID = selectedObjects(i);
        indexObject = find([objects.noradID]==ID);
        newObjects(i) = objects(indexObject);
        
        indexBC = find(medianBCestimates(:,1)==ID);
        BCestimates(i) = medianBCestimates(indexBC,2);
    end
    
    objects = newObjects;
    
    nop = length(objects);
    
    objectIDlabels = cell(1, nop);
    for i=1:nop
        objectIDlabels(i) = {num2str(objects(i).noradID)};
    end
    
    %% Check self-consistency of TLE data
    [objectDataSorted, covMEEerrors] = checkSelfConsistencyTLEs(objects, mu, jdate0TLEs, objectIDlabels);

    %% Generate observations from TLE data
    obsEpochs = jd0:dt/86400:jdf;
    [meeMeas] = generateObservationsMEE(objects,obsEpochs,GM_kms);
    initialState = meeMeas(:,1);
    
    if plotFigures
        % Plot orbital elements
        figure;
        for i=1:nop
            subplot(2,3,1); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,[objects(i).satrecs.a],'.'); hold on;
            subplot(2,3,2); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,[objects(i).satrecs.ecco],'.'); hold on;
            subplot(2,3,3); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,rad2deg([objects(i).satrecs.inclo]),'.'); hold on;
            subplot(2,3,4); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,rad2deg([objects(i).satrecs.nodeo]),'.'); hold on;
            subplot(2,3,5); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,rad2deg([objects(i).satrecs.argpo]),'.'); hold on;
            subplot(2,3,6); plot([objects(i).satrecs.jdsatepoch]-jdate0TLEs,([objects(i).satrecs.bstar])/median([objects(i).satrecs.bstar]),'.-'); hold on;
        end
        subplot(2,3,1); xlabel('Days since t_0'); ylabel('a [Earth radii]');legend(objectIDlabels,'Location','northeast');
        subplot(2,3,2); xlabel('Days since t_0'); ylabel('e [-]');legend(objectIDlabels,'Location','northeast');
        subplot(2,3,3); xlabel('Days since t_0'); ylabel('i [deg]'); legend(objectIDlabels,'Location','northeast');
        subplot(2,3,4); xlabel('Days since t_0'); ylabel('\Omega [deg]');legend(objectIDlabels,'Location','northeast');
        subplot(2,3,5); xlabel('Days since t_0'); ylabel('\omega [deg]');legend(objectIDlabels,'Location','northeast');
        subplot(2,3,6); xlabel('Days since t_0'); ylabel('Bstar');legend(objectIDlabels,'Location','northeast');
    end
    
    
    %% Load atmosphere data
    rt = 2*r;
    
    switch DMDmodel
        case 'MSISE2008'
            TA = load('DMDc_NRLMSISE_2008.mat');
            TA = TA.ROMdata;
            % Converting the dynamic and input matrices from discrete to continuous time
            [PhiC,Uh,q] = C2D(TA,rt,r); % q is number of inputs
            % Setup of ROM Modal Interpolation
            n_slt = 24;
            n_lat = 20;
            n_alt = 31;
            sltm=linspace(0,24,n_slt);
            latm=linspace(-90,90,n_lat);
            altm=linspace(100,700,n_alt);
            
            [Inp2] = Comp_Inputs_Var_Celestrak(jd0,jdf+20,SWmatDaily,SWmatMonthlyPred);
            DenS_Mean = TA.DenS_Mean;
            maxAtmAlt = 700;
            
            clear TA;
        case 'TIEGCM2008'
            TA = load('ROM_TIEGCM_2008.mat');
            [PhiC,Uh] = C2Dtiegcm2008(TA,r);
            sltm = TA.localsolartime;   n_slt = length(sltm);
            latm = TA.latitude';        n_lat = length(latm);
            altm = TA.altitude;         n_alt = length(altm);
            
            [Inp2] = Comp_Inputs_Var_Celestrak_TIEGCM(jd0,jdf+20,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);
            % TODO: FIX!!!
            Inp2([2 3],:)=Inp2([3 2],:); % Swap Doy and hour
            
            DenS_Mean = TA.densityDataLogMean;
            maxAtmAlt = 500;
            
            clear TA;
        case 'TIEGCM_1997_2008'
            TA = load('DMDc_1997_2008_Var.mat');
            % Required SW inputs: [F107; Kp; UT; doy]
            
            % Converting the dynamic and input matrices from discrete to continuous time
            [PhiC,Uh,q] = C2D(TA,rt,r); % q is number of inputs
            
            [InpTemp] = Comp_Inputs_Var_Celestrak_TIEGCM(jd0,jdf+20,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);
            Inp2 = [InpTemp(1,:); InpTemp(5,:); InpTemp(6,:); InpTemp(3,:); InpTemp(2,:)];
            
            % Setup of ROM Modal Interpolation
            n_slt = 24;
            n_lat = 20;
            n_alt = 16;
            sltm=linspace(0,24,n_slt);
            latm=linspace(-90,90,n_lat);
            altm=linspace(100,450,n_alt);
            DenS_Mean = TA.DenS_Mean;
            maxAtmAlt = 500;
            
            clear TA TI InpTemp;
        case 'TIE_GCM_ROM'
            TA = load(fullfile('..','HS-DMDc_TIEGCM-2019','TIE_GCM_ROM.mat'));
            
            % Converting the dynamic and input matrices from discrete to continuous time
            [PhiC,Uh,q] = C2D(TA,rt,r); % q is number of inputs
            PhiC = double(PhiC);
            Uh = double(Uh);
            % Required SW inputs: [F107; Kp; UT; doy]
            [InpTemp] = Comp_Inputs_Var_Celestrak_TIEGCM(jd0,jdf+20,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);
            Inp2 = [InpTemp(1,:); InpTemp(5,:); InpTemp(6,:); InpTemp(3,:); InpTemp(2,:)];
            
            % Setup of ROM Modal Interpolation
            n_slt = 24;
            n_lat = 20;
            n_alt = 16;
            sltm=linspace(0,24,n_slt);
            latm=linspace(-90,90,n_lat);
            altm=linspace(100,450,n_alt);
            DenS_Mean = zeros(n_slt*n_lat*n_alt,1);
            maxAtmAlt = 500;
            
            clear TA TI InpTemp;
        case 'JB2008_1999_2010'
            TA = load('JB2008_1999_2010_ROM_r100.mat');
            
            % Converting the dynamic and input matrices from discrete to continuous time
            [PhiC,Uh,Qrom] = C2D_JB2008(TA,r);
            
            [eopdata,SOLdata,DTCdata] = loadJB2008SWdata();
            [Inp2] = Comp_Inputs_JB2008(jd0,jdf+20,eopdata,SOLdata,DTCdata);
            
            % Setup of ROM Modal Interpolation
            sltm=TA.localSolarTimes;
            latm=TA.latitudes;
            altm=TA.altitudes;
            n_slt = length(sltm);
            n_lat = length(latm);
            n_alt = length(altm);
            
            DenS_Mean = TA.densityDataMeanLog;
            maxAtmAlt = 800;
            
            clear TA;
        case 'TIEGCM_1997_2008_new'
            TA = load('TIEGCM_1997_2008_ROM_r100.mat');
            
            % Converting the dynamic and input matrices from discrete to continuous time
            [PhiC,Uh,Qrom] = C2D_TIEGCM_new(TA,r);
            
            TIEGCM_SWdata = TA.SWdataFull;
            [Inp2] = Comp_Inputs_TIEGCM(jd0,jdf,TIEGCM_SWdata);
            
            % Setup of ROM Modal Interpolation
            sltm = TA.localSolarTimes;
            latm = TA.latitudes;
            altm = TA.altitudes;
            n_slt = length(sltm);
            n_lat = length(latm);
            n_alt = length(altm);
            
            DenS_Mean = TA.densityDataMeanLog;
            maxAtmAlt = 500;
            
            clear TA;
        case 'NRLMSISE_1997_2008'
            TA = load('NRLMSISE_1997_2008_ROM_r100_betterSWdata.mat');
            
            % Converting the dynamic and input matrices from discrete to continuous time
            [PhiC,Uh,Qrom] = C2D_NRLMSISE_betterSWdata(TA,r);
            
            [Inp2] = Comp_Inputs_NRLMSISE_1997_2008_betterSWdata(jd0,jdf,SWmatDaily,SWmatMonthlyPred);

            % Setup of ROM Modal Interpolation
            sltm = TA.localSolarTimes;
            latm = TA.latitudes;
            altm = TA.altitudes;
            n_slt = length(sltm);
            n_lat = length(latm);
            n_alt = length(altm);
            
            DenS_Mean = TA.densityDataMeanLog;
            maxAtmAlt = 700;
            
            clear TA;
        otherwise
            warning('No valid DMDc model selected!')
    end
    %%
    [SLTm,LATm,ALTm]=ndgrid(sltm,latm,altm);
    
    F_U{r} = [];
    for i = 1:r
        Uhr = reshape(Uh(:,i),n_slt,n_lat,n_alt); % i-th left singular vector on grid
        F_U{i} = griddedInterpolant(SLTm,LATm,ALTm,Uhr,'linear','linear'); % Create interpolant of Uhr
    end
    Mr = reshape(DenS_Mean,n_slt,n_lat,n_alt);
    M_U = griddedInterpolant(SLTm,LATm,ALTm,Mr,'linear','linear');
    AC = PhiC(1:r,1:r)/3600;
    BC = PhiC(1:r,r+1:end)/3600;
    
    %% Generate initial state guesses
    x0state = initialState;
    
    svs = 7;    %Size of state vector for each object [3xpos,3xvel,1xBC]
    x0g = zeros(svs*nop+r,1);   % initial state guess
    for i = 1:nop
        x0g(svs*(i-1)+1:svs*(i-1)+6,1) = x0state(6*(i-1)+1:6*(i-1)+6,1);
        x0g(svs*i) = BCestimates(i) * 1000;
    end
    
    % Compute the Initial Atmosphere State from MSIS
    UT = hr*3600+mn*60+sc;
    
    % Space weather data
    [ f107A, f107, ap ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jd0 );
    sltx = reshape(SLTm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    latx = reshape(LATm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    altx = reshape(ALTm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    
    % Density
    [eopdata,SOLdata,DTCdata] = loadJB2008SWdata();
    Den = zeros(numel(sltx),1);
    for i = 1:numel(sltx)
        lon = 15*(sltx(i)-UT/3600);
        
        Den_JB2008(i,1) = getDensityJB2008llajd(lon,latx(i),altx(i),jd0) * 1e-9;
    end
    
    % Initial guess for ROM
    z0_M = Uh'*(log10(Den_JB2008)-DenS_Mean); % JB2008 initialization
    
    x0g(end-r+1:end,1) = z0_M;
    
    clear TA TI
    
    %% Measurements and covariance
    
    Meas = meeMeas;
    
    RM = [];
    for i = 1:nop
            RMfactor = max(objects(i).satrecs(1).ecco/0.004,1);
            RM = [RM; [max(4*objects(i).satrecs(1).ecco,0.0023); RMfactor*3.0e-10; RMfactor*3.0e-10; 1.e-9; 1.e-9; 1e-8]];
    end
    RM = diag(RM);
    
    %% Optimization for estimating the process noise and Initial Covariance
    close all
    Pv = zeros(svs*nop+r,1); % state covariance
    Qv = zeros(svs*nop+r,1); % process coveriance
    for i = 1:nop
        Pv(svs*(i-1)+1) = RM(6*(i-1)+1,6*(i-1)+1);
        Pv(svs*(i-1)+2) = RM(6*(i-1)+2,6*(i-1)+2);
        Pv(svs*(i-1)+3) = RM(6*(i-1)+3,6*(i-1)+3);
        Pv(svs*(i-1)+4) = RM(6*(i-1)+4,6*(i-1)+4);
        Pv(svs*(i-1)+5) = RM(6*(i-1)+5,6*(i-1)+5);
        Pv(svs*(i-1)+6) = RM(6*(i-1)+6,6*(i-1)+6);
        Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.005)^2; % 0.5% BC error (std)
        if ismember(objects(i).noradID,[7337;8744;12138;12388;22;60;63;165;229;2611;14483;22875;23853;25769;26929;26996;614;750;1370;1808;2016;2129;2153;2622;3553;4221;4330;6073;20774;23278;25233;26405;27391;27392])
            Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.005)^2; % 0.5% BC error (std) NEW2
        else
            Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.1)^2; % 10% BC error (std)
        end
        Qv(svs*(i-1)+1) = 1.5e-8; %7.2e-4;
        Qv(svs*(i-1)+2) = 2e-14; %6.9e-12;
        Qv(svs*(i-1)+3) = 2e-14; %6.9e-12;
        Qv(svs*(i-1)+4) = 1e-14; %1.5e-12;
        Qv(svs*(i-1)+5) = 1e-14; %1.5e-12;
        Qv(svs*(i-1)+6) = 1e-12; %3.3e-12;
        Qv(svs*(i-1)+7) = 1e-16; % 1-std error is 1e-6 per 1 hour NEW
    end
    Pv(end-r+1:end) = (5e0)*ones(r,1);
    switch DMDmodel
        case 'TIEGCM_1997_2008_new'
            Qv(end-r+1:end) = 1.0*diag(Qrom); %1.0
        case 'JB2008_1999_2010'
            Qv(end-r+1:end) = 1.0*diag(Qrom); %1.0
        case 'NRLMSISE_1997_2008'
            Qv(end-r+1:end) = 1.0*diag(Qrom); %0.5/1.0
        otherwise
            Qv(end-r+1:end) = [0.06;0.03;0.15;0.06;0.017;0.0036;0.006;0.0012;0.0023;0.0035].^2; % These value are actually cov of ROM 1-hr prediction error computed during development, but used here as std: cov = measuredCov^2
    end
    Pv(end-r+1) = 2e1;
    P = diag(Pv);
    Q = diag(Qv);
    
    % Estimated state
    X_est = x0g;
    
    % Set state propagation and measurement functions
    et0  = cspice_str2et(strcat([num2str(jed2date(jd0),'%d %d %d %d %d %.10f') 'UTC']));
    stateFnc = @(xx,t0,tf) propagateStateMEE_FullGravDrag_New(xx,t0,tf,AC,BC,Inp2,r,nop,svs,F_U,M_U,maxAtmAlt,et0,jd0);
    measurementFcn = @(xx) fullmee2mee(xx,nop,svs);
    
    % Run Unscented Kalman filter estimation
    [X_est,Pv] = UKF(X_est,Meas,time,stateFnc,measurementFcn,P,RM,Q);
    Pv(:,1) = diag(P); % Add initial covariance to covariance history
    
    
    %% Plot results
    if plotFigures
        
        % Compute estimated position and velocity at observation epochs
        X_est_pv = X_est;
        for k = 1:nop
            for j=1:size(X_est,2)
                [pos,vel] = ep2pv(X_est((k-1)*svs+1:(k-1)*svs+6,j),GM_kms);
                X_est_pv((k-1)*svs+1:(k-1)*svs+3,j) = pos;
                X_est_pv((k-1)*svs+4:(k-1)*svs+6,j) = vel;
            end
        end
        
        % Plot estimated ROM modes and corresponding uncertainty
        ROMplot2 = figure;
        for i = 1:r
            subplot(ceil(r/4),4,i)
            plot(time/3600,X_est(end-r+i,:),'Linewidth',1); hold on;
            plot(time/3600,X_est(end-r+i,:) + 3*Pv(end-r+i,:).^0.5,'--k','Linewidth',1);
            plot(time/3600,X_est(end-r+i,:) - 3*Pv(end-r+i,:).^0.5,'--k','Linewidth',1);
            xlabel('Time, hrs');ylabel('z');
            title(sprintf('Mode %.0f',i));
            axis tight;
            set(gca,'fontsize', 14);
        end
        
        % Plot uncertainty in estimated ROM modes
        ROMcovplot = figure;
        for i = 1:r
            subplot(ceil(r/4),4,i)
            plot(time/3600,3*Pv(end-r+i,:).^0.5,'k');
            xlabel('Time, hrs');ylabel('z 3\sigma');
            title(sprintf('Mode %.0f',i));
        end
        
        % Plot estimated BC and corresponding uncertainty
        BCplot2 = figure;
        for i = 1:nop
            subplot(ceil(nop/2),2,i)
            plot(time/3600,X_est(svs*i,:)/1000,'Linewidth',1); hold on;
            plot(time/3600,X_est(svs*i,:)/1000 + 3*Pv(svs*i,:).^0.5/1000,'--k','Linewidth',1);
            plot(time/3600,X_est(svs*i,:)/1000 - 3*Pv(svs*i,:).^0.5/1000,'--k','Linewidth',1);
            xlabel('Time, hrs');ylabel('BC [m^2/kg]');
            title(sprintf('Orbit %.0f, BC=%.2f',objects(i).noradID,X_est(svs*i,end)));
            axis tight;
            set(gca,'fontsize', 14);
        end
        
        % Plot uncertainty in estimated BC
        BCcovplot = figure;
        for i = 1:nop
            subplot(ceil(nop/2),2,i)
            plot(time/3600,3*Pv(svs*i,:).^0.5./X_est(svs*i,:)*100,'k');
            xlabel('Time, hrs');ylabel('BC 3\sigma [%]');
            title(sprintf('Orbit %.0f, BC=%.2f',objects(i).noradID,X_est(svs*i,end)));
        end
        
        % Plot difference between observed and estimated orbital elements
        meePlot = figure;
        meeDiffFull = [];
        for i = 1:nop
            meeDiff = meeMeas((i-1)*6+1:i*6,:) - X_est((i-1)*svs+1:(i-1)*svs+6,:);
            meeDiff(6,:) = wrapToPi(meeDiff(6,:));
            for j=1:6
                subplot(2,3,j); hold on;
                plot(time/3600,meeDiff(j,:));
                plot(time/3600,3*sqrt(RM((i-1)*6+j,(i-1)*6+j))*ones(length(time),1),'k');
                plot(time/3600,-3*sqrt(RM((i-1)*6+j,(i-1)*6+j))*ones(length(time),1),'k');
            end
            meeDiffFull = [meeDiffFull meeDiff];
        end
        covMeasErrors = var(meeDiffFull,0,2)';
        subplot(2,3,1); xlabel('Time [hours]'); ylabel('p [km]'); legend(objectIDlabels);
        subplot(2,3,2); xlabel('Time [hours]'); ylabel('f [-]');
        subplot(2,3,3); xlabel('Time [hours]'); ylabel('g [-]');
        subplot(2,3,4); xlabel('Time [hours]'); ylabel('h [-]');
        subplot(2,3,5); xlabel('Time [hours]'); ylabel('k [-]');
        subplot(2,3,6); xlabel('Time [hours]'); ylabel('L [rad]');
        
        % Plot uncertainty in estimated orbits
        meeCovPlot = figure;
        for i = 1:nop
            for j=1:6
                subplot(2,3,j); hold on;
                plot(time/3600,Pv((i-1)*svs+j,:).^0.5);
            end
        end
        subplot(2,3,1); xlabel('Time [hours]'); ylabel('\sigma_p [km]'); legend(objectIDlabels);
        subplot(2,3,2); xlabel('Time [hours]'); ylabel('\sigma_f [-]');
        subplot(2,3,3); xlabel('Time [hours]'); ylabel('\sigma_g [-]');
        subplot(2,3,4); xlabel('Time [hours]'); ylabel('\sigma_h [-]');
        subplot(2,3,5); xlabel('Time [hours]'); ylabel('\sigma_k [-]');
        subplot(2,3,6); xlabel('Time [hours]'); ylabel('\sigma_L [rad]');
        
        % Plot position errors
        posPlot = figure;
        for k = 1:nop
            xx_pv_est = zeros(6,size(X_est,2));
            xx_pv_meas = zeros(6,size(meeMeas,2));
            for j=1:size(X_est,2)
                [pos,vel] = ep2pv(X_est((k-1)*svs+1:(k-1)*svs+6,j),GM_kms);
                xx_pv_est(1:3,j) = pos;
                xx_pv_est(4:6,j) = vel;
                [pos,vel] = ep2pv(meeMeas((k-1)*6+1:(k-1)*6+6,j),GM_kms);
                xx_pv_meas(1:3,j) = pos;
                xx_pv_meas(4:6,j) = vel;
            end
            posErrors = sqrt(sum( (xx_pv_est(1:3,:)-xx_pv_meas(1:3,:)) .^2,1));
            subplot(ceil(nop/2),2,k)
            plot(time/3600,posErrors); hold on;
            xlabel('Time, hrs');ylabel('Position error [km]');
            title(sprintf('Orbit %.0f, mean= %.2f',objects(k).noradID,mean(posErrors)));
        end
        
        % Plot uncertainty on grid
        slt_plot = 0:0.5:24;
        lat_plot = -90:4.5:90;
        [SLT,LAT]=ndgrid(slt_plot,lat_plot);
        
        sltx = reshape(SLT,length(slt_plot)*length(lat_plot),1);
        latx = reshape(LAT,length(slt_plot)*length(lat_plot),1);
        H_SL = zeros(numel(sltx),r);
        
        densCovPlot2 = figure;
        set(gcf,'Color','w');
        heights = 450:-50:350;
        nofHeights = length(heights);
        for i = 1:nofHeights
            height = heights(i);
            for ij = 1:numel(sltx)
                for jk = 1:r
                    H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height);
                end
            end
            Pyy = diag(H_SL(:,:) * diag(Pv(end-r+1:end,end)) * H_SL(:,:)');
            Pyyr = reshape(Pyy,length(slt_plot),length(lat_plot));
            Pyy1 = 100*Pyyr.^0.5*log(10);
            max1 = max(max(Pyy1)); min1 = min(min(Pyy1));
            subplot(nofHeights,1,i)
            contourf(slt_plot,lat_plot,Pyy1',100,'LineStyle','none');
            h = colorbar;caxis([min1 max1]);hold on;
            yticks([-90 0 90]); ylabel('Latitude [deg]');
            title(sprintf('Altitude = %.0f',height));
            if i == ceil(nofHeights/2); ylabel(h,'1\sigma error [%]','FontSize',14); end
            if i < nofHeights; set(gca,'XTickLabel',[]);end
            if i == nofHeights; xlabel('Local solar time'); xticks([0 6 12 18 24]); end
            set(gca,'FontSize',14);
        end
        
    end
    
catch errMsg
    rethrow(errMsg);
end

%------------- END OF CODE --------------