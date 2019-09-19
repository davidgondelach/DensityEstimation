% function runDensityEstimationTLE(yr,mth,dy,nofDays,DMDmodel,r,selectedObjects,varargin)
% 
% if nargin > 7
%     plotFigures = varargin{1};
% else
%     plotFigures = false;
% end
% if nargin > 8
%     gpsFitAndCompare = varargin{2};
% else
%     gpsFitAndCompare = false;
% end

try
    
    % %     clearvars;
    % %     clearvars -global;
    %
    %     addpath( 'AstroFunctions' );
    %     addpath( 'nrlmsise_matlab' );
    %     spicePath = fullfile('/Users','davidgondelach','Documents','mice');
    %     addpath( fullfile(spicePath,'src','mice') );
    %     addpath( fullfile(spicePath,'lib') );
    %
    %     % Results directory
    %     resultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/';
    global mu GM resultsDirPath GPSdataPath
    GM_kms = GM*1e-9;
    
    %% Settings
    % Start date
    % yr = 2018;
    % mth= 1;
    % dy = 1;
    % % yr = 2018;
    % % mth= 1;
    % % dy = 01;
    hr=0; mn=0; sc=0;
    doy = day(datetime(yr,mth,dy),'dayofyear');
    % % Number of days
    % nofDays = 10;
    
    %     jd0 = juliandate(yr,mth,dy,0,0,0);
    jd0 = juliandate(datetime(yr,mth,dy,0,0,0));
    %     jdf = juliandate(yr,mth,dy+nofDays,0,0,0);
    jdf = juliandate(datetime(yr,mth,dy+nofDays,0,0,0));
    
    % Time Interval for measurements
    dt = 3600;
    tf = (jdf-jd0)*24*60*60;
    time = [0:dt:tf]'; m=length(time);
    
    
    %     %% load standard kernels and reference frames
    %     % Clear cspice memory
    %     cspice_kclear;
    %     % Load SPK, PCK, LSK kernels
    %     kernelpath  = fullfile('Data','kernel.txt');
    %     cspice_furnsh( kernelpath );
    %
    %     global gravconst GM
    %     gravconst   = 6.67259e-20; % [km^3/kg/s^2]
    %
    %     gravitymodel     = 'EGM2008';
    %     gravmodeldegree  = 20;
    %     loadGravityModel( gravitymodel, gravmodeldegree );
    %     GM_kms = GM*1e-9;
    %
    SWpath = fullfile('Data','SW-All.txt');
    [ SWmatDaily, SWmatMonthlyPred ] = inputSWnrlmsise( SWpath );
    [ SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM ] = inputSWtiegcm( SWpath );
    %
    %     % Load Earth orientation parameters (needed for TEME to ECI transformation)
    %     global EOPMat
    %     EOPpath = fullfile('Data','EOP-All.txt');
    %     [ EOPMat ] = inputEOP_Celestrak( EOPpath );
    %
    %     %% Load TLEs
    %     % Setup the SGP4 propagator.
    %     global tumin mu radiusearthkm xke j2 j3 j4 j3oj2 opsmode whichconst
    %     opsmode = 'i';
    %     whichconst = 72;
    %     [tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );
    % xpdotp   =  1440.0 / (2.0*pi); % Conversion factor between SGP4 and TLE mean motion units [rev/day]/[rad/min]
    
    % selectedObjects = [16111;38710;35867;13770;15354;39452;11822;24547]; % 10 May 2017
    % selectedObjects = [ 1843;10363;11056;11156;11216;11269;24645;25018]; % 2001 -25738
    % selectedObjects = [16111;38710;35867;13770;15354;39452;11822;24547;41970]; % 10 May 2017
    
    
    maxAlt = 10000; % km
    [yrf, mthf, dyf, ~, ~, ~] = datevec(jdf+30-1721058.5);
    
    downloadTLEs = false;
    if downloadTLEs
        [objects] = downloadTLEsForEstimation(yr, mth, 1, yrf, mthf, dyf, maxAlt, selectedObjects);
    else
        [objects] = getTLEsForEstimation(yr, mth, 1, yrf, mthf, dyf, selectedObjects);
    end
    jdate0TLEs = juliandate(datetime(yr,mth,1,0,0,0));
    %     [objects] = filterTLEs(objects);
%     for i=length(objects):-1:1
%         % Remove corrected TLEs
%         isCorrectedTLE = [diff([objects(i).satrecs.jdsatepoch]) == 0, false];
%         if any(isCorrectedTLE)
%             objects(i).satrecs(isCorrectedTLE) = [];
%         end
%     end
%%
% addpath('/Users/davidgondelach/Documents/PhD/HEO-reentry/TLE-analysis');
%     for i=length(objects):-1:1
%         bstars = [objects(i).satrecs.bstar];
%         jdatesss = [objects(i).satrecs.jdsatepoch];
%         [ subsetBetterPerigeeRadii, subsetNonOutlierIndices ] = filterDataUsingMovingMedian( [objects(i).satrecs.bstar], 21, length(bstars), 10 );
%         figure;
%         plot(jdatesss,bstars,'r.');
%         hold on
%         plot(jdatesss(subsetNonOutlierIndices),bstars(subsetNonOutlierIndices),'b.');
%         title(num2str(objects(i).noradID));
%     end

    
    %%
    % Load BC estimates
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
        %     medianBCestimates = medianBCestimates20172018;
    % selectedObjects = [35867,38998,39452,38710,15354,13770,36801,43020]; % Jan 2018
    % selectedObjects = [35867,38998,39452,38710,15354,13770,25735,43020]; % Jan 2018
    % selectedObjects = [16111;38710;35867;13770;15354;28822;11822;24547]; % May 2017
    % selectedObjects = [16111;38710;35867;13770;15354;39452;11822;24547]; % 10 May 2017
    % selectedObjects = [16111;38710;35867;13770;15354;39452;11822;24547;24487;27115;38998;39029;39454]; % 1 May 2017 - Extended set
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
    
    %%
%     objects(end).satrecs(45)=[];
%     objects(end).satrecs(46)=[];
    [objectDataSorted, covMEEerrors] = checkSelfConsistencyTLEs(objects, mu, jdate0TLEs, objectIDlabels);

    %%
    % TODO: filter out non-consistent TLEs
    % TODO: filter TLEs with negative BC
    
    % Generate observations
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
%     if plotFigures
%         % Plot orbital elements
%         global radiusearthkm;
%         for i=1:nop
%             figure;
%             subplot(2,2,1); plot([objects(i).satrecs.a]*radiusearthkm,'.'); hold on; xlabel('Days since t_0'); ylabel('a [km]');
%             subplot(2,2,2); plot([objects(i).satrecs.a].*(1-[objects(i).satrecs.ecco])*radiusearthkm,'.'); hold on; xlabel('Days since t_0'); ylabel('r_p [km]');
%             subplot(2,2,3); plot(rad2deg([objects(i).satrecs.inclo]),'.'); hold on; xlabel('Days since t_0'); ylabel('i [deg]'); title(sprintf('Orbit %.0f',objects(i).noradID));
%             subplot(2,2,4); plot(rad2deg([objects(i).satrecs.bstar]),'.-'); hold on; xlabel('Days since t_0'); ylabel('Bstar');
%         end
%     end
    % % Plot observations
    % figure;
    % for i=1:nop
    %     subplot(2,3,1); plot(obsEpochs-obsEpochs(1),meeMeas(6*(i-1)+1,:),'.'); hold on;
    %     subplot(2,3,2); plot(obsEpochs-obsEpochs(1),meeMeas(6*(i-1)+2,:),'.'); hold on;
    %     subplot(2,3,3); plot(obsEpochs-obsEpochs(1),meeMeas(6*(i-1)+3,:),'.'); hold on;
    %     subplot(2,3,4); plot(obsEpochs-obsEpochs(1),meeMeas(6*(i-1)+4,:),'.'); hold on;
    %     subplot(2,3,5); plot(obsEpochs-obsEpochs(1),meeMeas(6*(i-1)+5,:),'.'); hold on;
    %     subplot(2,3,6); plot(obsEpochs-obsEpochs(1),meeMeas(6*(i-1)+6,:),'.'); hold on;
    % end
    % subplot(2,3,1); xlabel('Days since t_0'); ylabel('p [km]');
    % subplot(2,3,2); xlabel('Days since t_0'); ylabel('f [-]');
    % subplot(2,3,3); xlabel('Days since t_0'); ylabel('g [-]');
    % subplot(2,3,4); xlabel('Days since t_0'); ylabel('h [-]');
    % subplot(2,3,5); xlabel('Days since t_0'); ylabel('k [-]');
    % subplot(2,3,6); xlabel('Days since t_0'); ylabel('L [rad]');
    
    
    %% Load atmosphere data
    % rt = 20; r = 10;
    rt = 2*r;
    % DMDmodel = 'TIEGCM_1997_2008';
    
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
            % figure;
            % yyaxis left;  plot(Inp2(1,:)-jd0,Inp2(5,:)); ylabel('Daily F10.7');
            % yyaxis right; plot(Inp2(1,:)-jd0,Inp2(6,:)); ylabel('Ap');
            % xlabel('Days since t_0');
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
            
%             [InpTemp] = Comp_Inputs_Var_Celestrak_TIEGCM(jd0,jdf+20,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);
%             Inp2 = InpTemp(1:6,:); % jd,doy,hr,f107a,f107,kp
%             Inp2(7,1:end-1) = InpTemp(5,2:end); %future f107
%             Inp2(8,:) = InpTemp(7,:); %future Kp
%             Inp2(9,:) = Inp2(6,:).^2; %kp^2
%             Inp2(10,:) = Inp2(8,:).^2; %future kp^2
%             Inp2(11,:) = Inp2(6,:).*Inp2(5,:); %kp*f107
%             Inp2(12,:) = Inp2(8,:).*Inp2(7,:); %future kp*f107
            
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
            
%             [Inp2] = Comp_Inputs_NRLMSISE_1997_2008(jd0,jdf,SWmatDaily,SWmatMonthlyPred);
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
        
%         [output] = nrlmsise(altx(i),latx(i),lon,yr,doy,UT,sltx(i),...
%             f107A,f107,ap);
%         Den(i,1) = output.d(6);
        
        Den_JB2008(i,1) = getDensityJB2008llajd(lon,latx(i),altx(i),jd0) * 1e-9;
    end
    
    % Initial guess for ROM
%     z0_M = Uh'*(log10(Den)-DenS_Mean); % NRLMSISE-00 initialization
    z0_M = Uh'*(log10(Den_JB2008)-DenS_Mean); % JB2008 initialization
    
    x0g(end-r+1:end,1) = z0_M;
    
    clear TA TI
    
    %% Measurements and covariance
    
    Meas = meeMeas;
    % RM = diag(covMEEerrors);
    % RM = diag(repmat([0.045; 2.0e-5; 2.0e-5; 2.0e-5; 2.0e-5; 0.000125].^2,nop,1));
    % RM = diag(repmat([0.00246373105635138,3.04847081581516e-10,3.08919216470434e-10,5.73036093238983e-10,5.02344714375781e-10,1.33643157176182e-08],nop,1)); %Meas cov for 10d TLE test
    % RM = diag(repmat([0.0025;4.e-10;4.e-10;6.e-10;6.e-10;1.6e-08],nop,1)); %Meas cov for 10d TLE test
    % RM = diag(repmat([0.0025;4.e-10;4.e-10;6.e-10;6.e-10;1.6e-08],nop,1)); %Meas cov for 10d TLE test
%     RM = diag(repmat([0.0023;3.0e-10;3.0e-10;6.e-10;6.e-10;8.5e-09],nop,1)); % [0.00224124396855864,3.12130369073767e-10,2.89813559884430e-10,5.98683050296782e-10,4.82391313403328e-10,8.32582770502863e-09]
    
    RM = [];
    for i = 1:nop
%         if ismember(objects(i).noradID,[07337,08744,12138,12388])
%             RM = [RM; [10^2*0.0023; 3^2*3.0e-10; 3^2*3.0e-10; 6.e-10; 6.e-10; 8.5e-09]];
%         elseif yr >= 2017
%             RM = [RM; [0.0023; 3.0e-10; 3.0e-10; 6.e-10; 6.e-10; 8.5e-09]];
%         else
%             RM = [RM; [max(objects(i).satrecs(1).ecco,0.0023); 3.0e-10; 3.0e-10; 6.e-10; 6.e-10; 6.25e-08]];
%         else
            RMfactor = max(objects(i).satrecs(1).ecco/0.004,1);
            RM = [RM; [max(4*objects(i).satrecs(1).ecco,0.0023); RMfactor*3.0e-10; RMfactor*3.0e-10; 1.e-9; 1.e-9; 1e-8]];
%             RM = [RM; [max(4*objects(i).satrecs(1).ecco,0.0023); RMfactor*3.0e-10; RMfactor*3.0e-10; 1.e-9; 1.e-9; 4e-8]];
%         end
    end
    RM = diag(RM);
    % RM = 4*RM; %2001
    
    %% Optimization for estimating the process noise and Initial Covariance
    close all
    Pv = zeros(svs*nop+r,1); % state covariance
    Qv = zeros(svs*nop+r,1); % process coveriance
    for i = 1:nop
        %     Pv(svs*(i-1)+1:svs*(i-1)+3) = 1e-4;
        %     Pv(svs*(i-1)+4:svs*(i-1)+6) = 1e-5;
        %     Pv(svs*(i-1)+1) = covMEEerrors(6*(i-1)+1);
        %     Pv(svs*(i-1)+2) = covMEEerrors(6*(i-1)+2);
        %     Pv(svs*(i-1)+3) = covMEEerrors(6*(i-1)+3);
        %     Pv(svs*(i-1)+4) = covMEEerrors(6*(i-1)+4);
        %     Pv(svs*(i-1)+5) = covMEEerrors(6*(i-1)+5);
        %     Pv(svs*(i-1)+6) = covMEEerrors(6*(i-1)+6);
        Pv(svs*(i-1)+1) = RM(6*(i-1)+1,6*(i-1)+1);
        Pv(svs*(i-1)+2) = RM(6*(i-1)+2,6*(i-1)+2);
        Pv(svs*(i-1)+3) = RM(6*(i-1)+3,6*(i-1)+3);
        Pv(svs*(i-1)+4) = RM(6*(i-1)+4,6*(i-1)+4);
        Pv(svs*(i-1)+5) = RM(6*(i-1)+5,6*(i-1)+5);
        Pv(svs*(i-1)+6) = RM(6*(i-1)+6,6*(i-1)+6);
        %     Pv(svs*(i-1)+7) = 7e0;
        Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.005)^2; % 10% BC error (std)
        if ismember(objects(i).noradID,[7337;8744;12138;12388;22;60;63;165;229;2611;14483;22875;23853;25769;26929;26996;614;750;1370;1808;2016;2129;2153;2622;3553;4221;4330;6073;20774;23278;25233;26405;27391;27392])
%             Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.03)^2; % 3% BC error (std)
%             Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.001)^2; % 0.1% BC error (std) NEW
            Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.005)^2; % 0.5% BC error (std) NEW2
        else
            Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.1)^2; % 10% BC error (std)
        end
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
%         Qv(svs*(i-1)+7) = 1e-12; % 1-std error is 1e-6 per 1 hour
        Qv(svs*(i-1)+7) = 1e-16; % 1-std error is 1e-6 per 1 hour NEW
        %     [1.63233439347216e-08;2.41366395075001e-14;1.90396329360770e-14;9.08604764627490e-15;1.04736563945798e-14;1.12575443917672e-12]
    end
    Pv(end-r+1:end) = (5e0)*ones(r,1);
%     Qv(end-r+1:end) = 1e-14*ones(r,1);
    switch DMDmodel
        case 'TIEGCM_1997_2008_new'
%             Qv(end-r+1:end) = [0.026;0.016;0.011;0.032;0.024;0.0062;0.017;0.0069;0.0043;0.020]; % Cov of ROM 1-hr prediction error computed during development
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
    
    % Run Unscented Kalman filter
    [X_est,Pv] = UKF(X_est,Meas,time,stateFnc,measurementFcn,P,RM,Q);
    
    %% Save and plot results
    nowTimeStr = datestr(now,'yymmddHHMMSS');
    filenameBase = [resultsDirPath 'ukf_rom_tle_' DMDmodel '_'];
    testCaseName = [sprintf('%04d',yr), sprintf('%02d',mth), sprintf('%02d',dy), '_', num2str(nofDays) 'd_', num2str(nop), 'obj_'];
    
    Pv(:,1) = diag(P);
    save([filenameBase 'workspace_' testCaseName nowTimeStr]);
    
    %%
    if plotFigures
        
        X_est_pv = X_est;
        for k = 1:nop
            for j=1:size(X_est,2)
                [pos,vel] = ep2pv(X_est((k-1)*svs+1:(k-1)*svs+6,j),GM_kms);
                X_est_pv((k-1)*svs+1:(k-1)*svs+3,j) = pos;
                X_est_pv((k-1)*svs+4:(k-1)*svs+6,j) = vel;
            end
        end
        %
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
        savefig(ROMplot2,[filenameBase 'ROMmodes_' testCaseName nowTimeStr '.fig']);
        
        ROMcovplot = figure;
        for i = 1:r
            subplot(ceil(r/4),4,i)
            plot(time/3600,3*Pv(end-r+i,:).^0.5,'k');
            xlabel('Time, hrs');ylabel('z 3\sigma');
            title(sprintf('Mode %.0f',i));
        end
        savefig(ROMcovplot,[filenameBase 'ROMmodesCov_' testCaseName nowTimeStr '.fig']);
        
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
        savefig(BCplot2,[filenameBase 'BC_' testCaseName nowTimeStr '.fig']);
        
        BCcovplot = figure;
        for i = 1:nop
            subplot(ceil(nop/2),2,i)
            plot(time/3600,3*Pv(svs*i,:).^0.5./X_est(svs*i,:)*100,'k');
            xlabel('Time, hrs');ylabel('BC 3\sigma [%]');
            title(sprintf('Orbit %.0f, BC=%.2f',objects(i).noradID,X_est(svs*i,end)));
        end
        savefig(BCcovplot,[filenameBase 'BCcov_' testCaseName nowTimeStr '.fig']);
        
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
        savefig(meePlot,[filenameBase 'MEEestVmeas_' testCaseName nowTimeStr '.fig']);
        
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
        savefig(meeCovPlot,[filenameBase 'MEEcov_' testCaseName nowTimeStr '.fig']);
        
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
        savefig(posPlot,[filenameBase 'posErr_' testCaseName nowTimeStr '.fig']);
        
%         long_est = zeros(m,nop);lat_est = zeros(m,nop);height_est = zeros(m,nop);
%         for i = 1:nop
%             % Longitude, Latitude and Height of Satellite
%             [long_est(:,i),lat_est(:,i),height_est(:,i)]=gc2gd(X_est_pv(svs*(i-1)+1:svs*(i-1)+3,:)',yr,mth,dy,hr,mn,sc,dt,tf,1);
%         end
        
        % lst_obs = interp1(Inp2(1,:),Inp2(3,:),obsEpochs);
        % rho_est = zeros(size(long_est,1),nop);
        % slt_est = zeros(size(long_est,1),size(long_est,2));
        % densPlot1 = figure;
        % for i = 1:nop
        %     long_est(long_est(:,i)>180,i) = long_est(long_est(:,i)>180,i) - 360;
        %     slt_est(:,i) = lst_obs'+long_est(:,i)/15;
        %     slt_est(slt_est(:,i)>24,i) = slt_est(slt_est(:,i)>24,i)-24; slt_est(slt_est(:,i)<0,i) = slt_est(slt_est(:,i)<0,i)+24;
        %     UhI = zeros(size(long_est,1),r);
        %     for j = 1:r
        %         UhI(:,j) = F_U{j}(slt_est(:,i),lat_est(:,i),height_est(:,i));
        %     end
        %     MI_est = M_U(slt_est(:,i),lat_est(:,i),height_est(:,i));
        %     rho_est(:,i) = 10.^(sum(UhI.*X_est(end-r+1:end,:)',2)+MI_est);
        %     subplot(ceil(nop/2),2,i);
        %     plot(time/3600,rho_est(:,i));
        %     title(sprintf('Orbit %.0f',objects(k).noradID));
        %     xlabel('Time, hrs');ylabel('Density');
        % end
        
        %% Uncertainty on grid
        slt_plot = 0:0.5:24;
        lat_plot = -90:4.5:90;
        [SLT,LAT]=ndgrid(slt_plot,lat_plot);
        
        sltx = reshape(SLT,length(slt_plot)*length(lat_plot),1);
        latx = reshape(LAT,length(slt_plot)*length(lat_plot),1);
        H_SL = zeros(numel(sltx),r);
        
%         densCovPlot = figure;
%         set(gcf,'Color','w');
%         heights = 450:-10:350;
%         nofHeights = length(heights);
%         for i = 1:nofHeights
%             height = heights(i);
%             for ij = 1:numel(sltx)
%                 for jk = 1:r
%                     H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height);
%                 end
%             end
%             Pyy = diag(H_SL(:,:) * diag(Pv(end-r+1:end,1)) * H_SL(:,:)');
%             Pyyr = reshape(Pyy,n_slt,n_lat);
%             Pyy1 = 100*Pyyr.^0.5*log(10);
%             max1 = max(max(Pyy1)); min1 = min(min(Pyy1));
%             subplot(nofHeights,3,3*i-2)
%             contourf(sltm,latm,Pyy1',100,'LineStyle','none');
%             h = colorbar;caxis([min1 max1]);hold on;%ylabel(h,'1\sigma Error (%)','FontSize',20,'Fontweight','bold');hold on;
%             yticks([-90 0 90]);
%             title(sprintf('Altitude = %.0f',height));
%             if i < nofHeights; set(gca,'XTickLabel',[]);end
%             if i == floor(nofHeights/2); ylabel('Latitude'); end
%             if i == nofHeights; xlabel('Local Time'); xticks([0 6 12 18 24]); end
%             %     set(gca,'FontSize',24,'Fontweight','bold');
%             
%             for ij = 1:numel(sltx)
%                 for jk = 1:r
%                     H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height);
%                 end
%             end
%             PvLen = size(Pv,2);
%             Pyy = diag(H_SL(:,:) * diag(Pv(end-r+1:end,ceil(PvLen))) * H_SL(:,:)');
%             Pyyr = reshape(Pyy,n_slt,n_lat);
%             Pyy1 = 100*Pyyr.^0.5*log(10);
%             max1 = max(max(Pyy1)); min1 = min(min(Pyy1));
%             subplot(nofHeights,3,3*i-1)
%             contourf(sltm,latm,Pyy1',100,'LineStyle','none');
%             h = colorbar;caxis([min1 max1]);hold on;%ylabel(h,'1\sigma Error (%)','FontSize',20,'Fontweight','bold');hold on;
%             yticks([-90 0 90]);
%             title(sprintf('Altitude = %.0f',height));
%             if i < nofHeights; set(gca,'XTickLabel',[]);end
%             %     ylabel('Latitude');xlabel('Local Time');
%             if i == nofHeights; xlabel('Local Time'); xticks([0 6 12 18 24]); end
%             %     set(gca,'FontSize',24,'Fontweight','bold');
%             
%             for ij = 1:numel(sltx)
%                 for jk = 1:r
%                     H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height);
%                 end
%             end
%             Pyy = diag(H_SL(:,:) * diag(Pv(end-r+1:end,end)) * H_SL(:,:)');
%             Pyyr = reshape(Pyy,n_slt,n_lat);
%             Pyy1 = 100*Pyyr.^0.5*log(10);
%             max1 = max(max(Pyy1)); min1 = min(min(Pyy1));
%             subplot(nofHeights,3,3*i)
%             contourf(sltm,latm,Pyy1',100,'LineStyle','none');
%             h = colorbar;caxis([min1 max1]);hold on;
%             yticks([-90 0 90]);
%             title(sprintf('Altitude = %.0f',height));
%             if i == ceil(nofHeights/2); ylabel(h,'1\sigma Error (%)','FontSize',20,'Fontweight','bold'); end
%             if i < nofHeights; set(gca,'XTickLabel',[]);end
%             %     ylabel('Latitude');xlabel('Local Time');
%             if i == nofHeights; xlabel('Local Time'); xticks([0 6 12 18 24]); end
%             %     set(gca,'FontSize',24,'Fontweight','bold');
%         end
%         savefig(densCovPlot,[filenameBase 'densCov_' testCaseName nowTimeStr '.fig']);
        
        %
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
        savefig(densCovPlot2,[filenameBase 'densCovLessData_' testCaseName nowTimeStr '.fig']);
        
%         % Final Uncertainty per mode per height
%         set(gcf,'Color','w');
%         heights = 450:-50:350;
%         nofHeights = length(heights);
%         romModeCovPlot = figure;
%         for rom = 1:r
%             for i = 1:nofHeights
%                 height = heights(i);
%                 
%                 for ij = 1:numel(sltx)
%                     for jk = 1:r
%                         H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height);
%                     end
%                 end
%                 Pyy = diag(H_SL(:,rom) * diag(Pv(end-r+rom,end)) * H_SL(:,rom)');
%                 Pyyr = reshape(Pyy,n_slt,n_lat);
%                 Pyy1 = 100*Pyyr.^0.5*log(10);
%                 max1 = max(max(Pyy1)); min1 = min(min(Pyy1));
%                 subplot(r,nofHeights,(rom-1)*nofHeights+i)
%                 contourf(sltm,latm,Pyy1',100,'LineStyle','none');
%                 h = colorbar;caxis([min1 max1]);hold on;
%                 yticks([-90 0 90]);
%                 title(sprintf('Altitude = %.0f',height));
%                 if rom == ceil(r/2) && i == nofHeights; ylabel(h,'1\sigma Error (%)','FontSize',20,'Fontweight','bold'); end
%                 if rom < r; set(gca,'XTickLabel',[]);end
%                 
%                 if rom == r; xlabel('Local Time'); xticks([0 6 12 18 24]); end
%                 
%             end
%         end
%         savefig(romModeCovPlot,[filenameBase 'ROMmodeCov_' testCaseName nowTimeStr '.fig']);
        
    end
    
    if gpsFitAndCompare
        %% ### START COMPARISON WITH GPS DATA ###
        % 1) Fit orbit to GPS data and 2) compare prediction with GPS data
        %
        startDay = dy+nofDays;
        
        etf_ukf = et0 + time(end);
        
        % noradID, planetID, BCestimate
        testObjects(1,:) = {42874, '1056', 0.0154841705917428};
        testObjects(2,:) = {42880, '1047', 0.0149184720940771};
        testObjects(3,:) = {41609, '0E0E', 0.0285353336673982};
        testObjects(4,:) = {42861, '0F4A', 0.0352047559783136};
        testObjects(5,:) = {41970, '100B', 0.0403420104279107};
        
        % planetID = '0E0E'; % May 2017
        % planetID = '0F4A'; % Jan 2018
        % planetID = '100B'; % May 2017
        % planetID = '1056'; % NORAD='42874', BCmedian = 0.0154841705917428
        % planetID = '1047'; % NORAD='42880', BCmedian = 0.0149184720940771
        
        testID = 5;
        planetID = char(testObjects(testID,2));
        BC0 = cell2mat(testObjects(testID,3));
        
        % Load GPS data
        GPSobsWindowDays = 2;
        GPSobsWindowSeconds = GPSobsWindowDays*86400;
        %     [xx_GPS_obs,et_gps_obs,jdate_gps_obs] = getGPSdata(GPSdataPath,planetID,yr,mth,floor(startDay-1),yr,mth,ceil(startDay+GPSobsWindowDays+1));
        %
        %     % Remove GPS data outside observation window
        %     firstGPSindex   = find(et_gps_obs>etf_ukf,1);
        %     lastGPSindex    = find(et_gps_obs<et_gps_obs(firstGPSindex)+GPSobsWindowSeconds,1,'last');
        %     xx_GPS_obs      = xx_GPS_obs(:,firstGPSindex:lastGPSindex);
        %     et_gps_obs      = et_gps_obs(firstGPSindex:lastGPSindex);
        %     jdate_gps_obs   = jdate_gps_obs(firstGPSindex:lastGPSindex);
        
        [xx_GPS_obs,et_gps_obs,jdate_gps_obs] = getGPSdataET(GPSdataPath,planetID,etf_ukf,GPSobsWindowSeconds);
        
        % Initial state guess and observations
        x0 = [xx_GPS_obs(:,1); BC0];
        
        Obs = xx_GPS_obs(1:3,:);
        
        et0_gps  = et_gps_obs(1);
        time_obs = et_gps_obs-et0_gps;
        
        %
        includeSRPMoonSunPert = 1;
        %     stateFnc = @(xx,tt) propagateStatePosVel_FullGravDrag_NRLMSISE(xx,tt,et0_gps,includeSRPMoonSunPert);
        %     measurementFcn = @(xx) pvBC2p(xx);
        %     resScaling = [1;1;1];
        %     %
        %     residualFun = @(xx) computeResiduals(xx,time_obs,Obs,stateFnc,measurementFcn,resScaling);
        %     [x0_lsq_msise,residuals_lsq_msise] = LSQ(x0,residualFun);
        [x0_lsq_msise,residuals_lsq_msise] = performOD_GPS_NRLMSISE(x0,et0_gps,Obs,time_obs,includeSRPMoonSunPert);
        
        %%
        romState = X_est(end-r+1:end,end);
        useOtherROMforPrediction = false;
        if useOtherROMforPrediction
            TA = load('DMDc_1997_2008_Var.mat');
            % Required SW inputs: [F107; Kp; UT; doy]
            % Converting the dynamic and input matrices from discrete to continuous time
            rt_pred = 100; r_pred = 50;
            [PhiC_pred,Uh_pred,q_pred] = C2D(TA,rt_pred,r_pred);
            F_U_pred{r_pred} = [];
            for i = 1:r_pred
                Uhr_pred = reshape(Uh_pred(:,i),n_slt,n_lat,n_alt); % i-th left singular vector on grid
                F_U_pred{i} = griddedInterpolant(SLTm,LATm,ALTm,Uhr_pred,'linear','linear'); % Create interpolant of Uhr
            end
            AC_pred = PhiC_pred(1:r_pred,1:r_pred)/3600;
            BC_pred = PhiC_pred(1:r_pred,r_pred+1:end)/3600;
            
            romState = X_est(end-r+1:end,end);
            romState_pred = Uh_pred' * Uh * romState;
        else
            r_pred = r;
            F_U_pred = F_U;
            AC_pred = AC;
            BC_pred = BC;
            romState_pred = romState;
        end
        
        %     %
        %     tf = et_gps_obs(end) - etf_ukf;
        %     jdate0str  = cspice_et2utc( etf_ukf, 'J', 12 );
        %     jd0       = str2double(jdate0str(4:end));
        %     jdatefstr  = cspice_et2utc( et_gps_obs(end), 'J', 12 );
        %     jdf       = str2double(jdatefstr(4:end));
        %
        %     datetime0_gps = datetime(jdate_gps_obs(1),'convertfrom','juliandate');
        %     doy = day(datetime0_gps,'dayofyear');
        %
        %     jd20170101 = juliandate(datetime(2017,1,1));
        %     jd20171231 = juliandate(datetime(2017,12,31));
        %     jd20180101 = juliandate(datetime(2018,1,1));
        %     jd20181231 = juliandate(datetime(2018,12,31));
        %
        %     switch DMDmodel
        %         case 'MSISE2008'
        %             [Inputs] = Comp_Inputs_Var_Celestrak(jd20170101,jd20181231,SWmatDaily,SWmatMonthlyPred);
        %         case 'TIEGCM2008'
        %             [Inputs] = Comp_Inputs_Var_Celestrak_TIEGCM(jd20170101,jd20181231,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);
        %             % TODO: FIX!!!
        %             Inputs([2 3],:)=Inp2([3 2],:); % Swap Doy and hour
        %         case {'TIEGCM_1997_2008','TIE_GCM_ROM'}
        %             [InpTemp] = Comp_Inputs_Var_Celestrak_TIEGCM(jd20170101,jd20181231,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);
        %             Inputs = [InpTemp(1,:); InpTemp(5,:); InpTemp(6,:); InpTemp(3,:); InpTemp(2,:)];
        %         otherwise
        %             warning('No valid DMDc model selected!')
        %     end
        %     % figure;
        %     % yyaxis left;  plot(Inputs(1,:)-jd0,Inputs(5,:)); ylabel('Daily F10.7');
        %     % yyaxis right; plot(Inputs(1,:)-jd0,Inputs(6,:)); ylabel('Ap');
        %     % xlabel('Days since t_0');
        %
        %     % det = et_gps(end) - etf_ukf;
        %     opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
        %     [romTimeOut,romStateOut]=ode113(@(t,x) ROM_MSISE_ODE(t,x,AC_pred,BC_pred,Inputs,jd0),[0 tf],romState_pred,opts);
        %
        %     romStateTime = [romStateOut, romTimeOut+etf_ukf]; % ROM modes and ephemeris time
        [romStateTime] = predictROMdensity(romState_pred,etf_ukf,et_gps_obs(end),AC_pred,BC_pred,DMDmodel,SWmatDaily,SWmatMonthlyPred,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);
        
        stateFnc = @(x,t) propagateStatePosVel_FullGravDrag_ROM(x,t,et0_gps,romStateTime,r_pred,F_U_pred,M_U,includeSRPMoonSunPert);
        measurementFcn = @(xx) pvBC2p(xx);
        resScaling = [1;1;1];
        
        residualFun = @(xx) computeResiduals(xx,time_obs,Obs,stateFnc,measurementFcn,resScaling);
        [x0_lsq_rom,residuals_lsq_rom] = LSQ(x0,residualFun);
        %%
        
        
        if plotFigures
            OrbitFitPlot = figure; hold on;
            plot(time_obs/3600,sqrt(sum(reshape(residuals_lsq_rom,3,[]).^2,1)));
            plot(time_obs/3600,sqrt(sum(reshape(residuals_lsq_msise,3,[]).^2,1)));
            xlabel('Time [hr]'); ylabel('Position residual [km]'); set(gca,'FontSize',14);
            legend('ROM','NRLMSISE');
            savefig(OrbitFitPlot,[filenameBase testCaseName nowTimeStr '_' planetID '_GPSorbitFit_ROMvMSISE.fig']);
        end
        
        %
        sampleTime      = [0:5:et_gps_obs(end)-et0_gps];
        [xxfull_rom]    = stateFnc(x0_lsq_rom,sampleTime);
        [xxfull_msise]  = propagateStatePosVel_FullGravDrag_NRLMSISE(x0_lsq_msise,sampleTime,et0_gps,includeSRPMoonSunPert);
        
        rho_rom = zeros(1,size(xxfull_rom,2));
        rho_msise = zeros(1,size(xxfull_msise,2));
        for j=1:size(xxfull_rom,2)
            rho_rom(j) = getDensityROM(xxfull_rom(1:3,j)',sampleTime(j)+et0_gps,romStateTime,r_pred,F_U_pred,M_U);
        end
        for j=1:size(xxfull_msise,2)
            rho_msise(j) = getDensityNRLMSISE(xxfull_rom(1:3,j)',sampleTime(j)+et0_gps);
        end
        
        if plotFigures
            OrbitFitDensityPlot = figure;
            plot(sampleTime/3600,rho_rom*1e-9); hold on;
            plot(sampleTime/3600,rho_msise*1e-9);
            xlabel('Time [hr]'); ylabel('Density [kg/m^3]'); set(gca,'FontSize',14);
            legend('ROM','NRLMSISE');
            savefig(OrbitFitDensityPlot,[filenameBase testCaseName nowTimeStr '_' planetID '_GPSorbitFit_ROMvMSISEdensityAlongOrbit.fig']);
        end
        
        % Load GPS data
        predictionWindowDays = 10;
        predictionWindowSeconds = predictionWindowDays*86400;
        %     [xx_GPS_pred,et_gps_pred,jdate_gps_pred] = getGPSdata(GPSdataPath,planetID,yr,mth,floor(startDay-1),yr,mth,ceil(startDay+predictionWindowDays+1));
        %
        %     % Remove GPS data outside observation window
        %     firstGPSindex   = find(et_gps_pred>etf_ukf,1);
        %     lastGPSindex    = find(et_gps_pred<et_gps_pred(firstGPSindex)+predictionWindowSeconds,1,'last');
        %     xx_GPS_pred      = xx_GPS_pred(:,firstGPSindex:lastGPSindex);
        %     et_gps_pred      = et_gps_pred(firstGPSindex:lastGPSindex);
        %     jdate_gps_pred   = jdate_gps_pred(firstGPSindex:lastGPSindex);
        [xx_GPS_pred,et_gps_pred,jdate_gps_pred] = getGPSdataET(GPSdataPath,planetID,etf_ukf,predictionWindowSeconds);
        
        %
        % [romTimeOut_pred,romStateOut_pred]=ode113(@(t,x) ROM_MSISE_ODE(t,x,AC,BC,Inputs,jd0),[0 et_gps_pred(end) - etf_ukf],romState,opts);
        %     opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
        %     [romTimeOut_pred,romStateOut_pred]=ode113(@(t,x) ROM_MSISE_ODE(t,x,AC_pred,BC_pred,Inputs,jd0),[0 et_gps_pred(end) - etf_ukf],romState_pred,opts);
        %     romStateTime_pred2 = [romStateOut_pred, romTimeOut_pred+etf_ukf]; % ROM modes and ephemeris time
        [romStateTime_pred] = predictROMdensity(romState_pred,etf_ukf,et_gps_pred(end),AC_pred,BC_pred,DMDmodel,SWmatDaily,SWmatMonthlyPred,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);
        
        stateFnc = @(t,x) propagateStatePosVel_FullGravDrag_ROM(t,x,et0_gps,romStateTime_pred,r_pred,F_U_pred,M_U,includeSRPMoonSunPert);
        
        sampleTime_pred  = et_gps_pred-et0_gps;
        [xx_rom_pred]    = stateFnc(x0_lsq_rom,sampleTime_pred);
        [xx_msise_pred]  = propagateStatePosVel_FullGravDrag_NRLMSISE(x0_lsq_msise,sampleTime_pred,et0_gps,includeSRPMoonSunPert);
        
        % Convert state to measurement space
        [pos_rom_pred] = measurementFcn(xx_rom_pred);
        residuals_rom_pred = xx_GPS_pred(1:3,:)-pos_rom_pred;
        residuals_rom_pred = reshape(residuals_rom_pred,3,[]);
        
        [pos_msise_pred] = measurementFcn(xx_msise_pred);
        residuals_msise_pred = xx_GPS_pred(1:3,:)-pos_msise_pred;
        residuals_msise_pred = reshape(residuals_msise_pred,3,[]);
        
        for i=1:size(xx_GPS_pred,2)
            [cart2rtnMatrix] = computeCart2RTNMatrix(xx_GPS_pred(1:3,i), xx_GPS_pred(4:6,i));
            rr_Diff_rom_RTN(:,i) = cart2rtnMatrix*residuals_rom_pred(:,i);
            rr_Diff_msise_RTN(:,i) = cart2rtnMatrix*residuals_msise_pred(:,i);
        end
        
        if plotFigures
            OrbitPredictionPlot = figure; hold on;
            plot(sampleTime_pred/3600,sqrt(sum(reshape(residuals_rom_pred,3,[]).^2,1)));
            % plot(sampleTime_pred/3600,rr_Diff_rom_RTN(1,:));
            % plot(sampleTime_pred/3600,rr_Diff_rom_RTN(2,:));
            % plot(sampleTime_pred/3600,rr_Diff_rom_RTN(3,:));
            plot(sampleTime_pred/3600,sqrt(sum(reshape(residuals_msise_pred,3,[]).^2,1)));
            % plot(sampleTime_pred/3600,rr_Diff_msise_RTN(1,:));
            % plot(sampleTime_pred/3600,rr_Diff_msise_RTN(2,:));
            % plot(sampleTime_pred/3600,rr_Diff_msise_RTN(3,:));
            plot(ones(1,2)*(et_gps_obs(end)-et0_gps)/3600, [0 10], 'k--');
            xlabel('Time [hr]'); ylabel('Position error [km]'); set(gca,'FontSize',14);
            legend('ROM','NRLMSISE','End of fitting window','Location','NorthWest');
            title(['Position error w.r.t. GPS data (object ',planetID,')']);
            savefig(OrbitPredictionPlot,[filenameBase testCaseName nowTimeStr '_' planetID '_GPSorbitFit_ROMvMSISE_prediction.fig']);
        end
        
        %
        sampleTime_pred_full = [0:10:et_gps_pred(end)-et0_gps];
        [xxfull_rom_pred]    = stateFnc(x0_lsq_rom,sampleTime_pred_full);
        [xxfull_msise_pred]  = propagateStatePosVel_FullGravDrag_NRLMSISE(x0_lsq_msise,sampleTime_pred_full,et0_gps,includeSRPMoonSunPert);
        
        rho_rom_pred = zeros(1,size(xxfull_rom_pred,2));
        rho_msise_pred = zeros(1,size(xxfull_msise_pred,2));
        for j=1:size(xxfull_rom_pred,2)
            rho_rom_pred(j) = getDensityROM(xxfull_rom_pred(1:3,j)',sampleTime_pred_full(j)+et0_gps,romStateTime_pred,r,F_U,M_U);
        end
        for j=1:size(xxfull_msise_pred,2)
            rho_msise_pred(j) = getDensityNRLMSISE(xxfull_msise_pred(1:3,j)',sampleTime_pred_full(j)+et0_gps);
        end
        
        if plotFigures
            OrbitPredDensityPlot = figure;
            plot(sampleTime_pred_full/3600,rho_rom_pred*1e-9); hold on;
            plot(sampleTime_pred_full/3600,rho_msise_pred*1e-9);
            xlabel('Time [hr]'); ylabel('Density [kg/m^3]'); set(gca,'FontSize',14);
            legend('ROM','NRLMSISE');
            savefig(OrbitPredDensityPlot,[filenameBase testCaseName nowTimeStr '_' planetID '_GPSorbitFit_ROMvMSISEdensityAlongOrbitPred.fig']);
        end
        
        save([filenameBase 'workspace_' testCaseName nowTimeStr '_' planetID '_GPSorbitFit'], ...
            'gravmodeldegree','planetID','GPSobsWindowDays','xx_GPS_obs','et_gps_obs','jdate_gps_obs', ...
            'x0','Obs','et0_gps','time_obs','resScaling','x0_lsq_msise','residuals_lsq_msise', ...
            'jd0','jdf','doy','romState','romStateTime','x0_lsq_rom','residuals_lsq_rom','predictionWindowDays','xx_GPS_pred', ...
            'et_gps_pred','jdate_gps_pred','romStateTime_pred','sampleTime_pred','xx_rom_pred', ...
            'xx_msise_pred','residuals_rom_pred','residuals_msise_pred','sampleTime_pred_full');
        
    end
    
    
catch errMsg
    rethrow(errMsg);
end

% end