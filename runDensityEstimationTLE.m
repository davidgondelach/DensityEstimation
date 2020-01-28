function runDensityEstimationTLE(yr,mth,dy,nofDays,ROMmodel,r,selectedObjects,varargin)
%runDensityEstimationTLE - Estimate thermospheric density using TLE data.
% 
%  runDensityEstimationTLE(yr,mth,dy,nofDays,ROMmodel,r,selectedObjects) -
%     Estimate thermospheric density using TLE data: 
%     1) load TLE data, 2) load ballistic coefficient data, 3) generate 
%     observations from TLE data, 4) generate reduced-order density model,
%     5) initialize reduced-order density state, 6) set initial orbital
%     state guesses, 7) set measurement and process noise, 8) set initial  
%     state covariance, 9) run unscented Kalman filter to simultaneously
%     estimate the orbital states, ballistic coefficients and thermospheric
%     density from TLE-derived orbit observations.
% 
%  runDensityEstimationTLE(yr,mth,dy,nofDays,ROMmodel,r,selectedObjects,plotFigures)
%     Estimate thermospheric density using TLE data and plot results.
% 
%     yr                - start date: year
%     mth               - start date: month
%     dy                - start date: day
%     nofDays           - number of days of estimation window
%     ROMmodel          - name of reduced-order density model
%     r                 - dimension of reduced-order density model
%     selectedObjects   - NORAD IDs of objects used for estimation
%     plotFigures       - boolean: if true then plot results
%
% 
%     Copyright (C) 2019 by David Gondelach
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
%  Author: David Gondelach
%  Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
%  email: davidgondelach@gmail.com
%  Sep 2019; Last revision: 21-Jan-2020
%
%  Reference:
%  D.J. Gondelach and R. Linares, "Real-Time Thermospheric Density
%  Estimation Via Two-Line-Element Data Assimilation", 2019,
%  https://arxiv.org/abs/1910.00695
% 

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
    
    %% Datetime
    % Always start at midnight
    hr=0; mn=0; sc=0; 
    
    jd0 = juliandate(datetime(yr,mth,dy,0,0,0));
    jdf = juliandate(datetime(yr,mth,dy+nofDays,0,0,0));
    
    % Time Interval for measurements
    dt = 3600;
    tf = (jdf-jd0)*24*60*60;
    time = [0:dt:tf]'; m=length(time);
    
    %% Load space weather data
    SWpath = fullfile('Data','SW-All.txt');
    [ SWmatDaily, SWmatMonthlyPred ] = inputSWnrlmsise( SWpath );
    [ SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM ] = inputSWtiegcm( SWpath );
    [eopdata,SOLdata,DTCdata] = loadJB2008SWdata();
    
    %% Get TLE data
    maxAlt = 10000; % Maximum altitude of apogee [km], TLEs for objects with higher apogee will not be downloaded
    jdate0TLEs = juliandate(datetime(yr,mth,1,0,0,0));      % Start date of TLE collection window 
    [yrf, mthf, dyf, ~, ~, ~] = datevec(jdf+30-1721058.5);  % End date of TLE collection window 
    
    % Download or read TLE data
    downloadTLEs = false;
    if downloadTLEs
        username = "[USERNAME]"; % *** Specify your www.space-track.org username here! ***
        password = "[PASSWORD]"; % *** Specify your www.space-track.org password here! ***
        [objects] = downloadTLEsForEstimation(username, password, yr, mth, 1, yrf, mthf, dyf, maxAlt, selectedObjects);
    else
        [objects] = getTLEsForEstimation(yr, mth, 1, yrf, mthf, dyf, selectedObjects);
    end

    
    %% Load BC estimates
    % Ballistic coefficient data: NORAD ID and BC
    BCfilePath = fullfile('Data','BCdata.txt');
    [BCdata] = loadBCdata( BCfilePath );
    
    % Extract orbital data and BC data for selected objects
    for i=1:length(selectedObjects)
        % NORAD ID
        ID = selectedObjects(i);
        % Orbital data (duplicate to ensure data is available)
        newObjects(i) = objects([objects.noradID]==ID);
        % Ballistic coefficient
        BCestimates(i) = BCdata(BCdata(:,1)==ID,2);
    end
    objects = newObjects;
    
    % Number of objects
    nop = length(objects);
    
    % Cell with object ID text strings
    objectIDlabels = cell(1, nop);
    for i=1:nop
        objectIDlabels(i) = {num2str(objects(i).noradID)};
    end

    %% Generate observations from TLE data
    obsEpochs = jd0:dt/86400:jdf;
    [meeMeas] = generateObservationsMEE(objects,obsEpochs,GM_kms);
    
    if plotFigures
        % Plot orbital elements and bstar of TLEs
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
    
    
    %% Load reduced-order density models
    [AC,BC,Uh,F_U,Dens_Mean,M_U,SLTm,LATm,ALTm,maxAtmAlt,SWinputs,Qrom] = generateROMdensityModel(ROMmodel,r,jd0,jdf);
    
    
    %% Generate initial state guess
    % Use orbital state from TLE as initial orbital state guess
    x0orbits = meeMeas(:,1);
    
    % Size of state vector for each object [3xpos,3xvel,1xBC]
    svs = 7;
    
    % Initial state guess: Orbits, BCs and reduced order density
    x0g = zeros(svs*nop+r,1);
    for i = 1:nop
        % Orbital state guesses
        x0g(svs*(i-1)+1:svs*(i-1)+6,1) = x0orbits(6*(i-1)+1:6*(i-1)+6,1);
        % Ballistic coefficient guesses
        x0g(svs*i) = BCestimates(i) * 1000;
    end
    
    % Compute the initial atmosphere state from JB2008 density model
    % Seconds of day in UTC
    UT = hr*3600+mn*60+sc;
    
    % Grid points of ROM model
    sltx = reshape(SLTm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    latx = reshape(LATm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    altx = reshape(ALTm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    
    % Density at grid points according to JB2008 density model
    Dens_JB2008 = zeros(numel(sltx),1);
    for i = 1:numel(sltx)
        % Geographical longitude
        lon = 15*(sltx(i)-UT/3600);
        % Density from JB2008 density model
        Dens_JB2008(i,1) = getDensityJB2008llajd(lon,latx(i),altx(i),jd0) * 1e-9;
    end
    
    % Initialize state for ROM using JB2008
    z0_M = Uh'*(log10(Dens_JB2008)-Dens_Mean); % JB2008 initialization
    
    % Add initial ROM state to initial state guess
    x0g(end-r+1:end,1) = z0_M;
    
    
    %% Measurements and covariance
    % TLE-derived orbit observations in modified equinoctial elements
    Meas = meeMeas;
    
    % Measurement noise
    RM = [];
    for i = 1:nop
        RMfactor = max(objects(i).satrecs(1).ecco/0.004,1);
        RM = [RM; [max(4*objects(i).satrecs(1).ecco,0.0023); RMfactor*3.0e-10; RMfactor*3.0e-10; 1.e-9; 1.e-9; 1e-8]];
    end
    RM = diag(RM); % Convert to matrix with variances on the diagonal
    
    %% Process noise Q and initial state covariance P
    
    Pv = zeros(svs*nop+r,1); % state variance
    Qv = zeros(svs*nop+r,1); % process variance
    for i = 1:nop
        % Initial variance for orbital state in MEE (equal to measurement noise)
        Pv(svs*(i-1)+1) = RM(6*(i-1)+1,6*(i-1)+1);
        Pv(svs*(i-1)+2) = RM(6*(i-1)+2,6*(i-1)+2);
        Pv(svs*(i-1)+3) = RM(6*(i-1)+3,6*(i-1)+3);
        Pv(svs*(i-1)+4) = RM(6*(i-1)+4,6*(i-1)+4);
        Pv(svs*(i-1)+5) = RM(6*(i-1)+5,6*(i-1)+5);
        Pv(svs*(i-1)+6) = RM(6*(i-1)+6,6*(i-1)+6);
        
        % Initial variance for ballistic coefficient
        Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.01)^2; % 1% BC 1-sigma error
        
        % Process noise for orbital state in MEE
        Qv(svs*(i-1)+1) = 1.5e-8; 
        Qv(svs*(i-1)+2) = 2e-14; 
        Qv(svs*(i-1)+3) = 2e-14; 
        Qv(svs*(i-1)+4) = 1e-14; 
        Qv(svs*(i-1)+5) = 1e-14; 
        Qv(svs*(i-1)+6) = 1e-12;
        
        % Process noise for ballistic coefficient
        Qv(svs*(i-1)+7) = 1e-16; % 1-sigma error: 1e-8 per 1 hour
    end
    
    % Initial variance for reduced-order density state
    Pv(end-r+1:end) = (5e0)*ones(r,1);
    Pv(end-r+1) = 2e1; % First mode
    
    % Process noise for reduced-order density state
    % Use variance of ROM model 1-hour prediction error w.r.t. the training data
    Qv(end-r+1:end) = diag(Qrom);
    
    % Initial state covariance and process noise matrices
    P = diag(Pv); % Convert to matrix with variances on the diagonal
    Q = diag(Qv); % Convert to matrix with variances on the diagonal
    
    
    %% Density estimation
    % State estimate
    X_est = x0g; % Initial state guess
    
    % Initial Ephemeris Time (ET): ET is the number of seconds past the 
    % epoch of the J2000 reference frame in the time system known as 
    % Barycentric Dynamical Time (TDB).
    et0  = cspice_str2et(strcat([num2str(jed2date(jd0),'%d %d %d %d %d %.10f') 'UTC']));
    
    % Set state propagation and measurement functions
    stateFnc = @(xx,t0,tf) propagateState_MeeBcRom(xx,t0,tf,AC,BC,SWinputs,r,nop,svs,F_U,M_U,maxAtmAlt,et0,jd0);
    measurementFcn = @(xx) fullmee2mee(xx,nop,svs); % Function to convert from state with MEE,BC,ROM to only MEE
    
    % Run Unscented Kalman filter estimation
    [X_est,Pv] = UKF(X_est,Meas,time,stateFnc,measurementFcn,P,RM,Q);
    Pv(:,1) = diag(P); % Add initial state variance to state variance history
    
    
    %% Plot results
    if plotFigures
        
        % Plot estimated ROM modes and corresponding uncertainty
        figure;
        for i = 1:r
            subplot(ceil(r/4),4,i)
            plot(time/3600,X_est(end-r+i,:),'Linewidth',1); hold on;
            plot(time/3600,X_est(end-r+i,:) + 3*Pv(end-r+i,:).^0.5,'--k','Linewidth',1);
            plot(time/3600,X_est(end-r+i,:) - 3*Pv(end-r+i,:).^0.5,'--k','Linewidth',1);
            xlabel('Time, hrs');ylabel(sprintf('z_%.0f',i));legend('Estimate','3\sigma');
            title(sprintf('Mode %.0f',i));
            axis tight;
            set(gca,'fontsize', 14);
        end
        
        % Plot uncertainty in estimated ROM modes
        figure;
        for i = 1:r
            subplot(ceil(r/4),4,i)
            plot(time/3600,3*Pv(end-r+i,:).^0.5,'k');
            xlabel('Time, hrs'); ylabel(['z_{' num2str(i) '} 3\sigma']);
            title(sprintf('Mode %.0f',i));
        end
        
        % Plot estimated ballistic coefficients and corresponding uncertainty
        figure;
        for i = 1:nop
            subplot(ceil(nop/2),2,i)
            plot(time/3600,X_est(svs*i,:)/1000,'Linewidth',1); hold on;
            plot(time/3600,X_est(svs*i,:)/1000 + 3*Pv(svs*i,:).^0.5/1000,'--k','Linewidth',1);
            plot(time/3600,X_est(svs*i,:)/1000 - 3*Pv(svs*i,:).^0.5/1000,'--k','Linewidth',1);
            xlabel('Time, hrs');ylabel('BC [m^2/kg]');legend('Estimate','3\sigma');
            title(sprintf('Orbit %.0f, BC=%.2f',objects(i).noradID,X_est(svs*i,end)));
            axis tight;
            set(gca,'fontsize', 14);
        end
        
        % Plot uncertainty in estimated ballistic coefficients
        figure;
        for i = 1:nop
            subplot(ceil(nop/2),2,i)
            plot(time/3600,3*Pv(svs*i,:).^0.5./X_est(svs*i,:)*100,'k');
            xlabel('Time, hrs');ylabel('BC 3\sigma [%]');
            title(sprintf('Orbit %.0f, BC=%.2f',objects(i).noradID,X_est(svs*i,end)));
        end
        
        % Plot uncertainty in estimated equinoctial orbital elements
        figure;
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
        
        % Plot position errors w.r.t. TLE measurements
        figure;
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
        
        % Plot uncertainty in estimated density on local solar time v latitude grid
        slt_plot = 0:0.5:24;
        lat_plot = -90:4.5:90;
        [SLT,LAT]=ndgrid(slt_plot,lat_plot);
        sltx = reshape(SLT,length(slt_plot)*length(lat_plot),1);
        latx = reshape(LAT,length(slt_plot)*length(lat_plot),1);
        H_SL = zeros(numel(sltx),r);
        
        figure;
        set(gcf,'Color','w');
        heights = 500:-100:300; % Plot for different altitudes
        nofHeights = length(heights);
        for i = 1:nofHeights
            height = heights(i);
            for ij = 1:numel(sltx)
                for jk = 1:r
                    H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height); % ROM to grid transformation matrix
                end
            end
            Pyy = diag(H_SL(:,:) * diag(Pv(end-r+1:end,end)) * H_SL(:,:)'); % ROM covariance mapped to local solar time and latitude
            Pyyr = reshape(Pyy,length(slt_plot),length(lat_plot)); % Convert to matrix
            Pyy1 = 100*Pyyr.^0.5*log(10); % Convert covariance of log density to standard deviation of density
            % Plot estimated standard deviation (1-sigma error) of density
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
