function [densityMeas,sltMeas,latMeas,altMeas,UhMeas,MMeas] = generateDensitiesObservations(data,et0,time,r,F_U,M_U)

% Smooth density data
densities = movmean(data.densities,5); % Original data is spaced by approx 47s on average

% Get ephemeris time corresponding to density data
et_dens = zeros(length(data.seconds),1);
for i=1:length(data.seconds)
    secondsGPS = data.seconds(i);
    secondsTDT = secondsGPS + 19 + 32.184; % TDT = TAI + 32.184, TAI = GPS time + 19
    [hrTDT,minTDT,secTDT] = sec2hms( secondsTDT );
    yr = data.years(i) + 2000;
    doy_ = data.doys(i);
    if hrTDT >= 24
        hrTDT = hrTDT - 24;
        doy_ = doy_+1;
    end
    dv = datevec(datenum(yr, 1, doy_));
    mth = dv(2);
    dy = dv(3);
    et_dens(i)  = cspice_str2et(strcat([num2str([yr mth dy hrTDT minTDT secTDT],'%d %d %d %d %d %.10f') 'TDT']));
end

% Interpolate data at measurement epochs
et = et0 + time;
densityMeas = interp1(et_dens,densities,et);
sltMeas = interp1(et_dens,data.localtimes,et);
latMeas = interp1(et_dens,data.latitudes,et);
altMeas = interp1(et_dens,data.altitudes,et);

UhMeas = zeros(length(et),r);
for j = 1:r
    UhMeas(:,j) = F_U{j}(sltMeas,latMeas,altMeas);
end
MMeas = M_U(sltMeas,latMeas,altMeas);

end