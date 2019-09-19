function [rho_rom,varargout] = getDensitiesRom(data,etROM,romStates,r,F_U,M_U,varargin)

computeCov = nargin > 6;
if computeCov
    romCovs = varargin{1};
end

n = length(data.longitudes);
rho_rom = zeros(n,1);
if computeCov
    rho_rom_std = zeros(n,1);
end

for i=1:n
    
    lon = data.longitudes(i);
    slt = data.localtimes(i);
    lat = data.latitudes(i);
    alt = data.altitudes(i);
    
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
    et  = cspice_str2et(strcat([num2str([yr mth dy hrTDT minTDT secTDT],'%d %d %d %d %d %.10f') 'TDT']));
    
    % ROM density
    romState = interp1(etROM, romStates', et);
    slt(slt>24) = slt(slt>24)-24; slt(slt<0) = slt(slt<0)+24;
    UhI = zeros(1,r);
    for j = 1:r
        UhI(:,j) = F_U{j}(slt,lat,alt);
    end
    MI = M_U(slt,lat,alt);
    rho_rom(i) = 10.^(UhI*romState'+MI');
    % ROM covariance
    if computeCov
        romCov = interp1(etROM, romCovs', et);
        Pyy = diag(UhI * diag(romCov) * UhI');
        rho_rom_std(i) = Pyy.^0.5*log(10)*rho_rom(i); % true std
    end
end

if computeCov
    varargout{1} = rho_rom_std;
end

end