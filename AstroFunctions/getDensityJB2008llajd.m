function [rho] = getDensityJB2008llajd(lon,lat,alt,jdate)

persistent eopdata SOLdata DTCdata
if isempty(DTCdata)
    [eopdata,SOLdata,DTCdata] = loadJB2008SWdata();
end

% Date and time
[yy, mm, dd, hh, mnmn, ss] = datevec(jdate-1721058.5);
[doy] = dayofyear(yy,mm,dd);

% Space weather data
[MJD,GWRAS,SUN,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC] = computeJB2000SWinputs(yy,doy,hh,mnmn,ss,SOLdata,DTCdata,eopdata);

XLON = deg2rad(lon); % Lon
SAT(1) = mod(GWRAS + XLON, 2*pi);
SAT(2) = deg2rad(lat); % Lat
SAT(3) = alt;
[~,rho] = JB2008(MJD,SUN,SAT,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC);
rho = rho * 1e9; % to kg/km^3

end

