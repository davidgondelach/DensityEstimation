function [MJD,GWRAS,SUN,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC] = computeJB2000SWinputs(year,doy,hour,minute,sec,SOLdata,DTCdata,eopdata)
% Input: Datetime in UTC


% global const %PC

% SAT_Const
% constants
% load DE430Coeff.mat
% PC = DE430Coeff;
% 
% % read Earth orientation parameters
% fid = fopen('eop19620101.txt','r');
% %  ----------------------------------------------------------------------------------------------------
% % |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% % |(0h UTC)           "         "          s          s          "        "          "         "     s 
% %  ----------------------------------------------------------------------------------------------------
% eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
% fclose(fid);
% 
% % read space weather data
% fid = fopen('SOLFSMY.txt','r');
% %  ------------------------------------------------------------------------
% % | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
% %  ------------------------------------------------------------------------
% SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
% fclose(fid);

% % % % % year = 2001;
% % % % % doy = 200;
% year = 2004;
% doy = 32;
[month,day,~,~,~] = days2mdh(year,doy);
% hour = 12;
% minute = 34;
% sec = 23;
MJD = Mjday(year,month,day,hour,minute,sec);

% READ SOLAR INDICES
% USE 1 DAY LAG FOR F10 AND S10 FOR JB2008
JD = floor(MJD-1+2400000.5);
i = find(JD==SOLdata(3,:),1,'first');
SOL = SOLdata(:,i);
F10 = SOL(4);
F10B = SOL(5);
S10 = SOL(6);
S10B = SOL(7);

% USE 2 DAY LAG FOR M10 FOR JB2008
SOL = SOLdata(:,i-1);
XM10 = SOL(8);
XM10B = SOL(9);

% USE 5 DAY LAG FOR Y10 FOR JB2008
SOL = SOLdata(:,i-4);
Y10 = SOL(10);
Y10B = SOL(11);

% % READ GEOMAGNETIC STORM DTC VALUE
% fid = fopen('DTCFILE.txt','r');
% %  ------------------------------------------------------------------------
% % | YYYY DDD   DTC1 to DTC24
% %  ------------------------------------------------------------------------
% DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
% fclose(fid);

doy = finddays(year,month,day,hour,minute,sec);
i = find(year==DTCdata(1,:) & floor(doy)==DTCdata(2,:),1,'first');
DTC = DTCdata(:,i);
ii = floor(hour)+3;
DSTDTC1 = DTCdata(ii,i); %DTC(ii);
if ii >=26 % if hour >= 23
    DSTDTC2 = DTCdata(3,i+1); % Take first hour next day
else
    DSTDTC2 = DTCdata(ii+1,i); % Take next hour same day
end
DSTDTC = interp1([0 3600],[DSTDTC1 DSTDTC2],[minute*60+sec]) + 0.5;

% CONVERT POINT OF INTEREST LOCATION (RADIANS AND KM)
% CONVERT LONGITUDE TO RA
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
[DJMJD0, DATE] = iauCal2jd(year, month, day);
TIME = (60*(60*hour+minute)+sec)/86400;
UTC = DATE+TIME;
TT = UTC+TT_UTC/86400;
TUT = TIME+UT1_UTC/86400;
UT1 = DATE+TUT;
GWRAS = iauGmst06(DJMJD0, UT1, DJMJD0, TT);
% % % % % XLON = 60*const.Rad;
% % % % % SAT(1) = mod(GWRAS + XLON, 2*pi);
% % % % % SAT(2) = -70*const.Rad;
% % % % % SAT(3) = 400;
% XLON = lon*const.Rad; % Lon
% SAT(1) = mod(GWRAS + XLON, 2*pi);
% SAT(2) = lat*const.Rad; % Lat
% SAT(3) = ht;

% SET Sun's right ascension and declination (RADIANS)
% Difference between ephemeris time and universal time
% JD = MJD_UTC+2400000.5;
% [year, month, day, hour, minute, sec] = invjday(JD);
% days = finddays(year, month, day, hour, minute, sec);
% ET_UT = ETminUT(year+days/365.25);
% MJD_ET = MJD_UTC+ET_UT/86400;
% [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
%  r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE430(MJD_ET);

% MJD_TDB = Mjday_TDB(TT);
% [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
%  r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE430(MJD_TDB);
% ra_Sun2  = atan2(r_Sun(2), r_Sun(1));
% dec_Sun2 = atan2(r_Sun(3), sqrt(r_Sun(1)^2+r_Sun(2)^2));

% et  = cspice_str2et(strcat([num2str(jed2date(MJD_TDB+2400000.5),'%d %d %d %d %d %.10f') 'TDB']));
et  = cspice_str2et(strcat([num2str(jed2date(UTC+2400000.5),'%d %d %d %d %d %.10f') 'UTC']));
rr_sun = cspice_spkezr('Sun',et,'J2000','NONE', 'Earth');
rr_sun = rr_sun(1:3,1);
ra_Sun  = atan2(rr_sun(2), rr_sun(1));
dec_Sun = atan2(rr_sun(3), sqrt(rr_sun(1)^2+rr_sun(2)^2));

SUN(1)  = ra_Sun;
SUN(2)  = dec_Sun;

end