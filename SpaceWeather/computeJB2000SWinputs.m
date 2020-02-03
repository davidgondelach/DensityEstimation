function [MJD,GWRAS,SUN,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC] = computeJB2000SWinputs(year,doy,hour,minute,sec,SOLdata,DTCdata,eopdata)
% Input: Datetime in UTC and space weather data
% Output: space weather proxies in format for JB2008 atmosphere model
%
% This code is licensed under the GNU General Public License version 3.
%
% Based on code by M. Mahooti, 2018
% See https://www.mathworks.com/matlabcentral/fileexchange/56163-jacchia-bowman-atmospheric-density-model
%
% Modified by: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020


%------------- BEGIN CODE --------------

[month,day,~,~,~] = days2mdh(year,doy);
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

% GEOMAGNETIC STORM DTC VALUE
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

et  = cspice_str2et(strcat([num2str(jed2date(UTC+2400000.5),'%d %d %d %d %d %.10f') 'UTC']));
rr_sun = cspice_spkezr('Sun',et,'J2000','NONE', 'Earth');
rr_sun = rr_sun(1:3,1);
ra_Sun  = atan2(rr_sun(2), rr_sun(1));
dec_Sun = atan2(rr_sun(3), sqrt(rr_sun(1)^2+rr_sun(2)^2));

SUN(1)  = ra_Sun;
SUN(2)  = dec_Sun;

end

%------------- END OF CODE --------------
