function [long,lat,alt,alp,gst]=gc2gd(r,yr,mth,day,hr,min,sec,dt,tf,flag)
%function [long,lat,alt,alp,gst]=gc2gd(r,yr,mth,day,hr,min,sec,dt,tf,flag);
%
% This function converts from Geocentric to Geodetic quantities.
%
%  The inputs are:
%       r = Geocentric altitude (km)
%      yr = year, e.g. 1995
%     mth = month, e.g. Jan=1, Feb=2, etc.
%     day = day, e.g. 1-31
%      hr = hour, e.g. 0-23
%     min = minutes, e.g. 0-59
%     sec = seconds, e.g. 0-59
%      dt = sampling interval (sec)
%      tf = run time (sec)
%    flag = 1 to make longitude from -360 to 360 degrees
%
%  The outputs are:
%    long = longitude (deg)
%     lat = latitude (deg)
%     alt = Geodetic altitude (km)
%     alp = right ascension (deg)
%     gst = Sidereal time (deg)


% Time vector
if dt == 0
    t = 0;
else
    t=[0:dt:tf]';
end

% Flattening of Earth and Earth radius constants
f=1/298.257;
req=6378.14;

% Magnitude of the position
rmag=(r(:,1).^2+r(:,2).^2+r(:,3).^2).^(0.5);

% Altitude
delta=asin(r(:,3)./rmag);
alt=rmag-req*(1-f*sin(delta).^2-(f^2/2)*(sin(2*delta).^2).*(req./rmag-0.25));

% Latitude
sinfd=(req./rmag).*(f*sin(2*delta)+f*f*sin(4*delta).*(req./rmag-0.25));
lat=(delta+asin(sinfd))*180/pi;

% Sideral time at Greenwich
jdate=julian(yr,mth,day,hr,min,sec+t);
tdays=jdate-2415020;
jcent=tdays/36525;
ut=((sec+t)/60/60+min/60+hr)*360/24;
gst=99.6910+36000.7689*jcent+0.0004*jcent.*jcent+ut;

% Longitude
alp=atan2(r(:,2),r(:,1));
alp=unwrap(alp)*180/pi;
long=alp-gst;

% Make longitude from -360 to 360 degrees
ll=long(1);
if (ll<0);
    for k=1:1000,
        if (long(1)>0), break; end
        long=long+360;
    end
else
    for k=1:1000,
        if (long(1)<360), break; end
        long=long-360;
    end
end

if (flag==1),
    for k=1:1000,
        ii=find(long>360);
        if (isempty(ii)==1), break; end
        long=[long(1:ii(1)-1);long(ii(1):length(lat))-360];
    end
    for k=1:1000,
        ii=find(long<-360);
        if (isempty(ii)==1), break; end
        long=[long(1:ii(1)-1);long(ii(1):length(lat))+360];
    end
end
