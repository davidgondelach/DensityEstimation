function jd = julian(yr,mth,day,hr,min,sec);
%function jd = julian(yr,mth,day,hr,min,sec);
%
% This function computes the Julian date from calendar date.
%
%  The inputs are:
%     yr  = year, e.g. 1995
%     mth = month, e.g. Jan=1, Feb=2, etc.
%     day = day, e.g. 1-31
%     hr  = hour, e.g. 0-23
%     min = minutes, e.g. 0-59
%     sec = seconds, e.g. 0-59
%
%  The output is:
%      jd = Julian date
%         =  number of mean solar days (UT) since 4713 BC Jan 1, 12h UT


% compute day-of-year
doy = day + 31*(mth-1) - fix(2.2+0.4*mth).*(mth>2) + (mth>2).*(~rem(yr,4));

% compute Julian date
jd = 2415020 + (yr-1900)*365 + fix((yr-1901)/4) + (doy-1) + 0.5 + ...
       (3600*hr + 60*min + sec) / 86400;
