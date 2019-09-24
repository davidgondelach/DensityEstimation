function [doy] = dayofyear(year,month,day)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d = datenum(year,month,day);
doy = d - datenum(year,1,1) + 1;
end

