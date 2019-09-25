function [doy] = dayofyear(year,month,day)

d = datenum(year,month,day);
doy = d - datenum(year,1,1) + 1;

end

