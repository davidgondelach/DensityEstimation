function [ value, isterminal, direction ] = isdecayed( t, x, noo, svs )
%ISDECAYED - Check any object has decayed, i.e. is below 120 km altitude.

for i = 1:noo
    pos = x(svs*(i-1)+1:svs*(i-1)+3,:);
    objAlt(i) = altitude(pos');
end

value = min(objAlt) - 120;
isterminal = 1;
direction = 0;

end

