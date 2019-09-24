function [ value, isterminal, direction ] = isdecayed( t, x, noo, svs )
%ISDECAYED Summary of this function goes here
%   Detailed explanation goes here

% km = 1000;

for i = 1:noo
    pos = x(svs*(i-1)+1:svs*(i-1)+3,:);
%     lla = ecef2lla(km*pos');
%     objAlt(i) = lla(:,3)/km;
    objAlt(i) = altitude(pos');
end

value = min(objAlt) - 120;
isterminal = 1;
direction = 0;

end

