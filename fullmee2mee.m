function [mee] = fullmee2mee(Xp,nop,svs)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% mu = 398600.4415; Re = 6378.1363; J2 = 1.082626925638815*1e-3;

mee = zeros(nop*6,size(Xp,2));
for k = 1:nop
    mee(6*(k-1)+1:6*k,:) = Xp(svs*(k-1)+1:svs*(k-1)+6,:);
end

end

