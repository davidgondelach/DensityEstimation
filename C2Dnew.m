function [PhiC] = C2Dnew(TA)
% Converts the input and dynamic matrices from discrete (1hr model development resolution)
% to continuous and back to discrete time of propagation
%   Detailed explanation goes here

r = TA.r;
q = TA.q;
dth = 1;    %discrete time dt of the ROM in hours
Phi = [TA.A TA.B;zeros(q,r) eye(q)];
PhiC = logm(Phi)/dth;
% PhiD = expm(PhiC*dt); 
% aAD = PhiD(1:r,1:r); aBD = PhiD(1:r,r+1:end);

end

