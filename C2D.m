function [PhiC,Uh,q] = C2D(TA,rt,r)
% Converts the input and dynamic matrices from discrete (1hr model development resolution)
% to continuous and back to discrete time of propagation
%   Detailed explanation goes here

n=size(TA.U,1);q=size(TA.Ups,1);
Dii  = diag(1./diag(TA.Di));
Uh = TA.U(:,1:r); % left singular vectors
aA = Uh'*TA.DSOm*TA.Ui(:,1:rt)*Dii(1:rt,1:rt)*TA.Ui(1:n,1:rt)'*Uh;
aB = Uh'*TA.DSOm*TA.Ui(:,1:rt)*Dii(1:rt,1:rt)*TA.Ui(n+1:n+q,1:rt)';
dth = 1;    %discrete time dt of the ROM in hours
Phi = [aA aB;zeros(q,r) eye(q)];
PhiC = logm(Phi)/dth;
% PhiD = expm(PhiC*dt); 
% aAD = PhiD(1:r,1:r); aBD = PhiD(1:r,r+1:end);

end

