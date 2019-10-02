function [PhiC,Uh,Qrom] = generateROM_JB2008(TA,r)
%generateROM_JB2008 - Compute reduced-order dynamic density model based on 
% JB2008 density data

% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 24-Sep-2019

%------------- BEGIN CODE --------------

%% Reduced-order data
Uh = TA.U100(:,1:r);

%% ROM timesnaps
X1 = TA.densityDataLogVarROM100(1:r,1:end-1);
X2 = TA.densityDataLogVarROM100(1:r,2:end);

%% Space weather inputs: [doy; UThrs; F10; F10B; S10; S10B; XM10; XM10B; Y10; Y10B; DSTDTC; GWRAS; SUN(1); SUN(2)]
U1 = TA.SWdataFull(1:end-1,:)';
% Add future values
U1(15,:) = TA.SWdataFull(2:end,11)'; % DSTDTC
U1(16,:) = TA.SWdataFull(2:end,3)'; % F10
U1(17,:) = TA.SWdataFull(2:end,5)'; % S10
U1(18,:) = TA.SWdataFull(2:end,7)'; % XM10
U1(19,:) = TA.SWdataFull(2:end,9)'; % Y10
% Add quadratic DSTDTC
U1(20,:) = transpose(TA.SWdataFull(1:end-1,11).^2); % DSTDTC^2
U1(21,:) = transpose(TA.SWdataFull(2:end,11).^2); % DSTDTC^2
% % Add mixed terms
U1(22,:) = TA.SWdataFull(1:end-1,11)'.*TA.SWdataFull(1:end-1,3)';
U1(23,:) = TA.SWdataFull(2:end,11)'.*TA.SWdataFull(2:end,3)';

q = size(U1,1);

%% DMDc

% X2 = A*X1 + B*U1 = [A B]*[X1;U1] = Phi*Om
Om = [X1;U1];

% Phi = X2*pinv(Om)
Phi = (Om'\X2')';

A = Phi(1:r,1:r);
B = Phi(1:r,r+1:end);

dth = 1;    %discrete time dt of the ROM in hours
Phi = [A B;zeros(q,r) eye(q)];
PhiC = logm(Phi)/dth;

%% Covariance
X2Pred = A*X1 + B*U1;
errPred = X2Pred-X2;
Qrom = cov(errPred');

end

%------------- END OF CODE --------------
