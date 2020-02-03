function [PhiC,Uh,Qrom] = generateROM_TIEGCM(TA,r)
%generateROM_TIEGCM - Compute reduced-order dynamic density model based on 
% TIE-GCM density data
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020
%
%  Reference:
%  D.J. Gondelach and R. Linares, "Real-Time Thermospheric Density
%  Estimation Via Two-Line-Element Data Assimilation", Space Weather, 2020
%  https://doi.org/10.1029/2019SW002356 or https://arxiv.org/abs/1910.00695
% 

%------------- BEGIN CODE --------------

%% Reduced-order data
Uh = TA.U100(:,1:r);

%% ROM timesnaps
X1 = TA.densityDataLogVarROM100(1:r,1:end-1);
X2 = TA.densityDataLogVarROM100(1:r,2:end);

%% Space weather inputs: [doy; UThrs; F10; F10a; Kp]
U1 = TA.SWdataFull(1:end-1,:)';
% Add future values F10 and Kp
U1(6,:) = TA.SWdataFull(2:end,3)';
U1(7,:) = TA.SWdataFull(2:end,5)';
% Add quadratic Kp
U1(8,:) = transpose(TA.SWdataFull(1:end-1,5).^2);
U1(9,:) = transpose(TA.SWdataFull(2:end,5).^2);
% % Add mixed terms F10*Kp
U1(10,:) = TA.SWdataFull(1:end-1,5)'.*TA.SWdataFull(1:end-1,3)';
U1(11,:) = TA.SWdataFull(2:end,5)'.*TA.SWdataFull(2:end,3)';

q = size(U1,1);

%% DMDc

% X2 = A*X1 + B*U1 = [A B]*[X1;U1] = Phi*Om
Om = [X1;U1];

% Phi = X2*Om^+
Phi = (Om'\X2')';

% Discrete-time dynamic and input matrix
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
