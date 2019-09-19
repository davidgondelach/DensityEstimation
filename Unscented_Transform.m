function [Wm,Wc,L,lambda] = Unscented_Transform(x0f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

L=length(x0f);
alpha = 1;
beta = 2;
kappa = 3 - L; %0;
lambda = alpha^2*(L + kappa) - L;

W0m = lambda/(L + lambda);
W0c = lambda/(L + lambda)+(1-alpha^2+beta);
Wim = 1/(2*(L + lambda));

Wm = [W0m Wim+zeros(1,2*L)];
Wc = [W0c Wim+zeros(1,2*L)];

end

