function [x,residuals] = LSQ(x0,residualFun)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

options = optimoptions(@lsqnonlin,'TolX',1e-7,'TolFun',1e-7,'Maxiter', 20,'Display','iter');

[x, resnorm, residuals, exitflag, outputLSQ, ~, Jacobian] = lsqnonlin(residualFun, x0, [-Inf(6,1);0], [Inf(6,1);Inf], options);

end

