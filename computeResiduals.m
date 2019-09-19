function [residuals] = computeResiduals(X_est,time,Meas,stateFnc,measurementFcn,resScaling)

% Compute state at observation times
[Xp] = stateFnc(X_est,time);

% Convert state to measurement space
[Ym] = measurementFcn(Xp);

residuals = Meas-Ym;
residuals = residuals./resScaling;
residuals = reshape(residuals,1,[]);

end