function [tout,posDiff] = computePositionDifference(t1,pos1,t2,pos2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Interpolate COE
if length(t1) < length(t2)
    dataInterp(1,:) = interp1(t2, pos2(1,:), t1,'spline');
    dataInterp(2,:) = interp1(t2, pos2(2,:), t1,'spline');
    dataInterp(3,:) = interp1(t2, pos2(3,:), t1,'spline');
    difference = pos1 - dataInterp;
    
    tout = t1;
else
    dataInterp(1,:) = interp1(t1, pos1(1,:), t2,'spline');
    dataInterp(2,:) = interp1(t1, pos1(2,:), t2,'spline');
    dataInterp(3,:) = interp1(t1, pos1(3,:), t2,'spline');
    difference = pos2 - dataInterp;
    
    tout = t2;
end

posDiff = sqrt(sum(difference.^2,1));

end

