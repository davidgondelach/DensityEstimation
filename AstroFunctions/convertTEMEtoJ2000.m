function [reci, veci] = convertTEMEtoJ2000(rteme, vteme, jdate)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global EOPMat;
[ ~, ~, dut1, ~, ddpsi, ddeps, dat ] = computeEOP_Celestrak( EOPMat, jdate );

date = jed2date(jdate);
year = date(1);
mon = date(2);
day = date(3);
hr = date(4);
min = date(5);
sec = date(6);

timezone = 0;
[~, ~, ~, ~, ~, ~, ttt, ~, ~, ~, ~, ~, ~, ~, ~ ] ...
         = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );


[reci, veci, ~] = teme2eciNew(rteme, vteme, zeros(3,1), ttt, ddpsi, ddeps);

end

