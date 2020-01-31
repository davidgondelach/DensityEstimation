function [ xp, yp, dut1, lod, ddpsi, ddeps, dat ] = computeEOP_Celestrak( EOPMat, jdate )
%READEOP Reads Earth Orientation Parameters
%   Inptus:
%     EOPfilename
%
%   Outputs:
%     xp:    [arcsec]
%     yp:    [arcsec]
%     dut1:  [s]
%     lod:   [s]
%     ddpsi: [rad]
%     ddeps: [rad]
%     
%     

%% Output initialization
xp = 0;
yp = 0;
dut1 = 0;
lod = 0;
ddpsi = 0;
ddeps = 0;
dat = 26; % s, value for 1992 til July 1

%% Set reference date
% Julian Date of 1962 January 1 00:00:00
jdate0 = 2437665.5;

%% File processing

row = floor(jdate-jdate0)+1;

if row<=size(EOPMat,1)
    
    % Read PM-x
    xp = EOPMat(row,1);
    
    % Read PM-y
    yp = EOPMat(row,2);
    
    % Read UT1-UTC
    dut1 = EOPMat(row,3);
    
    % Read length of day
    lod = EOPMat(row,4);
    
    % Read dPsi
    ddpsi = EOPMat(row,5);
    
    % Read dEps
    ddeps = EOPMat(row,6);
    
end


%% Read DAT file

% Leap seconds (Data from UTC-TAI.history)
DAT = [       10,       11,       12,       13,       14,       15,       16,       17,       18,       19,       20,       21,       22,       23,       24,       25,       26,       27,       28,       29,       30,       31,       32,       33,       34,       35,       36,       37;
       2441317.5,2441499.5,2441683.5,2442048.5,2442413.5,2442778.5,2443144.5,2443509.5,2443874.5,2444239.5,2444786.5,2445151.5,2445516.5,2446247.5,2447161.5,2447892.5,2448257.5,2448804.5,2449169.5,2449534.5,2450083.5,2450630.5,2451179.5,2453736.5,2454832.5,2456109.5,2457204.5,2457754.5;
       2441499.5,2441683.5,2442048.5,2442413.5,2442778.5,2443144.5,2443509.5,2443874.5,2444239.5,2444786.5,2445151.5,2445516.5,2446247.5,2447161.5,2447892.5,2448257.5,2448804.5,2449169.5,2449534.5,2450083.5,2450630.5,2451179.5,2453736.5,2454832.5,2456109.5,2457204.5,2457754.5, Inf];
i=1;

while ~((jdate<DAT(3,i))&&(jdate>=DAT(2,i)))
    i = i+1;
end

dat = DAT(1,i);
   

end

