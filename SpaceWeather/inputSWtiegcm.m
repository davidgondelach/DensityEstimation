function [ SWmatDaily, SWmatMonthlyPred ] = inputSWtiegcm( swfName )
%READSW reads space weather file from CelesTrack
% [  ] = READSW(SWFNAME, JDATE, UTHR)
%
% Inputs for INPUTSW are:
% SWFNAME   :a string that contains space weather name
%
% Outputs for INPUTSW are:
% SWMATDAILY : 
%              matrix for F10.7Daily, F10.7Average, magnetic index
%              Daily observed and predicted AP (8)
%              from 1 Jan 2000 to end of Daily predicted
%
% SWMATMONTHLYPRED : 
%              matrix for Monthly predicted F10.7Daily, F10.7Average
%              Magnetic index and AP (8) from 1 Jan 2000 to end of 
%              Daily predicted
%
%

%% File processing

fid = fopen(swfName,'r+');

% Get rid of initial lines
for i=1:17
    fgetl(fid);
end

% Read number of observed points
str = fgetl(fid);
n_daily_obs = str2double(str(21:25));

% get rid of BEGIN OBSERVED
fgetl(fid);

SWaux = zeros(n_daily_obs, 11);

for i = 1:n_daily_obs
    
    str = fgetl(fid);
    
    SWaux(i, 1) = str2double(str(94:98)); % F10.7 Daily
    
    SWaux(i, 2) = str2double(str(102:106)); % F10.7 Average
    
%     SWaux(i, 3) = str2double(str(80:82)); % Daily Magnetic index
    
    SWaux(i, 4:11) = str2num([str(19:21),str(22:24),str(25:27),str(28:30),...
                         str(31:33),str(34:36),str(37:39),str(40:42)]) ./ 10; % Daily 3h Kp
                     
    if SWaux(i, 1) == 0
        SWaux(i, 1) = SWaux(i, 2);
    end
end

for i=1:3
    str=fgetl(fid);
end

pdt_pnt = str2double(str(28:29));

SWmatDaily = zeros( n_daily_obs + pdt_pnt, 11);
SWmatDaily(1:n_daily_obs, :) = SWaux;

clear SWaux;

% get rid of BEGIN DAILY_PREDICTED
fgetl(fid);

for i = n_daily_obs+1:n_daily_obs+pdt_pnt
    
    str = fgetl(fid);
    
    SWmatDaily(i, 1) =  str2double(str(94:98)); % F10.7 Daily
    
    SWmatDaily(i, 2) = str2double(str(102:106)); % F10.7 Average
    
%     SWmatDaily(i, 3) = str2double(str(80:82)); % Daily Magnetic index
    
    SWmatDaily(i, 4:11) = str2num([str(19:21),str(22:24),str(25:27),str(28:30),...
                         str(31:33),str(34:36),str(37:39),str(40:42)]) ./ 10; % Daily 3h Kp
    
end

for i=1:3
    str=fgetl(fid);
end

mpd_pnt = str2double(str(30:31));

SWmatMonthlyPred = zeros(mpd_pnt, 2);

% get rid of BEGIN MONTHLY_PREDICTED
fgetl(fid);

for i=1:mpd_pnt
    str=fgetl(fid);
    
    SWmatMonthlyPred(i, 1) = str2double(str(94:98)); % F10.7 Daily
    
    SWmatMonthlyPred(i, 2) = str2double(str(102:106)); % F10.7 Average
    
    % Daily Magnetic indeces are not available.

end

fclose(fid);

end