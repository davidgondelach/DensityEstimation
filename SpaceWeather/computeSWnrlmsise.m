function [ f107Average, f107Daily, magneticIndex ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jdate, varargin )
%COMPUTESWNRLMSISE reads space weather file from CelesTrack
% [  ] = COMPUTESWNRLMSISE(SWMATDAILY, SWMATMONTHLYPRED, JDATE, UTHR)
%
% Inputs for COMPUTESWNRLMSISE are:
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
% JDATE      : Julian Date
%
% USETODAYSF107 :
%              boolean: if true then use today's F10.7 value else use
%              yesterday's F10.7 value
%

%% Outputs initializations (sets default values for atmosnrlmsise00)
magneticIndex = 4*ones(1,7);
f107Daily = 150;
f107Average = 150;

%% Internal variables definition

% Julian Date of 1957 10 01 00:00:00
jdate1957 = 2436112.5;

%% Determine UT hour
% [UTyr,UTmo,~,UThr,~,~] = invjday(jdate);
% date = datetime(jdate,'ConvertFrom','juliandate');
% [UTyr, UTmo, ~, UThr, ~, ~] = datevec(date);
[UTyr, UTmo, ~, UThr, ~, ~] = datevec(jdate-1721058.5);

%% File processing
auxMI = zeros(1,32);

row = floor(jdate-jdate1957)+1;

if nargin>3 && varargin{1}
    rowf107Daily = row; % Use today's F10.7 flux
else
    rowf107Daily = row-1; % Use yesterday's F10.7 flux (needed for NRLMSISE-00)
end

if row <= 1
    
    error('No space weather data available for this date: \n JD %s',num2str(jdate));

elseif row <= size(SWmatDaily,1)
    
    f107Daily =  SWmatDaily(rowf107Daily, 1); % Daily F10.7 flux for previous day
    
    f107Average = SWmatDaily(row, 2); % 81 day average of F10.7 flux (centered on day)
    
    magneticIndex(1) = SWmatDaily(row, 3); % Daily
    column = ceil((24-UThr)/3);
    
    for i=1:4
        
        auxMI(i*8-7:i*8) = SWmatDaily(row-(4-i), 4:11);

    end
        
    auxMI = fliplr(auxMI);
    
    magneticIndex(2:5) = auxMI(column:column+3); % Now:3hr:Now-9hr
    
    magneticIndex(6) = mean(auxMI(column+4:column+11));
    magneticIndex(7) = mean(auxMI(column+12:column+19));
    
else
    
    % Determine UTyr and UTmo of Monthly prediction beginning
    [UTyrMP, UTmoMP, ~, ~, ~, ~] = invjday(jdate1957 + size(SWmatDaily,1));
    
     dUTyr = UTyr-UTyrMP;
     dUTmo = UTmo-UTmoMP;
        
     dmon = dUTmo + dUTyr*12;
        
     if dmon<0
         dmon=0;
     end
       
     if dmon<size(SWmatMonthlyPred,1);
        
        f107Daily = SWmatMonthlyPred(dmon+1,1);
        f107Average = SWmatMonthlyPred(dmon+1,2);
         
     end
        
end


end

