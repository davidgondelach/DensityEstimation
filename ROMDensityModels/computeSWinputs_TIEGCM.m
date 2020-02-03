function [Inputs] = computeSWinputs_TIEGCM(jd0,jdf,SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM)
%computeSWinputs_NRLMSISE - Compute space weather inputs for ROM-TIEGCM model
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

% Output hourly space weather
tt = jd0:1/24:jdf;
nofPoints = length(tt);

Inputs = zeros(7,nofPoints);
for i=1:nofPoints
    % Date and time
    jdate = tt(i);
    [yy, mm, dd, hh, mnmn, ss] = datevec(jdate-1721058.5);
    [doy] = dayofyear(yy,mm,dd);
    UThrs = hh + mnmn/60 + ss/3600;
    
    % Get TIE-GCM space weather data
    [ f107Average, f107Daily, kp ] = computeSWnrlmsise( SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM, jdate, true );

    % [jdate; doy; UThrs; F107; F107a; Kp]
    Inputs(1,i) = jdate;
    Inputs(2,i) = doy; 
    Inputs(3,i) = UThrs; 
    Inputs(4,i) = f107Daily; 
    Inputs(5,i) = f107Average;
    Inputs(6,i) = kp(2); % 3-hourly Kp
end

% Smooth space weather data
Inputs(4,:) = movmean(Inputs(4,:),[12 11]); % F10: average over 24h
Inputs(5,:) = movmean(Inputs(5,:),[12 11]); % F10: average over 24h
Inputs(6,:) = movmean(Inputs(6,:),3); % kp daily: average over 3h

% Add future (now+1hr) values F10 and Kp
Inputs(7,1:end-1) = Inputs(4,2:end); % F10 (now+1hr)
Inputs(8,1:end-1) = Inputs(6,2:end); % Kp (now+1hr)
[ ~, f107Daily, kp ] = computeSWnrlmsise( SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM, jdf+1/24 );
Inputs(7,end) = f107Daily; % F10 (now+1hr)
Inputs(8,end) = kp(2); % Kp (now+1hr)

% Add quadratic Kp
Inputs(9,:) = Inputs(6,:).^2; % Kp^2 (now)
Inputs(10,:) = Inputs(8,:).^2; % Kp^2 (now+1hr)
% Add mixed terms F10*Kp
Inputs(11,:) = Inputs(6,:).*Inputs(4,:); % Kp*F10 (now)
Inputs(12,:) = Inputs(8,:).*Inputs(7,:); % Kp*F10 (now+1hr)

end

%------------- END OF CODE --------------
