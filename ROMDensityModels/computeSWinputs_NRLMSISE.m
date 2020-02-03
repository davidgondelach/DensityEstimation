function [Inputs] = computeSWinputs_NRLMSISE(jd0,jdf,SWmatDaily,SWmatMonthlyPred)
%computeSWinputs_NRLMSISE - Compute space weather inputs for ROM-NRLMSISE model
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
    
    % Get NRLMSISE-00 space weather data
    [ f107Average, f107Daily, ap ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jdate );
    
    % [jdate; doy; UThrs; F10a; F10; ap]
    Inputs(1,i) = jdate;
    Inputs(2,i) = doy; 
    Inputs(3,i) = UThrs; 
    Inputs(4,i) = f107Average; 
    Inputs(5,i) = f107Daily; 
    Inputs(6:12,i) = ap'; % 7 ap indeces
end

% Smooth space weather data
Inputs(4,:) = movmean(Inputs(4,:),[12 11]); % F10: average over 24h
Inputs(5,:) = movmean(Inputs(5,:),[12 11]); % F10: average over 24h
Inputs(6,:) = movmean(Inputs(6,:),[12 11]); % Ap daily: average over 24h
Inputs(7:12,:) = movmean(Inputs(7:12,:)',3)'; % Ap: average over 3h

% Add future values F10 and Kp
Inputs(13:21,1:end-1) = Inputs(4:12,2:end);
[ f107Average, f107Daily, ap ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jdf+1/24 );
Inputs(13,end) = f107Average;
Inputs(14,end) = f107Daily;
Inputs(15:21,end) = ap';
% % Add quadratic Kp terms
Inputs(22:30,:) = Inputs(4:12,:).^2;
Inputs(31:39,:) = Inputs(13:21,:).^2;
% % Add mixed terms: F10*Kp
Inputs(40,:) = Inputs(5,:).*Inputs(7,:);
Inputs(41,:) = Inputs(14,:).*Inputs(16,:);

end

%------------- END OF CODE --------------
