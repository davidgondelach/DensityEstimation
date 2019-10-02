function [Inputs] = computeSWinputs_JB2008(jd0,jdf,eopdata,SOLdata,DTCdata)
%computeSWinputs_JB2008 - Compute space weather inputs for ROM-JB2008 model

% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 24-Sep-2019

%------------- BEGIN CODE --------------

tt = jd0:1/24:jdf;
nofPoints = length(tt);

Inputs = zeros(15,nofPoints);
for i=1:nofPoints
    jdate = tt(i);
    [yyUTC, mmUTC, ddUTC, hhUTC, mnmnUTC, ssUTC] = datevec(jdate-1721058.5);
    doyUTC = day(datetime(yyUTC, mmUTC, ddUTC),'dayofyear');
    UThrs = hhUTC + mnmnUTC/60 + ssUTC/3600;
    % JB2008 space weather
    [~,GWRAS,SUN,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC] = computeJB2000SWinputs(yyUTC,doyUTC,hhUTC,mnmnUTC,ssUTC,SOLdata,DTCdata,eopdata);
    
    % [jdate; doy; UThrs; F10; F10B; S10; S10B; XM10; XM10B; Y10; Y10B; DSTDTC; GWRAS; SUN(1); SUN(2)]
    Inputs(1:15,i) = [jdate; doyUTC; UThrs; F10; F10B; S10; S10B; XM10; XM10B; Y10; Y10B; DSTDTC; GWRAS; SUN(1); SUN(2)];
end
% Smooth DSTDTC
Inputs(12,:) = movmean(Inputs(12,:),12);
% Add future (now+1hr) SW data
Inputs(16,1:end-1) = Inputs(12,2:end); % DSTDTC
Inputs(17,1:end-1) = Inputs(4,2:end); % F10
Inputs(18,1:end-1) = Inputs(6,2:end); % S10
Inputs(19,1:end-1) = Inputs(8,2:end); % XM10
Inputs(20,1:end-1) = Inputs(10,2:end); % Y10

% Add future (now+1hr) SW data for last epoch
[yyUTC, mmUTC, ddUTC, hhUTC, mnmnUTC, ssUTC] = datevec(jdf+1/24-1721058.5);
doyUTC = day(datetime(yyUTC, mmUTC, ddUTC),'dayofyear');
% JB2008 space weather
[~,~,~,F10,~,S10,~,XM10,~,Y10,~,DSTDTC] = computeJB2000SWinputs(yyUTC,doyUTC,hhUTC,mnmnUTC,ssUTC,SOLdata,DTCdata,eopdata);
Inputs(16,end) = DSTDTC; % DSTDTC
Inputs(17,end) = F10; % F10
Inputs(18,end) = S10; % S10
Inputs(19,end) = XM10; % XM10
Inputs(20,end) = Y10; % Y10

% Add quadratic DSTDTC
Inputs(21,:) = Inputs(12,:).^2; % DSTDTC^2 (now)
Inputs(22,:) = Inputs(16,:).^2; % DSTDTC^2 (now+1hr)
% Add mixed terms
Inputs(23,:) = Inputs(12,:).*Inputs(4,:); % DSTDTC*F10 (now)
Inputs(24,:) = Inputs(16,:).*Inputs(17,:); % DSTDTC*F10 (now+1hr)

end

%------------- END OF CODE --------------
