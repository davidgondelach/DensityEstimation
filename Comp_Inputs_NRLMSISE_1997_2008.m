function [Inputs] = Comp_Inputs_NRLMSISE_1997_2008(jd0,jdf,SWmatDaily,SWmatMonthlyPred)

% Compute the Inputs to take the state to meaurement epoch and for
% computing the true trajectory

tt = jd0:1/24:jdf;
nofPoints = length(tt);

Inputs = zeros(7,nofPoints);
for i=1:nofPoints
    jdate = tt(i);
    [yy, mm, dd, hh, mnmn, ss] = datevec(jdate-1721058.5);
    [doy] = dayofyear(yy,mm,dd);
    UThrs = hh + mnmn/60 + ss/3600;
    
    [ f107Average, f107Daily, ap ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jdate, true );

    Inputs(1,i) = jdate;
    Inputs(2,i) = doy; 
    Inputs(3,i) = UThrs; 
    Inputs(4,i) = f107Average; 
    Inputs(5,i) = f107Daily; 
    Inputs(6,i) = ap(2);
end
% Add future (now+1hr) values F10 and Ap
Inputs(7,1:end-1) = Inputs(5,2:end); % F10 (now+1hr)
Inputs(8,1:end-1) = Inputs(6,2:end); % Ap (now+1hr)
[ ~, f107Daily, ap ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jdf+1/24, true );
Inputs(7,end) = f107Daily;
Inputs(8,end) = ap(2);

% Add quadratic Ap
Inputs(9,:) = Inputs(6,:).^2; % Ap^2 (now)
Inputs(10,:) = Inputs(8,:).^2; % Ap^2 (now+1hr)
% % Add mixed terms F10*Ap
Inputs(11,:) = Inputs(6,:).*Inputs(5,:); % Ap*F10 (now)
Inputs(12,:) = Inputs(8,:).*Inputs(7,:); % Ap*F10 (now+1hr)


end

