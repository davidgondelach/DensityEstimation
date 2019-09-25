function [Inputs] = Comp_Inputs_NRLMSISE_1997_2008(jd0,jdf,SWmatDaily,SWmatMonthlyPred)

tt = jd0:1/24:jdf;
nofPoints = length(tt);

Inputs = zeros(7,nofPoints);
for i=1:nofPoints
    jdate = tt(i);
    [yy, mm, dd, hh, mnmn, ss] = datevec(jdate-1721058.5);
    [doy] = dayofyear(yy,mm,dd);
    UThrs = hh + mnmn/60 + ss/3600;
    
    [ f107Average, f107Daily, ap ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jdate );

    Inputs(1,i) = jdate;
    Inputs(2,i) = doy; 
    Inputs(3,i) = UThrs; 
    Inputs(4,i) = f107Average; 
    Inputs(5,i) = f107Daily; 
    Inputs(6:12,i) = ap';
end

% Add future values F10 and Kp
Inputs(13:21,1:end-1) = Inputs(4:12,2:end);
[ f107Average, f107Daily, ap ] = computeSWnrlmsise( SWmatDaily, SWmatMonthlyPred, jdf+1/24 );
Inputs(13,end) = f107Average;
Inputs(14,end) = f107Daily;
Inputs(15:21,end) = ap';
% % Add quadratic Kp
Inputs(22:30,:) = Inputs(4:12,:).^2;
Inputs(31:39,:) = Inputs(13:21,:).^2;
% % Add mixed terms F10*Kp
Inputs(40,:) = Inputs(5,:).*Inputs(7,:);
Inputs(41,:) = Inputs(14,:).*Inputs(16,:);

end

