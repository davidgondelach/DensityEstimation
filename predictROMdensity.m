function [romStateTime] = predictROMdensity(romState_pred,et0,etf,AC_pred,BC_pred,DMDmodel,SWdata1,SWdata2,SWdata3,SWdata4)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

tf = etf - et0;
jdate0str  = cspice_et2utc( et0, 'J', 12 );
jd0       = str2double(jdate0str(4:end));
jdatefstr  = cspice_et2utc( etf, 'J', 12 );
jdf       = str2double(jdatefstr(4:end));

jd0inputs = floor(jd0-1);
jdfinputs = ceil(jdf+1);

switch DMDmodel
    case 'MSISE2008'
        [Inputs] = Comp_Inputs_Var_Celestrak(jd0inputs,jdfinputs,SWdata1,SWdata2);
    case 'TIEGCM2008'
        [Inputs] = Comp_Inputs_Var_Celestrak_TIEGCM(jd0inputs,jdfinputs,SWdata3,SWdata4);
        % TODO: FIX!!!
        Inputs([2 3],:)=Inp2([3 2],:); % Swap Doy and hour
    case {'TIEGCM_1997_2008','TIE_GCM_ROM'}
        [InpTemp] = Comp_Inputs_Var_Celestrak_TIEGCM(jd0inputs,jdfinputs,SWdata3,SWdata4);
        Inputs = [InpTemp(1,:); InpTemp(5,:); InpTemp(6,:); InpTemp(3,:); InpTemp(2,:)];
    case 'JB2008_1999_2010'
        [Inputs] = Comp_Inputs_JB2008(jd0inputs,jdfinputs,SWdata1,SWdata2,SWdata3);
%         % Smooth data
%         U1 = Inputs(2:end,:);
%         U1(3,:) = movmean(U1(3,:),24);
%         U1(5,:) = movmean(U1(5,:),24);
%         U1(7,:) = movmean(U1(7,:),24);
%         U1(9,:) = movmean(U1(9,:),24);
%         U1(11,:) = movmean(U1(11,:),3);
%         % Smooth data
%         U1(15,:) = movmean(U1(15,:),3);
%         U1(16,:) = movmean(U1(16,:),24);
%         U1(17,:) = movmean(U1(17,:),24);
%         U1(18,:) = movmean(U1(18,:),24);
%         U1(19,:) = movmean(U1(19,:),24);
%         U1(20,:) = movmean(U1(20,:),3);
%         U1(21,:) = movmean(U1(21,:),3);
%         U1(22,:) = movmean(U1(22,:),3);
%         U1(23,:) = movmean(U1(23,:),3);
%         Inputs(2:end,:) = U1;
    case 'NRLMSISE_1997_2008'
        [Inputs] = Comp_Inputs_NRLMSISE_1997_2008_betterSWdata(jd0,jdf,SWdata1,SWdata2);
    case 'TIEGCM_1997_2008_new'
        [Inputs] = Comp_Inputs_TIEGCM(jd0,jdf,SWdata1);
    otherwise
        warning('No valid DMDc model selected!')
end

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[romTimeOut,romStateOut]=ode113(@(t,x) ROM_MSISE_ODE(t,x,AC_pred,BC_pred,Inputs,jd0),[0 tf],romState_pred,opts);
romStateTime = [romStateOut, romTimeOut+et0]; % ROM modes and ephemeris time

end

