function [AC,BC,Uh,F_U,Dens_Mean,M_U,SLTm,LATm,ALTm,maxAtmAlt,SWinputs,Qrom] = generateROMdensityModel(DMDmodel,r,jd0,jdf)
%generateROMdensityModel - Generate reduced-order density model

% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 24-Sep-2019

%------------- BEGIN CODE --------------

switch DMDmodel
    case 'JB2008_1999_2010'
        TA = load('JB2008_1999_2010_ROM_r100.mat');
        
        % Compute reduced-order dynamic density model:
        % PhiC contains continuous-time dynamic and input matrices
        % Uh contains the POD spatial modes
        % Qrom is the covariance matrix of ROM prediction error
        [PhiC,Uh,Qrom] = generateROM_JB2008(TA,r);
        
        % Compute the space weather inputs in the estimation period
        [eopdata,SOLdata,DTCdata] = loadJB2008SWdata();
        [SWinputs] = computeSWinputs_JB2008(jd0,jdf+20,eopdata,SOLdata,DTCdata);
        
        % Maximum altitude of ROM density model
        maxAtmAlt = 800;
        
    case 'TIEGCM_1997_2008'
        TA = load('TIEGCM_1997_2008_ROM_r100.mat');
        
        % Compute reduced-order dynamic density model:
        % PhiC contains continuous-time dynamic and input matrices
        % Uh contains the POD spatial modes
        % Qrom is the covariance matrix of ROM prediction error
        [PhiC,Uh,Qrom] = generateROM_TIEGCM(TA,r);
        
        % Compute the space weather inputs in the estimation period
        TIEGCM_SWdata = TA.SWdataFull;
        [SWinputs] = computeSWinputs_TIEGCM(jd0,jdf,TIEGCM_SWdata);
        
        % Maximum altitude of ROM density model
        maxAtmAlt = 500;
        
    case 'NRLMSISE_1997_2008'
        TA = load('NRLMSISE_1997_2008_ROM_r100.mat');
        
        % Compute reduced-order dynamic density model:
        % PhiC contains continuous-time dynamic and input matrices
        % Uh contains the POD spatial modes
        % Qrom is the covariance matrix of ROM prediction error
        [PhiC,Uh,Qrom] = generateROM_NRLMSISE(TA,r);
        
        % Compute the space weather inputs in the estimation period
        SWpath = fullfile('Data','SW-All.txt');
        [ SWmatDaily, SWmatMonthlyPred ] = inputSWnrlmsise( SWpath );
        [SWinputs] = computeSWinputs_NRLMSISE(jd0,jdf,SWmatDaily,SWmatMonthlyPred);
        
        % Maximum altitude of ROM density model
        maxAtmAlt = 800;
        
    otherwise
        warning('No valid DMDc model selected!')
end
        
% Setup of ROM Modal Interpolation
sltm = TA.localSolarTimes;
latm = TA.latitudes;
altm = TA.altitudes;
n_slt = length(sltm);
n_lat = length(latm);
n_alt = length(altm);

% Mean density
Dens_Mean = TA.densityDataMeanLog;
        
% Generate full 3D grid in local solar time, latitude and altitude
[SLTm,LATm,ALTm]=ndgrid(sltm,latm,altm);

% Generate interpolant for each POD spatial mode in Uh
F_U{r} = [];
for i = 1:r
    Uhr = reshape(Uh(:,i),n_slt,n_lat,n_alt); % i-th left singular vector on grid
    F_U{i} = griddedInterpolant(SLTm,LATm,ALTm,Uhr,'linear','linear'); % Create interpolant of Uhr
end
% Generate interpolant for the mean density in DenS_Mean
Mr = reshape(Dens_Mean,n_slt,n_lat,n_alt);
M_U = griddedInterpolant(SLTm,LATm,ALTm,Mr,'linear','linear');
% Compute dynamic and input matrices
AC = PhiC(1:r,1:r)/3600;
BC = PhiC(1:r,r+1:end)/3600;

end

%------------- END OF CODE --------------
