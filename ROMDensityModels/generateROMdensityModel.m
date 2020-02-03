function [AC,BC,Uh,F_U,Dens_Mean,M_U,SLTm,LATm,ALTm,maxAtmAlt,SWinputs,Qrom] = generateROMdensityModel(ROMmodel,r,jd0,jdf)
%generateROMdensityModel - Generate reduced-order density model
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

switch ROMmodel
    case 'JB2008_1999_2010'
        TA = load('JB2008_1999_2010_ROM_r100.mat');
        
        % Compute reduced-order dynamic density model:
        % PhiC contains continuous-time dynamic and input matrices
        % Uh contains the POD spatial modes
        % Qrom is the covariance matrix of ROM prediction error
        [PhiC,Uh,Qrom] = generateROM_JB2008(TA,r);
        
        % Compute the space weather inputs in the estimation period
        [eopdata,SOLdata,DTCdata] = loadJB2008SWdata(); % Read JB2008 space weather inputs from file
        [SWinputs] = computeSWinputs_JB2008(jd0,jdf+1,eopdata,SOLdata,DTCdata);
        
        % Maximum altitude of ROM-JB2008 density model (for higher altitudes
        % density is set to zero)
        maxAtmAlt = 800;
        
    case 'TIEGCM_1997_2008'
        TA = load('TIEGCM_1997_2008_ROM_r100.mat');
        
        % Compute reduced-order dynamic density model:
        % PhiC contains continuous-time dynamic and input matrices
        % Uh contains the POD spatial modes
        % Qrom is the covariance matrix of ROM prediction error
        [PhiC,Uh,Qrom] = generateROM_TIEGCM(TA,r);
        
        % Compute the space weather inputs in the estimation period
        SWpath = fullfile('Data','SW-All.txt');
        [ SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM ] = inputSWtiegcm( SWpath ); % Read F10.7 and Kp values from file
        [SWinputs] = computeSWinputs_TIEGCM(jd0,jdf+1,SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM);
        
        % Maximum altitude of ROM-TIEGCM density model (for higher altitudes
        % density is set to zero)
        maxAtmAlt = 500; % Density is computed up to 500km through extrapolation (ROM-TIEGCM model itself goes up to 450km)
        
    case 'NRLMSISE_1997_2008'
        TA = load('NRLMSISE_1997_2008_ROM_r100.mat');
        
        % Compute reduced-order dynamic density model:
        % PhiC contains continuous-time dynamic and input matrices
        % Uh contains the POD spatial modes
        % Qrom is the covariance matrix of ROM prediction error
        [PhiC,Uh,Qrom] = generateROM_NRLMSISE(TA,r);
        
        % Compute the space weather inputs in the estimation period
        SWpath = fullfile('Data','SW-All.txt');
        [ SWmatDaily, SWmatMonthlyPred ] = inputSWnrlmsise( SWpath ); % Read F10.7 and ap values from file
        [SWinputs] = computeSWinputs_NRLMSISE(jd0,jdf+1,SWmatDaily,SWmatMonthlyPred);
        
        % Maximum altitude of ROM-NRLMSISE density model (for higher altitudes
        % density is set to zero)
        maxAtmAlt = 800;
        
    otherwise
        warning('No valid ROM model selected!')
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
