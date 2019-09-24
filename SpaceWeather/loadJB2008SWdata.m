function [eopdata,SOLdata,DTCdata] = loadJB2008SWdata()

global const
SAT_Const
constants
% DE430Coeff = load('DE430Coeff.mat');
% PC = DE430Coeff.DE430Coeff;

% read Earth orientation parameters
fid = fopen('Data/eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

% read space weather data
fid = fopen('Data/SOLFSMY.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% READ GEOMAGNETIC STORM DTC VALUE
fid = fopen('Data/DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);

end

