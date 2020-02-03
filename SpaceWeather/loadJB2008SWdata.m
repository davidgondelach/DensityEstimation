function [eopdata,SOLdata,DTCdata] = loadJB2008SWdata()
% LOADJB2008SWDATA - Read space weather data and Earth orientation 
% parameters needed for JB2008 density model from file
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

global const

SAT_Const;
constants;

% Read Earth orientation parameters
[ eopdata ] = inputEOP_Celestrak_Full( 'Data/EOP-All.txt' );

% Read space weather data: solar activity indices
SOLdata = readSOLFSMY('Data/SOLFSMY.txt')';

% Read geomagnetic storm DTC values
DTCdata = readDTCFILE('Data/DTCFILE.txt')';

end

%------------- END OF CODE --------------
