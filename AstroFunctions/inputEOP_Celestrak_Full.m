function [ EOPMat ] = inputEOP_Celestrak_Full( EOPfilename )
%INPUTEOP_Celestrak - Read Earth Orientation Parameters from text file. For
% reading: https://www.celestrak.com/SpaceData/EOP-All.txt
% 
%   Inputs:
%     EOPfilename
%
%   Outputs:
% 	  Year
%     Month (01-12)
% 	  Day
% 	  Modified Julian Date (Julian Date at 0h UT minus 2400000.5)
% 	  x (arc seconds)
% 	  y (arc seconds)
% 	  UT1-UTC (seconds)
% 	  Length of Day (seconds)
% 	  ??? (arc seconds)
% 	  ??? (arc seconds)
% 	  ?X (arc seconds)
% 	  ?Y (arc seconds)
% 	  Delta Atomic Time, TAI-UTC (seconds)
%     

%% File processing

% Open file
fid = fopen(EOPfilename, 'r');

% Skip initial lines
for i=1:34
    fgetl(fid);
end

% Read number of observed points
str = fgetl(fid);
nofObs = str2double(str(21:25));

fgetl(fid); % Go to next line

%% Read Earth Orientation Parameters
EOPMat = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 nofObs]);

%% Close file
fclose(fid);

end
