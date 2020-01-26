function [ EOPMat ] = inputEOP_Celestrak( EOPfilename )
%INPUTEOP_Celestrak - Read Earth Orientation Parameters from text file. For
% reading: https://www.celestrak.com/SpaceData/EOP-All.txt
% 
%   Inputs:
%     EOPfilename
%
%   Outputs:
%     xp:    [arcsec]
%     yp:    [arcsec]
%     dut1:  [s]
%     lod:   [s]
%     ddpsi: [rad]
%     ddeps: [rad]
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

% Initialize output
EOPMat = zeros(nofObs,6);

%% Read Earth Orientation Parameters
for ind=1:nofObs
    
    str = fgetl(fid);
    
    % Read PM-x
    EOPMat(ind,1) = str2double(str(18:26)); % arcsec
    
    % Read PM-y
    EOPMat(ind,2) = str2double(str(28:36)); % arcsec
    
    % Read UT1-UTC
    EOPMat(ind,3) = str2double(str(38:47));
    
    % Read length of day
    EOPMat(ind,4) = str2double(str(49:58));
    
    % Read dPsi
    EOPMat(ind,5) = str2double(str(60:68));  % arcsec
    
    % Read dEps
    EOPMat(ind,6) = str2double(str(70:78)); % arcsec
    
end

EOPMat(isnan(EOPMat)) = 0;
EOPMat(:,1:2) = EOPMat(:,1:2)/3600*pi/180; % rad
EOPMat(:,5:6) = EOPMat(:,5:6)/3600*pi/180; % rad

%% Close file
fclose(fid);

end
