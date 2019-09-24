function [ EOPMat ] = inputEOP_Celestrak( EOPfilename )
%READEOP Reads Earth Orientation Parameters
%   Inptus:
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
%     

%% Output initialization

% Read number of lines
nLines = countLines(EOPfilename);

EOPMat = zeros(nLines,6);

%% File processing

% Open file
fid = fopen(EOPfilename, 'r');

% Get rid of initial lines
for i=1:34
    fgetl(fid);
end
str = fgetl(fid);
n_daily_obs = str2double(str(21:25));
str = fgetl(fid);

for ind=1:n_daily_obs
    
    str = fgetl(fid);
    
    % Read PM-x
    EOPMat(ind,1) = str2double(str(18:26)); % arcsec
%     if isnan(EOPMat(ind,1))
%         EOPMat(ind,1) = 0;
%     end
    
    % Read PM-y
    EOPMat(ind,2) = str2double(str(28:36)); % arcsec
%     if isnan(EOPMat(ind,2))
%         EOPMat(ind,2) = 0;
%     end
    
    % Read UT1-UTC
    EOPMat(ind,3) = str2double(str(38:47));
%     if isnan(EOPMat(ind,3))
%         EOPMat(ind,3) = 0;
%     end
    
    % Read length of day
    EOPMat(ind,4) = str2double(str(49:58));
%     if isnan(EOPMat(ind,4))
%         EOPMat(ind,4) = 0;
%     end
    
    % Read dPsi
    EOPMat(ind,5) = str2double(str(60:68));  % arcsec
%     if isnan(EOPMat(ind,5))
%         EOPMat(ind,5) = 0;
%     end
    
    % Read dEps
    EOPMat(ind,6) = str2double(str(70:78)); % arcsec
%     if isnan(EOPMat(ind,6))
%         EOPMat(ind,6) = 0;
%     end
    
end

EOPMat(isnan(EOPMat)) = 0;
EOPMat(:,1:2) = EOPMat(:,1:2)/3600*pi/180; % rad
EOPMat(:,5:6) = EOPMat(:,5:6)/3600*pi/180; % rad

%% Close file
fclose(fid);

end


function count = countLines(fname)

fh = fopen(fname, 'rt');
assert(fh ~= -1, 'Could not read: %s', fname);
x = onCleanup(@() fclose(fh));
count = 0;
while ~feof(fh)
    count = count + sum( fread( fh, 16384, 'char' ) == char(10) );
end

end
