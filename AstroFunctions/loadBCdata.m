function [BCdata] = loadBCdata( BCfilePath )
% Read ballistic coefficient data from file

% Open ballistic coefficient data file
fid = fopen(BCfilePath,'r');

% Skip first two lines
for i=1:2
    fgetl(fid);
end

% Read data
%  ----------------
% |  ID      BC    |
% |        m^2/kg  |
%  ----------------
BCdata = fscanf(fid,'%d %f',[2 inf]);
fclose(fid);

BCdata = BCdata';

end

