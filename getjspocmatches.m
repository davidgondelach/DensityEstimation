function [matches] = getjspocmatches()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% File processing

fid = fopen('/Users/davidgondelach/Dropbox (MIT)/Research-ROM/GPSdata/jspoc_matches.txt','r+');

% Get rid of initial lines
for i=1:9
    fgetl(fid);
end

matches.PlanetObjectID = [];
matches.NoradID = [];

n = 0;
tline = fgetl(fid);
while ischar(tline)
    n = n+1;
    matches(n).PlanetObjectID = tline(1:4);
    matches(n).NoradID = str2double(tline(23:27));
    tline = fgetl(fid);
end

end

