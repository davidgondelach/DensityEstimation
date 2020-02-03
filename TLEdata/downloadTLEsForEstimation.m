function [objects] = downloadTLEsForEstimation(username, password, startYear, startMonth, startDay, endYear, endMonth, endDay, maxAlt, selectedObjects)
%downloadTLEsForEstimation - Download TLE data
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

curlCmd = 'curl -o TLEdata/estimationObjects.tle https://www.space-track.org/ajaxauth/login -d ';
usernameCmd = strcat('"identity=',username);
passwordCmd = strcat('&password=',password);
tleQuery = '&query=https://www.space-track.org/basicspacedata/query/class/tle/';
epochQuery = strcat('EPOCH/',int2str(startYear), '-', sprintf('%02d',startMonth), '-', sprintf('%02d',startDay), '--', sprintf('%02d',endYear), '-', sprintf('%02d',endMonth), '-', sprintf('%02d',endDay), '/');
apogeeQuery = strcat('APOGEE/%3C',num2str(maxAlt),'/');

if ~isempty(selectedObjects)
    objectsQuery = strcat('NORAD_CAT_ID/',sprintf('%05d',selectedObjects(1)));
    for i=2:length(selectedObjects)
        objectsQuery = strcat(objectsQuery,',%20',sprintf('%05d',selectedObjects(i)));
    end
    objectsQuery = strcat(objectsQuery,'/');
else
    objectsQuery = '';
end
orderingQuery = 'orderby/NORAD_CAT_ID%20asc/format/tle/"';
TLEdwnlcmd = strcat(curlCmd, usernameCmd, passwordCmd, tleQuery, epochQuery, apogeeQuery, objectsQuery, orderingQuery);
system(TLEdwnlcmd);

filename = fullfile('TLEdata','estimationObjects.tle');
[objects] = getTLEs(filename);

end

%------------- END OF CODE --------------
