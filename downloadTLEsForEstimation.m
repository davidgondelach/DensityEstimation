function [objects] = downloadTLEsForEstimation(startYear, startMonth, startDay, endYear, endMonth, endDay, maxAlt, selectedObjects)

curlCmd = "curl -o TLEdata/estimationObjects.tle https://www.space-track.org/ajaxauth/login -d 'identity=davidgondelach@gmail.com&password=VanWanroystraat8&query=https://www.space-track.org/basicspacedata/query/class/tle/";
epochQuery = strcat("EPOCH/",int2str(startYear), '-', sprintf('%02d',startMonth), '-', sprintf('%02d',startDay), '--', sprintf('%02d',endYear), '-', sprintf('%02d',endMonth), '-', sprintf('%02d',endDay), '/');
apogeeQuery = strcat("APOGEE/%3C",num2str(maxAlt),'/');

if ~isempty(selectedObjects)
    objectsQuery = strcat("NORAD_CAT_ID/",sprintf('%05d',selectedObjects(1)));
    for i=2:length(selectedObjects)
        objectsQuery = strcat(objectsQuery,",%20",sprintf('%05d',selectedObjects(i)));
    end
    objectsQuery = strcat(objectsQuery,"/");
else
    objectsQuery = "";
end
orderingQuery = "orderby/NORAD_CAT_ID%20asc/format/tle/'";
TLEdwnlcmd = strcat(curlCmd, epochQuery, apogeeQuery, objectsQuery, orderingQuery);
system(TLEdwnlcmd);

filename = fullfile('TLEdata','estimationObjects.tle');
[objects] = getTLEs(filename);

end

