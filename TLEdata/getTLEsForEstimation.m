function [objects] = getTLEsForEstimation(startYear, startMonth, startDay, endYear, endMonth, endDay, selectedObjects, getTLEsFromSingleFile)
%getTLEsForEstimation - Read TLE data from file
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

relativeDir = 'TLEdata';

selectedObjects = sort(selectedObjects);

if getTLEsFromSingleFile
    % Read TLEs from "estimationObjects.tle"
    filename = fullfile(relativeDir,'estimationObjects.tle');
    if ~isfile(filename)
        error('TLE file estimationObjects.tle was not found in %s directory.', relativeDir);
    end
    [objects] = getTLEs(filename);
else
    % Read TLEs from individual files named "[NORADID].tle"
    objects = struct('noradID',{},'satrecs',{});
    for i=1:length(selectedObjects)
        noradIDstr = sprintf('%05d',selectedObjects(i));
        filepath = fullfile(relativeDir,strcat(noradIDstr,'.tle'));
        if ~isfile(filepath)
            error('TLE file %s.tle was not found in %s directory.', noradIDstr, relativeDir);
        end
        [object] = getTLEs(filepath);
        objects = [objects,object];
    end
end

% Filter TLEs on date
jdStart = juliandate(datetime(startYear, startMonth, startDay,0,0,0));
jdEnd   = juliandate(datetime(endYear, endMonth, endDay,0,0,0));
for i=1:length(selectedObjects)
    % NORAD ID
    ID = selectedObjects(i);
    % Find object
    index = find([objects.noradID]==ID);
    if isempty(index)
        error('No TLEs found for object %.0f.', ID);
    end
    % Get object's TLE data
    object = objects(index);
    
    % Filter TLEs on date
    firstTLE = find([object.satrecs.jdsatepoch]>=jdStart,1,'first');
    lastTLE  = find([object.satrecs.jdsatepoch]<=jdEnd,1,'last');
    if isempty(firstTLE) || isempty(lastTLE)
        error('No TLEs found for object %.0f between %.0f-%.0f-%.0f and %.0f-%.0f-%.0f.', ID, startDay, startMonth, startYear, endDay, endMonth, endYear);
    end
    objects(index).satrecs = object.satrecs(firstTLE:lastTLE);
end

end

%------------- END OF CODE --------------
