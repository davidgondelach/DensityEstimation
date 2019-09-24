function [objects] = getTLEsForEstimation(startYear, startMonth, startDay, endYear, endMonth, endDay, selectedObjects)

relativeDir = 'TLEdata';

selectedObjects = sort(selectedObjects);

objects = struct('noradID',{},'satrecs',{});
for i=1:length(selectedObjects)
    noradIDstr = sprintf('%05d',selectedObjects(i));
    try % Try to read a ready TLEsStruct from the hard drive to avoid time-consmuing parsing.
        matFile = fullfile(relativeDir,strcat(noradIDstr,'.mat'));
        matData = load( matFile ); % We already have a TLEsStruct saved there.
        object = matData.object;
    catch
        filepath = fullfile(relativeDir,strcat(noradIDstr,'.tle'));
        [object] = getTLEs(filepath);
        save(matFile,'object');
    end
    
    jdStart = juliandate(datetime(startYear, startMonth, startDay,0,0,0));
    jdEnd   = juliandate(datetime(endYear, endMonth, endDay,0,0,0));
    firstTLE = find([object.satrecs.jdsatepoch]>=jdStart,1,'first');
    lastTLE  = find([object.satrecs.jdsatepoch]<=jdEnd,1,'last');
    object.satrecs = object.satrecs(firstTLE:lastTLE);
    objects = [objects,object];
end

end

