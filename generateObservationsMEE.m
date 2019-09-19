function [meeObs] = generateObservationsMEE(objects,obsEpochs,GM_kms)
%generateObservationsMEE Generate observations in MEE at specified
%observation epochs

nofObjects = length(objects);
nofObs = length(obsEpochs);
meeObs = zeros(6*nofObjects,nofObs);
for i=1:nofObjects
    for j=1:nofObs
        % Observation epoch
        obsEpoch = obsEpochs(j);
        % Find nearest newer TLE
        satrecIndex = find([objects(i).satrecs.jdsatepoch]>=obsEpoch,1,'first');
        diffObsTLEEpochMinutes = (obsEpoch - objects(i).satrecs(satrecIndex).jdsatepoch) * 24*60;
        % Compute SGP4 state at epoch
        [~, rtemeObs ,vtemeObs] = sgp4( objects(i).satrecs(satrecIndex), diffObsTLEEpochMinutes );
        % Convert to J2000
        [rj2000, vj2000] = convertTEMEtoJ2000(rtemeObs', vtemeObs', obsEpoch);
        meeObs(6*(i-1)+1:6*i,j) = pv2ep(rj2000,vj2000,GM_kms)';
    end
end

end