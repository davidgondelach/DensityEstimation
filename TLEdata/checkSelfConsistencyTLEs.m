function [objectDataSorted, covMEEerrors] = checkSelfConsistencyTLEs(objects, mu, jdate0TLEs, objectIDlabels)
%checkSelfConsistencyTLEs Compute and plot difference in position and in
%modified equinoctial elements taking from two different TLEs

noo = length(objects);

covMEEerrors = zeros(noo*6,1);
objectData = zeros(noo,9);

figure;
% Self consistency check
for i=1:noo
    clear diffR diffV epoch diffEpochInMin mee1 mee2 diffMee
    for j=length(objects(i).satrecs):-1:2
        index1 = j;
        index2 = j-1;
        
        diffEpochMinutes = (objects(i).satrecs(index1).jdsatepoch - objects(i).satrecs(index2).jdsatepoch) * 24 * 60;
        
        % Compute difference in position
        [~, rteme1 ,vteme1] = sgp4( objects(i).satrecs(index1), -diffEpochMinutes );
        [~, rteme2 ,vteme2] = sgp4( objects(i).satrecs(index2), 0 );
        diffR(j-1) = norm(rteme1 - rteme2);
        
        % Compute difference in MEE
        mee1 = pv2ep(rteme1,vteme1,mu)';
        mee2 = pv2ep(rteme2,vteme2,mu)';
        diffMee(j-1,:) = mee1 - mee2;
    end
    
    % Collect object data
    % [noradID, average position diff, # of TLEs,
    % a,e,i,RAAN,apogee,perigee]
    objectData(i,1) = objects(i).noradID;
    objectData(i,2) = mean(diffR);
    objectData(i,3) = length(objects(i).satrecs);
    objectData(i,4) = objects(i).satrecs(end).a;
    objectData(i,5) = objects(i).satrecs(end).ecco;
    objectData(i,6) = objects(i).satrecs(end).inclo;
    objectData(i,7) = objects(i).satrecs(end).nodeo;
    objectData(i,8) = objects(i).satrecs(end).alta;
    objectData(i,9) = objects(i).satrecs(end).altp;
    
    % Compute statistics of diff in MEE
    covMEEerrors((i-1)*6+1:i*6) = var(diffMee,0,1)';
    stdMEEerrors = std(diffMee,0,1)';

    % Plot difference in MEE
    diffMee(1:end,6) = wrapToPi(diffMee(1:end,6));
    diffMeeTimes = [objects(i).satrecs(1:end-1).jdsatepoch] - jdate0TLEs;
    
    subplot(2,3,1); plot(diffMeeTimes,diffMee(1:end,1)); hold on; %plot(diffMeeTimes([1,end]),3*stdMEEerrors(1)*ones(1,2),'--k'); plot(diffMeeTimes([1,end]),-3*stdMEEerrors(1)*ones(1,2),'--k');
    subplot(2,3,2); plot(diffMeeTimes,diffMee(1:end,2)); hold on; %plot(diffMeeTimes([1,end]),3*stdMEEerrors(2)*ones(1,2),'--k'); plot(diffMeeTimes([1,end]),-3*stdMEEerrors(2)*ones(1,2),'--k');
    subplot(2,3,3); plot(diffMeeTimes,diffMee(1:end,3)); hold on; %plot(diffMeeTimes([1,end]),3*stdMEEerrors(3)*ones(1,2),'--k'); plot(diffMeeTimes([1,end]),-3*stdMEEerrors(3)*ones(1,2),'--k');
    subplot(2,3,4); plot(diffMeeTimes,diffMee(1:end,4)); hold on; %plot(diffMeeTimes([1,end]),3*stdMEEerrors(4)*ones(1,2),'--k'); plot(diffMeeTimes([1,end]),-3*stdMEEerrors(4)*ones(1,2),'--k');
    subplot(2,3,5); plot(diffMeeTimes,diffMee(1:end,5)); hold on; %plot(diffMeeTimes([1,end]),3*stdMEEerrors(5)*ones(1,2),'--k'); plot(diffMeeTimes([1,end]),-3*stdMEEerrors(5)*ones(1,2),'--k');
    subplot(2,3,6); plot(diffMeeTimes,diffMee(1:end,6)); hold on; %plot(diffMeeTimes([1,end]),3*stdMEEerrors(6)*ones(1,2),'--k'); plot(diffMeeTimes([1,end]),-3*stdMEEerrors(6)*ones(1,2),'--k');
end
xlabeltext1 = 'Time [days]';
subplot(2,3,1);  xlabel(xlabeltext1); ylabel('\Delta p [km]');
subplot(2,3,2);  xlabel(xlabeltext1); ylabel('\Delta f [-]');
subplot(2,3,3);  xlabel(xlabeltext1); ylabel('\Delta g [-]'); legend(objectIDlabels,'Location','northeast');
subplot(2,3,4);  xlabel(xlabeltext1); ylabel('\Delta h [-]');
subplot(2,3,5);  xlabel(xlabeltext1); ylabel('\Delta k [-]');
subplot(2,3,6);  xlabel(xlabeltext1); ylabel('\Delta L [rad]');

objectDataSorted = sortrows(objectData,2);

end

