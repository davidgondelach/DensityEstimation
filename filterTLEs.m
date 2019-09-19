function [objects] = filterTLEs(objects)
%filterTLEs Filter out corrected TLEs and remove objects with increasing
%semi-major axis

for i=length(objects):-1:1
    
    % Remove corrected TLEs
    isCorrectedTLE = [diff([objects(i).satrecs.jdsatepoch]) == 0, false];
    if any(isCorrectedTLE)
        objects(i).satrecs(isCorrectedTLE) = [];
    end
    
    % Discard objects with less than 5 TLEs, more than 2-days TLE
    % separation and increasing semi-major axis
    if length(objects(i).satrecs) < 5 ... % At least 5 TLE
    || any(diff([objects(i).satrecs.jdsatepoch]) > 2) ... % Less than 2 days separation
    || objects(i).satrecs(end).no < objects(i).satrecs(1).no % No increase in semi-major axis, i.e. no maneuvers or strong SRP
        objects(i) = [];
        continue;
    end
    
    % Discard maneuvering objects (i.e. with increasing semi-major axis)
    if any(diff([objects(i).satrecs.no]) < 0)
        negNindeces = 2:length(objects(i).satrecs);
        negNindeces = negNindeces(diff([objects(i).satrecs.no]) < 0);
        for j=1:length(negNindeces)
            negNindex = negNindeces(j);
            if negNindex >=4 && negNindex <= length(objects(i).satrecs)-2 
                if objects(i).satrecs(negNindex+1).no < objects(i).satrecs(negNindex-2).no
                    if objects(i).satrecs(negNindex+2).no < objects(i).satrecs(negNindex-3).no
                        figure;
                        subplot(2,1,1);plot([objects(i).satrecs.no])
                        subplot(2,1,2);plot([objects(i).satrecs.a])
                        objects(i) = [];
                        break;
                    end
                end
            end
        end     
    end
end

end

