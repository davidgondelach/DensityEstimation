function [TLEstruct] = getTLEs(filename)
%getTLEs - Read TLE data from file and sort on data
%
% This code is licensed under the GNU General Public License version 3.
%
% Based on code by Aleksander Lidtke, University of Southampton, UK, July 2015
%
% Modified by: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

global whichconst; % The gravity constant that we're using.

finput = fullfile(filename); % Open the file for this object. Use OS independent path concatenation.
fid = fopen(finput);

%% Check the number of TLEs in the file.
fseek(fid, 0, 'eof');
fileSize = ftell(fid);
frewind(fid); % Put the file pointer at the beginning so that we read the whole file. 
data = fread(fid, fileSize, 'uint8');
frewind(fid);
% Count number of line-feeds and increase by one for EOF.
numLines = sum(data == 10) + 1;

if mod(numLines,2)==1
    numLines = numLines-1;
end

%% Parse the TLEs from the file by creating satrecs first.
k = 1; % Counter of TLEs read.

longstr1 = fgetl(fid); % Read the first line

while ischar(longstr1) && numel(longstr1)~=0 % Read all the TLEs for this object. Last line in the file is end of line char, don't want to try to convert it to a TLE.
    longstr2 = fgetl(fid); % Read the second line

    [satrec] = twoline2rv_edit( whichconst, longstr1, longstr2); % Initialise the satrec.
    satrecs(k) = satrec; % Parse this TLE by creating the satrec.
    EpochsJD(k) = satrec.jdsatepoch; % Need the epochs of the TLEs to sort them before appending to the final struct.
    noradID(k) = satrec.satnum;

    longstr1 = fgetl(fid); % Proceed to the next TLE.
    k=k+1;
end

fclose(fid); % Done with the file.


noradIDs = unique(noradID);
for i=1:length(noradIDs)
    object = noradIDs(i);
    objectIndeces = find(noradID==object);
    TLEstruct(i).noradID = object;
    
    objectEpochsJD = EpochsJD(objectIndeces);
    if ~issorted(objectEpochsJD)
        [~,SortedObjectEpochsJDIndices] = sort(objectEpochsJD);
        for j=1:length(objectEpochsJD)
            TLEstruct(i).satrecs(j) = satrecs(objectIndeces(SortedObjectEpochsJDIndices(j)));
        end
    else
        TLEstruct(i).satrecs = satrecs(objectIndeces);
    end
    
end

end

%------------- END OF CODE --------------
