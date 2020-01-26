function SOLFSMY = readSOLFSMY(filename, startRow, endRow)
%READSOLFSMY Import numeric data from a SOLFSMY text file as a matrix.
%   SOLFSMY = readSOLFSMY(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   SOLFSMY = readSOLFSMY(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   SOLFSMY = readSOLFSMY('SOLFSMY.txt', 1, 7913);
%


%% Initialize variables.
if nargin<=2
    startRow = 5;
    endRow = inf;
end

%% Format for each line of text:
%   column1: YEAR: double (%f)
%	column2: DOY: double (%f)
%   column3: JulianDay: double (%f)
%	column4: F10: double (%f)
%   column5: F81c: double (%f)
%	column6: S10: double (%f)
%   column7: S81c: double (%f)
%	column8: M10: double (%f)
%   column9: M81c: double (%f)
%	column10: Y10: double (%f)
%   column11: Y81c: double (%f)
formatSpec = '%6f%4f%12f%6f%6f%6f%6f%6f%6f%6f%6f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Create output variable
SOLFSMY = [dataArray{1:end-1}];
