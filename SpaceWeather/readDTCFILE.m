function DTCFILE = readDTCFILE(filename, startRow, endRow)
%READDTCFILE Import numeric data from a DTCFILE text file as a matrix.
%   DTCFILE = readDTCFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   DTCFILE = readDTCFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   DTCFILE = readDTCFILE('DTCFILE.txt', 1, 7912);
%

%% Initialize variables.
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format for each line of text:
%   column2: YEAR: double (%f)
%	column3: DOY: double (%f)
%   column4: DTC1: double (%f)
%	column5: DTC2: double (%f)
%   column6: DTC3: double (%f)
%	column7: DTC4: double (%f)
%   column8: DTC5: double (%f)
%	column9: DTC6: double (%f)
%   column10: DTC7: double (%f)
%	column11: DTC8: double (%f)
%   column12: DTC9: double (%f)
%	column13: DTC10: double (%f)
%   column14: DTC11: double (%f)
%	column15: DTC12: double (%f)
%   column16: DTC13: double (%f)
%	column17: DTC14: double (%f)
%   column18: DTC15: double (%f)
%	column19: DTC16: double (%f)
%   column20: DTC17: double (%f)
%	column21: DTC18: double (%f)
%   column22: DTC19: double (%f)
%	column23: DTC20: double (%f)
%   column24: DTC21: double (%f)
%	column25: DTC22: double (%f)
%   column26: DTC23: double (%f)
%	column27: DTC24: double (%f)
formatSpec = '%*3s%5f%4f%5f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%4f%f%[^\n\r]';

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
DTCFILE = [dataArray{1:end-1}];
