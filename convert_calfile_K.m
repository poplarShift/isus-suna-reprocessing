%% Import data from text file.
% Auto-generated by MATLAB on 2015/06/13 10:20:57

%% Initialize variables.
filename = 'ISUS132K.CAL';
delimiter = ',';
startRow = 26;
endRow = 281;

%% Format string for each line of text:
%   column2: double (%f)
%	column3: double (%f)
%   column4: double (%f)
%	column5: double (%f)
%   column6: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
lambda = dataArray{:, 1};
ENO3 = dataArray{:, 2};
ESW = dataArray{:, 3};
T_ESW = dataArray{:, 4};
ref_DIW = dataArray{:, 5};


%% Clear temporary variables
clearvars filename delimiter startRow endRow formatSpec fileID dataArray H ans;

save('ISUS132K_CAL.mat')
