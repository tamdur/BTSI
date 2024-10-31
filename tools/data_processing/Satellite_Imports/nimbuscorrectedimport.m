%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/TeddyA/Drive/Research/Data/Solar/instruments/PMOD_corrected/HF_corr_accessed_22_7_7.dat
%


%% Initialize variables.
filename = 'HF_corr_accessed_22_7_7.dat';
startRow = 8;

%% Format for each line of text:
%   columns are the following: year, month, day, julian day, number of
%   orbits with data, TSI
formatSpec = '%12d%10f%10f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
NIMBUSCORR.jd = dataArray{:, 1};
NIMBUSCORR.TSI = dataArray{:, 2};
NIMBUSCORR.TSIUNC = dataArray{:, 3};


%% Clear temporary variables
clearvars filename formatSpec fileID dataArray ans;