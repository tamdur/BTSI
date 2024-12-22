function [xOut,offset] = meanmin(dateIn,xIn,cycleInt)
%Subract away the mean of values for a specified solar cycle minimum
% INPUTS:
%           dateIn: Nx1 datetime vector for xIn
%           xIn: Nx1 vector of time series values passed in
%           cycleInt: Cycle minimum to select, passed in as integer
% OUTPUTS:
%           xOut: xIn, demeaned for cycleInt

load('lowspots_cutoff5_24_12_9.mat') %Load values from getlowspotindex


if length(dateIn) ~= length(lowAll)
    error(['Update the imported getlowspotindex results to be ' ...
            'consistent with input record'])
end
if cycleInt == 21
    offset = mean(xIn(low21,:));
elseif cycleInt == 22
    offset = mean(xIn(low22,:));
elseif cycleInt == 23
    offset = mean(xIn(low23,:));
elseif cycleInt == 24
    offset = mean(xIn(low24,:));
else
    error('Please specify an integer between 21 and 24 for cycleInt')
end
xOut = xIn - offset;
end

