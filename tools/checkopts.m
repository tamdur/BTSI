function optsOut = checkopts(optsIn)
%CHECKOPTS Create default settings for all Gibbs sampling options that have
%not been explicitly specified by the user


if ~isfield(optsIn,'randomizeChain') %Load default observation array, otherwise load provided one
    optsIn.randomizeChain=false;
end
if ~isfield(optsIn,'excludeFliers')
    optsIn.excludeFliers=true;%true to remove outlier observations from examined dataset
end
if ~isfield(optsIn,'reps')
    optsIn.reps=10500; %Total steps of Gibbs Sampler
end
if ~isfield(optsIn,'burn')
    optsIn.burn=500; %Steps excluded in burnin period
end
if ~isfield(optsIn,'logContributions')
    optsIn.logContributions=false; %Set to true to include record of innovation contributions
end
if ~isfield(optsIn,'dispProgress')
    optsIn.dispProgress=false; %Set to true to display progress of every 100 iterations
end
if ~isfield(optsIn,'lags')
    optsIn.lags=2; %Set process to AR(2)
end
if ~isfield(optsIn,'normalize')
    optsIn.normalize=true; %Set to true to normalize data within. For getpriors call
end
if ~isfield(optsIn,'magDependent')
    optsIn.magDependent=true; %Set to true to normalize data within. For getpriors call
end
if ~isfield(optsIn,'satOnly')
    optsIn.satOnly=false; %Set to true to not use any proxy data
end
if ~isfield(optsIn,'proxyModel')
    optsIn.proxyModel=false; %Set to true to calculate BTSI using only SORCE, proxies
end
if ~isfield(optsIn,'satOnlyDrift')
    optsIn.satOnlyDrift=false; %true to turn on satellite-only model with drift turned on
end
if ~isfield(optsIn,'virgo')
    optsIn.virgo=false; %true to set VIRGO to reference level satellite
end
optsIn.cmpStYr=1978;
optsOut=optsIn;
end

