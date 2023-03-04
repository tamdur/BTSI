%Wrapper to run a Gibbs sampling of BTSI
%
% Ted Amdur
% 2023/03/04

mainExperiment=1; %Run full experiment that results in published results and figures


if mainExperiment
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=1500; %Total length of chain, including burn-in
    opts.excludeFliers=1;%1 to remove outlier observations from examined dataset
    opts.satOnly=0;%1 to ignore proxy data and only use satellites (no drift calculation), 0 otherwise
    opts.proxyModel=0;%1 to use same datasets as NRLTSI2, 0 otherwise
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.excludeFliers=false;
    opts.logContributions=true;
    opts.normalize=false;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='chain_output/ar2_23_03_04.mat';
    runchain_22_04_25([],[],[],opts);
end