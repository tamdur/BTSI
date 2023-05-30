%Wrapper to run a Gibbs sampling of BTSI
%
% Ted Amdur
% 2023/03/04

clearvars
mainExperiment=0; %Run full experiment that results in published results and figures
pmodCorrections=0; %Run with Frohlich corrections
NRLTSIcomp=1; %Run with just spots, faculae, SORCE
virgo=0; %Run with VIRGO as reference satellite
satOnly=0; %Run satellite only with no drift
satOnlyDrift=0; %Run satellite only WITH drift

if mainExperiment
    opts.burn = 1000; %Number of burn-in reps assumed for chain length analysis
    opts.reps=11000; %Total length of chain, including burn-in
    opts.excludeFliers=true;%1 to remove outlier observations from examined dataset
    opts.satOnly=false;%1 to ignore proxy data and only use satellites (no drift calculation), 0 otherwise
    opts.proxyModel=0;%1 to use same datasets as NRLTSI2, 0 otherwise
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=true;
    opts.normalize=true;
    opts.magDependent=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='chain_output/ar2_23_05_11_long.mat';
    runchain_22_04_25([],[],[],[],opts);
end
if pmodCorrections
    opts.burn = 1000; %Number of burn-in reps assumed for chain length analysis
    opts.reps=2000; %Total length of chain, including burn-in
    opts.excludeFliers=true;%1 to remove outlier observations from examined dataset
    opts.satOnly=false;%1 to ignore proxy data and only use satellites (no drift calculation), 0 otherwise
    opts.proxyModel=0;%1 to use same datasets as NRLTSI2, 0 otherwise
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=false;
    opts.normalize=true;
    opts.magDependent=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='chain_output/ar2_23_05_10_pmodcorrections.mat';
    opts.obsmatrix='obs_23_05_10_pmod.mat';
    load(opts.obsmatrix)
    runchain_22_04_25(valM,oM,dateM,colLabels,opts);
end
if NRLTSIcomp
    opts.burn = 1000; %Number of burn-in reps assumed for chain length analysis
    opts.reps=2000; %Total length of chain, including burn-in
    opts.excludeFliers=true;%1 to remove outlier observations from examined dataset
    opts.satOnly=false;%1 to ignore proxy data and only use satellites (no drift calculation), 0 otherwise
    opts.proxyModel=1;%1 to use same datasets as NRLTSI2, 0 otherwise
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=true;
    opts.logContributions=false;
    opts.normalize=true;
    opts.magDependent=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='chain_output/ar2_23_05_10_nrltsicomp.mat';
    runchain_22_04_25([],[],[],[],opts);
end
if virgo
    opts.burn = 1000; %Number of burn-in reps assumed for chain length analysis
    opts.reps=2000; %Total length of chain, including burn-in
    opts.excludeFliers=true;%1 to remove outlier observations from examined dataset
    opts.satOnly=false;%1 to ignore proxy data and only use satellites (no drift calculation), 0 otherwise
    opts.proxyModel=0;%1 to use same datasets as NRLTSI2, 0 otherwise
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=false;
    opts.normalize=true;
    opts.magDependent=true;
    opts.virgo=true; %Use VIRGO baseline
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='chain_output/ar2_23_05_10_virgo.mat';
    runchain_22_04_25([],[],[],[],opts);
end
if satOnly
    opts.burn = 1000; %Number of burn-in reps assumed for chain length analysis
    opts.reps=2000; %Total length of chain, including burn-in
    opts.excludeFliers=true;%1 to remove outlier observations from examined dataset
    opts.satOnly=false;%1 to ignore proxy data and only use satellites (no drift calculation), 0 otherwise
    opts.proxyModel=0;%1 to use same datasets as NRLTSI2, 0 otherwise
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=false;
    opts.normalize=true;
    opts.magDependent=true;
    opts.satOnly=true; %Use VIRGO baseline
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='chain_output/ar2_23_05_10_satonly.mat';
    runchain_22_04_25([],[],[],[],opts);
end
if satOnlyDrift
    opts.burn = 1000; %Number of burn-in reps assumed for chain length analysis
    opts.reps=2000; %Total length of chain, including burn-in
    opts.excludeFliers=true;%1 to remove outlier observations from examined dataset
    opts.satOnly=false;%1 to ignore proxy data and only use satellites (no drift calculation), 0 otherwise
    opts.proxyModel=0;%1 to use same datasets as NRLTSI2, 0 otherwise
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=false;
    opts.normalize=true;
    opts.magDependent=true;
    opts.satOnlyDrift=true; %Use VIRGO baseline
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='chain_output/ar2_23_05_10_satonlydrift.mat';
    runchain_22_04_25([],[],[],[],opts);
end
