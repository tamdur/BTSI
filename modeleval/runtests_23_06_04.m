%Run tests of observer predictions in parallel computing environment
% Ted Amdur
% 10/25/22

% %Necessary cluster operations:
folder = fileparts(which(mfilename));
addpath(genpath(folder));
parpool('local',str2num(getenv('SLURM_CPUS_PER_TASK')))
usePMOD=1;

rng(1)
if ~usePMOD
    obsmatrix='obs_23_03_27.mat';
    load(obsmatrix); %From makeobsmatrix.m
    load excludeMask_23_03_27.mat %from exclude_fliers_22_04_26.m
else
    obsmatrix='obs_23_05_10_pmod.mat';
    load(obsmatrix); %From makeobsmatrix.m
    load excludeMask_PMOD_23_06_09.mat %from exclude_fliers_22_04_26.m
end
valM(excludeMask) = NaN;
oM(excludeMask) = false;
fracEx=0.5; %Exclude this fraction of results
reps=1000; %Length of kept chain
tN=2; %Number of experiments per observer

%Set options to be the same as the mainExperiment in btsi_23_03_04.m
opts.burn = 1000; %Number of burn-in reps assumed for chain length analysis
opts.reps=opts.burn+reps; %Total length of chain, including burn-in
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
if ~usePMOD
    opts.obsmatrix='obs_23_03_27.mat, half of an observer removed';
else
    opts.obsmatrix='obs_23_03_27.mat, half of an observer removed';
end

    
parfor col=1:length(colLabels)
    tic;
    len=ceil(fracEx.*sum(oM(:,col)))+1; %The number of observations to be excised
    iXM=NaN(tN,len); %Record of the observations excised from each experiment (rows)
    AM=zeros(tN,3,reps); %observation model draws for each experiment
    sigYM=zeros(tN,reps); %noise draws for each experiment
    pM=NaN(tN,len,reps); %Predictions for missing obs
    xM=NaN(tN,len,reps); %TSI predictions for missing obs times
    for ii=1:tN
        [valE,oE,iX]=exciseobs(valM,oM,col,fracEx); %Cut out a random half from an observer
        len=length(iX);
        %[obsPrior] = estimateobspriors(valE,oE,colLabels,1);
        
        %Run BTSI for valM with a set of observations excised
        [xAll,sigY,~,~,~,A,tau,outDat] = runchain_22_04_25(valE,oE,dateM,colLabels,opts);
        scaling=outDat.scaling(col);
        offset=outDat.offset(col);
        
        %[A,xAll,sigY,t] = runevalchain_22_10_25(valE,oE,colLabels,obsPrior);
         %pred=[ones(len,1) xAll(iX,ii) tau(iX,col)]; %make predictor matrix
         
         
         pO=zeros(len,reps); %Hold all predictions from this experiment
         for iN=1:reps
             pred=[ones(len,1) xAll(iX,iN) tau(iX,col)]; %make predictor matrix
             p=(squeeze(A(col,:,iN))*pred'.*scaling)'; %Predict from observation model
             p=p+randn(len,1).*sqrt(sigY(col,iN)).*scaling; %Add in assessed uncertainty
             pO(:,iN)=p;
         end
         iXM(ii,:)=iX';
         AM(ii,:,:)=squeeze(A(col,:,:));
         sigYM(ii,:)=sigY(col,:);
         pM(ii,:,:)=pO;
         xM(ii,:,:)=xAll(iX,:);
    end
    oTest(col).iX=iXM;
    oTest(col).A=AM;
    oTest(col).sigY=sigYM;
    oTest(col).p=pM;
    oTest(col).x=xM;
    oTest(col).tRun=toc;
end
save('modelevalPMOD_23_06_09.mat','oTest')