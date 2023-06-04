function [xAll,sigY,sigX,theta,a,A,tau,outDat] = runchain_22_04_25(valM,oM,dateM,colLabels,opts)
%runchain_22_04_25.m PRODUCE GIBBS SAMPLING OUTPUT OF BAYESIAN KALMAN
%FILTER MODEL USING SATLELLITE AND PROXY DATA.
%
% Ted Amdur 04/25/22
%Based upon runchain_22_03_23b.m that shows innovations with time, also
%allows for different ar order models. The following is an
%extensively-modified version of the code from Ex. 4 of Chapter 4 of
%Applied Bayesian Econometrics for Central Bankers, by Andrew Blake and
%Haroon Mumtaz.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Headers to modify
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outDat.script=mfilename; %Save name of script
if ~exist('valM','var') || isempty(valM) %Load default observation array, otherwise load provided one
    obsmatrix='obs_23_03_27'; %Load data array, with colLabels corresponding to observer source for each column
    load(obsmatrix); %From makeobsmatrix.m
else
     obsmatrix=opts.obsmatrix;
end
    
%Create default settings if not specified
opts = checkopts(opts);
%Develop a set of monthly observations from satellite and proxy
%observations
dateS=getdates;



%Twelve columns of valM correspond to the following observers:
%     "ACRIM1/SMM"
%     "ACRIM2/UARS"
%     "ACRIM3"
%     "BremenMgII"
%     "ERBE/ERBS"
%     "NIMBUS-7/HF"
%     "PREMOS/PICARD"
%     "SILSO"
%     "SORCE/TIM"
%     "TCTE"
%     "TSIS-1"
%     "VIRGO/SOHO"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of main script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.excludeFliers %Code to remove outliers using a past run of BTSI
    load excludeMask_23_03_27.mat %from exclude_fliers_22_04_26.m
    valM(excludeMask) = NaN;
    oM(excludeMask) = false;
end
opts.valM=valM;
opts.oM=oM;
opts.dateM=dateM;
if opts.satOnly
    %Satellite-only variation
    oindex=[1 1 1 1 1 1 1 1 0 1 1 1]; %oindex=1 for observers with varying offset, 0 for fixed
    tindex=[0 0 0 0 0 0 0 0 0 0 0 0]; %tindex= no satellite drift
    %tindex=[1 1 1 0 1 1 1 0 1 1 1 1]; %tindex= satellite drift
    pindex=[0 0 0 1 0 0 0 1 0 0 0 0]; %sindex=1 for observers with non-identity scaling to TSI, 0 otherwise
    oM(:,logical(pindex))=false; %Eliminate use of observations from proxies
elseif opts.satOnlyDrift
    %Satellite-only variation
    oindex=[1 1 1 1 1 1 1 1 0 1 1 1]; %oindex=1 for observers with varying offset, 0 for fixed
    tindex=[1 1 1 0 1 1 1 0 1 1 1 1]; %tindex= satellite drift
    pindex=[0 0 0 1 0 0 0 1 0 0 0 0]; %sindex=1 for observers with non-identity scaling to TSI, 0 otherwise
    oM(:,logical(pindex))=false; %Eliminate use of observations from proxies
elseif opts.proxyModel
    %NRLTSI2 version
    oindex=[1 1 1 1 1 1 1 1 0 1 1 1]; %oindex=1 for observers with varying offset, 0 for fixed
    tindex=[0 0 0 0 0 0 0 0 0 0 0 0]; %tindex=1 for observers with time dependent drift, 0 otherwise
    pindex=[0 0 0 1 0 0 0 1 0 0 0 0]; %sindex=1 for observers with non-identity scaling to TSI, 0 otherwise
    oM(:,[1 2 3 5 6 7 10 11 12])=false; %Eliminate use of observations from eliminated sources
    dateI=dateM.Year >= 2003 & dateM.Year <= 2014;
    oM(~dateI,9)=false;
elseif opts.virgo
    %Set SOHO/VIRGO as reference satellite
    oindex=[1 1 1 1 1 1 1 1 1 1 1 0]; %oindex=1 for observers with varying offset, 0 for fixed
    tindex=[1 1 1 0 1 1 1 0 1 1 1 1]; %tindex=1 for observers with time dependent drift, 0 otherwise
    pindex=[0 0 0 1 0 0 0 1 0 0 0 0]; %sindex=1 for observers with non-identity scaling to TSI, 0 otherwise
else
    %Specify priors for H coefficients
    oindex=[1 1 1 1 1 1 1 1 0 1 1 1]; %oindex=1 for observers with varying offset, 0 for fixed
    tindex=[1 1 1 0 1 1 1 0 1 1 1 1]; %tindex=1 for observers with time dependent drift, 0 otherwise
    pindex=[0 0 0 1 0 0 0 1 0 0 0 0]; %sindex=1 for observers with non-identity scaling to TSI, 0 otherwise
end
satindex=~pindex;
%Remove removed observers
keepI=sum(oM,1)~=0;
valM=valM(:,keepI);oM=oM(:,keepI);oindex=oindex(keepI);tindex=tindex(keepI);
pindex=pindex(keepI);satindex=satindex(keepI);colLabels=colLabels(keepI);

if isfield(opts,'normalize') && opts.normalize %Normalize proxy variables for matrix calculations
    scaling=ones(size(valM,2),1);
    offset=zeros(size(valM,2),1);
    proxies=find(pindex);
    for ii=1:length(proxies)
        scaling(proxies(ii))=nanstd(valM(oM(:,proxies(ii)),proxies(ii)));
        offset(proxies(ii))=nanmean(valM(oM(:,proxies(ii)),proxies(ii)));
        valM(:,proxies(ii))=(valM(:,proxies(ii))-offset(proxies(ii)))./scaling(proxies(ii));
    end
    outDat.offset=offset;
    outDat.scaling=scaling;
end


tsi.data = valM;
tsi.dateM=dateM;
tsi.oM=oM;
tsi.T=size(tsi.data,1);
tsi.L=opts.lags;  %number of lags in the VAR
tsi.N=1;  %Number of factors/variables in the transition equation
tsi.NN=size(tsi.data,2);% Number of observers

tsi.cx = ones(tsi.T,1); %ones vector in first row of z for offset
tau =repmat(linspace(0,tsi.T./120,tsi.T)',[1 tsi.NN]); %Make time rows for t
for ii=1:tsi.NN
    TM=mean(tau(oM(:,ii),ii));
    if isnan(TM)
        TM=0;
    end
    tau(:,ii)=tau(:,ii)-TM;
end
tsi.tau=tau;
tsi.tDependence=true; %Set to true to use drift predictor
tsi.magDependent=opts.magDependent; %True to use heteroscedastic TSI innovation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 1: establish starting values and priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load priors for observation model
[tsi.H0, tsi.Hsig, T0, th0,tsi.eta0,tsi.Xsig] = getpriors(tsi.NN,oindex,tindex,pindex,valM,oM,...
    satindex,colLabels,opts,[],dateM);

if opts.proxyModel
    th0(3)=0.2688; %Give SORCE/TIM the th0 estimate one gets from whole model
    tsi.H0(1:2,2)=[1.8588;1.7844];
    tsi.Hsig(1:2,2)=ones(2,1).*1E-5;
end
%get an intial guess for the process
x0=zeros(tsi.T,1);
guess1 = nanmean(tsi.data(:,~pindex),2); %mean of satellite observations
x0(~isnan(guess1))=guess1(~isnan(guess1));
if opts.randomizeChain
    rng('shuffle')
    x0=x0+randn(length(x0),1).*(nanstd(x0)./5);
end
tsi.x0=x0;
tsi.X0=[x0(1) zeros(1,tsi.L-1)];%state vector X[t-1|t-1] of x and nth order lagged x
tsi.ns=size(tsi.X0,2);
tsi.V00=eye(tsi.ns);  %v[t-1|t-1]
tsi.rmat=1E9.*ones(tsi.NN,1); %arbitrary starting value for the variance of process x
tsi.Sigma=eye(tsi.N);  %arbitrary starting value for the variance of transition model errors

%Prepare arrays to be saved
NS=opts.reps-opts.burn;
outDat.contributionChain=zeros(NS,tsi.T,tsi.NN);
outDat.m=zeros(NS,1);
outDat.b=zeros(NS,1);
xAll=zeros(tsi.T,NS);
A=zeros(tsi.NN,3,NS);
a=zeros(tsi.L,NS);
sigY=zeros(tsi.NN,NS);
sigX=zeros(tsi.T,NS);
theta=zeros(tsi.NN,NS);

mm=1;%index for saved chain post-burnin
tic
for m=1:opts.reps %Gibbs sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 2: sample loadings that compose H through Bayesian regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hload,errorN,infI]=coeffsample(tsi,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 3: sample variance of the observers from inverse gamma distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tsi.rmat,th] = epsiloninvgamma(T0,th0,infI,errorN,tsi.NN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 4: sample estimates of autoregressive parameters alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tsi.alpha,tsi.X,tsi.Y]=arsample(tsi,opts.cmpStYr);
%sample VAR covariance for time-dependent X uncertainty
[tsi.Sigma,mSigma,bSigma]=epsilonmagdependent(tsi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 5: Run Kalman filter to estimate x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tsi.hload=hload;
opts.nonLin=false;
[tsi.x0,contributionChain,tsi.x2,tsi.F,tsi.H]=carterkohn(tsi,opts,m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 7: If burnin completed, store state and observation model estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.dispProgress
    if mod(m,100) == 0 %display progress of chain
        disp([num2str(m) ' reps completed'])
    end
end
if m>opts.burn
    %save
    xAll(:,mm)=tsi.x2(:,1); %Reconstructed TSI
    for ii=1:tsi.NN
        A(ii,:,mm)=tsi.H(ii,[1 2 ii+tsi.L+1])'; %observation coefficients
    end
    a(:,mm)=tsi.F(1,:); %VAR autoregressive coefficients [lag1 lag2]
    sigY(:,mm)=tsi.rmat; %observer noise estimate
    sigX(:,mm) = tsi.Sigma; %TSI noise estimate
    theta(:,mm)=th; %Theta parameter estimate
    if exist('mSigma','var') && exist('bSigma','var')
        outDat.m(mm)=mSigma;
        outDat.b(mm)=bSigma;
    end
    if opts.logContributions
        outDat.contributionChain(mm,:,:)=contributionChain; %Innovation contributions
    end
    mm=mm+1;
end
    
end
if opts.dispProgress
   outDat.tElapsed=toc
end
opts.colLabels=colLabels;
outDat.H0=tsi.H0;outDat.Hsig=tsi.Hsig;outDat.T0=T0;outDat.th0=th0;
outDat.Xprior=tsi.eta0;outDat.Xsig=tsi.Xsig;
outDat.oindex=oindex;outDat.tindex=tindex;outDat.sindex=pindex;
outDat.satindex=satindex;
outDat.obsmatrix=obsmatrix;
outDat.opts=opts; %Save the input info
outDat.runDate=datetime;
if isfield(opts,'saveFile')
    save(opts.saveFile,'xAll','sigY','sigX','theta','a','A','tau','outDat','-v7.3')
end


end



