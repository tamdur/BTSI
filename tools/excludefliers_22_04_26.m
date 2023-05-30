%For a given run of runchain, identify observations that are more than X
%standard deviations from the expected values, as well as intervals at 
%beginning and end of satellite record that are suspicious, create a NaN 
%mask that allows for the exclusion of these values
%
% Ted Amdur
% 4/26/22, modified from excludefliers_22_2_02.m to use new runchain script

clearvars

cutoff = 0.15; %Fraction of full satellite record at beginning and end to exclude
pVal = 0.01; %Pvalue for which ends are removed
excludeSig = 3; %number of standard deviations to set for exclusion criteria
excludeVals.cutoff = cutoff; excludeVals.pVal = pVal; excludeVals.excludeSig=excludeSig;

%Load data, with colLabels corresponding to observer source for each column
load ar2_23_03_27_all.mat %Load a runchain output with no excision
load(outDat.obsmatrix)  %From makeobsmatrix.m
reps=size(A,3);
sindex=outDat.satindex;
nObs=sum(sindex);
T=size(xAll,1);
valSat=outDat.opts.valM(:,sindex);
%Get the output of predictions for each observer
samples = 1000; %Number of drawn samples to make CI
[ym,y5,y95,yAll] = estimatekalmanciy(A(sindex,:,:),xAll,sigY(sindex,:),tau(:,sindex)); 
for ii = 1:nObs
    rm(:,ii) = valSat(:,ii)- ym(:,ii);
end

%Get the standard error for each observer
SE = sqrt(mean(sigY(sindex,:),2));

excludeMask = false(size(valM));
sI=find(sindex);
for ii = 1:length(sI) %Skip the proxies
    %First go through the beginning of each record
    satI = find(oM(:,sI(ii)));
    lowI = satI(satI<quantile(satI,cutoff)); highI = satI(satI>quantile(satI,1-cutoff));
    pLow = []; pHigh = [];
    %Remove early interval if below cutoff
    for iL = 1:length(lowI)
        [~,p] = ztest(rm(lowI(1:iL),ii),0,SE(ii));
        pLow(iL) = p;
    end
    if min(pLow) < pVal
        [~,imin] = min(pLow);
        excludeMask(lowI(1:imin),sI(ii)) = true;
    end
    
    %Remove late interval if below cutoff
    for iL = 1:length(highI)
        [~,p] = ztest(rm(highI(end-iL+1:end),ii),0,SE(ii));
        pHigh(iL) = p;
    end
    if min(pHigh) < pVal
        [~,imin] = min(pHigh);
        excludeMask(highI(end-imin+1:end),sI(ii)) = true;
    end   
end

%Then, for each observer, find all observations that are more than X standard
%deviations from expected
for ii = 1:nObs
    excludeMask(abs(rm(:,ii)) > (excludeSig.*SE(ii)),ii) = true;
end



save('excludeMask_23_03_27.mat','excludeMask','excludeVals')

function [ym,y25,y975,yAll] = estimatekalmanciy(A,xAll,sigY,t)
%Return a confidence interval for observation variable given hidden process
%
%   INPUTS:
%           sIn: sIn object returned by gibbstests.m 
%           A: A parameter object returned by gibbstests.m
%           xAll: Estimates of hidden process
%           sigY: A parameter object for sigma of observers returned by gibbstests.m
%           rI: the realization from a multi-realization simulation to be
%               used


%values to be shared by all iterations
reps=size(A,3);
nObs=size(A,1);
T=size(xAll,1);
yAll = NaN(nObs,T,reps);
for ii=1:reps
    Xi=[ones(T,1) xAll(:,ii) ones(T,1)];
    for iN = 1:nObs
        Xi(:,3)=t(:,iN);
        yAll(iN,:,ii) = Xi*squeeze(A(iN,:,ii))'+ randn(T,1).*sqrt(sigY(iN,ii));
    end
end

ym = mean(yAll,3)';
y25 = quantile(yAll,0.025,3)';
y975 = quantile(yAll,0.975,3)';

end

