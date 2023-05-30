%Supplementary plots for  draft of TSI reconstruction manuscript. Built to take the 
%ar2 formulation
%
% Ted Amdur
% 9/20/22, updated from plotssupp_22_07_21.m

clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
btsiCompare=0; %Set to 1 to Plot BTSI from default vs different BTSI formulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS
btsiCompareTable=1; %Calculate solar constant and amplitude trends (w/ 95%CI) for alts
radForcingEffect=0; %1 to calculate global mean surface temperature effect from 
                    %the proposed degree of additional solar radiative
                    %forcing from BTSI
table1=0; %1 to calculate cycle trends for other TSI reconstructions
satireCompare=0;%1 to examine relationship between SATIRE-S and BTSI, 0 otherwise
satireTrendAnalysis=0; %1 to study modulation to SATIRE facular scaling
table1Supplement=0;%1 to calculate table 1 over whatever duration the product uses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fSize = 16;

excludeFliers=1;

load ar2_22_07_12.mat
dateStruct = getdates;
startDate = datejd(dateStruct.all(1));
endDate = datejd(dateStruct.all(2));
dates = [juliandate(datetime(1980,2,1)) juliandate(datetime(2021,12,31))]; %Use full interval
%dates=[juliandate(datetime(1980,2,1)) juliandate(datetime(2015,10,1))]; %Use shared interval
%Load data, with colLabels corresponding to observer source for each column
load obs_22_7_12.mat %From makeobsmatrix.m
valMAll = valM;dateMAll=dateM;
if excludeFliers %excludeFliers
    load excludeMask22_6_29.mat
    valM(excludeMask) = NaN;
    oM(excludeMask) = false;
    dateM=dateMAll;
end

%Reorder to match labels
%Order for proxies is sunspots first then MgII
%Order for sats is chronological from first observation: HF, ACRIM1, ERBE, ACRIM2,
%VIRGO/SOHO, ACRIM3, TIM/SORCE, PREMOS/PICARD, TCTE, TSIS-1
%Order as proxies followed by satellites in chronological order
lI=[8;4;6;1;5;2;12;3;9;7;10;11];
%sO=[7;3;6;4;12;5;9;8;10;11];
%Reorder arrays used in this plotting script
A=A(lI,:,:);
colLabels=colLabels(lI);
offsets=offsets(lI);
oM=oM(:,lI);
sigY=sigY(lI,:);
t=t(:,lI);
valM=valM(:,lI); valMAll=valMAll(:,lI);


%Colormap for plots
c1 = [51; 136; 68; 17; 153; 221; 204; 136; 170];
c2 = [34; 204; 170; 119; 153; 204; 102; 34; 68];
c3 = [136; 238; 153; 51; 51; 119; 119; 85; 153];
c = [c1 c2 c3]; c = c./255;
c(10,:) = c(1,:);

if btsiCompare
    showTrend = 1;
    smoothWindow = 6; %set smoothing (months)
    figure('Position',[10 10 1000 1200])
    
    cColor = get(gca,'colororder');
    %Get CI of our estimate,plot
    tsix = prctile(xAll',[.5 5 50 95 99.5])';
    for iS = 1:size(tsix,2)
        tsix(:,iS) = smoothPH(tsix(:,iS),smoothWindow);
    end
    [~,tsiO] = meaninterval(dateM,tsix(:,3),1990,2010);
    tsix = tsix-tsiO;
    x2 = [dateM', fliplr(dateM')];
    fill(x2,[tsix(:,1)',fliplr(tsix(:,end)')], [.85 .85 .85],'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[tsix(:,2)',fliplr(tsix(:,4)')], [.75 .75 .75],'FaceAlpha',0.5,...
        'LineStyle','none');
    
    %Get AR(1) reconstruction
    load ar1_22_07_12.mat
    x1 = prctile(xAll',50)';
    x1 = smoothPH(x1,smoothWindow);
    [~,x10] = meaninterval(dateMAll,x1,1990,2010);
    x1 = x1-x10; 
    
    %Get AR(3) reconstruction
    load ar3_22_07_12.mat
    x3 = prctile(xAll',50)';
    x3 = smoothPH(x3,smoothWindow);
    [~,x30] = meaninterval(dateMAll,x3,1990,2010);
    x3 = x3-x30; 
    
    %Get reconstruction keeping satellite artifacts
    load ar2_22_07_12_all.mat
    xa = prctile(xAll',50)';
    xa = smoothPH(xa,smoothWindow);
    [~,xa0] = meaninterval(dateMAll,xa,1990,2010);
    xa = xa-xa0; 
    
    %Get reconstruction scaling proxies with 3 most recent cycles
    load ar2_22_07_12_3cycle.mat
    x3c = mean(xAll,2);
    x3c = smoothPH(x3c,smoothWindow);
    [~,x3c0] = meaninterval(dateMAll,x3c,1990,2010);
    x3c = x3c-x3c0; 
    
    %Get reconstruction using PMOD corrections
    load ar2_22_07_12_pmodcorrections.mat
    xp = mean(xAll,2);
    xp = smoothPH(xp,smoothWindow);
    [~,xp0] = meaninterval(dateMAll,xp,1990,2010);
    xp = xp-xp0; 
    
    %Get reconstruction using VIRGO baseline
    load ar2_22_07_12_VIRGO.mat
    xv = mean(xAll,2);
    xv = smoothPH(xv,smoothWindow);
    [~,xv0] = meaninterval(dateMAll,xv,1990,2010);
    xv = xv-xv0; 
    
    ind = 1; %Index our contribution
    h(ind) = plot(dateM,tsix(:,3),'LineWidth',3,'Color','k');
    legendtxt(ind) = "BTSI";
    ind = ind + 1;
    hold on
    h(ind)=plot(dateMAll,x1,'--','LineWidth',1.5);
    legendtxt(ind)="BTSI AR(1)";
    ind=ind+1;
    hold on
    h(ind)=plot(dateMAll,x3,'--','LineWidth',1.5);
    legendtxt(ind)="BTSI AR(3)";
    ind=ind+1;
    h(ind)=plot(dateMAll,xa,'--','LineWidth',1.5);
    legendtxt(ind)="BTSI no artifact removal";
    ind=ind+1;
    h(ind)=plot(dateMAll,x3c,'--','LineWidth',1.5);
    legendtxt(ind)="BTSI 3 cycle regression";
    ind=ind+1;
    h(ind)=plot(dateMAll,xp,'--','LineWidth',1.5);
    legendtxt(ind)="BTSI PMOD satellite corrections";
    ind=ind+1;
    h(ind)=plot(dateMAll,xv,'--','LineWidth',1.5);
    legendtxt(ind)="BTSI VIRGO baseline";

    legend(h,legendtxt)
    legend boxoff
    set(gca,'FontSize',fSize)
    xlabel('Year')
    ylabel('TSI Anomaly (W/m^{2})')
    xlim([datetime(1978,1,1) datetime(2022,1,1)])
    ylim([-0.8 1.6])
    saveas(gcf,'plots/tsialternatives_22_07_12.png')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS
if btsiCompareTable
    %Produce table of TSI behavior for different alternatives, including
    %50th, 2.5th, and 97.5th percentile estimates for amplitude trend and
    %solar constant trend
    model=[];min=[];amp=[];total=[];
    warning('off','MATLAB:rankDeficientMatrix') %Some realizations for quant_reg are rank deficient
    intI=dateM>=datejd(dates(1))&dateM<=datejd(dates(2));
    
    
    
    %First entry is default
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"BTSI Default");
    xM=xMLR; %Save for proxy-model version at end
    %Get AR(1) reconstruction
    load ar1_22_07_12.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"AR(1)");
    
    %Get AR(3) reconstruction
    load ar3_22_07_12.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"AR(3)");
    
    %Get reconstruction keeping satellite artifacts
    load ar2_22_07_12_all.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"No artifact removal");
    
    %Get reconstruction scaling proxies with 3 most recent cycles
    load ar2_22_07_12_3cycle.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"3 cycle regression");
    
    %Get reconstruction using PMOD corrections
    load ar2_22_07_12_pmodcorrections.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"PMOD satellite corrections");
    
    
    %Get reconstruction using VIRGO baseline
    load ar2_22_07_12_VIRGO.mat
    
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"VIRGO baseline");
    
    %Get reconstruction using satellite-only
    load ar2_22_07_12_satonly.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"Satellite only");
    
    %Get reconstruction using proxy-model
    eI=sum(isnan(xM),2)==0;
    xMM=xM(eI,:);dateMM=dateM(eI);
    intII=intI(eI);
    [min95,amp95,total95]=getaltstats(dateMM(intII),xMM(intII,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"Proxy model");
    
    %Make table from results
    compare=table(model,amp,min,total);
end
if radForcingEffect
    %First, load the global mean surface temperature response from CMIP6
    %ensemble and HADCRUT5
    t=ncread('fig_3_4_panel_a.nc','time');
    tas=ncread('fig_3_4_panel_a.nc','tas');
    hadcrut=tas(:,61);
    cmip6ensemble=tas(:,60);
    
    %Then, load annual TSI time series from HEPPA-SOLARIS and BTSI
    [BTSI,yrBTSI] = monthlytoyearly(dateM,xAll);
    
    load oTSI_22_09_20.mat
    [HS,yrHS]=monthlytoyearly(oTSI(1).mthdatetime,oTSI(1).mthtsifull);
    
    %Subset the period of record for all obs
    hadcrut=hadcrut(t>1977 & t<2021); cmip6ensemble=cmip6ensemble(t>1977 & t<2021);
    t=t(t>1977 & t<2021);
    BTSI=BTSI(yrBTSI>1977 & yrBTSI<2021,:); yrBTSI=yrBTSI(yrBTSI>1977 & yrBTSI<2021);
    HS=HS(yrHS>1977 & yrHS<2021);yrHS=yrHS(yrHS>1977 & yrHS<2021);
    
    HS=HS-mean(HS); BTSI=BTSI-mean(BTSI(:));
    cmip2=(mean(BTSI,2)-HS).*0.1+cmip6ensemble;
end
if table1
    load oTSI_22_09_20.mat
    
    %Get anomaly relative to 1990-2010 mean
    datePMOD=oTSI(7).mthdatetime;
    tsiPMOD=oTSI(7).mthtsi;
    datePMOD=dateshift(datePMOD,'start','month');
    tsiPMOD=meaninterval(datePMOD,tsiPMOD,1990,2010);
    
    dateSolid=oTSI(9).mthdatetime;
    tsiSolid=oTSI(9).mthtsi;
    dateSolid=dateshift(dateSolid,'start','month');
    tsiSolid=meaninterval(dateSolid,tsiSolid,1990,2010);
    
    dateROB=oTSI(8).mthdatetime;
    tsiROB=oTSI(8).mthtsi;
    dateROB=dateshift(dateROB,'start','month');
    tsiROB=meaninterval(dateROB,tsiROB,1990,2010);
    
    dateNRL=oTSI(4).mthdatetime;
    tsiNRL=oTSI(4).mthtsi;
    dateNRL=dateshift(dateNRL,'start','month');
    tsiNRL=double(meaninterval(dateNRL,tsiNRL,1990,2010));
    
    dateSAT=oTSI(5).mthdatetime;
    tsiSAT=oTSI(5).mthtsi;
    dateSAT=dateshift(dateSAT,'start','month');
    tsiSAT=meaninterval(dateSAT,tsiSAT,1990,2010);
    
    dateAR6=oTSI(1).mthdatetime;
    tsiAR6=double(oTSI(1).mthtsi);
    dateAR6=dateshift(dateAR6,'start','month');
    tsiAR6=meaninterval(dateAR6,tsiAR6,1990,2010);
    
    %Get stats on each cycle minimum
    xms=mean(xAll,2); %TSI series used for minima calculation
    for ii = 1:size(dateStruct.cycles,1)-1
        cycleI = dateM > datejd(dateStruct.cycles(ii,1)+2*365) & dateM < ...
            datejd(dateStruct.cycles(ii,2)+2*365);
        minI = xms < quantile(xms(cycleI),0.10); %Lowest Nth percentile of obs is selected
        minCycleI{ii} = and(minI,cycleI);
        mcPMOD{ii}=ismember(datePMOD,dateM(minCycleI{ii}));
        mcSolid{ii}=ismember(dateSolid,dateM(minCycleI{ii}));
        mcNRL{ii}=ismember(dateNRL,dateM(minCycleI{ii}));
        mcSAT{ii}=ismember(dateSAT,dateM(minCycleI{ii}));
        mcROB{ii}=ismember(dateROB,dateM(minCycleI{ii}));
    end
    minCycleI=minCycleI';
    meanPMOD=setstats(tsiPMOD,mcPMOD');
    
    meanSolid=setstats(tsiSolid,mcSolid');
    
    meanNRL=setstats(tsiNRL,mcNRL');
    
    meanSAT=setstats(tsiSAT,mcSAT');
    
    meanROB=setstats(tsiROB,mcROB');
    %BTSI
    for ii=1:size(xAll,2)
        meanBTSI(ii,:)=setstats(meaninterval(dateM,xAll(:,ii),1990,2010),minCycleI);
    end
    meanMeanBTSI=mean(meanBTSI,1);
    
    %Get trend stats on each reconstruction over time they're all extent
    [a,m,t]=reconstructiontrends(datePMOD,tsiPMOD,dates);
    PMOD=[meanPMOD(end)-meanPMOD(end-1);a;m;t];
    
    [a,m,t]=reconstructiontrends(dateSolid,tsiSolid,dates);
    Solid=[meanSolid(end)-meanSolid(end-1);a;m;t];
    
    [a,m,t]=reconstructiontrends(dateROB,tsiROB,dates);
    ROB=[meanROB(end)-meanROB(end-1);a;m;t];
    
    [a,m,t]=reconstructiontrends(dateNRL,tsiNRL,dates);
    NRL=[meanNRL(end)-meanROB(end-1);a;m;t];
    
    [a,m,t]=reconstructiontrends(dateSAT,tsiSAT,dates);
    SAT=[meanSAT(end)-meanSAT(end-1);a;m;t];
    
    [a,m,t]=reconstructiontrends(dateAR6,tsiAR6,dates);
    AR6=[NaN;a;m;t];
    
    %BTSI
    for ii=1:1000
        interval=floor(size(xAll,2)./1000);
        [a,m,t]=reconstructiontrends(dateM,xAll(:,ii.*interval),dates);
        aa(ii)=a;mm(ii)=m;tt(ii)=t;
    end
    
    BTSI=[meanMeanBTSI(end)-meanMeanBTSI(end-1);mean(aa);mean(mm);mean(tt)];
    
    reconstructions=table(PMOD,Solid,ROB,NRL,SAT,AR6,BTSI);
end
if satireCompare
    plotFlag1=0;
    plotFlag2=0;
    plotFlag3=1;
    dates1 = [datetime(1978,11,1) datetime(1999,1,31)]; %Period without magnetograms
    dates2 = [datetime(1999,2,1) datetime(2021,11,31)]; %Period with magnetograms
    load oTSI_22_10_27.mat
    tsiS = oTSI(5).mthtsi;
    dateS=oTSI(5).mthdatetime;
    dateS=dateshift(dateS,'start','month');
    tsiS1=tsiS(dateS>=dates1(1) & dateS<dates1(2));
    dateS1=dateS(dateS>=dates1(1) & dateS<dates1(2));
    tsiS2=tsiS(dateS>=dates2(1) & dateS<dates2(2));
    dateS2=dateS(dateS>=dates2(1) & dateS<dates2(2));
    
%     
%     [~,tsiSO] = meaninterval(dateS,tsiS,2000,2010);
    
    btsi=mean(xAll,2);
    btsi1=btsi(dateM>=dates1(1) & dateM<dates1(2));
    dateM1=dateM(dateM>=dates1(1) & dateM<dates1(2));
    btsi2=btsi(dateM>=dates2(1) & dateM<dates2(2));
    dateM2=dateM(dateM>=dates2(1) & dateM<dates2(2));
%     [~,tsiBO] = meaninterval(oTSI(5).mthdatetime,tsiA,1990,2010);

    if plotFlag1
        figure2
        plot(dateM1,btsi1-mean(btsi1))
        hold on
        plot(dateS1,tsiS1-mean(tsiS1))
        legend('BTSI','SATIRE-S')
        set(gca,'FontSize',16)
    end
    if plotFlag2
        figure2
        plot(dateM2,btsi2-mean(btsi2))
        hold on
        plot(dateS2,tsiS2-mean(tsiS2))
        legend('BTSI','SATIRE-S')
        set(gca,'FontSize',16)
    end
    if plotFlag3
        figure2
        h(1)=plot((btsi1-mean(btsi1)),(tsiS1-mean(tsiS1)),'.');
        hold on
        h(2)=plot((btsi2-mean(btsi2)),(tsiS2-mean(tsiS2)),'.');
        hold on
        line([-.6 1.41],[-.6,1.41],'Color','k')
        xlim([-.75 1.45])
        ylim([-.75 1.45])
        xlabel('BTSI Anomaly (W/m^{2})')
        ylabel('SATIRE-S Anomaly (W/m^{2})')
        legendtxt(1)=string(['1978-1999: r=' num2str(corr(btsi1,tsiS1),2)]);
        legendtxt(2)=string(['1999-2021: r=' num2str(corr(btsi2,tsiS2),2)]);
        legend(h,legendtxt,'Location','NorthWest')
        legend boxoff
        set(gca,'FontSize',16)
    end
    
%     hold on %Plot SATIRE-S again during period of high-quality MDI/HMI data
%     dateMDIHMI=oTSI(5).mthdatetime(oTSI(5).mthdatetime>datetime(1999,02,02));
%     h(ind) = plot(dateMDIHMI,tsiA(oTSI(5).mthdatetime>datetime(1999,02,02))-tsiSO,'LineWidth',2.5,...
%         'Color',c(4,:));
end
if satireTrendAnalysis
    %download satire contribution data and full SATIRE product
    load satirecontributions.mat
    load satiretsi.mat
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pre-Processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Aug. 3, 2016 is missing in the contributions dataset - shift
    %everything up one day to account for this
    sC=satirecontributions;s=satiretsi;
    postI=sC.JD>juliandate(datetime(2016,8,3));
    sC.JD(postI)=sC.JD(postI)-1;
    
    %Remove days when facular contribution is negative or sunspot
    %contribution is positive
    iExclude=sC.fac < 0 | sC.spots > 0;
    sC=sC(~iExclude,:);
    
    %Remove obvious fliers/artifacts using normalized anomaly from running
    %mean
    fac20=smoothPH(sC.fac,20);spot20=smoothPH(sC.spots,20);
    nFac=normPH(sC.fac-fac20);nSpot=normPH(sC.spots-spot20);
    exclude5Sig=abs(nFac)>5 | abs(nSpot)>5;
    sC=sC(~exclude5Sig,:);
    %64 days are thus excluded (out of 10,562 original)
    
    %Assess which period has sufficient monthly coverage
    [sCTSI,mthDate,nObs] = dailytomonthly(sC.JD,sC.TSI);
    [sCFac,~,~] = dailytomonthly(sC.JD,sC.fac);
    [sCSpot,~,~] = dailytomonthly(sC.JD,sC.spots);
    %Based upon this, the period 3/1999 through 11/2021 has 93% daily coverage
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Regression analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Over the interval 3/1999 to 11/2021, regress BTSI over facular and
    %spot contribution
    iB=dateM>=datetime(1999,3,1)&dateM<datetime(2021,12,1);
    tsi=mean(xAll(iB,:),2);date=dateM(iB);
    iS=mthDate>=datetime(1999,3,1)&mthDate<datetime(2021,12,1);
    sCFacSub=sCFac(iS);sCSpotSub=sCSpot(iS);
    dateS=mthDate(iS);
    
    b=regress(tsi,[ones(length(tsi),1) sCFacSub sCSpotSub]);
    
    
    %Determine effect of this modification upon SATIRE record
    iV=nObs>10;
    [a,m,t]=reconstructiontrends(mthDate(iV),sCFac(iV)+sCSpot(iV),dates);
    [a2,m2,t2]=reconstructiontrends(mthDate(iV),sCFac(iV).*b(2)+sCSpot(iV),dates);
    
    %plot result of analysis
    plotI=0;
    if plotI
        figure
        plot(mthDate(iV),sCFac(iV)+sCSpot(iV),'.')
        hold on
        plot(mthDate(iV),sCFac(iV).*b(2)+sCSpot(iV),'.')
        hold on
        xlabel('Year')
        ylabel('TSI relative to quiet Sun (W/m^{2})')
        legend('Original','Scaled facular contribution')
        set(gca,'FontSize',16)
    end
end
if table1Supplement
    %Table 1, for the full duration of each reconstruction
    load oTSI_22_09_20.mat
    
    %Get anomaly relative to 1990-2010 mean
    datePMOD=oTSI(7).mthdatetime;
    tsiPMOD=oTSI(7).mthtsi;
    tsiPMOD=meaninterval(datePMOD,tsiPMOD,1990,2010);
    
    dateSolid=oTSI(2).mthdatetime;
    tsiSolid=oTSI(2).mthtsi;
    tsiSolid=meaninterval(dateSolid,tsiSolid,1990,2010);
    
    dateNRL=oTSI(4).mthdatetime;
    tsiNRL=oTSI(4).mthtsi;
    tsiNRL=double(meaninterval(dateNRL,tsiNRL,1990,2010));
    
    dateSAT=oTSI(5).mthdatetime;
    tsiSAT=oTSI(5).mthtsi;
    tsiSAT=meaninterval(dateSAT,tsiSAT,1990,2010);
    
    dateAR6=oTSI(1).mthdatetime;
    tsiAR6=double(oTSI(1).mthtsi);
    tsiAR6=meaninterval(dateAR6,tsiAR6,1990,2010);
    
    %Get trend stats on each reconstruction over time they're all extent
    [a,m,t]=reconstructiontrends(datePMOD,tsiPMOD,dates);
    PMOD=[a;m;t];
    
    [a,m,t]=reconstructiontrends(dateSolid,tsiSolid,dates);
    Solid=[a;m;t];
    
    [a,m,t]=reconstructiontrends(dateNRL,tsiNRL,dates);
    NRL=[a;m;t];
    
    [a,m,t]=reconstructiontrends(dateSAT,tsiSAT,dates);
    SAT=[a;m;t];
    
    [a,m,t]=reconstructiontrends(dateAR6,tsiAR6,dates);
    AR6=[a;m;t];
    
    %BTSI
    for ii=1:1000
        interval=floor(size(xAll,2)./1000);
        [a,m,t]=reconstructiontrends(dateM,xAll(:,ii.*interval),dates);
        aa(ii)=a;mm(ii)=m;tt(ii)=t;
    end
    
    BTSI=[mean(aa);mean(mm);mean(tt)];
    
    reconstructions=table(PMOD,Solid,NRL,SAT,AR6,BTSI);
end

function [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,modelName)
model=[model;modelName];
min=[min;min95];
amp=[amp;amp95];
total=[total;total95];
end
function [cycleAvg]=cyclestats(x,dates)
    %Return average of input x over the julian day intervals provided
    for ii=1:size(dates,1)
        cycleAvg(ii)=nanmean(x(cycI));
    end
end
function [setAvg]=setstats(x,dateI)
    %Return average of input x over the indexed values provided
    for ii=1:size(dateI,1)
        try
            setAvg(ii)=mean(x(dateI{ii}));
        catch
            setAvg(ii)=NaN;
        end
    end
end
function [min95,amp95,total95]=getaltstats(dateM,xAll)
%RETURN TREND AND 95% CI FOR REALIZATIONS OF TSI
% OUTPUT:
%           min95: trend of solar constant in form [2.5,50,97.5] percentile
%           amp95: trend of amplitude in form [2.5,50,97.5] percentile
for ii=1:size(xAll,2)
    b5(:,ii)=quant_reg([ones(size(dateM,1),1) juliandate(dateM)],xAll(:,ii),0.05);
    b95(:,ii)=quant_reg([ones(size(dateM,1),1) juliandate(dateM)],xAll(:,ii),0.95);
    b50(:,ii)=quant_reg([ones(size(dateM,1),1) juliandate(dateM)],xAll(:,ii),0.5);
end
min95=prctile(b5(2,:),[2.5 50 97.5]).*(365.25*10);
amp95=prctile(b95(2,:)-b5(2,:),[2.5 50 97.5]).*(365.25*10);
total95=prctile(b50(2,:),[2.5 50 97.5]).*(365.25*10);
end
function [ampTrend,minTrend,totTrend]=reconstructiontrends(dateM,x,dates)
    %Return trends given date vector and reconstruction vector, within the
    %dates interval
    dI=dateM>=datejd(dates(1))&dateM<=datejd(dates(2));
    dateM=dateM(dI);x=x(dI);
    dateM=dateM(~isnan(x));x=x(~isnan(x));
    b5=quant_reg([ones(size(dateM,1),1) juliandate(dateM)],x,0.05);
    b95=quant_reg([ones(size(dateM,1),1) juliandate(dateM)],x,0.95);
    b50=quant_reg([ones(size(dateM,1),1) juliandate(dateM)],x,0.50);
    ampTrend=(b95(2)-b5(2)).*(365.25*10);
    minTrend=b5(2).*(365.25*10);
    totTrend=b50(2).*(365.25*10);
end


