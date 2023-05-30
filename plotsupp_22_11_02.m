%Supplementary plots for  draft of TSI reconstruction manuscript.
%
% Ted Amdur
% 11/02/22, updated from plotssupp_22_09_20.m

clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each conditional, set to 1 to plot/calculate that figure/table, 0
% otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
satireCompare=0;%Examine relationship between SATIRE-S and BTSI, 0 otherwise
btsiCompare=0; %Plot BTSI from default vs different BTSI formulations
obsContributions=0; %Plot the relative contribution of each observer to BTSI over time
satsatcomp=1; %Plot overlapping satellite observation model predictions
priorsposteriors=0; %Plot three figures for: prior/posterior offset, prior/posterior drift, prior/posterior standard error
plotData=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE CALCULATIONS
table1=0; %Produce values for Table 1
table2=0; %Produce values for Table 2
btsiCompareTable=0; %Calculate solar constant and amplitude trends (w/ 95%CI) for alts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER CALCULATIONS
uncBTSI=0; %Calculate and plot the uncertainty in BTSI
solarvsanthroforcing = 0;%Calculate global mean surface temperature effect from 
                    %the proposed degree of additional solar radiative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fSize = 12;

excludeFliers=1;

load ar2_23_05_11_long.mat
dateStruct = getdates;
startDate = datejd(dateStruct.all(1));
endDate = datejd(dateStruct.all(2));
dates = [juliandate(datetime(1980,1,1)) juliandate(datetime(2021,12,31))]; %Use full interval
%dates = [juliandate(datetime(1978,11,1)) juliandate(datetime(2021,12,31))]; %T Response Interval
%dates=[juliandate(datetime(1980,2,1)) juliandate(datetime(2015,10,1))]; %Use shared interval
%Load data, with colLabels corresponding to observer source for each column
load(outDat.obsmatrix); %From makeobsmatrix.m
valMAll = valM;dateMAll=dateM;
if excludeFliers %excludeFliers
    load excludeMask_22_11_03.mat
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
tau=tau(:,lI);
valM=valM(:,lI); valMAll=valMAll(:,lI);
scaling=outDat.scaling(lI);


%Colormap for plots
c1 = [51; 136; 68; 17; 153; 221; 204; 136; 170];
c2 = [34; 204; 170; 119; 153; 204; 102; 34; 68];
c3 = [136; 238; 153; 51; 51; 119; 119; 85; 153];
c = [c1 c2 c3]; c = c./255;
c(10,:) = c(1,:);
cmat=        [ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    btsi=mean(xAll,2);
    btsi1=btsi(dateM>=dates1(1) & dateM<dates1(2));
    dateM1=dateM(dateM>=dates1(1) & dateM<dates1(2));
    btsi2=btsi(dateM>=dates2(1) & dateM<dates2(2));
    dateM2=dateM(dateM>=dates2(1) & dateM<dates2(2));
    
    %Plot figure
    figure2('Position',[10 10 1000 800])
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
    saveas(gcf,'plots/satirebtsicompare_22_11_07.png')
end
if btsiCompare
    showTrend = 1;
    smoothWindow = 6; %set smoothing (months)
    figure2('Position',[10 10 1000 1200])
    
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
    load ar1_22_11_04.mat
    x1 = prctile(xAll',50)';
    x1 = smoothPH(x1,smoothWindow);
    [~,x10] = meaninterval(dateMAll,x1,1990,2010);
    x1 = x1-x10; 
    
    %Get AR(3) reconstruction
    load ar3_22_11_04.mat
    x3 = prctile(xAll',50)';
    x3 = smoothPH(x3,smoothWindow);
    [~,x30] = meaninterval(dateMAll,x3,1990,2010);
    x3 = x3-x30; 
    
    %Get reconstruction keeping satellite artifacts
    load ar2_22_11_03_all.mat
    xa = prctile(xAll',50)';
    xa = smoothPH(xa,smoothWindow);
    [~,xa0] = meaninterval(dateMAll,xa,1990,2010);
    xa = xa-xa0; 
    
    %Get reconstruction scaling proxies with 3 most recent cycles
    load ar2_22_11_03_3cycle.mat
    x3c = mean(xAll,2);
    x3c = smoothPH(x3c,smoothWindow);
    [~,x3c0] = meaninterval(dateMAll,x3c,1990,2010);
    x3c = x3c-x3c0; 
    
    %Get reconstruction using PMOD corrections
    load ar2_22_11_04_pmodcorrections.mat
    xp = mean(xAll,2);
    xp = smoothPH(xp,smoothWindow);
    [~,xp0] = meaninterval(dateMAll,xp,1990,2010);
    xp = xp-xp0; 
    
    %Get reconstruction using VIRGO baseline
    load ar2_22_11_04_VIRGO.mat
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
    saveas(gcf,'plots/tsialternatives_22_11_07.png')
end
if obsContributions
     %Plot the relative contribution of each observer to the estimate of TSI
    smoothWindow = 36;
    conChain=outDat.contributionChain;
    conChain=conChain(:,:,lI);%Reorient to agree with formatting of this plotting function
    cn=squeeze(mean(abs(conChain),1));
    for ii=1:size(valM,2)
        cn(:,ii)=movmean(cn(:,ii),smoothWindow,'omitnan');
        cn(~oM(:,ii),ii)=NaN;
    end
    cn=cn./nansum(cn,2);
    
    %Plot
    figure2('Position',[10 10 1600 700])
    hold on
    ind=1;
    cp=colororder;
    for ii=1:2
        h(ind)=plot(dateM,cn(:,ind),'--','LineWidth',4,'Color',cp(ii,:));
        hold on
        ind=ind+1;
    end
    for ii=1:10
        h(ind)=plot(dateM,cn(:,ind),'LineWidth',3,'Color',c(ii,:));
        hold on
        ind=ind+1;
    end
    legend(h,colLabels,'Location','NorthWest','NumColumns',2)
    legend boxoff
    xlabel('Year')
    ylabel('Fractional contribution')
    set(gca,'FontSize',fSize)
    xlim([datetime(1979,1,1) datetime(2021,10,31)])
    %ylim([0 1])
    saveas(gcf,'plots/obscontribution_23_05_13.png')
    
end
if satsatcomp
    %Order is chronological from first observation: HF, ACRIM1, ERBE, ACRIM2,
    %VIRGO/SOHO, ACRIM3, TIM/SORCE, PREMOS/PICARD, TCTE, TSIS-1
    sO=[3;4;5;6;7;8;9;10;11;12];
    ind=1;
    for iS1 = 1:10
        for iS2 = 1:10
            if iS2 > iS1
                iOverlap = ~isnan(valM(:,sO(iS1))) & ~isnan(valM(:,sO(iS2)));
                iOverlapAll = ~isnan(valMAll(:,sO(iS1))) & ~isnan(valMAll(:,sO(iS2)));
                if sum(iOverlap) > 48
                    overlap(ind).sat1 = sO(iS1);
                    overlap(ind).sat2 = sO(iS2);
                    overlap(ind).ind = iOverlap;
                    overlap(ind).indall = iOverlapAll;
                    ind = ind + 1;
                end
            end
        end
    end
    figure2('Position',[10 10 2300 1300])
    tt = tiledlayout(7,6,'TileSpacing','Compact','Padding','Compact');
    ind=1;
    indD=1;
    cols=[1 3 4 5 6 7];
    rows=[2 3 4 5 6 7 9];
    deleteAxes=[];%Tiles that are empty and should be deleted
    for ii=1:42
        nexttile
        col=mod(ind-1,6)+1;
        row=ceil(ind/6);
        validPlot=0;
        for iO=1:length(overlap)
            if find(sO==overlap(iO).sat1)==cols(col) && find(sO==overlap(iO).sat2)==rows(row)
                s1 = overlap(iO).sat1; s2 = overlap(iO).sat2;
                oInd = overlap(iO).ind;
                oIndAll = overlap(iO).indall;
                xm=mean(xAll,2);
                residual2 = valM(oInd,s2)-xm(oInd);residual2=residual2-mean(residual2);
                residual1 = valM(oInd,s1)-xm(oInd);residual1=residual1-mean(residual1);
                plot(dateM(oInd),residual1,'.','MarkerSize',10,'Color',cmat(1,:));
                hold on
                plot(dateM(oInd),residual2,'.','MarkerSize',10,'Color',cmat(2,:));
                
                %Plot inferred trends for each satellite
                [t1,o1]=returntrend(A,tau(oInd,:),s1);
                [t2,o2]=returntrend(A,tau(oInd,:),s2);
                t1mm=mean(mean(t1,2)+mean(o1));
                t1m=mean(t1+o1,2)-t1mm;
                t1025=quantile(t1+o1,0.025,2)-t1mm;
                t1975=quantile(t1+o1,0.975,2)-t1mm;
                t2mm=mean(mean(t2,2)+mean(o2));
                t2m=mean(t2+o2,2)-t2mm;
                t2025=quantile(t2+o2,0.025,2)-t2mm;
                t2975=quantile(t2+o2,0.975,2)-t2mm;
                hold on
                h(1)=plot(dateM(oInd),t1m,'Color',cmat(1,:),'LineWidth',3.5);
                hold on
                plot(dateM(oInd),t1025,'Color',cmat(1,:),'LineWidth',2)
                hold on
                plot(dateM(oInd),t1975,'Color',cmat(1,:),'LineWidth',2)
                hold on
                h(2)=plot(dateM(oInd),t2m,'Color',cmat(2,:),'LineWidth',3.5);
                hold on
                plot(dateM(oInd),t2025,'Color',cmat(2,:),'LineWidth',2)
                hold on
                plot(dateM(oInd),t2975,'Color',cmat(2,:),'LineWidth',2)
                lgd=legend(h,colLabels(s1),colLabels(s2),'Location','SouthEast');
                lgd.FontSize=10;
                legend boxoff
                indD=indD+1;
                set(gca,'FontSize',fSize)
                if s1 < 6
                    ylim([-0.6 0.6])
                else 
                    ylim([-0.2 0.2])
                end
%                 rangeSet=[residual2;residual1];
%                 ylim([quantile(rangeSet,0.01),quantile(rangeSet,0.99)])
                validPlot=1;
            end      
        end
        if ~validPlot
            deleteAxes=[deleteAxes;indD];
        end
        ind=ind+1;
        indD=indD+1;
    end
    delete(tt.Children(53-deleteAxes))
    saveas(gcf,'plots/satresidual_23_05_13.png')
end
if priorsposteriors
    datesave='23_04_26'; %Date for figure name
    %Pull output from simulation
    reps = size(A,3);
    sC = squeeze(A(:,2,:));
    bC = squeeze(A(:,1,:));
    mC = squeeze(A(:,3,:));
    aP = a;
    sP = squeeze(sC(1:2,:))';
    oP = squeeze(bC(1:12,:))';
    lP = squeeze(mC(3:end,:))';
    rP = sigY;
    H0=outDat.H0(lI,:);
    Hsig=outDat.Hsig(lI,:);
    T0=outDat.T0(lI);
    th0=outDat.th0(lI);
    %Make the data that goes in the varpdf fields for each of 3 plots
    
    %First, do the offset variables
    varNames1 = ["o3";"o4";"o5";"o6";"o7";"o8";...
        "o10";"o11";"o12"];
    th1M1 = H0([3:8 10:12],1);
    th2M1 = sqrt(Hsig([3:8 10:12],1));
    vals1 = oP(:,3:end);
    %Then, do the drift variables
    varNames2 = ["l3";"l4";"l5";"l6";"l7";"l8";"l9";"l10";...
        "l11";"l12"];
    th1M2 = H0(3:end,3);
    th2M2 = sqrt(Hsig(3:end,3));
    vals2 = lP;
    %Last, do the noise variables
    varNames3=[];
    for ii = 3:12
        varNames3 = [varNames3;strcat("\epsilon_{",colLabels(ii),"}")];
    end
%     strcat("\sigma_{";"yN4";"yN5";"yN6";"yN7";"yN8";"yN9";...
%         "yN10";"yN11";"yN12"];
    th1M3 = T0(3:end);
    th2M3 = th0(3:end);
    vals3 = rP(3:end,:)'; 
    vals3(vals3 <= 0) = NaN;
    vals3 = vals3.^0.5;%Note 0.5 is choice to show std v var
    distType3 = repmat("invgamma",[10 1]);

    %------------------------------------------------------------------
    % First, plot estimated offsets 
    %------------------------------------------------------------------
    figure2('Position',[10 10 1600 600])
    %subplot('position',[.09 .7 .85 .27])
    offsetsI = [3:8 10:12]';
    numPlots = length(varNames1);
    pDim = ceil(sqrt(numPlots));
    for ii = 1:numPlots
        subplot(pDim,pDim,ii)
        hold on
        th1 = th1M1(ii) + offsets(offsetsI(ii));
        th2 = th2M1(ii);
        %Plot the prior distribution
        x = vals1;
        xL = min(x);xH = max(x);
        xPrior = randraw('norm', [th1 th2], [1e5 1] );
        xL = min([xPrior; vals1(:,ii)+ offsets(offsetsI(ii))]);
        xH = max([xPrior;vals1(:,ii)+ offsets(offsetsI(ii))]);
        pEdges = linspace(xL,xH,1000);
        [yP,xP] = histcounts(xPrior,pEdges); 
        [xP,yP,~] = histtoplot(yP,xP,50);
        %Plot the prior
        plot(xP,yP);
        %Plot the posterior distribution
        [yPost,xPost] = histcounts(vals1(:,ii)+ offsets(offsetsI(ii)),pEdges); 
        [xPost,yPost,~] = histtoplot(yPost,xPost,10);
        hold on
        plot(xPost,yPost)
        xlim([quantile(xPrior,0.001) quantile(xPrior,0.999)])
        ylim([0 2.*max(yP)])
        title(colLabels(offsetsI(ii)))
        xlabel("W/m^{2}")
        set(gca,'ytick',[])
        set(gca,'FontSize',fSize)
    end
    saveas(gcf,['plots/priorposterior1_' datesave '.png'])
    %------------------------------------------------------------------
    % Next, plot estimated linear drifts 
    %------------------------------------------------------------------
    figure2('Position',[10 10 1600 600])
    offsetsI = [3:12]';
    numPlots = length(varNames2);
    for ii = 1:numPlots
        subplot(3,4,ii)
        hold on
        th1 = th1M2(ii);
        th2 = th2M2(ii);
        %Plot the prior distribution
        x = vals2;
        xL = min(x);xH = max(x);
        xPrior = randraw('norm', [th1 th2], [1e5 1] );
        xL = min([xPrior; vals2(:,ii)]);
        xH = max([xPrior; vals2(:,ii)]);
        %Plot the prior
        pEdges = linspace(xL,xH,1000);
        [yP,xP] = histcounts(xPrior,pEdges);
        [xP,yP,~] = histtoplot(yP,xP,50);
        plot(xP,yP);
        %Plot the posterior distribution
        [yPost,xPost] = histcounts(vals2(:,ii),pEdges); 
        [xPost,yPost,~] = histtoplot(yPost,xPost,20);
        hold on 
        plot(xPost,yPost)
        if ii==6
            line([-.255 -.255],[0 14],'Color','k','LineWidth',2)
        end
        if ii==7
            line([.261 .261],[0 20],'Color','k','LineWidth',2)
        end
        title(colLabels(offsetsI(ii)))
        xlim([quantile(vals2(:,ii),0.001) quantile(vals2(:,ii),0.999)])
        xlabel("W/m^{2}/decade")
        set(gca,'ytick',[])
        set(gca,'FontSize',fSize)
    end
    saveas(gcf,['plots/priorposterior2_' datesave '.png'])
    %------------------------------------------------------------------
    % Last, plot noise estimates
    %------------------------------------------------------------------
    figure2('Position',[10 10 1600 600])
    offI = [29:38]';
    numPlots = length(varNames3);
    for ii = 1:numPlots
        subplot(3,4,ii)
        hold on
        th1 = th1M3(ii);
        th2 = th2M3(ii);
        %Plot the prior distribution
        x = vals3;
        xL = min(x);xH = max(x);
        %From randraw.m, form of drawing is y = randraw('gamma', [1 3 2], [1e5 1] );
        %where inputs are location, shape, and scale parameters,
        %respectively. Second input is number of draws
        if strcmp(distType3(ii),"normal")
            xPrior = randraw('norm', [th1 th2], [1e5 1] );
            xPrior(xPrior < 0) = NaN;
        else
            xPrior = drawgamma(0, th1, th2, 1e5);
        end
        xPrior = xPrior.^0.5;
        xL = max([min([xL(ii); min(xPrior)]); 1e-5]);
        xH = max([xH(ii); xPrior]);
        %Plot the prior
        pEdges = logspace(log10(xL),log10(xH),1000);
        [yP,xP] = histcounts(xPrior,pEdges);
        [xP,yP,~] = histtoplot(yP,xP,50);
        plot(xP,yP);
        %Plot the posterior distribution
        [yPost,xPost] = histcounts(vals3(:,ii),pEdges); 
        [xPost,yPost,~] = histtoplot(yPost,xPost,50);
        hold on
        plot(xPost,yPost)
        title(varNames3(ii))
        xmin = max([min([quantile(vals3(:,ii),0.001); quantile(xPrior,0.001)]); 1e-5]);
        xmax = max([quantile(vals3(:,ii),0.999); quantile(xPrior,0.999)]);
        xlim([xmin xmax])
        set(gca,'ytick',[])
        set(gca,'Xtick',linspace(xmin,xmax,4))
        xtickformat('%.2f')
        xlabel("W/m^{2}")
        set(gca,'FontSize',fSize)
    end
    saveas(gcf,['plots/priorposterior3_' datesave '.png'])
end
if plotData 
    showTrend = 1;
    smoothWindow = 6; %set smoothing (months)
    %Create x-axis points for cycle demarcation
    c21=datejd([dateStruct.cycles(1,:), fliplr(dateStruct.cycles(1,:))]);
    c23=datejd([dateStruct.cycles(3,:), fliplr(dateStruct.cycles(3,:))]);
    c25=datejd([dateStruct.cycles(5,:), fliplr(dateStruct.cycles(5,:))]);
    
    figure2('Position',[10 10 1000 1200])
    subplot('position',[.09 .53 .82 .45]) %Plot of satellite observations
    
     fill(c21,[1360 1360 1374 1374],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c23,[1360 1360 1374 1374],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c25,[1360 1360 1374 1374],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    ind = 1;
    for ii = 3:numel(lI) %Iterate over satellite observations
        if showTrend
            hold on
            [trends,offsets2]= returntrend(A,tau,ii);
            tM = mean(trends,2) + mean(offsets2)+offsets(ii);
            t995 = quantile(trends,.995,2)+quantile(offsets2,.995)+offsets(ii);
            t005 = quantile(trends,.005,2)+quantile(offsets2,.005)+offsets(ii);
            x2 = [dateM(oM(:,ii))', flip(dateM(oM(:,ii)))'];
            fill(x2,[t995(oM(:,ii))',flip(t005(oM(:,ii)))'], ...
                [1 .85 .85],'FaceAlpha',0.5,'LineStyle','none');
            hold on
            plot(dateM(oM(:,ii)),tM(oM(:,ii)),'Color','r')
            hold on
        end
        plot(dateMAll,valMAll(:,ii)+offsets(ii),'o','MarkerSize',4,...
            'Color',c(ind,:));
        hold on
        hh(ind) = plot(dateM,valM(:,ii)+offsets(ii),'.','MarkerSize',10,...
            'Color',c(ind,:));
        ind = ind + 1;
    end
    xlim([datetime(1978,1,1) datetime(2022,1,1)])
    ylim([1360 1374])
    legend(hh,colLabels(3:end),'NumColumns',2)
    legend boxoff
    xlabel('Year')
    ylabel('TSI (W/m^{2})')
    set(gca,'FontSize',fSize)
   
    
    subplot('position',[.09 .125 .82 .33]) %Plot of proxy observations
     yyaxis left
    fill(c21,[0 0 350 350],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c23,[0 0 350 350],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c25,[0 0 350 350],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    h(1)=plot(dateM,valM(:,1)+offsets(1),'.','MarkerSize',10);
    ylabel('Sunspot number')
    yyaxis right
    h(2)=plot(dateM,valM(:,2)+offsets(2),'.','MarkerSize',10);
    legend(h,'Silso sunspot number','Mg-II')
    legend boxoff
    ylabel('Mg-II index')
    xlabel('Year')
    xlim([datetime(1978,1,1) datetime(2022,1,1)])
    yyaxis left
    ylim([0 350])
    set(gca,'FontSize',fSize)
    
   
    
    
  
    saveas(gcf,'plots/tsicompare_23_2_28.png')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if table1
    load oTSI_23_02_01.mat
    
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
    
    dateM=dateshift(dateM,'start','month');
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
    [a,m,tau]=reconstructiontrends(datePMOD,tsiPMOD,dates);
    CPMDF=[meanPMOD(end)-meanPMOD(end-1);a;m;tau];
    
    [a,m,tau]=reconstructiontrends(dateSolid,tsiSolid,dates);
    Comm_Con=[meanSolid(end)-meanSolid(end-1);a;m;tau];
    
    [a,m,tau]=reconstructiontrends(dateROB,tsiROB,dates);
    ROB=[meanROB(end)-meanROB(end-1);a;m;tau];
    
    [a,m,tau]=reconstructiontrends(dateNRL,tsiNRL,dates);
    NRL=[meanNRL(end)-meanROB(end-1);a;m;tau];
    
    [a,m,tau]=reconstructiontrends(dateSAT,tsiSAT,dates);
    SAT=[meanSAT(end)-meanSAT(end-1);a;m;tau];
    
    [a,m,tau]=reconstructiontrends(dateAR6,tsiAR6,dates);
    AR6=[NaN;a;m;tau];
    
    %BTSI
    for ii=1:1000
        interval=floor(size(xAll,2)./1000);
        [a,m,tau]=reconstructiontrends(dateM,xAll(:,ii.*interval),dates);
        aa(ii)=a;mm(ii)=m;tt(ii)=tau;
    end
    
    BTSI=[meanMeanBTSI(end)-meanMeanBTSI(end-1);prctile(aa,50);prctile(mm,50);prctile(tt,50)];
    BTSI25=[prctile(meanBTSI(:,end)-meanBTSI(:,end-1),2.5);prctile(aa,2.5);prctile(mm,2.5);prctile(tt,2.5)];
    BTSI975=[prctile(meanBTSI(:,end)-meanBTSI(:,end-1),97.5);prctile(aa,97.5);prctile(mm,97.5);prctile(tt,97.5)];
    
    reconstructions=table(CPMDF,Comm_Con,ROB,NRL,SAT,AR6,BTSI,BTSI25,BTSI975)
end
if table2
    clear tabout
    offsets([3;4;5;6;7;8;9;10;11;12])=offsets([3;4;5;6;7;8;9;10;11;12])-offsets(9);
    table2order=[3;4;5;6;7;8;9;10;11;12;1;2];
    T0=outDat.T0(lI(table2order));
    th0=outDat.th0(lI(table2order))';
    e0=sqrt(th0./(T0'-1));
    for ii=1:length(colLabels)
        ind=table2order(ii);
        Asub=squeeze(A(ind,:,:));
        tabout(:,ii)=[offsets(ind);prctile(Asub(1,:).*scaling(ind),[2.5 50 97.5])'+offsets(ind);...
            outDat.H0(ind,2).*scaling(ind);prctile(Asub(2,:).*scaling(ind),[2.5 50 97.5])';...
            prctile(Asub(3,:),[2.5 50 97.5])';...
            e0(ii)./outDat.H0(ind,2);prctile(sqrt(sigY(ind,:))./Asub(2,:),[2.5 50 97.5])'];
    end
    rows=["a_0";"mean 2.5";"mean 50";"mean 97.5";"b_0";...
        "scaling 2.5";"scaling 50";"scaling 97.5";...
        "drift 2.5";"drift 50";"drift 97.5";"\epsilon_0";"error 2.5";"error 50";"error 97.5"];
    displayTable=array2table(tabout,'VariableNames',colLabels(table2order),'RowNames',rows)
end
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
    load ar1_22_11_04.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"AR(1)");
    
    %Get AR(3) reconstruction
    load ar3_22_11_04.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"AR(3)");
    
    %Get reconstruction keeping satellite artifacts
    load ar2_22_11_03_all.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"No artifact removal");
    
    %Get reconstruction scaling proxies with 3 most recent cycles
    load ar2_22_11_03_3cycle.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"3 cycle regression");
    
    %Get reconstruction using PMOD corrections
    load ar2_22_11_04_pmodcorrections.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"PMOD satellite corrections");
    
    
    %Get reconstruction using VIRGO baseline
    load ar2_22_11_04_VIRGO.mat
    
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"VIRGO baseline");
    
    %Get reconstruction using satellite-only, drift allowed
    load ar2_22_11_04_satdrift.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"Sat. only w/ drift");
    
    %Get reconstruction using satellite-only, no drift allowed
    load ar2_22_11_04_satnodrift.mat
    [min95,amp95,total95]=getaltstats(dateM(intI),xAll(intI,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"Sat. only w/o drift");
    
    %Get reconstruction using proxy-model
    eI=sum(isnan(xM),2)==0;
    xMM=xM(eI,:);dateMM=dateM(eI);
    intII=intI(eI);
    [min95,amp95,total95]=getaltstats(dateMM(intII),xMM(intII,:));
    [min,amp,total,model]=appendaltstats(min,amp,total,model,min95,amp95,total95,"Proxy model");
    
    %Make table from results
    compare=table(model,amp,min,total);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if uncBTSI
    alpha=0.95; %CI width for analysis
    tsix = quantile(xAll,[(1-alpha)/2,1-(1-alpha)/2],2);
    figure
    plot(dateM,tsix(:,2)-tsix(:,1))
end
if solarvsanthroforcing
    stYr=1978;endYr=2021;
    %Calculate contribution from BTSI solar forcing
    iValid = dateM.Year>=stYr&dateM.Year<=endYr;
    xAll=xAll(iValid,:)'+offsets(9);
    dateM=dateM(iValid);
    X = [ones(size(dateM,1),1) juliandate(dateM)];
    for ii=1:size(xAll,1)
        b(:,ii)=regress(xAll(ii,:)',[ones(size(dateM,1),1) juliandate(dateM)]);
        yB(:,ii)=[ones(size(dateM,1),1) juliandate(dateM)]*b(:,ii);
    end
    fBTSI= b(2,:)*(juliandate(dateM(end))-juliandate(dateM(1))); %BTSI TSI change 1978-2021
    %IPCC AR6 conversion to ERF (Ch.7 pg. 46):
    EFF=0.72*(1-0.29)./4;
    erfBTSI=fBTSI.*EFF;
    
    %Calculate SOLARIS-HEPPA forcing
    load oTSI_22_10_27.mat
    dateHEPPA=oTSI(9).mthdatetime;
    tsi=oTSI(9).mthtsi;
    iHEPPA=dateHEPPA>=dateM(1)&dateHEPPA<=dateM(end);
    dH=dateHEPPA(iHEPPA);
    tH=tsi(iHEPPA);
    b=regress(tH,[ones(size(dH,1),1) juliandate(dH)]);
    fSOLID = b(2)*(juliandate(dH(end))-juliandate(dH(1)));
    erfS=fSOLID.*EFF;
    
    %Calculate ERF from Anthropogenic forcing, via loadallforcing.m and
    %loadanth.m
    load forcingSet.mat
    iSt=f(2).yr==stYr;
    iEnd=f(2).yr==endYr;
    erfB=f(2).preds(iEnd,:)-f(2).preds(iSt,:);
    
    %Convert to temperature response using 1.8 K per 3.7 W/m^{2} TCR from
     %Ribes et al. 2021
     TS=erfS.*(1.8/3.7);
     TB=erfBTSI.*(1.8/3.7);
end

%  if radForcingEffect
%      load oTSI_22_09_20.mat
%      dateSolid=oTSI(9).mthdatetime;
%      tsiSolid=oTSI(9).mthtsi;
%      dateSolid=dateshift(dateSolid,'start','month');
%      tsiSolid=meaninterval(dateSolid,tsiSolid,1990,2010);
% %      sInd=find(dateSolid >= datejd(dates(1)) & dateSolid <= datejd(dates(end)));
% %      bSolid=regress(tsiSolid(sInd),[ones(length(sInd),1) (juliandate(dateSolid(sInd))-mean(juliandate(dateSolid(sInd))))./(365.25*10)]);
%      [~,~,tSolid]=reconstructiontrends(dateSolid,tsiSolid,dates);
% %      bInd=find(dateM >= datejd(dates(1)) & dateM <= datejd(dates(end)));
% %      bBTSI=regress(mean(xAll(bInd,:),2),[ones(length(bInd),1) (juliandate(dateM(bInd))-mean(juliandate(dateM(bInd))))./(365.25*10)]);
%      [~,~,tBTSI]=reconstructiontrends(dateM,mean(xAll,2),dates);
%      %Convert to ERF using IPCC AR6 Conversion factor of 1/7.8 (Pg. 957)
%      erfS=tSolid.*(diff(dates)./(365.25*10))./7.8;
%      erfB=tBTSI.*(diff(dates)./(365.25*10))./7.8;
%      
%      %Convert to temperature response using 1.8 K per 3.7 W/m^{2} TCR from
%      %Ribes et al. 2021
%      TS=erfS.*(1.8/3.7);
%      TB=erfB.*(1.8/3.7);
%  end


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


