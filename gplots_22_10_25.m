%Plots for draft of TSI reconstruction manuscript. Built to take the 
%ar2 formulation
%
% Ted Amdur
% 10/25/22 

clearvars



tsiComparison = 1; %Set to 1 to plot Figure 1 of manuscript
cycleMin = 0; %Set to 1 to plot Figure 2 of manuscript
tsiSunspots=0; %Set to 1 to plot Figure 3 of manuscript

fSize = 16;

load ar2_22_11_03.mat %from runchain_22_04_25.m
load(outDat.obsmatrix); %From makeobsmatrix.m
excludeFliers=outDat.excludeFliers;
dateStruct = getdates;
startDate = datejd(dateStruct.all(1));
endDate = datejd(dateStruct.all(2));
%dates = [juliandate(datetime(1978,11,1)) juliandate(datetime(2021,12,1))]; %Use full interval
dates=[juliandate(datetime(1980,2,1)) juliandate(datetime(2021,12,31))]; %Use shared interval
%Load data, with colLabels corresponding to observer source for each column
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
t=t(:,lI);
valM=valM(:,lI); valMAll=valMAll(:,lI);


%Colormap for plots
c1 = [51; 136; 68; 17; 153; 221; 204; 136; 170];
c2 = [34; 204; 170; 119; 153; 204; 102; 34; 68];
c3 = [136; 238; 153; 51; 51; 119; 119; 85; 153];
c = [c1 c2 c3]; c = c./255;
c(10,:) = c(1,:);

if tsiComparison 
    showTrend = 0;
    smoothWindow = 6; %set smoothing (months)
    %Create x-axis points for cycle demarcation
    c21=datejd([dateStruct.cycles(1,:), fliplr(dateStruct.cycles(1,:))]);
    c23=datejd([dateStruct.cycles(3,:), fliplr(dateStruct.cycles(3,:))]);
    c25=datejd([dateStruct.cycles(5,:), fliplr(dateStruct.cycles(5,:))]);
    
    figure2('Position',[10 10 1000 1200])
    subplot('position',[.09 .85 .85 .13]) %Plot of proxy observations
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
    xlim([datetime(1978,1,1) datetime(2022,1,1)])
    yyaxis left
    ylim([0 350])
    set(gca,'FontSize',fSize)
    
    subplot('position',[.09 .575 .85 .25]) %Plot of satellite observations
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
            [trends,offsets2]= returntrend(A,t,ii);
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
    ylabel('TSI (W/m^{2})')
    set(gca,'FontSize',fSize)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot of reconstructions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot('position',[.09 .04 .85 .50]) 
    fill(c21,[-1.2 -1.2 1.8 1.8],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c23,[-1.2 -1.2 1.8 1.8],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c25,[-1.2 -1.2 1.8 1.8],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
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
    
%     %Get satellite only reconstruction
    load ar2_22_11_04_satnodrift.mat
    xms = mean(xAll,2);
    xms = smoothPH(xms,smoothWindow);
    [~,xmsO] = meaninterval(dateMAll,xms,1990,2010);
    xms = xms-xmsO; 
    
%     %Get MLR model reconstruction
    xmmlr=mean(xMLR,2);
    xmmlr=smoothPH(xmmlr,smoothWindow);
    [~,xmmlrO] = meaninterval(dateMAll,xmmlr,1990,2010);
    xmmlr = xmmlr-xmmlrO;
    
    %Plot other TSI reconstructions
    ind = 1;
    load oTSI_22_10_27.mat
    %Plot PMOD
    hold on
    tsiA = smoothPH(oTSI(7).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(7).mthdatetime,tsiA,1990,2010);
    h(ind) = plot(oTSI(7).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(7,:));
    legendtxt(ind) = string(oTSI(7).product);
    ind = ind + 1;
    
    %Plot SOLID
    hold on
    tsiA = smoothPH(oTSI(9).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(9).mthdatetime,tsiA,1990,2010);
    h(ind) = plot(oTSI(9).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(8,:));
    legendtxt(ind) = "Comm.-Consensus";
    ind = ind + 1;
    
    %Plot NRLTSI2
    hold on
    tsiA = smoothPH(oTSI(4).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(4).mthdatetime,tsiA,1990,2010);
    h(ind) = plot(oTSI(4).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(5,:));
    legendtxt(ind) = string(oTSI(4).product);
    ind = ind + 1;
    

    %Plot SATIRE-S during whole interval
    hold on
    tsiA = smoothPH(oTSI(5).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(5).mthdatetime,tsiA,1990,2010);
    plot(oTSI(5).mthdatetime,tsiA-tsiAO,'--','LineWidth',2.5,...
        'Color',c(4,:));
    legendtxt(ind)=string(oTSI(5).product);
    hold on %Plot SATIRE-S again during period of high-quality MDI/HMI data
    dateMDIHMI=oTSI(5).mthdatetime(oTSI(5).mthdatetime>datetime(1999,02,02));
    h(ind) = plot(dateMDIHMI,tsiA(oTSI(5).mthdatetime>datetime(1999,02,02))-tsiAO,'LineWidth',2.5,...
        'Color',c(4,:));
    
    
    
    ind = ind + 1; %Index our contribution
    hold on
    h(ind) = plot(dateM,tsix(:,3),'LineWidth',3.5,'Color','k');
    legendtxt(ind) = "BTSI";
    ind = ind + 1;
    hold on
    h(ind)=plot(dateMAll,xms,'LineWidth',2.5,'Color',c(1,:));
    legendtxt(ind)="BTSI sat. only";
    ind=ind+1;
    hold on
    h(ind)=plot(dateMAll,xmmlr,'LineWidth',2.5,'Color',c(2,:));
    legendtxt(ind)="BTSI proxy model";

    legend(h,legendtxt,'NumColumns',2)
    legend boxoff
    set(gca,'FontSize',fSize)
    xlabel('Year')
    ylabel('TSI anomaly from 1990-2010 mean (W/m^{2})')
    xlim([datetime(1980,2,1) datetime(2022,1,1)])
    ylim([-0.9 1.25])
    text(datejd(dateStruct.cycles(1,1))+years(4.5),-0.8,'Cycle 21','FontSize',14)
    text(datejd(dateStruct.cycles(2,1))+years(3.25),-0.8,'Cycle 22','FontSize',14)
    text(datejd(dateStruct.cycles(3,1))+years(4.75),-0.8,'Cycle 23','FontSize',14)
    text(datejd(dateStruct.cycles(4,1))+years(4.25),-0.8,'Cycle 24','FontSize',14)
    saveas(gcf,'plots/tsicompare_23_2_28.png')
end
if cycleMin
    %First, calculate the percentage of time that each simulation spends
    %inside of the 95% confidence interval

    xAll=xAll'+offsets(9);xms=mean(xAll,1);
    figure2('Position',[10 10 1500 900])
    subplot('position',[.09 .31 .85 .68])
    trendInd=dateM>=datejd(dates(1))&dateM<datejd(dates(2));
    X = [ones(size(dateM(trendInd),1),1) juliandate(dateM(trendInd))];
    warning('off','MATLAB:rankDeficientMatrix') %Some realizations for quant_reg are rank deficient
    for ii=1:size(xAll,1)
        b5(:,ii)=quant_reg(X,xAll(ii,trendInd)',0.05);
        b95(:,ii)=quant_reg(X,xAll(ii,trendInd)',0.95);
        yhat5(:,ii) = X*b5(:,ii);
        yhat95(:,ii) = X*b95(:,ii);
    end
    hold on
    yplot5=prctile(yhat5',[.5 5 50 95 99.5])';
    yplot95=prctile(yhat95',[.5 5 50 95 99.5])';
    x2 = [dateM(trendInd)', fliplr(dateM(trendInd)')];
    fill(x2,[yplot5(:,1)',fliplr(yplot5(:,5)')],c(2,:).*1.05,'FaceAlpha',...
        0.5,'LineStyle','none');
    fill(x2,[yplot5(:,2)',fliplr(yplot5(:,4)')],c(2,:).*.9,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    plot(dateM(trendInd),yplot5(:,3),'Color',c(2,:).*.9,'LineWidth',2)
    hold on
    fill(x2,[yplot95(:,1)',fliplr(yplot95(:,5)')],c(8,:).*1.05,'FaceAlpha',...
        0.5,'LineStyle','none');
    fill(x2,[yplot95(:,2)',fliplr(yplot95(:,4)')],c(8,:).*.9,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    plot(dateM(trendInd),yplot95(:,3),'Color',c(8,:),'LineWidth',2)
    hold on
    plot(dateM,xAll','.','Color',[0.65 0.65 0.65])
    hold on
    plot(dateM,xms,'Color','k','LineWidth',2)
    xlabel('Year')
    ylabel('TSI (W/m^{2})')
    set(gca,'FontSize',fSize)
    box off
    subplot('position',[.09 .07 .85 .18])
    histogram(b5(2,:).*(365.25*10),linspace(-.40,0.02,85),'FaceColor',c(2,:))
    hold on
    histogram(b95(2,:).*(365.25*10),linspace(-.40,0.02,85),'FaceColor',c(8,:))
    hold on
    histogram((b95(2,:)-b5(2,:)).*(365.25*10),linspace(-.40,0.02,85),'FaceColor',[0.6 0.6 0.6])
    ylabel('Realizations')
    xlabel('W/m^{2}/decade')
    xlim([-0.35 0])
    set(gca,'FontSize',fSize)
    saveas(gcf,'plots/mincycle_22_11_07.png')
end
if tsiSunspots
    xAll=xAll-nanmean(xAll(:));xMLR=xMLR-nanmean(xMLR(:));
    xm = nanmean(xAll,2);
    xmlrm=nanmean(xMLR,2);
    yM=valM(:,1)+offsets(1);
    nanI=isnan(yM);
    y=yM(~nanI);x=xm(~nanI);
    %Determine piecewise regression for tsi-sunspot relationship
    %Code derived from following post:
    %https://www.mathworks.com/matlabcentral/answers/344059-piece-wise-linear-fitting
    %Assume the curve is defined in terms of two linear segments, with break at b1
    plusfun = @(x) max(x,0);
    model = @(P,x) P(1) + P(2)*x + (P(3)-P(2))*plusfun(x-P(4));
    %Parameter set P includes: offset, slope 1, slope 2, break value
    lb=[-50,-500,-500,min(x)+1E-8];
    ub=[50,500,500,max(x)-1E-8];
    P0=(lb+ub)./2; %naive, unbiased initial search location
    [Pfit,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(model,P0,x,y,lb,ub);
    conf = nlparci(Pfit,residual,'jacobian',jacobian);
    [Ypred,delta] = nlpredci(model,x,Pfit,residual,'Jacobian',jacobian);
    modelpred=model(Pfit,x);
    lowpred=modelpred+delta;
    highpred=modelpred-delta;
    [linfit,linint]=regress(y,[ones(length(x),1) x]);
    [~,sI]=sort(x); %Plot line in order on the x-axis
    
    figure2('Position',[10 10 1000 1000])
    x2=[x(sI); flipud(x(sI))];
    
    %Plot CI from BTSI for sunspots
    spotPred=zeros(size(x,1),size(xAll,2))';
    for iR=1:size(xAll,2)
        spotPred(iR,:)=squeeze(A(1,[2 3],iR))*[x';t(~nanI,1)']+offsets(1);
    end
    linPreds=prctile(spotPred,[2.5 50 97.5])';
    linlow=linPreds(:,1);
    linhigh=linPreds(:,3);
    
    fill(x2,[linlow(sI); flipud(linhigh(sI))], [.85 .85 .85],'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[lowpred(sI); flipud(highpred(sI))],[1 .85 .85],'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    
    
    
    cColor = colororder;
    for ii = 1:size(dateStruct.cycles,1)
        cycleI = dateM > datejd(dateStruct.cycles(ii,1)) & dateM < datejd(dateStruct.cycles(ii,2));
        hold on
        h(ii) = plot(xm(cycleI),valM(cycleI,1)+offsets(1),'.','MarkerSize',15,'Color',c(2*ii-1,:));
    end
    hold on
    h(ii+1)=plot(x(sI),linPreds(sI,2),'Color','k','LineWidth',1.5);
    hold on
    h(ii+2)=plot(x(sI),modelpred(sI),'Color','r','LineWidth',1.5);
    legend(h,'Cycle 21','Cycle 22','Cycle 23','Cycle 24','Cycle 25','BTSI prediction','Piecewise prediction','Location','NorthWest')
    legend boxoff
    hold on
    xlim([-.8 1.35])
    ylim([0 300])
    xlabel('BTSI anomaly (W/m^{2})')
    ylabel('SILSO sunspot number')
    set(gca,'FontSize',fSize)
    saveas(gcf,'plots/tsispots_22_11_07.png')
    
end





