function [xMLR,b] = makexmlr(dateM,sindex,oindex,valM,xAll,offsets)
%MAKEXMLR 


xMLR=NaN(size(valM,1),size(xAll,2));
for ii=1:size(xAll,2)
    dateI=dateM.Year >= 2003 & dateM.Year <= 2014;
    pI=find(sindex);
    p=[ones(size(valM,1),1) valM(:,pI(1))-min(valM(:,pI(1)),[],'omitnan')...
        valM(:,pI(2))-min(valM(:,pI(2)),[],'omitnan')];
%     pall=[ones(size(valM,1),1) valMAll(:,pI(1))-min(valMAll(:,pI(1)),[],'omitnan')...
%         valMAll(:,pI(2))-min(valMAll(:,pI(2)),[],'omitnan')];
    [b,bint]=regress(xAll(dateI,ii)+offsets(~oindex), p(dateI,:),0.6827);
    b=b+randn(3,1).*((bint(:,2)-bint(:,1))/2);%include uncertainty in coefficient estimate
    xMLR(:,ii)=p*b;
end
end

