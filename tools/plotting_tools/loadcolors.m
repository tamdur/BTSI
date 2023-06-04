function [cAll,cCI1,cCI2] = loadcolors
%LOADCOLORS Load colorscheme for BTSI-Long plots

%First, color scheme for labels, up to 10
c1 = [51; 136; 68; 17; 153; 221; 204; 136; 170];
c2 = [34; 204; 170; 119; 153; 204; 102; 34; 68];
c3 = [136; 238; 153; 51; 51; 119; 119; 85; 153];
c = [c1 c2 c3]; c = c./255;
c(10,:) = c(1,:);
cAll=c;

%Next, CI color scheme 1 (95% then 80%)
c95=[.65 .65 .65];
c80=[.4 .4 .4];
cCI1=[c95;c80];
    
    
%Next, CI color scheme 2 (95% then 80%)
c95=[250 188 188];
c80=[250 142 142];
cCI2=[c95;c80];
cCI2=cCI2./255;

end

