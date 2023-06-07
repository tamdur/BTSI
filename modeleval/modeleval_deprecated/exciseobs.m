function [valE,oE,iX]=exciseobs(valM,oM,col,fracEx)

L=sum(oM(:,col)); %Number of observations
oI=find(oM(:,col)); %indices of observations
rI=ceil(rand.*L); %select a random start index to cut out
nani=mod(rI:rI+ceil(fracEx.*L),L)+1; %The indices to excise
%nani=nani(nani~=0); %Remove zero indexing
iX=oI(nani); %Map excisions to the locations of observations
valM(iX,col)=NaN; %Remove excisions from valM
oM(iX,col)=false; %Remove excisions from oM
valE=valM;
oE=oM;
end