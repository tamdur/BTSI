function [low21,low22,low23,low24,lowAll] = getlowspotindex(dateM,valM,offsets,colLabels,sx)
%GETLOWSPOTINDEX Return indices for months with low solar activity, subset
%by solar cycle

%first, select for datetimes when SILSO spot count is less than sx
    sptI=find(strcmp(colLabels,"SILSO")); %Get column of SILSO observations
    lowAll=(valM(:,sptI)+offsets(sptI))<sx;
    dateS=getdates;
    low24=and(lowAll,dateM.Year<2025 & dateM.Year > 2015);
    low23=and(lowAll,dateM.Year<2015 & dateM.Year > 2003);
    low22=and(lowAll,dateM.Year<2000 & dateM.Year > 1990);
    low21=and(lowAll,dateM.Year<1990 & dateM.Year > 1980);
end

