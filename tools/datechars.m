function [yrchar,mthchar,daychar,pthDate] = datechars
%DATECHARS Return two digit character arrays for year, month, day, given
%current date.
% Ted Amdur
% 2023/04/22

currD=datetime;
yr=num2str(year(datetime));
mth=num2str(month(datetime));
if month(datetime)<10
    mth=['0' mth];
end
dy=num2str(day(datetime));
if day(datetime)<10
    dy=['0' dy];
end
yrchar=yr(3:4);
mthchar=mth;
daychar=dy;
pthDate=[yrchar '_' mthchar '_' daychar];
end

