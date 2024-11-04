function yearfraction = yearfraction(startDate,endDate)
%YEARFRACTION returns a fraction, in years, based on the number of days 
%between dates StartDate and EndDate using the given day-count Basis.
endDateMinus1 = datetime(endDate.Year,startDate.Month,startDate.Day);
yd = endDate.Year - startDate.Year;
yearfraction = yd + caldays(between(endDateMinus1,endDate,'days'))./365;
end

