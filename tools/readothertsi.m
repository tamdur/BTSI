%Load alternative TSI reconstructions
%Ted Amdur
%12/18/21, updated 07/13/22 to use new PMOD
%Updated 23/02/01 to show both SOLID
%Updated 24/11/03 to use updated files within BTSI repos

clearvars
%Read SOLARIS-HEPPA
filePath = 'solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA-3-2_gn_185001-229912.nc';
tsi = ncread(filePath,'tsi');
t = ncread(filePath,'time');
time = datetime(1850,1,1) + days(t);
%Remove extrapolations not based in observation
tsiFull = tsi;
tsi(time.Year >= 2015) = NaN;

oTSI(1).product = 'SOLARIS-HEPPA';
oTSI(1).path = filePath;
oTSI(1).version = '3.2';
oTSI(1).mthtsi = tsi;
oTSI(1).mthtsifull = tsiFull;
oTSI(1).mthdatetime = time;

% Read Dudok de wit et al.
importdudokdewit
oTSI(2).product = 'SOLID';
oTSI(2).path = 'dudok_dewit_accessed_22_1_18.dat';
oTSI(2).version = '1.1';
oTSI(2).tsi = TSIcomposite.TSIc;
oTSI(2).tsiUnc = TSIcomposite.dTSIc;
oTSI(2).datetime = datejd(TSIcomposite.julian_dat);

[mthVal,mthDate] = dailytomonthly(TSIcomposite.julian_dat,TSIcomposite.TSIc);
mthVal2 = dailytomonthly(TSIcomposite.julian_dat,TSIcomposite.dTSIc);
oTSI(2).mthtsi = mthVal;
oTSI(2).mthtsiUnc = mthVal2;
oTSI(2).mthdatetime = mthDate;

%Read PMOD Composite
importPMOD
tsi = PMOD.TSI_new_;
tsi(tsi<1000) = NaN;
t = PMOD.epoch1980;
time = datetime(1980,1,1) + days(t);

oTSI(3).product = 'PMOD';
oTSI(3).path = 'PMOD_accessed_22_1_18.dat';
oTSI(3).version = '42_65_1805';
oTSI(3).tsi = tsi;
oTSI(3).datetime = time;

[mthVal,mthDate] = dailytomonthly(juliandate(time),tsi);
oTSI(3).mthtsi = mthVal;
oTSI(3).mthdatetime = mthDate;

%Read NRLTSI2
nrlPth = 'tsi_v02r01_monthly_s188201_e202312_c20240223_accessed_24_10_23.nc';
tsi = ncread(nrlPth,'TSI');
time = datetime(1610,1,1+ncread(nrlPth,'time'));
oTSI(4).product = 'NRLTSI2';
oTSI(4).path = nrlPth;
oTSI(4).version = 'v02r01';
oTSI(4).mthtsi = tsi;
oTSI(4).mthtsiUnc=ncread(nrlPth,'TSI_UNC');
oTSI(4).mthdatetime = time;

%Read SATIRE
importsatire
time = satiretsi.JD;
tsi = satiretsi.TSI;
tsiUnc=satiretsi.TSI_H-satiretsi.TSI_L;

oTSI(5).product = 'SATIRE-S';
oTSI(5).path = 'SATIRE-S_TSI_latest_accessed_24_10_23.txt';
oTSI(5).version = '20220923';
oTSI(5).tsi = tsi;
oTSI(5).datetime = datejd(time);

[mthVal,mthDate] = dailytomonthly(time,tsi);
mthVal2 = dailytomonthly(time,tsiUnc);
oTSI(5).mthtsi = mthVal;
oTSI(5).mthdatetime = mthDate;
oTSI(5).mthtsiUnc=mthVal2;

%Read ACRIM
acrimcompimport
importdudokdewit
oTSI(6).product = 'ACRIM composite';
oTSI(6).path = 'nnaa3Acrimcomposite_accessed_22_02_11.txt';
oTSI(6).version = 'NNAA3';
oTSI(6).tsi = acrimcomp.tsi;
oTSI(6).datetime = datetime(floor(acrimcomp.yearf),1,1)+...
    years(acrimcomp.yearf-floor(acrimcomp.yearf));
[mthVal,mthDate] = dailytomonthly(juliandate(oTSI(6).datetime),oTSI(6).tsi);
oTSI(6).mthtsi = mthVal;
oTSI(6).mthdatetime = mthDate;

%Read 2022 iteration of PMOD (CPMDF)
importPMOD2022
oTSI(7).product='PMOD CPMDF';
oTSI(7).path='MergedPMOD_NobaselineScaleCycle23_JPM_September2024_accessed_24_10_23.txt';
oTSI(7).version='CPMDF2';
oTSI(7).tsi=PMOD2022.TSIAfterWavelet;
oTSI(7).datetime=datejd(PMOD2022.JD);
[mthVal,mthDate] = dailytomonthly(juliandate(oTSI(7).datetime),oTSI(7).tsi);
mthVal2 = dailytomonthly(juliandate(oTSI(7).datetime),PMOD2022.TSIUnc2);
oTSI(7).mthtsi = mthVal;
oTSI(7).mthdatetime = mthDate;
oTSI(7).mthtsiUnc=mthVal2;

%Read ROB TSI composite
importROB
oTSI(8).product='ROB';
oTSI(8).path='ROB_TSI_composite_latest_accessed_24_10_23.txt';
oTSI(9).source='https://www.sidc.be/observations/space-based-timelines/tsi';
oTSI(8).version='2022_02_28';
oTSI(8).tsi=ROB.TSI;
oTSI(8).datetime=datejd(ROB.JD);
[mthVal,mthDate] = dailytomonthly(juliandate(oTSI(8).datetime),oTSI(8).tsi);
oTSI(8).mthtsi = mthVal;
oTSI(8).mthdatetime = mthDate;

% Read Dudok de wit et al. updated
importdudokdewitupdated
oTSI(9).product='SOLID (Kopp Updates)';
oTSI(9).path='TSI_Composite-SIST_accessed_24_10_23.txt';
oTSI(9).source='https://spot.colorado.edu/~koppg/TSI/TSI_Composite-SIST.txt';
oTSI(9).version='2023_08_11';
oTSI(9).tsi=SOLID2.tsi;
oTSI(9).datetime=datejd(SOLID2.JD);
[mthVal,mthDate] = dailytomonthly(juliandate(oTSI(9).datetime),oTSI(9).tsi);
oTSI(9).mthtsi = mthVal;
oTSI(9).mthdatetime = mthDate;

% Read Dudok de wit et al. uncorrected
importdudokdewit
oTSI(10).product = 'SOLID Uncorrected';
oTSI(10).path = 'dudok_dewit_accessed_22_1_18.dat';
oTSI(10).version = '1.1';
oTSI(10).tsi = TSIcomposite.TSIc;
oTSI(10).tsiUnc = TSIcomposite.eTSIo;
oTSI(10).datetime = datejd(TSIcomposite.julian_dat);

[mthVal,mthDate] = dailytomonthly(TSIcomposite.julian_dat,TSIcomposite.TSIc);
mthVal2 = dailytomonthly(TSIcomposite.julian_dat,TSIcomposite.TSIc);
oTSI(10).mthtsi = mthVal;
oTSI(10).mthtsiUnc = mthVal2;
oTSI(10).mthdatetime = mthDate;

%Last, load AR7 recontruction
filePath = 'multiple_input4MIPs_solar_CMIP_SOLARIS-HEPPA-CMIP-4-4_gn_185001-202312_accessed_24_11_03.nc';
tsi = ncread(filePath,'tsi');
t = ncread(filePath,'time');
time = datetime(1850,1,1) + days(t);
%Remove extrapolations not based in observation
tsiFull = tsi;
tsi(time.Year >= 2024) = NaN;

oTSI(11).product = 'SOLARIS-HEPPA-CMIP-4-4';
oTSI(11).path = filePath;
oTSI(11).version = '4.4';
oTSI(11).mthtsi = tsi;
oTSI(11).mthtsifull = tsiFull;
oTSI(11).mthdatetime = time;

save('oTSI_24_11_03.mat','oTSI')


