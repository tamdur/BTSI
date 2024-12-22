%Create netcdf file of completed BTSI file 
%
% Ted Amdur
% December 16, 2024

clearvars

load('ar2_24_11_08_long.mat'); %Load BTSI output
load('obs_24_12_11.mat') %Load source data

%Index BTSI to SORCE/TIM
xAll=xAll+offsets(9);xms=mean(xAll,2);

% Create a new NetCDF file
ncid = netcdf.create('../BTSIv1.nc', 'CLOBBER'); 

% Define dimensions
dimid_time = netcdf.defDim(ncid, 'time', 550);
dimid_realization = netcdf.defDim(ncid, 'realization', 10000);

% Define variables
varid_time = netcdf.defVar(ncid, 'time', 'double', dimid_time);
netcdf.putAtt(ncid, varid_time, 'units', 'Julian Date'); 
netcdf.putAtt(ncid, varid_time, 'long_name', 'Time at monthly resolution, expressed as julian date at the center of the respective month');

varid_realization = netcdf.defVar(ncid, 'realization', 'int', dimid_realization);
netcdf.putAtt(ncid, varid_realization, 'long_name', 'BTSI realization number');

varid_TSI = netcdf.defVar(ncid, 'TSI', 'double', [dimid_time dimid_realization]);
netcdf.putAtt(ncid, varid_TSI, 'units', 'W/m^2');
netcdf.putAtt(ncid, varid_TSI, 'long_name', 'Total solar irradiance, 10,000 realizations of BTSI indexed to SORCE/TIM mean TSI'); 

varid_meanTSI = netcdf.defVar(ncid, 'Mean TSI', 'double', dimid_time);
netcdf.putAtt(ncid, varid_meanTSI, 'units', 'W/m^2');
netcdf.putAtt(ncid, varid_meanTSI, 'long_name', 'Mean of 10,000 BTSI total solar irradiance realizations indexed to SORCE/TIM mean TSI'); 

% Add global attributes (metadata)
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'Monthly TSI calculated from BTSI model: mean and 10,000 realizations');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'source', 'BTSI v1');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'institution', 'Aon Impact Forecasting, Harvard University Department of Earth and Planetary Sciences');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'author', 'Ted Amdur');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'author email', 'tamdur4@gmail.com');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'creation_date', string(datetime));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'summary', 'This dataset contains 10,000 realizations of total irradiance as a function of time using the Bayesian TSI (BTSI) approach to combine satellite and proxy sources. TSI values assume the accuracy of mean TSI as estimated by SORCE/TIM.');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'time_coverage_start', '1978-11-01');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'repository', 'https://github.com/tamdur/BTSI');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'time_resolution', 'monthly');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'time_coverage_end', '2024-08-31');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'date_modified', '2024-11-08');
% ... add other global attributes ...

% End define mode
netcdf.endDef(ncid);

% Write data to the NetCDF file
netcdf.putVar(ncid, varid_time, juliandate(dateM)); % Convert datetime to string
netcdf.putVar(ncid, varid_realization, 1:10000); 
netcdf.putVar(ncid, varid_TSI, xAll);
netcdf.putVar(ncid, varid_meanTSI, xms);

% Close the NetCDF file
netcdf.close(ncid);