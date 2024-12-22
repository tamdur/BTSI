# BTSI
## A Bayesian hierarchical model for assimilating direct and indirect observations of total solar irradiance (TSI), producing joint estimate of TSI, TSI uncertainty, and observation error.

## Getting started
BTSIv1.nc contains the BTSI output for the satellite era, last updated December 22, 2024.

### Running analysis with BTSI code
Download BTSI to a Matlab directory. Running the MCMC model requires running runchain_22_04_25.m from the main project directory, a formatted array of observations (obsmatrix), prior hyperparameters for the observation model, and an observation mask array for observations that are excised. 

A completed observation mask array and set of prior hyperparameters are located within the /mat_files directory as excludeMask_24_12_10.mat and obspriors_22_06_23.mat, respectively. An array of observers, updated as of November 3 2024, is located in /mat_files as obs_24_12_11.mat. 

To run runchain_22_04_25.m, first open the runchain_22_04_25.m script and edit any desired headers under "Headers to modify", including the name and path of the output, the observation array to be used, and any alternative model specifications one would like to use.

## Manuscript plots
gplots_24_11_22.m and plotsupp_24_11_22.m reproduce the figures and tables cited in Amdur and Huybers 24 and shown in the plots folder. These scripts assume access to specific BTSI model outputs performed by the lead author. To request the mat files needed to run these scripts, please contact the corresponding author in Amdur and Huybers 24.

## Updating the MCMC model with new data

### Create a formatted set of observation text files
- Update getdates.m such that dateS.all runs through updated time series
- Update importmgiiupdated to use revised mgii file
- Update makesilso to use revised sunspot file
- Update makeinstruments to use updated satellite observations
- Last, run makeobsfiles.m

### Create a formatted array of observations (obsmatrix)
- Run makeobsmatrix.m, paying attention to whether it is referring to a directory with PMOD-corrected satellite files or not

## Updating the TSI reconstructions used for comparison to BTSI
- Update readothertsi.m to include all updated TSI reconstructions, then save the file in mat_files directory

### Generate an observation mask array for observations that are excised
- Update getpriors.m to use the updated set of oTSI data produced by readothertsi.m
- Update runchain_22_04_25.m to use the updated obsmatrix
- Run btsi_23_03_24.m with the noExclude option to create a base run for exclusion calculation.


