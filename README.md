# ArtsApp
Repository for storing source code for the automatic processing of occurrence data, environmental covariates, and automate the distribution modelling of various species.  Consists of two main folders:

1. __Source__ General source code containing the functions called by the main processes.
2. __ServerSideProcesses__ Scripts that perform the processing of the data and the distribution modelling.

## Server-side Processes
These scripts are designed to be ran directly on the server.  They all have the a section at the top of the file under the "Set Global Variables" heading where the input and output location are specified and parameters controlling their behaviour.

### covariateUpdate.R
This script downloads the climate data from WORLDCLIM, processes it (clipping out the region of interest), and creates the relevant GeoTIFF and RDS files.  It has the following options at the top of the script, which you can change to suit which way you want the script to behave:

1. __appLoc__ The location where the Source and ServerSideProcesses directories reside on the local server.
2. __climateZipURL__ The URL of the zip file of the climate data.
3. __climateZipSubLoc__ A vector of sub-folders containing the climate data in the zip file.  In WORLDCLIM 2.0 they are all stored in the base of the zip file, so this is no longer required.
4. __climateLayerExt__ A file extension for the climate layers.
5. __outputLoc__ The location to store the outputs.
6. __outputSubLoc__ The subfolder in the output location to store the outputs.  Warning: all files in this directory will be overwritten.
7. __clipRegionLoc__ The directory where the shapefile denoting the clipping region is stored.
8. __clipRegionLayerName__ The layer name of the shapefile denoting the clipping region.
9. __climateVariableDescriptions__ A list containing the file names of the climate variables to use (without the file suffix) as a names attribute and a description of the environmental variable.  The description is what will be used in figures of the environmental variable (response functions and histograms).

### occurrenceDataUpdate.R
This script downloads occurrence data from GBIF and processes them: keeps only records in the clipping region and removes those records with poor spatial information or outside the dates of interest.  It has the following options at the top of the script, which you can change to suit which way you want the script to behave:

1. __appLoc__ The location where the Source and ServerSideProcesses directories reside on the local server.
2. __speciesInfoLoc__ The location of a text file file containing the information of the species you wish to collect data for.  The file must be tab-delimited text file with the first row giving column headers.  There must be at least 2 columns in the text file: 'taxonKey' which gives the unique taxon key used in GBIF for the species, and 'shortName' which gives a short name used to name files related to the species (the will be fewer problems if there is no whitespace in the short name).
3. __templateGridRDSLoc_ The location of the rds file containing the climate data.  This is created in the covariateUpdate.R script and is placed in the output folder specified in that script.  The climate data is used as spatial template for choosing which occurrence points to retain and build outputs for.
4. __outputLoc__ The location to store the outputs.
5. __outputSubLoc__ The subfolder in the output location to store the outputs.  Warning: all files in this directory will be overwritten.
6. __gbifLoginCredentials__ A list containing the GBIF login credentials.  The lsit should have three components:
  6.1. __username__ The GBIF username
  6.2. __password__ The GBIF password
  6.3. __email__ The email used in GBIF registration (it used often used as a notification email that the occurrence data is ready for download)
7. __gbifMaxTime__ The maximum amount of time (in seconds) to wait for GBIF to prepare the download.
8. __gbifAttemptTime__ The amount of time (in seconds) to wait between queries to the GBIF server to see if the download is ready.
9. __yearLimits__ A vector of years specifying the range of years to take observation from (this is given by range(yearLimits, na.rm = TRUE).  For example: 'yearLimits <- c(1960, Inf)' would select all occurrence records between 1960 and the current day.