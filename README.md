# ArtsApp
Repository for storing source code for the automatic processing of occurrence data, environmental covariates, and automate the distribution modelling of various species.  Consists of two main folders:

1. __Source__ General source code containing the functions called by the main processes.
2. __ServerSideProcesses__ Scripts that perform the processing of the data and the distribution modelling.

## Server-side Processes
These scripts are designed to be ran directly on the server.  They all have the a section at the top of the file under the "Set Global Variables" heading where the input and output location are specified and parameters controlling their behaviour.

### covariateUpdate.R
This script downloads the climate data from WORLDCLIM, processes it (clipping out the region of interest), and creates the relevant GeoTIFF and RDS files.  It has the following options at the top of the script, which you can change to suit the way you want the script to behave:

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
This script downloads occurrence data from GBIF and processes them: keeps only records in the clipping region and removes those records with poor spatial information or outside the dates of interest.  It has the following options at the top of the script, which you can change to suit the way you want the script to behave:

1. __appLoc__ The location where the Source and ServerSideProcesses directories reside on the local server.
2. __speciesInfoLoc__ The location of a text file file containing the information of the species you wish to collect data for.  The file must be tab-delimited text file with the first row giving column headers.  There must be at least 2 columns in the text file: 'taxonKey' which gives the unique taxon key used in GBIF for the species, and 'shortName' which gives a short name used to name files related to the species (the will be fewer problems if there is no whitespace in the short name).
3. __templateGridRDSLoc__ The location of the rds file containing the climate data.  This is created in the covariateUpdate.R script and is placed in the output folder specified in that script.  The climate data is used as spatial template for choosing which occurrence points to retain and build outputs for.
4. __outputLoc__ The location to store the outputs.
5. __outputSubLoc__ The subfolder in the output location to store the outputs.  Warning: all files in this directory will be overwritten.
6. __gbifLoginCredentials__ A list containing the GBIF login credentials.  The list should have three components:
    1. __username__ The GBIF username.
    2. __password__ The GBIF password.
    3. __email__ The email used in GBIF registration (it used often used as a notification email that the occurrence data is ready for download).
7. __gbifMaxTime__ The maximum amount of time (in seconds) to wait for GBIF to prepare the download.
8. __gbifAttemptTime__ The amount of time (in seconds) to wait between queries to the GBIF server to see if the download is ready.
9. __yearLimits__ A vector of years specifying the range of years to take observation from (this is given by range(yearLimits, na.rm = TRUE).  For example: 'yearLimits <- c(1960, Inf)' would select all occurrence records between 1960 and the current day.

### runDistributionModels.R
This script runs the species distribution model for the different species and then produces a series of outputs for each modelled species.  It has the following options at the top of the script, which you can change to suit the way you want the script to behave:

1. __appLoc__ The location where the Source and ServerSideProcesses directories reside on the local server.
2. __occupancyLoc__ The location of the rds file containing the processed occurrence records.  This file is produced by the occurrenceDataUpdate.R script.
3. __climateLoc__ The location of the rds file containing the climate data.  This is created in the covariateUpdate.R script and is placed in the output folder specified in that script.
4. __meshBoundaryLoc__ The location of the shapefile containing the boundary information of the region of interest.
5. __meshBoundaryName__ The name of the shapefile layer containing the boundary information of the region of interest.
6. __outputLoc__ The location to store the outputs.
7. __outputSubLoc__ The subfolder in the output location to store the outputs.  Warning: all files in this directory will be overwritten.
8. __spdeAlpha__ The fractional operator order of the underlying SPDE model controlling the spatial random effects.  See [this tutorial](https://folk.ntnu.no/fuglstad/Lund2016/Session6/inla-spde-howto.pdf) for more information on SPDE modelling.
9. __meshSpecification__ A list containing the parameters that set the properties of the spatial mesh over which to define the spatial random effects.  See [this tutorial](https://www.math.ntnu.no/inla/r-inla.org/tutorials/spde/spde-tutorial.pdf) for more information on defining appropriate meshes for SPDE models.
10. __responseDensity__ The density of samples on the response curve plots.  Lower values represent quicker computation times but coarser resolution on the resultant figures.

The script produces a number of files.  For each species there is the following set of files:
1. __[Species]\_MeanEst.tif__ A GeoTIFF.  The mean estimate of the probability of occurrence.
2. __[Species]\_LowerEst.tif__ A GeoTIFF.  The lower estimate (2.5th percentile) of the probability of occurrence.
3. __[Species]\_UpperEst.tif__ A GeoTIFF.  The upper estimate (97.5th percentile) of the probability of occurrence.
4. __[Species]\_UncertaintyEst.tif__ A GeoTIFF.  The width of the prediction interval (UpperEst - LowerEst)
5. __[Species]\_MeanLinearPred.tif__ A GeoTIFF.  The mean estimate of the linear predictor (logit of the probability of occurrence).
6. __[Species]\_MeanClimatePred.tif__ A GeoTIFF.  The mean 'climate-only' part of the prediction with no spatial random effects.
7. __[Species]\_MeanSpatialPred.tif__ A GeoTIFF.  The mean spatial random effect (i.e. the non-climatic effects) on the species distribution.
8. __[Species]\_ModelPredictions.rds__ An RDS file containing the INLA model object and some extra information that can be used for the construction of the response curves.
9. __[Species]\_[Climate].pdf__ A series of PDFs containing the response curves for each the climate variables.

In addition the following files are produced:
1. __spatialEffectsMesh__ ESRI shapefile files containing the spatial random effects Delaunay triangulation.
2. __modelOutputSummary.rds__ An RDS file containing the model summary information and the graphical objects used in plotting the species response curves.