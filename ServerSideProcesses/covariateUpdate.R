## 1. ------ SET GLOBAL VARIABLES ------
# Set the location of the application
appLoc <- "E:/Dropbox (Personal)/Work/ArtsApp/ArtApp_GitHub"
# The URL of the climate layers
climateZipURL <- "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_10m_bio.zip"
# The sub-folder path to the climate data in the unpacked zip folder
climateZipSubLoc <- ""
# The extension of the climate layers
climateLayerExt <- "tif"
# The folder to put the processed outputs in (and to hold temporary files during
# processing)
outputLoc <- "E:/Temp"
# An optional sub-folder to store the processed outputs in (and to hold temporary
# files during processing)
outputSubLoc <- "WORLDCLIMData"
# Set the location of the clipping region shapefile
clipRegionLoc <- "E:/Dropbox (Personal)/Work/BiodiversitetVeinett/Data/GADMNorway"
# Set the name of the clipping region layer
clipRegionLayerName <- "NOR_adm0"
# Set the names and description of the climate variables to import
climateVariableDescriptions <- list(
	"wc2.0_bio_10m_01" = expression("Annual Mean Temperature (d" * degree * "C)"),
	"wc2.0_bio_10m_04" = expression("Annual Temperature Range (d" * degree * "C)"),
	"wc2.0_bio_10m_12" = "Annual Precipitation (mm)",
	"wc2.0_bio_10m_15" = "Coefficient of Precipitation Variation (%)"
)
## 2. ------ IMPORT LIBRARIES AND SOURCE CODE ------
library(ggplot2)																											# Import plotting functions
library(rgdal)																												# Import the spatial input/output libraries
source(paste(appLoc, "/Source/progressIndicators.R", sep = ""))				# Import the progress indicators
source(paste(appLoc, "/Source/spatialGridDataClip.R", sep = ""))			# Import the spatial clipping functions
source(paste(appLoc, "/Source/processWORLDCLIMData.R", sep = ""))			# Import the WORLDCLIM processing functions

## 3. ------ DOWNLOAD AND PROCESS WORLDCLIM DATA ------
# Retrieve the clipping region
clipRegion <- readOGR(clipRegionLoc, clipRegionLayerName)
# Set a progress indicator
curProgress <- txtProgressBar(title = "Downloading and processing WORLDCLIM data", style = 3)
# Start the processing of the climate data
processedClimateData <- tryCatch(processWORLDCLIMData(names(climateVariableDescriptions), climateZipURL, outputSubLoc, climateZipSubLoc,
	climateLayerExt, TRUE, FALSE, outputLoc, clipPoly = clipRegion, progressUpdater = createProgressUpdater(curProgress)), error = function(err) {
		# If an error occurs half-way through processing then delete the temporary directory
		if(dir.exists(paste(outputLoc, outputSubLoc, sep = "/"))) {
			unlink(paste(outputLoc, outputSubLoc, sep = "/"), recursive = TRUE)
		}
		stop(paste("error encountered during processing of WORLDCLIM data:", err, sep = " "))
})
# Save the outputs
saveRDS(list(climateData = processedClimateData, climateDescription = climateVariableDescriptions), file = paste(outputLoc, "/", outputSubLoc, "/ProcessedClimateData.rds", sep = ""))
# Close the progress bar
close(curProgress)