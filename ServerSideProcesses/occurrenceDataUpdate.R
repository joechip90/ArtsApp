## 1. ------ SET GLOBAL VARIABLES ------
# Set the location of the application
appLoc <- "E:/Dropbox (Personal)/Work/ArtsApp/ArtApp_GitHub"
# Set the location of the species information csv file.  This file must have at least the following two columns:
#		> taxonKey - The taxon ID (as used in GBIF) denoting the taxon of interest
#		> shortName - A description of the taxon to be used in file output names
speciesInfoLoc <- "E:/Dropbox (Personal)/Work/ArtsApp/SpeciesList.txt"
# Set the location of the template grid RDS file
templateGridRDSLoc <- "E:/Temp/WORLDCLIMData/ProcessedClimateData.rds"
# The folder to put the processed outputs in (and to hold temporary files during
# processing)
outputLoc <- "E:/Temp"
# An optional sub-folder to store the processed outputs in (and to hold temporary
# files during processing)
outputSubLoc <- "OccurrenceData"
# GBIF login credentials
gbifLoginCredentials <- list(
	username = "joechip90",               # The GBIF username  
	password = "badt00th",                # The GBIF password
	email = "joechip90@googlemail.com"    # The GBIF email
)
# Maximum time (in seconds) to wait before failure when attempting to download occurrence data
gbifMaxTime <- 14400
# Time (in seconds) to wait between checks to see if GBIF data is ready to download
gbifAttemptTime <- 300
# The range of years to download (inclusive)
yearLimits <- c(1960, Inf)
## 2. ------ IMPORT LIBRARIES AND SOURCE CODE ------
library(rgdal)																												# Import the spatial input/output libraries
library(rgbif)																												# Import the GBIF download libraries
source(paste(appLoc, "/Source/progressIndicators.R", sep = ""))				# Import the progress indicators
source(paste(appLoc, "/Source/processOccupancyData.R", sep = ""))			# Import the occurrence data processing functions
## 3. ------ DOWNLOAD AND PROCESS OCCUPANCY DATA ------
# Import the species information from the species information file
speciesInfo <- read.table(speciesInfoLoc, header = TRUE)
rownames(speciesInfo) <- speciesInfo$taxonKey                         # Set the row names to be the taxon key
# Import the template grid from the climate data
templateGrid <- readRDS(templateGridRDSLoc)$climateData
# Set a progress indicator
curProgress <- txtProgressBar(title = "Downloading and processing GBIF occurrence data", style = 3)
# Process the occurrence data
processedOccurrenceData <- processOccupancyData(speciesInfo$taxonKey, templateGrid, paste(outputLoc, outputSubLoc, sep = "/"),
	yearLimits, gbifLoginCredentials$username, gbifLoginCredentials$password, gbifLoginCredentials$email,
	gbifAttemptTime, gbifMaxTime, createProgressUpdater(curProgress))
# Close the progress bar
close(curProgress)
# Save the outputs of the occurrence data
fullOutputLoc <- paste(outputLoc, outputSubLoc, sep = "/")
if(dir.exists(fullOutputLoc)) {
	# Delete the file if it already exists
	unlink(fullOutputLoc, recursive = TRUE)
}
if(exists("processedOccurrenceData")) {
	# Create the output folder
	dir.create(fullOutputLoc)
	lapply(X = speciesInfo$taxonKey, FUN = function(curID, outputLoc, occurrenceData, shortNames) {
		# Retrieve the current short name for the species
		curName <- shortNames[as.character(curID)]
		# Retrieve the current occurrences for the species
		curOccurrences <- occurrenceData$rawOccurrences[occurrenceData$rawOccurrences@data$taxonKey == curID, ]
		# Write a shapefile containing the relevant occurrences
		writeOGR(curOccurrences, outputLoc, curName, "ESRI Shapefile")
		# Write a GeoTIFF containing the occurrence density
		writeGDAL(occurrenceData$griddedOccurrences[as.character(curID)], paste(outputLoc, "/", curName, ".tif", sep = ""), "GTiff")
	}, outputLoc = fullOutputLoc, occurrenceData = processedOccurrenceData, shortNames = setNames(as.character(speciesInfo$shortName), speciesInfo$taxonKey))
}
# Save the entire processed occurrence data
saveRDS(append(processedOccurrenceData, list(speciesInfo = speciesInfo)), file = paste(outputLoc, "/", outputSubLoc, "/ProcessedOccurrenceData.rds", sep = ""))
