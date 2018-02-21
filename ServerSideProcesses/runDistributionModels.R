## 1. ------ SET GLOBAL VARIABLES ------
# Set the location of the application
appLoc <- "E:/Dropbox (Personal)/Work/ArtsApp/ArtApp_GitHub"
# Set the location of the occupancy data (rds file)
occupancyLoc <- "E:/Temp/OccurrenceData/ProcessedOccurrenceData.rds"
# Set the location of the climate covariates (rds file)
climateLoc <- "E:/Temp/WORLDCLIMData/ProcessedClimateData.rds"
# Set the location of the mesh boundary shapefile
meshBoundaryLoc <- "E:/Dropbox (Personal)/Work/BiodiversitetVeinett/Data/GADMNorway"
meshBoundaryName <- "NOR_adm0"
# The folder to put the processed outputs in (and to hold temporary files during
# processing)
outputLoc <- "E:/Temp"
# An optional sub-folder to store the processed outputs in (and to hold temporary
# files during processing)
outputSubLoc <- "SDMOutput"
# Set the SPDE model fractional operator order
spdeAlpha <- 1.5
# Set the mesh specification parameters for the spatial random effect
meshSpecification <- list(
	offset = c(2.0, 4.0),
	cutoff = 0.05,            # Set the minimum vertex distance for the mesh points
	max.edge = c(0.2, 0.6)    # Set the mesh boundary buffering options
)
# Set the point density of the response curve plots
responseDensity <- 100
## 2. ------ IMPORT LIBRARIES AND SOURCE CODE ------
library(grid)                                                         # Import the grid graphics library
library(ggplot2)																											# Import the graphical plotting function libraries
library(ggthemes)                                                     # Import the extra plotting themes
library(rgdal)																												# Import the spatial input/output libraries
library(INLA)																													# Import the nested Laplace approximation libraries
source(paste(appLoc, "/Source/progressIndicators.R", sep = ""))				# Import the progress indicators
source(paste(appLoc, "/Source/runINLAModels_SpatialGLM.R", sep = ""))	# Import the distribution modelling functions
source(paste(appLoc, "/Source/inlaMeshToSpatialPolygons.R", sep = ""))# Import the mesh conversion functions
## 3. ------ PERFORM DISTRIBUTION MODELLING ------
# Import the processed occurrence data
processedOccurrenceData <- readRDS(occupancyLoc)
# Import the processed climate data
processedCovariateData <- readRDS(climateLoc)
# Import the mesh boundary information
meshBoundary <- readOGR(meshBoundaryLoc, meshBoundaryName)
# Set a progress indicator
curProgress <- txtProgressBar(title = "Running species distribution models", style = 3)
# Remove any previously stored model outputs
fullOutputLoc <- paste(outputLoc, outputSubLoc, sep = "/")
if(dir.exists(fullOutputLoc)) {
	# Delete the file if it already exists
	unlink(fullOutputLoc, recursive = TRUE)
}
if(dir.create(fullOutputLoc)) {
	# Run the model outputs
	modelOutputs <- tryCatch(runINLASpatialGLM(processedOccurrenceData$griddedOccurrences, processedCovariateData$climateData, spdeAlpha, meshSpecification,
		meshBoundary, createProgressUpdater(curProgress), responseDensity, fullOutputLoc, TRUE, TRUE), error = function(err) {
			if(dir.exists(fullOutputLoc)) {
				# Delete the file if it already exists
				unlink(fullOutputLoc, recursive = TRUE)
			}
			stop(paste("error encountered during species distribution modelling:", err, sep = " "))
	})
	# Save the mesh object to the output folder (for some reason this needs to a SpatialPolygonsDataFrame object so we have to make a
	# dummy data frame to attach to the output mesh object in order to be able to store it)
	writeOGR(SpatialPolygonsDataFrame(modelOutputs$spatialMesh, data.frame(delaunayID = names(modelOutputs$spatialMesh), row.names = names(modelOutputs$spatialMesh))), fullOutputLoc, "spatialEffectsMesh", "ESRI Shapefile")
	# Create figures of the response curves - one for each species and climate covariate
	responsePlots <- apply(X = as.matrix(processedOccurrenceData$speciesInfo), FUN = function(curInfoRow, occurrencePoints, climateData, fullOutputLoc) {
		# Retrieve the point occurrence records that belong to the current species
		curOccurrencePoints <- occurrencePoints[occurrencePoints$taxonKey == as.integer(curInfoRow["taxonKey"]), ]
		# Retrieve a data frame of climate values at observation points
		gridIndeces <- unique(getGridIndex(coordinates(curOccurrencePoints), getGridTopology(climateData$climateData), all.inside = FALSE))
		gridIndeces <- gridIndeces[!is.na(gridIndeces)]
		curClimateValues <- climateData$climateData@data[gridIndeces, ]
		# Retrieve the response curve information
		curResponseCurve <- readRDS(paste(fullOutputLoc, "/Species_", curInfoRow["taxonKey"], "_ModelPredictions.rds", sep = ""))$responsePredictions
		# Create response curve plots
		outputPlots <- lapply(X = colnames(climateData$climateData@data), FUN = function(curClimateName, curClimateValues, curResponseCurve,
			curClimateDescription, allClimateValues, fullOutputLoc, speciesName) {
				# Create a histogram for the entire climate space
				allHist <- hist(allClimateValues[, curClimateName], plot = FALSE)
				curHist <- hist(curClimateValues[, curClimateName], plot = FALSE, breaks = allHist$breaks)
				# Calculate the proportion of occupied available habitat
				propAvailable <- curHist$counts / allHist$counts
				propAvailable <- ifelse(is.finite(propAvailable) & !is.na(propAvailable), propAvailable, 0.0)
				# Calculate the cell widths for the columns
				cellMids <- allHist$mids
				cellWidths <- (allHist$breaks[2:length(allHist$breaks)] - allHist$breaks[1:(length(allHist$breaks) - 1)]) * 0.98
				# Create the ggplot objects to plot the proportional occupation histogram and the climate-only prediction
				occPropGrob <- ggplot(data.frame(cellMids = cellMids, propAvailable = propAvailable, cellWidths = cellWidths),
					aes(cellMids, propAvailable)) + geom_col(width = cellWidths) + theme_tufte() + geom_rangeframe() +
					labs(x = curClimateDescription[[curClimateName]], y = "Proportion of occupied cells")
				climatePredGrob <- ggplot(curResponseCurve[[curClimateName]], aes(covarVal, meanEst)) + geom_rangeframe() +
					geom_ribbon(aes(ymin = lowerEst, ymax = upperEst), fill = "grey70") + geom_line(size = 1) + theme_tufte() +
					theme(axis.title.x = element_blank()) + labs(y = "Probability of occurrence")
				# Draw the grid objects to an output file
				pdf(paste(fullOutputLoc, "/", speciesName, "_", curClimateName, ".pdf", sep = ""), width = 5, height = 8)
				grid.draw(rbind(ggplotGrob(climatePredGrob), ggplotGrob(occPropGrob), size = "last"))
				dev.off()
				# Set the output graphical object outputs
				list(
					proportionOccupied = occPropGrob,
					probabilityOccurrence = climatePredGrob)
		}, curClimateValues = curClimateValues, curResponseCurve = curResponseCurve, curClimateDescription = climateData$climateDescription,
			allClimateValues = climateData$climateData@data, fullOutputLoc = fullOutputLoc, speciesName = curInfoRow["shortName"])
		# Set the names of the climate response plots
		names(outputPlots) <- colnames(climateData$climateData@data)
		outputPlots
	}, occurrencePoints = processedOccurrenceData$rawOccurrences, climateData = processedCovariateData, fullOutputLoc = fullOutputLoc, MARGIN = 1)
	# Set the names of the species climate response plots
	names(responsePlots) <- processedOccurrenceData$speciesInfo$shortName
	# Save the model summary information
	saveRDS(append(modelOutputs, list(responsePlots = responsePlots)), file = paste(fullOutputLoc, "/modelOutputSummary.rds", sep = ""))
	# Rename the files that use the species IDs with the species short names
	apply(X = as.matrix(processedOccurrenceData$speciesInfo), FUN = function(curRow, fullOutputLoc) {
		# List the files that are the relevant to the current species
		curFiles <- list.files(fullOutputLoc, paste("^Species_", curRow["taxonKey"], "_*", sep = ""))
		# Rename those files
		sapply(X = curFiles, FUN = function(curFile, fullOutputLoc, curName) {
			file.rename(paste(fullOutputLoc, "/", curFile, sep = ""), paste(fullOutputLoc, "/", gsub("^.*_", curName, curFile, perl = TRUE), sep = ""))
		}, fullOutputLoc = fullOutputLoc, curName = curRow["shortName"])
	}, fullOutputLoc = fullOutputLoc, MARGIN = 1)
} else {
	stop("unable to create output directory")
}
# Close the progress bar
close(curProgress)
