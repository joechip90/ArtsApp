## 1. ------ SET GLOBAL VARIABLES ------
# This script relies on you setting an environment variable "WORKSPACE_ARTSAPP_ELEVATION" to set a workspace
# location and also setting "WORKSPACE_ARTSAPP_GITHUB" to assign the local position of the ArtsApp source files.
# You can set environmental variables in R using the usethis::edit_r_environ function
# Set the location of the application
appLoc <- Sys.getenv("WORKSPACE_ARTSAPP_GITHUB")
# Set a workspace location for the input data files
workspaceLoc <- Sys.getenv("WORKSPACE_ARTSAPP_ELEVATION")
# Name of the data file
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
# Set a directory for the output to be stored
outputLoc <- "C:/Temp/ElevationRangeAnalysisTest"
# Number of Monte Carlo simulations used to estimate elevational attributes
elevMC <- 10000
## 2. ------ IMPORT LIBRARIES AND SOURCE CODE ------
library(grid)                                                       			# Import the grid graphics library
library(ggplot2)																													# Import the graphical plotting function libraries
library(ggthemes)                                                     		# Import the extra plotting themes
library(rgdal)																														# Import the spatial input/output libraries
library(raster)																														# Import the raster handling routines
library(INLA)																															# Import the nested Laplace approximation libraries
library(elevatr)																													# Import the elevation data import library
source(paste(appLoc, "/Source/progressIndicators.R", sep = ""))						# Import the progress indicators
source(paste(appLoc, "/Source/spatialGridDataClip.R", sep = ""))					# Import the spatial grid manipulation functions
source(paste(appLoc, "/Source/processWORLDCLIMData.R", sep = ""))					# Import the WORLDCLIM data processing functions
source(paste(appLoc, "/Source/runINLAModels_SpatialGLM.R", sep = ""))			# Import the distribution modelling functions
source(paste(appLoc, "/Source/inlaMeshToSpatialPolygons.R", sep = ""))   	# Import the mesh conversion functions
source(paste(appLoc, "/Source/elevationalRangeAnalysis.R", sep = ""))			# Import the elevational range analysis
# The URL of the climate layers
climateZipURL <- "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_5m_bio.zip"
# The sub-folder path to the climate data in the unpacked zip folder
climateZipSubLoc <- ""
# The extension of the climate layers
climateLayerExt <- "tif"
# The folder to put the processed outputs in (and to hold temporary files during
# processing)
temporaryWorkspaceLoc <- tempdir()
# An optional sub-folder to store the processed outputs in (and to hold temporary
# files during processing)
temporaryWORLDCLIMSubLoc <- "WORLDCLIMData"
# Set the names and description of the climate variables to import
climateVariableDescriptions <- list(
	"wc2.0_bio_5m_01" = expression("Annual Mean Temperature (d" * degree * "C)"),
	"wc2.0_bio_5m_04" = expression("Annual Temperature Range (d" * degree * "C)"),
	"wc2.0_bio_5m_12" = "Annual Precipitation (mm)",
	"wc2.0_bio_5m_15" = "Coefficient of Precipitation Variation (%)"
)
# Set the CRS string of the occurrence data
crsOccurrence <- CRS("+init=epsg:4326")
# Create the output folder for the analysis
if(dir.exists(outputLoc)) {
	unlink(outputLoc, recursive = TRUE)
}
dir.create(outputLoc)
## 3. ------ PERFORM ELEVATIONAL RANGE ANALYSIS ------
# Import the occurrence data
rawOccurrenceData <- read.table(paste(workspaceLoc, "ExampleOccurrenceData_gbif_sp_costaRica.csv", sep = "/"), header = TRUE)
occurrenceData <- SpatialPointsDataFrame(as.matrix(rawOccurrenceData[, c("longitude", "latitude")]), data = rawOccurrenceData[, c("species", "date", "basisOfRecord")], proj4string = crsOccurrence)
# Set the limits of the analysis to be the extent of the occurrence data
xExtent <- bbox(occurrenceData)[1, ] + c(-0.05, 0.05) * diff(bbox(occurrenceData)[1, ])
yExtent <- bbox(occurrenceData)[2, ] + c(-0.05, 0.05) * diff(bbox(occurrenceData)[2, ])
# Download and process the bioclimate layers
bioclimateGrid <- processWORLDCLIMData(names(climateVariableDescriptions), climateZipURL, temporaryWORLDCLIMSubLoc,
	climateZipSubLoc, climateLayerExt, TRUE, TRUE, temporaryWorkspaceLoc, progressUpdater = NULL,
	clipX = xExtent, clipY = yExtent, clipCRS = crsOccurrence)
# Download and process the elevation layer
elevationGrid <- as(get_elev_raster(stack(bioclimateGrid), 5), "SpatialGridDataFrame")
# Change the coordinate system of the occurrence data to match the coordinate system of the bioclimate data
occcurrenceData <- spTransform(occurrenceData, CRS(proj4string(bioclimateGrid)))
# Retrieve the list of species in the occurrence data
speciesNames <- unique(as.character(occurrenceData@data$species))
# Reshape the spatial points object to spatial grid object with '1' being a presence in a cell and '0' being an absence in a cell
occurrenceGrid <- SpatialGridDataFrame(
	as(bioclimateGrid, "SpatialGrid"),
	as.data.frame(setNames(lapply(
		X = speciesNames, FUN = function(curSpecies, inGridTopo, occurrenceData) {
			# Initialise an output full of absences
			curOutput <- rep(0, prod(inGridTopo@cells.dim))
			# set all cells equal to 1 that have coordinates of the current species within them
			curOutput[getGridIndex(coordinates(occurrenceData)[as.character(occurrenceData@data$species) == curSpecies, , drop = FALSE], inGridTopo)] <- 1
			curOutput
		}, inGridTopo = getGridTopology(bioclimateGrid), occurrenceData = occurrenceData),
	gsub(" ", "_", speciesNames, fixed = TRUE))),
	proj4string = CRS(proj4string(bioclimateGrid))
)
# Blank out those cells that don't have bioclimate data (i.e. fall in the sea, outside region of interest)
occurrenceGrid@data[apply(X = as.matrix(bioclimateGrid@data), FUN = anyNA, MARGIN = 1), ] <- NA
curProgress <- NULL
# Run the elevational range analysis function (with a progress indicator)
{
	# Set a progress indicator
	#curProgress <- txtProgressBar(title = "Running species elevational range analysis", style = 3)
	# Run the elevational range analysis
	rangeOutput <- runElevationalRangeAnalysis(
		occurrenceGrid[1:3],
		bioclimateGrid,
		elevationGrid,
		spdeAlpha,
		meshSpecification,
		outFolder = outputLoc,
		responseDensity = responseDensity,
		progressUpdater = curProgress,
		elevMC = elevMC)
	# Close the progress indicator
	#close(curProgess)
}
