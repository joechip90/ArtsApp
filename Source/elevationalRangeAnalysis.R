## 1. ------ FUNCTION TO PERFORM ELEVATIONAL RANGE ANALYSIS ------
#' @title Produce Elevational Range Analysis
#' 
#' @description 
#' This function performs a generalised linear regression model with a spatial SPDE random effect using INLA on
#' provided occurrence data and a set of linear covariates.  This model then calculate the statistics relating
#' to the elevational range of each of the species included in the analysis.
#' 
#' @param occurrenceGrid A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing the gridded
#' occurrence data.  All values in the associated \code{data.frame} greater than zero will be treated as a one
#' for the purposes of the model.
#' @param bioclimateGrid A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing the climate data.
#' @param elevationGrid A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing the elevation data.
#' @param alpha A \code{numeric} scalar containing the SPDE fractional operater order.
#' @param meshParameters A \code{list} containing the parameters to use in the construction of the spatial mesh
#' to approximate the random effects.  These parameters are passed to the \code{\link[INLA::inla.mesh.2d]{inla.mesh.2d}}
#' and \code{\link[INLA::inla.nonconvex.hull]{inla.nonconvex.hull}} functions.
#' @param meshBoundary A \code{\link[sp::SpatialPolygons]{SpatialPolygons}} object containing the boundary of the
#' mesh.  If this is \code{NULL} then instead construct a border polygon using the \code{\link[INLA::inla.nonconvex.hull]{inla.nonconvex.hull}}
#' function and parameters passed to \code{meshParameters}.
#' @param progressUpdater A function created using \code{\link[createProgressUpdater]{createProgressUpdater}}
#' that can be used to update a progress inficator.  \code{NULL} means that no progress indication
#' will be provided.
#' @param responseDensity An \code{integer} scalar.  The number of sampling locations to use when building the
#' response curves.
#' @param outFolder A character scalar denoting the folder to place the SDM output.
#' @param createGeoTIFF A logical scalar.  If \code{TRUE} then a series of GeoTIFF files are created
#' in the \code{outFolder} directory with the model predictions for each species.
#' @param createRObFile A logical scalar.  If \code{TRUE} then an R object file (with extension .rds) is
#' created in the \code{outFolder} with the fitted model object for each species.
#' @param inlaVerbose A logical scalar.  If \code{TRUE}, \code{\link[INLA::inla]{inla}} is run in verbose mode.
#' @param inlaKeep A logical scalar.  If \code{TRUE}, \code{\link[INLA::inla]{inla}} retains working files.  These are
#' stored in \code{outFolder}.
#' @param inlaDebug A logical scalar.  If \code{TRUE}, \code{\link[INLA::inla]{inla}} produces debugging information.
#' @param elevMC An integer scalar.  Determines the number of Monte-Carlo iterations to perform when calculating the
#' distribution of elevational trait statistics.
#'
#' @return A \code{list} containing two elements:
#' \itemize{
#'  \item{modelSummaries}{A list of \code{\link[INLA::inla]{inla}} model objects for each species in the \code{occurrenceData}
#'  \code{data.frame}.}
#'  \item{spatialMesh}{A \code{\link[sp::SpatialPolygons]{SpatialPolygons}} object containing the spatial mesh used in the
#'  \code{\link[INLA::inla]{inla}} function.}
#' }
#' 
#' @examples
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nina.no}
#' @seealso \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} \code{\link[INLA::inla.mesh.2d]{inla.mesh.2d}}
#' \code{\link[INLA::inla.nonconvex.hull]{inla.nonconvex.hull}}
#' @export
#'
runElevationalRangeAnalysis <- function(occurrenceGrid, bioclimateGrid, elevationGrid, alpha = 1.5, meshParameters = list(), meshBoundary = NULL, progressUpdater = NULL, responseDensity = 100,
	outFolder = tempdir(), createGeoTIFF = TRUE, createRObFile = TRUE, inlaVerbose = FALSE, inlaKeep = FALSE, inlaDebug = FALSE, elevMC = 10000) {
	### 1.1 ==== Sanity check the function inputs ====
	# Import the occurrence data SpatialPoints
	inOccurrenceData <- tryCatch(as(occurrenceGrid, "SpatialGridDataFrame"), error = function(err) {
		stop("error encountered importing occurrence data: ", err)
	})
	# Import the bioclimate as SpatialGridDataFrame
	inBioclimateGrid <- tryCatch(as(bioclimateGrid, "SpatialGridDataFrame"), error = function(err) {
		stop("error encountered importing bioclimatic variables: ", err)
	})
	# Import the elevation as SpatialGridDataFrame
	inElevationGrid <- tryCatch(as(elevationGrid, "SpatialGridDataFrame"), error = function(err) {
		stop("error encountered importing elevation data: ", err)
	})
	# Retrieve the elevational values at the relevant points on the bioclimate grid
	queryPoints <- spTransform(SpatialPoints(coordinates(inBioclimateGrid), CRS(proj4string(inBioclimateGrid))), CRS(proj4string(inElevationGrid)))
	elevationValues <- inElevationGrid@data[getGridIndex(coordinates(queryPoints), getGridTopology(inElevationGrid)), 1]
	### 1.2 === Run the species distribution models ====
	# Run the species distribution models (using the spatial version with mesh parameters and ridge regularisation)
	sdmOutputs <- runINLASpatialGLM(inOccurrenceData, inBioclimateGrid, alpha, meshParameters, meshBoundary, progressUpdater, responseDensity, outFolder, createGeoTIFF, createRObFile, inlaVerbose, inlaKeep, inlaDebug)
	### 1.3 === Retrieve the elevational limits for each species ====
	createElevationTraits <- function(curSpecies, elevMC, elevationValues, predictType, occurrenceData) {
		# Retrieve the species predictions
		curSpeciesOutput <- readRDS(paste(outFolder, "Species_", curSpecies, "_ModelPredictions.rds", sep = ""))
		# Calculate the probability of occurrences at each cell by applying the inverse logit to the prediction
		probVals <- 1.0 / (1.0 + exp(-curSpeciesOutput$spatialPredictions@data[, predictType]))
		# Sample a random set of occurrence using the prediction
		randomOcc <- t(sapply(X = probVals, FUN = function(curProbVal, elevMC) {
			runif(elevMC, 0.0, 1.0) <= curProbVal
		}, elevMC = elevMC))
		# Ensure that cells that have an observation are set to TRUE
		randomOcc[occurrenceData@data[, curSpecies], ] <- TRUE
		# Create the lowest, highest, and range of elevations in the sampled grid
		elevRanges <- t(apply(X = randomOcc, FUN = function(curOccVec, elevationValues) {
			# Find the elevations that the species has been observed at
			obsElevation <- elevationValues[isTRUE(curOccVec)]
			# Find the range of elevations
			obsRange <- range(obsElevation, na.rm = TRUE)
			c(obsRange[1], obsRange[2], diff(obsRange))
		}, MARGIN = 1, elevationValues = elevationValues))
		setNames(c(
			mean(elevRanges[, 1]),
			sd(elevRange[, 1]),
			quantile(elevRanges[, 1], probs = c(0.05, 0.5)),
			mean(elevRanges[, 2]),
			sd(elevRange[, 2]),
			quantile(elevRanges[, 2], probs = c(0.95, 0.5)),
			mean(elevRanges[, 3]),
			sd(elevRange[, 3]),
			quantile(elevRanges[, 3], probs = c(0.025, 0.975, 0.5)),
		), c(
			"estMinElevationMean",
			"estMinElevationSD",
			"estMinElevation5thPercentile", "estMinElevationMedian",
			"estMaxElevationMean",
			"estMaxElevationSD",
			"estMaxElevation95thPercentile", "estMaxElevationMedian",
			"estRangeElevationMean",
			"estRangeElevationSD",
			"estRangeElevation2.5thPercentile", "estRangeElevation97.5thPercentile", "estRangeElevationMedian",
		))
	}
	# Create predictions of the elevational range of each of the species dependent upon the modelled distribution
	realisedElevation <- sapply(X = names(inOccurrenceData), FUN = createElevationTraits, elevMC = elevMC, elevationValues = elevationValues, predictType = "meanLinearPred", occurrenceData = inOccurrenceData)
	fundamentalElevation <- sapply(X = names(inOccurrenceData), FUN = createElevationTraits, elevMC = elevMC, elevationValues = elevationValues, predictType = "meanClimatePred", occurrenceData = inOccurrenceData)
	rownames(realisedElevation) <- names(inOccurrenceData)
	colnames(realisedElevation) <- paste(colnames(realisedElevation), "fullPrediction", sep = "_")
	rownames(fundamentalElevation) <- names(inOccurrenceData)
	colnames(fundamentalElevation) <- paste(colnames(fundamentalElevation), "climateOnly", sep = "_")
	# Observed elevational range
	knownElevation <- t(apply(X = as.matrix(inOccurrenceData@data), FUN = function(curObs, elevationValues) {
		range(elevationValues[isTrue(curObs)], na.rm = TRUE)
	}, MARGIN = 2, elevationValues = elevationValues))
	colnames(knownElevation) <- c("observedMinElevation", "observedMaxElevation")
	rownames(knownElevation) <- names(inOccurrenceData)
	# Link all the elevation tables together
	elevationAttributes <- cbind(knownbElevation, realisedElevation, fundamentalElevation)
	# Produce a table in the output folder
	write.csv2(as.data.frame(elevationAttributes), paste(outFolder, "/ElevationAttributes.csv", sep = ""))
	append(sdmOutputs, list(
		elevationAttributes = elevationAttributes
	))
}