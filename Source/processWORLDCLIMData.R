## 1. ------ FUNCTION TO DOWNLOAD AND PROCESS WORLDCLIM DATA ------
#' @title Download and Process Data from WORLDCLIM
#' 
#' @description 
#' This function connects to the server holding the WORLDCLIM data, downloads it, and then
#' processes it into a single \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}}
#' object.  The temporary download files are then deleted.
#' 
#' @param climNames A character vector containing the names of climate variables that are required.
#' @param climURL A character scalar containing the URL where the bioclimate data are kept.
#' @param climFolder A character scalar denoting the sub-folder that the WORLDCLIM data will be placed
#' locally (possibly temporaryily).
#' @param climSubFolder A character scalar denoting the sub-folder where the WORLDCLIM data is stored
#' in the downloaded WORLDCLIM data compressed file.
#' @param climExt A character scalar denoting the extension for the climate data.
#' @param createGeoTIFF A logical scalar.  If \code{TRUE} then a series of GeoTIFF files are created
#' in the \code{outFolder} directory for each layer.
#' @param createRObFile A logical scalar.  If \code{TRUE} then an R object file (with extension .rds) is
#' created in the \code{outFolder} durectory for each layer.
#' @param outFolder A character scalar denoting the folder to place the uncompressed and processed
#' WORLDCLIM data.
#' @param clipX A numeric vector.  \code{range(clipX, na.rm =TRUE)} returns the x-coordinates to clip
#' the imported data too (retains any cells whose centres fall within this range).  Not used if
#' \code{clipPoly} is not \code{NULL}.
#' @param clipY A numeric vector.  \code{range(clipY, na.rm = TRUE)} returns the y-coordinates to clip
#' the imported data too (retains any cells whose centres fall within this range).  Not used if
#' \code{clipPoly} is not \code{NULL}.
#' @param clipCRS PROJ4 string defining the coordinate-reference system of the \code{clipX} and
#' \code{clipY} parameters.  Not used if \code{clipPoly} is not \code{NULL}.
#' @param clipPoly A \code{\link[sp:SpatialPolygons]{SpatialPolygons}} or
#' \code{\link[sp:SpatialPolygonsDataFrame]{SpatialPolygonsDataFrame}} used to clip the selected
#' climate layers.
#' @param progressUpdater A function created using \code{\link[createProgressUpdater]{createProgressUpdater}}
#' that can be used to update a progress inficator.  \code{NULL} means that no progress indication
#' will be provided.
#' 
#' @return An object of type \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing
#' the imported and clipped raster data.
#' 
#' @examples
#' library(ggplot2)																											# Import plotting functions
#' library(rgdal)																												# Import the spatial input/output libraries
#' # The URL of the climate layers
#' climateZipURL <- "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_10m_bio.zip"
#' # The sub-folder path to the climate data in the unpacked zip folder
#' climateZipSubLoc <- ""
#' # The extension of the climate layers
#' climateLayerExt <- "txt"
#' # The folder to put the processed outputs in (and to hold temporary files during
#' # processing)
#' outputLoc <- "temp"
#' # An optional sub-folder to store the processed outputs in (and to hold temporary
#' # files during processing)
#' outputSubLoc <- "WORLDCLIMData"
#' 
#' # Set the names and description of the climate variables to import
#' climateVariableDescriptions <- list(
#' "wc2.0_bio_10m_01" = expression("Annual Mean Temperature (d" * degree * "C)"),
#' "wc2.0_bio_10m_04" = expression("Annual Temperature Range (d" * degree * "C)"),
#' "wc2.0_bio_10m_12" = "Annual Precipitation (mm)",
#' "wc2.0_bio_10m_15" = "Coefficient of Precipitation Variation (%)"
#' )
#' 
#' # Set a progress indicator
#' curProgress <- txtProgressBar(title = "Downloading and processing WORLDCLIM data", style = 3)
#' # Start the processing
#' processedOutput <- processWORLDCLIMData(names(climateVariableDescriptions), climateZipURL, outputSubLoc, climateZipSubLoc,
#'     climateLayerExt, TRUE, TRUE, outputLoc, progressUpdater = createProgressUpdater(curProgress),
#'     clipX = c(3.223898, 31.261007), clipY = c(57.215897, 71.737400), clipCRS = CRS("+init=epsg:3857"))
#' # Close the progress bar
#' close(curProgress)
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[sp:SpatialPolygonsDataFrame]{SpatialPolygonsDataFrame}}
#' @export
#'
processWORLDCLIMData <- function(climNames, climURL, climFolder = "WORLDCLIMData", climSubFolder = "", climExt = "", createGeoTIFF = TRUE,
	createRObFile = TRUE, outFolder = tempdir(), clipX = c(-Inf, Inf), clipY = c(-Inf, Inf), clipCRS = CRS("+init=epsg:4326"), clipPoly = NULL,
	progressUpdater = NULL) {
	### 1.1 ==== Sanity check the function inputs ====
	updateProgressBar(progressUpdater, value = 0.0, details = "processing inputs")
	# Sanity check the climate file extension
	inClimExt <- tryCatch(as.character(climExt), error = function(err) {
		stop(paste("invalid entry for the climate file extension:", err, sep = " "))
	})
	if(length(inClimExt) <= 0) {
		stop("invalid entry for the climate file extension: zero vector length")
	} else if(length(inClimExt) > 1) {
		warning("climate file extension vector length greater than one: only the first element will be used")
		inClimExt <- inClimExt[1]
	}
	if(any(is.na(inClimExt))) {
		stop("invalid entry for the climate file extensions: NA values present")
	}
	# Sanity check the output generation flags
	inCreateGeoTIFF <- tryCatch(as.logical(createGeoTIFF), error = function(err) {
		stop(paste("invalid entry for the GeoTIFF creation flag:", err, sep = " "))
	})
	if(length(inCreateGeoTIFF) <= 0) {
		stop("invalid entry for the GeoTIFF creation flag: zero vector length")
	} else if(length(inCreateGeoTIFF) > 1) {
		warning("GeoTIFF creation flag vector length greater than one: only the first element will be used")
		inCreateGeoTIFF <- inCreateGeoTIFF[1]
	}
	if(any(is.na(inCreateGeoTIFF))) {
		stop("invalid entry for the GeoTIFF creation flag: NA values present")
	}
	inCreateRObFile <- tryCatch(as.logical(createRObFile), error = function(err) {
		stop(paste("invalid entry for the R object file creation flag:", err, sep = " "))
	})
	if(length(inCreateRObFile) <= 0) {
		stop("invalid entry for the R object file creation flag: zero vector length")
	} else if(length(inCreateRObFile) > 1) {
		warning("R object creation flag vector length greater than one: only the first element will be used")
		inCreateRObFile <- inCreateRObFile[1]
	}
	if(any(is.na(inCreateRObFile))) {
		stop("invalid entry for the R object creation flag: NA values present")
	}
	# Sanity check the input names for the climate variables
	inClimNames <- tryCatch(as.character(climNames), error = function(err) {
		stop(paste("invalid entry for the names of the climate variables:", err, sep = " "))
	})
	if(length(inClimNames) <= 0) {
		stop("invalid entry for the names of the climate variables: zero vector length")
	}
	if(any(is.na(inClimNames))) {
		stop("invalid entry for the names of the climate variables: NA values present")
	}
	# Sanity check the input location for the climate URL
	inClimURL <- tryCatch(as.character(climURL), error = function(err) {
		stop(paste("invalid entry for the location of the climate data:", err, sep = " "))
	})
	if(length(inClimURL) <= 0) {
		stop("invalid entry for the location of the climate data: zero vector length")
	} else if(length(inClimURL) > 1) {
		warning("length of vector specifying the location of the climate data is greater than one: only the first element will be used")
		inClimURL <- inClimURL[1]
	}
	if(is.na(inClimURL)) {
		stop("invalid entry for the location of the climate data: value is NA")
	}
	# Sanity check the input of the clip polygon and convert it accordingly
	inClipPoly <- clipPoly
	if(is.null(inClipPoly)) {
		if(!is.null(clipX) && !is.null(clipY) && !is.null(clipCRS)) {
			# Set the clip extent from the inputs
			inClipRange <- tryCatch(cbind(range(clipX, na.rm = TRUE), range(clipY, na.rm = TRUE)), error = function(err) {
				stop(paste("invalid values for the clipping range provided:", err, sep = " "))
			})
			# Respecify the clip coordinates as SpatialPoints
			inClipPoly <- tryCatch(SpatialPoints(inClipRange, proj4string = clipCRS), error = function(err) {
				stop(paste("error thrown whilst converting clipping region to spatial points:", err, sep = " "))
			})
		} else {
			stop("if no clipping polygon is provided then the clipping coordinates must be non-NULL")
		}
	}
	if(is(inClipPoly, "SpatialPoints") || is(inClipPoly, "SpatialPointsDataFrame")) {
		# Input is a spatial points object so use the limits of the coordinates as the spatial extent
		clipPolyCoords <- coordinates(inClipPoly)
		clipPolyCoordsX <- range(clipPolyCoords[, 1])
		clipPolyCoordsY <- range(clipPolyCoords[, 2])
		# Convert these coordinates into a spatial polygon object
		inClipPoly <- tryCatch(SpatialPolygons(list(Polygons(list(Polygon(cbind(
			clipPolyCoordsX[c(1, 1, 2, 2, 1)], clipPolyCoordsY[c(1, 2, 2, 1, 1)]
		))), ID = "borderFrame")), proj4string = CRS(proj4string(inClipPoly))), error = function(err) {
			stop(paste("error thrown whilst converting clipping region to spatial polygons:", err, sep = " "))
		})
	} else if(is(inClipPoly, "SpatialPolygons") || is(inClipPoly, "SpatialPolygonsDataFrame")) {
		# Input is a SpatialPolygons object so use that directly
		inClipPoly <- as(inClipPoly, "SpatialPolygons")
	} else {
		stop("invalid input for the clipping extent: invalid input data type")
	}
	# Sanity check the input of the location of the climate folder
	inClimFolder <- tryCatch(as.character(climFolder), error = function(err) {
		stop(paste("invalid entry for the folder specification:", err, sep = " "))
	})
	if(length(inClimFolder) <= 0) {
		stop("invalid entry for the folder specification: zero vector length")
	} else if(length(inClimFolder) > 1) {
		warning("length of vector specifying the location to the place the climate data is greater than one: only the first element will be used")
		inClimFolder <- inClimFolder[1]
	}
	if(is.na(inClimFolder)) {
		stop("invalid entry for the folder specification: NA values present")
	}
	# Sanity check the input of the location of the climate subfolders
	inClimSubFolder <- tryCatch(as.character(climSubFolder), error = function(err) {
		stop(paste("invalid entry for the subfolder specification:", err, sep = " "))
	})
	if(length(inClimSubFolder) <= 0) {
		stop("invalid entry for the subfolder specification: zero vector length")
	}
	if(any(is.na(inClimSubFolder))) {
		stop("invalid entry for the subfolder specification: NA values present")
	}
	inClimSubFolder <- inClimSubFolder[0:(length(inClimNames) - 1) %% length(inClimSubFolder) + 1]
	# Sanity check the location of the output folder
	inOutFolder <- tryCatch(as.character(outFolder), error = function(err) {
		stop(paste("invalid entry for the output folder location:", err, sep = " "))
	})
	if(length(inOutFolder) <= 0) {
		stop("invalid entry for the output folder location: zero vector length")
	} else if(length(inOutFolder) > 1) {
		warning("length of vector specifying the location of the ouput folder is greater than one: only the first element will be used")
		inOutFolder <- inOutFolder[1]
	}
	### 1.2 ==== Download the WORLDCLIM data ====
	updateProgressBar(progressUpdater, value = 0.2, details = "downloading data")
	# Download the WORLDCLIM data
	if(file.exists(paste(inOutFolder, "/", inClimFolder, ".zip", sep = ""))) {
		# If the file already exists then delete any previous copies of it
		if(unlink(paste(inOutFolder, "/", inClimFolder, ".zip", sep = "")) != 0) {
			stop("unable to delete existing climate data")
		}
	}
	if(download.file(inClimURL, paste(inOutFolder, "/", inClimFolder, ".zip", sep = "")) != 0) {
		stop("unable to download WORLDCLIM data from server")
	}
	# Unzip the WORLDCLIM data
	if(dir.exists(paste(inOutFolder, "/", inClimFolder, sep = ""))) {
		if(unlink(paste(inOutFolder, "/", inClimFolder, sep = ""), recursive = TRUE) != 0) {
			stop("unable to delete existing climate data")
		}
	}
	unzip(paste(inOutFolder, "/", inClimFolder, ".zip", sep = ""),
		exdir = paste(inOutFolder, "/", inClimFolder, sep = ""))
	### 1.3 ==== Import and process the downloaded data ====
	# Open the bioclimate data and read it into a list
	climLocs <- paste(inOutFolder, inClimFolder, sep = "/")
	climLocs <- ifelse(inClimSubFolder == "", climLocs, paste(climLocs, inClimSubFolder, sep = "/"))
	climLocs <- paste(climLocs, climNames, sep = "/")
	if(inClimExt != "") {
		climLocs <- paste(climLocs, inClimExt, sep = ".")
	}
	bioclimList <- lapply(X = 1:length(climNames), FUN = function(curClimIndex, climLocs, inClipPoly, progressUpdater) {
		curClimLoc <- climLocs[curClimIndex]
		updateProgressBar(progressUpdater, value = 0.2 + 0.6 * (curClimIndex - 1) / length(climLocs), details = paste("importing", curClimLoc, sep = " "))
		# Retrieve the spatial information from the raster layer
		rasterInfo <- tryCatch(GDALinfo(curClimLoc), error = function(err) {
			stop(paste("error thrown whilst extracting GDAL information:", err, sep = " "))
		})
		# Retrieve the projection of the input raster and transform the input polygon if it is of another projection
		rasterProjInfo <- attr(rasterInfo, "projection")
		curClipPoly <- inClipPoly
		if(proj4string(curClipPoly) != rasterProjInfo) {
			curClipPoly <- tryCatch(spTransform(curClipPoly, CRS(rasterProjInfo)), error = function(err) {
				stop(paste("error thrown during coordinate projection:", err, sep = " "))
			})
		}
		longLimits <- range(fortify(curClipPoly)[, 1], na.rm = TRUE)
		latLimits <- range(fortify(curClipPoly)[, 2], na.rm = TRUE)
		# Retrieve the maximum and minimum cells of the target region
		lowerXCell <- max(floor((longLimits[1] - rasterInfo["ll.x"]) / rasterInfo["res.x"]), 0)
		lowerYCell <- max(floor((latLimits[1] - rasterInfo["ll.y"]) / rasterInfo["res.y"]), 0)
		upperXCell <- min(ceiling((longLimits[2] - rasterInfo["ll.x"]) / rasterInfo["res.x"]), rasterInfo["columns"])
		upperYCell <- min(ceiling((latLimits[2] - rasterInfo["ll.y"]) / rasterInfo["res.y"]), rasterInfo["rows"])
		# Ensure that the target region intersects with the climate data
		if(upperXCell < lowerXCell || upperYCell < lowerYCell){
			stop("bioclimate data region does not intersect with region of interest")
		}
		# Read the cropped raster data
		outVal <- tryCatch(readGDAL(curClimLoc, offset = c(rasterInfo["rows"] - upperYCell, lowerXCell),
			region.dim = c(upperYCell - lowerYCell, upperXCell - lowerXCell)), error = function(err) {
				stop(paste("error thrown whist important raster information:", err, sep = " "))
			})
		outVal
	}, climLocs = climLocs, inClipPoly = inClipPoly, progressUpdater = progressUpdater)
	# Test to ensure all the climate layers have the same spatial properties
	if(length(bioclimList) > 1) {
		if(!all(sapply(X = bioclimList[2:length(bioclimList)], FUN = function(curSpace, refSpace) {
			identical(as(curSpace, "SpatialGrid"), refSpace)
		}, refSpace = as(bioclimList[[1]], "SpatialGrid")))) {
			stop("bioclimate data does not all share the same spatial properties")
		}
	}
	# Reduce the bioclimate down to a data frame
	bioclimDataFrame <- as.data.frame(do.call(cbind, lapply(X = bioclimList, FUN = function(curClim) {
		curClim$band1
	})))
	colnames(bioclimDataFrame) <- inClimNames
	# Create a new SpatialGridDataFrame containing the entire dataset
	bioclimData <- SpatialGridDataFrame(getGridTopology(bioclimList[[1]]), bioclimDataFrame, CRS(proj4string(bioclimList[[1]])))
	### 1.4 ==== Remove temporary files ====
	# Delete the downloaded and unpacked data
	if(unlink(paste(inOutFolder, "/", inClimFolder, ".zip", sep = "")) != 0) {
		warning("unable to delete downloaded climate data")
	}
	if(unlink(paste(inOutFolder, "/", inClimFolder, sep = ""), recursive = TRUE) != 0) {
		warning("unable to delete downloaded climate data")
	}
	### 1.5 ==== Clip the data to the required spatial extent ====
	updateProgressBar(progressUpdater, value = 0.9, details = "clipping the data to the required spatial extent")
	# Clip the bioclimate data to the relevant region
	outValues <- spatialGridDataClip(bioclimData, clipPoly = inClipPoly)
	### 1.6 ==== Store the outputs ====
	updateProgressBar(progressUpdater, value = 1.0, details = "saving the outputs")
	if(inCreateGeoTIFF || inCreateRObFile) {
		if(dir.create(paste(inOutFolder, inClimFolder, sep = "/"))) {
			# Save the R object file if it is requested
			if(inCreateRObFile) {
				saveRDS(outValues, file = paste(inOutFolder, "/", inClimFolder, "/clippedClimateData.rds", sep = ""))
			}
			# Save the GeoTIFF if they are requested
			if(inCreateGeoTIFF) {
				lapply(X = colnames(outValues@data), FUN = function(curVar, outValues, outPath) {
					writeGDAL(outValues[curVar], paste(outPath, "/", curVar, ".tif", sep = ""), drivername = "GTiff")
				}, outValues = outValues, outPath = paste(inOutFolder, inClimFolder, sep = "/"))
			}
		} else {
			stop("unable to create folder to store outputs")
		}
	}
	outValues
}