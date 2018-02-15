## 1. ------ FUNCTION TO CLIP SPATIAL GRID DATA FRAMES ------
#' @title Import Spatial Grid Data and Clip to Given Extent
#' 
#' @description
#' This function imports a series of \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}}
#' objects or loads the spatially-referenced raster data from \code{rgdal}-compliant files, clips the
#' data to the required spatial extent, and stores the resulting data in a single
#' \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}}.
#' 
#' @param spatGrid Either a \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} object, list
#' of \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} objects, or a character vector
#' containing the locations of \code{rgdal}-compliant files (they will be imported using the
#' \code{\link[rgdal::readGDAL]{readGDAL}} function).
#' @param clipX A numeric vector.  \code{range(clipX, na.rm =TRUE)} returns the x-coordinates to clip
#' the imported data too (retains any cells whose centres fall within this range).  Not used if
#' \code{clipPoly} is not \code{NULL}.
#' @param clipY A numeric vector.  \code{range(clipY, na.rm = TRUE)} returns the y-coordinates to clip
#' the imported data too (retains any cells whose centres fall within this range).  Not used if
#' \code{clipPoly} is not \code{NULL}.
#' @param clipCRS PROJ4 string defining the coordinate-reference system of the \code{clipX} and
#' \code{clipY} parameters.  Not used if \code{clipPoly} is not \code{NULL}.
#' @param clipPoly A \code{\link[sp:SpatialPolygons]{SpatialPolygons}} or
#' \code{\link[sp:SpatialPolygonsDataFrame]{SpatialPolygonsDataFrame}} used to clip the elements
#' of \code{spatGrid}.
#' 
#' @return An object of type \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing
#' the imported and clipped raster data.
#' 
#' @examples
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[INLA::inla]{inla}}
#' @export
#'
spatialGridDataClip <- function(spatGrid, clipX = c(-Inf, Inf), clipY = c(-Inf, Inf), clipCRS = CRS("+init=epsg:4326"), clipPoly = NULL) {
	outObject <- NULL
	### 1.1 ==== Sanity test the inputs ====
	inSpatInfo <- NULL
	if(is.character(spatGrid)) {
		# The input spatial grid is a character vector so interpret those as file names to be passed
		# the rgdal import function
		if(length(spatGrid) > 0) {
			# Get the GDAL information for the first element in the vector
			inGDALinfo <- GDALinfo(spatGrid[1])
			inSpatInfo <- SpatialGrid(grid = GridTopology(
				cellcentre.offset = setNames(c(
					inGDALinfo[["ll.x"]] + 0.5 * inGDALinfo[["res.x"]],
					inGDALinfo[["ll.y"]] + 0.5 * inGDALinfo[["res.y"]]), c("x", "y")),
				cellsize = setNames(c(inGDALinfo[["res.x"]], inGDALinfo[["res.y"]]), c("x", "y")),
				cells.dim = setNames(c(inGDALinfo[["columns"]], inGDALinfo[["rows"]]), c("x", "y"))
			), proj4string = CRS(attr(inGDALinfo, "projection")))
		} else {
			stop("invalid input for the spatial grid: character vector length zero")
		}
	} else if(is(spatGrid, "SpatialGridDataFrame")) {
		# The input spatial grid is sp object type so retrieve the spatial information from that
		inSpatInfo <- as(spatGrid, "SpatialGrid")
	} else if(is.list(spatGrid)) {
		# The input spatial grid is a list then retrieve the spatial information from the first element
		if(length(spatGrid) > 0 && is(spatGrid[[1]], "SpatialGridDataFrame")) {
			inSpatInfo <- as(spatGrid[[1]], "SpatialGrid")
		} else {
			stop("invalid input for the spatial grid: element of list is not a SpatialGrid")
		}
	} else {
		stop("invalid input for the spatial grid: invalid input data type")
	}
	xlim <- summary(inSpatInfo)$bbox[1, ]
	ylim <- summary(inSpatInfo)$bbox[2, ]
	inClipPoly <- clipPoly
	if(is.null(inClipPoly)) {
		if(!is.null(clipX) && !is.null(clipY) && !is.null(clipCRS)) {
			# Set the clip extent from the inputs
			inClipRange <- tryCatch(cbind(range(climX, na.rm = TRUE), range(climY, na.rm = TRUE)), error = function(err) {
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
	### 1.2 ==== Process the clipping region ====
	# Transform the coordinate system of the clipping polygon to that of the spatial information
	inClipPoly <- spTransform(inClipPoly, CRS = CRS(proj4string(inSpatInfo)))
	# Contract the limits of the spatial extent
	clipPolyCoordsX <- inClipPoly@bbox[1, ]
	clipPolyCoordsY <- inClipPoly@bbox[2, ]
	xlim <- c(max(xlim[1], clipPolyCoordsX[1]), min(xlim[2], clipPolyCoordsX[2]))
	ylim <- c(max(ylim[1], clipPolyCoordsY[1]), min(ylim[2], clipPolyCoordsY[2]))
	# Retrieve the maximum and minimum cells to extract data from
	lowerXCell <- max(floor((xlim[1] - attr(inSpatInfo, "bbox")[1, 1]) / attr(attr(inSpatInfo, "grid"), "cellsize")[1]), 0)
	lowerYCell <- max(floor((ylim[1] - attr(inSpatInfo, "bbox")[2, 1]) / attr(attr(inSpatInfo, "grid"), "cellsize")[2]), 0)
	upperXCell <- min(ceiling((xlim[2] - attr(inSpatInfo, "bbox")[1, 1]) / attr(attr(inSpatInfo, "grid"), "cellsize")[1]), attr(attr(inSpatInfo, "grid"), "cells.dim")[1])
	upperYCell <- min(ceiling((ylim[2] - attr(inSpatInfo, "bbox")[2, 1]) / attr(attr(inSpatInfo, "grid"), "cellsize")[2]), attr(attr(inSpatInfo, "grid"), "cells.dim")[2])
	### 1.3 ==== Set the output spatial grid ====
	# Create a new SpatialGrid that the output will be stored within
	outSpatGrid <- SpatialGrid(GridTopology(
		cellcentre.offset = attr(attr(inSpatInfo, "grid"), "cellcentre.offset") +
			(c(lowerXCell, lowerYCell) + 0.5) * attr(attr(inSpatInfo, "grid"), "cellsize"),
		cellsize = attr(attr(inSpatInfo, "grid"), "cellsize"),
		cells.dim = c(upperXCell - lowerXCell, upperYCell - lowerYCell)), proj4string = CRS(proj4string(inSpatInfo)))
	### 1.4 === Import and process the input spatial grid ====
	# Condense the spatial information down into one SpatialGridDataFrame
	inSpatGrid <- spatGrid
	## 1.4.1 Load file if the input spatial grid parameter is a character ----
	if(is.character(inSpatGrid)) {
		# Import the rasters but selecting only the relevant cells that are required
		inSpatGrid <- lapply(X = spatGrid, FUN = function(curFileLoc, inSpatInfo, inOffset, inRegionDim) {
			readGDAL(curFileLoc, offset = inOffset, region.dim = inRegionDim)
		}, inSpatInfo = inSpatInfo,
			inOffset = c(attr(attr(inSpatInfo, "grid"), "cells.dim")[2] - upperYCell, lowerXCell),
			inRegionDim = c(upperYCell - lowerYCell, upperXCell - lowerXCell))
	}
	## 1.4.2 Process inputs if the input spatial grid parameter is a list ----
	if(is.list(inSpatGrid)) {
		# The input is a list of rasters: combine all the data from all of the rasters in
		# the list into one data frame
		firstGrid <- inSpatGrid[[1]]
		inSpatGridData <- do.call(cbind, lapply(X = inSpatGrid, FUN = function(curGrid, firstGrid) {
			inCurGrid <- tryCatch(as(curGrid, "SpatialGridDataFrame"), error = function(err) {
				stop("input spatial grid list element is not a SpatialGridDataFrame object")
			})
			if(!identical(as(firstGrid, "SpatialGrid"), as(inCurGrid, "SpatialGrid"))) {
				stop("input spatial grids do not have the same spatial extent, resolution, or projection")
			}
			inCurGrid@data
		}, firstGrid = firstGrid))
		inSpatGrid <- SpatialGridDataFrame(as(firstGrid, "SpatialGrid"), inSpatGridData)
	}
	## 1.4.3 Ensure that the input spatial grid is a SpatialGridDataFrame ----
	if(!is(inSpatGrid, "SpatialGridDataFrame")) {
		stop("input spatial grid is not a valid type")
	} else {
		### 1.5 ==== Import and process the input spatial grid ====
		# Look up the relevant indeces in the inSpatGrid and copy those across to the output spatial grid
		outSpatGridData <- inSpatGrid@data[apply(X = coordinates(outSpatGrid), FUN = function(curCoords, llcorner, cellsize, ncells) {
			outIndex <- NA
			# Calculate the cell indeces of the output coordinates in the X and Y direction
			cellIndeces <- floor((curCoords - llcorner) / cellsize)
			if(all(cellIndeces >= c(0, 0)) && all(cellIndeces < ncells)) {
				# Covert those indeces to an index on the flattened data structure
				outIndex <- (ncells[2] - cellIndeces[2]) * ncells[1] + cellIndeces[1] + 1
			}
			outIndex
		}, MARGIN = 1, llcorner = attr(inSpatGrid, "bbox")[, 1], cellsize = attr(attr(inSpatGrid, "grid"), "cellsize"), ncells = attr(attr(inSpatGrid, "grid"), "cells.dim")), ]
		if(is.null(dim(outSpatGridData)) || length(dim(outSpatGridData)) <= 1) {
			dim(outSpatGridData) <- c(length(outSpatGridData), 1)
		}
		# Retrieve the rows that fall within the area of interest and set all other rows to NA
		removeRow <- is.na(over(as(outSpatGrid, "SpatialPoints"), inClipPoly))
		outSpatGridData[removeRow, ] <- NA
		colnames(outSpatGridData) <- colnames(inSpatGrid@data)
		# Return the spatial grid data frame
		outObject <- SpatialGridDataFrame(getGridTopology(outSpatGrid), data = as.data.frame(outSpatGridData), proj4string = CRS(proj4string(outSpatGrid)))
	}
	outObject
}
