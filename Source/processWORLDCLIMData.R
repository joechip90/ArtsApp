processWORLDCLIMData <- function(climNames, climURL, latLimits, longLimits, climFolder = "WORLDCLIMData", climSubFolder = "", outFolder = tempdir()) {
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
	# Sanity check the latitudinal limits
	inLatLimits <- tryCatch(range(as.double(latLimits), na.rm = TRUE), error = function(err) {
		stop(paste("invalid entry for the latitudinal limits:", err, sep = " "))
	})
	if(any(is.na(inLatLimits))) {
		inLatLimits <- c(-Inf, Inf)
	}
	if(inLatLimits[1] == inLatLimits[2]) {
		stop("invalid entry for the latitudinal limits: maximum and minimum values are not different")
	}
	# Sanity check the longitudinal limits
	inLongLimits <- tryCatch(range(as.double(longLimits), na.rm = TRUE), error = function(err) {
		stop(paste("invalid entry for the longitudinal limits:", err, sep = " "))
	})
	if(any(is.na(inLongLimits))) {
		inLongLimits <- c(-Inf, Inf)
	}
	if(inLongLimits[1] == inLongLimits[2]) {
		stop("invalid entry for the longitudinal limits: maximum and minimum values are not different")
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
	# Open the bioclimate data and read it into a list
	bioclimList <- lapply(X = paste(inOutFolder, inClimFolder, inClimSubFolder, climNames, sep = "/"), FUN = function(curClimLoc, latLimits, longLimits) {
		# Retrieve the spatial information from the raster layer
		rasterInfo <- GDALinfo(curClimLoc)
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
		readGDAL(curClimLoc, offset = c(rasterInfo["rows"] - upperYCell, lowerXCell),
			region.dim = c(upperYCell - lowerYCell, upperXCell - lowerXCell))
	}, latLimits = inLatLimits, longLimits = inLongLimits)
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
	# Delete the downloaded and unpacked data
	if(unlink(paste(inOutFolder, "/", inClimFolder, ".zip", sep = "")) != 0) {
		stop("unable to delete downloaded climate data")
	}
	if(unlink(paste(inOutFolder, "/", inClimFolder, sep = ""), recursive = TRUE) != 0) {
		stop("unable to delete downloaded climate data")
	}
	# Return the bioclimate data
	bioclimData
}