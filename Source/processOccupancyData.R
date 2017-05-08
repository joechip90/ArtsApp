processOccupancyData <- function(speciesIDs, templateGrid, outFolder = paste(tempdir(), "occData", sep = "/"),
	yearLimits = c(1960, Inf), GBIFUsername = Sys.getenv("GBIFUser", unset = NA),
	GBIFPasword = Sys.getenv("GBIFPass", unset = NA), GBIFEmail = Sys.getEnv("GBIFEmail", unset = NA),
	GBIFattemptGBIFTime = 1800, maxGBIFTime = 14400) {
	# Sanity check the inputs for the species IDs
	inSpeciesIDs <- tryCatch(as.character(speciesIDs), error = function(err) {
		stop(paste("invalid entry for species IDs:", err, sep = " "))
	})
	if(length(inSpeciesIDs) <= 0) {
		stop("invalid entry for species IDs: zero-length specification vector")
	}
	if(anyNA(inSpeciesIDs)) {
		stop("invalid entry for species IDs: NAs present in specification vector")
	}
	# Sanity check the inputs for templatre grid
	inTemplateGrid <- tryCatch(as(templateGrid, "SpatialGrid"), error = function(err) {
		stop(paste("invalid entry for the template grid:", err, sep = " "))
	})
	# Sanity check the inputs for the output folder
	inOutFolder <- tryCatch(as.character(outFolder), error = function(err) {
		stop(paste("invalid entry for the output folder:", err, sep = " "))
	})
	if(length(inOutFolder) <= 0) {
		stop("invalid entry for the output folder: zero-length specification vector")
	} else if(length(inOutFolder) > 1) {
		warning("output folder specification vector length greater than one: only the first element will be used")
		inOutFolder <- inOutFolder[1]
	}
	if(is.na(inOutFolder)) {
		stop("invalid entry for the output folder: NAs present in specification vector")
	}
	# Sanity check the year limits for the occurrence data
	inYearLimits <- tryCatch(range(as.double(yearLimits), na.rm = TRUE), error = function(err) {
		stop(paste("invalid entry for the year limits:", err, sep = " "))
	})
	if(anyNA(inYearLimits)) {
		stop("invalid entry for the year limits: NAs present in the specification vector")
	}
	# Sanity check the GBIF download attempt interval
	inAttemptGBIFTime <- tryCatch(as.integer(attemptGBIFTime), error = function(err) {
		stop(paste("invalid entry for the GBIF download attempt interval:", err, sep = " "))
	})
	if(length(inAttemptGBIFTime) <= 0) {
		stop("invalid entry for the GBIF download attempt interval: zero-length specification vector")
	} else if(length(inAttemptGBIFTime) > 1) {
		warning("GBIF download attempt interval specification vector length greater than one: only the first element will be used")
		inAttemptGBIFTime <- inAttemptGBIFTime[1]
	}
	if(is.na(inAttemptGBIFTime)) {
		stop("invalid entry for the GBIF download attempt interval: NAs present in specification vector")
	} else if(inAttemptGBIFTime < 1) {
		stop("invalid entry for the GBIF download attempt interval: interval cannot be less than one")
	}
	# Sanity check the GBIF maximum download attempt time
	inMaxGBIFTime <- tryCatch(as.integer(maxGBIFTime), error = function(err) {
		stop(paste("invalid entry for the GBIF maximum download attempt time:", err, sep = " "))
	})
	if(length(inMaxGBIFTime) <= 0) {
		stop("invalid entry for the GBIF maximum download attempt time: zero-length specification vector")
	} else if(length(inMaxGBIFTime) > 1) {
		warning("GBIF maximum download attempt time specification vector length greater than one: only the first element will be used")
		inMaxGBIFTime <- inMaxGBIFTime[1]
	}
	if(is.na(inMaxGBIFTime)) {
		stop("invalid entry for the GBIF maximum download attempt time: NAs present in specification vector")
	} else if(inMaxGBIFTime < 1) {
		stop("invalid entry for the GBIF maximum download attempt time: attempt time cannot be less than one")
	}
	# Sanity check the GBIF username and password
	gbifCredentials <- tryCatch(c(
		as.character(GBIFUsername),
		as.character(GBIFPassword),
		as.character(GBIFEmail)), error = function(err) {
		stop(paste("invalid entry for the GBIF credentials:", err, sep = " "))
	})
	if(length(gbifCredentials) != 3) {
		stop("invalid entry for the GBIF credentials: username, password, and email must have a length of 1")
	}
	if(anyNA(gbifCredentials)) {
		stop("invalid entry for the GBIF credentials: NAs present in the credentials")
	}
	# Retrieve the latitude and longitude limits from the template grid
	latLimits <- as.double(attr(inTemplateGrid, "bbox")[2, ])
	longLimits <- as.double(attr(inTemplateGrid, "bbox")[1, ])
	# Set the extra queries relating to the year limits
	yearQuery <- c()
	if(is.finite(inYearLimits[1])) {
		yearQuery <- c(yearQuery, paste("\t\t\t{\"type\":\"greaterThanOrEquals\", \"key\":\"YEAR\", \"value\":\"", as.character(as.integer(inYearLimits[1])), "\"},", sep = ""))
	}
	if(is.finite(inYearLimits[2])) {
		yearQuery <- c(yearQuery, paste("\t\t\t{\"type\":\"lessThanOrEquals\", \"key\":\"YEAR\", \"value\":\"", as.character(as.integer(inYearLimits[2])), "\"},", sep = ""))
	}
	# Produce a query for the species records from GBIF
	downloadQuery <- paste(c("{",
		paste("\t\"creator\":\"", gbifCredentials[1], "\",", sep = ""),
		paste("\t\"notification_address\":[\"", gbifCredentials[3], "\"],", sep = ""),
		"\t\"predicate\":{",
		"\t\t\"type\":\"and\",",
		"\t\t\"predicates\":[",
		"\t\t\t{\"type\":\"or\", \"predicates\":[",
		paste("\t\t\t\t{\"type\":\"equals\", \"key\":\"TAXON_KEY\", \"value\":\"", inSpeciesIDs, "\"}", sep = "", collapse = ",\n"),
		"\t\t\t]},",
		"\t\t\t{\"type\":\"equals\", \"key\":\"HAS_COORDINATE\", \"value\":\"TRUE\"},",
		"\t\t\t{\"type\":\"equals\", \"key\":\"HAS_GEOSPATIAL_ISSUE\", \"value\":\"FALSE\"},",
		paste("\t\t\t{\"type\":\"greaterThanOrEquals\", \"key\":\"DECIMAL_LONGITUDE\", \"value\":\"", as.character(longLimits[1]), "\"},", sep = ""),
		paste("\t\t\t{\"type\":\"lessThanOrEquals\", \"key\":\"DECIMAL_LONGITUDE\", \"value\":\"", as.character(longLimits[2]), "\"},", sep = ""),
		paste("\t\t\t{\"type\":\"greaterThanOrEquals\", \"key\":\"DECIMAL_LATITUDE\", \"value\":\"", as.character(latLimits[1]), "\"},", sep = ""),
		paste("\t\t\t{\"type\":\"lessThanOrEquals\", \"key\":\"DECIMAL_LATITUDE\", \"value\":\"", as.character(latLimits[2]), "\"}", sep = ""),
		yearQuery,
		"\t\t]",
		"\t}",
	"}\n"), collapse = "\n")
	# Retrieve a download ID from GBIF
	downloadID <- occ_download(body = downloadQuery, user = gbifCredentials[1], pwd = gbifCredentials[2], email = gbifCredentials[3])
	# Wait for GBIF to process the download request
	totalTime <- 0
	while(totalTime < inMaxGBIFTime & (
		occ_download_meta(downloadID)$status == "RUNNING" ||
		occ_download_meta(downloadID)$status == "PREPARING"
	)) {
		Sys.sleep(min(inMaxGBIFTime - totalTime, inAttemptGBIFTime))
		totalTime <- totalTime + inAttemptGBIFTime
	}
	if(occ_download_meta(downloadID)$status != "SUCCESS") {
		stop("unable to acquire occurrence data from GBIF")
	}
	# Delete any existing copy of the data
	if(dir.exists(inOutFolder)) {
		if(unlink(inOutFolder, recursive = TRUE) != 0) {
			stop("unable to delete existing occupancy data")
		}
	}
	# Download and import the occurrence data into R
	occDataRaw <- tryCatch(
		as.data.frame(occ_download_import(occ_download_get(as.character(downloadID), path = inOutFolder), downloadID, path = inOutFolder)),
	error = function(err) {
		stop(paste("error encountered during download of occurrence data:", err, sep = " "))
	})
	# Get the occurrence count information
	occCounts <- t(apply(X = coordinates(templateGrid), FUN = function(curCoords, xRes, yRes, specIDs, occData) {
		sapply(X = specIDs, FUN = function(curSpec, xRes, yRes, curCoords, occData) {
			sum(ifelse(
				as.character(curSpec) == as.character(occData$taxonKey) &
				curCoords[1] - xRes / 2.0 <= occData$xCoord &
				curCoords[1] + xRes / 2.0 > occData$xCoord &
				curCoords[2] - yRes / 2.0 <= occData$yCoord &
				curCoords[2] + yRes / 2.0 > occData$yCoord,
			1, 0))
		}, xRes = xRes, yRes = yRes, curCoords = curCoords,
			occData = occData[!is.na(occData$xCoord) & !is.na(occData$yCoord) & !is.na(occData$taxonKey), ])
	}, MARGIN = 1, xRes = templateGrid@grid@cellsize[1], yRes = templateGrid@grid@cellsize[2], specIDs = inSpeciesIDs,
		occData = data.frame(xCoord = occDataRaw$decimalLongitude, yCoord = occDataRaw$decimalLatitude, taxonKey = occDataRaw$taxonKey)))
	dim(occCounts) <- c(nrow(coordinates(templateGrid)), length(inSpeciesIDs))
	colnames(occCounts) <- paste("numRecords", inSpeciesIDs, sep = "_")
	# Clear up the folder after processing
	if(unlink(inOutFolder, recursive = TRUE) != 0) {
		warning("unable to delete downloaded occupancy data")
	}
	# Return the gridded and raw occurrence data
	list(
		rawOccurrences = occDataRaw,
		griddedOccurrences = SpatialGridDataFrame(getGridTopology(templateGrid), as.data.frame(occCounts), CRS(proj4string(templateGrid)))
	)
}