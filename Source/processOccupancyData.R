## 1. ------ FUNCTION TO DOWNLOAD AND PROCESS SPECIES OCCUPANCY DATA ------
#' @title Download and Process Data from GBIF
#' 
#' @description 
#' This function connects to the GBIF server and downloads the relevant
#' occurrence data and processes them into two formats: one
#' \code{\link[sp::SpatialPointsDataFrame]{SpatialPointsDataFrame}}
#' object that contains the raw occurrence data but with the projection
#' information added, and one \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}}
#' containing the observation density for each downloaded species.
#' 
#' @param speciesIDs An integer vector containing the taxanomic ID codes of the species you wish to
#' download from GBIF.
#' @param templateGrid A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} object to use as
#' a spatial template for the downloaded points.  Only records that fall on cells with non-\code{NA}
#' values in the first layer of \code{templateGrid} will be included.
#' @param outFolder Scalar character giving the location to store the temporary outputs (warning: contents
#' will be deleted after processing).
#' @param yearLimits Integer vector where \code{range(yearLimits, na.rm = TRUE)} specifies the lower and
#' upper ranges of the years to use when selecting GBIF data.
#' @param GBIFUsername Scalar character conatining the username of the GBIF user.
#' @param GBIFPassword Scalar character containing the password of the GBIF user.
#' @param GBIFEmail Scalar character containing the email address of the GBIF user.
#' @param GBIFAttemptTime Scalar integer denoting the number of seconds to wait between each
#' connection to GBIF to query the status of the processed download.
#' @param maxGBIFTime Scalar integer denoting the maximum number of seconds to wait before the connection
#' to GBIF fails.
#' @param progressUpdater A function created using \code{\link[createProgressUpdater]{createProgressUpdater}}
#' that can be used to update a progress inficator.  \code{NULL} means that no progress indication
#' will be provided.
#' 
#' @return A list containing two elements:
#' \describe{
#'     \item{\code{rawOccurrences}}{A \code{\link[sp::SpatialPointsDataFrame]{SpatialPointsDataFrame}}
#'     object containing the occurrence data downloaded with GBIF.}
#'     \item{\code{griddedOccurrences}}{A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}}
#'     object containing the gridded versions of the occurrence data with the density of observations
#'     for each species in \code{speciesIDs}.}
#' }
#' 
#' @examples
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}}
#' \code{\link[sp::SpatialPointsDataFrame]{SpatialPointsDataFrame}}
#' @export
#'
processOccupancyData <- function(speciesIDs, templateGrid, outFolder = paste(tempdir(), "occData", sep = "/"),
	yearLimits = c(1960, Inf), GBIFUsername = Sys.getenv("GBIFUser", unset = NA),
	GBIFPassword = Sys.getenv("GBIFPass", unset = NA), GBIFEmail = Sys.getenv("GBIFEmail", unset = NA),
	GBIFAttemptTime = 1800, maxGBIFTime = 14400, progressUpdater = NULL) {
	### 1.1 ==== Sanity check the function inputs ====
	updateProgressBar(progressUpdater, value = 0.0, details = "processing inputs")
	# Set the GBIF projection information
	GBIFcrsString <- "+proj=longlat +ellps=WGS84 +no_defs"
	GBIFcrs <- CRS(GBIFcrsString)
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
	inTemplateGrid <- templateGrid
	if(class(inTemplateGrid) != "SpatialGrid") {
		inTemplateGrid <- tryCatch(as(templateGrid, "SpatialGridDataFrame"), error = function(err) {
			stop(paste("invalid entry for the template grid:", err, sep = " "))
		})
	}
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
	inAttemptGBIFTime <- tryCatch(as.integer(GBIFAttemptTime), error = function(err) {
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
	### 1.2 ==== Create GBIF credentials ====
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
	updateProgressBar(progressUpdater, value = 0.4, details = "contacting GBIF database")
	### 1.3 ==== Create GBIF query ====
	# Retrieve the latitude and longitude limits from the template grid
	latLongCoords <- cbind(
		attr(inTemplateGrid, "bbox")[1, ][c(1, 1, 2, 2)],
		attr(inTemplateGrid, "bbox")[2, ][c(1, 2, 1, 2)])
	rownames(latLongCoords) <- NULL
	latLongLimits <- SpatialPoints(latLongCoords, proj4string = CRS(proj4string(inTemplateGrid)))
	# Project the coordinates if neccessary
	if(proj4string(latLongLimits) != GBIFcrsString) {
		latLongLimits <- spTransform(latLongLimits, GBIFcrs)
	}
	latLimits <- range(coordinates(latLongLimits)[, 2], na.rm = TRUE)
	longLimits <- range(coordinates(latLongLimits)[, 1], na.rm = TRUE)
	# Set the extra queries relating to the year limits
	yearQuery <- c()
	if(is.finite(inYearLimits[1])) {
		yearQuery <- c(yearQuery, paste("\t\t\t{\"type\":\"greaterThanOrEquals\", \"key\":\"YEAR\", \"value\":\"", as.character(as.integer(inYearLimits[1])), "\"}", sep = ""))
	}
	if(is.finite(inYearLimits[2])) {
		yearQuery <- c(yearQuery, paste("\t\t\t{\"type\":\"lessThanOrEquals\", \"key\":\"YEAR\", \"value\":\"", as.character(as.integer(inYearLimits[2])), "\"}", sep = ""))
	}
	if(length(yearQuery) > 1) {
		yearQuery[1:(length(yearQuery) - 1)] <- paste(yearQuery[1:(length(yearQuery) - 1)], ",", sep = "")
	}
	if(length(yearQuery) > 0) {
		yearQuery[1] <- paste(",", yearQuery[1], sep = "")
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
	### 1.4 ==== Download GBIF occurrence data ====
	# Retrieve a download ID from GBIF
	downloadID <- tryCatch(occ_download(body = downloadQuery, user = gbifCredentials[1], pwd = gbifCredentials[2], email = gbifCredentials[3]), error = function(err) {
		stop(paste("error encountered during initialisation of GBIF download:", err, "(check login credentials)", sep = " "))
	})
	# Wait for GBIF to process the download request
	totalTime <- 0
	while(totalTime < inMaxGBIFTime & (
		occ_download_meta(downloadID)$status == "RUNNING" ||
		occ_download_meta(downloadID)$status == "PREPARING"
	)) {
		Sys.sleep(min(inMaxGBIFTime - totalTime, inAttemptGBIFTime))
		totalTime <- totalTime + inAttemptGBIFTime
	}
	updateProgressBar(progressUpdater, value = 0.6, details = "downloading occurrence data")
	if(!(occ_download_meta(downloadID)$status %in% c("SUCCEEDED", "SUCCESS"))) {
		stop("unable to acquire occurrence data from GBIF")
	}
	# Delete any existing copy of the data
	if(dir.exists(inOutFolder)) {
		if(unlink(inOutFolder, recursive = TRUE) != 0) {
			stop("unable to delete existing occupancy data")
		}
	}
	if(!dir.create(inOutFolder)) {
		stop("unable to create output folder")
	}
	# Download and import the occurrence data into R
	occDataRaw <- tryCatch(
		as.data.frame(occ_download_import(occ_download_get(as.character(downloadID), path = inOutFolder), downloadID, path = inOutFolder)),
	error = function(err) {
		stop(paste("error encountered during download of occurrence data:", err, sep = " "))
	})
	updateProgressBar(progressUpdater, value = 0.8, details = "processing occurrence data")
	### 1.5 ==== Process occurrence data ====
	# Get the occurrence count information
	occCounts <- t(apply(X = coordinates(templateGrid), FUN = function(curCoords, xRes, yRes, specIDs, occData) {
		outVals <- sapply(X = specIDs, FUN = function(curSpec, xRes, yRes, curCoords, occData) {
			sum(ifelse(
				as.character(curSpec) == as.character(occData$taxonKey) &
				curCoords[1] - xRes / 2.0 <= occData$xCoord &
				curCoords[1] + xRes / 2.0 > occData$xCoord &
				curCoords[2] - yRes / 2.0 <= occData$yCoord &
				curCoords[2] + yRes / 2.0 > occData$yCoord,
			1, 0))
		}, xRes = xRes, yRes = yRes, curCoords = curCoords,
			occData = occData[!is.na(occData$xCoord) & !is.na(occData$yCoord) & !is.na(occData$taxonKey), ])
	}, MARGIN = 1, xRes = inTemplateGrid@grid@cellsize[1], yRes = inTemplateGrid@grid@cellsize[2], specIDs = inSpeciesIDs,
		occData = data.frame(xCoord = occDataRaw$decimalLongitude, yCoord = occDataRaw$decimalLatitude, taxonKey = occDataRaw$taxonKey)))
	dim(occCounts) <- c(nrow(coordinates(inTemplateGrid)), length(inSpeciesIDs))
	colnames(occCounts) <- as.character(inSpeciesIDs)
	# Define which cells in the template grid have output values
	haveData <- rep(TRUE, nrow(occCounts))
	if(is(inTemplateGrid, "SpatialGridDataFrame")) {
		# If the template grid has a data frame attached then retrieve any rows that have data
		haveData <- apply(X = as.matrix(templateGrid@data), FUN = function(curRow) {
			!all(is.na(curRow))
		}, MARGIN = 1)
	}
	# Set the relevant rows of the occurrence counts to NA
	occCounts[!haveData, ] <- NA
	# Clear up the folder after processing
	if(unlink(inOutFolder, recursive = TRUE) != 0) {
		warning("unable to delete downloaded occupancy data")
	}
	# Process the raw data
	rawOccurrences <- SpatialPointsDataFrame(
		# Assign coordinates based on the decimal longitude and latitude values
		cbind(occDataRaw$decimalLongitude, occDataRaw$decimalLatitude),
		# Assign all points that are not the coordinates to be attributes in the attribute table
		occDataRaw[, !(colnames(occDataRaw) %in% c("decimalLongitude", "decimalLatitude"))],
		proj4string = GBIFcrs
	)
	# Transform the raw occurrences if they are on a different projection than the template grid
	if(proj4string(templateGrid) != GBIFcrsString) {
		rawOccurrences <- spTransform(rawOccurrences, CRS(proj4string(templateGrid)))
	}
	# Remove points that do not fall on cells with data
	rawOccurrences <- rawOccurrences[haveData[getGridIndex(coordinates(rawOccurrences), getGridTopology(inTemplateGrid), all.inside = FALSE)], ]
	# Return the gridded and raw occurrence data
	list(
		rawOccurrences = rawOccurrences,
		griddedOccurrences = SpatialGridDataFrame(getGridTopology(templateGrid), as.data.frame(occCounts), CRS(proj4string(templateGrid)))
	)
}