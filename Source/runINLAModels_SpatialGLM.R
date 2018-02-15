## 1. ------ FUNCTION TO MODEL OCCUPANCY USING INLA SPATIAL GLM ------
runINLASpatialGLM <- function(occurrenceData, climateData, alpha = 1.5, meshParameters = list(), progressUpdater = NULL) {
	### 1.1 ==== Sanity check the function inputs ====
	updateProgressBar(progressUpdater, value = 0.0, details = "processing inputs")
	# Check the occurrence data
	inOccurrenceData <- tryCatch(as(occurrenceData, "SpatialGridDataFrame"), error = function(err) {
		stop(paste("invalid input for the occurrence data:", err, sep = " "))
	})
	# Check the climate data
	inClimateData <- tryCatch(as(climateData, "SpatialGridDataFrame"), error = function(err) {
		stop(paste("invalid input for the environmental covariates:", err, sep = " "))
	})
	# Ensure that they have the same grid topology and coordinate reference system
	if(!identical(as(inOccurrenceData, "SpatialGrid"), as(inClimateData, "SpatialGrid"))) {
		stop("occurrence data and environmental covariate data do not share the same grid topology and/or coordinate projection system")
	}
	# Test the value of the fractional operator order
	inAlpha <- tryCatch(as.double(alpha), error = function(err) {
		stop(paste("invalid input for the fractional operator order:", err, sep = " "))
	})
	if(is.null(inAlpha) || length(inAlpha) <= 0) {
		stop("invalid input for the fractional operator order: vector has zero length")
	} else if(length(inAlpha) > 1) {
		warning("fractional operator order specification vector length greater than one: only the first element will be used")
		inAlpha <- inAlpha[1]
	}
	if(anyNA(inAlpha)) {
		stop("invalid input for the fractional operator order: NA values detected")
	}
	# Check the parameters to generate the spatial mesh
	inMeshParameters <- tryCatch(as.list(meshParameters), error = function(err) {
		stop(paste("invalid input for the mesh parameters:", err, sep = " "))
	})
	# Set the mesh parameters from the input list and use default values if none supplied
	inOffset <- meshParameters$offset
	inN <- meshParameters$n
	inMaxEdge <- meshParameters$max.edge
	inMinAngle <- meshParameters$min.angle
	inCutoff <- meshParameters$cutoff
	inMaxNStrict <- meshParameters$max.n.strict
	inMaxN <- meshParameters$max.n
	inConvex <- meshParameters$convex
	if(is.null(inConvex)) {
		inConvex <- -0.15
	}
	inConcave <- meshParameters$concave
	if(is.null(inConcave)) {
		inConcave <- inConvex
	}
	inResolution <- meshParameters$resolution
	if(is.null(inResolution)) {
		inResolution <- 40
	}
	inEps <- meshParameters$eps
	### 1.2 ==== Generate the spatial random effects mesh ====
	updateProgressBar(progressUpdater, value = 0.2, details = "generating mesh for spatial random effect")
	# Determine which rows to use
	useRow <- !apply(X = as.matrix(inClimateData@data), FUN = anyNA, MARGIN = 1)
	# Create a bounding hull from the climate data coordinates
	boundHull <- inla.nonconvex.hull(points = coordinates(inClimateData)[useRow, ], convex = inConvex, concave = inConcave, resolution = inResolution, eps = inEps, crs = CRS(proj4string(inClimateData)))
	# Create a mesh to define the spatial random effects on
	effectsMesh <- inla.mesh.2d(boundary = boundHull, n = inN, offset = inOffset, max.edge = inMaxEdge, min.angle = inMinAngle, cutoff = inCutoff, max.n.strict = inMaxNStrict, max.n = inMaxN, crs = CRS(proj4string(inClimateData)))
	# Create a projection matrix associated with the mesh
	effectsProjMat <- inla.spde.make.A(mesh = effectsMesh, loc = coordinates(inClimateData)[useRow, ])
	### 1.3 ==== Formulate the model specification ====
	
}