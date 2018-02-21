## 1. ------ FUNCTION TO CONVERT INLA MESH OBJECTS TO SPATIAL POLYGONS ------
#' @title Converts the Delaunay Triangulation in INLA Mesh Objects to SpatialPolygons
#' 
#' @description 
#' Converts the constructed Delaunay triangulation from a two-dimensional INLA \link[INLA::inla.mesh.create]{mesh} object
#' into a \code{\link[sp::SpatialPolygons]{SpatialPolygons}} object.  This can be useful if you want to save the mesh
#' in standard vector file formats or perform further processing on the mesh.
#' 
#' @param meshInput An \code{\link[INLA::inla.mesh.create]{inla.mesh}} object.  This mesh must be two-dimensional.
#' 
#' @return A \code{\link[sp::SpatialPolygons]{SpatialPolygons}} object containing the Delaunay triangulation from
#' the input mesh.  Note that inner and outer boundary information is not present in the output.
#' 
#' @examples
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[sp::SpatialPolygons]{SpatialPolygons}} \code{\link[sp::SpatialPolygons]{SpatialPolygons}}
#' @export
#'
inlaMeshToSpatialPolygons <- function(meshInput) {
	# Sanity test the input
	inMeshInput <- tryCatch(as(meshInput, "inla.mesh"), error = function(err) {
		stop(paste("invalid input mesh:", err, sep = " "))
	})
	# Ensure the mesh has the correct manifold type
	if(inMeshInput$manifold != "R2") {
		stop("invalid input mesh: mesh must be two-dimensional")
	}
	# Retrieve the Delaunay triangles as polygon objects
	delaunayPolygons <- lapply(X = 1:nrow(inMeshInput$graph$tv), FUN = function(curRowIndex, pointCoords, pointIndeces) {
		# Retrieve the vertex coordinates of the current triangle
		curIndeces <- pointIndeces[curRowIndex, ]
		# Construct a Polygons object to contain the triangle
		Polygons(list(Polygon(
			pointCoords[c(curIndeces, curIndeces[1]), ], hole = FALSE
		)), ID = paste("Delaunay", curRowIndex, sep = "_"))
	}, pointCoords = inMeshInput$loc[, c(1, 2)], pointIndeces = inMeshInput$graph$tv)
	# Convert the Polygons into a SpatialPolygons object
	SpatialPolygons(delaunayPolygons, proj4string = inMeshInput$crs)
}