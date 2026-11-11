# Estimate point elevations by equal-weight least squares.
# Based on Condit's solve.topo(), using sparse QR for the observation equations.
# Requires Matrix and igraph. Returns x, y, z sorted by x and then y.
# Each row is one measurement of dz = z_j - z_i.
# A point is identified by its exact (x, y) coordinates.
estimate_point_elevations <- function(survey, reference.xy = NULL, reference.z = 0)
{
	##################################################
	# READ AND VALIDATE INPUTS
	##################################################

	# Select columns by name so their order in the input does not matter.
	columns <- c("x_i", "y_i", "x_j", "y_j", "dz")
	if(!is.data.frame(survey) && !is.matrix(survey))
		stop("survey must be a data frame or a matrix.")
	if(!all(columns %in% colnames(survey)) ||
	   anyDuplicated(colnames(survey)[colnames(survey) %in% columns]))
		stop("survey must contain one of each column: x_i, y_i, x_j, y_j, dz.")
	survey <- as.data.frame(survey)[columns]

	complete <- complete.cases(survey)
	if(!any(complete)) stop("No complete elevation differences.")
	if(!all(vapply(survey, function(x) is.numeric(x) && is.null(dim(x)), logical(1))) ||
	   any(is.infinite(as.matrix(survey))))
		stop("Survey columns must be numeric, with finite values or NA.")
	if(!is.numeric(reference.z) || length(reference.z) != 1L || !is.finite(reference.z))
		stop("reference.z must be one finite number.")


	##################################################
	# SELECT USABLE OBSERVATIONS
	##################################################

	# Exclude incomplete rows and observations whose endpoints coincide.
	same.point <- complete & survey$x_i == survey$x_j & survey$y_i == survey$y_j
	keep <- complete & !same.point
	if(!any(keep)) stop("No usable elevation differences.")
	if(any(!keep))
		warning(sum(!keep), " incomplete or same-point observation(s) omitted.", call. = FALSE)
	survey <- survey[keep, , drop = FALSE]


	##################################################
	# IDENTIFY POINTS AND CHECK CONNECTIONS
	##################################################

	# Stack the starting and ending coordinates, then sort by x and y.
	n <- nrow(survey)
	xy <- data.frame(x = c(survey$x_i, survey$x_j), y = c(survey$y_i, survey$y_j))
	sorted <- order(xy$x, xy$y)
	xy <- xy[sorted, , drop = FALSE]

	# Identify distinct points by comparing adjacent coordinate pairs exactly.
	last <- nrow(xy)
	new.point <- c(TRUE, xy$x[-1L] != xy$x[-last] | xy$y[-1L] != xy$y[-last])

	# Map point IDs back to the original endpoint order to preserve alignment with dz.
	id <- integer(2L * n)
	id[sorted] <- cumsum(new.point)

	# Keep one output row per distinct point, ordered by x and y.
	points <- xy[new.point, , drop = FALSE]
	row.names(points) <- NULL

	from <- id[seq_len(n)]
	to <- id[n + seq_len(n)]

	# Check that usable measurements connect all retained points.
	graph <- igraph::graph_from_edgelist(cbind(from, to), directed = FALSE)
	if(!igraph::is_connected(graph))
		stop("Usable observations form disconnected groups. One datum cannot tie them together.")


	##################################################
	# SET THE REFERENCE POINT
	##################################################

	# Default to the first point in (x, y) order, independently of input row order.
	anchor <- 1L

	# If reference coordinates are supplied, locate that exact point.
	if(!is.null(reference.xy))
	{
		if(!is.numeric(reference.xy) || length(reference.xy) != 2L ||
		   any(!is.finite(reference.xy)))
			stop("reference.xy must contain two finite numbers.")
		anchor <- which(points$x == reference.xy[1L] & points$y == reference.xy[2L])
		if(!length(anchor)) stop("reference.xy must occur in a usable observation.")
	}


	##################################################
	# FIT ELEVATIONS AND RETURN RESULTS
	##################################################

	# Every measurement, including repeats, contributes one equation for z_j - z_i.
	observation <- seq_len(n)
	A <- Matrix::sparseMatrix(i = c(observation, observation),
	                          j = c(from, to),
	                          x = rep(c(-1, 1), each = n),
	                          dims = c(n, nrow(points)))

	# Fix the reference point at zero by removing its column, then fit by sparse QR.
	fit <- Matrix::qr(A[, -anchor, drop = FALSE])

	# Shift every fitted elevation by the reference elevation.
	points$z <- rep(reference.z, nrow(points))
	points$z[-anchor] <- as.numeric(Matrix::qr.coef(fit, survey$dz)) + reference.z

	# Return one {x, y, z} row per point in the usable observations.
	points
}
