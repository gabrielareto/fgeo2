# Recover point coordinates with an affine transformation.
#
# known.xy.of.controls contains known plot coordinates. The corresponding
# clicked controls must have the same row order. The result uses the units and
# coordinate system of the known controls. At least three non-collinear
# controls are required.

apply_affine_transformation <- function(known.xy.of.controls,
											clicked.xy.of.controls,
											clicked.xy.of.points)
{
	inputs <- list(
		known.xy.of.controls,
		clicked.xy.of.controls,
		clicked.xy.of.points
	)
	valid <- sapply(inputs, function(x)
		is.matrix(x) && is.numeric(x) &&
		ncol(x) == 2 && all(is.finite(x)))

	if(any(!valid))
		stop("All coordinate inputs must be numeric matrices with two finite columns.")

	number.controls <- nrow(known.xy.of.controls)
	if(number.controls < 3 ||
	   number.controls != nrow(clicked.xy.of.controls))
		stop("Affine control matrices must contain the same three or more points.")

	# Fit the transformation from known to clicked coordinates, then invert it.
	coefficients <- qr.solve(
		cbind(1, known.xy.of.controls),
		clicked.xy.of.controls
	)

	estimated.xy.of.points <- cbind(
		clicked.xy.of.points[, 1] - coefficients[1, 1],
		clicked.xy.of.points[, 2] - coefficients[1, 2]
	) %*% solve(coefficients[-1, , drop = FALSE])

	colnames(estimated.xy.of.points) <- c("x", "y")
	rownames(estimated.xy.of.points) <- rownames(clicked.xy.of.points)
	estimated.xy.of.points
}


# Recover point coordinates with a homography.
#
# known.xy.of.controls contains known plot coordinates. The corresponding
# clicked controls must have the same row order. The result uses the units and
# coordinate system of the known controls. At least four non-degenerate
# controls are required.

apply_homography_transformation <- function(known.xy.of.controls,
											 clicked.xy.of.controls,
											 clicked.xy.of.points)
{
	inputs <- list(
		known.xy.of.controls,
		clicked.xy.of.controls,
		clicked.xy.of.points
	)
	valid <- sapply(inputs, function(x)
		is.matrix(x) && is.numeric(x) &&
		ncol(x) == 2 && all(is.finite(x)))

	if(any(!valid))
		stop("All coordinate inputs must be numeric matrices with two finite columns.")

	number.controls <- nrow(known.xy.of.controls)
	if(number.controls < 4 ||
	   number.controls != nrow(clicked.xy.of.controls))
		stop("Homography control matrices must contain the same four or more points.")

	x <- clicked.xy.of.controls[, 1]
	y <- clicked.xy.of.controls[, 2]
	gx <- known.xy.of.controls[, 1]
	gy <- known.xy.of.controls[, 2]

	# Each control contributes two equations for eight coefficients.
	equations <- rbind(
		cbind(x, y, 1, 0, 0, 0, -gx*x, -gx*y),
		cbind(0, 0, 0, x, y, 1, -gy*x, -gy*y)
	)
	coefficients <- qr.solve(equations, c(gx, gy))

	# Transform the points and convert from homogeneous coordinates.
	homography <- matrix(c(coefficients, 1), nrow = 3, byrow = TRUE)
	estimated.xy.of.points <- cbind(clicked.xy.of.points, 1) %*%
		t(homography)
	estimated.xy.of.points <- estimated.xy.of.points[, 1:2, drop = FALSE] /
		estimated.xy.of.points[, 3]

	colnames(estimated.xy.of.points) <- c("x", "y")
	rownames(estimated.xy.of.points) <- rownames(clicked.xy.of.points)
	estimated.xy.of.points
}
