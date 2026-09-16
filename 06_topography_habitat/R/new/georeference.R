# Recover point coordinates with an affine transformation.
#
# known.xy.of.controls contains the known plot coordinates of the controls.
# clicked.xy.of.controls contains their coordinates in the image, in the same
# row order. clicked.xy.of.points contains the image coordinates to recover.
# The result is an x-y matrix in the units and coordinate system of the known
# controls. At least three non-collinear controls are required.

apply_affine_transformation <- function(known.xy.of.controls,
											clicked.xy.of.controls,
											clicked.xy.of.points)
{
	known.xy.of.controls <- as.matrix(known.xy.of.controls)
	clicked.xy.of.controls <- as.matrix(clicked.xy.of.controls)
	clicked.xy.of.points <- as.matrix(clicked.xy.of.points)

	inputs <- list(
		known.xy.of.controls = known.xy.of.controls,
		clicked.xy.of.controls = clicked.xy.of.controls,
		clicked.xy.of.points = clicked.xy.of.points
	)

	if(any(!vapply(inputs, is.numeric, logical(1))) ||
	   any(vapply(inputs, ncol, integer(1)) != 2L))
		stop("All coordinate inputs must have exactly two numeric columns.")

	if(any(!vapply(inputs, function(x) all(is.finite(x)), logical(1))))
		stop("Coordinates cannot contain NA, NaN, Inf or -Inf.")

	if(nrow(known.xy.of.controls) != nrow(clicked.xy.of.controls))
		stop("The two control matrices must have the same number of rows.")

	if(nrow(known.xy.of.controls) < 3L)
		stop("An affine transformation requires at least three controls.")

	if(nrow(clicked.xy.of.points) < 1L)
		stop("At least one point coordinate is required.")

	known.center <- colMeans(known.xy.of.controls)
	clicked.center <- colMeans(clicked.xy.of.controls)
	known.centered <- sweep(known.xy.of.controls, 2, known.center)
	clicked.centered <- sweep(clicked.xy.of.controls, 2, clicked.center)

	if(qr(known.centered)$rank < 2L)
		stop("The known controls must not be collinear.")

	# Fit known-to-clicked coordinates because the known locations are fixed
	# and the clicks contain measurement error. Centering improves stability
	# when plot coordinates have a large global offset.
	coefficients <- qr.solve(known.centered, clicked.centered)

	if(qr(coefficients)$rank < 2L)
		stop("The clicked controls do not define an invertible transformation.")

	estimated.xy.of.points <- sweep(
		clicked.xy.of.points,
		2,
		clicked.center
	) %*% solve(coefficients)
	estimated.xy.of.points <- sweep(
		estimated.xy.of.points,
		2,
		known.center,
		"+"
	)

	colnames(estimated.xy.of.points) <- c("x", "y")
	rownames(estimated.xy.of.points) <- rownames(clicked.xy.of.points)
	estimated.xy.of.points
}


# Recover point coordinates with a planar projective transformation.
#
# known.xy.of.controls contains the known plot coordinates of the controls.
# clicked.xy.of.controls contains their coordinates in the image, in the same
# row order. clicked.xy.of.points contains the image coordinates to recover.
# The result is an x-y matrix in the units and coordinate system of the known
# controls. At least four non-degenerate controls are required.

apply_homography_transformation <- function(known.xy.of.controls,
											 clicked.xy.of.controls,
											 clicked.xy.of.points)
{
	known.xy.of.controls <- as.matrix(known.xy.of.controls)
	clicked.xy.of.controls <- as.matrix(clicked.xy.of.controls)
	clicked.xy.of.points <- as.matrix(clicked.xy.of.points)

	inputs <- list(
		known.xy.of.controls = known.xy.of.controls,
		clicked.xy.of.controls = clicked.xy.of.controls,
		clicked.xy.of.points = clicked.xy.of.points
	)

	if(any(!vapply(inputs, is.numeric, logical(1))) ||
	   any(vapply(inputs, ncol, integer(1)) != 2L))
		stop("All coordinate inputs must have exactly two numeric columns.")

	if(any(!vapply(inputs, function(x) all(is.finite(x)), logical(1))))
		stop("Coordinates cannot contain NA, NaN, Inf or -Inf.")

	if(nrow(known.xy.of.controls) != nrow(clicked.xy.of.controls))
		stop("The two control matrices must have the same number of rows.")

	if(nrow(known.xy.of.controls) < 4L)
		stop("A homography requires at least four controls.")

	if(nrow(clicked.xy.of.points) < 1L)
		stop("At least one point coordinate is required.")

	# Center and scale both coordinate systems before fitting. This leaves the
	# result unchanged but improves numerical stability for global coordinates.
	known.center <- colMeans(known.xy.of.controls)
	clicked.center <- colMeans(clicked.xy.of.controls)
	known.centered <- sweep(known.xy.of.controls, 2, known.center)
	clicked.centered <- sweep(clicked.xy.of.controls, 2, clicked.center)
	known.scale <- sqrt(mean(known.centered^2))
	clicked.scale <- sqrt(mean(clicked.centered^2))

	if(qr(known.centered)$rank < 2L || qr(clicked.centered)$rank < 2L)
		stop("Known and clicked controls must each span two dimensions.")

	if(!is.finite(known.scale) || known.scale == 0 ||
	   !is.finite(clicked.scale) || clicked.scale == 0)
		stop("The controls do not span two coordinate systems.")

	known.scaled <- known.centered / known.scale
	clicked.scaled <- clicked.centered / clicked.scale
	x <- clicked.scaled[, 1]
	y <- clicked.scaled[, 2]
	gx <- known.scaled[, 1]
	gy <- known.scaled[, 2]

	# Each control contributes two equations for eight coefficients.
	equations <- rbind(
		cbind(x, y, 1, 0, 0, 0, -gx*x, -gx*y),
		cbind(0, 0, 0, x, y, 1, -gy*x, -gy*y)
	)

	if(qr(equations)$rank < 8L)
		stop("The controls do not define a unique homography.")

	coefficients <- qr.solve(equations, c(gx, gy))
	homography <- matrix(c(coefficients, 1), nrow = 3, byrow = TRUE)
	control.denominator <- (cbind(clicked.scaled, 1) %*%
		t(homography))[, 3]

	if(any(!is.finite(control.denominator) |
	       abs(control.denominator) < sqrt(.Machine$double.eps)) ||
	   length(unique(sign(control.denominator))) > 1L)
		stop("Projective horizon crosses controls; check control identities.")

	points.scaled <- sweep(
		clicked.xy.of.points,
		2,
		clicked.center
	) / clicked.scale
	points.homogeneous <- cbind(points.scaled, 1) %*% t(homography)
	denominator <- points.homogeneous[, 3]

	if(any(!is.finite(denominator) |
	       abs(denominator) < sqrt(.Machine$double.eps)) ||
	   any(sign(denominator) != sign(control.denominator[1])))
		stop("At least one point is on or across the projective horizon.")

	estimated.xy.of.points <- points.homogeneous[, 1:2, drop = FALSE] /
		denominator
	estimated.xy.of.points <- estimated.xy.of.points * known.scale
	estimated.xy.of.points <- sweep(
		estimated.xy.of.points,
		2,
		known.center,
		"+"
	)

	if(any(!is.finite(estimated.xy.of.points)))
		stop("The homography produced non-finite coordinates.")

	colnames(estimated.xy.of.points) <- c("x", "y")
	rownames(estimated.xy.of.points) <- rownames(clicked.xy.of.points)
	estimated.xy.of.points
}
