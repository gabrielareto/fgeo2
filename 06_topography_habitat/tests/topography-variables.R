# Standalone checks for derive_topography().
# Run from the repository root or from 06_topography_habitat/tests.
# Core checks need terra. Whitebox checks run only when its program is installed.

path.to.fgeo2 = if(file.exists("fgeo2.Rproj")) "." else "../.."
source(file.path(path.to.fgeo2, "06_topography_habitat/R/new/topography-variables.R"))

local({
	library(terra)

	##################################################
	# ANALYTIC PLANE AND UNCHANGED DEM INPUT
	##################################################

	# On a 5 m grid, this plane has slope atan(sqrt(.1^2 + .2^2)).
	# Its downhill direction is southwest, TPI is zero, and TRI is 0.875 m.
	plane.dem <- rast(xmin = 0, xmax = 100, ymin = 0, ymax = 100,
	                  resolution = 5, crs = "local")
	cell.xy <- xyFromCell(plane.dem, seq_len(ncell(plane.dem)))
	values(plane.dem) <- 100 + 0.1 * cell.xy[, 1] + 0.2 * cell.xy[, 2]
	original.values <- values(plane.dem)
	original.names <- names(plane.dem)
	variables = c("slope", "aspect", "northness", "eastness", "tpi", "tri")
	plane.layers <- derive_topography(dem = plane.dem, variables = variables,
	                                  position.radius = 10)
	expected.slope <- atan(sqrt(0.1^2 + 0.2^2)) * 180 / pi
	expected.aspect <- (atan2(-0.1, -0.2) * 180 / pi) %% 360
	centre <- cellFromRowCol(plane.dem, 10, 10)
	stopifnot(
		all(is.finite(values(plane.layers)[centre, ])),
		identical(names(plane.layers), c("elevation_m", "slope_deg", "aspect_deg",
		                                "northness", "eastness", "tpi_m", "tri_m")),
		compareGeom(plane.layers, plane.dem, stopOnError = FALSE),
		identical(names(plane.dem), original.names),
		identical(values(plane.dem), original.values),
		max(abs(values(plane.layers[["elevation_m"]]) - original.values)) == 0,
		max(abs(values(plane.layers[["slope_deg"]]) - expected.slope), na.rm = TRUE) < 1e-8,
		max(abs(values(plane.layers[["aspect_deg"]]) - expected.aspect), na.rm = TRUE) < 1e-8,
		max(abs(values(plane.layers[["northness"]]) + 0.2 / sqrt(0.05)), na.rm = TRUE) < 1e-8,
		max(abs(values(plane.layers[["eastness"]]) + 0.1 / sqrt(0.05)), na.rm = TRUE) < 1e-8,
		max(abs(values(plane.layers[["tpi_m"]])), na.rm = TRUE) < 1e-8,
		max(abs(values(plane.layers[["tri_m"]]) - 0.875), na.rm = TRUE) < 1e-8
	)

	# A flat surface has no meaningful aspect, but zero slope/TPI/TRI.
	flat.dem <- plane.dem
	values(flat.dem) <- 100
	flat <- derive_topography(dem = flat.dem, variables = variables, position.radius = 10)
	stopifnot(all(is.na(values(flat[["aspect_deg"]]))),
	          all(is.na(values(flat[["northness"]]))),
	          all(is.na(values(flat[["eastness"]]))),
	          max(abs(values(flat[["slope_deg"]])), na.rm = TRUE) == 0,
	          max(abs(values(flat[["tpi_m"]])), na.rm = TRUE) == 0,
	          max(abs(values(flat[["tri_m"]])), na.rm = TRUE) == 0)

	##################################################
	# POINT INPUTS, ORDER, DUPLICATES AND MISSING ROWS
	##################################################

	xyz <- expand.grid(x = seq(0, 100, 20), y = seq(0, 100, 20))
	xyz$z <- 100 + 0.1 * xyz$x + 0.2 * xyz$y
	points.dem <- derive_topography(xyz = xyz, variables = "elevation")
	reordered <- derive_topography(xyz = xyz[nrow(xyz):1, c("z", "y", "x")],
	                                variables = "elevation")
	duplicated <- derive_topography(xyz = rbind(xyz, xyz[1:3, ]), variables = "elevation")
	unnamed.matrix <- unname(as.matrix(xyz[c("x", "y", "z")]))
	matrix.dem <- derive_topography(xyz = unnamed.matrix, variables = "elevation")
	stopifnot(all(res(points.dem) == 5),
	          all(as.vector(ext(points.dem)) == c(0, 100, 0, 100)),
	          identical(names(points.dem), "elevation_m"))
	for(other in list(reordered, duplicated, matrix.dem))
		stopifnot(compareGeom(points.dem, other, stopOnError = FALSE),
		          isTRUE(all.equal(values(points.dem), values(other), tolerance = 1e-10)))

	warning.seen <- FALSE
	missing.row <- withCallingHandlers(
		derive_topography(xyz = rbind(xyz, c(NA_real_, NA_real_, NA_real_)),
		                   variables = "elevation"),
		warning = function(w)
		{
			warning.seen <<- TRUE
			invokeRestart("muffleWarning")
		}
	)
	stopifnot(warning.seen,
	          isTRUE(all.equal(values(points.dem), values(missing.row), tolerance = 1e-10)))

	# A complete regular grid can use bilinear interpolation without interp.
	bilinear <- derive_topography(xyz = xyz, variables = "elevation", interpolation = "bilinear")
	stopifnot(all(is.finite(values(bilinear))),
	          max(abs(values(bilinear) - (100 + 0.1 * cell.xy[, 1] + 0.2 * cell.xy[, 2]))) < 1e-5)

	# Unavailable support stays missing when the DEM is smoothed.
	hull.dem <- plane.dem
	values(hull.dem)[cell.xy[, 1] + cell.xy[, 2] > 100] <- NA
	smoothed <- derive_topography(dem = hull.dem, variables = "elevation", smoothing.cells = 3)
	stopifnot(identical(unname(is.na(values(hull.dem))), unname(is.na(values(smoothed)))))

	##################################################
	# INVALID OR AMBIGUOUS INPUTS
	##################################################

	geographic.dem <- plane.dem
	crs(geographic.dem) <- "EPSG:4326"
	feet.dem <- plane.dem
	crs(feet.dem) <- "EPSG:2263"
	unknown.crs <- plane.dem
	crs(unknown.crs) <- ""
	rectangular.dem <- rast(xmin = 0, xmax = 100, ymin = 0, ymax = 100,
	                        nrows = 10, ncols = 20, crs = "local")
	values(rectangular.dem) <- 100
	conflicting <- rbind(xyz, transform(xyz[1, ], z = z + 1))
	infinite <- xyz
	infinite$z[1] <- Inf

	invalid.calls <- list(
		quote(derive_topography()),
		quote(derive_topography(xyz = xyz, dem = plane.dem)),
		quote(derive_topography(dem = c(plane.dem, plane.dem))),
		quote(derive_topography(dem = geographic.dem, variables = "slope")),
		quote(derive_topography(dem = feet.dem, variables = "slope")),
		quote(derive_topography(dem = unknown.crs, variables = "slope")),
		quote(derive_topography(dem = rectangular.dem, variables = "slope")),
		quote(derive_topography(dem = plane.dem, variables = "misspelled")),
		quote(derive_topography(dem = plane.dem, variables = c("all", "slope"))),
		quote(derive_topography(dem = plane.dem, resolution = 10)),
		quote(derive_topography(dem = plane.dem, extent = c(0, 100, 0, 100))),
		quote(derive_topography(dem = plane.dem, crs = "local")),
		quote(derive_topography(xyz = xyz, variables = "elevation", resolution = 0)),
		quote(derive_topography(xyz = xyz, variables = "elevation", extent = c(0, 102, 0, 100))),
		quote(derive_topography(xyz = xyz, variables = "elevation", interpolation = "misspelled")),
		quote(derive_topography(dem = plane.dem, variables = "slope", smoothing.cells = 2)),
		quote(derive_topography(dem = plane.dem, variables = "slope", neighbors = 6)),
		quote(derive_topography(dem = plane.dem, variables = "slope", neighbors = "8")),
		quote(derive_topography(xyz = xyz, variables = "elevation", crs = "EPSG:4326")),
		quote(derive_topography(dem = plane.dem, variables = "twi", flow.preparation = "misspelled")),
		quote(derive_topography(xyz = conflicting, variables = "elevation")),
		quote(derive_topography(xyz = infinite, variables = "elevation")),
		quote(derive_topography(xyz = xyz[1, ] * NA_real_, variables = "elevation")),
		quote(derive_topography(xyz = data.frame(x = 0:4, y = 0:4, z = 0:4), variables = "elevation"))
	)
	for(call in invalid.calls)
		stopifnot(inherits(suppressWarnings(try(eval(call), silent = TRUE)), "try-error"))

	cat("Core terrain checks passed.\n")

	##################################################
	# OPTIONAL WHITEBOX, LIFECYCLE AND SAVED RASTERS
	##################################################

	has.whitebox <- requireNamespace("whitebox", quietly = TRUE) &&
	                file.exists(whitebox::wbt_exe_path(shell_quote = FALSE))
	if(has.whitebox)
	{
		# Local-only calculations must not depend on an available Whitebox executable.
		old.exe <- whitebox::wbt_exe_path(shell_quote = FALSE)
		tryCatch(
			{
				whitebox::wbt_exe_path(exe_path = tempfile("unavailable-whitebox-"))
				local.only <- derive_topography(dem = plane.dem, variables = "slope")
				stopifnot(identical(names(local.only), c("elevation_m", "slope_deg")))
			},
			finally = whitebox::wbt_exe_path(exe_path = old.exe)
		)

		# Zero-slope TWI is undefined; the function must retain NA rather than
		# silently replacing zero slope with an arbitrary small positive number.
		flat.twi <- derive_topography(dem = flat.dem, variables = "twi")
		stopifnot(all(is.na(values(flat.twi[["twi"]]))),
		          all(values(flat.twi[["elevation_m"]]) == 100))
		filled.twi <- derive_topography(dem = plane.dem, variables = "twi", flow.preparation = "fill")
		stopifnot(any(is.finite(values(filled.twi[["twi"]]))))

		records <- read.csv(file.path(path.to.fgeo2, "data_examples/topography-survey.csv"))
		survey <- data.frame(x_i = records$from_x_m, y_i = records$from_y_m,
		                     x_j = records$to_x_m, y_j = records$to_y_m,
		                     dz = records$height_difference_m)
		source(file.path(path.to.fgeo2, "06_topography_habitat/R/new/topography-solver.R"))
		elevations <- estimate_point_elevations(survey, reference.xy = c(0, 0), reference.z = 371)
		work.before <- list.files(tempdir(), pattern = "^fgeo2-topography-work-", full.names = TRUE)
		all.layers <- derive_topography(xyz = elevations, smoothing.cells = 3)
		stopifnot(identical(work.before,
		          list.files(tempdir(), pattern = "^fgeo2-topography-work-", full.names = TRUE)))
		stopifnot(nlyr(all.layers) == 13L, nrow(all.layers) == 48L, ncol(all.layers) == 48L,
		          identical(names(all.layers)[1], "elevation_m"),
		          min(values(all.layers[["drainage_distance_m"]]), na.rm = TRUE) >= 0,
		          min(values(all.layers[["height_above_drainage_m"]]), na.rm = TRUE) >= 0)
		for(variable in c("slope", "tpi", "twi", "height_above_drainage"))
		{
			selected <- derive_topography(xyz = elevations, variables = variable, smoothing.cells = 3)
			stopifnot(nlyr(selected) == 2L,
			          isTRUE(all.equal(values(selected), values(all.layers[[names(selected)]]), tolerance = 1e-8)))
		}

		# A later call must not overwrite or invalidate the earlier returned rasters.
		before <- values(all.layers)
		shifted <- elevations
		shifted$z <- shifted$z + 100
		again <- derive_topography(xyz = shifted, smoothing.cells = 3)
		gc()
		stopifnot(isTRUE(all.equal(before, values(all.layers), tolerance = 0)),
		          max(abs(values(again[["elevation_m"]]) - before[, 1] - 100)) < 1e-8)
		paths <- sources(all.layers)
		stopifnot(all(file.exists(paths[nzchar(paths)])))

		# Cleanup also runs when an impossible channel threshold stops the call.
		no.channels <- try(derive_topography(xyz = elevations, variables = "drainage_distance",
		                                     channel.area.m2 = 1e12), silent = TRUE)
		stopifnot(inherits(no.channels, "try-error"),
		          identical(work.before,
		          list.files(tempdir(), pattern = "^fgeo2-topography-work-", full.names = TRUE)))

		# An already draining plane can skip conditioning entirely.
		unconditioned <- derive_topography(dem = plane.dem, variables = "contributing_area",
		                                   flow.preparation = "none")
		stopifnot(any(is.finite(values(unconditioned[["contributing_area_m2"]]))))

		file <- tempfile(fileext = ".tif")
		writeRaster(all.layers, file, datatype = "FLT8S", overwrite = TRUE)
		restored <- rast(file)
		stopifnot(compareGeom(all.layers, restored, stopOnError = FALSE),
		          identical(names(all.layers), names(restored)),
		          identical(units(all.layers), units(restored)),
		          isTRUE(all.equal(values(all.layers), values(restored), tolerance = 1e-10)))
		unlink(file)
		cat("Whitebox, variable selection, temporary-file lifecycle and GeoTIFF checks passed.\n")
	}
	if(!has.whitebox) cat("Whitebox checks skipped: executable not installed.\n")
})
