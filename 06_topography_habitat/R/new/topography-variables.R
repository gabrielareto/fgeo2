# Derive topographical raster layers from point elevations or an existing DEM.
#
# INPUT AND OUTPUT
# Supply exactly one input:
# - xyz: a data frame or matrix containing x, y and z in metres. An unnamed
#   three-column matrix is interpreted in that order. Extra columns are ignored.
#   Incomplete rows are omitted with a warning. Identical records are collapsed;
#   conflicting elevations at the same coordinates must be resolved by the user.
# - dem: one SpatRaster elevation layer, with square cells and a metre-based
#   planar CRS. Its grid is retained. Interpolation is skipped.
#
# Returns a standard SpatRaster: elevation first, then the requested variables.
# Names and units are retained. No plots or permanent files are produced.
# WhiteboxTools uses a private temporary directory. Requested results are copied
# to a separate R temporary GeoTIFF before that directory is removed. Save the
# returned raster with terra::writeRaster() to retain it across R sessions.
#
# VARIABLES
# variables = "all" returns these thirteen layers. Alternatively supply a vector
# of the keys on the left; elevation is always included, and dependencies are
# calculated only as needed. Unrequested intermediate layers are not returned.
# elevation              -> elevation_m              : elevation (m)
# slope                  -> slope_deg                : inclination (degrees)
# aspect                 -> aspect_deg               : downhill direction (degrees)
# northness, eastness    -> northness, eastness        : cos/sin of aspect
# tpi                    -> tpi_m                    : focal minus neighbour mean (m)
# tri                    -> tri_m                    : mean absolute elevation
#                                                      difference to 8 neighbours (m)
# plan_curvature         -> plan_curvature_per_m      : positive = divergence (1/m)
# profile_curvature      -> profile_curvature_per_m   : positive = downhill
#                                                      flow acceleration (1/m)
# contributing_area      -> contributing_area_m2      : D-infinity upslope area (m2)
# twi                    -> twi                      : log(SCA / tan(slope))
# drainage_distance      -> drainage_distance_m       : D8 downslope path length (m)
# height_above_drainage  -> height_above_drainage_m    : height above the channel
#                                                      reached along D8 flow (m)
#
# PARAMETERS
# resolution: square output cell size for xyz (m). Default 5; not survey accuracy.
# extent: calculation extent c(xmin,xmax,ymin,ymax) for xyz. Width and height must
#   be divisible by resolution. NULL uses the point bounds, rounded outwards to
#   multiples of resolution. Keep surrounding terrain for hydrological analyses.
# crs: actual CRS of xyz; "local" means an arbitrary Cartesian plane in metres.
#   This labels the supplied coordinates; it does not transform them.
# interpolation: "idw" (terra), "linear" or "akima" (optional interp package),
#   or "bilinear" for a complete regular xyz grid. No silent method fallback.
# search.radius, idw.power, idw.smooth: IDW search distance (m), distance power,
#   and distance smoothing (m). min.points/max.points limit the IDW neighbours.
# smoothing.cells: odd moving-mean width in cells. Default 1 preserves the DEM.
#   Missing cells remain missing. Interpolated cells outside the point hull are
#   masked. A supplied DEM retains its original support, including internal NAs.
# neighbors: 4 or 8 neighbours for slope/aspect. TRI always uses 8 neighbours.
# position.radius: square TPI neighbourhood radius (m), rounded up to whole cells;
#   the focal cell is excluded. This is not the legacy CTFS quadrat convexity.
# flow.preparation: "breach" (default, with filling where breaching fails),
#   "fill" (raise closed depressions), or "none" for an already prepared DEM.
#   Preparation affects a separate routing copy, not the returned elevation.
# breach.distance.cells: maximum least-cost breaching search distance in cells.
# channel.area.m2: D8 contributing-area threshold defining the modelled channels.
#   Default 1000 is an example starting value, not a site-calibrated threshold.
#
# METHODS AND LIMITS
# terra is required. whitebox and its separately installed executable are needed
# only for curvature/hydrology; interp only for linear/Akima interpolation.
# All elevations and horizontal distances must be in metres. Aspect is relative
# to grid north; geographic exposure requires a north-aligned/georeferenced grid.
# Curvature uses the WhiteboxTools 2.4.0 definitions (5x5 fit, log=FALSE).
# TPI and curvature retain NA where their full neighbourhood is unavailable;
# flat-cell aspect and zero-slope TWI are undefined. Hydrological variables use
# only the supplied domain and cannot recover missing upstream terrain.
# TWI is a wetness proxy; drainage and HAND are not observed streams/water levels.
# For definitions, ecological references and parameter comparisons, see
# ../../tutorials/derive-topographical-variables.qmd.
#
# EXAMPLES
# topography <- derive_topography(xyz = elevations, resolution = 5)
# local.layers <- derive_topography(xyz = elevations, variables = c("slope", "tpi"))
# from.dem <- derive_topography(dem = my.dem, variables = c("slope", "tri"))

derive_topography <- function(xyz = NULL, dem = NULL, variables = "all",
                              resolution = 5, extent = NULL, crs = "local",
                              interpolation = "idw", search.radius = 60,
                              idw.power = 2, idw.smooth = 0,
                              min.points = 3, max.points = 12,
                              smoothing.cells = 1, neighbors = 8,
                              position.radius = 20, flow.preparation = "breach",
                              breach.distance.cells = 20, channel.area.m2 = 1000)
{
	##################################################
	# CHECK INPUTS AND SELECT CALCULATIONS
	##################################################

	if(!requireNamespace("terra", quietly = TRUE)) stop("Install the terra package.")
	if(is.null(xyz) == is.null(dem)) stop("Supply exactly one of xyz or dem.")
	layer.names <- c(elevation = "elevation_m", slope = "slope_deg", aspect = "aspect_deg",
	                 northness = "northness", eastness = "eastness", tpi = "tpi_m",
	                 tri = "tri_m", plan_curvature = "plan_curvature_per_m",
	                 profile_curvature = "profile_curvature_per_m",
	                 contributing_area = "contributing_area_m2", twi = "twi",
	                 drainage_distance = "drainage_distance_m",
	                 height_above_drainage = "height_above_drainage_m")
	layer.units <- c("m", "degrees", "degrees", "1", "1", "m", "m",
	                 "1/m", "1/m", "m2", "1", "m", "m")
	names(layer.units) <- names(layer.names)
	if(!is.character(variables) || !length(variables) || anyNA(variables))
		stop("variables must be 'all' or a character vector of variable keys.")
	if(identical(variables, "all")) variables <- names(layer.names)
	if(any(!variables %in% names(layer.names)))
		stop("Unknown variable(s): ", paste(setdiff(variables, names(layer.names)), collapse = ", "))
	variables <- unique(c("elevation", variables))

	need.aspect <- any(variables %in% c("aspect", "northness", "eastness"))
	need.slope <- need.aspect || any(variables %in% c("slope", "twi"))
	need.curvature <- any(variables %in% c("plan_curvature", "profile_curvature"))
	need.sca <- any(variables %in% c("contributing_area", "twi"))
	need.channels <- any(variables %in% c("drainage_distance", "height_above_drainage"))
	need.flow <- need.sca || need.channels
	need.whitebox <- need.curvature || need.flow

	if(!is.numeric(smoothing.cells) || is.complex(smoothing.cells) || length(smoothing.cells) != 1L ||
	   !is.finite(smoothing.cells) || smoothing.cells <= 0 || smoothing.cells %% 2 != 1)
		stop("smoothing.cells must be a positive odd integer; 1 means no smoothing.")
	if(need.slope && (!is.numeric(neighbors) || is.complex(neighbors) ||
	   length(neighbors) != 1L || is.na(neighbors) || !neighbors %in% c(4, 8)))
		stop("neighbors must be 4 or 8.")
	if("tpi" %in% variables && (!is.numeric(position.radius) || is.complex(position.radius) || length(position.radius) != 1L ||
	   !is.finite(position.radius) || position.radius <= 0))
		stop("position.radius must be one positive finite distance in metres.")
	if(need.whitebox)
	{
		if(!requireNamespace("whitebox", quietly = TRUE)) stop("Install the whitebox package.")
		executable <- whitebox::wbt_exe_path(shell_quote = FALSE)
		if(length(executable) != 1L || !file.exists(executable))
			stop("WhiteboxTools is not installed/configured. Run whitebox::install_whitebox() or whitebox::wbt_init().")
	}
	if(need.flow && (!is.character(flow.preparation) || length(flow.preparation) != 1L ||
	   is.na(flow.preparation) || !flow.preparation %in% c("breach", "fill", "none")))
		stop("flow.preparation must be 'breach', 'fill' or 'none'.")
	if(need.flow && flow.preparation == "breach" &&
	   (!is.numeric(breach.distance.cells) || is.complex(breach.distance.cells) || length(breach.distance.cells) != 1L ||
	   !is.finite(breach.distance.cells) || breach.distance.cells < 1 ||
	   breach.distance.cells != floor(breach.distance.cells)))
		stop("breach.distance.cells must be a positive integer.")
	if(need.channels && (!is.numeric(channel.area.m2) || is.complex(channel.area.m2) || length(channel.area.m2) != 1L ||
	   !is.finite(channel.area.m2) || channel.area.m2 <= 0))
		stop("channel.area.m2 must be one positive finite area in square metres.")


	##################################################
	# BUILD A DEM OR RETAIN THE SUPPLIED DEM
	##################################################

	if(!is.null(dem))
	{
		if(!inherits(dem, "SpatRaster") || terra::nlyr(dem) != 1L || !terra::hasValues(dem))
			stop("dem must be one SpatRaster layer containing elevations.")
		if(!missing(resolution) || !is.null(extent) || !missing(crs))
			stop("resolution, extent and crs apply to xyz only. Adjust a supplied DEM explicitly before this call.")
		surface <- dem
	}

	if(!is.null(xyz))
	{
		# Normalize a coordinate table without requiring any census-specific fields.
		if(!is.data.frame(xyz) && !is.matrix(xyz)) stop("xyz must be a data frame or matrix.")
		if(is.matrix(xyz) && is.null(colnames(xyz)) && ncol(xyz) == 3L)
			colnames(xyz) <- c("x", "y", "z")
		columns <- c("x", "y", "z")
		if(!all(columns %in% colnames(xyz)) || anyDuplicated(colnames(xyz)[colnames(xyz) %in% columns]))
			stop("xyz must contain one of each column: x, y and z.")
		xyz <- as.data.frame(xyz)[columns]
		if(!all(sapply(xyz, function(x) is.numeric(x) && !is.complex(x) && is.null(dim(x)))) ||
		   any(is.infinite(as.matrix(xyz)))) stop("xyz columns must be numeric, with finite values or NA.")
		complete <- complete.cases(xyz)
		if(any(!complete)) warning(sum(!complete), " incomplete xyz row(s) omitted.", call. = FALSE)
		xyz <- unique(xyz[complete, , drop = FALSE])
		if(anyDuplicated(xyz[c("x", "y")])) stop("Conflicting elevations at the same x,y coordinates.")
		if(nrow(xyz) < 3L || length(grDevices::chull(xyz$x, xyz$y)) < 3L)
			stop("xyz needs at least three distinct, non-collinear points.")
		xyz <- xyz[order(xyz$x, xyz$y), , drop = FALSE]

		# The output grid is explicit. Automatic bounds are rounded outwards.
		if(!is.numeric(resolution) || is.complex(resolution) || length(resolution) != 1L ||
		   !is.finite(resolution) || resolution <= 0) stop("resolution must be one positive finite cell size.")
		if(is.null(extent))
			extent <- c(floor(min(xyz$x) / resolution), ceiling(max(xyz$x) / resolution),
			            floor(min(xyz$y) / resolution), ceiling(max(xyz$y) / resolution)) * resolution
		if(!is.numeric(extent) || is.complex(extent) || length(extent) != 4L || any(!is.finite(extent)) ||
		   extent[1] >= extent[2] || extent[3] >= extent[4])
			stop("extent must be finite c(xmin, xmax, ymin, ymax), with positive width and height.")
		grid.dim <- c(extent[2] - extent[1], extent[4] - extent[3]) / resolution
		if(any(abs(grid.dim - round(grid.dim)) > 1e-7))
			stop("The extent width and height must be divisible by resolution.")
		surface <- terra::rast(xmin = extent[1], xmax = extent[2], ymin = extent[3], ymax = extent[4],
		                       resolution = resolution, crs = crs)
		linear.units <- terra::linearUnits(surface)
		if(!is.finite(linear.units) || abs(linear.units - 1) > 1e-8)
			stop("xyz must use a metre-based planar CRS, or crs = 'local'.")

		# Interpolation choices are explicit; no method is substituted after failure.
		if(!is.character(interpolation) || length(interpolation) != 1L || is.na(interpolation) ||
		   !interpolation %in% c("idw", "linear", "akima", "bilinear"))
			stop("interpolation must be 'idw', 'linear', 'akima' or 'bilinear'.")
		if(interpolation == "idw")
		{
			parameters <- list(search.radius = search.radius, idw.power = idw.power,
			                   idw.smooth = idw.smooth, min.points = min.points, max.points = max.points)
			for(parameter in names(parameters))
			{
				value <- parameters[[parameter]]
				if(!is.numeric(value) || is.complex(value) || length(value) != 1L || !is.finite(value) || value < 0)
					stop(parameter, " must be one finite, nonnegative number.")
			}
			if(search.radius <= 0 || idw.power <= 0 || min.points < 1 || max.points < min.points ||
			   min.points != floor(min.points) || max.points != floor(max.points))
				stop("IDW requires positive radius/power and integer point limits with 1 <= min.points <= max.points.")
			surface <- terra::interpIDW(surface, as.matrix(xyz), radius = search.radius,
			                             power = idw.power, smooth = idw.smooth,
			                             minPoints = min.points, maxPoints = max.points)
		}
		if(interpolation %in% c("linear", "akima"))
		{
			if(!requireNamespace("interp", quietly = TRUE)) stop("Install interp for linear/Akima interpolation.")
			if(nrow(xyz) < 4L) stop("interp requires at least four elevation points.")
			cell.xy <- terra::xyFromCell(surface, seq_len(terra::ncell(surface)))
			terra::values(surface) <- interp::interp(xyz$x, xyz$y, xyz$z,
				xo = cell.xy[, 1], yo = cell.xy[, 2], method = interpolation, output = "points")$z
		}
		if(interpolation == "bilinear")
		{
			x <- sort(unique(xyz$x))
			y <- sort(unique(xyz$y))
			if(nrow(xyz) != length(x) * length(y) ||
			   max(abs(diff(x) - diff(x)[1])) > 1e-7 * diff(x)[1] ||
			   max(abs(diff(y) - diff(y)[1])) > 1e-7 * diff(y)[1])
				stop("Bilinear interpolation requires a complete, regularly spaced xyz grid.")
			grid.dem <- terra::rast(xyz, type = "xyz", crs = crs)
			surface <- terra::resample(grid.dem, surface, method = "bilinear")
		}
		points <- terra::vect(xyz, geom = c("x", "y"), crs = terra::crs(surface))
		surface <- terra::mask(surface, terra::convHull(points), touches = FALSE)
	}

	# Every method below uses a square, planar raster with metre units.
	linear.units <- terra::linearUnits(surface)
	if(!is.finite(linear.units) || abs(linear.units - 1) > 1e-8)
		stop("Use a planar CRS with metre units (or crs = 'local' for local xyz coordinates).")
	cell.size <- terra::res(surface)
	if(abs(cell.size[1] - cell.size[2]) > 1e-8 * max(cell.size))
		stop("The DEM must have square cells; resample it explicitly before this call.")
	if(any(terra::global(is.infinite(surface), "sum", na.rm = TRUE)[, 1] > 0))
		stop("DEM elevations cannot contain Inf or -Inf.")
	if(terra::global(!is.na(surface), "sum", na.rm = TRUE)[1, 1] == 0)
		stop("The DEM has no usable elevations. Check the extent and interpolation settings.")

	# Smoothing is optional and never fills cells that were already missing.
	if(smoothing.cells > 1)
		surface <- terra::focal(surface, w = smoothing.cells, fun = "mean",
		                        na.rm = TRUE, na.policy = "omit")
	layers <- list(elevation = surface)


	##################################################
	# CALCULATE LOCAL TERRAIN VARIABLES
	##################################################

	if(need.slope)
		layers$slope <- terra::terrain(surface, v = "slope", unit = "degrees", neighbors = neighbors)
	if(need.aspect)
	{
		layers$aspect <- terra::terrain(surface, v = "aspect", unit = "degrees", neighbors = neighbors)
		layers$aspect[layers$slope == 0] <- NA
		if("northness" %in% variables) layers$northness <- cos(layers$aspect * pi / 180)
		if("eastness" %in% variables) layers$eastness <- sin(layers$aspect * pi / 180)
	}
	if("tpi" %in% variables)
	{
		radius.cells <- ceiling(position.radius / cell.size[1])
		weights <- matrix(1, 2 * radius.cells + 1, 2 * radius.cells + 1)
		weights[radius.cells + 1, radius.cells + 1] <- NA
		layers$tpi <- surface - terra::focal(surface, w = weights, fun = "mean", na.rm = FALSE)
	}
	if("tri" %in% variables) layers$tri <- terra::terrain(surface, v = "TRI")


	##################################################
	# CALCULATE CURVATURE AND FLOW-BASED VARIABLES
	##################################################

	if(need.whitebox)
	{
		# Isolate external-tool files. Returned layers are detached before cleanup.
		work.dir <- tempfile("fgeo2-topography-work-")
		dir.create(work.dir)
		on.exit(unlink(work.dir, recursive = TRUE, force = TRUE), add = TRUE)
		work.dir <- normalizePath(work.dir, winslash = "/")
		dem.file <- file.path(work.dir, "elevation.tif")
		terra::writeRaster(surface, dem.file, overwrite = TRUE)

		if(need.curvature)
		{
			# Keep curvature on the physical DEM, not the drainage-conditioned copy.
			support <- terra::focal(!is.na(surface), w = 5, fun = "sum", fillvalue = 0)
			if("plan_curvature" %in% variables)
			{
				file <- file.path(work.dir, "plan-curvature.tif")
				tool.output <- whitebox::wbt_plan_curvature(dem.file, file, log = FALSE, zfactor = 1, verbose_mode = FALSE)
				if(!file.exists(file)) stop("WhiteboxTools plan curvature failed. ", paste(tool.output, collapse = "\n"))
				layers$plan_curvature <- terra::rast(file)
				layers$plan_curvature[support < 25] <- NA
			}
			if("profile_curvature" %in% variables)
			{
				file <- file.path(work.dir, "profile-curvature.tif")
				tool.output <- whitebox::wbt_profile_curvature(dem.file, file, log = FALSE, zfactor = 1, verbose_mode = FALSE)
				if(!file.exists(file)) stop("WhiteboxTools profile curvature failed. ", paste(tool.output, collapse = "\n"))
				layers$profile_curvature <- terra::rast(file)
				layers$profile_curvature[support < 25] <- NA
			}
		}
		if(need.flow)
		{
			# Prepare a separate copy. Real depressions may require another choice.
			hydro.file <- dem.file
			if(flow.preparation != "none") hydro.file <- file.path(work.dir, "elevation-for-flow.tif")
			if(flow.preparation == "breach")
				tool.output <- whitebox::wbt_breach_depressions_least_cost(dem.file, hydro.file,
					dist = breach.distance.cells, fill = TRUE, verbose_mode = FALSE)
			if(flow.preparation == "fill")
				tool.output <- whitebox::wbt_fill_depressions(dem.file, hydro.file,
					fix_flats = TRUE, verbose_mode = FALSE)
			if(!file.exists(hydro.file)) stop("WhiteboxTools drainage preparation failed. ", paste(tool.output, collapse = "\n"))
			if(need.sca)
			{
				sca.file <- file.path(work.dir, "specific-contributing-area.tif")
				tool.output <- whitebox::wbt_d_inf_flow_accumulation(hydro.file, sca.file,
					out_type = "Specific Contributing Area", log = FALSE, verbose_mode = FALSE)
				if(!file.exists(sca.file)) stop("WhiteboxTools D-infinity accumulation failed. ", paste(tool.output, collapse = "\n"))
				if("contributing_area" %in% variables)
					layers$contributing_area <- terra::rast(sca.file) * cell.size[1]
				if("twi" %in% variables)
				{
					slope.file <- file.path(work.dir, "slope.tif")
					file <- file.path(work.dir, "twi.tif")
					terra::writeRaster(layers$slope, slope.file, overwrite = TRUE)
					tool.output <- whitebox::wbt_wetness_index(sca.file, slope.file, file, verbose_mode = FALSE)
					if(!file.exists(file)) stop("WhiteboxTools wetness index failed. ", paste(tool.output, collapse = "\n"))
					layers$twi <- terra::rast(file)
				}
			}
			if(need.channels)
			{
				# D8 supplies a single-path network for both distance and HAND.
				area.file <- file.path(work.dir, "d8-area.tif")
				streams.file <- file.path(work.dir, "drainage.tif")
				tool.output <- whitebox::wbt_d8_flow_accumulation(hydro.file, area.file,
					out_type = "catchment area", log = FALSE, verbose_mode = FALSE)
				if(!file.exists(area.file)) stop("WhiteboxTools D8 accumulation failed. ", paste(tool.output, collapse = "\n"))
				tool.output <- whitebox::wbt_extract_streams(area.file, streams.file,
					threshold = channel.area.m2, zero_background = TRUE, verbose_mode = FALSE)
				if(!file.exists(streams.file)) stop("WhiteboxTools channel extraction failed. ", paste(tool.output, collapse = "\n"))
				stream.maximum <- terra::global(terra::rast(streams.file), "max", na.rm = TRUE)[1, 1]
				if(!is.finite(stream.maximum) || stream.maximum <= 0)
					stop("No channels at channel.area.m2. Review the threshold or use a larger drainage domain.")
				if("drainage_distance" %in% variables)
				{
					file <- file.path(work.dir, "distance.tif")
					tool.output <- whitebox::wbt_downslope_distance_to_stream(hydro.file, streams.file, file,
						dinf = FALSE, verbose_mode = FALSE)
					if(!file.exists(file)) stop("WhiteboxTools drainage distance failed. ", paste(tool.output, collapse = "\n"))
					layers$drainage_distance <- terra::rast(file)
				}
				if("height_above_drainage" %in% variables)
				{
					file <- file.path(work.dir, "hand.tif")
					tool.output <- whitebox::wbt_elevation_above_stream(hydro.file, streams.file, file, verbose_mode = FALSE)
					if(!file.exists(file)) stop("WhiteboxTools height above drainage failed. ", paste(tool.output, collapse = "\n"))
					layers$height_above_drainage <- terra::rast(file)
				}
			}
		}
	}


	##################################################
	# RETURN STANDARD RASTER LAYERS WITH NAMES AND UNITS
	##################################################

	out <- terra::rast(layers[variables])
	names(out) <- unname(layer.names[variables])
	terra::units(out) <- unname(layer.units[variables])
	if(need.whitebox)
	{
		# This file lives outside work.dir and remains valid for the R session.
		result.file <- tempfile("fgeo2-topography-result-", fileext = ".tif")
		out <- terra::writeRaster(out, result.file, datatype = "FLT8S", overwrite = TRUE)
	}
	out
}
