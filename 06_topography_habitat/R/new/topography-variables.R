# Derive topographical raster layers from point elevations or an existing DEM.
#
# INPUT AND OUTPUT
# Supply exactly one input:
# - xyz: a data frame or matrix containing x, y and z in metres. An unnamed
#   three-column matrix is interpreted in that order. Extra columns are ignored.
#   Incomplete rows are omitted with a warning. Identical records are collapsed;
#   conflicting elevations at the same coordinates must be resolved by the user.
#   Default linear interpolation requires at least four points (interp package).
# - dem: one SpatRaster elevation layer, with square cells and a metre-based
#   planar CRS. Its grid is retained. Interpolation is skipped.
#   Interpolation and previous smoothing are already in the cell elevations;
#   derivative settings must be supplied again. Use smoothing.cells = 1 to avoid
#   smoothing a previously processed DEM twice.
#
# Returns a standard SpatRaster: elevation first, then the requested variables.
# Names and units are retained. No plots or permanent files are produced.
# With fill.missing = TRUE, attr(result, "filled.cells") is a same-grid raster
# marking substitutions in this call (1 = filled, 0 = unchanged). Save it separately;
# ordinary raster transformations and GeoTIFF export may drop R attributes.
# WhiteboxTools uses a private temporary directory. Requested results are copied
# to a separate R temporary GeoTIFF before that directory is removed. Save the
# returned raster with terra::writeRaster() to retain it across R sessions.
#
# VARIABLES
# variables = "all" returns these thirteen layers. Alternatively supply a vector
# of the keys on the left; elevation is always included, and dependencies are
# calculated only as needed. Unrequested intermediate layers are not returned.
# The CTFS quantities have raster counterparts in elevation, slope, TPI (relative
# position, analogous to quadrat convexity) and aspect (downhill direction).
# Their spatial support and algorithms differ; curvature is distinct from TPI.
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
# interpolation: "linear" (default, interp package), "idw" (terra), "akima"
#   (interp), or "bilinear" for a complete regular xyz grid. No silent fallback.
#   Linear interpolation preserves planes and does not create new extrema.
#   IDW can imprint the sampling grid on slope/curvature; Akima can overshoot.
# search.radius, idw.power, idw.smooth: IDW search distance (m), distance power,
#   and distance smoothing (m). min.points/max.points limit the IDW neighbours.
# smoothing.cells: odd moving-mean width in cells. Default 1 preserves the DEM.
#   Smoothing preserves missing cells; nearest-value completion happens only
#   after all derivatives and flow calculations, when fill.missing = TRUE.
# neighbors: 4 or 8 neighbours for slope/aspect. TRI always uses 8 neighbours.
# position.radius: square TPI neighbourhood radius (m), rounded up to whole cells;
#   the focal cell is excluded. This is not the legacy CTFS quadrat convexity.
# flow.preparation: "breach" (default, with filling where breaching fails),
#   "fill" (raise closed depressions), or "none" for an already prepared DEM.
#   Preparation affects a separate routing copy, not the returned elevation.
# breach.distance.cells: maximum least-cost breaching search distance in cells.
# channel.area.m2: D8 contributing-area threshold defining the modelled channels.
#   Default 1000 is an example starting value, not a site-calibrated threshold.
# fill.missing: TRUE (default) replaces non-finite output cells with the nearest
#   finite estimate of that variable, anywhere inside the existing raster extent.
#   This includes gaps, borders and cells outside the point hull. Valid estimates
#   are unchanged. Aspect is filled once; northness/eastness use that direction.
#   FALSE retains the unfilled calculations. A layer with no finite donor raises
#   an error with TRUE; this can occur for aspect/TWI on entirely flat terrain.
#
# METHODS AND LIMITS
# terra is required. interp is needed for default linear/Akima point input;
# DEM, IDW and bilinear input do not need it. whitebox and its separately
# installed executable are needed only for curvature/hydrology.
# All elevations and horizontal distances must be in metres. Aspect is relative
# to grid north; geographic exposure requires a north-aligned/georeferenced grid.
# Curvature uses the WhiteboxTools 2.4.0 definitions (5x5 fit, log=FALSE).
# Before completion, TPI/curvature require their full neighbourhood; flat-cell
# aspect and zero-slope TWI are undefined. Filled values are spatial proxies,
# including for drainage distance and HAND; no new flow path is calculated.
# Completion does not expand the grid or recover missing upstream terrain.
# A missing channel network still requires a revised domain or channel threshold.
# TWI is a wetness proxy; drainage and HAND are not observed streams/water levels.
# Raster cells and D8 flow paths remain discrete. Smoother-looking maps do not
# recover terrain detail absent from the survey. Do not smooth the derived
# drainage layers just to hide their grid. Supply surrounding terrain for edges.
# For definitions, ecological references and parameter comparisons, see
# ../../tutorials/derive-topographical-variables.qmd.
#
# EXAMPLES
# topography <- derive_topography(xyz = elevations, resolution = 5)
# local.layers <- derive_topography(xyz = elevations, variables = c("slope", "tpi"))
# from.dem <- derive_topography(dem = my.dem, variables = c("slope", "tri"))
# unfilled <- derive_topography(xyz = elevations, fill.missing = FALSE)

derive_topography <- function(xyz = NULL, dem = NULL, variables = "all",
                              resolution = 5, extent = NULL, crs = "local",
                              interpolation = "linear", search.radius = 60,
                              idw.power = 2, idw.smooth = 0,
                              min.points = 3, max.points = 12,
                              smoothing.cells = 1, neighbors = 8,
                              position.radius = 20, flow.preparation = "breach",
                              breach.distance.cells = 20, channel.area.m2 = 1000,
                              fill.missing = TRUE)
{
	##################################################
	# CHECK INPUT TYPES AND REQUESTED VARIABLES
	##################################################

	if(!requireNamespace("terra", quietly = TRUE)) stop("Install the terra package.")

	if(is.null(xyz) == is.null(dem)) stop("Supply exactly one of xyz or dem.")

	if(!is.logical(fill.missing) || length(fill.missing) != 1L || is.na(fill.missing))
		stop("fill.missing must be TRUE or FALSE.")

	from.xyz <- !is.null(xyz)
	from.dem <- !is.null(dem)

	# Keep the output names and units together, independent of input type.
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

	# Some outputs require intermediate calculations that are not returned.
	need.northness <- "northness" %in% variables
	need.eastness <- "eastness" %in% variables
	need.aspect <- any(variables %in% c("aspect", "northness", "eastness"))
	need.slope <- need.aspect || any(variables %in% c("slope", "twi"))
	need.tpi <- "tpi" %in% variables
	need.tri <- "tri" %in% variables
	need.plan <- "plan_curvature" %in% variables
	need.profile <- "profile_curvature" %in% variables
	need.curvature <- need.plan || need.profile
	need.area <- "contributing_area" %in% variables
	need.twi <- "twi" %in% variables
	need.sca <- need.area || need.twi
	need.distance <- "drainage_distance" %in% variables
	need.hand <- "height_above_drainage" %in% variables
	need.channels <- need.distance || need.hand
	need.flow <- need.sca || need.channels
	need.whitebox <- need.curvature || need.flow


	##################################################
	# CHECK CALCULATION PARAMETERS AND DEPENDENCIES
	##################################################

	# Validate only the parameters needed by the selected calculations.
	if(!is.numeric(smoothing.cells) || is.complex(smoothing.cells) || length(smoothing.cells) != 1L ||
	   !is.finite(smoothing.cells) || smoothing.cells <= 0 || smoothing.cells %% 2 != 1)
		stop("smoothing.cells must be a positive odd integer; 1 means no smoothing.")

	if(need.slope && (!is.numeric(neighbors) || is.complex(neighbors) ||
	   length(neighbors) != 1L || is.na(neighbors) || !neighbors %in% c(4, 8)))
		stop("neighbors must be 4 or 8.")

	if(need.tpi && (!is.numeric(position.radius) || is.complex(position.radius) || length(position.radius) != 1L ||
	   !is.finite(position.radius) || position.radius <= 0))
		stop("position.radius must be one positive finite distance in metres.")

	# Local variables do not require WhiteboxTools or its executable.
	if(need.whitebox && !requireNamespace("whitebox", quietly = TRUE))
		stop("Install the whitebox package.")

	if(need.whitebox) executable <- whitebox::wbt_exe_path(shell_quote = FALSE)

	if(need.whitebox && (length(executable) != 1L || !file.exists(executable)))
		stop("WhiteboxTools is not installed/configured. Run whitebox::install_whitebox() or whitebox::wbt_init().")

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

	# These cases are separate operations, not nested alternatives.
	prepare.flow <- need.flow && flow.preparation != "none"
	breach.flow <- need.flow && flow.preparation == "breach"
	fill.flow <- need.flow && flow.preparation == "fill"


	##################################################
	# RETAIN A SUPPLIED DEM OR PREPARE ELEVATION POINTS
	##################################################

	# A supplied DEM keeps its grid. Changes to that grid must be explicit.
	if(from.dem && (!inherits(dem, "SpatRaster") || terra::nlyr(dem) != 1L || !terra::hasValues(dem)))
		stop("dem must be one SpatRaster layer containing elevations.")

	if(from.dem && (!missing(resolution) || !is.null(extent) || !missing(crs)))
		stop("resolution, extent and crs apply to xyz only. Adjust a supplied DEM explicitly before this call.")

	if(from.dem) surface <- dem

	# Normalize point coordinates without requiring any census-specific fields.
	if(from.xyz && !is.data.frame(xyz) && !is.matrix(xyz))
		stop("xyz must be a data frame or matrix.")

	if(from.xyz && is.matrix(xyz) && is.null(colnames(xyz)) && ncol(xyz) == 3L)
		colnames(xyz) <- c("x", "y", "z")

	columns <- c("x", "y", "z")

	if(from.xyz && (!all(columns %in% colnames(xyz)) ||
	   anyDuplicated(colnames(xyz)[colnames(xyz) %in% columns])))
		stop("xyz must contain one of each column: x, y and z.")

	if(from.xyz) xyz <- as.data.frame(xyz)[columns]

	if(from.xyz && (!all(sapply(xyz, function(x) is.numeric(x) && !is.complex(x) && is.null(dim(x)))) ||
	   any(is.infinite(as.matrix(xyz)))))
		stop("xyz columns must be numeric, with finite values or NA.")

	# Missing rows and identical repeats are handled before interpolation.
	if(from.xyz) complete <- complete.cases(xyz)

	if(from.xyz && any(!complete))
		warning(sum(!complete), " incomplete xyz row(s) omitted.", call. = FALSE)

	if(from.xyz) xyz <- unique(xyz[complete, , drop = FALSE])

	if(from.xyz && anyDuplicated(xyz[c("x", "y")]))
		stop("Conflicting elevations at the same x,y coordinates.")

	if(from.xyz && (nrow(xyz) < 3L || length(grDevices::chull(xyz$x, xyz$y)) < 3L))
		stop("xyz needs at least three distinct, non-collinear points.")

	if(from.xyz) xyz <- xyz[order(xyz$x, xyz$y), , drop = FALSE]


	##################################################
	# BUILD THE OUTPUT GRID FOR POINT INPUTS
	##################################################

	# Automatic bounds are rounded outwards; an explicit extent is retained.
	if(from.xyz && (!is.numeric(resolution) || is.complex(resolution) || length(resolution) != 1L ||
	   !is.finite(resolution) || resolution <= 0))
		stop("resolution must be one positive finite cell size.")

	if(from.xyz && is.null(extent))
		extent <- c(floor(min(xyz$x) / resolution), ceiling(max(xyz$x) / resolution),
		            floor(min(xyz$y) / resolution), ceiling(max(xyz$y) / resolution)) * resolution

	if(from.xyz && (!is.numeric(extent) || is.complex(extent) || length(extent) != 4L || any(!is.finite(extent)) ||
	   extent[1] >= extent[2] || extent[3] >= extent[4]))
		stop("extent must be finite c(xmin, xmax, ymin, ymax), with positive width and height.")

	if(from.xyz)
	  grid.dim <- c(extent[2] - extent[1], extent[4] - extent[3]) / resolution

	if(from.xyz && any(abs(grid.dim - round(grid.dim)) > 1e-7))
		stop("The extent width and height must be divisible by resolution.")

	if(from.xyz)
	{
		surface <- terra::rast(xmin = extent[1], xmax = extent[2], ymin = extent[3], ymax = extent[4],
		                       resolution = resolution, crs = crs)
		linear.units <- terra::linearUnits(surface)
	}

	if(from.xyz && (!is.finite(linear.units) || abs(linear.units - 1) > 1e-8))
		stop("xyz must use a metre-based planar CRS, or crs = 'local'.")


	##################################################
	# INTERPOLATE ELEVATION POINTS ONTO THE GRID
	##################################################

	# Only point input uses interpolation. No method is substituted after failure.
	if(from.xyz && (!is.character(interpolation) || length(interpolation) != 1L || is.na(interpolation) ||
	   !interpolation %in% c("idw", "linear", "akima", "bilinear")))
		stop("interpolation must be 'idw', 'linear', 'akima' or 'bilinear'.")

	use.idw <- from.xyz && interpolation == "idw"
	use.interp <- from.xyz && interpolation %in% c("linear", "akima")
	use.bilinear <- from.xyz && interpolation == "bilinear"

	# Inverse-distance weighting uses an explicit search and point count.
	parameters <- list()

	if(use.idw)
		parameters <- list(search.radius = search.radius, idw.power = idw.power,
		                   idw.smooth = idw.smooth, min.points = min.points, max.points = max.points)

	for(parameter in names(parameters))
	{
		value <- parameters[[parameter]]

		if(!is.numeric(value) || is.complex(value) || length(value) != 1L || !is.finite(value) || value < 0)
			stop(parameter, " must be one finite, nonnegative number.")
	}

	if(use.idw && (search.radius <= 0 || idw.power <= 0 || min.points < 1 || max.points < min.points ||
	   min.points != floor(min.points) || max.points != floor(max.points)))
		stop("IDW requires positive radius/power and integer point limits with 1 <= min.points <= max.points.")

	if(use.idw)
		surface <- terra::interpIDW(surface, as.matrix(xyz), radius = search.radius,
		                             power = idw.power, smooth = idw.smooth,
		                             minPoints = min.points, maxPoints = max.points)

	# Linear and Akima interpolation evaluate the same raster-cell centres.
	if(use.interp && !requireNamespace("interp", quietly = TRUE))
		stop("Install interp for linear/Akima interpolation.")

	if(use.interp && nrow(xyz) < 4L)
	  stop("interp requires at least four elevation points.")

	if(use.interp)
	{
		cell.xy <- terra::xyFromCell(surface, seq_len(terra::ncell(surface)))
		terra::values(surface) <- interp::interp(xyz$x, xyz$y, xyz$z,
			xo = cell.xy[, 1], yo = cell.xy[, 2], method = interpolation, output = "points")$z
	}

	# Bilinear interpolation is restricted to a complete regular input grid.
	if(use.bilinear)
	{
		x <- sort(unique(xyz$x))
		y <- sort(unique(xyz$y))
	}

	if(use.bilinear && (nrow(xyz) != length(x) * length(y) ||
	   max(abs(diff(x) - diff(x)[1])) > 1e-7 * diff(x)[1] ||
	   max(abs(diff(y) - diff(y)[1])) > 1e-7 * diff(y)[1]))
		stop("Bilinear interpolation requires a complete, regularly spaced xyz grid.")

	if(use.bilinear)
	{
		grid.dem <- terra::rast(xyz, type = "xyz", crs = crs)
		surface <- terra::resample(grid.dem, surface, method = "bilinear")
	}

	# Do not extrapolate beyond the horizontal support of the measured points.
	if(from.xyz)
	{
		points <- terra::vect(xyz, geom = c("x", "y"), crs = terra::crs(surface))
		surface <- terra::mask(surface, terra::convHull(points), touches = FALSE)
	}


	##################################################
	# CHECK AND OPTIONALLY SMOOTH THE DEM
	##################################################

	# Both input paths now use square raster cells and planar metre coordinates.
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

	# A moving mean never fills cells that were already missing.
	if(smoothing.cells > 1)
		surface <- terra::focal(surface, w = smoothing.cells, fun = "mean",
		                        na.rm = TRUE, na.policy = "omit")

	layers <- list(elevation = surface)


	##################################################
	# CALCULATE LOCAL TERRAIN VARIABLES
	##################################################

	# Inclination and downhill orientation share the chosen 4/8-neighbour method.
	if(need.slope)
		layers$slope <- terra::terrain(surface, v = "slope", unit = "degrees", neighbors = neighbors)

	if(need.aspect)
	{
		layers$aspect <- terra::terrain(surface, v = "aspect", unit = "degrees", neighbors = neighbors)
		layers$aspect[layers$slope == 0] <- NA
	}

	if(need.northness)
	  layers$northness <- cos(layers$aspect * pi / 180)

	if(need.eastness)
	  layers$eastness <- sin(layers$aspect * pi / 180)

	# TPI compares the focal cell with its stated neighbourhood, excluding itself.
	if(need.tpi)
	{
		radius.cells <- ceiling(position.radius / cell.size[1])
		weights <- matrix(1, 2 * radius.cells + 1, 2 * radius.cells + 1)
		weights[radius.cells + 1, radius.cells + 1] <- NA
		layers$tpi <- surface - terra::focal(surface, w = weights, fun = "mean", na.rm = FALSE)
	}

	# TRI always uses the eight neighbours, independently of the slope setting.
	if(need.tri)
	  layers$tri <- terra::terrain(surface, v = "TRI")


	##################################################
	# PREPARE TEMPORARY FILES FOR WHITEBOXTOOLS
	##################################################

	# Isolate external-tool files. Returned layers are detached before cleanup.
	if(need.whitebox)
	{
		work.dir <- tempfile("fgeo2-topography-work-")
		dir.create(work.dir)
		on.exit(unlink(work.dir, recursive = TRUE, force = TRUE), add = TRUE)
		work.dir <- normalizePath(work.dir, winslash = "/")
		dem.file <- file.path(work.dir, "elevation.tif")

		# Retain elevation precision before second derivatives and flow routing.
		terra::writeRaster(surface, dem.file, datatype = "FLT8S", overwrite = TRUE)
	}


	##################################################
	# CALCULATE PLAN AND PROFILE CURVATURE
	##################################################

	# Curvature uses the physical DEM and requires a full 5 by 5 neighbourhood.
	if(need.curvature)
		support <- terra::focal(!is.na(surface), w = 5, fun = "sum", fillvalue = 0)

	if(need.plan)
	{
		plan.file <- file.path(work.dir, "plan-curvature.tif")
		tool.output <- whitebox::wbt_plan_curvature(dem.file, plan.file,
			log = FALSE, zfactor = 1, verbose_mode = FALSE)
	}

	if(need.plan && !file.exists(plan.file))
		stop("WhiteboxTools plan curvature failed. ", paste(tool.output, collapse = "\n"))

	if(need.plan)
	{
		layers$plan_curvature <- terra::rast(plan.file)
		layers$plan_curvature[support < 25] <- NA
	}

	if(need.profile)
	{
		profile.file <- file.path(work.dir, "profile-curvature.tif")
		tool.output <- whitebox::wbt_profile_curvature(dem.file, profile.file,
			log = FALSE, zfactor = 1, verbose_mode = FALSE)
	}

	if(need.profile && !file.exists(profile.file))
		stop("WhiteboxTools profile curvature failed. ", paste(tool.output, collapse = "\n"))

	if(need.profile)
	{
		layers$profile_curvature <- terra::rast(profile.file)
		layers$profile_curvature[support < 25] <- NA
	}


	##################################################
	# PREPARE A SEPARATE SURFACE FOR FLOW ROUTING
	##################################################

	# Breaching/filling change only the routing copy, never the local terrain DEM.
	if(need.flow)
	  hydro.file <- dem.file

	if(prepare.flow)
	  hydro.file <- file.path(work.dir, "elevation-for-flow.tif")

	if(breach.flow)
		tool.output <- whitebox::wbt_breach_depressions_least_cost(dem.file, hydro.file,
			dist = breach.distance.cells, fill = TRUE, verbose_mode = FALSE)

	if(fill.flow)
		tool.output <- whitebox::wbt_fill_depressions(dem.file, hydro.file,
			fix_flats = TRUE, verbose_mode = FALSE)

	if(need.flow && !file.exists(hydro.file))
		stop("WhiteboxTools drainage preparation failed. ", paste(tool.output, collapse = "\n"))


	##################################################
	# CALCULATE CONTRIBUTING AREA AND WETNESS
	##################################################

	# D-infinity may distribute flow between neighbouring directions.
	if(need.sca)
	{
		sca.file <- file.path(work.dir, "specific-contributing-area.tif")
		tool.output <- whitebox::wbt_d_inf_flow_accumulation(hydro.file, sca.file,
			out_type = "Specific Contributing Area", log = FALSE, verbose_mode = FALSE)
	}

	if(need.sca && !file.exists(sca.file))
		stop("WhiteboxTools D-infinity accumulation failed. ", paste(tool.output, collapse = "\n"))

	if(need.area)
	  layers$contributing_area <- terra::rast(sca.file) * cell.size[1]

	# TWI combines unlogged specific contributing area with slope in degrees.
	if(need.twi)
	{
		slope.file <- file.path(work.dir, "slope.tif")
		twi.file <- file.path(work.dir, "twi.tif")
		terra::writeRaster(layers$slope, slope.file, datatype = "FLT8S", overwrite = TRUE)
		tool.output <- whitebox::wbt_wetness_index(sca.file, slope.file, twi.file, verbose_mode = FALSE)
	}

	if(need.twi && !file.exists(twi.file))
		stop("WhiteboxTools wetness index failed. ", paste(tool.output, collapse = "\n"))

	if(need.twi) layers$twi <- terra::rast(twi.file)


	##################################################
	# DEFINE CHANNELS AND CALCULATE DRAINAGE POSITION
	##################################################

	# D8 supplies one path for both distance to drainage and height above drainage.
	if(need.channels)
	{
		area.file <- file.path(work.dir, "d8-area.tif")
		streams.file <- file.path(work.dir, "drainage.tif")
		tool.output <- whitebox::wbt_d8_flow_accumulation(hydro.file, area.file,
			out_type = "catchment area", log = FALSE, verbose_mode = FALSE)
	}

	if(need.channels && !file.exists(area.file))
		stop("WhiteboxTools D8 accumulation failed. ", paste(tool.output, collapse = "\n"))

	if(need.channels)
		tool.output <- whitebox::wbt_extract_streams(area.file, streams.file,
			threshold = channel.area.m2, zero_background = TRUE, verbose_mode = FALSE)

	if(need.channels && !file.exists(streams.file))
		stop("WhiteboxTools channel extraction failed. ", paste(tool.output, collapse = "\n"))

	if(need.channels)
	  stream.maximum <- terra::global(terra::rast(streams.file), "max", na.rm = TRUE)[1, 1]

	if(need.channels && (!is.finite(stream.maximum) || stream.maximum <= 0))
		stop("No channels at channel.area.m2. Review the threshold or use a larger drainage domain.")

	if(need.distance)
	{
		distance.file <- file.path(work.dir, "distance.tif")
		tool.output <- whitebox::wbt_downslope_distance_to_stream(hydro.file, streams.file, distance.file,
			dinf = FALSE, verbose_mode = FALSE)
	}

	if(need.distance && !file.exists(distance.file))
		stop("WhiteboxTools drainage distance failed. ", paste(tool.output, collapse = "\n"))

	if(need.distance)
	  layers$drainage_distance <- terra::rast(distance.file)

	if(need.hand)
	{
		hand.file <- file.path(work.dir, "hand.tif")
		tool.output <- whitebox::wbt_elevation_above_stream(hydro.file, streams.file, hand.file, verbose_mode = FALSE)
	}

	if(need.hand && !file.exists(hand.file))
		stop("WhiteboxTools height above drainage failed. ", paste(tool.output, collapse = "\n"))

	if(need.hand)
	  layers$height_above_drainage <- terra::rast(hand.file)


	##################################################
	# OPTIONALLY COMPLETE THE CALCULATED LAYERS
	##################################################

	filled.cells <- NULL
	if(fill.missing)
	{
		# Record the original gaps before replacing anything. Completion is last:
		# substituted elevation/slope must not feed back into terrain or routing.
		filled.cells <- terra::rast(lapply(layers[variables], function(x) !is.finite(x)))
		completion.variables <- setdiff(unique(c(if(need.aspect) "aspect", variables)),
		                                c("northness", "eastness"))
		# A raster diagonal includes every possible donor within this extent.
		fill.radius <- sqrt(sum((terra::res(surface) * c(ncol(surface), nrow(surface)))^2))

		for(variable in completion.variables)
		{
			layer <- layers[[variable]]
			missing <- !is.finite(layer)
			if(terra::global(missing, "sum", na.rm = TRUE)[1, 1] == 0) next
			if(terra::global(is.finite(layer), "sum", na.rm = TRUE)[1, 1] == 0)
				stop("Cannot fill ", layer.names[[variable]], ": no finite estimates. ",
				     "Use fill.missing = FALSE or revise the inputs/calculation settings.")

			# Borrow one nearest value, without averaging directions or smoothing
			# valid estimates. Non-finite values cannot serve as donors.
			layer[missing] <- NA
			donors <- terra::as.points(layer, na.rm = TRUE)
			nearest <- terra::interpNear(layer, donors, field = names(layer)[1],
			                             radius = fill.radius, interpolate = FALSE)
			layers[[variable]] <- terra::cover(layer, nearest)
		}

		# Complete the direction once, even if only a component was requested.
		if(need.northness) layers$northness <- cos(layers$aspect * pi / 180)
		if(need.eastness) layers$eastness <- sin(layers$aspect * pi / 180)
	}


	##################################################
	# RETURN STANDARD RASTER LAYERS WITH NAMES AND UNITS
	##################################################

	out <- terra::rast(layers[variables])
	names(out) <- unname(layer.names[variables])
	terra::units(out) <- unname(layer.units[variables])

	# This file lives outside work.dir and remains valid for the R session.
	if(need.whitebox)
	{
		result.file <- tempfile("fgeo2-topography-result-", fileext = ".tif")
		out <- terra::writeRaster(out, result.file, datatype = "FLT8S", overwrite = TRUE)
	}

	if(fill.missing)
	{
		names(filled.cells) <- names(out)
		terra::units(filled.cells) <- rep("1", terra::nlyr(out))
		if(need.whitebox)
		{
			# Keep provenance readable after the private Whitebox folder is removed.
			mask.file <- tempfile("fgeo2-topography-filled-", fileext = ".tif")
			filled.cells <- terra::writeRaster(filled.cells, mask.file, datatype = "INT1U",
			                                  overwrite = TRUE)
		}
		attr(out, "filled.cells") <- filled.cells
	}

	out
}
