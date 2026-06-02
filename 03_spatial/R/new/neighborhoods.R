# Code to study neighborhoods
# Gabriel Arellano
# gabriel.arellano.torres@gmail.com

get_distances_up_to_r <- function(query = NULL, reference = NULL, radius) {
  
  # Fast calculation of distances up to a given radius.
  # Note: if this function breaks, then you will have to find an
  # alternative that will calculate distances fast and return
  # a sparse distance matrix of the form {i, j, distance between i and j}.
  
  # Parameters:
  # query: the set of focal points, around which we are going to look
  # reference: the set of potential neighbors
  # radius: the size of the local neighborhood around the focal point
  
  # Checks and warnings:
  if(is.null(query)) stop("You must provide 'query'")
  if(is.null(reference)) stop("You must provide 'reference', even if reference = query")
  if(ncol(query) != 2) stop("'query' must have two columns with numerical coordinates")
  if(ncol(reference) != 2) stop("'reference' must have two columns with numerical coordinates; possibly reference = query")
  
  query <- as.data.frame(query)
  reference <- as.data.frame(reference)
  if(!is.numeric(query[,1]) | !is.numeric(query[,2]) | !is.numeric(reference[,1]) | !is.numeric(reference[,2]))
    stop("'query' and 'reference' must have numerical columns with coordinates")
  
  # To handle NAs, we have to record the original
  # indexes, then remove NAs, then come back
  # to original indexes before returning output.
  query$original_i <- 1:nrow(query)
  reference$original_j <- 1:nrow(reference)
  query <- na.omit(query)
  reference <- na.omit(reference)
  original_i <- query$original_i
  original_j <- reference$original_j
  
  # Reshape input, keep coordinates only
  query <- as.matrix(query[,c(1, 2)])
  reference <- as.matrix(reference[,c(1, 2)])
  storage.mode(query) <- "double"
  storage.mode(reference) <- "double"
  
  # Calculate distances
  nn <- dbscan::frNN(
    x = reference,
    query = query,
    eps = radius,
    sort = FALSE,
    approx = 0 # precision parameter, rarely needed to update
  )
  
  # Reorganize output
  k <- lengths(nn$id)
  d <- data.table::data.table(
    query_i = rep(seq_along(k), k),
    reference_j = unlist(nn$id, use.names = FALSE),
    dist = unlist(nn$dist, use.names = FALSE)
  )
  
  d <- d[d$dist <= radius,]
  
  # Recover the original indexes:
  d$query_i <- original_i[d$query_i]
  d$reference_j <- original_j[d$reference_j]
  d
}

# small example
if(FALSE)
{
  query <- data.frame(x = runif(3), y = runif(3))
  reference <- data.frame(x = runif(10), y = runif(10))
  query[2,] <- c(NA, NA)
  reference[2,] <- c(NA, NA)

  get_distances_up_to_r(query, reference, radius = 100)
}



get_boundary <- function(boundary = NULL, window = NULL, xy = NULL, q = 5)
{
  # Minimal for spatial objects
  if (!is.null(boundary)) {
    if (inherits(boundary, "sfg")) {
      boundary <- sf::st_sfc(boundary)
    } else {
      boundary <- sf::st_geometry(boundary)
    }
  }
  
  # From "window" to proper sf object
  if (is.null(boundary) && !is.null(window)) {
    xmin <- window[1]
    ymin <- window[2]
    xmax <- window[3]
    ymax <- window[4]
    
    boundary <- sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(xmin, ymin),
        c(xmax, ymin),
        c(xmax, ymax),
        c(xmin, ymax),
        c(xmin, ymin)
      )))
    )
  }
  
  # If just the points are provided, guess the boundary
  if (is.null(boundary) && is.null(window)) {
    
    xy <- as.matrix(xy)
    stopifnot(ncol(xy) == 2L, nrow(xy) >= 2L)
    storage.mode(xy) <- "double"
    
    knn <- FNN::get.knn(xy, k = 1)
    buffer <- q * median(knn$nn.dist[, 1], na.rm = TRUE)
    delta <- buffer / 10
    
    cell_x <- floor(xy[, 1] / delta)
    cell_y <- floor(xy[, 2] / delta)
    cell_id <- paste(cell_x, cell_y, sep = "_")
    xy <- xy[!duplicated(cell_id), , drop = FALSE]
    
    pts <- sf::st_as_sf(
      as.data.frame(xy),
      coords = 1:2
    )
    
    boundary <- sf::st_buffer(pts, dist = buffer + delta)
    boundary <- sf::st_union(boundary)
    boundary <- sf::st_simplify(boundary, dTolerance = delta, preserveTopology = TRUE)
  }
  
  if (is.null(boundary)) {
    stop("Provide one of: boundary, window, or xy.")
  }
  
  # Final polygon cleanup
  boundary <- sf::st_geometry(boundary)
  boundary <- sf::st_make_valid(boundary)
  
  gtype <- sf::st_geometry_type(boundary)
  if (any(gtype == "GEOMETRYCOLLECTION")) 
    boundary <- sf::st_collection_extract(boundary, "POLYGON")
  
  if (length(boundary) == 0 || all(sf::st_is_empty(boundary))) 
    stop("Boundary does not contain polygon geometry.")
  
  boundary <- sf::st_union(boundary)
  boundary <- sf::st_make_valid(boundary)
  boundary
}


# Clasify points based on a given boundary.
# The boundary can be provided as a sfc_POLYGON of arbitrary shape
# or a rectangular window (xmin, ymin, xmax, ymax). If it is not
# provided it will be assumed rectangular around the points provided.

core_vs_edge <- function(xy,
                         boundary = NULL,
                         window = NULL,
                         radius) {
  
  # Parameters:
  # xy: a matrix or data.frame with x and y coordinates in 1st and 2nd column respectively
  # boundary: boundary defined as a sfc_POLYGON object, for boundaries of any shape
  # window: alternative to "boundary" for rectangular boundaries, as c(xmin, ymin, xmax, ymax)
  # radius: the distance to the boundary that separates core from edge points
  # assume_inside: will accelerate if we know all xy points are inside the boundary
  
  # Confirm 2-columns matrix for the point coordinates
  xy <- as.matrix(xy)
  stopifnot(ncol(xy) == 2L)
  storage.mode(xy) <- "double"
  
  # Confirm or guess boundary:
  boundary <- get_boundary(boundary = boundary, window = window, xy = xy)
  
  # Determine which points are inside:
  pts <- sf::st_as_sf(
    as.data.frame(xy),
    coords = 1:2,
    crs = sf::st_crs(boundary),
    remove = FALSE
  )
  inside <- lengths(sf::st_intersects(pts, boundary)) > 0
  
  # Densify boundary once, then use get_distances_up_to_r()
  boundary <- sf::st_boundary(boundary)
  boundary <- sf::st_segmentize(boundary, dfMaxLength = radius / 10) 
  
  boundary_xy <- sf::st_coordinates(boundary)[, c("X", "Y"), drop = FALSE]
  boundary_xy <- unique(boundary_xy)
  storage.mode(boundary_xy) <- "double"
  
  # Identify points near the boundary:
  near <- rep(FALSE, nrow(xy))
  ii <- which(inside)
  if (length(ii)) {
    d <- get_distances_up_to_r(
      query = xy[ii, , drop = FALSE],
      reference = boundary_xy,
      radius = 1.05 * radius # small safety margin
    )
    
    near[ii[unique(d$query_i)]] <- TRUE
  }
  
  # Build output into three classes:
  location <- rep(NA, length(inside))
  location[!inside] <- "out"
  location[inside & near] <- "edge"
  location[inside & !near] <- "core"
  location
}


# The neighborhood analyses are scale-explicit, looking at circles or rings of
# increasing sizes. The user may not have a predefined set of scales at which
# to measure the metrics of interest. This function implements useful defaults,
# given a set of query points:

get_radii <- function(query, rmax = NA, m = 50, case = c("circles", "rings"))
{
  # Default case: circles
  if(length(case) != 1) case = "circles"
  
  # Maximum distance: same defaults as Ripley's K function in spatstat
  if(is.na(rmax))
  {
    a = diff(range(query[,1], na.rm = TRUE))
    b = diff(range(query[,2], na.rm = TRUE))
    lambda = nrow(query) / (a * b)
    rmax = min(c(a/4, b/4, sqrt(1000 / (pi * lambda))))
  }
  
  # Default for circles: simple equal increments in radii
  radii.c <- seq(from = rmax/m, to = rmax, by = rmax/m)
  
  # Default for rings: rings of equal area
  # (It returns the outer radii: the inner radii is 0 and the same sequence but shifted)
  radii.r <- rmax*sqrt((1:m)/m)
  
  # Output
  if(case == "circles") radii <- radii.c
  if(case == "rings")   radii <- radii.r
  radii
}





# An important piece of information to control for edge effects are
# the areas of the neighborhoods at different scales, both for the
# "core" points (easy) and the "edge" points (difficult). The code
# can handle irregular shapes, estimating areas numerically within
# a tolerance to error, defined by the user. To do so, it creates
# dense grids (but not excessively dense) within the valid boundary,
# and then counts points in the grid within the provided radii.
# The function below integrates useful defaults and rules of thumb.

get_area_of_neighborhoods <- function(query = NULL,
                                      boundary = NULL,
                                      window = NULL,
                                      radii = NULL,
                                      rmax = NA,
                                      m = 50,
                                      case = c("circles", "rings"),
                                      tol = 0.05)
{
  
  # query: The set of focal or query points.
  #        2-columns table with x and y coordinates in 1st and 2nd column, respectively.
  # boundary: boundary defined as a sfc_POLYGON object, for boundaries of any shape
  # window: rectangular, alternative to "boundary", as c(xmin, ymin, xmax, ymax)
  # radii: defines the scales of the circular neighborhoods.
  # rmax: maximum radius for neighborhoods, if radii not provided
  # m: number of scales to evaluate, if radii not provided
  # case: circular vs. ring-shaped neighborhoods: will affect area and default radii
  # tol: the accepted error in area, 5% by default.
  
  # Default case: circles
  if(length(case) != 1) case = "circles"
  
  # Confirm or guess boundary
  boundary <- get_boundary(boundary = boundary, window = window, xy = query)
  bb <- sf::st_bbox(boundary)
  xmin = bb["xmin"]
  xmax = bb["xmax"]
  ymin = bb["ymin"]
  ymax = bb["ymax"]
  
  # Define radii, if the user has not provided them yet:
  if(is.null(radii))
    radii <- get_radii(query = query, rmax = rmax, m = m, case = case)
  m = length(radii)
  rmax = max(radii)
  
  # The areas of neighborhoods at different scales will 
  # repeat similar calculations, and will use different
  # grids for each scale. The alternative is to commit
  # to a very dense grid with low error in the worst case
  # scenario, but costly calculations when it is not really
  # needed to achieve a given tolerance. 
  
  # This runs a loop across scales. Before the loop,
  # it creates empty objects to store results.
  areas <- matrix(NA, nrow = nrow(query), ncol = m) # point x scale matrix
  steps <- rep(NA, m)
  for(i in 1:m)
  {
    # The target scale
    radius = radii[i]
    
    # Define the resolution of the grid.
    # Our approach is the Gauss circle problem,
    # and the associated error is not perfectly
    # known, but relatively easy to estimate numerically.
    # There are useful rules of thumb, empirical results, etc.
    
    # This is a useful heuristic for the step associated
    # with a given error in the case "points within circles",
    # when the points come from a regular grid:
    if(case == "circles")
      step = radius * (1.772 * tol)^(2/3)
    
    # The "points within rings" follows a similar
    # logic, but adds more penalization, because
    # of the greater amount of perimeter.
    # Note that outer rings, if of equal area,
    # are thinner and thinner, and perimeter:area
    # ratio grows fast, thus demanding denser and
    # denser grids to constrain error within tolerance.
    if(case == "rings")
    {
      if(i == 1) r0 = 0
      if(i >= 2) r0 = radii[i-1]
      r1 = radius
      r_area = sqrt(r1^2 - r0^2)          # area-equivalent circle radius
      B_ring = (r1 + r0) / r_area         # boundary inflation vs. equal-area circle
      step = r_area * (1.772 * tol)^(2/3) / B_ring^(1/3)
    }
    
    # Because the regular grid is the "easy" case,
    # we need to improve the resolution so it
    # returns the same error around the (possibly complex) boundary.
    # This is a conservative heuristic derived from Steiner / parallel-body argument:
    K = 1 + 4 / pi
    step = step / K
    steps[i] <- step
    
    # Put all the areas together: initialize with full circles,
    # then overwrite with the areas calculated for the edge cases:
    A0 = pi*radius^2
    A <- rep(A0, nrow(query))
    ceo_query <- core_vs_edge(xy = query, boundary = boundary, radius = radius)
    query_edge <- query[ceo_query == "edge",,drop = FALSE]
    
    if(nrow(query_edge) > 0)
    {
      # Note the 2r distance needed for the grid to cover the entire
      # neighborhood of query points that are at r distance from the edge.
      grid <- expand.grid(x = seq(from = xmin - step/2, to = xmax + step/2, by = step),
                          y = seq(from = ymin - step/2, to = ymax + step/2, by = step))
      
      ceo_grid  <- core_vs_edge(xy = grid, boundary = boundary, radius = 2*radius)
      grid_edge <- grid[ceo_grid == "edge",,drop = FALSE]
      
      # Use the neighbor-identifying workhorse to identify points in the grid that
      # are inside the neighborhoods. Their number is a proxy for neighborhood area.
      d <- get_distances_up_to_r(query = query_edge,
                                  reference = grid_edge,
                                  radius = radius)
      
      A[ceo_query == "edge"] <- pmin(A0,
                                     tabulate(d$query_i, nbins = nrow(query_edge)) * step^2)
    }
    
    areas[,i] <- A
    #print(paste0(i, "/", m))
  }
  
  # The ring case is solved from the circles case,
  # just by sequential substraction:
  if(case == "rings")
  {
    large <- areas[, 2:m, drop = FALSE]
    small <- areas[, 1:(m-1), drop = FALSE]
    areas <- cbind(areas[,1], large - small) # the inner ring is a circle, zero substraction
  }
  
  # Output:
  if(case == "circles")
    scales <- data.frame(radius = radii,
                         grid_step = steps,
                         grid_points_per_circle = pi * radii^2 / steps^2,
                         error = tol)
  
  if(case == "rings")
  {
    inner_radii <- c(0, radii[-m])
    scales <- data.frame(inner_radius = inner_radii,
                         outer_radius = radii,
                         grid_step = steps,
                         grid_points_per_ring = pi * (radii^2 - inner_radii^2) / steps^2,
                         error = tol)
  }
  
  list(areas = areas, case = case, scales = scales, boundary = boundary)
}








