# <function> <name> readelevdata </name> <description> A function to read a table of elevation data into a matrix appropriate
# for mapping. It can read a text table (if text==TRUE) or a dataframe.  In either case, there must be x then y coordinates,
# followed by elevation, at every corner across the plot, using grid of gridsize. The output is a matrix, as needed for R's
# contour function.

# </description> <arguments> </arguments> <sample> </sample> <source>
readelevdata <- function(elevfile, gridsize = 5, text = FALSE) {
    if (text) 
        elevdata <- read.table(elevfile, header = T, sep = "\t") else elevdata <- elevfile
    
    xdim <- max(elevdata$x)
    ydim <- max(elevdata$y)
    
    elevmat <- matrix(elevdata$elev, nrow = 1 + ydim/gridsize, ncol = 1 + xdim/gridsize, byrow = F)
    
    return(elevmat)
}
# </source> </function>

# <function> <name> elev.to.list </name> <description> A function which reads a dataframe with x,y,elevation for a given grid
# size and converts to a list

# note the names of the columns in the dataframe must be: x y elev

# </description> <arguments> object 1 in list is: the input dataframe object 2 in list is: matrix of elevation value sin the
# orientation of the plot

# </arguments> <sample> </sample> <source>
elev.to.list <- function(elevfile, gridsize = 5) {
    xdim <- max(elevfile$x)
    ydim <- max(elevfile$y)
    
    elevmat <- matrix(elevfile$elev, nrow = 1 + ydim/gridsize, ncol = 1 + xdim/gridsize, byrow = F)
    
    return(list(col = elevfile, mat = elevmat))
}



# </source> </function> <function> <name> allquadratslopes </name> <description> Calculating slope of all quadrats in a plot
# calcslope takes 3 elevations and finds the slope of the plane through them quadslope divides the 4 corners of a quadrat into
# 4 different groups of 3 stakes, takes the slope of each, then averages these were first written in C++ see slopeelev.cpp for
# more on the geometry allquadratslopes goes through all 20x20 m quadrats in a plot and finds the slope, mean elevation, and
# convexity of each convexity is the mean elevation of one 20x20 m quadrat relative (minus) the mean of its immediate neighbors

# Helene Muller-Landau added a section to correct convexity in edge quadrats

# </description> <arguments> </arguments> <sample> </sample> <source>
allquadratslopes <- function(elev, gridsize = 20, plotdim = c(1000, 500), edgecorrect = TRUE) {
    rw <- cl <- 0
    on.exit(cat(rw, " ", cl, "\n"))
    
    quadrat.locations <- make_quadrats_fast(
      outer.xmin = 0,
      outer.ymin = 0,
      outer.xmax = plotdim[1],
      outer.ymax = plotdim[2],
      quadrat.side.x = gridsize,
      quadrat.side.y = gridsize
    )
    
    totalquads <- nrow(quadrat.locations)
    
    x.breaks <- sort(unique(c(
      quadrat.locations$xmin,
      quadrat.locations$xmax
    )))
    
    y.breaks <- sort(unique(c(
      quadrat.locations$ymin,
      quadrat.locations$ymax
    )))
    
    columns <- length(x.breaks)
    rows <- length(y.breaks)
    
    cat("Calculating topographic indices for ", totalquads, " quadrats\n")
    
    elevdata <- elev$col[elev$col$x%%gridsize == 0 & elev$col$y%%gridsize == 0, ]
    elevmat <- matrix(elevdata$elev, nrow = rows, ncol = columns, byrow = F)
    
    meanelev <- convex <- convex2 <- slope <- numeric()
    corner <- sideht <- numeric()
    
    # Mean elevation and slope of each quadrat.
    for(i in seq_len(totalquads))
    {
      q <- quadrat.locations[i, ]
      
      column.i <- match(q$xmin, x.breaks)
      row.i <- match(q$ymin, y.breaks)
      
      corner[1] <- elevmat[row.i, column.i]
      corner[2] <- elevmat[row.i + 1, column.i]
      corner[3] <- elevmat[row.i + 1, column.i + 1]
      corner[4] <- elevmat[row.i, column.i + 1]
      
      meanelev[i] <- mean(corner)
      slope[i] <- quadslope(corner, gridsize = gridsize)[1]
      
      if(i %% 1000 == 0)
        cat("Finding elevation and slope of quadrat ", i, "\n")
    }
    
    if(!requireNamespace("FNN", quietly = TRUE))
      stop("install package 'FNN' to calculate neighboring quadrats")
    
    if(totalquads < 2)
      stop("at least two quadrats are required to calculate convexity")
    
    centroids <- cbind(
      x = (quadrat.locations$xmin + quadrat.locations$xmax) / 2,
      y = (quadrat.locations$ymin + quadrat.locations$ymax) / 2
    )
    
    k <- min(8, totalquads - 1)
    knn <- FNN::get.knn(centroids, k = k)
    
    maximum.distance <- sqrt(2) * gridsize
    
    # Convexity relative to quadrats sharing an edge or corner.
    for(i in seq_len(totalquads))
    {
      neighbor.quads <- knn$nn.index[i, ]
      neighbor.distances <- knn$nn.dist[i, ]
      
      neighbor.quads <- neighbor.quads[
        neighbor.distances <=
          maximum.distance + sqrt(.Machine$double.eps)
      ]
      
      convex[i] <- meanelev[i] - mean(meanelev[neighbor.quads])
      
      if(i %% 1000 == 0)
        cat("Finding convexity of quadrat ", i, "\n")
    }
    
    # correcting convexity in edge quadrats, based on center of the 20x20 rather than surrounding 20x20s. This requires that the
    # elev$mat has an elevation at the middle of every grid cell.
    if(edgecorrect)
    {
      edge.quadrats <- which(
        quadrat.locations$xmin == 0 |
          quadrat.locations$ymin == 0 |
          quadrat.locations$xmax == plotdim[1] |
          quadrat.locations$ymax == plotdim[2]
      )
      
      for(i in edge.quadrats)
      {
        midx <- (
          quadrat.locations$xmin[i] +
            quadrat.locations$xmax[i]
        ) / 2
        
        midy <- (
          quadrat.locations$ymin[i] +
            quadrat.locations$ymax[i]
        ) / 2
        
        middle.i <- which(
          elev$col$x == midx &
            elev$col$y == midy
        )
        
        if(length(middle.i) != 1)
          stop("one centre elevation is required for every edge quadrat")
        
        midelev <- elev$col$elev[middle.i]
        convex[i] <- midelev - meanelev[i]
      }
    }
    
    return(data.frame(meanelev = meanelev, convex = convex, slope = slope))
}



# </source> </function> <function> <name> quadslope </name> <description> Given the elevation at four corners of a square of
# side=gridsize, this estimates the slope in degrees of the terrain over that square. The slope is calculated for each of the 4
# planes defined by omitting one of the points.  The 4 slopes are averaged. Returns both the mean slope and the sqrt of the
# variance of the 4 different slopes.  </description> <arguments> cornerelev: vector of 4 elevations <br> gridsize: the side of
# the square </arguments> <sample> </sample> <source>
quadslope <- function(cornerelev, gridsize = 20) {
    slope <- numeric(4)
    z <- numeric(3)
    
    for (j in 1:4) {
        post <- 1
        
        for (k in (j + 1):(j + 3)) {
            if (k > 4) 
                m <- k%%4 else m <- k
            
            z[post] <- cornerelev[m]
            post <- post + 1
        }
        
        slope[j] <- calcslope(z, gridsize)
    }
    
    return(c(mean(slope), sqrt(var(slope))))
}

# </source> </function> <function> <name> calcslope </name> <description> Given the z-coordinate of 3 points in 3D space whose
# x-y coordinates are corners of a right, isoceles triangle whose short side = gridsize, this calculates the slope in degrees
# of the plane through the points.  The second point is the one between the two short sides.  </description> <arguments> z:
# numeric vector of length 3 <br> gridsize: numeric scalar </arguments> <sample> </sample> <source>
calcslope <- function(z, gridsize = 20) {
    z2 <- z[3] - z[2]
    z1 <- z[1] - z[2]
    
    if (z2 == 0) {
        if (z1 == 0) 
            return(0) else denom <- sqrt(1 + (gridsize/z1)^2)
        
        theta1 <- acos((-gridsize/z1)/denom)
        theta2 <- acos((gridsize/z1)/denom)
    } else {
        denom <- sqrt(1 + (z1/z2)^2 + (gridsize/z2)^2)
        
        theta1 <- acos((-gridsize/z2)/denom)
        theta2 <- acos((gridsize/z2)/denom)
    }
    
    if (theta1 <= theta2) 
        return(180 * theta1/pi) else return(180 * theta2/pi)
}



# </source> </function> <function> <name> calc.gradient </name> <description> Calculate flow using Seibert & McGlynn algorithm.
# Takes a 3x3 matrix of elevations and works on central point; also requires grid size (usually 20 m)
# z=matrix(c(c(268.7,275.9,283.2),c(275.9,282.8,290.0),c(283.2,290.0,297)),nrow=3,byrow=TRUE)

# This calculates the gradient for the 8 triangular facets around the center point, following Seibert & McGlynn The output is a
# data.frame of direction and slope for the 8 facets, starting with the lower left and moving clockwise

# </description> <arguments> </arguments> <sample> </sample> <source>
calc.gradient <- function(elev, grid = 20) {
    elev.ctr <- elev[2, 2]
    elev.corner <- c(elev[1, 1], elev[1, 3], elev[3, 3], elev[3, 1])
    elev.side <- c(elev[2, 1], elev[2, 3], elev[2, 3], elev[1, 2])
    
    x.ctr <- y.ctr <- grid
    x.corner <- c(0, 0, 2 * grid, 2 * grid)
    y.corner <- c(0, 2 * grid, 2 * grid, 0)
    x.side <- c(0, grid, 2 * grid, grid)
    y.side <- c(grid, 2 * grid, grid, 0)
    # browser()
    
    corner <- side <- matrix(ncol = 2, nrow = 4)
    
    for (i in 1:4) {
        z1 <- elev.corner[i] - elev.ctr
        z2 <- elev.side[i] - elev.ctr
        
        x1 <- (x.corner[i] - x.ctr)
        x2 <- (x.side[i] - x.ctr)
        
        y1 <- (y.corner[i] - y.ctr)
        y2 <- (y.side[i] - y.ctr)
        
        corner[i, ] <- calc.directionslope(z1, z2, x1, x2, y1, y2)
    }
    
    for (i in 1:4) {
        z2 <- elev.corner[i] - elev.ctr
        z1 <- elev.side[i] - elev.ctr
        
        x2 <- x.corner[i] - x.ctr
        x1 <- x.side[i] - x.ctr
        
        y2 <- y.corner[i] - y.ctr
        y1 <- y.side[i] - y.ctr
        
        side[i, ] <- calc.directionslope(z1, z2, x1, x2, y1, y2)
    }
    
    result <- matrix(ncol = 2, nrow = 8)
    colnames(result) <- c("direction", "slope")
    result[c(1, 3, 5, 7), ] <- corner
    result[c(2, 4, 6, 8), ] <- side
    
    return(data.frame(result))
}



# </source> </function> <function> <name> calc.directionslope </name> <description> This runs equations 1-3 of Seibert &
# McGlynn

# </description> <arguments> </arguments> <sample> </sample> <source>
calc.directionslope <- function(z1, z2, x1, x2, y1, y2) {
    normal <- c(z1 * y2 - z2 * y1, z1 * x2 - z2 * x1, y1 * x2 - y2 * x1)
    # browser()
    
    if (normal[1] == 0) {
        if (normal[2] >= 0) 
            direction <- 0 else direction <- pi
    } else if (normal[1] > 0) 
        direction <- pi/2 - atan(normal[2]/normal[1]) else direction <- 3 * pi/2 - atan(normal[2]/normal[1])
    
    slope <- (-1) * tan(acos(normal[3]/sqrt(sum(normal^2))))
    
    return(c(direction, slope))
}


# </source> </function>
