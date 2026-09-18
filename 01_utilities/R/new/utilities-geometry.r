# <function> <name> perpendicular.distance </name> <description> Distance from a point to a line (so it's
# the perpendicular distance); m and b are slope and intercept; x and y are coordinates. If both b,m and x,y are vectors, they
# must all be same length.  Note check for infinite slope, meaning that the intercept b is the x-intercept.

# </description> <display>true</display> <update>true</update> <arguments> b: y-intercept m: line slope </arguments> <sample>
# </sample> <source>
perpendicular.distance <- function(b, m, x, y) {
    len <- ifelse(length(b) > length(x), length(b), length(x))
    minf <- is.infinite(m)
    result <- rep(NA, len)
    
    result <- sqrt((y - m * x - b)^2/(1 + m^2))
    result[minf] <- abs(x - b)[minf]
    
    return(result)
}
# </source> </function>

# <function> <name> insideRectangle </name> <description> Checks a vector of coordinates x, y to return which are inside a
# rectangle. For a much more general function for checking whether points are inside polygons, use the function inout() in the
# package splancs.  </description> <arguments> </arguments> <sample> </sample> <source>
insideRectangle <- function(x, y, xrange, yrange) {
    inX <- x >= xrange[1] & x < xrange[2] & !is.na(x)
    inY <- y >= yrange[1] & y < yrange[2] & !is.na(y)
    return(inX & inY)
}
# </source> </function>


# <function> <name> are.ptsinside </name> <description> Checks many points (dataframe pt with x and y) against a single quadrat
# whose corners are given by as xlo, ylo, xhi, yhi.  It returns a logical vector, TRUE for the points inside. This is same as
# insideRectangle, but accepting input as a matrix pts and a single vector of the four corners of the rectange.  </description>
# <arguments> </arguments> <sample> </sample> <source>
are.ptsinside <- function(pts, coord) {
    return(insideRectangle(x = pts[, 1], y = pts[, 2],
        xrange = coord[c(1, 3)], yrange = coord[c(2, 4)]))
}

# </source> </function>

# <function> <name> ispt.inside </name> <description> Check a single pt (x and y) against a large number of quadrats whose
# corners are given by the rows of coord, xlo, ylo, xhi, yhi.  It returns the fraction of quadrats which the point falls
# inside. This is exactly like are.ptsinside() but allows there to be many rectangles, defined by a dataframe coord.
# </description> <arguments> </arguments> <sample> </sample> <source>
ispt.inside <- function(pt, coord) {
    norect <- dim(coord)[1]
    inside <- logical(norect)
    
    for (i in seq_len(norect)) inside[i] <- are.ptsinside(t(pt), drop(as.matrix(coord[i, ])))
    
    return(mean(inside))
}
# </source> </function>
