# Functions to handle quadrats, splits, local coordinates, and global coordinates
# Gabriel Arellano
# gabriel.arellano.torres@gmail.com


##################################################################################
### FUNCTIONS THAT RETURN QUADRAT LOCATIONS
##################################################################################

# Auxiliary functions to make different types of sequences of letters:

# 1) First n labels: a, b, ..., z, aa, ab, ..., az, ba, ...
letter_seq <- function(n) {
  f <- function(i) {
    x <- ""
    while (i > 0) {
      i <- i - 1
      x <- paste0(letters[i %% 26 + 1], x)
      i <- i %/% 26
    }
    x
  }
  vapply(seq_len(n), f, character(1))
}

# 2) Labels with exactly k characters:
# k = 2 gives aa, ab, ..., az, ba, bb, ..., bz, ...
letter_seq_fixed <- function(n, k) {
  grid <- expand.grid(rep(list(letters), k), stringsAsFactors = FALSE)
  do.call(paste0, rev(grid))[1:n]
}


# AUXILIARY TO ENSURE NON-OVERLAPPING QUADRATS
 
make_quadrats_non_overlapping <- function(quadrats, k = 8)
{
  # Required structure
  needed <- c("quadrat_id", "xmin", "ymin", "xmax", "ymax")
  missing <- setdiff(needed, names(quadrats))
  
  if(length(missing) > 0)
    stop("quadrats is missing: ", paste(missing, collapse = ", "))
  
  if(!requireNamespace("FNN", quietly = TRUE))
    stop("install package 'FNN' to use make_quadrats_non_overlapping()")
  
  if(nrow(quadrats) < 2)
    return(quadrats)
  
  # Original centroids define relative quadrat positions.
  cx <- (quadrats$xmin + quadrats$xmax) / 2
  cy <- (quadrats$ymin + quadrats$ymax) / 2
  xy <- cbind(cx, cy)
  
  # Find local neighbors only.
  k <- min(k, nrow(quadrats) - 1)
  nn <- FNN::get.knn(xy, k = k)$nn.index
  
  # Unique neighbor pairs.
  pairs <- data.frame(
    i = rep(seq_len(nrow(quadrats)), each = k),
    j = as.vector(t(nn))
  )
  
  pairs <- data.frame(
    i = pmin(pairs$i, pairs$j),
    j = pmax(pairs$i, pairs$j)
  )
  
  pairs <- unique(pairs)
  
  for(p in seq_len(nrow(pairs)))
  {
    i <- pairs$i[p]
    j <- pairs$j[p]
    
    # Check whether rectangles overlap. Touching is OK.
    overlap.x <- quadrats$xmin[i] < quadrats$xmax[j] &
      quadrats$xmax[i] > quadrats$xmin[j]
    
    overlap.y <- quadrats$ymin[i] < quadrats$ymax[j] &
      quadrats$ymax[i] > quadrats$ymin[j]
    
    if(!overlap.x | !overlap.y)
      next
    
    # Relative position of centroids.
    dx <- cx[j] - cx[i]
    dy <- cy[j] - cy[i]
    
    if(dx == 0 & dy == 0)
    {
      warning("two quadrats have the same centroid; overlap could not be resolved")
      next
    }
    
    # Separate along the dominant centroid-separation axis.
    horizontal <- abs(dx) >= abs(dy)
    
    # j is mostly to the right of i.
    if(horizontal & dx > 0)
    {
      cut <- (cx[i] + cx[j]) / 2
      quadrats$xmax[i] <- min(quadrats$xmax[i], cut)
      quadrats$xmin[j] <- max(quadrats$xmin[j], cut)
    }
    
    # j is mostly to the left of i.
    if(horizontal & dx < 0)
    {
      cut <- (cx[i] + cx[j]) / 2
      quadrats$xmin[i] <- max(quadrats$xmin[i], cut)
      quadrats$xmax[j] <- min(quadrats$xmax[j], cut)
    }
    
    # j is mostly above i.
    if(!horizontal & dy > 0)
    {
      cut <- (cy[i] + cy[j]) / 2
      quadrats$ymax[i] <- min(quadrats$ymax[i], cut)
      quadrats$ymin[j] <- max(quadrats$ymin[j], cut)
    }
    
    # j is mostly below i.
    if(!horizontal & dy < 0)
    {
      cut <- (cy[i] + cy[j]) / 2
      quadrats$ymin[i] <- max(quadrats$ymin[i], cut)
      quadrats$ymax[j] <- min(quadrats$ymax[j], cut)
    }
  }
  
  # Sanity check.
  bad <- quadrats$xmax <= quadrats$xmin |
    quadrats$ymax <= quadrats$ymin
  
  if(any(bad))
    warning(sum(bad), " quadrat(s) ended with invalid dimensions after overlap adjustment")
  
  quadrats
}



# USING PLOT & QUADRAT DIMENSIONS, AND DESCRIPTORS OF CONVENTION
get_quadrat_locations_from_id_convention <- function(
  # PARAMETERS
  # Dimensions of the large unit
  outer.xmin = NA,
  outer.ymin = NA,
  outer.xmax = NA,
  outer.ymax = NA,
  
  # Dimensions of the smaller sub-units
  quadrat.side.x = NA,
  quadrat.side.y = NA,
  
  # How the ID is built
  built = c("col+row", "row+col", "sequence"),
  example.part.for.col = "1", # c("A", "1", "01", "000001")
  example.part.for.row = "1", # c("A", "1", "01", "000001")
  example.part.for.sequence = "1", # c("0001", "1", "A")
  prefix = "",
  suffix = "",
  separator = "", # c("-", ",", " ")
  
  # Even more flexibility, overwriting:
  seq.of.row.parts.along.one.col = NA, # seq()
  seq.of.col.parts.along.one.row = NA, # seq()
  entire.sequence   = NA, # seq()
  
  # How the sequences will run
  start.col.at = NA,
  start.row.at = start.col.at,
  start.sequence.at = 1,
  along.columns.first.then.next.column = TRUE,
  along.rows.first.then.next.row = FALSE,
  changing.directions = FALSE # south to north in first column, north to south in second, etc
  )

{
  
  # Checks, stops, warnings, and defaults
  
  if(is.na(quadrat.side.x) | is.null(quadrat.side.x))
    if(is.na(quadrat.side.y) | is.null(quadrat.side.y))
      stop("specify quadrat size through quadrat.side.x and quadrat.side.y")
  
  if(is.na(outer.xmax) | is.null(outer.xmax) | is.na(outer.ymax) | is.null(outer.ymax))
    stop("specify dimensions through outer.xmax and outer.ymax")
  
  if(is.na(outer.xmin) | is.null(outer.xmin))
  {
    outer.xmin = 0
    warning("outer.xmin = 0 assumed by default")
  }
    
  if(is.na(outer.ymin) | is.null(outer.ymin))
  {
    outer.ymin = 0
    warning("outer.ymin = 0 assumed by default")
  }
  
  if(length(built) != 1)
  {
    built = "sequence"
    warning("built = 'sequence' assumed by default")
  }
  
  if(is.na(quadrat.side.x) | is.null(quadrat.side.x))
  {
    quadrat.side.x <- quadrat.side.y
    warning("quadrat.side.x assumed equal to quadrat.side.y")
  }
  
  if(is.na(quadrat.side.y) | is.null(quadrat.side.y))
  {
    quadrat.side.y <- quadrat.side.x
    warning("quadrat.side.y assumed equal to quadrat.side.x")
  }
  
  if(is.na(start.col.at) | is.null(start.col.at))
  {
    start.col.at = as.numeric(example.part.for.col)
    warning("start.col.at taken from example.part.for.col by default")
  }
  
  if(is.na(start.row.at) | is.null(start.row.at))
  {
    start.row.at = as.numeric(example.part.for.row)
    warning("start.row.at taken from example.part.for.row by default")
  }
  
  
  # Generating sub-units within units:
  # (Small leftovers preferred over unusually large quadrats)
  xmins <- seq(from = outer.xmin, to = outer.xmax, by = quadrat.side.x)
  xmins <- xmins[xmins < outer.xmax]
  ymins <- seq(from = outer.ymin, to = outer.ymax, by = quadrat.side.y)
  ymins <- ymins[ymins < outer.ymax]
  
  # Add corner coordinates and col+row indexes
  grid.exact  <- expand.grid(xmin = xmins, ymin = ymins)
  grid.exact$xmax <- pmin(outer.xmax, grid.exact$xmin + quadrat.side.x)
  grid.exact$ymax <- pmin(outer.ymax, grid.exact$ymin + quadrat.side.y)
  grid.colrow <- expand.grid(col = 1:length(xmins), row =  1:length(ymins))
  grid <- data.frame(grid.exact, grid.colrow)
  head(grid)
  
  # Reorder the table if we are going to add
  # one single simple sequence for all sub-units:
  if(built == "sequence")
  {
    # avoid conflicts
    miss1 = is.na(along.columns.first.then.next.column) | is.null(along.columns.first.then.next.column)
    miss2 = is.na(along.rows.first.then.next.row) | is.null(along.rows.first.then.next.row)
    if(miss1 & miss2) "clarify whether labelling is along columns first, or along rows first"
    
    if(along.columns.first.then.next.column == along.rows.first.then.next.row)
      stop("clarify whether labelling is along columns first, or along rows first")
    
    # re-order
    if(along.columns.first.then.next.column)
      o <- order(grid$col, grid$row)
    
    if(along.rows.first.then.next.row)
      o <- order(grid$row, grid$col)
    
    grid <- grid[o,]
    
    # change directions if necessary, for even columns or rows
    if(changing.directions)
    {
      
      if(along.columns.first.then.next.column)
      {
        for(i in seq(from = 2, to = max(grid$col), by = 2))
        {
          chunk <- which(grid$col == i)
          grid[chunk,] <- grid[rev(chunk),]
        }
      }
      
      if(along.rows.first.then.next.row)
      {
        for(i in seq(from = 2, to = max(grid$row), by = 2))
        {
          chunk <- which(grid$row == i)
          grid[chunk,] <- grid[rev(chunk),]
        }
      }
      
    } # end of if changing directions
  } # done reordering the table if we are going to add sequence
  
  # Identify whether the sequences are numbers or letters:
  if(built %in% c("col+row", "row+col"))
  {
    #use.numbers.for.col1 = is.numeric(as.numeric(example.part.for.col)) # would fail: is.numeric(as.numeric("a")) = TRUE
    #use.numbers.for.col2 = is.numeric(as.numeric(start.col.at))         # would fail: is.numeric(as.numeric("a")) = TRUE
    use.numbers.for.col1 = !is.na(suppressWarnings(as.numeric(example.part.for.col)))
    use.numbers.for.col2 = !is.na(suppressWarnings(as.numeric(start.col.at)))
    
    if(use.numbers.for.col1 != use.numbers.for.col2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for 'col' part")
    
    #use.numbers.for.row1 = is.numeric(as.numeric(example.part.for.row))  # would fail: is.numeric(as.numeric("a")) = TRUE
    #use.numbers.for.row2 = is.numeric(as.numeric(start.row.at))          # would fail: is.numeric(as.numeric("a")) = TRUE
    use.numbers.for.row1 = !is.na(suppressWarnings(as.numeric(example.part.for.row)))
    use.numbers.for.row2 = !is.na(suppressWarnings(as.numeric(start.col.row)))
    
    if(use.numbers.for.row1 != use.numbers.for.row2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for 'row' part")
    
    use.numbers.for.col = use.numbers.for.col1
    use.numbers.for.row = use.numbers.for.row1
    
    use.letters.for.col1 = grepl('([a-zA-Z])', example.part.for.col)
    use.letters.for.col2 = grepl('([a-zA-Z])', start.col.at)
    if(use.letters.for.col1 != use.letters.for.col2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for 'col' part")
    
    use.letters.for.row1 = grepl('([a-zA-Z])', example.part.for.row)
    use.letters.for.row2 = grepl('([a-zA-Z])', start.row.at)
    if(use.letters.for.row1 != use.letters.for.row2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for 'row' part")
    
    use.letters.for.col = use.letters.for.col1
    use.letters.for.row = use.letters.for.row1
    
    if(use.numbers.for.col == use.letters.for.col)
      stop("clarify between numbers and letters for 'col' part")
    if(use.numbers.for.row == use.letters.for.row)
      stop("clarify between numbers and letters for 'row' part")
  }
  
  if(built == "sequence")
  {
    #use.numbers.for.sequence1 = is.numeric(as.numeric(example.part.for.sequence)) # would fail: is.numeric(as.numeric("a")) = TRUE
    #use.numbers.for.sequence2 = is.numeric(as.numeric(start.sequence.at))         # would fail: is.numeric(as.numeric("a")) = TRUE
    use.numbers.for.sequence1 = !is.na(suppressWarnings(as.numeric(example.part.for.sequence)))
    use.numbers.for.sequence2 = !is.na(suppressWarnings(as.numeric(start.sequence.at)))
    
    if(use.numbers.for.sequence1 != use.numbers.for.sequence2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for the sequence")
    
    use.numbers.for.sequence = use.numbers.for.sequence1
    
    use.letters.for.sequence1 = grepl('([a-zA-Z])', example.part.for.sequence)
    use.letters.for.sequence2 = grepl('([a-zA-Z])', start.sequence.at)
    if(use.letters.for.sequence1 != use.letters.for.sequence2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for the sequence")
    
    use.letters.for.sequence = use.letters.for.sequence1
    
    if(use.numbers.for.sequence == use.letters.for.sequence)
      stop("clarify between numbers and letters for the sequence")
  }
  
  # Build the sequences for the col and row parts:
  if(built %in% c("col+row", "row+col"))
  {
    # Get the basic sequences first, then expand
    if(use.numbers.for.col)
    {
      # checks
      if(as.numeric(example.part.for.col) == 0 & start.col.at > 0)
        warning(paste0("start.col.at = ", start.col.at, " takes priority, although example suggests start = 0"))
      
      # initialize
      basic.col.seq <- seq(from = as.numeric(start.col.at), length.out = length(unique(grid$col)), by = 1)
      
      # add zeros if needed
      add <- pmax(0, nchar(example.part.for.col) - nchar(basic.col.seq))
      basic.col.seq <- paste0(strrep("0", add), basic.col.seq)
    } # end of basic sequence to build col parts with numbers
    
    if(use.numbers.for.row)
    {
      # checks
      if(as.numeric(example.part.for.row) == 0 & start.row.at > 0)
        warning(paste0("start.row.at = ", start.row.at, " takes priority, although example suggests start = 0"))
      
      # initialize
      basic.row.seq <- seq(from = as.numeric(start.row.at), length.out = length(unique(grid$row)), by = 1)
      
      # add zeros if needed
      add <- pmax(0, nchar(example.part.for.row) - nchar(basic.row.seq))
      basic.row.seq <- paste0(strrep("0", add), basic.row.seq)
    } # end of basic sequence to build row parts with numbers
    
    if(use.letters.for.col)
    {
      l = FALSE
      u = FALSE
      if(substr(example.part.for.col, 1, 1) %in% letters) l = TRUE
      if(substr(example.part.for.col, 1, 1) %in% LETTERS) u = TRUE
      if(!l & !u) stop("the letter in 'example.part.for.col' cannot be recognized")
      n = length(unique(grid$col))
      
      if(nchar(example.part.for.col) == 1)
        basic.col.seq <- letter_seq(n)
      
      if(nchar(example.part.for.col) >  1)
      {
        k = ceiling(log(n, base = 26))
        if(nchar(example.part.for.col) != k)
          stop("unclear example part for 'col': you need a different number of characters")
        
        basic.col.seq <- letter_seq_fixed(n = n, k = k)
      }
      
      if(u) basic.col.seq <- toupper(basic.col.seq)  
    } # end of basic sequence to build col parts with letters
    
    if(use.letters.for.row)
    {
      l = FALSE
      u = FALSE
      if(substr(example.part.for.row, 1, 1) %in% letters) l = TRUE
      if(substr(example.part.for.row, 1, 1) %in% LETTERS) u = TRUE
      if(!l & !u) stop("the letter in 'example.part.for.row' cannot be recognized")
      n = length(unique(grid$row))
      
      if(nchar(example.part.for.row) == 1)
        basic.row.seq <- letter_seq(n)
      
      if(nchar(example.part.for.row) >  1)
      {
        k = ceiling(log(n, base = 26))
        if(nchar(example.part.for.row) != k)
          stop("unclear example part for 'row': you need a different number of characters")
        
        basic.row.seq <- letter_seq_fixed(n = n, k = k)
      }
      
      if(u) basic.row.seq <- toupper(basic.row.seq)  
    } # end of basic sequence to build row parts with letters
    
    # Expand and align with the grid object:
    grid$col.part <- basic.col.seq[grid$col]
    grid$row.part <- basic.row.seq[grid$row]
  }
  
  # Build the entire sequence, if that's the convention:
  if(built == "sequence")
  {
    if(use.numbers.for.sequence)
    {
      # initialize, then add zeros if needed
      s <- seq(from = as.numeric(start.sequence.at), length.out = nrow(grid), by = 1)
      add <- pmax(0, nchar(example.part.for.sequence) - nchar(s))
      grid$entire <- paste0(strrep("0", add), s)
    }
    
    if(use.letters.for.sequence)
    {
      l = FALSE
      u = FALSE
      if(substr(example.part.for.sequence, 1, 1) %in% letters) l = TRUE
      if(substr(example.part.for.sequence, 1, 1) %in% LETTERS) u = TRUE
      if(!l & !u) stop("the letter in 'example.part.for.sequence' cannot be recognized")
      n = nrow(grid)
      
      if(nchar(example.part.for.sequence) == 1)
        s <- letter_seq(n)
      
      if(nchar(example.part.for.sequence) >  1)
      {
        k = ceiling(log(n, base = 26))
        if(nchar(example.part.for.sequence) != k)
          stop("unclear example.part.for.sequence: you need a different number of characters")
        
        s <- letter_seq_fixed(n = n, k = k)
      }
      
      if(u) s <- toupper(s)  
      grid$entire <- s
      
    }
    
  }
  
  
  # If sequences for the col parts or row parts are provided, use them:
  if(!is.na(seq.of.col.parts.along.one.row) & !is.null(seq.of.col.parts.along.one.row))
  {
    if(length(seq.of.col.parts.along.one.row) != length(unique(grid$col)))
      stop("the sequence provided for 'col' parts does not match number of columns -- perhaps it's our 'leftover' convention")
    
    grid$col.part <- seq.of.col.parts.along.one.row[grid$col]
  }
  
  if(!is.na(seq.of.row.parts.along.one.col) & !is.null(seq.of.row.parts.along.one.col))
  {
    if(length(seq.of.row.parts.along.one.col) != length(unique(grid$row)))
      stop("the sequence provided for 'row' parts does not match number of rows -- perhaps it's our 'leftover' convention")
    
    grid$row.part <- seq.of.row.parts.along.one.col[grid$row]
  }
  
  
  # Combine col parts and row parts:
  if(built == "col+row")
    grid$entire <- paste0(grid$col.part, separator, grid$row.part)
  
  if(built == "row+col")
    grid$entire <- paste0(grid$row.part, separator, grid$col.part)
  
  
  # Add prefix and suffix, if there are any:
  if(!is.na(prefix) & !is.null(prefix))
    grid$entire <- paste0(prefix, grid$entire)
  
  if(!is.na(suffix) & !is.null(suffix))
    grid$entire <- paste0(suffix, grid$entire)
  
  
  # If an entire sequence was provided, then use it as it is!
  if(!is.na(entire.sequence) & !is.null(entire.sequence))
  {
    warning("you provided the entire sequence, so that is being used")
    grid$entire <- entire.sequence
  }
  
  
  # Other warnings:
  if(any(grid$xmax - grid$xmin != quadrat.side.x))
    warning("the quadrats don't fit perfectly along x: quadrats at the end are smaller")
  
  if(any(grid$ymax - grid$ymin != quadrat.side.y))
    warning("the quadrats don't fit perfectly along y: quadrats at the end are smaller")
  
  
  # Trim and rename the output:
  out <- grid[,c("entire", "xmin", "ymin", "xmax", "ymax"), drop = FALSE]
  colnames(out)[1] <- "quadrat_id"
  rownames(out) <- out$quadrat_id

  # Sort in the same way, for consistency:
  o <- order(out$xmin, out$ymin)
  out <- out[o,,drop = FALSE]
  out
}
  



# USING A SET OF POINTS WITH COORDINATES
get_quadrat_locations_from_outer_coordinates <- function(
    quadrat.id = NA, # the sequence of quadrat membership of a set of points
    outer.x = NA, # the x coordinates of a set of points
    outer.y = NA, # the y coordinates of a set of points
    outer.xmin = NA, # the xmin of the area in which the points are located, usually 0
    outer.ymin = NA, # the ymin of the area in which the points are located, usually 0
    outer.xmax = NA, # the xmax of the area in which the points are located, ~ plot dimension along x
    outer.ymax = NA,  # the ymax of the area in which the points are located, ~ plot dimension along y
    regular = TRUE # will assume standard regular grids, complete and non-overlapping
)
{
  # Defaults for outer dimensions
  if(is.na(outer.xmin) | is.null(outer.xmin)) outer.xmin = floor(min(outer.x, na.rm = TRUE))
  if(is.na(outer.ymin) | is.null(outer.ymin)) outer.ymin = floor(min(outer.y, na.rm = TRUE))
  if(is.na(outer.xmax) | is.null(outer.xmax)) outer.xmax = ceiling(max(outer.x, na.rm = TRUE))
  if(is.na(outer.ymax) | is.null(outer.ymax)) outer.ymax = ceiling(max(outer.y, na.rm = TRUE))
  
  # Estimate the locations of the quadrats based on the locations of the points
  xmins <- aggregate(outer.x, by = list(quadrat.id), FUN = min, na.rm = TRUE)
  ymins <- aggregate(outer.y, by = list(quadrat.id), FUN = min, na.rm = TRUE)
  xmaxs <- aggregate(outer.x, by = list(quadrat.id), FUN = max, na.rm = TRUE)
  ymaxs <- aggregate(outer.y, by = list(quadrat.id), FUN = max, na.rm = TRUE)
  
  colnames(xmins) <- c("quadrat_id", "xmin")
  colnames(xmaxs) <- c("quadrat_id", "xmax")
  colnames(ymins) <- c("quadrat_id", "ymin")
  colnames(ymaxs) <- c("quadrat_id", "ymax")
  
  locations <- purrr::reduce(.x = list(xmins, ymins, xmaxs, ymaxs),
                             merge, by = c("quadrat_id"), all = TRUE)
  
  rownames(locations) <- locations$quadrat_id
  
  # The quadrats estimated from the locations, as they 
  # adjust exactly to the observations, will be smaller
  # than the true quadrats in general. This is an adjustment
  # to increase quadrat dimensions, under the assumption
  # of random distribution of points. 
  nq <- table(na.omit(data.frame(quadrat.id, outer.x, outer.y))[, 1])
  adj <- (nq[as.character(locations$quadrat_id)] + 1) / (nq[as.character(locations$quadrat_id)] - 1)
  adj[adj > 2] <- 2 # capped  
  buf.x <- (locations$xmax - locations$xmin) * (adj - 1) / 2
  buf.y <- (locations$ymax - locations$ymin) * (adj - 1) / 2
  locations$xmin <- locations$xmin - buf.x
  locations$ymin <- locations$ymin - buf.y
  locations$xmax <- locations$xmax + buf.x
  locations$ymax <- locations$ymax + buf.y
  
  # The previous step reduces bias in quadrat size
  # at the cost of generating some overlapping quadrats.
  # This is mostly theoretical and won't affect in almost
  # any application, but it can be an issue in certain
  # non-standard cases. If this step causes problems, cancel these lines.
  # Note that this won't run if regular = TRUE, because in that 
  # case the quadrats will be non-overlapping by definition.
  if(!regular) 
  {
    locations <- make_quadrats_non_overlapping(locations, k = 8)
    rownames(locations) <- locations$quadrat_id
  }
  
  # Estimate a standard, regular, complete grid,
  # but possibly rectangular quadrats.
  if(regular)
  {
    # Estimate grid steps in x and y dimensions
    quadrat.side.x = median((locations$xmax - locations$xmin), na.rm = TRUE)
    quadrat.side.y = median((locations$ymax - locations$ymin), na.rm = TRUE)
    
    quadrat.side.x = (outer.xmax - outer.xmin) /
      round((outer.xmax - outer.xmin) / quadrat.side.x)
    
    quadrat.side.y = (outer.ymax - outer.ymin) /
      round((outer.ymax - outer.ymin) / quadrat.side.y)
    
    # Build a regular grid from scratch.
    # For coherence with other functions: small leftovers preferred.
    xmins <- seq(from = outer.xmin, to = outer.xmax, by = quadrat.side.x)
    xmins <- xmins[xmins < outer.xmax]
    
    ymins <- seq(from = outer.ymin, to = outer.ymax, by = quadrat.side.y)
    ymins <- ymins[ymins < outer.ymax]
    
    grid <- expand.grid(xmin = xmins, ymin = ymins)
    grid$xmax <- pmin(outer.xmax, grid$xmin + quadrat.side.x)
    grid$ymax <- pmin(outer.ymax, grid$ymin + quadrat.side.y)
    
    # Then, compare the regular grid with the input data,
    # to figure out which ID corresponds to each grid cell.
    look.at <- which(!is.na(outer.x) & !is.na(outer.y) & !is.na(quadrat.id))
    
    ncol.grid <- length(xmins)
    nrow.grid <- length(ymins)
    
    col <- floor((outer.x[look.at] - outer.xmin) / quadrat.side.x) + 1
    row <- floor((outer.y[look.at] - outer.ymin) / quadrat.side.y) + 1
    
    # Include points that fall exactly on the ultimate right/top boundary.
    col[outer.x[look.at] == outer.xmax] <- ncol.grid
    row[outer.y[look.at] == outer.ymax] <- nrow.grid
    
    inside <-
      col >= 1 & col <= ncol.grid &
      row >= 1 & row <= nrow.grid
    
    grid.i <- col[inside] + (row[inside] - 1) * ncol.grid
    
    pairs <- data.frame(
      grid.i = grid.i,
      quadrat_id = as.character(quadrat.id[look.at][inside])
    )
    
    # Count observed quadrat IDs within each generated grid cell.
    counts <- aggregate(
      x = list(n = rep(1, nrow(pairs))),
      by = list(grid.i = pairs$grid.i, quadrat_id = pairs$quadrat_id),
      FUN = length
    )
    
    # For each grid cell, keep the most frequent observed quadrat ID.
    counts <- counts[order(counts$grid.i, -counts$n), ]
    best <- counts[!duplicated(counts$grid.i), ]
    
    most.likely.quadrat.ids <- rep(NA_character_, nrow(grid))
    most.likely.quadrat.ids[best$grid.i] <- best$quadrat_id
    
    if(any(is.na(most.likely.quadrat.ids)))
      warning("some grid cells had no observed points with non-NA quadrat IDs")
    
    # Add likely quadrat ID to the grid:
    grid <- data.frame(quadrat_id = most.likely.quadrat.ids, grid)
  }
  
  
  # Organize output, order in consistent way:
  if(regular)  out <- grid
  if(!regular) out <- locations
  
  o <- order(out$xmin, out$ymin)
  out <- out[o, , drop = FALSE]
  
  if(!anyNA(out$quadrat_id) && !anyDuplicated(out$quadrat_id))
    rownames(out) <- out$quadrat_id
  
  out
}



#########################################################################
####### FUNCTIONS THAT HANDLE SETS OF POINTS, GIVEN QUADRAT LOCATIONS
#########################################################################

# This function takes outer (~global) coordinates and a set of quadrats,
# and returns memberships into those quadrats + inner (~local) coordinates.

# WARNING: it is flexible by design, so it would assign points to
# random, unequal, overlapping, etc. quadrats, and some points may
# not be assigned to any quadrat, if that's the input!

assign_points_to_quadrats <- function(
    outer.x = NA,
    outer.y = NA,
    quadrat.locations = NA,
    boundary.rule = "left_bottom"
)
{
  
  something_missing <- function(x) {
    is.null(x) || length(x) == 0 || all(is.na(x))
  }
  
  # Check for boundary rule
  if(length(boundary.rule) != 1)
  {
    boundary.rule = "left_bottom"
    warning("'boundary.rule' assumed 'left_bottom' by default; only valid alternative is 'right_top'")
  }
  
  if(!boundary.rule %in% c("left_bottom", "right_top"))
  {
    boundary.rule = "left_bottom"
    warning("'boundary.rule' assumed 'left_bottom' by default; only valid alternative is 'right_top'")
  }
  
  # Checks: point coordinates
  if (missing(outer.x) || is.null(outer.x) || length(outer.x) == 0)
    stop("outer.x must be supplied")
  
  if (missing(outer.y) || is.null(outer.y) || length(outer.y) == 0)
    stop("outer.y must be supplied")
  
  if (length(outer.x) != length(outer.y))
    stop("outer.x and outer.y must have the same length")
  
  # All-NA vectors are allowed, even if R stored them as logical.
  if (all(is.na(outer.x)))
    outer.x <- as.numeric(outer.x)
  
  if (all(is.na(outer.y)))
    outer.y <- as.numeric(outer.y)
  
  if (!is.numeric(outer.x) || !is.numeric(outer.y))
    stop("outer.x and outer.y must be numeric")
  
  # Checks: quadrat locations
  if (something_missing(quadrat.locations))
    stop("quadrat.locations must be supplied")
  
  needed <- c("quadrat_id", "xmin", "ymin", "xmax", "ymax")
  missing <- setdiff(needed, names(quadrat.locations))
  
  if (length(missing) > 0)
    stop("quadrat.locations is missing: ", paste(missing, collapse = ", "))
  
  if (nrow(quadrat.locations) == 0)
    stop("quadrat.locations has no rows")
  
  if (any(is.na(quadrat.locations$quadrat_id)))
    stop("quadrat_id cannot contain NA")
  
  if (anyDuplicated(quadrat.locations$quadrat_id))
    stop("quadrat_id must be unique in quadrat.locations")
  
  coord.cols <- c("xmin", "ymin", "xmax", "ymax")
  
  if (!all(sapply(quadrat.locations[coord.cols], is.numeric)))
    stop("xmin, ymin, xmax, and ymax must be numeric")
  
  if (any(is.na(quadrat.locations[coord.cols])))
    stop("quadrat locations cannot contain NA")
  
  if (any(quadrat.locations$xmax <= quadrat.locations$xmin))
    stop("all quadrats must have xmax > xmin")
  
  if (any(quadrat.locations$ymax <= quadrat.locations$ymin))
    stop("all quadrats must have ymax > ymin")
  
  quadrat.locations$quadrat_id <- as.character(quadrat.locations$quadrat_id)
  
  # Original input positions
  original.i <- seq_along(outer.x)
  
  # Work only with complete coordinates.
  # These working copies are used only to determine quadrat assignment.
  keep <- !is.na(outer.x) & !is.na(outer.y)
  i <- original.i[keep]
  x <- outer.x[keep]
  y <- outer.y[keep]
  
  if (!any(keep))
    warning("no x-y coordinates were supplied: all points will receive NA quadrat assignments")
  
  # Outer extent represented by the quadrat-location table
  outer.xmin <- min(quadrat.locations$xmin)
  outer.ymin <- min(quadrat.locations$ymin)
  outer.xmax <- max(quadrat.locations$xmax)
  outer.ymax <- max(quadrat.locations$ymax)
  
  # Points exactly on the ultimate outer boundary may be excluded by the
  # half-open boundary rule. Move only those working coordinates slightly
  # inside, only for assignment. Final coordinates are taken from the
  # original input vectors: outer.x[assigned$i] and outer.y[assigned$i].
  quadrat.side.x <- min(quadrat.locations$xmax - quadrat.locations$xmin)
  quadrat.side.y <- min(quadrat.locations$ymax - quadrat.locations$ymin)
  small.adjustment.x <- quadrat.side.x / 1000
  small.adjustment.y <- quadrat.side.y / 1000
  
  if(boundary.rule == "left_bottom")
  {
    x[x == outer.xmax] <- outer.xmax - small.adjustment.x
    y[y == outer.ymax] <- outer.ymax - small.adjustment.y
  }
  
  if(boundary.rule == "right_top")
  {
    x[x == outer.xmin] <- outer.xmin + small.adjustment.x
    y[y == outer.ymin] <- outer.ymin + small.adjustment.y
  }
  
  # Empty assignment table. Coordinates are added after quadrat assignment.
  assigned <- data.frame(
    i = integer(),
    quadrat_id = character(),
    stringsAsFactors = FALSE
  )
  
  pieces <- list(assigned)
  
  # ---------------------------------------------------------------------------
  # Assignment
  # ---------------------------------------------------------------------------
  
  for (q in seq_len(nrow(quadrat.locations)))
  {
    quad <- quadrat.locations[q, ]
    
    if (boundary.rule == "left_bottom")
      inside <- x >= quad$xmin & x < quad$xmax & y >= quad$ymin & y < quad$ymax
    
    if (boundary.rule == "right_top")
      inside <- x > quad$xmin & x <= quad$xmax & y > quad$ymin & y <= quad$ymax
    
    j <- which(inside)
    
    if (length(j) > 0)
    {
      pieces[[length(pieces) + 1]] <- data.frame(
        i = i[j],
        quadrat_id = quad$quadrat_id,
        stringsAsFactors = FALSE
      )
    }
  }
  
  assigned <- do.call(rbind, pieces)
  
  # ---------------------------------------------------------------------------
  # Keep all unassigned original points as NA quadrat assignments
  # ---------------------------------------------------------------------------
  
  assigned.i <- unique(assigned$i)
  outside.i <- setdiff(original.i, assigned.i)
  
  if (length(outside.i) > 0)
  {
    outside <- data.frame(
      i = outside.i,
      quadrat_id = NA_character_,
      stringsAsFactors = FALSE
    )
    
    assigned <- rbind(assigned, outside)
  }
  
  # ---------------------------------------------------------------------------
  # Add coordinates and inner coordinates
  # ---------------------------------------------------------------------------
  
  assigned <- assigned[order(assigned$i), , drop = FALSE]
  rownames(assigned) <- NULL
  
  assigned$outer.x <- outer.x[assigned$i]
  assigned$outer.y <- outer.y[assigned$i]
  
  q <- match(assigned$quadrat_id, quadrat.locations$quadrat_id)
  
  assigned$inner.x <- assigned$outer.x - quadrat.locations$xmin[q]
  assigned$inner.y <- assigned$outer.y - quadrat.locations$ymin[q]
  
  assigned <- assigned[, c("i", "outer.x", "outer.y", "quadrat_id", "inner.x", "inner.y")]
  
  # ---------------------------------------------------------------------------
  # Warnings
  # ---------------------------------------------------------------------------
  
  if (length(outside.i) > 0)
  {
    warning(
      length(outside.i),
      " point(s) were not assigned to any quadrat; ",
      "quadrat_id, inner.x, and inner.y were set to NA"
    )
  }
  
  n.per.point <- table(assigned$i[!is.na(assigned$quadrat_id)])
  multiple.i <- names(n.per.point)[n.per.point > 1]
  
  if (length(multiple.i) > 0)
  {
    warning(
      length(multiple.i),
      " point(s) were assigned to more than one quadrat; ",
      "this is expected if quadrats overlap"
    )
  }
  
  assigned
}



# This function does the opposite: takes membership into quadrats, inner coordinates
# within those quadrats, and quadrat location, and compiles all those coordinates
# into outer coordinates, where all points are located relative to the coordinate
# system that describes the quadrat locations. 
# For example: to compile full-plot datasets from quadrat-level datasets.





#########################################################################
####### FUNCTIONS THAT HANDLE QUADRATS: 1-to-n and n-to-1
#########################################################################


