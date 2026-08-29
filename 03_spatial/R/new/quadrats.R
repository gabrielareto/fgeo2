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


# THE CORE FUNCTION TO SPLIT INTO QUADRATS (FAST, no IDs complications)
make_quadrats_fast <- function(
    outer.xmin = NA,
    outer.ymin = NA,
    outer.xmax = NA,
    outer.ymax = NA,
    quadrat.side.x = NA,
    quadrat.side.y = NA
)
{
  something_missing <- function(x) {
    is.null(x) || length(x) == 0 || all(is.na(x))
  }
  
  # Required dimensions
  if(something_missing(outer.xmax) || something_missing(outer.ymax))
    stop("specify dimensions through outer.xmax and outer.ymax")
  
  if(something_missing(quadrat.side.x) && something_missing(quadrat.side.y))
    stop("specify quadrat size through quadrat.side.x and/or quadrat.side.y")
  
  # Defaults
  if(something_missing(outer.xmin))
  {
    outer.xmin <- 0
    warning("outer.xmin = 0 assumed by default")
  }
  
  if(something_missing(outer.ymin))
  {
    outer.ymin <- 0
    warning("outer.ymin = 0 assumed by default")
  }
  
  if(something_missing(quadrat.side.x))
  {
    quadrat.side.x <- quadrat.side.y
    warning("quadrat.side.x assumed equal to quadrat.side.y")
  }
  
  if(something_missing(quadrat.side.y))
  {
    quadrat.side.y <- quadrat.side.x
    warning("quadrat.side.y assumed equal to quadrat.side.x")
  }
  
  # Basic checks
  if(length(outer.xmin) != 1 || !is.numeric(outer.xmin) || !is.finite(outer.xmin))
    stop("outer.xmin must be one finite numeric value")
  
  if(length(outer.ymin) != 1 || !is.numeric(outer.ymin) || !is.finite(outer.ymin))
    stop("outer.ymin must be one finite numeric value")
  
  if(length(outer.xmax) != 1 || !is.numeric(outer.xmax) || !is.finite(outer.xmax))
    stop("outer.xmax must be one finite numeric value")
  
  if(length(outer.ymax) != 1 || !is.numeric(outer.ymax) || !is.finite(outer.ymax))
    stop("outer.ymax must be one finite numeric value")
  
  if(length(quadrat.side.x) != 1 || !is.numeric(quadrat.side.x) || !is.finite(quadrat.side.x))
    stop("quadrat.side.x must be one finite numeric value")
  
  if(length(quadrat.side.y) != 1 || !is.numeric(quadrat.side.y) || !is.finite(quadrat.side.y))
    stop("quadrat.side.y must be one finite numeric value")
  
  if(outer.xmax <= outer.xmin)
    stop("outer.xmax must be greater than outer.xmin")
  
  if(outer.ymax <= outer.ymin)
    stop("outer.ymax must be greater than outer.ymin")
  
  if(quadrat.side.x <= 0 || quadrat.side.y <= 0)
    stop("quadrat sides must be greater than zero")
  
  # Split the outer extent. Small leftovers remain at the right and top edges.
  xmins <- seq(from = outer.xmin, to = outer.xmax, by = quadrat.side.x)
  xmins <- xmins[xmins < outer.xmax]
  
  ymins <- seq(from = outer.ymin, to = outer.ymax, by = quadrat.side.y)
  ymins <- ymins[ymins < outer.ymax]
  
  grid <- expand.grid(
    xmin = xmins,
    ymin = ymins,
    KEEP.OUT.ATTRS = FALSE
  )
  
  grid$xmax <- pmin(outer.xmax, grid$xmin + quadrat.side.x)
  grid$ymax <- pmin(outer.ymax, grid$ymin + quadrat.side.y)
  
  
  # Sort by columns: move along y first, then to the next x.
  # This matches the historical CTFS sequential quadrat indexes.
  o <- order(grid$xmin, grid$ymin)
  grid <- grid[o, , drop = FALSE]
  rownames(grid) <- NULL
  
  # Add very simple IDs
  grid <- data.frame(
    quadrat_id = as.character(seq_len(nrow(grid))),
    grid,
    stringsAsFactors = FALSE
  )
  
  grid
}


# USING PLOT & QUADRAT DIMENSIONS, AND DESCRIPTORS OF CONVENTION
make_quadrats <- function(
  # Dimensions of the large unit
  outer.xmin = NA,
  outer.ymin = NA,
  outer.xmax = NA,
  outer.ymax = NA,
  
  # Dimensions of the smaller sub-units
  quadrat.side.x = NA,
  quadrat.side.y = NA,
  
  # How the ID is built
  built = "sequence",
  example.part.for.col = "1",
  example.part.for.row = "1",
  example.part.for.sequence = "1",
  prefix = "",
  suffix = "",
  separator = "",
  
  # Explicit sequences, if supplied
  seq.of.row.parts.along.one.col = NA,
  seq.of.col.parts.along.one.row = NA,
  entire.sequence = NA,
  
  # How the sequences will run
  start.col.at = NA,
  start.row.at = start.col.at,
  start.sequence.at = 1,
  along.columns.first.then.next.column = TRUE,
  along.rows.first.then.next.row = FALSE,
  changing.directions = FALSE
)
{
  something_missing <- function(x) {
    is.null(x) || length(x) == 0 || all(is.na(x))
  }
  
  # Check ID convention
  if(length(built) != 1 ||
     is.na(built) ||
     !built %in% c("col+row", "row+col", "sequence"))
    stop("built must be 'col+row', 'row+col', or 'sequence'")
  
  # Resolve quadrat sides here because they are also used below.
  if(something_missing(quadrat.side.x) && something_missing(quadrat.side.y))
    stop("specify quadrat size through quadrat.side.x and/or quadrat.side.y")
  
  if(something_missing(quadrat.side.x))
  {
    quadrat.side.x <- quadrat.side.y
    warning("quadrat.side.x assumed equal to quadrat.side.y")
  }
  
  if(something_missing(quadrat.side.y))
  {
    quadrat.side.y <- quadrat.side.x
    warning("quadrat.side.y assumed equal to quadrat.side.x")
  }
  
  # Create the grid geometry.
  grid <- make_quadrats_fast(
    outer.xmin = outer.xmin,
    outer.ymin = outer.ymin,
    outer.xmax = outer.xmax,
    outer.ymax = outer.ymax,
    quadrat.side.x = quadrat.side.x,
    quadrat.side.y = quadrat.side.y
  )
  
  # The primitive IDs are replaced below.
  grid$quadrat_id <- NULL
  
  # Column and row positions within the grid.
  grid$col <- match(grid$xmin, sort(unique(grid$xmin)))
  grid$row <- match(grid$ymin, sort(unique(grid$ymin)))
  
  # Defaults for starting values.
  if(built %in% c("col+row", "row+col"))
  {
    if(something_missing(start.col.at))
    {
      start.col.at <- example.part.for.col
      warning("start.col.at taken from example.part.for.col by default")
    }
    
    if(something_missing(start.row.at))
    {
      start.row.at <- example.part.for.row
      warning("start.row.at taken from example.part.for.row by default")
    }
  }
  
  # Reorder before assigning one sequence to all quadrats.
  if(built == "sequence")
  {
    if(length(along.columns.first.then.next.column) != 1 ||
       is.na(along.columns.first.then.next.column) ||
       !is.logical(along.columns.first.then.next.column))
      stop("along.columns.first.then.next.column must be TRUE or FALSE")
    
    if(length(along.rows.first.then.next.row) != 1 ||
       is.na(along.rows.first.then.next.row) ||
       !is.logical(along.rows.first.then.next.row))
      stop("along.rows.first.then.next.row must be TRUE or FALSE")
    
    if(along.columns.first.then.next.column == along.rows.first.then.next.row)
      stop("clarify whether labelling is along columns first, or along rows first")
    
    if(length(changing.directions) != 1 ||
       is.na(changing.directions) ||
       !is.logical(changing.directions))
      stop("changing.directions must be TRUE or FALSE")
    
    if(along.columns.first.then.next.column)
      o <- order(grid$col, grid$row)
    
    if(along.rows.first.then.next.row)
      o <- order(grid$row, grid$col)
    
    grid <- grid[o, , drop = FALSE]
    
    if(changing.directions)
    {
      if(along.columns.first.then.next.column && max(grid$col) >= 2)
      {
        for(i in seq(from = 2, to = max(grid$col), by = 2))
        {
          chunk <- which(grid$col == i)
          grid[chunk, ] <- grid[rev(chunk), ]
        }
      }
      
      if(along.rows.first.then.next.row && max(grid$row) >= 2)
      {
        for(i in seq(from = 2, to = max(grid$row), by = 2))
        {
          chunk <- which(grid$row == i)
          grid[chunk, ] <- grid[rev(chunk), ]
        }
      }
    }
  }
  
  # Identify whether column and row parts use numbers or letters.
  if(built %in% c("col+row", "row+col"))
  {
    use.numbers.for.col1 <- !is.na(suppressWarnings(as.numeric(example.part.for.col)))
    use.numbers.for.col2 <- !is.na(suppressWarnings(as.numeric(start.col.at)))
    
    if(use.numbers.for.col1 != use.numbers.for.col2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for 'col' part")
    
    use.numbers.for.row1 <- !is.na(suppressWarnings(as.numeric(example.part.for.row)))
    use.numbers.for.row2 <- !is.na(suppressWarnings(as.numeric(start.row.at)))
    
    if(use.numbers.for.row1 != use.numbers.for.row2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for 'row' part")
    
    use.numbers.for.col <- use.numbers.for.col1
    use.numbers.for.row <- use.numbers.for.row1
    
    use.letters.for.col1 <- grepl("([a-zA-Z])", example.part.for.col)
    use.letters.for.col2 <- grepl("([a-zA-Z])", start.col.at)
    
    if(use.letters.for.col1 != use.letters.for.col2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for 'col' part")
    
    use.letters.for.row1 <- grepl("([a-zA-Z])", example.part.for.row)
    use.letters.for.row2 <- grepl("([a-zA-Z])", start.row.at)
    
    if(use.letters.for.row1 != use.letters.for.row2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for 'row' part")
    
    use.letters.for.col <- use.letters.for.col1
    use.letters.for.row <- use.letters.for.row1
    
    if(use.numbers.for.col == use.letters.for.col)
      stop("clarify between numbers and letters for 'col' part")
    
    if(use.numbers.for.row == use.letters.for.row)
      stop("clarify between numbers and letters for 'row' part")
  }
  
  # Identify whether a single sequence uses numbers or letters.
  if(built == "sequence")
  {
    use.numbers.for.sequence1 <- !is.na(suppressWarnings(as.numeric(example.part.for.sequence)))
    use.numbers.for.sequence2 <- !is.na(suppressWarnings(as.numeric(start.sequence.at)))
    
    if(use.numbers.for.sequence1 != use.numbers.for.sequence2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for the sequence")
    
    use.numbers.for.sequence <- use.numbers.for.sequence1
    
    use.letters.for.sequence1 <- grepl("([a-zA-Z])", example.part.for.sequence)
    use.letters.for.sequence2 <- grepl("([a-zA-Z])", start.sequence.at)
    
    if(use.letters.for.sequence1 != use.letters.for.sequence2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for the sequence")
    
    use.letters.for.sequence <- use.letters.for.sequence1
    
    if(use.numbers.for.sequence == use.letters.for.sequence)
      stop("clarify between numbers and letters for the sequence")
  }
  
  # Build column and row sequences.
  if(built %in% c("col+row", "row+col"))
  {
    if(use.numbers.for.col)
    {
      if(as.numeric(example.part.for.col) == 0 && as.numeric(start.col.at) > 0)
        warning(
          "start.col.at = ", start.col.at,
          " takes priority, although example suggests start = 0"
        )
      
      basic.col.seq <- seq(
        from = as.numeric(start.col.at),
        length.out = length(unique(grid$col)),
        by = 1
      )
      
      add <- pmax(0, nchar(example.part.for.col) - nchar(basic.col.seq))
      basic.col.seq <- paste0(strrep("0", add), basic.col.seq)
    }
    
    if(use.numbers.for.row)
    {
      if(as.numeric(example.part.for.row) == 0 && as.numeric(start.row.at) > 0)
        warning(
          "start.row.at = ", start.row.at,
          " takes priority, although example suggests start = 0"
        )
      
      basic.row.seq <- seq(
        from = as.numeric(start.row.at),
        length.out = length(unique(grid$row)),
        by = 1
      )
      
      add <- pmax(0, nchar(example.part.for.row) - nchar(basic.row.seq))
      basic.row.seq <- paste0(strrep("0", add), basic.row.seq)
    }
    
    if(use.letters.for.col)
    {
      lower <- substr(example.part.for.col, 1, 1) %in% letters
      upper <- substr(example.part.for.col, 1, 1) %in% LETTERS
      
      if(!lower && !upper)
        stop("the letter in 'example.part.for.col' cannot be recognized")
      
      n <- length(unique(grid$col))
      
      if(nchar(example.part.for.col) == 1)
        basic.col.seq <- letter_seq(n)
      
      if(nchar(example.part.for.col) > 1)
      {
        k <- nchar(example.part.for.col)
        
        if(n > 26^k)
          stop("the column labels need more characters")
        
        basic.col.seq <- letter_seq_fixed(n = n, k = k)
      }
      
      if(upper)
        basic.col.seq <- toupper(basic.col.seq)
    }
    
    if(use.letters.for.row)
    {
      lower <- substr(example.part.for.row, 1, 1) %in% letters
      upper <- substr(example.part.for.row, 1, 1) %in% LETTERS
      
      if(!lower && !upper)
        stop("the letter in 'example.part.for.row' cannot be recognized")
      
      n <- length(unique(grid$row))
      
      if(nchar(example.part.for.row) == 1)
        basic.row.seq <- letter_seq(n)
      
      if(nchar(example.part.for.row) > 1)
      {
        k <- nchar(example.part.for.row)
        
        if(n > 26^k)
          stop("the row labels need more characters")
        
        basic.row.seq <- letter_seq_fixed(n = n, k = k)
      }
      
      if(upper)
        basic.row.seq <- toupper(basic.row.seq)
    }
    
    grid$col.part <- basic.col.seq[grid$col]
    grid$row.part <- basic.row.seq[grid$row]
  }
  
  # Build one sequence for all quadrats.
  if(built == "sequence")
  {
    if(use.numbers.for.sequence)
    {
      s <- seq(
        from = as.numeric(start.sequence.at),
        length.out = nrow(grid),
        by = 1
      )
      
      add <- pmax(0, nchar(example.part.for.sequence) - nchar(s))
      grid$entire <- paste0(strrep("0", add), s)
    }
    
    if(use.letters.for.sequence)
    {
      lower <- substr(example.part.for.sequence, 1, 1) %in% letters
      upper <- substr(example.part.for.sequence, 1, 1) %in% LETTERS
      
      if(!lower && !upper)
        stop("the letter in 'example.part.for.sequence' cannot be recognized")
      
      n <- nrow(grid)
      
      if(nchar(example.part.for.sequence) == 1)
        s <- letter_seq(n)
      
      if(nchar(example.part.for.sequence) > 1)
      {
        k <- nchar(example.part.for.sequence)
        
        if(n > 26^k)
          stop("the sequence labels need more characters")
        
        s <- letter_seq_fixed(n = n, k = k)
      }
      
      if(upper)
        s <- toupper(s)
      
      grid$entire <- s
    }
  }
  
  # Explicit column and row sequences take priority.
  if(!something_missing(seq.of.col.parts.along.one.row))
  {
    if(length(seq.of.col.parts.along.one.row) != length(unique(grid$col)))
      stop(
        "the sequence provided for 'col' parts does not match number of columns -- ",
        "perhaps it's our 'leftover' convention"
      )
    
    grid$col.part <- seq.of.col.parts.along.one.row[grid$col]
  }
  
  if(!something_missing(seq.of.row.parts.along.one.col))
  {
    if(length(seq.of.row.parts.along.one.col) != length(unique(grid$row)))
      stop(
        "the sequence provided for 'row' parts does not match number of rows -- ",
        "perhaps it's our 'leftover' convention"
      )
    
    grid$row.part <- seq.of.row.parts.along.one.col[grid$row]
  }
  
  # Combine column and row parts.
  if(built == "col+row")
    grid$entire <- paste0(grid$col.part, separator, grid$row.part)
  
  if(built == "row+col")
    grid$entire <- paste0(grid$row.part, separator, grid$col.part)
  
  # Add prefix and suffix.
  if(!something_missing(prefix))
    grid$entire <- paste0(prefix, grid$entire)
  
  if(!something_missing(suffix))
    grid$entire <- paste0(grid$entire, suffix)
  
  # A complete supplied sequence has final priority.
  if(!something_missing(entire.sequence))
  {
    if(length(entire.sequence) != nrow(grid))
      stop("entire.sequence must have one value per quadrat")
    
    grid$entire <- as.character(entire.sequence)
  }
  
  # Warn about smaller leftover quadrats.
  widths <- grid$xmax - grid$xmin
  heights <- grid$ymax - grid$ymin
  
  if(!isTRUE(all.equal(widths, rep(quadrat.side.x, length(widths)))))
    warning("the quadrats don't fit perfectly along x: quadrats at the end are smaller")
  
  if(!isTRUE(all.equal(heights, rep(quadrat.side.y, length(heights)))))
    warning("the quadrats don't fit perfectly along y: quadrats at the end are smaller")
  
  # Trim, rename, and sort output.
  out <- grid[, c("entire", "xmin", "ymin", "xmax", "ymax"), drop = FALSE]
  names(out)[1] <- "quadrat_id"
  out$quadrat_id <- as.character(out$quadrat_id)
  
  if(anyNA(out$quadrat_id))
    stop("generated quadrat IDs contain NA")
  
  if(anyDuplicated(out$quadrat_id))
    stop("generated quadrat IDs are not unique")
  
  o <- order(out$xmin, out$ymin)
  out <- out[o, , drop = FALSE]
  
  rownames(out) <- out$quadrat_id
  
  out
}



# USING A SET OF POINTS WITH COORDINATES
get_quadrat_locations_from_outer_coordinates <- function(
    quadrat.id = NA,
    outer.x = NA,
    outer.y = NA,
    outer.xmin = NA,
    outer.ymin = NA,
    outer.xmax = NA,
    outer.ymax = NA,
    regular = TRUE
)
{
  something_missing <- function(x) {
    is.null(x) || length(x) == 0 || all(is.na(x))
  }
  
  # Checks
  if(something_missing(quadrat.id))
    stop("quadrat.id must be supplied")
  
  if(something_missing(outer.x))
    stop("outer.x must be supplied")
  
  if(something_missing(outer.y))
    stop("outer.y must be supplied")
  
  if(length(quadrat.id) != length(outer.x) ||
     length(quadrat.id) != length(outer.y))
    stop("quadrat.id, outer.x, and outer.y must have the same length")
  
  if(!is.numeric(outer.x) || !is.numeric(outer.y))
    stop("outer.x and outer.y must be numeric")
  
  if(any(is.infinite(outer.x)) || any(is.infinite(outer.y)))
    stop("outer.x and outer.y cannot contain Inf or -Inf")
  
  if(length(regular) != 1 || is.na(regular) || !is.logical(regular))
    stop("regular must be TRUE or FALSE")
  
  complete <- !is.na(quadrat.id) & !is.na(outer.x) & !is.na(outer.y)
  
  if(!any(complete))
    stop("no complete quadrat.id, outer.x, and outer.y observations were supplied")
  
  # Defaults for outer dimensions
  inferred.bounds <- c(
    outer.xmin = something_missing(outer.xmin),
    outer.ymin = something_missing(outer.ymin),
    outer.xmax = something_missing(outer.xmax),
    outer.ymax = something_missing(outer.ymax)
  )
  
  if(inferred.bounds[["outer.xmin"]])
    outer.xmin <- floor(min(outer.x, na.rm = TRUE))
  
  if(inferred.bounds[["outer.ymin"]])
    outer.ymin <- floor(min(outer.y, na.rm = TRUE))
  
  if(inferred.bounds[["outer.xmax"]])
    outer.xmax <- ceiling(max(outer.x, na.rm = TRUE))
  
  if(inferred.bounds[["outer.ymax"]])
    outer.ymax <- ceiling(max(outer.y, na.rm = TRUE))
  
  if(any(inferred.bounds))
    warning(
      "outer bounds inferred from observed coordinates: ",
      paste(names(inferred.bounds)[inferred.bounds], collapse = ", ")
    )
  
  if(length(outer.xmin) != 1 || !is.numeric(outer.xmin) || !is.finite(outer.xmin))
    stop("outer.xmin must be one finite numeric value")
  
  if(length(outer.ymin) != 1 || !is.numeric(outer.ymin) || !is.finite(outer.ymin))
    stop("outer.ymin must be one finite numeric value")
  
  if(length(outer.xmax) != 1 || !is.numeric(outer.xmax) || !is.finite(outer.xmax))
    stop("outer.xmax must be one finite numeric value")
  
  if(length(outer.ymax) != 1 || !is.numeric(outer.ymax) || !is.finite(outer.ymax))
    stop("outer.ymax must be one finite numeric value")
  
  if(outer.xmax <= outer.xmin)
    stop("outer.xmax must be greater than outer.xmin")
  
  if(outer.ymax <= outer.ymin)
    stop("outer.ymax must be greater than outer.ymin")
  
  # Use complete observations to estimate quadrat locations.
  dat <- data.frame(
    quadrat_id = as.character(quadrat.id[complete]),
    outer.x = outer.x[complete],
    outer.y = outer.y[complete],
    stringsAsFactors = FALSE
  )
  
  xmins <- aggregate(outer.x ~ quadrat_id, data = dat, FUN = min)
  ymins <- aggregate(outer.y ~ quadrat_id, data = dat, FUN = min)
  xmaxs <- aggregate(outer.x ~ quadrat_id, data = dat, FUN = max)
  ymaxs <- aggregate(outer.y ~ quadrat_id, data = dat, FUN = max)
  
  names(xmins)[2] <- "xmin"
  names(ymins)[2] <- "ymin"
  names(xmaxs)[2] <- "xmax"
  names(ymaxs)[2] <- "ymax"
  
  locations <- Reduce(
    function(x, y) merge(x, y, by = "quadrat_id", all = TRUE, sort = FALSE),
    list(xmins, ymins, xmaxs, ymaxs)
  )
  
  # Observed ranges underestimate true quadrat dimensions.
  # Increase them under the assumption of randomly distributed points.
  nq <- table(dat$quadrat_id)
  
  adj <- (nq[locations$quadrat_id] + 1) /
    (nq[locations$quadrat_id] - 1)
  
  adj[adj > 2] <- 2
  
  buf.x <- (locations$xmax - locations$xmin) * (adj - 1) / 2
  buf.y <- (locations$ymax - locations$ymin) * (adj - 1) / 2
  
  locations$xmin <- locations$xmin - buf.x
  locations$ymin <- locations$ymin - buf.y
  locations$xmax <- locations$xmax + buf.x
  locations$ymax <- locations$ymax + buf.y
  
  # For non-regular cases, remove local overlaps caused by the adjustment.
  if(!regular)
  {
    bad <- locations$xmax <= locations$xmin |
      locations$ymax <= locations$ymin
    
    if(any(bad))
      stop("could not infer positive dimensions for all quadrats")
    
    locations <- make_quadrats_non_overlapping(locations, k = 8)
    
    bad <- locations$xmax <= locations$xmin |
      locations$ymax <= locations$ymin
    
    if(any(bad))
      stop("some quadrats ended with invalid dimensions after overlap adjustment")
  }
  
  # Estimate a standard regular and complete grid.
  if(regular)
  {
    quadrat.side.x <- median(locations$xmax - locations$xmin, na.rm = TRUE)
    quadrat.side.y <- median(locations$ymax - locations$ymin, na.rm = TRUE)
    
    if(!is.finite(quadrat.side.x) || quadrat.side.x <= 0)
      stop("could not infer a positive quadrat side along x")
    
    if(!is.finite(quadrat.side.y) || quadrat.side.y <= 0)
      stop("could not infer a positive quadrat side along y")
    
    ncol.grid <- max(
      1,
      round((outer.xmax - outer.xmin) / quadrat.side.x)
    )
    
    nrow.grid <- max(
      1,
      round((outer.ymax - outer.ymin) / quadrat.side.y)
    )
    
    quadrat.side.x <- (outer.xmax - outer.xmin) / ncol.grid
    quadrat.side.y <- (outer.ymax - outer.ymin) / nrow.grid
    
    grid <- make_quadrats_fast(
      outer.xmin = outer.xmin,
      outer.ymin = outer.ymin,
      outer.xmax = outer.xmax,
      outer.ymax = outer.ymax,
      quadrat.side.x = quadrat.side.x,
      quadrat.side.y = quadrat.side.y
    )
    
    # Compare the regular grid with the input data to identify each grid cell.
    look.at <- which(complete)
    
    col <- floor((outer.x[look.at] - outer.xmin) / quadrat.side.x) + 1
    row <- floor((outer.y[look.at] - outer.ymin) / quadrat.side.y) + 1
    
    # Include points exactly on the ultimate right and top boundaries.
    col[outer.x[look.at] == outer.xmax] <- ncol.grid
    row[outer.y[look.at] == outer.ymax] <- nrow.grid
    
    inside <-
      col >= 1 & col <= ncol.grid &
      row >= 1 & row <= nrow.grid
    
    if(any(!inside))
      warning(
        sum(!inside),
        " complete observation(s) fell outside the outer extent and were ignored"
      )
    
    most.likely.quadrat.ids <- rep(NA_character_, nrow(grid))
    
    if(any(inside))
    {
      grid.i <- row[inside] + (col[inside] - 1) * nrow.grid
      
      pairs <- data.frame(
        grid.i = grid.i,
        quadrat_id = as.character(quadrat.id[look.at][inside]),
        stringsAsFactors = FALSE
      )
      
      # Count observed quadrat IDs within each generated grid cell.
      counts <- aggregate(
        x = list(n = rep(1, nrow(pairs))),
        by = list(
          grid.i = pairs$grid.i,
          quadrat_id = pairs$quadrat_id
        ),
        FUN = length
      )
      
      ids.per.cell <- table(counts$grid.i)
      conflicting.cells <- sum(ids.per.cell > 1)
      
      if(conflicting.cells > 0)
        warning(
          conflicting.cells,
          " grid cell(s) contained more than one quadrat ID; ",
          "the most frequent ID was retained"
        )
      
      # For each grid cell, keep the most frequent observed quadrat ID.
      counts <- counts[order(counts$grid.i, -counts$n), ]
      best <- counts[!duplicated(counts$grid.i), ]
      
      most.likely.quadrat.ids[best$grid.i] <- best$quadrat_id
    }
    
    if(anyNA(most.likely.quadrat.ids))
      warning("some grid cells had no observed points with non-NA quadrat IDs")
    
    grid$quadrat_id <- most.likely.quadrat.ids
  }
  
  # Organize output.
  if(regular)
    out <- grid
  
  if(!regular)
    out <- locations
  
  out <- out[, c("quadrat_id", "xmin", "ymin", "xmax", "ymax"), drop = FALSE]
  
  o <- order(out$xmin, out$ymin)
  out <- out[o, , drop = FALSE]
  
  if(anyDuplicated(out$quadrat_id[!is.na(out$quadrat_id)]))
    warning("some quadrat IDs occur in more than one inferred grid cell")
  
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
    point.id = NA,
    quadrat.locations = NA,
    boundary.rule = "left_bottom"
)
{
  
  something_missing <- function(x) {
    is.null(x) || length(x) == 0 || all(is.na(x))
  }
  
  # Check boundary rule.
  if(length(boundary.rule) != 1 ||
     is.na(boundary.rule) ||
     !boundary.rule %in% c("left_bottom", "right_top"))
    stop("boundary.rule must be 'left_bottom' or 'right_top'")
  
  # Check point coordinates.
  if(missing(outer.x) || is.null(outer.x) || length(outer.x) == 0)
    stop("outer.x must be supplied")
  
  if(missing(outer.y) || is.null(outer.y) || length(outer.y) == 0)
    stop("outer.y must be supplied")
  
  if(length(outer.x) != length(outer.y))
    stop("outer.x and outer.y must have the same length")
  
  # All-NA vectors are allowed, even if R stored them as logical.
  if(all(is.na(outer.x)))
    outer.x <- as.numeric(outer.x)
  
  if(all(is.na(outer.y)))
    outer.y <- as.numeric(outer.y)
  
  if(!is.numeric(outer.x) || !is.numeric(outer.y))
    stop("outer.x and outer.y must be numeric")
  
  if(any(is.infinite(outer.x)) || any(is.infinite(outer.y)))
    stop("outer.x and outer.y cannot contain Inf or -Inf")
  
  # Use sequential point IDs by default.
  if(something_missing(point.id))
    point.id <- as.character(seq_along(outer.x))
  
  if(length(point.id) != length(outer.x))
    stop("point.id, outer.x, and outer.y must have the same length")
  
  if(any(is.na(point.id)))
    stop("point.id cannot contain NA")
  
  if(anyDuplicated(point.id))
    stop("point.id must be unique")
  
  point.id <- as.character(point.id)
  
  # Check quadrat locations.
  if(something_missing(quadrat.locations))
    stop("quadrat.locations must be supplied")
  
  needed <- c("quadrat_id", "xmin", "ymin", "xmax", "ymax")
  missing <- setdiff(needed, names(quadrat.locations))
  
  if(length(missing) > 0)
    stop("quadrat.locations is missing: ", paste(missing, collapse = ", "))
  
  if(nrow(quadrat.locations) == 0)
    stop("quadrat.locations has no rows")
  
  if(any(is.na(quadrat.locations$quadrat_id)))
    stop("quadrat_id cannot contain NA")
  
  if(anyDuplicated(quadrat.locations$quadrat_id))
    stop("quadrat_id must be unique in quadrat.locations")
  
  coord.cols <- c("xmin", "ymin", "xmax", "ymax")
  
  if(!all(sapply(quadrat.locations[coord.cols], is.numeric)))
    stop("xmin, ymin, xmax, and ymax must be numeric")
  
  if(any(is.na(quadrat.locations[coord.cols])))
    stop("quadrat locations cannot contain NA")
  
  if(any(is.infinite(as.matrix(quadrat.locations[coord.cols]))))
    stop("quadrat locations cannot contain Inf or -Inf")
  
  if(any(quadrat.locations$xmax <= quadrat.locations$xmin))
    stop("all quadrats must have xmax > xmin")
  
  if(any(quadrat.locations$ymax <= quadrat.locations$ymin))
    stop("all quadrats must have ymax > ymin")
  
  quadrat.locations$quadrat_id <- as.character(quadrat.locations$quadrat_id)
  
  # Work only with complete coordinates.
  # x and y are used only to determine quadrat assignment.
  point.i <- which(!is.na(outer.x) & !is.na(outer.y))
  x <- outer.x[point.i]
  y <- outer.y[point.i]
  
  if(length(point.i) == 0)
    warning("no x-y coordinates were supplied: all points will receive NA quadrat assignments")
  
  # Outer extent represented by the quadrat-location table.
  outer.xmin <- min(quadrat.locations$xmin)
  outer.ymin <- min(quadrat.locations$ymin)
  outer.xmax <- max(quadrat.locations$xmax)
  outer.ymax <- max(quadrat.locations$ymax)
  
  # Points exactly on the ultimate excluded boundary are moved slightly
  # inward only to achieve complete assignment. Final coordinates always
  # come directly from the original outer.x and outer.y input vectors.
  small.adjustment.x <- min(quadrat.locations$xmax - quadrat.locations$xmin) / 1000
  small.adjustment.y <- min(quadrat.locations$ymax - quadrat.locations$ymin) / 1000
  
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
  
  # Empty assignment table.
  assigned <- data.frame(
    point.i = integer(),
    quadrat_id = character(),
    stringsAsFactors = FALSE
  )
  
  pieces <- list(assigned)
  
  # Assign points to quadrats.
  for(q in seq_len(nrow(quadrat.locations)))
  {
    quad <- quadrat.locations[q, ]
    
    if(boundary.rule == "left_bottom")
      inside <- x >= quad$xmin & x < quad$xmax &
        y >= quad$ymin & y < quad$ymax
    
    if(boundary.rule == "right_top")
      inside <- x > quad$xmin & x <= quad$xmax &
        y > quad$ymin & y <= quad$ymax
    
    j <- which(inside)
    
    if(length(j) > 0)
    {
      pieces[[length(pieces) + 1]] <- data.frame(
        point.i = point.i[j],
        quadrat_id = quad$quadrat_id,
        stringsAsFactors = FALSE
      )
    }
  }
  
  assigned <- do.call(rbind, pieces)
  
  # Keep unassigned points as NA quadrat assignments.
  assigned.i <- unique(assigned$point.i)
  outside.i <- setdiff(seq_along(outer.x), assigned.i)
  
  if(length(outside.i) > 0)
  {
    outside <- data.frame(
      point.i = outside.i,
      quadrat_id = NA_character_,
      stringsAsFactors = FALSE
    )
    
    assigned <- rbind(assigned, outside)
  }
  
  # Sort by original point order.
  assigned <- assigned[order(assigned$point.i), , drop = FALSE]
  rownames(assigned) <- NULL
  
  # Add IDs and coordinates from the original input.
  assigned$point_id <- point.id[assigned$point.i]
  assigned$outer.x <- outer.x[assigned$point.i]
  assigned$outer.y <- outer.y[assigned$point.i]
  
  q <- match(assigned$quadrat_id, quadrat.locations$quadrat_id)
  
  assigned$inner.x <- assigned$outer.x - quadrat.locations$xmin[q]
  assigned$inner.y <- assigned$outer.y - quadrat.locations$ymin[q]
  
  assigned <- assigned[, c(
    "point_id",
    "outer.x",
    "outer.y",
    "quadrat_id",
    "inner.x",
    "inner.y"
  )]
  
  # Warnings.
  if(length(outside.i) > 0)
  {
    warning(
      length(outside.i),
      " point(s) were not assigned to any quadrat; ",
      "quadrat_id, inner.x, and inner.y were set to NA"
    )
  }
  
  n.per.point <- table(assigned$point_id[!is.na(assigned$quadrat_id)])
  multiple.points <- names(n.per.point)[n.per.point > 1]
  
  if(length(multiple.points) > 0)
  {
    warning(
      length(multiple.points),
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
# Because overlapping quadrats are allowed, in general, it checks for consistency
# of the possibly multiple locations for the same point. That consistency is
# evaluated within a given tolerance. 

compile_points_from_quadrats <- function(
    inner.x = NA,
    inner.y = NA,
    point.id = NA,
    quadrat.id = NA,
    quadrat.locations = NA,
    tolerance = 1e-8
)
{
  
  something_missing <- function(x) {
    is.null(x) || length(x) == 0 || all(is.na(x))
  }
  
  # Checks: coordinates and quadrat membership
  if(missing(inner.x) || is.null(inner.x) || length(inner.x) == 0)
    stop("inner.x must be supplied")
  
  if(missing(inner.y) || is.null(inner.y) || length(inner.y) == 0)
    stop("inner.y must be supplied")
  
  if(missing(quadrat.id) || is.null(quadrat.id) || length(quadrat.id) == 0)
    stop("quadrat.id must be supplied")
  
  if(length(inner.x) != length(inner.y) ||
     length(inner.x) != length(quadrat.id))
    stop("inner.x, inner.y, and quadrat.id must have the same length")
  
  # All-NA vectors are allowed, even if R stored them as logical.
  if(all(is.na(inner.x)))
    inner.x <- as.numeric(inner.x)
  
  if(all(is.na(inner.y)))
    inner.y <- as.numeric(inner.y)
  
  if(!is.numeric(inner.x) || !is.numeric(inner.y))
    stop("inner.x and inner.y must be numeric")
  
  quadrat.id <- as.character(quadrat.id)
  
  # Use sequential point IDs by default.
  # Without explicit point IDs, each input row is treated as a separate point.
  if(something_missing(point.id))
    point.id <- as.character(seq_along(inner.x))
  
  if(length(point.id) != length(inner.x))
    stop("point.id, inner.x, inner.y, and quadrat.id must have the same length")
  
  if(any(is.na(point.id)))
    stop("point.id cannot contain NA")
  
  point.id <- as.character(point.id)
  
  # Check tolerance.
  if(length(tolerance) != 1 ||
     !is.numeric(tolerance) ||
     is.na(tolerance) ||
     tolerance <= 0)
    stop("tolerance must be one positive numeric value")
  
  # Checks: quadrat locations
  if(something_missing(quadrat.locations))
    stop("quadrat.locations must be supplied")
  
  needed <- c("quadrat_id", "xmin", "ymin", "xmax", "ymax")
  missing <- setdiff(needed, names(quadrat.locations))
  
  if(length(missing) > 0)
    stop("quadrat.locations is missing: ", paste(missing, collapse = ", "))
  
  if(nrow(quadrat.locations) == 0)
    stop("quadrat.locations has no rows")
  
  if(any(is.na(quadrat.locations$quadrat_id)))
    stop("quadrat_id cannot contain NA")
  
  if(anyDuplicated(quadrat.locations$quadrat_id))
    stop("quadrat_id must be unique in quadrat.locations")
  
  coord.cols <- c("xmin", "ymin", "xmax", "ymax")
  
  if(!all(sapply(quadrat.locations[coord.cols], is.numeric)))
    stop("xmin, ymin, xmax, and ymax must be numeric")
  
  if(any(is.na(quadrat.locations[coord.cols])))
    stop("quadrat locations cannot contain NA")
  
  if(any(quadrat.locations$xmax <= quadrat.locations$xmin))
    stop("all quadrats must have xmax > xmin")
  
  if(any(quadrat.locations$ymax <= quadrat.locations$ymin))
    stop("all quadrats must have ymax > ymin")
  
  quadrat.locations$quadrat_id <- as.character(quadrat.locations$quadrat_id)
  
  # Match each input row to its quadrat.
  q <- match(quadrat.id, quadrat.locations$quadrat_id)
  
  missing.quadrats <- unique(quadrat.id[!is.na(quadrat.id) & is.na(q)])
  
  if(length(missing.quadrats) > 0)
  {
    warning(
      length(missing.quadrats),
      " quadrat ID(s) were not found in quadrat.locations; ",
      "their outer coordinates will be NA"
    )
  }
  
  # Reconstruct outer coordinates for every [point x quadrat] combination.
  outer.x <- inner.x + quadrat.locations$xmin[q]
  outer.y <- inner.y + quadrat.locations$ymin[q]
  
  # Use only complete reconstructed locations.
  complete <- !is.na(outer.x) & !is.na(outer.y)
  
  reconstructed <- data.frame(
    point_id = point.id[complete],
    outer.x = outer.x[complete],
    outer.y = outer.y[complete],
    rounded.x = round(outer.x[complete] / tolerance),
    rounded.y = round(outer.y[complete] / tolerance),
    stringsAsFactors = FALSE
  )
  
  # Check whether repeated reconstructions of each point agree.
  rounded.locations <- unique(
    reconstructed[, c("point_id", "rounded.x", "rounded.y")]
  )
  
  n.locations <- table(rounded.locations$point_id)
  inconsistent.ids <- names(n.locations[n.locations > 1])
  
  if(length(inconsistent.ids) > 0)
    stop(
      "inconsistent outer coordinates for point_id: ",
      paste(inconsistent.ids, collapse = ", ")
    )
  
  # Use the centroid of all reconstructed locations for each point.
  centroids <- aggregate(
    cbind(outer.x, outer.y) ~ point_id,
    data = reconstructed,
    FUN = mean
  )
  
  # Return one row per point, preserving first-appearance order.
  out <- data.frame(
    point_id = unique(point.id),
    stringsAsFactors = FALSE
  )
  
  j <- match(out$point_id, centroids$point_id)
  out$outer.x <- centroids$outer.x[j]
  out$outer.y <- centroids$outer.y[j]
  
  out
}










