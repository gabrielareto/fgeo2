# Functions to handle quadrats, splits, local coordinates, and global coordinates
# Gabriel Arellano
# gabriel.arellano.torres@gmail.com




# Auxiliary, to make different types of sequences of letters:

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

# 2)  Labels with exactly k characters:
# k = 2 gives aa, ab, ..., az, ba, bb, ..., bz, ...
letter_seq_fixed <- function(n, k) {
  grid <- expand.grid(rep(list(letters), k), stringsAsFactors = FALSE)
  do.call(paste0, rev(grid))[1:n]
}




# Quadrat-IDing conventions

get_quadrat_locations_from_id_convention <- function(
  # PARAMETERS
  # Dimensions of the large unit
  xmin = NA,
  ymin = NA,
  xmax = NA,
  ymax = NA,
  
  # Dimensions of the smaller sub-units
  quadrat.side.x = NA,
  quadrat.side.y = NA,
  
  # How the ID is built
  built = c("col+row", "row+col", "sequence"),
  example.part.for.col = "01", # c("A", "1", "01", "000001")
  example.part.for.row = "01", # c("A", "1", "01", "000001")
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
  
  if(is.na(xmax) | is.null(ymax) | is.na(ymax) | is.null(ymax))
    stop("specify dimensions through xmax and ymax")
  
  if(is.na(xmin) | is.null(xmin))
  {
    xmin = 0
    warning("xmin = 0 assumed by default")
  }
    
  if(is.na(ymin) | is.null(ymin))
  {
    ymin = 0
    warning("ymin = 0 assumed by default")
  }
  
  if(length(built) != 1)
  {
    built = "col+row"
    warning("built = 'col+row' assumed by default")
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
  xmins <- seq(from = xmin, to = xmax, by = quadrat.side.x)
  xmins <- xmins[xmins < xmax]
  ymins <- seq(from = ymin, to = ymax, by = quadrat.side.y)
  ymins <- ymins[ymins < ymax]
  
  # Add corner coordinates to col+row indexes
  grid.exact  <- expand.grid(xmin = xmins, ymin = ymins)
  grid.exact$xmax <- pmin(xmax, grid.exact$xmin + quadrat.side.x)
  grid.exact$ymax <- pmin(ymax, grid.exact$ymin + quadrat.side.y)
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
    use.numbers.for.col1 = is.numeric(as.numeric(example.part.for.col))
    use.numbers.for.col2 = is.numeric(as.numeric(start.col.at))
    if(use.numbers.for.col1 != use.numbers.for.col2)
      stop("numbers vs. letters conflict between 'example' and 'start at' for 'col' part")
    
    use.numbers.for.row1 = is.numeric(as.numeric(example.part.for.row))
    use.numbers.for.row2 = is.numeric(as.numeric(start.row.at))
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
    use.numbers.for.sequence1 = is.numeric(as.numeric(example.part.for.sequence))
    use.numbers.for.sequence2 = is.numeric(as.numeric(start.sequence.at))
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
    warning("the quadrats don't fit perfectly along x")
  
  if(any(grid$ymax - grid$ymin != quadrat.side.y))
    warning("the quadrats don't fit perfectly along y")
  
  
  # Trim and rename the output:
  out <- grid[,c("entire", "xmin", "ymin", "xmax", "ymax"), drop = FALSE]
  colnames(out)[1] <- "quadrat_id"
  rownames(out) <- out$quadrat_id

  # Sort in the same way, for consistency:
  o <- order(out$xmin, out$ymin)
  out <- out[o,,drop = FALSE]
  out
}
  


qq <- get_quadrat_locations_from_id_convention(xmax = 100,
                                         ymax = 100,
                                         quadrat.side.x = 20,
                                         example.part.for.col = "00",
                                         example.part.for.row = "00")

head(qq)[,c(1,2,4,3,5), 10]


# para ir para atrás:
# partir por separador: se puede detectar
# quitar no alfanuméricos
# chequear letras vs. numeros
# si letras: posiciones en secuencias que se construyen con esas auxiliares
# si numeros: quitar ceros, muy directo
# si secuencia: resolver los órdenes

