# Useful code to evaluate code
# Gabriel Arellano
# gabriel.arellano.torres@gmail.com


################################################################
### FUNCTIONS TO LIST FUNCTIONS AND MAP DEPENDENCIES
################################################################

# Function to return all dependencies among
# all functions in the global environment.
get_all_dependencies <- function()
{
  f <- mget(ls(.GlobalEnv), .GlobalEnv)
  f <- f[vapply(f, is.function, TRUE)]
  n <- names(f)
  
  edges <- do.call(rbind, Map(\(nm, fun) {
    z <- intersect(codetools::findGlobals(fun, merge = FALSE)$functions, n)
    z <- setdiff(z, nm)
    if (length(z)) data.frame(from = z, to = nm)
  }, n, f))
  edges
}

# All downstream and upstream dependencies of
# a given function, based on a set of dependencies.

downstream_dependents <- function(x, edges) {
  out <- character()
  repeat {
    new <- setdiff(edges$to[edges$from %in% c(x, out)], out)
    if (!length(new)) break
    out <- c(out, new)
  }
  out
}

upstream_dependencies <- function(x, edges) {
  out <- character()
  repeat {
    new <- setdiff(edges$from[edges$to %in% c(x, out)], out)
    if (!length(new)) break
    out <- c(out, new)
  }
  out
}


# List the functions in a target script:
functions_in_script <- function(file) {
  x <- parse(file)
  
  is_fun_assign <- function(e)
    is.call(e) &&
    as.character(e[[1]]) %in% c("<-", "=", "<<-") &&
    is.call(e[[3]]) &&
    identical(e[[3]][[1]], as.name("function"))
  
  vapply(x[vapply(x, is_fun_assign, logical(1))], \(e) as.character(e[[2]]), "")
}


################################################################
### FUNCTIONS TO BE ABLE TO READ FROM THE HISTORY OF THE PROJECT
################################################################
# This will be necessary for reproducibility of the
# evaluation scripts, as the code itself will be
# updated.

# Maybe equivalent to gert::git_worktree_add() and gert::git_worktree_remove()

commit_from_date <- function(path, date, ref = "HEAD")
{
  path <- normalizePath(path)
  
  cutoff <- paste0(as.Date(date), " 23:59:59")
  
  commit_id <- system2(
    "git",
    c(
      "-C", shQuote(path),
      "rev-list", "-n", "1",
      paste0("--before=", shQuote(cutoff)),
      ref
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!is.null(attr(commit_id, "status"))) {
    stop(paste(commit_id, collapse = "\n"), call. = FALSE)
  }
  
  commit_id[1]
}


make_dir_to_point_in_history <- function(path, commit)
{
  path <- normalizePath(path)
  
  # Resolve commit / tag to full commit hash
  commit_id <- system2(
    "git",
    c(
      "-C", shQuote(path),
      "rev-parse", paste0(commit, "^{commit}")
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!is.null(attr(commit_id, "status"))) {
    stop(paste(commit_id, collapse = "\n"), call. = FALSE)
  }
  
  commit_id <- commit_id[1]
  
  # Temporary historical project path
  historical_path <- file.path(
    tempdir(),
    paste0(basename(path), "_at_", substr(commit_id, 1, 12))
  )
  
  # Create historical checkout
  status <- system2(
    "git",
    c(
      "-C", shQuote(path),
      "worktree", "add", "--detach",
      shQuote(historical_path),
      commit_id
    )
  )
  
  if (!identical(status, 0L)) {
    stop("git worktree add failed.", call. = FALSE)
  }
  
  historical_path
}



remove_dir_to_point_in_history <- function(path, historical_path)
{
  path <- normalizePath(path)
  
  status <- system2(
    "git",
    c(
      "-C", shQuote(path),
      "worktree", "remove", "--force",
      shQuote(historical_path)
    )
  )
  
  if (!identical(status, 0L)) {
    stop("git worktree remove failed.", call. = FALSE)
  }
  
  invisible(TRUE)
}