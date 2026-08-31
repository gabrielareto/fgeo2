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


# Newer, more powerful version that accepts specific packages too.
get_all_dependencies <- function(
    packages = character(),
    include_global = TRUE,
    internal = TRUE,
    resolve_imports = TRUE,
    resolve_unique_names = FALSE,
    envir = .GlobalEnv
) {
  make_id <- function(pkg, nm) paste(pkg, nm, sep = "::")
  key <- function(pkg, nm) paste(pkg, nm, sep = "\r")
  
  drop_calls <- c(
    "{", "(", "if", "for", "while", "repeat", "function", "return",
    "<-", "<<-", "=", "[", "[[", "$", "@", "::", ":::",
    "+", "-", "*", "/", "^", "%%", "%/%", ":", "~",
    "!", "!=", "==", "<", ">", "<=", ">=",
    "&", "&&", "|", "||"
  )
  
  get_env_functions <- function(env, pkg, nms) {
    vals <- lapply(nms, function(nm) {
      tryCatch(get(nm, envir = env, inherits = FALSE), error = function(e) NULL)
    })
    
    keep <- vapply(vals, is.function, logical(1))
    vals <- vals[keep]
    nms <- nms[keep]
    
    setNames(vals, make_id(pkg, nms))
  }
  
  get_global_functions <- function() {
    nms <- ls(envir, all.names = TRUE)
    get_env_functions(envir, ".GlobalEnv", nms)
  }
  
  get_package_functions <- function(pkg) {
    ns <- asNamespace(pkg)
    
    nms <- if (internal) {
      ls(ns, all.names = TRUE)
    } else {
      getNamespaceExports(pkg)
    }
    
    get_env_functions(ns, pkg, unique(nms))
  }
  
  get_qualified_calls <- function(fun) {
    out <- character()
    
    walk <- function(x) {
      if (!is.call(x)) return(invisible(NULL))
      
      head <- x[[1]]
      
      op <- if (is.symbol(head)) {
        as.character(head)
      } else {
        NA_character_
      }
      
      if (
        length(x) >= 3 &&
        !is.na(op) &&
        op %in% c("::", ":::")
      ) {
        pkg <- as.character(x[[2]])
        nm  <- as.character(x[[3]])
        
        if (length(pkg) == 1 && length(nm) == 1) {
          out <<- c(out, paste(pkg, nm, sep = "::"))
        }
      }
      
      for (i in seq_along(x)[-1]) {
        walk(x[[i]])
      }
      
      invisible(NULL)
    }
    
    if (is.function(fun) && typeof(fun) == "closure") {
      walk(body(fun))
    }
    
    unique(out)
  }
  
  get_import_map <- function(packages, known_ids) {
    out <- list()
    
    for (pkg in packages) {
      imports <- tryCatch(getNamespaceImports(pkg), error = function(e) list())
      
      for (ipkg in names(imports)) {
        x <- imports[[ipkg]]
        
        if (is.environment(x)) {
          local <- remote <- ls(x, all.names = TRUE)
        } else if (is.character(x)) {
          remote <- unname(x)
          local <- names(x)
          if (is.null(local) || any(local == "")) local <- remote
        } else {
          next
        }
        
        for (i in seq_along(local)) {
          id <- make_id(ipkg, remote[i])
          if (id %in% known_ids) {
            out[[key(pkg, local[i])]] <- id
          }
        }
      }
    }
    
    out
  }
  
  # Collect functions ---------------------------------------------------------
  
  funs <- list()
  
  if (include_global) {
    funs <- c(funs, get_global_functions())
  }
  
  for (pkg in packages) {
    funs <- c(funs, get_package_functions(pkg))
  }
  
  if (!length(funs)) {
    return(data.frame(
      from_id = character(),
      to_id = character(),
      from_package = character(),
      from_function = character(),
      to_package = character(),
      to_function = character()
    ))
  }
  
  ids <- names(funs)
  
  meta <- data.frame(
    id = ids,
    package = sub("::.*$", "", ids),
    function_name = sub("^.*::", "", ids),
    stringsAsFactors = FALSE
  )
  
  local_lookup <- split(meta$id, key(meta$package, meta$function_name))
  name_lookup  <- split(meta$id, meta$function_name)
  
  import_map <- if (resolve_imports) {
    get_import_map(packages, ids)
  } else {
    list()
  }
  
  # Build edges ---------------------------------------------------------------
  
  edge_list <- Map(function(to_id, fun) {
    to <- meta[match(to_id, meta$id), ]
    
    if (!is.function(fun) || typeof(fun) != "closure") return(NULL)
    
    calls <- tryCatch(
      codetools::findGlobals(fun, merge = FALSE)$functions,
      error = function(e) character()
    )
    
    calls <- unique(setdiff(calls, drop_calls))
    
    from_ids <- character()
    
    for (called in calls) {
      # 1. Same namespace first
      k_local <- key(to$package, called)
      hit <- local_lookup[[k_local]]
      
      if (!is.null(hit)) {
        from_ids <- c(from_ids, hit)
        next
      }
      
      # 2. Imported function
      k_import <- key(to$package, called)
      hit <- import_map[[k_import]]
      
      if (!is.null(hit)) {
        from_ids <- c(from_ids, hit)
        next
      }
      
      # 3. Optional weak fallback: unique name among all included functions
      if (resolve_unique_names) {
        hit <- name_lookup[[called]]
        if (!is.null(hit) && length(hit) == 1) {
          from_ids <- c(from_ids, hit)
        }
      }
    }
    
    # Explicit pkg::fun or pkg:::fun calls
    qualified <- get_qualified_calls(fun)
    qualified <- qualified[qualified %in% ids]
    
    from_ids <- unique(c(from_ids, qualified))
    from_ids <- setdiff(from_ids, to_id)
    
    if (!length(from_ids)) return(NULL)
    
    from <- meta[match(from_ids, meta$id), ]
    
    data.frame(
      from_id = from$id,
      to_id = to$id,
      from_package = from$package,
      from_function = from$function_name,
      to_package = to$package,
      to_function = to$function_name,
      stringsAsFactors = FALSE
    )
  }, ids, funs)
  
  edges <- do.call(rbind, edge_list)
  
  if (is.null(edges)) {
    edges <- data.frame(
      from_id = character(),
      to_id = character(),
      from_package = character(),
      from_function = character(),
      to_package = character(),
      to_function = character()
    )
  }
  
  rownames(edges) <- NULL
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

  if (!length(commit_id) || !nzchar(commit_id[1])) {
    stop(
      "No commit was found on or before ", as.Date(date), ". ",
      "If this is a shallow clone, run 'git fetch --unshallow' first.",
      call. = FALSE
    )
  }
  
  commit_id[1]
}


make_dir_to_point_in_history <- function(path, commit)
{
  path <- normalizePath(path)

  # Remove registrations left by interrupted evaluations.
  status <- system2(
    "git",
    c("-C", shQuote(path), "worktree", "prune")
  )

  if (!identical(status, 0L)) {
    stop("git worktree prune failed.", call. = FALSE)
  }
  
  # Resolve commit / tag to full commit hash
  commit_id <- system2(
    "git",
    c("-C", shQuote(path), "rev-parse", paste0(commit, "^{commit}")),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!is.null(attr(commit_id, "status"))) {
    stop(
      paste(commit_id, collapse = "\n"),
      "\nIf this is a shallow clone, run 'git fetch --unshallow' first.",
      call. = FALSE
    )
  }
  
  commit_id <- commit_id[1]
  
  # Temporary historical project path
  historical_path <- file.path(
    tempdir(),
    paste0(basename(path), "_at_", substr(commit_id, 1, 12))
  )
  
  # If it already exists, reuse it only if it points to the right commit
  if (dir.exists(historical_path)) {
    
    existing_commit <- system2(
      "git",
      c("-C", shQuote(historical_path), "rev-parse", "HEAD"),
      stdout = TRUE,
      stderr = TRUE
    )
    
    if (!is.null(attr(existing_commit, "status"))) {
      stop(
        "The historical path already exists but is not a valid Git worktree:\n",
        historical_path,
        call. = FALSE
      )
    }
    
    if (!identical(existing_commit[1], commit_id)) {
      stop(
        "The historical path already exists but points to a different commit:\n",
        historical_path,
        call. = FALSE
      )
    }
    
    warning(
      "Historical path already exists; reusing it:\n",
      historical_path,
      call. = FALSE
    )
    
    return(historical_path)
  }
  
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

