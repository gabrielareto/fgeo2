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
  
  get_syntax_calls <- function(fun, known_names) {
    local_dynamic <- character()
    callbacks <- character()
    qualified <- character()
    dynamic_functions <- c("get", "get0", "mget", "match.fun", "do.call")
    callback_arguments <- list(
      apply = 3L, lapply = 2L, sapply = 2L, vapply = 2L,
      tapply = 3L, eapply = 2L, rapply = 2L,
      Map = 1L, mapply = 1L, Reduce = 1L, Filter = 1L,
      Find = 1L, Position = 1L, Negate = 1L, Vectorize = 1L,
      do.call = 1L, match.fun = 1L, outer = 3L, sweep = 4L,
      by = 3L, optim = c(2L, 3L), optimize = 1L,
      nlm = 1L, nlminb = 2:4, uniroot = 1L, integrate = 1L
    )
    callback_argument_names <- c(
      "FUN", ".f", "fn", "gr", "objective", "gradient", "hessian"
    )

    literal_names <- function(x) {
      if (is.character(x)) return(x)

      if (
        is.call(x) &&
        identical(x[[1]], as.name("c"))
      ) {
        values <- as.list(x)[-1]
        if (length(values) && all(vapply(values, is.character, logical(1)))) {
          return(unlist(values, use.names = FALSE))
        }
      }

      character()
    }

    has_explicit_environment <- function(x, op) {
      args <- as.list(x)[-1]
      arg_names <- names(args)
      if (is.null(arg_names)) arg_names <- rep("", length(args))

      if ("envir" %in% arg_names || "pos" %in% arg_names) return(TRUE)

      environment_position <- switch(
        op,
        get = 2L,
        get0 = 2L,
        mget = 2L,
        do.call = 4L,
        Inf
      )

      length(args) >= environment_position &&
        !nzchar(arg_names[[environment_position]])
    }

    walk <- function(x) {
      if (is.pairlist(x)) {
        for (i in seq_along(x)) {
          if (identical(x[[i]], quote(expr = ))) next
          walk(x[[i]])
        }
        return(invisible(NULL))
      }

      if (is.expression(x)) {
        for (i in seq_along(x)) walk(x[[i]])
        return(invisible(NULL))
      }

      if (!is.call(x)) return(invisible(NULL))

      head <- x[[1]]
      op <- if (is.symbol(head)) as.character(head) else NA_character_

      called <- op
      if (
        is.call(head) &&
        length(head) >= 3 &&
        as.character(head[[1]]) %in% c("::", ":::")
      ) {
        called <- as.character(head[[3]])
      }

      if (
        length(x) >= 3 &&
        !is.na(op) &&
        op %in% c("::", ":::")
      ) {
        pkg <- as.character(x[[2]])
        nm <- as.character(x[[3]])

        if (length(pkg) == 1 && length(nm) == 1) {
          qualified <<- c(qualified, paste(pkg, nm, sep = "::"))
        }
      }

      if (
        length(x) >= 2 &&
        !is.na(called) &&
        called %in% dynamic_functions
      ) {
        dynamic_names <- literal_names(x[[2]])
        if (!has_explicit_environment(x, called)) {
          local_dynamic <<- c(local_dynamic, dynamic_names)
        }
      }

      args <- as.list(x)[-1]
      arg_names <- names(args)
      if (is.null(arg_names)) arg_names <- rep("", length(args))

      callback_positions <- which(arg_names %in% callback_argument_names)
      positional_callbacks <- if (!is.na(called)) {
        callback_arguments[[called]]
      } else {
        NULL
      }
      if (!is.null(positional_callbacks)) {
        positional_callbacks <- positional_callbacks[
          positional_callbacks <= length(args) &
            vapply(positional_callbacks, function(i) {
              all(!nzchar(arg_names[seq_len(i)]))
            }, logical(1))
        ]
      }

      positions <- c(callback_positions, positional_callbacks)
      if (!is.na(called) && called %in% known_names) {
        positions <- c(positions, seq_along(args))
      }

      positions <- unique(positions[positions <= length(args)])
      skip_dynamic <- !is.na(called) &&
        called %in% dynamic_functions &&
        has_explicit_environment(x, called)

      for (i in positions) {
        if (skip_dynamic) next

        if (is.symbol(args[[i]])) {
          callbacks <<- c(callbacks, as.character(args[[i]]))
        } else if (i %in% c(callback_positions, positional_callbacks)) {
          local_dynamic <<- c(local_dynamic, literal_names(args[[i]]))
        }
      }

      # Include the call head, where pkg::function is stored.
      for (i in seq_along(x)) walk(x[[i]])

      invisible(NULL)
    }

    if (is.function(fun) && typeof(fun) == "closure") {
      defaults <- formals(fun)
      for (i in seq_along(defaults)) {
        if (identical(defaults[[i]], quote(expr = ))) next
        if (is.symbol(defaults[[i]])) {
          callbacks <- c(callbacks, as.character(defaults[[i]]))
        }
      }

      walk(formals(fun))
      walk(body(fun))
    }

    list(
      local_dynamic = unique(local_dynamic),
      callbacks = unique(callbacks),
      qualified = unique(qualified)
    )
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

    syntax_calls <- get_syntax_calls(fun, unique(meta$function_name))
    
    globals <- tryCatch(
      codetools::findGlobals(fun, merge = FALSE),
      error = function(e) list(functions = character(), variables = character())
    )

    locals <- c(
      names(formals(fun)),
      tryCatch(
        codetools::findFuncLocals(formals(fun), body(fun)),
        error = function(e) character()
      )
    )
    
    calls <- c(
      globals$functions,
      intersect(syntax_calls$callbacks, globals$variables)
    )

    # Character names in get(), do.call(), match.fun(), and related calls.
    calls <- c(
      calls,
      setdiff(syntax_calls$local_dynamic, locals)
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
    qualified <- syntax_calls$qualified
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

# Source all R/new scripts from the latest commit on the default branch.
# load_github_snapshot()
#
# Source a specific script from a full commit SHA into the current environment.
# load_github_snapshot(commit = "a89081ed2e6da798f492d80b4fbb8c80618285b8",
#                      scripts = "03_spatial/R/new/quadrats.R")
#
# Use the latest commit through a date (UTC) and load into a separate environment.
# historical <- new.env()
# load_github_snapshot(commit.date = "2026-11-06", envir = historical)
#
# Save selected files under destination/<commit>/, preserving their paths.
# load_github_snapshot(scripts = "03_spatial/R/new/quadrats.R",
#                      destination = "~/fgeo2-snapshots")
#
# Use repo = "owner/name" for another public repository.
# Install jsonlite to look up a commit or discover script paths.

load_github_snapshot <- function(repo = "gabrielareto/fgeo2", commit = NULL,
                                 commit.date = NULL, scripts = NULL,
                                 destination = NULL, envir = parent.frame())
{
	if(!is.null(commit) && !is.null(commit.date))
		stop("Set either commit or commit.date, not both.")
	if(!is.null(commit) && (length(commit) != 1L || !grepl("^[0-9a-fA-F]{40}$", commit)))
		stop("commit must be a full 40-character SHA.")
	if(is.null(destination) && !is.environment(envir))
		stop("envir must be an R environment.")

	api.url <- paste0("https://api.github.com/repos/", repo)
	if(is.null(commit))
	{
		query <- paste0(api.url, "/commits?per_page=1")
		if(!is.null(commit.date))
		{
			day <- as.Date(commit.date)
			if(length(day) != 1L || is.na(day)) stop("Use one date in YYYY-MM-DD format.")
			query <- paste0(query, "&until=", day, "T23:59:59Z")
		}

		revisions <- jsonlite::fromJSON(query, simplifyVector = FALSE)
		if(!length(revisions)) stop("No matching commit found for ", repo)
		commit <- revisions[[1L]]$sha
	}
	commit <- tolower(commit)

	if(is.null(scripts))
	{
		tree.url <- paste0(api.url, "/git/trees/", commit, "?recursive=1")
		tree <- jsonlite::fromJSON(tree.url)
		if(isTRUE(tree$truncated)) stop("GitHub returned an incomplete file list.")

		paths <- tree$tree$path[tree$tree$type == "blob"]
		scripts <- sort(grep("(^|/)R/new/.*\\.[Rr]$", paths, value = TRUE))
	}

	if(!is.character(scripts) || !length(scripts) || anyNA(scripts) ||
	   any(!nzchar(scripts) | grepl("^/|^[A-Za-z]:|\\\\|(^|/)\\.\\.?(/|$)", scripts)))
		stop("scripts must contain one or more paths relative to the repository root.")

	message(repo, " commit: ", commit)
	raw.url <- paste0("https://raw.githubusercontent.com/", repo, "/", commit, "/")
	directory <- NULL
	if(!is.null(destination))
	{
		directory <- file.path(path.expand(destination), commit)
		temporary <- character()
		on.exit(unlink(temporary), add = TRUE)
	}

	for(file in scripts)
	{
		parts <- strsplit(file, "/", fixed = TRUE)[[1L]]
		url.path <- paste(vapply(parts, utils::URLencode, "",
		                         reserved = TRUE, repeated = TRUE), collapse = "/")
		url <- paste0(raw.url, url.path)

		if(is.null(directory))
		{
			source(url, local = envir, encoding = "UTF-8")
			next
		}

		path <- file.path(directory, file)
		dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
		temporary <- tempfile(tmpdir = dirname(path))
		status <- utils::download.file(url, temporary, mode = "wb")
		if(status != 0L || !file.rename(temporary, path))
			stop("Could not save ", file, " to ", path)
	}

	invisible(list(commit = commit, scripts = scripts, directory = directory))
}
