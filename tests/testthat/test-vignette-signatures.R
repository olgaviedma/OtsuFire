# Guard: vignette code must call the package with arguments that exist.
#
# WHY THIS TEST EXISTS
# The vignettes set `eval = FALSE` (their examples need full-scale annual
# mosaics) and `purl = FALSE` (so R CMD check does not extract the code
# either). The result is that vignette code is rendered but never executed
# and never parsed against the package, so a stale example survives every
# check indefinitely. Three defects reached the 2.0.0 release branch that
# way: `build_supervised_burned_config(scenario = ...)` after the rename to
# `run_label`, plus `change_index_mosaic()` and `validate_fire_maps()`
# documented with argument names that never existed.
#
# WHAT IT CHECKS
# Every call to an exported OtsuFire function inside a vignette code chunk is
# matched against that function's formals:
#   * a named argument that is not a formal fails the test;
#   * more supplied arguments than the function accepts fails the test.
# Both are skipped for functions taking `...`.
#
# WHAT IT CANNOT CHECK
# Argument ORDER in positional calls. `score_burned_patches(cfg, patches)`
# instead of `score_burned_patches(patches, cfg)` has the right arity and the
# right names, and only fails at run time. Prefer named arguments in vignette
# examples so that this guard can see them.
#
# It deliberately does NOT execute anything: it only parses.

vignette_dir <- testthat::test_path("..", "..", "vignettes")

# Extract the R code from ```{r ...} fenced chunks of an .Rmd file.
.of_vignette_chunk_code <- function(path) {
  ln <- readLines(path, warn = FALSE)
  starts <- grep("^```\\{r", ln)
  ends   <- grep("^```\\s*$", ln)
  code <- character(0)
  for (s in starts) {
    e <- ends[ends > s]
    if (!length(e)) next
    if (e[1] > s + 1) code <- c(code, ln[(s + 1):(e[1] - 1)])
  }
  code
}

# Collect "function(arg = ...)" problems from parsed code.
.of_bad_arguments <- function(parsed, exports) {
  problems <- character(0)
  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))
    fn <- x[[1]]
    if (is.name(fn)) {
      nm <- as.character(fn)
      if (nm %in% exports) {
        fmls <- names(formals(get(nm, envir = asNamespace("OtsuFire"))))
        given <- names(as.list(x))[-1]
        given <- given[!is.na(given) & nzchar(given)]
        if (!("..." %in% fmls)) {
          unknown <- setdiff(given, fmls)
          if (length(unknown)) {
            problems <<- c(problems, sprintf(
              "%s(): argument(s) not in the signature: %s",
              nm, paste(unknown, collapse = ", ")
            ))
          }
          n_supplied <- length(as.list(x)) - 1L
          if (n_supplied > length(fmls)) {
            problems <<- c(problems, sprintf(
              "%s(): %d arguments supplied but the function accepts %d",
              nm, n_supplied, length(fmls)
            ))
          }
        }
      }
    }
    for (i in seq_along(x)) {
      # An empty symbol is R's marker for a missing argument, as in `x[i, ]`.
      # It cannot be inspected at all: touching it raises "argument is
      # missing, with no default". Skip exactly that error and nothing else.
      tryCatch(
        walk(x[[i]]),
        error = function(e) {
          if (!grepl("is missing, with no default", conditionMessage(e), fixed = TRUE)) {
            stop(e)
          }
        }
      )
    }
  }
  for (ex in parsed) walk(ex)
  unique(problems)
}

test_that("vignette code calls the package with arguments that exist", {
  skip_if_not(dir.exists(vignette_dir), "vignettes/ not available")
  rmds <- list.files(vignette_dir, pattern = "[.]Rmd$", full.names = TRUE)
  skip_if(length(rmds) == 0L, "no vignettes found")

  exports <- getNamespaceExports("OtsuFire")

  for (rmd in rmds) {
    code <- .of_vignette_chunk_code(rmd)
    if (!length(code)) next

    parsed <- tryCatch(
      parse(text = paste(code, collapse = "\n")),
      error = function(e) {
        fail(sprintf("%s: vignette code does not parse: %s",
                     basename(rmd), conditionMessage(e)))
        NULL
      }
    )
    if (is.null(parsed)) next

    problems <- .of_bad_arguments(parsed, exports)
    expect_identical(
      problems, character(0),
      info = sprintf("%s\n  %s", basename(rmd), paste(problems, collapse = "\n  "))
    )
  }
})
