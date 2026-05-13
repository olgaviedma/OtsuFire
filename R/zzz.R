# Package load hooks. Kept minimal on purpose; heavy defaults belong in
# the relevant builder functions, not at load time.

.onLoad <- function(libname, pkgname) {
  # Option defaults can be set here once implementation starts.
  invisible()
}

.onAttach <- function(libname, pkgname) {
  # Intentionally silent.
  invisible()
}
