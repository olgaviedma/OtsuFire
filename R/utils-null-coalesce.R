# Scalar-safe null-coalescing operator. Internal only — NOT exported.
# Returns `y` when `x` is NULL or scalar NA; otherwise `x`. Semantics frozen
# by handoff Section 58b (F9 cleanup).

`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1L && is.na(x[[1L]]))) y else x
}
