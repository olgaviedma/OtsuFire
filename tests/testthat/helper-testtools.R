#' Get tool path from environment variable
#'
#' Retrieves the path to an external tool (e.g., GDAL or Python) using an environment variable.
#' If the variable is not defined or the file does not exist, it safely returns NULL and warns the user.
#'
#' @param var_name The name of the environment variable to look up.
#' @param must_exist Logical. If TRUE (default), will warn and return NULL if the file doesn't exist.
#' @return A character path if found, or NULL otherwise.
get_env_tool <- function(var_name, must_exist = TRUE) {
  path <- Sys.getenv(var_name, unset = NA_character_)
  
  if (is.na(path) || path == "") {
    warning(glue::glue("Environment variable '{var_name}' is not set."))
    return(NULL)
  }
  
  if (must_exist && !file.exists(path)) {
    warning(glue::glue("Path defined by '{var_name}' does not exist: {path}"))
    return(NULL)
  }
  
  return(path)
}
