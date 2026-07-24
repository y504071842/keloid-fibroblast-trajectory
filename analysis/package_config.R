script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

repo_root_env <- Sys.getenv("KFT_REPO_ROOT", unset = "")
repo_root <- if (nzchar(repo_root_env)) {
  normalizePath(repo_root_env, winslash = "/", mustWork = FALSE)
} else {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
}

package_root <- repo_root
audit_root <- Sys.getenv("KFT_RESULTS_ROOT", unset = file.path(repo_root, "results"))
input_root <- Sys.getenv("KFT_INPUT_ROOT", unset = file.path(repo_root, "data", "frozen_inputs"))
reference_root <- file.path(repo_root, "reference")

path_in_package <- function(...) file.path(repo_root, ...)

required_input <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop(
      label, " was not found: ", path,
      "\nSee data/README.md for acquisition and environment-variable instructions.",
      call. = FALSE
    )
  }
  path
}

external_object <- function(env_var, default_path, label) {
  required_input(Sys.getenv(env_var, unset = default_path), label)
}
