# Machine-specific paths. Copy to config.local.R and edit both lines.
# config.local.R is not tracked by git.
#
#   DATA_DIR     downloaded source files. Read only to this pipeline.
#   RESULTS_DIR  everything generated. Must differ from DATA_DIR.
#
# The pipeline stops if CHANGE_ME is still present, so an unedited copy
# fails loudly rather than silently falling back to a default.
#
# Use forward slashes on Windows. No trailing slash.
#   DATA_DIR    <- "D:/glioma/data"
#   RESULTS_DIR <- "D:/glioma/results"

DATA_DIR    <- "CHANGE_ME/path/to/data"
RESULTS_DIR <- "CHANGE_ME/path/to/results"
