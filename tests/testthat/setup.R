# Test setup file
# Clean up any dppackage output files that might be created during tests

# Function to clean up dppackage output files
cleanup_dppackage_files <- function() {
  # Remove any dppackage output files
  dppackage_files <- list.files(
    pattern = "dppackage.*\\.out$", 
    full.names = TRUE, 
    recursive = TRUE
  )
  if (length(dppackage_files) > 0) {
    unlink(dppackage_files)
  }
}

# Clean up on test start
cleanup_dppackage_files()

# Also clean up on test end
.onLoad <- function(libname, pkgname) {
  cleanup_dppackage_files()
}
