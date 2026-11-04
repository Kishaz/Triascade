# One-command launcher for CascadeScope.
# From this folder:  Rscript run.R   (or in R:  source("run.R"))
# Installs any missing runtime packages from CRAN, then starts the app.
pkgs <- c("shiny", "bslib", "plotly", "DT", "ggplot2", "patchwork", "readxl", "openxlsx")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  message("Installing missing packages: ", paste(miss, collapse = ", "))
  install.packages(miss, repos = "https://cloud.r-project.org")
}
shiny::runApp(".", launch.browser = TRUE)
