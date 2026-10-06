# Install CRAN packages used by the manuscript analysis scripts.
packages <- c(
  "MASS", "RColorBrewer", "Rtsne", "brglm2", "cluster", "colorspace",
  "dplyr", "forcats", "ggplot2", "ggrepel", "glue", "gt", "lmtest",
  "logistf", "openxlsx", "ordinal", "patchwork", "purrr", "readr",
  "readxl", "sandwich", "scales", "stringr", "tibble", "tidyr",
  "tidyverse", "uwot", "vegan"
)

installed <- rownames(installed.packages())
missing <- setdiff(packages, installed)

if (length(missing)) {
  install.packages(missing)
} else {
  message("All required CRAN packages are already installed.")
}
