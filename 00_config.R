# ============================================================
# 00_config.R
# Public project configuration for PTB-MNP analysis
# ============================================================
# Run scripts from the repository root. By default PROJECT_ROOT is the
# current working directory. Override with MNP_PROJECT_ROOT if needed.
#
# Participant-level and validation-source data are not distributed in this
# repository. See data/README.md for expected inputs and access notes.
# ============================================================

PROJECT_ROOT <- normalizePath(
  Sys.getenv("MNP_PROJECT_ROOT", unset = getwd()),
  mustWork = TRUE
)
setwd(PROJECT_ROOT)

# Restricted/model-ready inputs. These paths may be overridden with
# environment variables without editing tracked code.
DIR_MODEL_READY <- Sys.getenv(
  "MNP_MODEL_READY_DIR",
  unset = file.path("data", "restricted", "model_ready")
)

PERIBANK_XLSX <- Sys.getenv(
  "MNP_PERIBANK_XLSX",
  unset = file.path("data", "restricted", "peribank_metadata.xlsx")
)

SVI_PATH <- Sys.getenv(
  "MNP_SVI_PATH",
  unset = file.path("data", "external", "SVI_2022_US_ZCTA.csv")
)

SDI_PATH <- Sys.getenv(
  "MNP_SDI_PATH",
  unset = file.path("data", "external", "rgcsdi-2015-2019-zcta.csv")
)

# Outputs
DIR_TABLES <- file.path("outputs", "tables")
DIR_FIGS   <- file.path("outputs", "figs")
DIR_LOGS   <- file.path("outputs", "logs")
DIR_DIAG   <- file.path("outputs", "diagnostics")

for (d in c(DIR_TABLES, DIR_FIGS, DIR_LOGS, DIR_DIAG)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Canonical polymer order used across analyses.
POLYMERS <- c("PE", "SBR", "PP", "PVC", "N66", "PET", "N6", "PMMA", "PC", "PU", "ABS", "PS")

POLYMER_FULL_NAMES <- c(
  PE   = "Polyethylene",
  PP   = "Polypropylene",
  PVC  = "Polyvinyl chloride",
  PET  = "Polyethylene terephthalate",
  N66  = "Nylon-66",
  N6   = "Nylon-6",
  PU   = "Polyurethane",
  PMMA = "Polymethyl methacrylate",
  PC   = "Polycarbonate",
  ABS  = "Acrylonitrile-butadiene-styrene",
  PS   = "Polystyrene",
  SBR  = "Styrene-butadiene rubber"
)

DETECTION_EPS <- 0
LOG10_OFFSET  <- 1
MIN_DETECTED_N_FOR_CONTINUOUS_MODEL <- 10L
FORCE_PE_CONTROL_BLANK_ZERO <- TRUE
