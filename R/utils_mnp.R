# ============================================================
# R/utils_mnp.R
# Shared helpers for PTB-MNP analysis
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

safe_num <- function(x) {
  if (is.numeric(x)) return(x)
  if (is.logical(x)) return(as.numeric(x))
  suppressWarnings(readr::parse_number(as.character(x)))
}

safe_chr <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  dplyr::na_if(x, "")
}

safe_log10p1 <- function(x) {
  x <- safe_num(x)
  log10(x + 1)
}

safe_log2p1 <- function(x) {
  x <- safe_num(x)
  log2(x + 1)
}

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

pull_existing <- function(df, candidates, default = NA) {
  nm <- first_existing_col(df, candidates)
  if (is.na(nm)) return(rep(default, nrow(df)))
  df[[nm]]
}

as_01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x) || is.integer(x)) {
    x <- suppressWarnings(as.numeric(x))
    return(as.integer(ifelse(is.na(x), NA, ifelse(x != 0, 1, 0))))
  }
  v <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    v %in% c("1","yes","y","true","t","present","positive") ~ 1L,
    v %in% c("0","no","n","false","f","none","absent","negative") ~ 0L,
    TRUE ~ NA_integer_
  )
}

as_yesno_factor <- function(x) {
  z <- as_01(x)
  factor(ifelse(is.na(z), NA_character_, ifelse(z == 1, "Yes", "No")), levels = c("No", "Yes"))
}

pretty_p <- function(p, digits = 3) {
  out <- rep(NA_character_, length(p))
  ok <- !is.na(p)
  out[ok & p < 0.001] <- "<0.001"
  out[ok & p >= 0.001] <- formatC(p[ok & p >= 0.001], format = "f", digits = digits)
  out
}

fmt_n_pct <- function(n, N) {
  if (is.na(n) || is.na(N) || N == 0) return(NA_character_)
  sprintf("%d/%d (%.1f%%)", n, N, 100 * n / N)
}

median_iqr <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  q <- stats::quantile(x, c(0.25, 0.50, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), q[2], q[1], q[3])
}

mean_sd <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  sprintf(paste0("%.", digits, "f \u00B1 %.", digits, "f"), mean(x), stats::sd(x))
}

poly_lower <- function(x) tolower(as.character(x))

poly_col_ug_g <- function(poly) paste0("mnp_", tolower(poly), "_ug_g")
poly_col_sum_ug <- function(poly) paste0("mnp_", tolower(poly), "_sum_ug")
legacy_col_ug_g <- function(poly) paste0(toupper(poly), "_ug_per_g")
legacy_col_log10 <- function(poly) paste0(toupper(poly), "_log10ugpg")

poly_family <- function(exposure) {
  p <- tolower(exposure)
  p <- sub("_log10ugpg$", "", p)
  p <- sub("_ug_per_g$", "", p)
  p <- sub("^mnp_", "", p)
  p <- sub("_ug_g$", "", p)
  p <- trimws(p)

  dplyr::case_when(
    p %in% c("pp","pe") ~ "polyolefins",
    p %in% c("pet") ~ "polyesters",
    p %in% c("pu") ~ "polyurethanes",
    p %in% c("pvc") ~ "vinyls",
    p %in% c("pmma") ~ "methacrylates",
    p %in% c("pc") ~ "polycarbonates",
    p %in% c("n6","n66") ~ "polyamides",
    p %in% c("ps","abs") ~ "styrenics",
    p %in% c("sbr") ~ "elastomers",
    TRUE ~ "other"
  )
}

make_poly_palette <- function(poly_codes_upper, mode = c("familyish","distinct")) {
  mode <- match.arg(mode)
  if (!requireNamespace("RColorBrewer", quietly = TRUE) ||
      !requireNamespace("colorspace", quietly = TRUE)) {
    stop("Need RColorBrewer + colorspace. install.packages(c('RColorBrewer','colorspace'))", call. = FALSE)
  }

  paired12 <- RColorBrewer::brewer.pal(12, "Paired")
  family_pal <- c(
    polyolefins     = paired12[2],
    polyesters      = paired12[4],
    polyurethanes   = paired12[6],
    vinyls          = paired12[8],
    methacrylates   = paired12[10],
    polycarbonates  = paired12[12],
    polyamides      = paired12[1],
    styrenics       = paired12[3],
    elastomers      = paired12[5],
    other           = "grey70"
  )

  if (mode == "distinct") {
    cols <- paired12[seq_along(poly_codes_upper)]
    names(cols) <- poly_codes_upper
    return(cols)
  }

  make_shades <- function(base_col, k) {
    if (k <= 1) return(base_col)
    t <- seq(-0.22, 0.22, length.out = k)
    vapply(t, function(tt) {
      if (tt < 0) colorspace::darken(base_col, abs(tt)) else colorspace::lighten(base_col, tt)
    }, character(1))
  }

  lower <- tolower(poly_codes_upper)
  fam <- vapply(lower, poly_family, character(1))
  out <- setNames(rep(NA_character_, length(poly_codes_upper)), poly_codes_upper)

  for (f in unique(fam)) {
    idx <- which(fam == f)
    out[idx] <- make_shades(family_pal[[f]], length(idx))
  }
  out
}

write_csv_clean <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, path, na = "")
  invisible(path)
}


pad_zip <- function(x) {
  x_chr <- as.character(x)
  x_chr <- stringr::str_replace_all(x_chr, "[^0-9]", "")
  x_chr[x_chr == ""] <- NA_character_
  too_long <- !is.na(x_chr) & nchar(x_chr) > 5
  x_chr[too_long] <- stringr::str_sub(x_chr[too_long], 1, 5)
  stringr::str_pad(x_chr, width = 5, side = "left", pad = "0")
}

yes_no_like <- function(x) {
  v <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    v %in% c("yes","y","true","t","1","present","positive") ~ TRUE,
    v %in% c("no","n","false","f","0","none","absent","negative") ~ FALSE,
    TRUE ~ NA
  )
}

derive_ethnicity_race <- function(dat) {

  eth <- safe_chr(
    pull_existing(dat, c("Ethnicity", "ethnicity"),
                  default = NA_character_)
  )

  eth_clean <- stringr::str_squish(
    stringr::str_to_lower(eth)
  )

  # Explicitly distinguish Hispanic, non-Hispanic,
  # and unknown ethnicity.
  non_hisp <- !is.na(eth_clean) &
    stringr::str_detect(
      eth_clean,
      "non[- ]?hispanic|not hispanic|non[- ]?latino"
    )

  hisp <- !is.na(eth_clean) &
    !non_hisp &
    stringr::str_detect(
      eth_clean,
      "hispanic|latino|latina|latinx"
    )

  ethnicity_group <- dplyr::case_when(
    non_hisp ~ "Non-Hispanic",
    hisp ~ "Hispanic",
    TRUE ~ "Unknown ethnicity"
  )

  race_cols <- c(
    ai_an = "Race - American Indian/Alaska Native",
    asian = "Race - Asian",
    nhpi = "Race - Native Hawaiian or other Pacific Islander",
    black = "Race - Black or African American",
    white = "Race - White",
    unknown = "Race - Not reported/Unknown"
  )

  get_yes <- function(nm) {
    if (!nm %in% names(dat)) {
      return(rep(FALSE, nrow(dat)))
    }
    x <- yes_no_like(dat[[nm]])
    ifelse(is.na(x), FALSE, x)
  }

  race_yes <- tibble::tibble(
    ai_an = get_yes(race_cols[["ai_an"]]),
    asian = get_yes(race_cols[["asian"]]),
    nhpi = get_yes(race_cols[["nhpi"]]),
    black = get_yes(race_cols[["black"]]),
    white = get_yes(race_cols[["white"]]),
    unknown = get_yes(race_cols[["unknown"]])
  )

  n_races <- rowSums(
    race_yes[, c("ai_an", "asian", "nhpi",
                 "black", "white")],
    na.rm = TRUE
  )

  race_label <- dplyr::case_when(
    race_yes$unknown ~ "Unknown/Not reported",
    n_races > 1 ~ "Multiracial",
    race_yes$ai_an ~ "American Indian/Alaska Native",
    race_yes$asian ~ "Asian",
    race_yes$nhpi ~ "Native Hawaiian/Other Pacific Islander",
    race_yes$black ~ "Black or African American",
    race_yes$white ~ "White",
    TRUE ~ NA_character_
  )

  dplyr::case_when(
    ethnicity_group == "Hispanic" &
      is.na(race_label) ~ "Hispanic/Latino",

    ethnicity_group == "Hispanic" ~
      paste("Hispanic", race_label),

    ethnicity_group == "Non-Hispanic" &
      !is.na(race_label) ~
      paste("Non-Hispanic", race_label),

    ethnicity_group == "Unknown ethnicity" &
      !is.na(race_label) ~
      paste("Unknown ethnicity", race_label),

    TRUE ~ "Unknown ethnicity/race"
  )
}


derive_prior_ptb_from_parity <- function(parity) {
  # Parity often encoded as: term preterm abortions living
  # Example: "2 0 0 2" or "1 0 0 1"
  p <- stringr::str_squish(as.character(parity))
  nums <- stringr::str_extract_all(p, "\\d+")
  preterm_n <- purrr::map_dbl(nums, function(z) {
    if (length(z) >= 2) as.numeric(z[[2]]) else NA_real_
  })
  as.integer(ifelse(is.na(preterm_n), NA_integer_, ifelse(preterm_n > 0, 1L, 0L)))
}

collapse_hospital <- function(x) {
  v <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    stringr::str_detect(v, "ben\\s*taub|btgh|harris\\s*health") ~ "Ben Taub",
    stringr::str_detect(v, "tch|texas\\s*children|pavilion") ~ "TCH",
    is.na(v) | v == "" ~ NA_character_,
    TRUE ~ as.character(x)
  )
}
