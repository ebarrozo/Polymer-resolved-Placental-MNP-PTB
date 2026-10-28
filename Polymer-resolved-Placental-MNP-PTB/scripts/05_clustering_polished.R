#!/usr/bin/env Rscript
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 05_clustering.R
# CLR-based placental MNP composition clustering
# Main: UMAP + cluster enrichment + stacked composition/burden views
# Supplement: PCA, PCA scree, t-SNE, PERMANOVA/adonis tests
# ============================================================
# Input:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#   outputs/tables/mnp_ptb_cohort__polymer_detection_summary.csv
#
# Outputs:
#   Main-ish figures:
#     outputs/figs/fig_clustering__umap_by_delivery_class.png
#     outputs/figs/fig_clustering__umap_by_cluster.png
#     outputs/figs/fig_clustering__cluster_preterm_fraction.png
#     outputs/figs/fig_clustering__stacked_subject_burden_by_delivery_class.png
#     outputs/figs/fig_clustering__stacked_subject_burden_by_cluster.png
#     outputs/figs/fig_clustering__cluster_polymer_heatmap.png
#
#   Supplement figures:
#     outputs/figs/supp_clustering__pca_by_delivery_class.png
#     outputs/figs/supp_clustering__pca_by_cluster.png
#     outputs/figs/supp_clustering__pca_scree.png
#     outputs/figs/supp_clustering__tsne_by_delivery_class.png
#     outputs/figs/supp_clustering__tsne_by_cluster.png
#     outputs/figs/supp_clustering__cluster_burden_boxplot.png
#
#   Tables:
#     outputs/tables/clustering__polymer_matrix_used.csv
#     outputs/tables/clustering__sample_scores_and_clusters.csv
#     outputs/tables/clustering__cluster_summary.csv
#     outputs/tables/clustering__cluster_enrichment_delivery_class.csv
#     outputs/tables/clustering__composition_tests_delivery_class.csv
#     outputs/tables/clustering__k_selection_silhouette.csv
#     outputs/tables/clustering__clr_pca_outlier_audit.csv
#
# Notes:
#   - Clustering uses polymer composition, not absolute total burden.
#   - Sparse polymers are excluded from the clustering matrix by default.
#   - Initial PCA/CLR outliers are audited, then excluded from primary clustering by default.
#   - UMAP/t-SNE are run on CLR-scaled matrix after primary exclusions.
#   - If uwot or Rtsne is not installed, those plots are skipped with clear messages.
#   - PNG + CSV only.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(cluster)
})

source("00_config.R")
source(file.path("R", "utils_mnp.R"))
source(file.path("R", "figure_theme_mnp.R"))

# ----------------------------
# Settings
# ----------------------------
IN_CSV <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")
IN_DETECT <- file.path(DIR_TABLES, "mnp_ptb_cohort__polymer_detection_summary.csv")

K_FORCE <- NA_integer_       # set to 2L, 3L, etc. to override silhouette
K_MIN <- 2L
K_MAX <- 8L
KM_NSTART <- 100L

PSEUDO_PROP <- 1e-6
EXCLUDE_SPARSE_POLYMERS_FROM_CLUSTERING <- TRUE
EXCLUDE_TOTAL_ZERO_ROWS <- TRUE
EXCLUDE_CLR_PCA_OUTLIERS_FOR_PRIMARY_CLUSTERING <- TRUE
CLR_PCA_OUTLIER_Z <- 4

RUN_UMAP <- TRUE
RUN_TSNE <- TRUE
UMAP_N_NEIGHBORS <- 15L
UMAP_MIN_DIST <- 0.10
TSNE_PERPLEXITY_FRAC <- 0.15

BASE_SIZE <- 11
set.seed(1234)

# ----------------------------
# Outputs
# ----------------------------
OUT_MATRIX <- file.path(DIR_TABLES, "clustering__polymer_matrix_used_polished.csv")
OUT_SCORES <- file.path(DIR_TABLES, "clustering__sample_scores_and_clusters_polished.csv")
OUT_SUMMARY <- file.path(DIR_TABLES, "clustering__cluster_summary_polished.csv")
OUT_ENRICH <- file.path(DIR_TABLES, "clustering__cluster_enrichment_delivery_class_polished.csv")
OUT_TESTS <- file.path(DIR_TABLES, "clustering__composition_tests_delivery_class_polished.csv")
OUT_K <- file.path(DIR_TABLES, "clustering__k_selection_silhouette_polished.csv")
OUT_OUTLIERS <- file.path(DIR_TABLES, "clustering__clr_pca_outlier_audit_polished.csv")

OUT_UMAP_DELIVERY <- file.path(DIR_FIGS, "fig_clustering__umap_by_delivery_class_polished.png")
OUT_UMAP_CLUSTER <- file.path(DIR_FIGS, "fig_clustering__umap_by_cluster_polished.png")
OUT_PRETERM <- file.path(DIR_FIGS, "fig_clustering__cluster_preterm_fraction_polished.png")
OUT_STACK_DELIVERY <- file.path(DIR_FIGS, "fig_clustering__stacked_subject_burden_by_delivery_class_polished.png")
OUT_STACK_CLUSTER <- file.path(DIR_FIGS, "fig_clustering__stacked_subject_burden_by_cluster_polished.png")
OUT_HEATMAP <- file.path(DIR_FIGS, "fig_clustering__cluster_polymer_heatmap_polished.png")

OUT_SUPP_PCA_DELIVERY <- file.path(DIR_FIGS, "supp_clustering__pca_by_delivery_class_polished.png")
OUT_SUPP_PCA_CLUSTER <- file.path(DIR_FIGS, "supp_clustering__pca_by_cluster_polished.png")
OUT_SUPP_PCA_SCREE <- file.path(DIR_FIGS, "supp_clustering__pca_scree_polished.png")
OUT_SUPP_TSNE_DELIVERY <- file.path(DIR_FIGS, "supp_clustering__tsne_by_delivery_class_polished.png")
OUT_SUPP_TSNE_CLUSTER <- file.path(DIR_FIGS, "supp_clustering__tsne_by_cluster_polished.png")
OUT_SUPP_BURDEN <- file.path(DIR_FIGS, "supp_clustering__cluster_burden_boxplot_polished.png")

OUT_LOG <- file.path(DIR_LOGS, "05_clustering__run_log_polished.txt")

# ----------------------------
# Helpers
# ----------------------------
read_model_table <- function(path) {
  readr::read_csv(
    path,
    show_col_types = FALSE,
    guess_max = 10000,
    col_types = cols(
      SampleID = col_character(),
      PregnancyID = col_character(),
      `Pregnancy ID` = col_character(),
      home_zip_2nd_3rd_trim = col_character(),
      .default = col_guess()
    )
  ) %>%
    mutate(
      SampleID = as.character(SampleID),
      PregnancyID = as.character(PregnancyID),
      Delivery_class = factor(as.character(Delivery_class), levels = c("Term", "Preterm")),
      PTB01 = case_when(
        Delivery_class == "Preterm" ~ 1L,
        Delivery_class == "Term" ~ 0L,
        TRUE ~ NA_integer_
      ),
      GA = safe_num(GA),
      Cumulative_MNPs_ug_per_g = safe_num(Cumulative_MNPs_ug_per_g)
    )
}

clr_from_composition <- function(X, pseudo = 1e-6) {
  X[!is.finite(X)] <- 0
  row_tot <- rowSums(X, na.rm = TRUE)
  P <- sweep(X, 1, row_tot, "/")
  P[!is.finite(P)] <- 0
  P <- P + pseudo
  P <- sweep(P, 1, rowSums(P), "/")
  logP <- log(P)
  sweep(logP, 1, rowMeans(logP), "-")
}

avg_sil_width <- function(x_scaled, k, nstart = 100L) {
  km <- stats::kmeans(x_scaled, centers = k, nstart = nstart)
  sil <- cluster::silhouette(km$cluster, dist(x_scaled))
  mean(sil[, "sil_width"])
}

cluster_name_from_profile <- function(profile_named, top_n = 2) {
  vals <- as.numeric(profile_named)
  nms <- names(profile_named)
  nms <- str_remove(nms, "^mean_prop_")
  top <- nms[order(vals, decreasing = TRUE)][seq_len(min(top_n, length(vals)))]
  paste(top, collapse = "/")
}

safe_ellipse_data <- function(dat, group_col) {
  dat %>%
    group_by(.data[[group_col]]) %>%
    filter(n() >= 3) %>%
    ungroup()
}

format_p_title <- function(p) {
  if (is.na(p)) return("p=NA")
  if (p < 0.001) return("p<0.001")
  paste0("p=", formatC(p, format = "f", digits = 3))
}

# ----------------------------
# Load
# ----------------------------
df0 <- read_model_table(IN_CSV)
df <- df0 %>% filter(!is.na(Delivery_class))

if (!file.exists(IN_DETECT)) {
  warning("Detection summary not found. Sparse filtering will use observed >0 counts in model table.")
  detect_tbl <- tibble(
    Polymer = POLYMERS,
    n_detected_gt0 = vapply(POLYMERS, function(p) {
      col <- legacy_col_ug_g(p)
      sum(is.finite(safe_num(df[[col]])) & safe_num(df[[col]]) > 0)
    }, numeric(1)),
    continuous_model_ok = n_detected_gt0 >= MIN_DETECTED_N_FOR_CONTINUOUS_MODEL
  )
} else {
  detect_tbl <- readr::read_csv(IN_DETECT, show_col_types = FALSE) %>%
    mutate(Polymer = toupper(as.character(Polymer)))
}

poly_keep <- POLYMERS
if (EXCLUDE_SPARSE_POLYMERS_FROM_CLUSTERING && "continuous_model_ok" %in% names(detect_tbl)) {
  sparse <- detect_tbl %>%
    filter(Polymer %in% POLYMERS, continuous_model_ok %in% FALSE) %>%
    pull(Polymer)
  poly_keep <- setdiff(poly_keep, sparse)
}

if (length(poly_keep) < 3) {
  stop("Fewer than 3 polymers remain for clustering after sparse filtering: ", paste(poly_keep, collapse = ", "), call. = FALSE)
}

ug_cols <- legacy_col_ug_g(poly_keep)
missing_ug <- setdiff(ug_cols, names(df))
if (length(missing_ug) > 0) stop("Missing ug/g columns: ", paste(missing_ug, collapse = ", "), call. = FALSE)

X_raw0 <- df %>%
  select(all_of(ug_cols)) %>%
  mutate(across(everything(), safe_num)) %>%
  as.matrix()

colnames(X_raw0) <- poly_keep
X_raw0[!is.finite(X_raw0)] <- 0

row_total0 <- rowSums(X_raw0, na.rm = TRUE)
keep_rows <- is.finite(row_total0)
if (EXCLUDE_TOTAL_ZERO_ROWS) keep_rows <- keep_rows & row_total0 > 0

df_initial <- df[keep_rows, , drop = FALSE]
X_raw_initial <- X_raw0[keep_rows, , drop = FALSE]

if (nrow(X_raw_initial) < 10) stop("Too few rows remain for clustering.", call. = FALSE)

# ----------------------------
# Initial PCA outlier audit
# ----------------------------
P_initial <- sweep(X_raw_initial, 1, rowSums(X_raw_initial), "/")
P_initial[!is.finite(P_initial)] <- 0
X_clr_initial <- clr_from_composition(X_raw_initial, pseudo = PSEUDO_PROP)
X_scaled_initial <- scale(X_clr_initial)

pca_initial <- prcomp(X_scaled_initial, center = FALSE, scale. = FALSE)
pca_var_initial <- (pca_initial$sdev^2) / sum(pca_initial$sdev^2)

initial_scores <- as_tibble(pca_initial$x[, 1:min(5, ncol(pca_initial$x)), drop = FALSE], .name_repair = "minimal") %>%
  rename_with(~ paste0("PC", seq_along(.x))) %>%
  bind_cols(
    df_initial %>%
      transmute(
        SampleID,
        PregnancyID,
        Delivery_class,
        PTB01,
        GA,
        Cumulative_MNPs_ug_per_g
      )
  ) %>%
  mutate(
    PC1_z = as.numeric(scale(PC1)),
    PC2_z = as.numeric(scale(PC2)),
    clr_pca_outlier = abs(PC1_z) > CLR_PCA_OUTLIER_Z | abs(PC2_z) > CLR_PCA_OUTLIER_Z
  )

outlier_audit <- initial_scores %>%
  filter(clr_pca_outlier) %>%
  select(
    SampleID,
    PregnancyID,
    Delivery_class,
    GA,
    Cumulative_MNPs_ug_per_g,
    PC1,
    PC2,
    PC1_z,
    PC2_z
  )

write_csv_clean(outlier_audit, OUT_OUTLIERS)
message("CLR/PCA outliers flagged: ", nrow(outlier_audit))

# ----------------------------
# Primary clustering data after optional outlier exclusion
# ----------------------------
keep_primary <- rep(TRUE, nrow(df_initial))
if (EXCLUDE_CLR_PCA_OUTLIERS_FOR_PRIMARY_CLUSTERING) {
  keep_primary <- !initial_scores$clr_pca_outlier
}

df_clust <- df_initial[keep_primary, , drop = FALSE]
X_raw <- X_raw_initial[keep_primary, , drop = FALSE]

P <- sweep(X_raw, 1, rowSums(X_raw), "/")
P[!is.finite(P)] <- 0
X_clr <- clr_from_composition(X_raw, pseudo = PSEUDO_PROP)
X_scaled <- scale(X_clr)

# PCA on primary matrix
pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as_tibble(pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE], .name_repair = "minimal") %>%
  rename_with(~ paste0("PC", seq_along(.x))) %>%
  bind_cols(
    df_clust %>%
      transmute(
        SampleID,
        PregnancyID,
        Delivery_class,
        PTB01,
        GA,
        Cumulative_MNPs_ug_per_g,
        SDI_score = if ("SDI_score" %in% names(.)) safe_num(SDI_score) else NA_real_,
        SVI_score = if ("SVI_score" %in% names(.)) safe_num(SVI_score) else NA_real_,
        Hospital = if ("Hospital" %in% names(.)) as.character(Hospital) else NA_character_,
        Ethnicity_Race = if ("Ethnicity_Race" %in% names(.)) as.character(Ethnicity_Race) else NA_character_
      )
  )

# ----------------------------
# K selection and clustering
# ----------------------------
max_k_allowed <- min(K_MAX, nrow(X_scaled) - 1)
k_grid <- K_MIN:max_k_allowed

k_tbl <- purrr::map_dfr(k_grid, function(k) {
  tibble(k = k, avg_silhouette = avg_sil_width(X_scaled, k, nstart = KM_NSTART))
})

if (is.na(K_FORCE)) {
  K_CHOSEN <- k_tbl$k[which.max(k_tbl$avg_silhouette)]
} else {
  K_CHOSEN <- as.integer(K_FORCE)
}

km <- kmeans(X_scaled, centers = K_CHOSEN, nstart = KM_NSTART)

scores <- scores %>%
  mutate(
    Cluster = factor(km$cluster),
    Cluster_chr = as.character(Cluster)
  )

# ----------------------------
# Cluster profiles / names
# ----------------------------
prop_tbl <- as_tibble(P, .name_repair = "minimal") %>%
  setNames(poly_keep) %>%
  mutate(SampleID = scores$SampleID, Cluster = scores$Cluster)

cluster_profile <- prop_tbl %>%
  group_by(Cluster) %>%
  summarise(
    n_profile = n(),
    across(all_of(poly_keep), ~ mean(.x, na.rm = TRUE), .names = "mean_prop_{.col}"),
    .groups = "drop"
  )

prof_cols <- names(cluster_profile)[stringr::str_detect(names(cluster_profile), "^mean_prop_")]

cluster_names <- cluster_profile %>%
  dplyr::mutate(
    Cluster_Name = purrr::pmap_chr(
      dplyr::select(., dplyr::all_of(prof_cols)),
      function(...) {
        vals <- c(...)
        names(vals) <- prof_cols
        cluster_id <- as.character(dplyr::cur_data_all()$Cluster[1])
        # cluster_id is not reliable inside pmap, so rebuilt below after row index join
        cluster_name_from_profile(vals)
      }
    )
  ) %>%
  dplyr::mutate(
    Cluster_Name = paste0("Cluster ", as.character(Cluster), ": ", Cluster_Name)
  ) %>%
  dplyr::select(Cluster, Cluster_Name)

scores <- scores %>%
  left_join(cluster_names, by = "Cluster")

ellipse_delivery_df <- safe_ellipse_data(scores, "Delivery_class")
ellipse_cluster_df <- safe_ellipse_data(scores, "Cluster")

# ----------------------------
# UMAP and t-SNE
# ----------------------------
umap_ran <- FALSE
if (RUN_UMAP && requireNamespace("uwot", quietly = TRUE)) {
  set.seed(1234)
  um <- uwot::umap(
    X_scaled,
    n_neighbors = min(UMAP_N_NEIGHBORS, nrow(X_scaled) - 1),
    min_dist = UMAP_MIN_DIST,
    metric = "euclidean",
    scale = FALSE,
    ret_model = FALSE,
    verbose = FALSE
  )
  scores$UMAP1 <- um[, 1]
  scores$UMAP2 <- um[, 2]
  umap_ran <- TRUE
} else {
  message("UMAP skipped. Install package 'uwot' to generate UMAP panels.")
  scores$UMAP1 <- NA_real_
  scores$UMAP2 <- NA_real_
}

tsne_ran <- FALSE
if (RUN_TSNE && requireNamespace("Rtsne", quietly = TRUE)) {
  set.seed(1234)
  perplexity <- max(5, floor(nrow(X_scaled) * TSNE_PERPLEXITY_FRAC))
  perplexity <- min(perplexity, floor((nrow(X_scaled) - 1) / 3))
  ts <- Rtsne::Rtsne(
    X_scaled,
    dims = 2,
    perplexity = perplexity,
    check_duplicates = FALSE,
    pca = FALSE,
    verbose = FALSE
  )
  scores$tSNE1 <- ts$Y[, 1]
  scores$tSNE2 <- ts$Y[, 2]
  tsne_ran <- TRUE
} else {
  message("t-SNE skipped. Install package 'Rtsne' to generate t-SNE panels.")
  scores$tSNE1 <- NA_real_
  scores$tSNE2 <- NA_real_
}

# ----------------------------
# Cluster summary and enrichment
# ----------------------------
cluster_summary <- scores %>%
  group_by(Cluster, Cluster_Name) %>%
  summarise(
    n = n(),
    n_preterm = sum(Delivery_class == "Preterm", na.rm = TRUE),
    n_term = sum(Delivery_class == "Term", na.rm = TRUE),
    pct_preterm = 100 * n_preterm / n,
    GA_median = median(GA, na.rm = TRUE),
    GA_q1 = quantile(GA, 0.25, na.rm = TRUE, names = FALSE),
    GA_q3 = quantile(GA, 0.75, na.rm = TRUE, names = FALSE),
    total_burden_median_ug_g = median(Cumulative_MNPs_ug_per_g, na.rm = TRUE),
    total_burden_q1_ug_g = quantile(Cumulative_MNPs_ug_per_g, 0.25, na.rm = TRUE, names = FALSE),
    total_burden_q3_ug_g = quantile(Cumulative_MNPs_ug_per_g, 0.75, na.rm = TRUE, names = FALSE),
    SDI_median = median(SDI_score, na.rm = TRUE),
    SVI_median = median(SVI_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(cluster_profile, by = "Cluster")

cluster_enrich <- purrr::map_dfr(levels(scores$Cluster), function(cl) {
  in_idx <- scores$Cluster == cl
  tab <- table(
    Cluster = factor(ifelse(in_idx, "In", "Out"), levels = c("In", "Out")),
    Delivery = scores$Delivery_class
  )

  ft <- suppressWarnings(fisher.test(tab))

  tibble(
    Cluster = cl,
    Cluster_Name = unique(scores$Cluster_Name[scores$Cluster == cl]),
    n_in = sum(in_idx),
    preterm_in = sum(scores$Delivery_class[in_idx] == "Preterm", na.rm = TRUE),
    term_in = sum(scores$Delivery_class[in_idx] == "Term", na.rm = TRUE),
    preterm_out = sum(scores$Delivery_class[!in_idx] == "Preterm", na.rm = TRUE),
    term_out = sum(scores$Delivery_class[!in_idx] == "Term", na.rm = TRUE),
    pct_preterm_in = 100 * preterm_in / n_in,
    OR_preterm = unname(ft$estimate),
    p_value = unname(ft$p.value)
  )
}) %>%
  mutate(q_BH = p.adjust(p_value, method = "BH"))

# ----------------------------
# Composition tests: Term vs Preterm
# ----------------------------
composition_tests <- tibble(
  test = character(),
  matrix = character(),
  n = integer(),
  df = numeric(),
  statistic = numeric(),
  R2 = numeric(),
  p_value = numeric(),
  note = character()
)

if (requireNamespace("vegan", quietly = TRUE)) {
  ad <- vegan::adonis2(dist(X_scaled) ~ scores$Delivery_class, permutations = 999)
  composition_tests <- bind_rows(
    composition_tests,
    tibble(
      test = "PERMANOVA/adonis2",
      matrix = "CLR-scaled Euclidean distance",
      n = nrow(X_scaled),
      df = ad$Df[1],
      statistic = ad$F[1],
      R2 = ad$R2[1],
      p_value = ad$`Pr(>F)`[1],
      note = "Tests whether MNP composition differs between Term and Preterm."
    )
  )

  bd <- vegan::betadisper(dist(X_scaled), scores$Delivery_class)
  bd_perm <- vegan::permutest(bd, permutations = 999)
  composition_tests <- bind_rows(
    composition_tests,
    tibble(
      test = "PERMDISP/betadisper",
      matrix = "CLR-scaled Euclidean distance",
      n = nrow(X_scaled),
      df = bd_perm$tab$Df[1],
      statistic = bd_perm$tab$F[1],
      R2 = NA_real_,
      p_value = bd_perm$tab$`Pr(>F)`[1],
      note = "Tests whether within-group dispersion differs; helps interpret PERMANOVA."
    )
  )
} else {
  composition_tests <- bind_rows(
    composition_tests,
    tibble(
      test = "PERMANOVA/adonis2",
      matrix = "CLR-scaled Euclidean distance",
      n = nrow(X_scaled),
      df = NA_real_,
      statistic = NA_real_,
      R2 = NA_real_,
      p_value = NA_real_,
      note = "Skipped because package 'vegan' is not installed."
    )
  )
}

permanova_p <- composition_tests %>%
  filter(test == "PERMANOVA/adonis2") %>%
  pull(p_value) %>%
  .[1]

permanova_R2 <- composition_tests %>%
  filter(test == "PERMANOVA/adonis2") %>%
  pull(R2) %>%
  .[1]

# ----------------------------
# Write tables
# ----------------------------
matrix_out <- bind_cols(
  scores %>% select(SampleID, PregnancyID, Delivery_class, Cluster, Cluster_Name),
  as_tibble(X_raw, .name_repair = "minimal") %>% setNames(paste0(poly_keep, "_ug_per_g_used")),
  as_tibble(P, .name_repair = "minimal") %>% setNames(paste0(poly_keep, "_proportion")),
  as_tibble(X_clr, .name_repair = "minimal") %>% setNames(paste0(poly_keep, "_clr"))
)

write_csv_clean(matrix_out, OUT_MATRIX)
write_csv_clean(scores, OUT_SCORES)
write_csv_clean(cluster_summary, OUT_SUMMARY)
write_csv_clean(cluster_enrich, OUT_ENRICH)
write_csv_clean(composition_tests, OUT_TESTS)
write_csv_clean(k_tbl, OUT_K)

# Compatibility outputs for old scripts
saveRDS(scores, file.path(PROJECT_ROOT, "df_clust.rds"))
write_csv_clean(scores, file.path(DIR_TABLES, "df_clust_polished.csv"))

# ----------------------------
# Plot inputs
# ----------------------------
poly_pal <- make_poly_palette(poly_keep, mode = "familyish")
poly_full <- paste0(POLYMER_FULL_NAMES[poly_keep], " (", poly_keep, ")")

stack_long <- as_tibble(X_raw, .name_repair = "minimal") %>%
  setNames(poly_keep) %>%
  mutate(
    SampleID = scores$SampleID,
    Delivery_class = scores$Delivery_class,
    Cluster = scores$Cluster,
    Cluster_Name = scores$Cluster_Name,
    Cumulative_MNPs_ug_per_g = scores$Cumulative_MNPs_ug_per_g
  ) %>%
  pivot_longer(all_of(poly_keep), names_to = "Polymer", values_to = "ug_g") %>%
  mutate(
    Polymer = factor(Polymer, levels = poly_keep),
    PolymerLabel = factor(paste0(POLYMER_FULL_NAMES[as.character(Polymer)], " (", as.character(Polymer), ")"),
                          levels = poly_full)
  )

subject_order_delivery <- stack_long %>%
  distinct(SampleID, Delivery_class, Cumulative_MNPs_ug_per_g) %>%
  arrange(Delivery_class, desc(Cumulative_MNPs_ug_per_g)) %>%
  mutate(SampleID_ordered = factor(SampleID, levels = unique(SampleID)))

subject_order_cluster <- stack_long %>%
  distinct(SampleID, Cluster_Name, Cumulative_MNPs_ug_per_g) %>%
  arrange(Cluster_Name, desc(Cumulative_MNPs_ug_per_g)) %>%
  mutate(SampleID_ordered = factor(SampleID, levels = unique(SampleID)))

stack_delivery <- stack_long %>%
  left_join(subject_order_delivery %>% select(SampleID, SampleID_ordered), by = "SampleID")

stack_cluster <- stack_long %>%
  left_join(subject_order_cluster %>% select(SampleID, SampleID_ordered), by = "SampleID")

# ----------------------------
# Main dimensional reduction: UMAP
# ----------------------------
if (umap_ran) {
  umap_subtitle <- paste0(
    "PERMANOVA on CLR composition: R²=",
    ifelse(is.na(permanova_R2), "NA", formatC(permanova_R2, digits = 3, format = "f")),
    ", ",
    format_p_title(permanova_p)
  )

  p_umap_delivery <- ggplot(scores, aes(x = UMAP1, y = UMAP2, shape = Delivery_class)) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(
      title = "UMAP of Placental MNP Composition",
      subtitle = umap_subtitle,
      x = "UMAP1",
      y = "UMAP2",
      shape = "Delivery class"
    ) +
    theme_mnp(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  ggsave_mnp(OUT_UMAP_DELIVERY, p_umap_delivery, width = 7.2, height = 5.6, dpi = 300)

  p_umap_cluster <- ggplot(scores, aes(x = UMAP1, y = UMAP2, shape = Cluster)) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(
      title = "UMAP of Placental MNP Composition by Cluster",
      subtitle = paste0("K = ", K_CHOSEN, "; selected by average silhouette unless overridden"),
      x = "UMAP1",
      y = "UMAP2",
      shape = "Cluster"
    ) +
    theme_mnp(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  ggsave_mnp(OUT_UMAP_CLUSTER, p_umap_cluster, width = 7.2, height = 5.6, dpi = 300)
}

# ----------------------------
# Supplement: PCA and scree
# ----------------------------
pc_labs <- c(
  PC1 = paste0("PC1 (", scales::percent(pca_var[1], accuracy = 0.1), ")"),
  PC2 = paste0("PC2 (", scales::percent(pca_var[2], accuracy = 0.1), ")")
)

p_pca_delivery <- ggplot(scores, aes(x = PC1, y = PC2, shape = Delivery_class)) +
  geom_point(size = 2.5, alpha = 0.80) +
  stat_ellipse(
    data = ellipse_delivery_df,
    aes(group = Delivery_class),
    linewidth = 0.5,
    level = 0.68,
    na.rm = TRUE
  ) +
  labs(
    title = "CLR-PCA of Placental MNP Composition",
    subtitle = paste0("Polymers used: ", paste(poly_keep, collapse = ", ")),
    x = pc_labs["PC1"],
    y = pc_labs["PC2"],
    shape = "Delivery class"
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave_mnp(OUT_SUPP_PCA_DELIVERY, p_pca_delivery, width = 7.2, height = 5.6, dpi = 300)

p_pca_cluster <- ggplot(scores, aes(x = PC1, y = PC2, shape = Cluster)) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(
    data = ellipse_cluster_df,
    aes(group = Cluster),
    linewidth = 0.5,
    level = 0.68,
    na.rm = TRUE
  ) +
  labs(
    title = "CLR-PCA of Placental MNP Composition by K-Means Cluster",
    subtitle = paste0("K = ", K_CHOSEN),
    x = pc_labs["PC1"],
    y = pc_labs["PC2"],
    shape = "Cluster"
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave_mnp(OUT_SUPP_PCA_CLUSTER, p_pca_cluster, width = 7.2, height = 5.6, dpi = 300)

scree <- tibble(
  PC = seq_along(pca_var),
  variance_explained = pca_var,
  cumulative_variance = cumsum(pca_var)
) %>%
  filter(PC <= min(12, n()))

p_scree <- ggplot(scree, aes(x = PC, y = variance_explained)) +
  geom_col(width = 0.75, alpha = 0.70) +
  geom_line(aes(y = cumulative_variance), linewidth = 0.7) +
  geom_point(aes(y = cumulative_variance), size = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = scree$PC) +
  labs(
    title = "PCA Scree Plot for CLR-Transformed MNP Composition",
    subtitle = "Bars = variance explained; line = cumulative variance",
    x = "Principal component",
    y = "Variance explained"
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave_mnp(OUT_SUPP_PCA_SCREE, p_scree, width = 7.2, height = 5.2, dpi = 300)

# ----------------------------
# Supplement: t-SNE
# ----------------------------
if (tsne_ran) {
  p_tsne_delivery <- ggplot(scores, aes(x = tSNE1, y = tSNE2, shape = Delivery_class)) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(
      title = "t-SNE of Placental MNP Composition",
      x = "tSNE1",
      y = "tSNE2",
      shape = "Delivery class"
    ) +
    theme_mnp(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  ggsave_mnp(OUT_SUPP_TSNE_DELIVERY, p_tsne_delivery, width = 7.2, height = 5.6, dpi = 300)

  p_tsne_cluster <- ggplot(scores, aes(x = tSNE1, y = tSNE2, shape = Cluster)) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(
      title = "t-SNE of Placental MNP Composition by Cluster",
      x = "tSNE1",
      y = "tSNE2",
      shape = "Cluster"
    ) +
    theme_mnp(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  ggsave_mnp(OUT_SUPP_TSNE_CLUSTER, p_tsne_cluster, width = 7.2, height = 5.6, dpi = 300)
}

# ----------------------------
# Cluster heatmap
# ----------------------------
heat_df <- cluster_profile %>%
  select(Cluster, starts_with("mean_prop_")) %>%
  pivot_longer(starts_with("mean_prop_"), names_to = "Polymer", values_to = "mean_prop") %>%
  mutate(
    Polymer = str_remove(Polymer, "^mean_prop_"),
    Polymer = factor(Polymer, levels = poly_keep),
    Cluster = factor(Cluster)
  ) %>%
  left_join(cluster_names, by = "Cluster")

p_heat <- ggplot(heat_df, aes(x = Polymer, y = Cluster_Name, fill = mean_prop)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Mean Polymer Composition by Cluster",
    x = NULL,
    y = NULL,
    fill = "Mean proportion"
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave_mnp(OUT_HEATMAP, p_heat, width = 8.5, height = 4.8, dpi = 300)

# ----------------------------
# Cluster enrichment / preterm fraction
# ----------------------------
p_preterm_df <- cluster_summary %>%
  mutate(
    n = as.integer(n),
    n_preterm = as.integer(n_preterm),
    Cluster_Name = factor(Cluster_Name, levels = Cluster_Name[order(pct_preterm)]),
    preterm_label = paste0(n_preterm, "/", n)
  )

p_preterm <- ggplot(p_preterm_df, aes(x = Cluster_Name, y = pct_preterm)) +
  geom_col(width = 0.75, alpha = 0.75) +
  geom_text(aes(label = preterm_label), hjust = -0.10, size = 3.3) +
  coord_flip() +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    limits = c(0, max(100, max(p_preterm_df$pct_preterm, na.rm = TRUE) * 1.12))
  ) +
  labs(
    title = "Preterm Fraction by MNP Composition Cluster",
    subtitle = "Labels show preterm/total within cluster",
    x = NULL,
    y = "Preterm deliveries"
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave_mnp(OUT_PRETERM, p_preterm, width = 7.2, height = 4.8, dpi = 300)

# ----------------------------
# Stacked subject burden plots
# ----------------------------
p_stack_delivery <- ggplot(stack_delivery, aes(x = SampleID_ordered, y = ug_g, fill = Polymer)) +
  geom_col(width = 0.85) +
  facet_grid(. ~ Delivery_class, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = poly_pal, drop = FALSE) +
  scale_y_continuous(labels = label_number()) +
  labs(
    title = "Placental MNP Burden by Subject",
    subtitle = paste0(
      "Subjects ordered within delivery class by cumulative burden; PERMANOVA: R²=",
      ifelse(is.na(permanova_R2), "NA", formatC(permanova_R2, digits = 3, format = "f")),
      ", ",
      format_p_title(permanova_p)
    ),
    x = "Subject",
    y = expression("MNP burden ("*mu*"g/g tissue)"),
    fill = "Polymer"
  ) +
  theme_mnp(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave_mnp(OUT_STACK_DELIVERY, p_stack_delivery, width = 13.5, height = 6.2, dpi = 300)

p_stack_cluster <- ggplot(stack_cluster, aes(x = SampleID_ordered, y = ug_g, fill = Polymer)) +
  geom_col(width = 0.85) +
  facet_grid(. ~ Cluster_Name, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = poly_pal, drop = FALSE) +
  scale_y_continuous(labels = label_number()) +
  labs(
    title = "Placental MNP Burden by Composition Cluster",
    subtitle = "Subjects ordered within cluster by cumulative burden",
    x = "Subject",
    y = expression("MNP burden ("*mu*"g/g tissue)"),
    fill = "Polymer"
  ) +
  theme_mnp(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave_mnp(OUT_STACK_CLUSTER, p_stack_cluster, width = 13.5, height = 6.2, dpi = 300)

# Supplement burden boxplot
scores_plot <- scores %>%
  mutate(Cluster_Name = factor(Cluster_Name, levels = unique(cluster_summary$Cluster_Name)))

p_burden <- ggplot(scores_plot, aes(x = Cluster_Name, y = Cumulative_MNPs_ug_per_g)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.55, size = 1.6) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), labels = label_number()) +
  coord_flip() +
  labs(
    title = "Total MNP Burden by Composition Cluster",
    x = NULL,
    y = expression("Cumulative MNP burden ("*mu*"g/g tissue)")
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave_mnp(OUT_SUPP_BURDEN, p_burden, width = 7.2, height = 4.8, dpi = 300)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "05_clustering.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Input rows with Delivery_class: ", nrow(df)),
  paste0("Initial rows before CLR/PCA outlier exclusion: ", nrow(df_initial)),
  paste0("CLR/PCA outliers flagged: ", nrow(outlier_audit)),
  paste0("Rows used for primary clustering: ", nrow(scores)),
  paste0("Polymers used: ", paste(poly_keep, collapse = ", ")),
  paste0("Sparse filtering enabled: ", EXCLUDE_SPARSE_POLYMERS_FROM_CLUSTERING),
  paste0("Outliers excluded for primary clustering: ", EXCLUDE_CLR_PCA_OUTLIERS_FOR_PRIMARY_CLUSTERING),
  paste0("UMAP generated: ", umap_ran),
  paste0("t-SNE generated: ", tsne_ran),
  paste0("K chosen: ", K_CHOSEN),
  "",
  "K selection:",
  capture.output(print(k_tbl)),
  "",
  "Cluster summary:",
  capture.output(print(cluster_summary)),
  "",
  "Cluster enrichment:",
  capture.output(print(cluster_enrich)),
  "",
  "Composition tests:",
  capture.output(print(composition_tests)),
  "",
  "Wrote:",
  OUT_MATRIX,
  OUT_SCORES,
  OUT_SUMMARY,
  OUT_ENRICH,
  OUT_TESTS,
  OUT_K,
  OUT_OUTLIERS,
  OUT_UMAP_DELIVERY,
  OUT_UMAP_CLUSTER,
  OUT_PRETERM,
  OUT_STACK_DELIVERY,
  OUT_STACK_CLUSTER,
  OUT_HEATMAP,
  OUT_SUPP_PCA_DELIVERY,
  OUT_SUPP_PCA_CLUSTER,
  OUT_SUPP_PCA_SCREE,
  OUT_SUPP_TSNE_DELIVERY,
  OUT_SUPP_TSNE_CLUSTER,
  OUT_SUPP_BURDEN
)

cat(paste(run_log, collapse = "\n"), file = OUT_LOG)
message(paste(run_log, collapse = "\n"))
