# ============================================================
# 19_figure5_analytical_validation_panels.R
# Figure 5 panels for PTB-MNP placenta manuscript
#
# Creates four separate high-resolution PNGs:
#   Fig5a_validation_workflow.png
#   Fig5b_whole_process_recovery.png
#   Fig5c_PU_PET_QC_categories_by_delivery_class.png
#   Fig5d_PU_PET_analytical_robustness.png
#
# Run from PROJECT_ROOT:
#   source("scripts/19_figure5_analytical_validation_panels.R")
#
# Assumes 00_config.R and outputs from
# 18_placenta_spike_recovery_validation.R are present.
# ============================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(forcats)
})

# ------------------------------------------------------------
# USER OPTION
# ------------------------------------------------------------
# Control provenance confirmed by the analytical laboratory:
# control block A corresponds to the placenta validation run.
SPIKE_ONLY_CONTROL_BLOCK <- "A"

# ------------------------------------------------------------
# Locate latest validation output folder
# ------------------------------------------------------------
validation_dirs <- list.dirs(
  DIR_TABLES,
  full.names = TRUE,
  recursive = FALSE
)

validation_dirs <- validation_dirs[
  grepl("^18_placenta_spike_recovery_validation__", basename(validation_dirs))
]

if (length(validation_dirs) == 0) {
  stop(
    "No validation output directory found under ", DIR_TABLES,
    ". Run scripts/18_placenta_spike_recovery_validation.R first.",
    call. = FALSE
  )
}

# Prefer the most recently modified validation output folder.
validation_dir <- validation_dirs[
  order(file.info(validation_dirs)$mtime, decreasing = TRUE)
][1]

message("Using validation output directory: ", validation_dir)

fig_dir <- file.path(DIR_FIGS, "Figure5_analytical_validation")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

req <- c(
  "spike_only_process_recovery.csv",
  "placenta_paired_polymer_recovery.csv",
  "pu_pet_qc_category_counts_by_delivery_class.csv",
  "pu_pet_qc_binary_ptb_sensitivity_models.csv",
  "pu_pet_qc_bias_reduced_models.csv",
  "pu_pet_batchmatched_maxblank_subtracted_continuous_models.csv"
)

missing_req <- req[!file.exists(file.path(validation_dir, req))]
if (length(missing_req) > 0) {
  stop(
    "Required validation outputs are missing:\n  ",
    paste(missing_req, collapse = "\n  "),
    "\nRe-run the current version of scripts/18_placenta_spike_recovery_validation.R.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Shared plotting settings
# ------------------------------------------------------------
poly_order <- c("PE", "PP", "PS", "ABS", "SBR", "PMMA",
                "PC", "PVC", "PU", "PET", "N6", "N66")

theme_nc <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 9, hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(8, 8, 8, 8)
  )

save_panel <- function(plot, filename, width, height) {
  out <- file.path(fig_dir, filename)
  ggsave(
    filename = out,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  message("Wrote: ", out)
}

# ============================================================
# PANEL A — validation workflow schematic
# ============================================================

# Built entirely with ggplot2 so it exports as a standalone PNG
# and remains reproducible.

box_df <- tibble::tribble(
  ~xmin, ~xmax, ~ymin, ~ymax, ~label,
   0.2,  1.8,  3.0,  3.8, "Fresh placenta",
   2.4,  4.0,  4.0,  4.8, "Unspiked aliquot",
   2.4,  4.0,  2.0,  2.8, "Spiked aliquot",
      4.55,  7.45,  2.90,  3.90, "Digestion → centrifugation
→ secondary processing
→ transfer",
   7.8,  9.2,  3.0,  3.8, "Py-GC/MS",
   9.9, 12.3,  3.0,  3.8, "Paired whole-process\nrecovery",
   0.2,  1.8,  0.5,  1.3, "Spike-only control",
   2.8,  5.2,  0.5,  1.3, "Same workflow\nwithout tissue",
   6.0,  8.0,  0.5,  1.3, "Py-GC/MS",
   8.8, 12.3,  0.5,  1.3, "Baseline whole-process\nrecovery"
)

p5a <- ggplot() +
  geom_rect(
    data = box_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white", color = "black", linewidth = 0.6
  ) +
  geom_text(
    data = box_df %>% filter(!grepl("^Digestion", label)),
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 3.0, lineheight = 0.95
  ) +
  geom_text(
    data = box_df %>% filter(grepl("^Digestion", label)),
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 2.65, lineheight = 0.92
  ) +
  # Fresh placenta to paired aliquots
  geom_segment(
    aes(x = 1.8, y = 3.4, xend = 2.25, yend = 4.4),
    arrow = arrow(length = unit(0.11, "inches")), linewidth = 0.5
  ) +
  geom_segment(
    aes(x = 1.8, y = 3.4, xend = 2.25, yend = 2.4),
    arrow = arrow(length = unit(0.11, "inches")), linewidth = 0.5
  ) +
  # Paired aliquots to identical processing
  geom_segment(
    aes(x = 4.0, y = 4.4, xend = 4.42, yend = 3.65),
    arrow = arrow(length = unit(0.11, "inches")), linewidth = 0.5
  ) +
  geom_segment(
    aes(x = 4.0, y = 2.4, xend = 4.42, yend = 3.15),
    arrow = arrow(length = unit(0.11, "inches")), linewidth = 0.5
  ) +
  # Processing to Py-GC/MS to recovery
  geom_segment(
    aes(x = 7.45, y = 3.4, xend = 7.7, yend = 3.4),
    arrow = arrow(length = unit(0.11, "inches")), linewidth = 0.5
  ) +
  geom_segment(
    aes(x = 9.2, y = 3.4, xend = 9.8, yend = 3.4),
    arrow = arrow(length = unit(0.11, "inches")), linewidth = 0.5
  ) +
  # Spike-only lane
  geom_segment(
    aes(x = 1.8, y = 0.9, xend = 2.7, yend = 0.9),
    arrow = arrow(length = unit(0.11, "inches")), linewidth = 0.5
  ) +
  geom_segment(
    aes(x = 5.2, y = 0.9, xend = 5.9, yend = 0.9),
    arrow = arrow(length = unit(0.11, "inches")), linewidth = 0.5
  ) +
  geom_segment(
    aes(x = 8.0, y = 0.9, xend = 8.7, yend = 0.9),
    arrow = arrow(length = unit(0.11, "inches")), linewidth = 0.5
  ) +
  annotate(
    "text", x = 3.2, y = 1.55,
    label = "12-polymer standard added before digestion",
    size = 2.9, fontface = "italic"
  ) +
  coord_cartesian(xlim = c(0, 12.6), ylim = c(0.2, 5.0), clip = "off") +
  labs(title = "Validation design") +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.margin = margin(8, 8, 8, 8)
  )

save_panel(
  p5a,
  "Fig5a_validation_workflow.png",
  width = 10.2,
  height = 4.0
)

# ============================================================
# PANEL B — whole-process recovery with/without placenta
# ============================================================

spike_only <- read_csv(
  file.path(validation_dir, "spike_only_process_recovery.csv"),
  show_col_types = FALSE
)

placenta_rec <- read_csv(
  file.path(validation_dir, "placenta_paired_polymer_recovery.csv"),
  show_col_types = FALSE
) %>%
  filter(paired_unspiked_available)

spike_only <- spike_only %>%
  filter(ControlBlock == SPIKE_ONLY_CONTROL_BLOCK)

if (nrow(spike_only) == 0) {
  stop("No placenta-run spike-only controls found in block A.")
}

control_label <- "Spike-only control (placenta run)"

b_spike <- spike_only %>%
  transmute(
    Polymer,
    Sample = paste(ControlBlock, ControlLabel, sep = ": "),
    Recovery = process_recovery_pct,
    Group = control_label
  )

b_plac <- placenta_rec %>%
  transmute(
    Polymer,
    Sample = BaseSampleID,
    Recovery = apparent_recovery_pct,
    Group = "Placenta"
  )

b_dat <- bind_rows(b_spike, b_plac) %>%
  filter(Polymer %in% poly_order, is.finite(Recovery)) %>%
  mutate(
    Polymer = factor(Polymer, levels = poly_order),
    Group = factor(Group, levels = c(control_label, "Placenta"))
  )

b_sum <- b_dat %>%
  group_by(Polymer, Group) %>%
  summarise(
    median = median(Recovery, na.rm = TRUE),
    q1 = quantile(Recovery, 0.25, na.rm = TRUE),
    q3 = quantile(Recovery, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p5b <- ggplot(
  b_dat,
  aes(x = Polymer, y = Recovery, shape = Group)
) +
  geom_hline(yintercept = 100, linetype = 2, linewidth = 0.5) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.50
    ),
    size = 2.0,
    alpha = 0.70
  ) +
  geom_errorbar(
    data = b_sum,
    aes(x = Polymer, ymin = q1, ymax = q3, group = Group),
    position = position_dodge(width = 0.50),
    width = 0.15,
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = b_sum,
    aes(x = Polymer, y = median, shape = Group, group = Group),
    position = position_dodge(width = 0.50),
    size = 2.5,
    fill = "white",
    stroke = 0.9
  ) +
  scale_shape_manual(values = c(17, 21)) +
  scale_y_continuous(
    limits = c(min(-5, min(b_dat$Recovery, na.rm = TRUE)), 150),
    breaks = c(0, 25, 50, 75, 100, 125, 150),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(
    title = "Whole-process spike recovery with and without placental matrix",
    x = NULL,
    y = "Apparent recovery (%)",
    shape = NULL
  ) +
  theme_nc +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

save_panel(
  p5b,
  "Fig5b_whole_process_recovery.png",
  width = 8.5,
  height = 5.0
)

# ============================================================
# PANEL C — PU/PET analytical-confidence categories
# ============================================================

qc <- read_csv(
  file.path(validation_dir, "pu_pet_qc_category_counts_by_delivery_class.csv"),
  show_col_types = FALSE
) %>%
  filter(
    PTB %in% c("Term", "Preterm"),
    Polymer %in% c("PU", "PET")
  ) %>%
  mutate(
    PTB = factor(
      PTB,
      levels = c("Term", "Preterm"),
      labels = c("Term\n(n = 85)", "Preterm\n(n = 86)")
    ),
    Polymer = factor(
      Polymer,
      levels = c("PU", "PET"),
      labels = c(
        "Polyurethane (PU)",
        "Polyethylene terephthalate (PET)"
      )
    ),
    QC = factor(
      QC,
      levels = c(
        "not_detected",
        "below_LOB",
        "above_LOB_below_LOD",
        "above_LOD_below_LOQ",
        "above_LOQ"
      ),
      labels = c(
        "Not detected",
        "< LOB",
        "LOB–LOD",
        "LOD–LOQ",
        "> LOQ"
      )
    )
  )

# Nature-style ordered analytical-confidence gradient.
qc_cols <- c(
  "Not detected" = "#f0f0f0",
  "< LOB" = "#cccccc",
  "LOB–LOD" = "#9ecae1",
  "LOD–LOQ" = "#fdae6b",
  "> LOQ" = "#e6550d"
)

p5c <- ggplot(
  qc,
  aes(x = PTB, y = pct_within_group / 100, fill = QC)
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.25
  ) +
  facet_wrap(~ Polymer, nrow = 1) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.015))
  ) +
  scale_fill_manual(values = qc_cols, drop = FALSE) +
  labs(
    title = "PU and PET analytical confidence by delivery class",
    x = NULL,
    y = "Percent of placentae",
    fill = NULL
  ) +
  theme_nc +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.65, "cm")
  )

save_panel(
  p5c,
  "Fig5c_PU_PET_QC_categories_by_delivery_class.png",
  width = 7.2,
  height = 4.1
)

# ============================================================
# PANEL D — PU/PET targeted analytical robustness
# ============================================================

bin <- read_csv(
  file.path(validation_dir, "pu_pet_qc_binary_ptb_sensitivity_models.csv"),
  show_col_types = FALSE
)

firth <- read_csv(
  file.path(validation_dir, "pu_pet_qc_bias_reduced_models.csv"),
  show_col_types = FALSE
)

blank <- read_csv(
  file.path(validation_dir, "pu_pet_batchmatched_maxblank_subtracted_continuous_models.csv"),
  show_col_types = FALSE
)

d_cont <- blank %>%
  filter(
    Polymer %in% c("PU", "PET"),
    Sensitivity %in% c(
      "Original_recomputed_from_analysis_ready_table",
      "MaxBlank_subtracted_per_run"
    )
  ) %>%
  transmute(
    Polymer,
    Section = "Continuous burden models",
    Model = case_when(
      Sensitivity == "Original_recomputed_from_analysis_ready_table" ~
        "Primary continuous\n(per doubling)",
      Sensitivity == "MaxBlank_subtracted_per_run" ~
        "Max-blank subtracted\n(per doubling)"
    ),
    OR = OR_per_doubling_plus1,
    CI_lo,
    CI_hi,
    Family = "Continuous"
  )

d_bin <- bin %>%
  filter(
    Polymer %in% c("PU", "PET"),
    Threshold %in% c("at_or_above_LOD", "above_LOQ")
  ) %>%
  transmute(
    Polymer,
    Section = "Analytical-threshold models",
    Model = case_when(
      Threshold == "at_or_above_LOD" ~ "≥LOD vs <LOD",
      Threshold == "above_LOQ" ~ ">LOQ vs ≤LOQ"
    ),
    OR,
    CI_lo,
    CI_hi,
    Family = "Threshold"
  )

d_firth <- firth %>%
  filter(
    Polymer %in% c("PU", "PET"),
    Threshold == "above_LOQ",
    is.finite(OR)
  ) %>%
  transmute(
    Polymer,
    Section = "Analytical-threshold models",
    Model = ">LOQ vs ≤LOQ\n(Firth)",
    OR,
    CI_lo,
    CI_hi,
    Family = "Bias-reduced"
  )

model_levels <- rev(c(
  "Primary continuous\n(per doubling)",
  "Max-blank subtracted\n(per doubling)",
  "≥LOD vs <LOD",
  ">LOQ vs ≤LOQ",
  ">LOQ vs ≤LOQ\n(Firth)"
))

d_plot <- bind_rows(d_cont, d_bin, d_firth) %>%
  mutate(
    Polymer = factor(
      Polymer,
      levels = c("PU", "PET"),
      labels = c(
        "Polyurethane (PU)",
        "Polyethylene terephthalate (PET)"
      )
    ),
    Model = factor(Model, levels = model_levels),
    Family = factor(
      Family,
      levels = c("Continuous", "Threshold", "Bias-reduced")
    )
  )

robust_cols <- c(
  "Continuous" = "black",
  "Threshold" = "#2b8cbe",
  "Bias-reduced" = "#d7301f"
)

# Section-label dataframe for annotation in each facet.
section_labels <- tibble::tribble(
  ~Polymer, ~x, ~y, ~label,
  "Polyurethane (PU)", 0.74, 4.55, "Continuous burden models",
  "Polyurethane (PU)", 0.74, 2.55, "Analytical-threshold models",
  "Polyethylene terephthalate (PET)", 0.74, 4.55, "Continuous burden models",
  "Polyethylene terephthalate (PET)", 0.74, 2.55, "Analytical-threshold models"
) %>%
  mutate(
    Polymer = factor(
      Polymer,
      levels = c("Polyurethane (PU)", "Polyethylene terephthalate (PET)")
    )
  )

p5d <- ggplot(
  d_plot,
  aes(x = OR, y = Model, color = Family)
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.5,
    color = "grey40"
  ) +
  # Divider between continuous and threshold models
  geom_hline(
    yintercept = 3.5,
    linewidth = 0.45,
    color = "grey75"
  ) +
  geom_segment(
    aes(x = CI_lo, xend = CI_hi, yend = Model),
    linewidth = 0.8
  ) +
  geom_point(size = 2.7) +
  facet_wrap(~ Polymer, nrow = 1) +
  scale_x_log10(
    limits = c(0.7, 35),
    breaks = c(1, 2, 5, 10, 25),
    labels = c("1", "2", "5", "10", "25")
  ) +
  scale_color_manual(values = robust_cols) +
  labs(
    title = "PU and PET analytical robustness",
    x = "Odds ratio for preterm birth (95% CI)",
    y = NULL,
    color = NULL
  ) +
  theme_nc +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35)
  )

save_panel(
  p5d,
  "Fig5d_PU_PET_analytical_robustness.png",
  width = 8.2,
  height = 5.1
)

# ============================================================
# Optional PDF copies for vector-quality assembly
# ============================================================

ggsave(
  file.path(fig_dir, "Fig5a_validation_workflow.pdf"),
  p5a, width = 10.2, height = 4.0, units = "in", bg = "white"
)
ggsave(
  file.path(fig_dir, "Fig5b_whole_process_recovery.pdf"),
  p5b, width = 8.5, height = 5.0, units = "in", bg = "white"
)
ggsave(
  file.path(fig_dir, "Fig5c_PU_PET_QC_categories_by_delivery_class.pdf"),
  p5c, width = 7.2, height = 4.1, units = "in", bg = "white"
)
ggsave(
  file.path(fig_dir, "Fig5d_PU_PET_analytical_robustness.pdf"),
  p5d, width = 8.2, height = 5.1, units = "in", bg = "white"
)

# ------------------------------------------------------------
# Reproducibility note
# ------------------------------------------------------------
note <- c(
  paste0("Validation output directory: ", validation_dir),
  paste0(
    "Spike-only control block used in Fig5b: ",
    ifelse(
      is.na(SPIKE_ONLY_CONTROL_BLOCK),
      "pooled A+B (working figure; awaiting matrix mapping)",
      SPIKE_ONLY_CONTROL_BLOCK
    )
  ),
  "Fig5d uses complete-case adjusted PTB model outputs from the validation script.",
  "Threshold and Firth ORs are not directly comparable in scale to per-doubling continuous ORs; labels preserve the exposure definition."
)

writeLines(
  note,
  file.path(fig_dir, "Figure5_generation_notes.txt")
)

message("")
message("Figure 5 panels complete.")
message("Output directory: ", fig_dir)
message("Files:")
message("  Fig5a_validation_workflow.png")
message("  Fig5b_whole_process_recovery.png")
message("  Fig5c_PU_PET_QC_categories_by_delivery_class.png")
message("  Fig5d_PU_PET_analytical_robustness.png")
