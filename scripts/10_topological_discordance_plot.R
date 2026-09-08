# =============================================================================
# PUBLICATION-QUALITY TOPOLOGICAL DISCORDANCE BOXPLOT
# =============================================================================
#
# INPUT:
#   matched_bootstrap_replicates.csv
#
# REQUIRED COLUMNS:
#   family
#   comparison
#   grf
#
# DESIGN:
#   Three independent panels:
#
#       Upstream-CDS
#       CDS-Downstream
#       Upstream-Downstream
#
#   Each panel contains five independent boxplots:
#
#       IMP | NDM | KPC | TEM | SHV
#
#   Box:
#       Q1 to Q3
#
#   Black line inside box:
#       median
#
#   Standard whiskers:
#       1.5 × IQR
#
#   Black points:
#       individual bootstrap replicates
#
#   Black capped error bar:
#       2.5th–97.5th bootstrap percentile interval
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. PACKAGES
# -----------------------------------------------------------------------------

required_packages <- c("ggplot2", "dplyr", "readr")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing R packages: ", paste(missing, collapse = ", "),
    ". Install them from environment-r.yml before running this script."
  )
}

library(ggplot2)
library(dplyr)
library(readr)


# -----------------------------------------------------------------------------
# 2. INPUT FILE
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i >= length(args)) stop(sprintf("Missing value for %s", flag))
  args[[i + 1L]]
}

input_file <- arg_value("--input", "matched_bootstrap_replicates.csv")


# -----------------------------------------------------------------------------
# 3. OUTPUT DIRECTORY
# -----------------------------------------------------------------------------

output_dir <- arg_value("--output-dir", "results/topological_discordance_plot")


# -----------------------------------------------------------------------------
# 4. CREATE OUTPUT DIRECTORY
# -----------------------------------------------------------------------------

if (!dir.exists(output_dir)) {

  dir.create(
    output_dir,
    recursive = TRUE
  )

}


# -----------------------------------------------------------------------------
# 5. READ DATA
# -----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("READING BOOTSTRAP DATA\n")
cat("============================================================\n")

dat <- read_csv(
  input_file,
  show_col_types = FALSE
)


# -----------------------------------------------------------------------------
# 6. CHECK REQUIRED COLUMNS
# -----------------------------------------------------------------------------

required_columns <- c(
  "family",
  "comparison",
  "grf"
)

missing_columns <- setdiff(
  required_columns,
  colnames(dat)
)

if (length(missing_columns) > 0) {

  stop(
    paste0(
      "\nERROR: Required columns are missing:\n",
      paste(
        missing_columns,
        collapse = ", "
      ),
      "\n\nColumns found in your CSV:\n",
      paste(
        colnames(dat),
        collapse = ", "
      ),
      "\n"
    )
  )

}


# -----------------------------------------------------------------------------
# 7. CLEAN DATA
# -----------------------------------------------------------------------------

dat <- dat %>%

  filter(
    is.finite(grf)
  ) %>%

  mutate(

    # -------------------------------------------------------------------------
    # FAMILY ORDER
    # -------------------------------------------------------------------------

    family = factor(
      family,
      levels = c(
        "IMP",
        "NDM",
        "KPC",
        "TEM",
        "SHV"
      )
    ),

    # -------------------------------------------------------------------------
    # COMPARISON ORDER
    # -------------------------------------------------------------------------

    comparison = factor(
      comparison,
      levels = c(
        "upstream_vs_cds",
        "cds_vs_downstream",
        "upstream_vs_downstream"
      ),
      labels = c(
        "Upstream–CDS",
        "CDS–Downstream",
        "Upstream–Downstream"
      )
    )

  ) %>%

  filter(
    !is.na(family),
    !is.na(comparison)
  )


# -----------------------------------------------------------------------------
# 8. DATA SUMMARY
# -----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("DATA SUMMARY\n")
cat("============================================================\n")

cat(
  "Total bootstrap observations:",
  nrow(dat),
  "\n\n"
)

cat("Observations per family and comparison:\n\n")

print(
  dat %>%
    count(
      comparison,
      family
    )
)


# -----------------------------------------------------------------------------
# 9. CALCULATE SUMMARY STATISTICS
# -----------------------------------------------------------------------------
#
# median:
#   50th percentile
#
# lower_95:
#   2.5th percentile
#
# upper_95:
#   97.5th percentile
#
# These are calculated independently for every:
#
#   comparison × family
#
# combination.
#

summary_dat <- dat %>%

  group_by(
    comparison,
    family
  ) %>%

  summarise(

    median_grf =
      median(
        grf,
        na.rm = TRUE
      ),

    lower_95 =
      quantile(
        grf,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE
      ),

    upper_95 =
      quantile(
        grf,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE
      ),

    n =
      sum(
        is.finite(grf)
      ),

    .groups = "drop"

  )


# -----------------------------------------------------------------------------
# 10. FAMILY COLOURS
# -----------------------------------------------------------------------------
#
# Five clearly distinct colours.
#
# IMP = gold
# NDM = light blue
# KPC = green
# TEM = orange
# SHV = dark blue
#

family_colors <- c(

  "IMP" = "#E69F00",

  "NDM" = "#56B4E9",

  "KPC" = "#009E73",

  "TEM" = "#D55E00",

  "SHV" = "#0072B2"

)


# -----------------------------------------------------------------------------
# 11. CREATE THE PLOT
# -----------------------------------------------------------------------------

p <- ggplot(

  dat,

  aes(
    x = family,
    y = grf
  )

) +


  # ===========================================================================
  # A. 95% BOOTSTRAP INTERVAL
  # ===========================================================================
  #
  # This is drawn FIRST so that the boxplot remains visually dominant.
  #
  # Black vertical line:
  #
  #       97.5%
  #          ─
  #          │
  #          │
  #          │
  #          ─
  #        2.5%
  #
  # It is NOT the boxplot whisker.
  #
  # It represents the bootstrap percentile interval.
  #

  geom_errorbar(

    data = summary_dat,

    aes(
      x = family,
      ymin = lower_95,
      ymax = upper_95
    ),

    width = 0.18,

    linewidth = 0.90,

    colour = "black",

    inherit.aes = FALSE

  ) +


  # ===========================================================================
  # B. INDIVIDUAL BOOTSTRAP POINTS
  # ===========================================================================
  #
  # ALL POINTS ARE BLACK.
  #
  # The points are distributed horizontally within the box width.
  #
  # x has NO biological/quantitative meaning, so the horizontal jitter is
  # purely for visualization.
  #
  # jitter.width = 0.22 keeps points inside/near the box width.
  #

  geom_jitter(

    width = 0.22,

    height = 0,

    size = 0.85,

    alpha = 0.30,

    colour = "black",

    stroke = 0,

    show.legend = FALSE

  ) +


  # ===========================================================================
  # C. BOXPLOT
  # ===========================================================================
  #
  # Each family gets its own independent box.
  #
  # Box:
  #       Q1 -> Q3
  #
  # Black horizontal line:
  #       median
  #
  # Whiskers:
  #       1.5 × IQR
  #
  # Whisker caps:
  #       horizontal black lines
  #
  # Standard outlier symbols are suppressed because ALL individual bootstrap
  # observations are already shown as black points.
  #

  geom_boxplot(

    aes(
      fill = family
    ),

    width = 0.62,

    linewidth = 0.75,

    colour = "black",

    alpha = 0.82,

    outlier.shape = NA,

    show.legend = TRUE

  ) +


  # ===========================================================================
  # D. FAMILY COLOUR SCALE
  # ===========================================================================

  scale_fill_manual(

    values = family_colors,

    breaks = c(
      "IMP",
      "NDM",
      "KPC",
      "TEM",
      "SHV"
    ),

    drop = FALSE,

    name = NULL

  ) +


  # ===========================================================================
  # E. FACET INTO THREE INDIVIDUAL PANELS
  # ===========================================================================
  #
  # This is the key change.
  #
  # Instead of putting all three comparisons into one x-axis,
  # each comparison becomes its OWN PANEL.
  #
  # Each panel independently contains:
  #
  #       IMP  NDM  KPC  TEM  SHV
  #
  # The Y axis remains shared across all three panels.
  #

  facet_wrap(

    ~ comparison,

    nrow = 1,

    scales = "fixed"

  ) +


  # ===========================================================================
  # F. Y AXIS
  # ===========================================================================

  scale_y_continuous(

    breaks = seq(
      0.65,
      1.00,
      by = 0.05
    ),

    expand = expansion(

      mult = c(
        0.015,
        0.035
      )

    )

  ) +


  # ===========================================================================
  # G. AXIS LABELS
  # ===========================================================================

  labs(

    x = NULL,

    y = "Normalized Clustering Information Distance"

  ) +


  # ===========================================================================
  # H. PUBLICATION THEME
  # ===========================================================================

  theme_classic(

    base_size = 15

  ) +

  theme(

    # -------------------------------------------------------------------------
    # NO GRID
    # -------------------------------------------------------------------------

    panel.grid.major = element_blank(),

    panel.grid.minor = element_blank(),


    # -------------------------------------------------------------------------
    # PANEL STRIPS
    # -------------------------------------------------------------------------
    #
    # These are the three comparison labels.
    # They are NOT a figure title.
    #

    strip.background = element_blank(),

    strip.text = element_text(

      size = 15,

      face = "bold",

      colour = "black",

      margin = margin(
        b = 8
      )

    ),


    # -------------------------------------------------------------------------
    # PANEL SPACING
    # -------------------------------------------------------------------------

    panel.spacing = unit(

      1.25,

      "cm"

    ),


    # -------------------------------------------------------------------------
    # AXIS LINES
    # -------------------------------------------------------------------------

    axis.line = element_line(

      linewidth = 0.8,

      colour = "black"

    ),


    # -------------------------------------------------------------------------
    # AXIS TICKS
    # -------------------------------------------------------------------------

    axis.ticks = element_line(

      linewidth = 0.7,

      colour = "black"

    ),

    axis.ticks.length = unit(

      0.18,

      "cm"

    ),


    # -------------------------------------------------------------------------
    # X-AXIS FAMILY LABELS
    # -------------------------------------------------------------------------

    axis.text.x = element_text(

      size = 13,

      colour = "black",

      angle = 0,

      margin = margin(
        t = 8
      )

    ),


    # -------------------------------------------------------------------------
    # Y-AXIS LABELS
    # -------------------------------------------------------------------------

    axis.text.y = element_text(

      size = 13,

      colour = "black"

    ),


    # -------------------------------------------------------------------------
    # Y-AXIS TITLE
    # -------------------------------------------------------------------------

    axis.title.y = element_text(

      size = 15,

      face = "plain",

      colour = "black",

      margin = margin(
        r = 12
      )

    ),


    # -------------------------------------------------------------------------
    # LEGEND
    # -------------------------------------------------------------------------

    legend.position = "top",

    legend.direction = "horizontal",

    legend.justification = "center",

    legend.text = element_text(

      size = 13,

      colour = "black"

    ),

    legend.key.width = unit(

      0.75,

      "cm"

    ),

    legend.key.height = unit(

      0.45,

      "cm"

    ),

    legend.spacing.x = unit(

      0.25,

      "cm"

    ),

    legend.margin = margin(

      b = 8

    ),


    # -------------------------------------------------------------------------
    # NO FIGURE TITLE
    # -------------------------------------------------------------------------

    plot.title = element_blank(),

    plot.subtitle = element_blank(),

    plot.caption = element_blank(),


    # -------------------------------------------------------------------------
    # PLOT MARGINS
    # -------------------------------------------------------------------------

    plot.margin = margin(

      t = 8,

      r = 15,

      b = 15,

      l = 15

    )

  )


# -----------------------------------------------------------------------------
# 12. DISPLAY
# -----------------------------------------------------------------------------

print(p)


# -----------------------------------------------------------------------------
# 13. SAVE HIGH-RESOLUTION PNG
# -----------------------------------------------------------------------------

ggsave(

  filename = file.path(
    output_dir,
    "topological_discordance_boxplot.png"
  ),

  plot = p,

  width = 13,

  height = 7.5,

  units = "in",

  dpi = 600,

  bg = "white"

)


# -----------------------------------------------------------------------------
# 14. SAVE VECTOR PDF
# -----------------------------------------------------------------------------

ggsave(

  filename = file.path(
    output_dir,
    "topological_discordance_boxplot.pdf"
  ),

  plot = p,

  width = 13,

  height = 7.5,

  units = "in",

  device = cairo_pdf,

  bg = "white"

)


# -----------------------------------------------------------------------------
# 15. SAVE VECTOR SVG
# -----------------------------------------------------------------------------

ggsave(

  filename = file.path(
    output_dir,
    "topological_discordance_boxplot.svg"
  ),

  plot = p,

  width = 13,

  height = 7.5,

  units = "in",

  bg = "white"

)


# -----------------------------------------------------------------------------
# 16. SAVE SUMMARY TABLE
# -----------------------------------------------------------------------------

write.csv(

  summary_dat,

  file.path(
    output_dir,
    "bootstrap_summary_for_plot.csv"
  ),

  row.names = FALSE

)


# -----------------------------------------------------------------------------
# 17. FINISH
# -----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("PLOT GENERATION COMPLETE\n")
cat("============================================================\n")

cat("\nOutput folder:\n")

cat(
  normalizePath(
    output_dir
  )
)

cat("\n\nFiles created:\n")

cat(
  "  topological_discordance_boxplot.png\n"
)

cat(
  "  topological_discordance_boxplot.pdf\n"
)

cat(
  "  topological_discordance_boxplot.svg\n"
)

cat(
  "  bootstrap_summary_for_plot.csv\n"
)

cat("\n")