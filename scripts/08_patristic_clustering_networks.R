#!/usr/bin/env Rscript
## =============================================================================
## bla_phylo_cluster_workflow.R
##
## Phylogenetic cluster architecture and cross-region cluster correspondence
## for five Klebsiella pneumoniae beta-lactamase families (IMP, NDM, KPC, TEM, SHV)
##
## For each family, three IQ-TREE trees are analysed independently:
##   upstream (1 kb) | CDS | downstream (1 kb)
##
## METHOD (UNCHANGED from the previously validated working version -- this
## revision only rewrites Figure 1 and Figure 2 rendering code):
##   1. Read Newick trees, intersect tip labels shared across all 3 regions.
##   2. Patristic distances: D <- ape::cophenetic.phylo(tree)   [no normalization]
##   3. Clustering on the FULL patristic distance matrix:
##        n <= 300               -> PAM (cluster::pam, diss = TRUE)
##        n  > 300               -> hierarchical, average linkage
##                                   (hclust(as.dist(D), method = "average"))
##      TEM downstream (n ~ 565) and SHV (n ~ 1971) therefore ALWAYS use
##      hierarchical average-linkage clustering and NEVER PAM.
##   4. k chosen by maximum average silhouette width, k in 2..min(15, floor(n/15))
##        n <= 1000  -> silhouette evaluated on the FULL distance matrix
##        n  > 1000  -> silhouette evaluated on a fixed random subset of 1000
##                       sequences (set.seed(42)); final cluster labels for ALL
##                       sequences still come from the full hclust/cutree.
##   5. PCoA (ape::pcoa, Cailliez correction) of the SAME patristic distances,
##      used ONLY for 2D visualization. Clustering is never performed on PCoA
##      coordinates. % variance explained on axes 1-2 is also extracted here,
##      for axis labelling only.
##   6. Clusters renumbered deterministically per family x region, C1 = largest.
##   7. Cross-region correspondence (same taxa): ARI, NMI, normalized VI,
##      computed with base-R implementations (no aricode dependency).
##   8. Upstream -> CDS-variant -> Downstream transition counts (full,
##      sequence-level resolution) -- the exact data underlying Figure 2.
##      A separate, VISUALIZATION-ONLY copy of this table pools rare CDS
##      variants into "Other variants" for Figure 2 legibility; the full
##      resolution table is untouched and is what gets written to CSV.
##
## IMPORTANT INTERPRETATION NOTE (see caption text below and in console output):
##   This script produces "phylogenetic cluster architecture" and "cluster-
##   membership correspondence" results only. It does NOT test for, and must
##   NOT be described as evidence of, recombination, linkage disequilibrium,
##   ancestry, or horizontal gene transfer.
##
## USAGE:
##   Rscript bla_phylo_cluster_workflow.R
##
## REQUIRED PACKAGES (install once):
##   Rscript -e 'install.packages(c("ape","cluster","ggplot2","dplyr","tidyr","patchwork","ggrepel"), repos="https://cloud.r-project.org")'
## SVG export uses svglite if available; otherwise PDF/PNG are produced and SVG
## is silently skipped (no crash). `grid` is a base R package shipped with
## every R installation and requires no separate install.
## =============================================================================

suppressWarnings(suppressPackageStartupMessages({
  library(ape)
  library(cluster)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(grid)       # base R package (ships with R) -- used only for unit()
}))

HAVE_SVGLITE <- requireNamespace("svglite", quietly = TRUE)
HAVE_GGREPEL <- requireNamespace("ggrepel", quietly = TRUE)
if (!HAVE_GGREPEL) {
  stop("Package ggrepel is required for the publication-quality Figure 1 labels. Install it with: install.packages('ggrepel') or conda install -c conda-forge r-ggrepel")
}

## -----------------------------------------------------------------------
## 0. CONFIGURATION
## -----------------------------------------------------------------------
FAMILIES  <- c("imp", "kpc", "ndm", "tem", "shv")
FAMILY_LABELS <- c(imp = "IMP", ndm = "NDM", kpc = "KPC", tem = "TEM", shv = "SHV")
REGIONS   <- c("upstream", "cds", "downstream")
REGION_LABELS <- c(upstream = "Upstream", cds = "CDS", downstream = "Downstream")

## Scalar-lookup helpers that always return an UNNAMED length-1 string.
## (Indexing a named vector like FAMILY_LABELS["tem"] returns a named result,
##  which triggers spurious "row names were found from a short variable"
##  warnings when mixed into data.frame() with longer columns.)
fam_label <- function(fam) unname(FAMILY_LABELS[fam])
region_label <- function(region) unname(REGION_LABELS[region])

TREE_DIR      <- "trees"
METADATA_CSV  <- "metadata.csv"
OUT_ROOT      <- "bla_cluster_results"
DIR_CLUSTERING <- file.path(OUT_ROOT, "clustering")
DIR_COMPARISON <- file.path(OUT_ROOT, "comparison")
DIR_PCOA       <- file.path(OUT_ROOT, "pcoa")
DIR_TREECHECK  <- file.path(OUT_ROOT, "trees_checked")
DIR_FIGURES    <- file.path(OUT_ROOT, "figures")

for (d in c(OUT_ROOT, DIR_CLUSTERING, DIR_COMPARISON, DIR_PCOA, DIR_TREECHECK, DIR_FIGURES)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

PAM_MAX_N        <- 300     # n <= this -> PAM; n > this -> hierarchical average linkage
SIL_SUBSET_MAX_N <- 1000    # n > this -> evaluate silhouette on a fixed random subset
SIL_SUBSET_SIZE  <- 1000
KMAX_CAP         <- 15
KMAX_DIVISOR     <- 15
GLOBAL_SEED      <- 42

## Figure-2 display-only parameter: CDS variants occurring fewer than this many
## times in a family (summed across all upstream/downstream combinations) are
## shown pooled as "Other variants" in Figure 2 ONLY. This never changes the
## underlying sequence-level data or any CSV output.
DISPLAY_VARIANT_MIN_N <- 5

set.seed(GLOBAL_SEED)

CLUSTER_PALETTE <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#B15928",
  "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462", "#BC80BD"
)
get_palette <- function(n) rep(CLUSTER_PALETTE, length.out = n)

## Separate qualitative palette for CDS variants in Figure 2 (kept distinct
## from CLUSTER_PALETTE only so cluster colors in Fig 1 and variant colors in
## Fig 2 are never visually confused with one another).
VARIANT_PALETTE <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
  "#1F78B4", "#B15928", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#A6CEE3", "#33A02C", "#FB9A99", "#FDBF6F"
)
OTHER_GREY <- "#B0B0B0"   # reserved for "Other variants" / "Unknown" only

## Figure-2 node/label geometry (visualization only -- does not affect any
## flow weight, node height, or stacking order, all of which are computed
## from n_sequences exactly as before).
NODE_HALF_WIDTH       <- 0.09   # half-width of each stratum block (was 0.035)
INSIDE_LABEL_MIN_FRAC <- 0.10   # a node must span >= this fraction of its
                                 # column's total height to hold its label
                                 # INSIDE the block; smaller nodes get an
                                 # outside label with a short leader line,
                                 # or are omitted if no safe position exists

## -----------------------------------------------------------------------
## 1. METADATA (optional)
## -----------------------------------------------------------------------
metadata <- NULL
if (file.exists(METADATA_CSV)) {
  metadata <- tryCatch(read.csv(METADATA_CSV, stringsAsFactors = FALSE),
                        error = function(e) NULL)
  if (!is.null(metadata) && !all(c("tip_label", "family", "variant") %in% names(metadata))) {
    warning("metadata.csv found but missing required columns (tip_label, family, variant). Ignoring it.")
    metadata <- NULL
  }
}
if (is.null(metadata)) {
  message("WARNING: metadata.csv not found or invalid.")
  message("Variant identities will be extracted from CDS tip labels instead.")
}

## Extract a beta-lactamase variant string from a tip label. Handles forms such
## as blaKPC-2, KPC-2, blaKPC_2, KPC_2, bla_KPC_2, kpc-2, etc. Falls back to
## "Unknown" rather than crashing.
extract_variant <- function(tip_labels, family_prefix) {
  fp <- toupper(family_prefix)
  pattern <- paste0("(?i)(BLA[_-]?)?", fp, "[_-]?([0-9]+[A-Za-z]*)")
  m <- regmatches(tip_labels, regexpr(pattern, tip_labels, perl = TRUE))
  out <- rep("Unknown", length(tip_labels))
  has_match <- nchar(m) > 0
  if (any(has_match)) {
    num <- regmatches(m, regexpr("[0-9]+[A-Za-z]*$", m, perl = TRUE))
    out[has_match] <- paste0(fp, "-", num)
  }
  out
}

get_variants_for_family <- function(fam, tip_labels, cds_tip_labels) {
  if (!is.null(metadata)) {
    sub <- metadata[metadata$family == fam | toupper(metadata$family) == toupper(fam), ]
    if (nrow(sub) > 0) {
      v <- setNames(rep("Unknown", length(tip_labels)), tip_labels)
      matched <- intersect(tip_labels, sub$tip_label)
      if (length(matched) > 0) {
        v[matched] <- sub$variant[match(matched, sub$tip_label)]
      }
      missing <- setdiff(tip_labels, sub$tip_label)
      if (length(missing) > 0) {
        v[missing] <- extract_variant(missing, fam)
      }
      return(v)
    }
  }
  setNames(extract_variant(tip_labels, fam), tip_labels)
}

## -----------------------------------------------------------------------
## 2. TREE READING / CHECKING
## -----------------------------------------------------------------------
tree_check_rows <- list()

read_and_check_tree <- function(path, fam, region) {
  status <- "OK"
  n_tips <- NA_integer_
  tr <- NULL
  if (!file.exists(path)) {
    status <- "MISSING_FILE"
  } else {
    tr <- tryCatch(ape::read.tree(path), error = function(e) NULL)
    if (is.null(tr)) {
      status <- "UNREADABLE"
    } else if (is.null(tr$tip.label) || length(tr$tip.label) < 3) {
      status <- "TOO_FEW_TIPS"
      n_tips <- length(tr$tip.label)
    } else {
      n_tips <- length(tr$tip.label)
    }
  }
  tree_check_rows[[length(tree_check_rows) + 1]] <<- data.frame(
    family = fam, region = region, path = path, status = status, n_tips = n_tips
  )
  if (status != "OK") return(NULL)
  tr
}

## -----------------------------------------------------------------------
## 3. PATRISTIC DISTANCE + VALIDATION
## -----------------------------------------------------------------------
compute_patristic <- function(tree, keep_tips) {
  drop_these <- setdiff(tree$tip.label, keep_tips)
  if (length(drop_these) > 0) {
    tree <- ape::drop.tip(tree, drop_these)
  }
  if (is.null(tree$tip.label) || length(tree$tip.label) < 3) return(NULL)
  D <- ape::cophenetic.phylo(tree)
  D <- D[keep_tips[keep_tips %in% rownames(D)], keep_tips[keep_tips %in% colnames(D)], drop = FALSE]
  ## validation
  if (any(is.na(D)) || any(is.nan(D)) || any(is.infinite(D))) {
    warning("Non-finite values in patristic distance matrix; attempting to clean.")
    D[!is.finite(D)] <- NA
    if (any(is.na(D))) return(NULL)
  }
  diag(D) <- 0
  D <- (D + t(D)) / 2   # enforce numerical symmetry
  D
}

## -----------------------------------------------------------------------
## 4. PCoA (visualization only) -- also returns % variance explained on
##    axes 1-2, used only for axis labelling in Figure 1.
## -----------------------------------------------------------------------
run_pcoa <- function(D) {
  pco <- tryCatch(ape::pcoa(as.dist(D), correction = "cailliez"),
                   error = function(e) NULL)
  if (is.null(pco)) {
    pco <- ape::pcoa(as.dist(D))
  }
  vec <- pco$vectors.cor
  if (is.null(vec) || ncol(vec) < 2) vec <- pco$vectors
  if (ncol(vec) < 2) {
    ## pathological case: fewer than 2 usable axes -- pad with zeros
    vec <- cbind(vec, matrix(0, nrow(vec), 2 - ncol(vec)))
  }
  coords <- vec[, 1:2, drop = FALSE]
  colnames(coords) <- c("PCo1", "PCo2")
  rownames(coords) <- rownames(D)

  varexp <- pco$values$Rel_corr_eig
  if (is.null(varexp)) varexp <- pco$values$Relative_eig
  varexp_pct <- if (!is.null(varexp) && length(varexp) >= 2) {
    round(100 * as.numeric(varexp[1:2]), 1)
  } else {
    c(NA_real_, NA_real_)
  }
  list(coords = coords, varexp_pct = varexp_pct)
}

## -----------------------------------------------------------------------
## 5. CLUSTERING: PAM (n<=300) or hierarchical average linkage (n>300)
##    k chosen by maximum average silhouette width
##    *** UNCHANGED FROM THE WORKING VERSION ***
## -----------------------------------------------------------------------
compute_silhouette_avg <- function(labels_sub, d_sub) {
  if (length(unique(labels_sub)) < 2) return(NA_real_)
  sil <- tryCatch(cluster::silhouette(as.integer(as.factor(labels_sub)), d_sub),
                   error = function(e) NULL)
  if (is.null(sil)) return(NA_real_)
  mean(sil[, "sil_width"])
}

cluster_region <- function(D, fam, region) {
  n <- nrow(D)
  taxa <- rownames(D)
  dmat_full <- as.dist(D)
  kmax <- max(2, min(KMAX_CAP, floor(n / KMAX_DIVISOR)))
  if (kmax < 2) kmax <- 2

  method <- if (n <= PAM_MAX_N) "PAM" else "hierarchical average linkage"
  message(sprintf("Processing: %s / %s", fam_label(fam), region_label(region)))
  message(sprintf("  n = %d", n))
  message(sprintf("  method = %s", method))

  use_subset_for_sil <- n > SIL_SUBSET_MAX_N
  sil_idx <- seq_len(n)
  if (use_subset_for_sil) {
    set.seed(GLOBAL_SEED)
    sil_idx <- sort(sample(seq_len(n), size = min(SIL_SUBSET_SIZE, n)))
    message(sprintf("  silhouette evaluation = %d-sequence random subset (final labels use full data)",
                     length(sil_idx)))
  } else {
    message("  silhouette evaluation = full dataset")
  }
  sub_taxa <- taxa[sil_idx]
  D_sub <- D[sil_idx, sil_idx, drop = FALSE]
  dmat_sub <- as.dist(D_sub)

  hc_full <- NULL
  if (method == "hierarchical average linkage") {
    hc_full <- stats::hclust(dmat_full, method = "average")
  }

  sil_curve <- data.frame(k = integer(0), silhouette = numeric(0))
  fits_full <- list()   # k (character) -> full-data label vector (named)

  for (k in 2:kmax) {
    message(sprintf("  testing k = %d [%s]", k, method))

    if (method == "PAM") {
      fit_full <- cluster::pam(dmat_full, k = k, diss = TRUE)
      labels_full <- as.integer(fit_full$clustering)
      names(labels_full) <- taxa
      if (use_subset_for_sil) {
        labels_sub <- labels_full[sub_taxa]
        sil_val <- compute_silhouette_avg(labels_sub, dmat_sub)
      } else {
        sil_val <- compute_silhouette_avg(labels_full, dmat_full)
      }
    } else {
      labels_full <- as.integer(stats::cutree(hc_full, k = k))
      names(labels_full) <- taxa
      labels_sub <- labels_full[sub_taxa]
      sil_val <- compute_silhouette_avg(labels_sub, dmat_sub)
    }

    ## guard against degenerate/empty/singleton-only partitions
    tab <- table(labels_full)
    if (length(tab) < 2) sil_val <- NA_real_

    sil_curve <- rbind(sil_curve, data.frame(k = k, silhouette = sil_val))
    fits_full[[as.character(k)]] <- labels_full
  }

  valid <- !is.na(sil_curve$silhouette)
  if (!any(valid)) {
    ## fallback: k = 2 with whatever partition was produced
    best_k <- 2
    warning(sprintf("No valid silhouette values for %s/%s; defaulting to k=2.", fam, region))
  } else {
    best_k <- sil_curve$k[valid][which.max(sil_curve$silhouette[valid])]
  }
  best_labels <- fits_full[[as.character(best_k)]]
  best_sil <- sil_curve$silhouette[sil_curve$k == best_k][1]

  ## sanity checks
  if (length(best_labels) != n) stop(sprintf("Label length mismatch for %s/%s", fam, region))
  if (is.null(names(best_labels)) || length(names(best_labels)) != length(best_labels)) {
    stop(sprintf("Missing/incorrect names on cluster labels for %s/%s", fam, region))
  }

  message(sprintf("  selected k = %d", best_k))
  message(sprintf("  silhouette = %.4f", ifelse(is.na(best_sil), NA, best_sil)))

  ## deterministic renumbering: C1 = largest cluster, by decreasing size
  sizes <- sort(table(best_labels), decreasing = TRUE)
  remap <- setNames(seq_along(sizes), names(sizes))
  final_labels <- paste0("C", remap[as.character(best_labels)])
  names(final_labels) <- names(best_labels)

  list(
    labels = final_labels,
    k = best_k,
    silhouette = best_sil,
    method = method,
    sil_curve = sil_curve,
    n = n
  )
}

## -----------------------------------------------------------------------
## 6. Cluster summary + composition tables
##    *** UNCHANGED FROM THE WORKING VERSION ***
## -----------------------------------------------------------------------
build_cluster_summary <- function(fam, region, labels, variants) {
  taxa <- names(labels)
  df <- data.frame(
    family = fam_label(fam), region = region_label(region),
    tip_label = taxa, cluster = unname(labels[taxa]),
    variant = unname(variants[taxa]), stringsAsFactors = FALSE
  )

  n_total <- nrow(df)
  summ <- df %>%
    group_by(cluster) %>%
    summarise(
      n_sequences = n(),
      n_variants  = n_distinct(variant),
      .groups = "drop"
    ) %>%
    mutate(family = fam_label(fam), region = region_label(region),
           fraction = n_sequences / n_total) %>%
    select(family, region, cluster, n_sequences, fraction, n_variants)

  cluster_totals <- df %>% count(cluster, name = "cluster_n")
  dom <- df %>%
    count(cluster, variant, name = "n") %>%
    left_join(cluster_totals, by = "cluster") %>%
    group_by(cluster) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(cluster, dominant_variant = variant,
              dominant_variant_fraction = n / cluster_n)

  summ <- summ %>% left_join(dom, by = "cluster") %>%
    mutate(.cluster_num = as.integer(sub("^C", "", cluster))) %>%
    arrange(.cluster_num) %>%
    select(-.cluster_num)

  list(assignments = df, summary = summ)
}

## -----------------------------------------------------------------------
## 7. Base-R cross-partition agreement statistics: ARI, NMI, normalized VI
##    *** UNCHANGED FROM THE WORKING VERSION ***
## -----------------------------------------------------------------------
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  sum_comb <- function(x) sum(choose(x, 2))
  a_sums <- rowSums(tab); b_sums <- colSums(tab)
  index <- sum_comb(tab)
  expected <- sum_comb(a_sums) * sum_comb(b_sums) / choose(n, 2)
  max_index <- 0.5 * (sum_comb(a_sums) + sum_comb(b_sums))
  if (max_index == expected) return(1)
  (index - expected) / (max_index - expected)
}

normalized_mutual_information <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  Pij <- tab / n
  Pi <- rowSums(Pij); Pj <- colSums(Pij)
  Hi <- -sum(Pi[Pi > 0] * log(Pi[Pi > 0]))
  Hj <- -sum(Pj[Pj > 0] * log(Pj[Pj > 0]))
  outer_ij <- outer(Pi, Pj)
  nz <- Pij > 0
  MI <- sum(Pij[nz] * log(Pij[nz] / outer_ij[nz]))
  denom <- sqrt(Hi * Hj)
  if (denom == 0) return(0)
  MI / denom
}

normalized_variation_of_information <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  Pij <- tab / n
  Pi <- rowSums(Pij); Pj <- colSums(Pij)
  Hi <- -sum(Pi[Pi > 0] * log(Pi[Pi > 0]))
  Hj <- -sum(Pj[Pj > 0] * log(Pj[Pj > 0]))
  outer_ij <- outer(Pi, Pj)
  nz <- Pij > 0
  MI <- sum(Pij[nz] * log(Pij[nz] / outer_ij[nz]))
  VI <- Hi + Hj - 2 * MI
  if (n <= 1) return(0)
  VI / log(n)
}

## -----------------------------------------------------------------------
## 8. MAIN PER-FAMILY LOOP
##    *** UNCHANGED FROM THE WORKING VERSION. Figure 1 rendering uses
##        the full patristic matrices retained below; clustering logic is untouched. ***
## -----------------------------------------------------------------------
all_family_results <- list()      # fam -> region -> list(labels, k, silhouette, ...)
all_pcoa           <- list()      # fam -> region -> coords data.frame
all_pcoa_varexp    <- list()      # fam -> region -> c(PCo1_pct, PCo2_pct)
all_patristic      <- list()      # fam -> region -> full patristic distance matrix
final_summary_rows <- list()

for (fam in FAMILIES) {

  message("")
  message("============================================")
  message(sprintf("FAMILY: %s", fam_label(fam)))
  message("============================================")

  paths <- setNames(
    file.path(TREE_DIR, paste0(fam, "_", REGIONS, ".treefile")),
    REGIONS
  )
  message("Reading:")
  for (p in paths) message(sprintf("  %s", p))

  trees <- list()
  ok <- TRUE
  for (region in REGIONS) {
    tr <- read_and_check_tree(paths[[region]], fam, region)
    if (is.null(tr)) {
      warning(sprintf("Tree missing/unreadable for %s / %s (%s). Skipping family %s.",
                       fam, region, paths[[region]], fam_label(fam)))
      ok <- FALSE
    } else {
      trees[[region]] <- tr
    }
  }
  if (!ok) next

  shared_taxa <- Reduce(intersect, lapply(trees, function(t) t$tip.label))
  message(sprintf("Shared taxa across all 3 regions: %d", length(shared_taxa)))
  if (length(shared_taxa) < 10) {
    warning(sprintf("Fewer than 10 shared taxa for family %s; skipping.", fam_label(fam)))
    next
  }
  shared_taxa <- sort(shared_taxa)

  ## variant labels (from metadata if available, else parsed from CDS tips)
  variants <- get_variants_for_family(fam, shared_taxa, shared_taxa)

  region_results <- list()
  region_pcoa <- list()
  region_varexp <- list()
  region_cluster_data <- list()
  region_patristic <- list()

  region_ok <- TRUE
  for (region in REGIONS) {
    D <- compute_patristic(trees[[region]], shared_taxa)
    if (is.null(D) || nrow(D) < 10) {
      warning(sprintf("Invalid/too-small distance matrix for %s / %s. Skipping family %s.",
                       fam, region, fam_label(fam)))
      region_ok <- FALSE
      break
    }

    pco_result <- run_pcoa(D)
    coords <- pco_result$coords
    varexp_pct <- pco_result$varexp_pct

    cl <- cluster_region(D, fam, region)

    out <- build_cluster_summary(fam, region, cl$labels, variants)

    region_results[[region]] <- cl
    region_pcoa[[region]] <- data.frame(
      tip_label = rownames(coords), PCo1 = coords[, 1], PCo2 = coords[, 2],
      cluster = unname(cl$labels[rownames(coords)]),
      variant = unname(variants[rownames(coords)]),
      stringsAsFactors = FALSE
    )
    region_varexp[[region]] <- varexp_pct
    region_cluster_data[[region]] <- out
    region_patristic[[region]] <- D

    ## ---- write per-region CSVs (unchanged schemas) ----
    write.csv(out$assignments,
      file.path(DIR_CLUSTERING, sprintf("%s_%s_cluster_assignments.csv", fam, region)),
      row.names = FALSE)
    write.csv(out$summary,
      file.path(DIR_CLUSTERING, sprintf("%s_%s_cluster_summary.csv", fam, region)),
      row.names = FALSE)
    write.csv(cl$sil_curve,
      file.path(DIR_CLUSTERING, sprintf("%s_%s_silhouette.csv", fam, region)),
      row.names = FALSE)
    write.csv(data.frame(
        family = fam_label(fam), region = region_label(region),
        n_sequences = cl$n, method = cl$method,
        selected_k = cl$k, selected_silhouette = cl$silhouette
      ),
      file.path(DIR_CLUSTERING, sprintf("%s_%s_clustering_info.csv", fam, region)),
      row.names = FALSE)
    write.csv(region_pcoa[[region]],
      file.path(DIR_PCOA, sprintf("%s_%s_pcoa.csv", fam, region)),
      row.names = FALSE)

    final_summary_rows[[length(final_summary_rows) + 1]] <- data.frame(
      family = fam_label(fam), region = region_label(region),
      n = cl$n, method = cl$method, k = cl$k, silhouette = cl$silhouette
    )
  }
  if (!region_ok) next

  all_family_results[[fam]] <- region_results
  all_pcoa[[fam]] <- region_pcoa
  all_pcoa_varexp[[fam]] <- region_varexp
  all_patristic[[fam]] <- region_patristic

  ## ---- cross-region concordance (this family) ----
  lab_up <- region_results[["upstream"]]$labels[shared_taxa]
  lab_cds <- region_results[["cds"]]$labels[shared_taxa]
  lab_down <- region_results[["downstream"]]$labels[shared_taxa]

  concordance_pairs <- list(
    list(name = "CDS vs Upstream", a = lab_cds, b = lab_up),
    list(name = "CDS vs Downstream", a = lab_cds, b = lab_down),
    list(name = "Upstream vs Downstream", a = lab_up, b = lab_down)
  )
  fam_conc_rows <- list()
  for (cp in concordance_pairs) {
    fam_conc_rows[[length(fam_conc_rows) + 1]] <- data.frame(
      family = fam_label(fam),
      comparison = cp$name,
      n_taxa = length(shared_taxa),
      ARI = adjusted_rand_index(cp$a, cp$b),
      NMI = normalized_mutual_information(cp$a, cp$b),
      VI_normalized = normalized_variation_of_information(cp$a, cp$b)
    )
  }
  if (fam == FAMILIES[1]) {
    concordance_all <- bind_rows(fam_conc_rows)
  } else {
    concordance_all <- bind_rows(concordance_all, bind_rows(fam_conc_rows))
  }

  ## ---- transition counts: upstream_cluster -> cds_variant -> downstream_cluster ----
  ## (full sequence-level resolution -- this is the authoritative data written to CSV)
  trans_df <- data.frame(
    family = fam_label(fam),
    upstream_cluster = unname(lab_up),
    cds_variant = unname(variants[shared_taxa]),
    downstream_cluster = unname(lab_down),
    stringsAsFactors = FALSE
  )
  trans_counts <- trans_df %>%
    count(family, upstream_cluster, cds_variant, downstream_cluster, name = "n_sequences")

  if (fam == FAMILIES[1]) {
    transition_all <- trans_counts
  } else {
    transition_all <- bind_rows(transition_all, trans_counts)
  }
}

if (!exists("concordance_all")) concordance_all <- data.frame(
  family = character(0), comparison = character(0), n_taxa = integer(0),
  ARI = numeric(0), NMI = numeric(0), VI_normalized = numeric(0)
)
if (!exists("transition_all")) transition_all <- data.frame(
  family = character(0), upstream_cluster = character(0),
  cds_variant = character(0), downstream_cluster = character(0), n_sequences = integer(0)
)

## Full-resolution outputs -- UNCHANGED, this is the authoritative numerical data
write.csv(concordance_all, file.path(DIR_COMPARISON, "cross_region_cluster_concordance.csv"),
          row.names = FALSE)
write.csv(transition_all, file.path(DIR_COMPARISON, "cluster_transition_counts.csv"),
          row.names = FALSE)

tree_check_df <- bind_rows(tree_check_rows)
write.csv(tree_check_df, file.path(DIR_TREECHECK, "tree_check_summary.csv"), row.names = FALSE)

processed_families <- names(all_family_results)
if (length(processed_families) == 0) {
  stop("No families were successfully processed. Check trees_checked/tree_check_summary.csv.")
}

## =========================================================================

## =========================================================================
## =========================================================================
## FIGURE 1: PATRISTIC-DISTANCE CLUSTER-NETWORK ARCHITECTURE
##
## FINAL LAYOUT: 3 columns x 5 rows
##   columns = Upstream | CDS | Downstream
##   rows    = a) IMP | b) KPC | c) NDM | d) TEM | e) SHV
##
## Visual encoding:
##   * One node = one final cluster.
##   * Node AREA is proportional to cluster size (n sequences), using one
##     common size scale across the entire figure.
##   * Node colour = cluster identity (C1, C2, ...).
##   * Node positions start from classical MDS of the complete median
##     inter-cluster patristic-distance matrix, followed only by display-only
##     node separation.
##   * Edges are the minimum-spanning backbone of the complete cluster-distance
##     matrix. Thicker edges indicate smaller median inter-cluster distances.
##   * Every cluster is labelled OUTSIDE its node using arrowed leader lines.
##     Labels are plain text (NO label boxes) and ggrepel prevents text-text
##     and text-node overlap.
##   * No edge-distance numbers are drawn.
##   * No rooting or evolutionary direction is implied.
##   * No caption is embedded in Figure 1; figure text belongs in the manuscript.
##
## All calculations and numerical outputs remain unchanged.
## =========================================================================
message("")
message("============================================")
message("FIGURE 1: 3 x 5 patristic cluster networks")
message("============================================")

## ---- cluster-level median patristic distances -----------------------------
get_cluster_distance_data <- function(fam, region) {
  D <- all_patristic[[fam]][[region]]
  labs <- all_family_results[[fam]][[region]]$labels
  taxa <- names(labs)

  clusters <- unique(unname(labs))
  clusters <- clusters[order(as.integer(sub("^C", "", clusters)))]

  n_vec <- as.integer(table(factor(unname(labs), levels = clusters)))
  names(n_vec) <- clusters

  M <- matrix(0, nrow = length(clusters), ncol = length(clusters),
              dimnames = list(clusters, clusters))
  within <- setNames(numeric(length(clusters)), clusters)

  for (i in seq_along(clusters)) {
    ci <- clusters[i]
    mi <- taxa[unname(labs) == ci]

    if (length(mi) <= 1L) {
      within[ci] <- 0
    } else {
      block <- D[mi, mi, drop = FALSE]
      vals <- block[upper.tri(block)]
      within[ci] <- median(vals, na.rm = TRUE)
    }

    if (i < length(clusters)) {
      for (j in (i + 1L):length(clusters)) {
        cj <- clusters[j]
        mj <- taxa[unname(labs) == cj]
        vals <- D[mi, mj, drop = FALSE]
        md <- median(vals, na.rm = TRUE)
        if (!is.finite(md)) md <- 0
        M[i, j] <- md
        M[j, i] <- md
      }
    }
  }
  diag(M) <- 0

  list(D = D, labels = labs, clusters = clusters, n = n_vec,
       within = within, between = M)
}

## ---- minimum spanning tree using only base R ------------------------------
mst_from_distance <- function(M) {
  k <- nrow(M)
  if (k < 2L) {
    return(data.frame(i = integer(0), j = integer(0), distance = numeric(0)))
  }

  inds <- which(upper.tri(M), arr.ind = TRUE)
  if (nrow(inds) == 0L) {
    return(data.frame(i = integer(0), j = integer(0), distance = numeric(0)))
  }

  edges <- data.frame(
    i = inds[, 1], j = inds[, 2], distance = M[inds], stringsAsFactors = FALSE
  )
  edges <- edges[is.finite(edges$distance), , drop = FALSE]
  if (!nrow(edges)) {
    return(data.frame(i = integer(0), j = integer(0), distance = numeric(0)))
  }

  edges <- edges[order(edges$distance, edges$i, edges$j), , drop = FALSE]
  parent <- seq_len(k)

  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }

  chosen <- logical(nrow(edges))
  n_chosen <- 0L
  for (r in seq_len(nrow(edges))) {
    a <- find_root(edges$i[r])
    b <- find_root(edges$j[r])
    if (a != b) {
      parent[a] <- b
      chosen[r] <- TRUE
      n_chosen <- n_chosen + 1L
      if (n_chosen == k - 1L) break
    }
  }

  edges[chosen, , drop = FALSE]
}

## ---- classical MDS starting coordinates -----------------------------------
cluster_mds <- function(M) {
  k <- nrow(M)
  if (k == 1L) {
    xy <- matrix(c(0, 0), ncol = 2)
    rownames(xy) <- rownames(M)
    return(xy)
  }

  if (k == 2L) {
    xy <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
    rownames(xy) <- rownames(M)
    return(xy)
  }

  cs <- tryCatch(
    stats::cmdscale(as.dist(M), k = 2, eig = TRUE, add = TRUE),
    error = function(e) NULL
  )

  if (is.null(cs) || is.null(cs$points) || nrow(cs$points) != k) {
    theta <- seq(0, 2 * pi, length.out = k + 1L)[-(k + 1L)]
    xy <- cbind(cos(theta), sin(theta))
  } else {
    xy <- as.matrix(cs$points)
    if (ncol(xy) < 2L) xy <- cbind(xy[, 1], rep(0, k))
    xy <- xy[, 1:2, drop = FALSE]
  }

  for (j in 1:2) {
    rg <- range(xy[, j], finite = TRUE)
    if (!all(is.finite(rg)) || diff(rg) == 0) {
      xy[, j] <- 0
    } else {
      xy[, j] <- 2 * (xy[, j] - mean(rg)) / diff(rg)
    }
  }

  rownames(xy) <- rownames(M)
  xy
}

## ---- display-only node separation -----------------------------------------
## This changes plotting coordinates only. It deliberately uses generous
## separation so node glyphs themselves never overlap.
separate_nodes <- function(xy, n_values, max_iter = 3000L) {
  xy <- as.matrix(xy)
  xy0 <- xy
  k <- nrow(xy)
  if (k <= 1L) return(xy)

  span_x <- diff(range(xy[, 1], finite = TRUE))
  span_y <- diff(range(xy[, 2], finite = TRUE))
  span <- max(span_x, span_y)
  if (!is.finite(span) || span <= 0) span <- 1

  nmax <- max(n_values, na.rm = TRUE)
  if (!is.finite(nmax) || nmax <= 0) nmax <- 1

  ## Radii used ONLY for collision detection. Actual plotted node area is set
  ## later from n and one common scale across all panels.
  radii <- span * (0.045 + 0.075 * sqrt(n_values / nmax))
  target_gap <- span * 0.035

  for (iter in seq_len(max_iter)) {
    moved <- FALSE
    disp <- matrix(0, nrow = k, ncol = 2)

    for (i in seq_len(k - 1L)) {
      for (j in (i + 1L):k) {
        v <- xy[j, ] - xy[i, ]
        d <- sqrt(sum(v^2))
        if (!is.finite(d) || d < 1e-10) {
          ang <- ((i * 137L + j * 67L) %% 360) * pi / 180
          v <- c(cos(ang), sin(ang))
          d <- 1
        }

        required <- radii[i] + radii[j] + target_gap
        if (d < required) {
          push <- 0.68 * (required - d)
          u <- v / d
          disp[i, ] <- disp[i, ] - push * u
          disp[j, ] <- disp[j, ] + push * u
          moved <- TRUE
        }
      }
    }

    if (!moved) break
    xy <- xy + disp
    xy <- 0.993 * xy + 0.007 * xy0
  }

  for (j in 1:2) {
    rg <- range(xy[, j], finite = TRUE)
    if (!all(is.finite(rg)) || diff(rg) == 0) {
      xy[, j] <- 0
    } else {
      xy[, j] <- 1.75 * (xy[, j] - mean(rg)) / diff(rg)
    }
  }
  xy
}

## ---- build all network data first -----------------------------------------
## Doing this before plotting lets us derive ONE node-size scale for the whole
## figure. With this common scale, circle area is proportional to sequence n
## across all 15 panels, rather than being re-normalized within each panel.
network_data <- list()
for (fam in FAMILIES) {
  for (region in REGIONS) {
    if (!fam %in% processed_families) next

    z <- get_cluster_distance_data(fam, region)
    xy0 <- cluster_mds(z$between)
    xy <- separate_nodes(xy0, as.numeric(z$n[rownames(xy0)]))

    nodes <- data.frame(
      family = fam,
      region = region,
      cluster = rownames(xy),
      x = xy[, 1], y = xy[, 2],
      n = as.numeric(z$n[rownames(xy)]),
      stringsAsFactors = FALSE
    )

    nodes$label <- paste0(
      nodes$cluster, " (n=",
      format(nodes$n, big.mark = ",", scientific = FALSE, trim = TRUE), ")"
    )

    mst <- mst_from_distance(z$between)
    if (nrow(mst) > 0L) {
      edges <- data.frame(
        family = fam, region = region,
        cluster_1 = z$clusters[mst$i], cluster_2 = z$clusters[mst$j],
        x = nodes$x[mst$i], y = nodes$y[mst$i],
        xend = nodes$x[mst$j], yend = nodes$y[mst$j],
        distance = mst$distance,
        stringsAsFactors = FALSE
      )
      d_rng <- range(edges$distance, finite = TRUE)
      if (!all(is.finite(d_rng)) || diff(d_rng) == 0) {
        edges$linewidth <- 1.05
      } else {
        edges$linewidth <- 0.60 + 1.35 *
          (d_rng[2] - edges$distance) / diff(d_rng)
      }
    } else {
      edges <- data.frame(
        family = character(0), region = character(0),
        cluster_1 = character(0), cluster_2 = character(0),
        x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
        distance = numeric(0), linewidth = numeric(0), stringsAsFactors = FALSE
      )
    }

    ## Exact cluster-level quantities underlying each plotted network.
    write.csv(as.data.frame(z$between),
      file.path(DIR_CLUSTERING,
        sprintf("%s_%s_cluster_patristic_distance_matrix.csv", fam, region)),
      row.names = TRUE)
    write.csv(nodes[, c("cluster", "n")],
      file.path(DIR_CLUSTERING,
        sprintf("%s_%s_cluster_network_nodes.csv", fam, region)),
      row.names = FALSE)

    mst_out <- if (nrow(mst) > 0L) {
      data.frame(
        cluster_1 = z$clusters[mst$i],
        cluster_2 = z$clusters[mst$j],
        median_patristic_distance = mst$distance,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        cluster_1 = character(0), cluster_2 = character(0),
        median_patristic_distance = numeric(0), stringsAsFactors = FALSE
      )
    }
    write.csv(mst_out,
      file.path(DIR_CLUSTERING,
        sprintf("%s_%s_cluster_network_edges.csv", fam, region)),
      row.names = FALSE)

    network_data[[paste(fam, region, sep = "__")]] <- list(
      nodes = nodes, edges = edges
    )
  }
}

if (!length(network_data)) stop("No network data available for Figure 1.")

GLOBAL_MAX_CLUSTER_N <- max(
  unlist(lapply(network_data, function(z) z$nodes$n)), na.rm = TRUE
)
if (!is.finite(GLOBAL_MAX_CLUSTER_N) || GLOBAL_MAX_CLUSTER_N <= 0) {
  GLOBAL_MAX_CLUSTER_N <- 1
}

cluster_colours <- setNames(
  CLUSTER_PALETTE,
  paste0("C", seq_along(CLUSTER_PALETTE))
)

make_fig1_panel <- function(fam, region) {
  key <- paste(fam, region, sep = "__")
  nd <- network_data[[key]]
  nodes <- nd$nodes
  edges <- nd$edges
  cl <- all_family_results[[fam]][[region]]

  xr <- range(c(nodes$x, edges$x, edges$xend), finite = TRUE)
  yr <- range(c(nodes$y, edges$y, edges$yend), finite = TRUE)
  sx <- diff(xr); sy <- diff(yr)
  if (!is.finite(sx) || sx <= 0) sx <- 1
  if (!is.finite(sy) || sy <= 0) sy <- 1

  ## Enough surrounding space for labels, but still compact.
  pad_x <- 0.58 * max(sx, sy)
  pad_y <- 0.58 * max(sx, sy)

  fam_index <- match(fam, FAMILIES)
  fam_letter <- letters[fam_index]
  title_txt <- sprintf("%s) %s", fam_letter, fam_label(fam))
  subtitle_txt <- sprintf(
    "%s | n = %d | k = %d | silhouette = %s",
    region_label(region), cl$n, cl$k,
    ifelse(is.na(cl$silhouette), "NA", sprintf("%.3f", cl$silhouette))
  )

  ## One common size scale: plotted circle area is proportional to n.
  ## scale_size_area() preserves this relationship; the same limits are used
  ## for every panel, so n remains visually comparable across the figure.
  p <- ggplot() +
    geom_segment(
      data = edges,
      aes(x = x, y = y, xend = xend, yend = yend, linewidth = linewidth),
      colour = "grey62", lineend = "round", inherit.aes = FALSE
    ) +
    scale_linewidth_identity() +
    geom_point(
      data = nodes,
      aes(x = x, y = y, size = n, fill = cluster),
      shape = 21, colour = "grey10", stroke = 0.65, alpha = 0.98,
      inherit.aes = FALSE
    ) +
    scale_size_area(
      max_size = 13.0,
      limits = c(0, GLOBAL_MAX_CLUSTER_N),
      breaks = NULL,
      guide = "none"
    ) +
    scale_fill_manual(values = cluster_colours, guide = "none") +
    ## Plain text labels with a very subtle white knockout behind the glyphs.
    ## This is NOT a label box; it only prevents an underlying connector from
    ## visually running through the letterforms.
    ggrepel::geom_text_repel(
      data = nodes,
      aes(x = x, y = y, label = label),
      colour = "grey12",
      size = 3.05,
      bg.color = "white",
      bg.r = 0.16,
      fontface = "plain",
      box.padding = 1.05,
      point.padding = 1.05,
      force = 14,
      force_pull = 0.02,
      max.iter = 100000,
      max.time = 30,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.colour = "grey48",
      segment.size = 0.32,
      segment.alpha = 0.92,
      arrow = grid::arrow(
        length = grid::unit(0.10, "cm"),
        type = "closed"
      ),
      seed = 42 + match(fam, FAMILIES) * 100 + match(region, REGIONS),
      direction = "both",
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    labs(title = title_txt, subtitle = subtitle_txt) +
    coord_cartesian(
      xlim = c(xr[1] - pad_x, xr[2] + pad_x),
      ylim = c(yr[1] - pad_y, yr[2] + pad_y),
      clip = "off",
      expand = FALSE
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(
        size = 11.5, face = "plain", hjust = 0, colour = "grey10",
        margin = margin(b = 1.2)
      ),
      plot.subtitle = element_text(
        size = 7.6, hjust = 0, colour = "grey38",
        margin = margin(b = 2)
      ),
      plot.margin = margin(7, 9, 7, 7),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )

  p
}

fig1_panels <- list()
for (fam in FAMILIES) {
  for (region in REGIONS) {
    key <- paste(fam, region, sep = "__")
    if (fam %in% processed_families && key %in% names(network_data)) {
      fig1_panels[[length(fig1_panels) + 1L]] <- make_fig1_panel(fam, region)
    }
  }
}

## IMPORTANT: no caption and no facet-strip boxes. The row identity is carried
## directly by a) IMP, b) KPC, c) NDM, d) TEM, e) SHV in each panel title, while
## the genomic region appears as the subtitle.
figure1 <- patchwork::wrap_plots(
  fig1_panels,
  ncol = 3,
  nrow = 5,
  widths = rep(1, 3),
  heights = rep(1, 5),
  guides = "keep"
)

## Large publication canvas; not constrained to A4.
FIG1_WIDTH_IN  <- 18
FIG1_HEIGHT_IN <- 20

fig1_pdf <- file.path(DIR_FIGURES, "Figure1_patristic_cluster_networks_3x5.pdf")
fig1_png <- file.path(DIR_FIGURES, "Figure1_patristic_cluster_networks_3x5.png")
fig1_svg <- file.path(DIR_FIGURES, "Figure1_patristic_cluster_networks_3x5.svg")

pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
ggsave(
  fig1_pdf, figure1,
  width = FIG1_WIDTH_IN, height = FIG1_HEIGHT_IN,
  units = "in", device = pdf_device
)
ggsave(
  fig1_png, figure1,
  width = FIG1_WIDTH_IN, height = FIG1_HEIGHT_IN,
  units = "in", dpi = 500, bg = "white"
)
if (HAVE_SVGLITE) {
  ggsave(
    fig1_svg, figure1,
    width = FIG1_WIDTH_IN, height = FIG1_HEIGHT_IN,
    units = "in"
  )
} else {
  message("svglite not available: skipping SVG export for Figure 1 (PDF/PNG saved).")
}
message("Figure 1: 3 x 5 upstream/CDS/downstream patristic cluster networks saved (labels repel all text and are masked from connectors).")
message("  Connector thickness: thicker = smaller median inter-cluster patristic distance (within each panel).")

## FIGURE 2: Upstream cluster -> CDS variant -> Downstream cluster
##        Only labels that fit inside their own nodes are shown.
##        Small nodes are intentionally left unlabeled; no outside text/arrows.
## =========================================================================
message("Figure 2: building publication-ready cluster correspondence panels...")

## ---- build a per-family, DISPLAY-ONLY aggregated transition table --------
build_display_transitions <- function(trans_fam) {
  if (nrow(trans_fam) == 0) {
    return(trans_fam %>% mutate(cds_variant_display = character(0)))
  }
  variant_totals <- trans_fam %>%
    group_by(cds_variant) %>%
    summarise(total_n = sum(n_sequences), .groups = "drop")
  keep_variants <- variant_totals$cds_variant[
    variant_totals$total_n >= DISPLAY_VARIANT_MIN_N & variant_totals$cds_variant != "Unknown"
  ]
  trans_fam %>%
    mutate(cds_variant_display = case_when(
      cds_variant == "Unknown" ~ "Unknown",
      cds_variant %in% keep_variants ~ cds_variant,
      TRUE ~ "Other variants"
    ))
}

display_trans_by_family <- list()
for (fam in processed_families) {
  trans_fam <- transition_all %>% filter(family == fam_label(fam))
  display_trans_by_family[[fam]] <- build_display_transitions(trans_fam)
}

## ---- deterministic GLOBAL variant colour map (built once, used everywhere) ----
all_display_variants <- bind_rows(display_trans_by_family) %>%
  filter(!cds_variant_display %in% c("Unknown", "Other variants")) %>%
  group_by(cds_variant_display) %>%
  summarise(total_n = sum(n_sequences), .groups = "drop") %>%
  arrange(desc(total_n))

if (nrow(all_display_variants) > 0) {
  variant_color_map <- setNames(
    rep(VARIANT_PALETTE, length.out = nrow(all_display_variants)),
    all_display_variants$cds_variant_display
  )
} else {
  variant_color_map <- character(0)
}
variant_color_map[["Unknown"]] <- OTHER_GREY
variant_color_map[["Other variants"]] <- OTHER_GREY

order_for_stack <- function(heights, bottom_keys = character(0)) {
  nm <- names(heights)
  is_bottom <- nm %in% bottom_keys
  bottom_part <- heights[is_bottom]
  if (length(bottom_part) > 0) bottom_part <- bottom_part[order(names(bottom_part))]
  rest_part <- heights[!is_bottom]
  if (length(rest_part) > 0) rest_part <- sort(rest_part)
  c(bottom_part, rest_part)
}

make_node_positions <- function(heights) {
  heights[is.na(heights)] <- 0
  total <- sum(heights)
  if (total == 0) total <- 1
  gap_frac <- 0.15
  n_nodes <- length(heights)
  gap <- total * gap_frac / max(n_nodes - 1, 1)
  ytop <- numeric(n_nodes); ybot <- numeric(n_nodes)
  cursor <- 0
  for (i in seq_len(n_nodes)) {
    ytop[i] <- cursor + heights[i]
    ybot[i] <- cursor
    cursor <- cursor + heights[i] + gap
  }
  data.frame(node = names(heights), ybot = ybot, ytop = ytop,
             ymid = (ybot + ytop) / 2, gap = gap, stringsAsFactors = FALSE)
}

# Figure 2 label policy:
#   Put the node/allele text INSIDE the block whenever that block is large enough
#   to contain it.  There are NO outside labels and NO leader arrows.  Small or
#   narrow blocks are intentionally left unlabelled.
F2_INSIDE_MIN_FRAC <- 0.050

build_inside_only_labels <- function(pos, x_center) {
  total <- sum(pos$ytop - pos$ybot)
  if (!is.finite(total) || total <= 0) total <- 1

  ## Only label a node when its vertical height is sufficient for the text.
  ## A lower threshold allows C1/C2/allele labels to appear whenever the
  ## node has enough vertical room; genuinely tiny blocks remain blank.
  pos <- pos %>% mutate(
    height_frac = (ytop - ybot) / total,
    show_label = height_frac >= F2_INSIDE_MIN_FRAC,
    placement = ifelse(show_label, "inside", "suppress")
  )

  ## IMPORTANT: label_x is the x-coordinate of the column, not the node y-midpoint.
  ## The previous version accidentally assigned ymid to label_x, which could place
  ## the text away from the node and make useful labels disappear from view.
  pos$label_x <- x_center
  pos$label_y <- pos$ymid
  pos
}

build_flow_panel <- function(fam_key, trans_display) {
  fam_lab <- fam_label(fam_key)
  if (nrow(trans_display) == 0) {
    return(ggplot() + theme_void() + labs(title = fam_lab, subtitle = "No data"))
  }

  node_height <- function(data, col_vals) {
    tot <- data %>% group_by(.data[[col_vals]]) %>%
      summarise(n = sum(n_sequences), .groups = "drop")
    setNames(tot$n, tot[[col_vals]])
  }

  up_h   <- order_for_stack(node_height(trans_display, "upstream_cluster"))
  var_h  <- order_for_stack(node_height(trans_display, "cds_variant_display"),
                             bottom_keys = c("Other variants", "Unknown"))
  down_h <- order_for_stack(node_height(trans_display, "downstream_cluster"))

  pos_up   <- make_node_positions(up_h)
  pos_var  <- make_node_positions(var_h)
  pos_down <- make_node_positions(down_h)

  node_df <- bind_rows(
    pos_up   %>% mutate(x = 1),
    pos_var  %>% mutate(x = 2),
    pos_down %>% mutate(x = 3)
  )

  make_flow_polygons <- function(trans_sub, left_col, right_col, color_col,
                                 pos_left, pos_right, x_left, x_right) {
    left_order  <- pos_left$node[order(pos_left$ymid)]
    right_order <- pos_right$node[order(pos_right$ymid)]
    left_cursor  <- setNames(pos_left$ybot, pos_left$node)
    right_cursor <- setNames(pos_right$ybot, pos_right$node)

    trans_sub <- trans_sub %>%
      mutate(.lrank = match(.data[[left_col]], left_order),
             .rrank = match(.data[[right_col]], right_order)) %>%
      arrange(.lrank, .rrank)

    poly_list <- list()
    for (i in seq_len(nrow(trans_sub))) {
      lval <- as.character(trans_sub[[left_col]][i])
      rval <- as.character(trans_sub[[right_col]][i])
      cval <- as.character(trans_sub[[color_col]][i])
      h <- trans_sub$n_sequences[i]
      y0L <- left_cursor[lval]; y1L <- y0L + h
      y0R <- right_cursor[rval]; y1R <- y0R + h
      left_cursor[lval] <- y1L
      right_cursor[rval] <- y1R
      n_curve <- 24
      xs <- seq(x_left, x_right, length.out = n_curve)
      ease <- (1 - cos(seq(0, pi, length.out = n_curve))) / 2
      y_bottom <- y0L + (y0R - y0L) * ease
      y_top    <- y1L + (y1R - y1L) * ease
      poly_id <- paste(lval, rval, x_left, sep = "__")
      poly_list[[length(poly_list) + 1]] <- data.frame(
        x = c(xs, rev(xs)), y = c(y_bottom, rev(y_top)),
        group = poly_id, fill_key = cval, stringsAsFactors = FALSE
      )
    }
    bind_rows(poly_list)
  }

  agg_left <- trans_display %>%
    count(upstream_cluster, cds_variant_display, wt = n_sequences, name = "n_sequences")
  agg_right <- trans_display %>%
    count(cds_variant_display, downstream_cluster, wt = n_sequences, name = "n_sequences")

  flow1 <- make_flow_polygons(agg_left, "upstream_cluster", "cds_variant_display",
                              "cds_variant_display", pos_up, pos_var, 1, 2)
  flow2 <- make_flow_polygons(agg_right, "cds_variant_display", "downstream_cluster",
                              "cds_variant_display", pos_var, pos_down, 2, 3)

  # ONLY labels that physically fit INSIDE their own node are drawn.
  # Small nodes are intentionally left blank; no external text and no leader
  # arrows are used anywhere in Figure 2.
  pos_up_lab   <- build_inside_only_labels(pos_up,   x_center = 1)
  pos_var_lab  <- build_inside_only_labels(pos_var,  x_center = 2)
  pos_down_lab <- build_inside_only_labels(pos_down, x_center = 3)

  label_all <- bind_rows(pos_up_lab %>% mutate(x = 1),
                          pos_var_lab %>% mutate(x = 2),
                          pos_down_lab %>% mutate(x = 3)) %>%
    filter(show_label)

  p <- ggplot() +
    geom_polygon(data = flow1,
                 aes(x = x, y = y, group = group, fill = fill_key), alpha = 0.58) +
    geom_polygon(data = flow2,
                 aes(x = x, y = y, group = group, fill = fill_key), alpha = 0.44) +
    geom_rect(
      data = node_df,
      aes(xmin = x - NODE_HALF_WIDTH, xmax = x + NODE_HALF_WIDTH,
          ymin = ybot, ymax = ytop),
      fill = "grey25"
    )

  if (nrow(label_all) > 0) {
    p <- p + geom_text(
      data = label_all,
      aes(x = label_x, y = label_y, label = node),
      inherit.aes = FALSE,
      size = 1.62,
      fontface = "bold",
      colour = "white"
    )
  }

  y_min <- min(node_df$ybot, na.rm = TRUE)
  y_max <- max(node_df$ytop, na.rm = TRUE)
  y_pad <- max(0.035 * (y_max - y_min), 0.35)

  p +
    scale_fill_manual(values = variant_color_map, guide = "none") +
    scale_x_continuous(
      limits = c(0.42, 3.58),
      breaks = c(1, 2, 3),
      labels = c("Upstream\ncluster", "CDS\nvariant", "Downstream\ncluster"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(limits = c(y_min - y_pad, y_max + y_pad),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = fam_lab, x = NULL, y = NULL) +
    theme_void(base_size = 7) +
    theme(
      plot.title = element_text(size = 7.4, face = "bold", hjust = 0,
                                margin = margin(b = 1.5)),
      axis.text.x = element_text(size = 5.5, color = "grey25",
                                 lineheight = 0.92, margin = margin(t = 2)),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(2, 5, 4, 5)
    )
}

fig2_panels <- list()
for (fam in FAMILIES) {
  if (fam %in% processed_families) {
    fig2_panels[[fam]] <- build_flow_panel(fam, display_trans_by_family[[fam]])
  } else {
    fig2_panels[[fam]] <- ggplot() + theme_void() +
      labs(title = fam_label(fam), subtitle = "skipped (see tree_check_summary.csv)")
  }
}

fig2_layout <- fig2_panels[FAMILIES]

figure2 <- wrap_plots(fig2_layout, ncol = 1, nrow = 5,
                      heights = rep(1, 5), guides = "keep") +
  plot_annotation(
    theme = theme(
      plot.caption = element_text(size = 5.6, hjust = 0,
                                  lineheight = 1.05,
                                  margin = margin(t = 6, r = 2, b = 1, l = 2)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 3, r = 4, b = 4, l = 4)
    )
  )

FIG2_WIDTH_IN  <- 11
FIG2_HEIGHT_IN <- 16.5

ggsave(file.path(DIR_FIGURES, "Figure2_cluster_membership_correspondence.pdf"),
       figure2, width = FIG2_WIDTH_IN, height = FIG2_HEIGHT_IN, units = "in", device = pdf_device)
ggsave(file.path(DIR_FIGURES, "Figure2_cluster_membership_correspondence.png"),
       figure2, width = FIG2_WIDTH_IN, height = FIG2_HEIGHT_IN, units = "in", dpi = 400, bg = "white")
if (HAVE_SVGLITE) {
  ggsave(file.path(DIR_FIGURES, "Figure2_cluster_membership_correspondence.svg"),
         figure2, width = FIG2_WIDTH_IN, height = FIG2_HEIGHT_IN, units = "in")
} else {
  message("svglite not available: skipping SVG export for Figure 2 (PDF/PNG saved).")
}
message("Figure 2: saved (publication-scale PDF/PNG", ifelse(HAVE_SVGLITE, "/SVG", ""), ")")

## =========================================================================
## FINAL SUMMARY
## =========================================================================
message("")
message("============================================")
message("FINAL CLUSTERING SUMMARY")
message("============================================")
final_summary_df <- bind_rows(final_summary_rows)
if (nrow(final_summary_df) > 0) {
  for (i in seq_len(nrow(final_summary_df))) {
    r <- final_summary_df[i, ]
    message(sprintf("%-4s | %-10s | n=%-5d | %-28s | k=%-3d | silhouette=%s",
                     r$family, r$region, r$n, r$method, r$k,
                     ifelse(is.na(r$silhouette), "NA", sprintf("%.4f", r$silhouette))))
  }
}
write.csv(final_summary_df, file.path(OUT_ROOT, "final_clustering_summary.csv"), row.names = FALSE)

message("")
message("============================================")
message("FIGURES")
message("============================================")
message(sprintf("Figure 1: %s", file.path(DIR_FIGURES, "Figure1_patristic_cluster_networks_3x5.pdf")))
message(sprintf("Figure 2: %s", file.path(DIR_FIGURES, "Figure2_cluster_membership_correspondence.pdf")))

writeLines(capture.output(sessionInfo()), file.path(OUT_ROOT, "session_info.txt"))

message("")
message("DONE.")