#!/usr/bin/env Rscript

# =============================================================================
# TOPOLOGICAL DISCORDANCE ANALYSIS
# Upstream vs CDS vs Downstream
#
# Families:
#   IMP, NDM, KPC, TEM, SHV
#
# Current tree directory:
#   trees
#
# Expected filenames:
#   imp_upstream.treefile
#   imp_cds.treefile
#   imp_downstream.treefile
#
# Main analysis:
#   1. Tree discovery and classification
#   2. Tip-label normalization
#   3. Common-taxon audit
#   4. Tree pruning
#   5. Normalized Clustering Information Distance
#   6. Normalized RF distance
#   7. Quartet distance
#   8. Permutation null model
#   9. Matched-taxon bootstrap
#  10. Upstream-Downstream vs flank comparisons
#  11. Tanglegrams
#  12. Main quantitative plot
# =============================================================================


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i >= length(args)) stop(sprintf("Missing value for %s", flag))
  args[[i + 1L]]
}

CONFIG <- list(

  # ---------------------------------------------------------------------------
  # TREE DIRECTORY
  # ---------------------------------------------------------------------------
  tree_dir =
    arg_value("--tree-dir", "trees"),


  # ---------------------------------------------------------------------------
  # FAMILY PATTERNS
  # ---------------------------------------------------------------------------
  family_patterns = c(
    IMP = "(?i)imp",
    NDM = "(?i)ndm",
    KPC = "(?i)kpc",
    TEM = "(?i)tem",
    SHV = "(?i)shv"
  ),


  # ---------------------------------------------------------------------------
  # REGION PATTERNS
  #
  # IMPORTANT:
  # Do NOT use \\bcds\\b because "_" is a word character in regex.
  #
  # These patterns detect:
  #
  # imp_upstream.treefile
  # imp_cds.treefile
  # imp_downstream.treefile
  # ---------------------------------------------------------------------------
  region_patterns = c(

    upstream =
      "(?i)(^|[_\\.-])up(stream)?([_\\.-]|$)",

    cds =
      "(?i)(^|[_\\.-])cds([_\\.-]|$)",

    downstream =
      "(?i)(^|[_\\.-])down(stream)?([_\\.-]|$)"
  ),


  # ---------------------------------------------------------------------------
  # TREE EXTENSIONS
  # ---------------------------------------------------------------------------
  tree_file_ext =
    "\\.(treefile|contree|nwk|tre|tree)$",


  # ---------------------------------------------------------------------------
  # IF .contree AND .treefile BOTH EXIST
  # ---------------------------------------------------------------------------
  prefer_contree = TRUE,


  # ---------------------------------------------------------------------------
  # LOW SUPPORT COLLAPSE
  # ---------------------------------------------------------------------------
  support_collapse_threshold = 70,


  # ---------------------------------------------------------------------------
  # OUTPUT
  # ---------------------------------------------------------------------------
  out_dir =
    arg_value("--output-dir", "results/topological_discordance"),


  # ---------------------------------------------------------------------------
  # STATISTICS
  # ---------------------------------------------------------------------------
  n_null_perm = 1000,

  n_bootstrap = 200,

  matched_n_cap = 100,

  quartet_leaf_cap = 500,

  seed = 20260820,


  # ---------------------------------------------------------------------------
  # TANGLEGRAM
  # ---------------------------------------------------------------------------
  tanglegram_max_tips = 50,


  # ---------------------------------------------------------------------------
  # FAMILY ORDER
  # ---------------------------------------------------------------------------
  family_levels =
    c(
      "IMP",
      "NDM",
      "KPC",
      "TEM",
      "SHV"
    ),


  # ---------------------------------------------------------------------------
  # FAMILY COLOURS
  # ---------------------------------------------------------------------------
  palette = c(
    IMP = "#E69F00",
    NDM = "#56B4E9",
    KPC = "#009E73",
    TEM = "#D55E00",
    SHV = "#0072B2"
  )
)


# =============================================================================
# 1. PACKAGE CHECK
# =============================================================================

required_pkgs <- c(
  "ape",
  "phangorn",
  "TreeDist",
  "Quartet",
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "phytools"
)


missing_pkgs <-
  required_pkgs[
    !required_pkgs %in% rownames(installed.packages())
  ]


if (length(missing_pkgs) > 0) {

  stop(
    paste0(
      "\nMissing R packages:\n",
      paste(
        missing_pkgs,
        collapse = ", "
      ),
      "\n\nInstall these packages first.\n"
    )
  )
}


suppressPackageStartupMessages({

  library(ape)

  library(phangorn)

  library(TreeDist)

  library(Quartet)

  library(ggplot2)

  library(dplyr)

  library(tidyr)

  library(stringr)

  library(purrr)

  library(phytools)

})


set.seed(
  CONFIG$seed
)


dir.create(
  CONFIG$out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 2. CLASSIFY TREE FILES
# =============================================================================

classify_files <- function(cfg) {

  message("")
  message("============================================================")
  message("SEARCHING FOR TREE FILES")
  message("============================================================")
  message("")
  message(
    "Tree directory: ",
    cfg$tree_dir
  )
  message("")


  if (!dir.exists(cfg$tree_dir)) {

    stop(
      "\nTree directory does not exist:\n",
      cfg$tree_dir,
      "\n"
    )
  }


  files <-
    list.files(
      cfg$tree_dir,
      full.names = TRUE,
      recursive = TRUE
    )


  files <-
    files[
      grepl(
        cfg$tree_file_ext,
        files,
        perl = TRUE
      )
    ]


  if (length(files) == 0) {

    stop(
      "\nNo tree files found in:\n",
      cfg$tree_dir,
      "\n"
    )
  }


  message(
    "Found ",
    length(files),
    " tree files."
  )


  # ---------------------------------------------------------------------------
  # Classify one file
  # ---------------------------------------------------------------------------

  classify_one <- function(f) {

    base <-
      basename(f)


    # Family
    fam_hits <-
      names(cfg$family_patterns)[
        purrr::map_lgl(
          cfg$family_patterns,
          function(pattern) {

            grepl(
              pattern,
              base,
              perl = TRUE
            )
          }
        )
      ]


    # Region
    reg_hits <-
      names(cfg$region_patterns)[
        purrr::map_lgl(
          cfg$region_patterns,
          function(pattern) {

            grepl(
              pattern,
              base,
              perl = TRUE
            )
          }
        )
      ]


    tibble(

      path = f,

      filename = base,

      family =
        if (length(fam_hits) > 0)
          fam_hits[1]
      else
        NA_character_,

      region =
        if (length(reg_hits) > 0)
          reg_hits[1]
      else
        NA_character_,

      is_contree =
        grepl(
          "\\.contree$",
          f,
          perl = TRUE
        )
    )
  }


  manifest <-
    purrr::map_dfr(
      files,
      classify_one
    )


  # ---------------------------------------------------------------------------
  # PRINT ALL CLASSIFIED FILES
  # ---------------------------------------------------------------------------

  message("")
  message("Detected tree files:")
  message("")


  print(
    manifest %>%
      select(
        filename,
        family,
        region
      ) %>%
      arrange(
        family,
        region,
        filename
      )
  )


  # ---------------------------------------------------------------------------
  # UNCLASSIFIED FILES
  # ---------------------------------------------------------------------------

  unclassified <-
    manifest %>%
    filter(
      is.na(family) |
      is.na(region)
    )


  if (nrow(unclassified) > 0) {

    warning(
      paste0(
        "\nUnclassified tree files:\n",
        paste(
          unclassified$filename,
          collapse = "\n"
        )
      )
    )
  }


  manifest <-
    manifest %>%
    filter(
      !is.na(family),
      !is.na(region)
    )


  # ---------------------------------------------------------------------------
  # HANDLE DUPLICATES
  # ---------------------------------------------------------------------------

  manifest <-
    manifest %>%
    group_by(
      family,
      region
    ) %>%
    group_modify(
      function(d, key) {

        if (nrow(d) <= 1) {

          return(d)
        }


        if (
          cfg$prefer_contree &&
          any(d$is_contree)
        ) {

          return(
            d %>%
              filter(is_contree)
          )
        }


        warning(
          "Multiple files detected for ",
          key$family,
          " / ",
          key$region,
          ". Using first file: ",
          d$filename[1]
        )


        d[1, , drop = FALSE]
      }
    ) %>%
    ungroup()


  # ---------------------------------------------------------------------------
  # EXPECTED FAMILY x REGION COMBINATIONS
  # ---------------------------------------------------------------------------

  expected <-
    expand.grid(
      family =
        cfg$family_levels,

      region =
        c(
          "upstream",
          "cds",
          "downstream"
        ),

      stringsAsFactors = FALSE
    )


  missing_combo <-
    expected %>%
    anti_join(
      manifest,
      by =
        c(
          "family",
          "region"
        )
    )


  if (nrow(missing_combo) > 0) {

    warning(
      paste0(
        "\nMissing family/region combinations:\n",
        paste(
          sprintf(
            "  %s - %s",
            missing_combo$family,
            missing_combo$region
          ),
          collapse = "\n"
        )
      )
    )
  }


  # ---------------------------------------------------------------------------
  # SAVE MANIFEST
  # ---------------------------------------------------------------------------

  write.csv(
    manifest,
    file.path(
      cfg$out_dir,
      "tree_manifest.csv"
    ),
    row.names = FALSE
  )


  message("")
  message(
    "Tree classification completed."
  )
  message("")


  manifest
}


# =============================================================================
# 3. NORMALIZE TIP LABELS
# =============================================================================

normalize_tip_labels <- function(tr) {

  # Example:
  #
  # IMP-1|CP141568
  #
  # becomes:
  #
  # CP141568


  tr$tip.label <-
    sub(
      ".*\\|",
      "",
      tr$tip.label
    )


  tr$tip.label <-
    trimws(
      tr$tip.label
    )


  tr
}


# =============================================================================
# 4. COLLAPSE LOW SUPPORT
# =============================================================================

collapse_low_support <- function(
  tr,
  threshold
) {

  if (is.null(threshold)) {

    return(tr)
  }


  if (
    is.null(tr$node.label) ||
    length(tr$node.label) == 0
  ) {

    return(tr)
  }


  supp <-
    suppressWarnings(
      as.numeric(
        tr$node.label
      )
    )


  if (
    all(
      is.na(supp)
    )
  ) {

    return(tr)
  }


  finite_supp <-
    supp[
      is.finite(supp)
    ]


  # Convert 0-1 support to 0-100
  if (
    length(finite_supp) > 0 &&
    all(
      finite_supp <= 1
    )
  ) {

    supp <-
      supp * 100
  }


  node_numbers <-
    which(
      is.finite(supp) &
      supp < threshold
    ) +
    Ntip(tr)


  if (
    length(node_numbers) == 0
  ) {

    return(tr)
  }


  low_edges <-
    which(
      tr$edge[, 2] %in%
        node_numbers
    )


  if (
    length(low_edges) > 0
  ) {

    if (
      is.null(
        tr$edge.length
      )
    ) {

      tr$edge.length <-
        rep(
          1,
          nrow(
            tr$edge
          )
        )
    }


    tr$edge.length[
      low_edges
    ] <- 0
  }


  tr <-
    di2multi(
      tr,
      tol = 1e-8
    )


  tr
}


# =============================================================================
# 5. LOAD FAMILY TREES
# =============================================================================

load_family_trees <- function(
  manifest,
  family,
  cfg
) {

  sub <-
    manifest %>%
    filter(
      family == !!family
    )


  regions <-
    c(
      "upstream",
      "cds",
      "downstream"
    )


  trs <-
    setNames(
      vector(
        "list",
        length(regions)
      ),
      regions
    )


  for (r in regions) {

    row <-
      sub %>%
      filter(
        region == r
      )


    if (
      nrow(row) == 0
    ) {

      trs[[r]] <- NULL

      next
    }


    message(
      "Loading ",
      family,
      " ",
      r,
      ": ",
      basename(
        row$path[1]
      )
    )


    tr <-
      ape::read.tree(
        row$path[1]
      )


    tr <-
      normalize_tip_labels(
        tr
      )


    # -------------------------------------------------------------------------
    # Duplicate labels
    # -------------------------------------------------------------------------

    dup <-
      duplicated(
        tr$tip.label
      )


    if (
      any(dup)
    ) {

      warning(
        family,
        " ",
        r,
        ": duplicate accession labels detected."
      )


      # Keep first occurrence
      keep <-
        !duplicated(
          tr$tip.label
        )


      tr <-
        ape::keep.tip(
          tr,
          tr$tip.label[keep]
        )
    }


    tr <-
      ape::unroot(
        tr
      )


    tr <-
      collapse_low_support(
        tr,
        cfg$support_collapse_threshold
      )


    trs[[r]] <-
      tr
  }


  trs
}


# =============================================================================
# 6. AUDIT AND PRUNE
# =============================================================================

audit_and_prune <- function(
  trs,
  family
) {

  labels <-
    purrr::map(
      trs,
      function(tr) {

        if (
          is.null(tr)
        ) {

          return(
            character(0)
          )
        }


        tr$tip.label
      }
    )


  nonempty <-
    labels[
      purrr::map_lgl(
        labels,
        function(x) {

          length(x) > 0
        }
      )
    ]


  if (
    length(nonempty) < 3
  ) {

    return(
      list(
        audit = tibble(),
        pruned = trs,
        n_common = 0
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Common taxa across ALL THREE
  # ---------------------------------------------------------------------------

  all_present <-
    Reduce(
      intersect,
      nonempty
    )


  # ---------------------------------------------------------------------------
  # Pairwise audit
  # ---------------------------------------------------------------------------

  audit_rows <-
    list()


  regions <-
    names(trs)


  for (
    i in seq_along(regions)
  ) {

    for (
      j in seq_along(regions)
    ) {

      if (
        j <= i
      )
        next


      r1 <-
        regions[i]

      r2 <-
        regions[j]


      if (
        is.null(
          trs[[r1]]
        ) ||
        is.null(
          trs[[r2]]
        )
      )
        next


      shared <-
        intersect(
          labels[[r1]],
          labels[[r2]]
        )


      only1 <-
        setdiff(
          labels[[r1]],
          labels[[r2]]
        )


      only2 <-
        setdiff(
          labels[[r2]],
          labels[[r1]]
        )


      # -----------------------------------------------------------------------
      # CORRECT R LIST INDEXING
      # -----------------------------------------------------------------------

    audit_rows[[paste(r1, r2, sep = "_vs_")]] <- tibble(
  family = family,

  comparison =
    paste0(
      r1,
      "_vs_",
      r2
    ),

  n_region1 =
    length(
      labels[[r1]]
    ),

  n_region2 =
    length(
      labels[[r2]]
    ),

  n_shared =
    length(
      shared
    ),

  n_only_region1 =
    length(
      only1
    ),

  n_only_region2 =
    length(
      only2
    )
)
    }
  }


  audit <-
    bind_rows(
      audit_rows
    )


  # ---------------------------------------------------------------------------
  # PRUNE ALL TREES TO COMMON TAXA
  # ---------------------------------------------------------------------------

  pruned <-
    purrr::map(
      trs,
      function(tr) {

        if (
          is.null(tr)
        )
          return(NULL)


        if (
          length(all_present) < 4
        )
          return(NULL)


        ape::keep.tip(
          tr,
          all_present
        )
      }
    )


  # ---------------------------------------------------------------------------
  # FORCE SAME TAXON ORDER
  # ---------------------------------------------------------------------------

  if (
    !is.null(
      pruned$upstream
    )
  ) {

    tip_order <-
      pruned$upstream$tip.label


    pruned <-
      purrr::map(
        pruned,
        function(tr) {

          if (
            is.null(tr)
          )
            return(NULL)


          tr
        }
      )
  }


  list(
    audit = audit,

    pruned = pruned,

    n_common =
      length(
        all_present
      )
  )
}


# =============================================================================
# 7. NORMALIZED CLUSTERING INFORMATION DISTANCE
# =============================================================================

grf_distance <- function(
  tr1,
  tr2
) {

  result <-
    tryCatch(

      TreeDist::ClusteringInfoDistance(
        tr1,
        tr2,
        normalize = TRUE
      ),

      error = function(e) {

        warning(
          "ClusteringInfoDistance failed: ",
          conditionMessage(e)
        )

        NA_real_
      }
    )


  as.numeric(
    result
  )
}


# =============================================================================
# 8. NORMALIZED RF DISTANCE
# =============================================================================

nrf_distance <- function(
  tr1,
  tr2
) {

  n <-
    Ntip(
      tr1
    )


  if (
    n <= 3
  )
    return(NA_real_)


  raw <-
    tryCatch(

      phangorn::RF.dist(
        tr1,
        tr2,
        normalize = FALSE
      ),

      error = function(e) {

        warning(
          "RF.dist failed: ",
          conditionMessage(e)
        )

        NA_real_
      }
    )


  maxrf <-
    2 *
    (
      n - 3
    )


  if (
    maxrf <= 0
  )
    return(NA_real_)


  as.numeric(
    raw
  ) /
    maxrf
}


# =============================================================================
# 9. QUARTET DISTANCE
# =============================================================================

quartet_distance_safe <- function(
  tr1,
  tr2,
  cap
) {

  n <-
    Ntip(
      tr1
    )


  if (
    n < 4 ||
    n > cap
  )
    return(NA_real_)


  result <-
    tryCatch(

      Quartet::QuartetDistance(
        tr1,
        tr2
      ),

      error = function(e) {

        NA_real_
      }
    )


  as.numeric(
    result
  )
}


# =============================================================================
# 10. PAIRWISE DISTANCES
# =============================================================================

all_pairwise_distances <- function(
  pruned,
  cfg
) {

  combos <-
    list(

      c(
        "upstream",
        "cds"
      ),

      c(
        "cds",
        "downstream"
      ),

      c(
        "upstream",
        "downstream"
      )
    )


  purrr::map_dfr(
    combos,

    function(cmb) {

      r1 <-
        cmb[1]

      r2 <-
        cmb[2]


      if (
        is.null(
          pruned[[r1]]
        ) ||
        is.null(
          pruned[[r2]]
        )
      ) {

        return(NULL)
      }


      tibble(

        comparison =
          paste0(
            r1,
            "_vs_",
            r2
          ),

        grf =
          grf_distance(
            pruned[[r1]],
            pruned[[r2]]
          ),

        nrf =
          nrf_distance(
            pruned[[r1]],
            pruned[[r2]]
          ),

        quartet =
          quartet_distance_safe(
            pruned[[r1]],
            pruned[[r2]],
            cfg$quartet_leaf_cap
          )
      )
    }
  )
}


# =============================================================================
# 11. PERMUTATION NULL MODEL
# =============================================================================

permutation_null <- function(
  tr1,
  tr2,
  n_perm,
  metric_fun
) {

  observed <-
    metric_fun(
      tr1,
      tr2
    )


  labels <-
    tr1$tip.label


  null_values <-
    numeric(
      n_perm
    )


  for (
    i in seq_len(
      n_perm
    )
  ) {

    shuffled <-
      tr2


    shuffled$tip.label <-
      sample(
        labels
      )


    null_values[i] <-
      metric_fun(
        tr1,
        shuffled
      )
  }


  null_values <-
    null_values[
      is.finite(
        null_values
      )
    ]


  if (
    length(null_values) < 2
  ) {

    return(
      list(
        observed = observed,
        null_mean = NA_real_,
        null_sd = NA_real_,
        z = NA_real_,
        percentile = NA_real_
      )
    )
  }


  null_mean <-
    mean(
      null_values
    )


  null_sd <-
    sd(
      null_values
    )


  z <-
    if (
      is.finite(null_sd) &&
      null_sd > 0
    ) {

      (
        observed -
          null_mean
      ) /
        null_sd

    } else {

      NA_real_
    }


  percentile <-
    mean(
      null_values <= observed
    )


  list(

    observed =
      observed,

    null_mean =
      null_mean,

    null_sd =
      null_sd,

    z =
      z,

    percentile =
      percentile
  )
}


# =============================================================================
# 12. STRATIFIED TAXON SAMPLING
# =============================================================================

stratified_taxon_sample <- function(
  reference_tree,
  n_target,
  k_clades = NULL
) {

  n <-
    Ntip(
      reference_tree
    )


  if (
    n <= n_target
  ) {

    return(
      reference_tree$tip.label
    )
  }


  if (
    is.null(k_clades)
  ) {

    k_clades <-
      min(
        30,
        max(
          4,
          ceiling(
            n_target / 3
          )
        )
      )
  }


  cop <-
    ape::cophenetic.phylo(
      reference_tree
    )


  hc <-
    hclust(
      as.dist(cop),
      method = "average"
    )


  clades <-
    cutree(
      hc,
      k = k_clades
    )


  df <-
    tibble(

      tip =
        reference_tree$tip.label,

      clade =
        clades[
          reference_tree$tip.label
        ]
    )


  sampled <-
    df %>%
    group_by(
      clade
    ) %>%
    group_modify(
      function(.x, .y) {

        n_take <-
          max(
            1,
            round(
              n_target *
                nrow(.x) /
                n
            )
          )


        n_take <-
          min(
            n_take,
            nrow(.x)
          )


        slice_sample(
          .x,
          n = n_take
        )
      }
    ) %>%
    ungroup() %>%
    pull(tip) %>%
    unique()


  if (
    length(sampled) > n_target
  ) {

    sampled <-
      sample(
        sampled,
        n_target
      )
  }


  sampled
}


# =============================================================================
# 13. MATCHED TAXON BOOTSTRAP
# =============================================================================

matched_bootstrap <- function(
  pruned,
  cfg
) {

  ref_tree <-
    pruned$cds


  if (
    is.null(ref_tree)
  ) {

    warning(
      "CDS tree missing. Bootstrap skipped."
    )


    return(
      tibble()
    )
  }


  purrr::map_dfr(
    seq_len(
      cfg$n_bootstrap
    ),

    function(b) {

      tips <-
        stratified_taxon_sample(
          ref_tree,
          cfg$matched_n_cap
        )


      if (
        length(tips) < 8
      )
        return(NULL)


      sub <-
        purrr::map(
          pruned,

          function(tr) {

            if (
              is.null(tr)
            )
              return(NULL)


            ape::keep.tip(
              tr,
              tips
            )
          }
        )


      d <-
        all_pairwise_distances(
          sub,
          cfg
        )


      d$rep <-
        b


      d
    }
  )
}


# =============================================================================
# 14. CENTRAL HYPOTHESIS
# =============================================================================

test_central_hypothesis <- function(
  boot_df,
  family
) {

  if (
    nrow(boot_df) == 0
  )
    return(NULL)


  wide <-
    boot_df %>%
    select(
      rep,
      comparison,
      grf
    ) %>%
    pivot_wider(
      names_from =
        comparison,

      values_from =
        grf
    )


  required <-
    c(
      "upstream_vs_cds",
      "cds_vs_downstream",
      "upstream_vs_downstream"
    )


  if (
    !all(
      required %in%
        names(wide)
    )
  ) {

    warning(
      "Incomplete pairwise comparisons for ",
      family
    )


    return(NULL)
  }


  wide <-
    wide %>%
    mutate(

      flank_mean =
        (
          upstream_vs_cds +
            cds_vs_downstream
        ) /
        2,

      diff =
        upstream_vs_downstream -
        flank_mean
    )


  valid <-
    wide$diff[
      is.finite(
        wide$diff
      )
    ]


  if (
    length(valid) < 2
  )
    return(NULL)


  wt <-
    suppressWarnings(
      wilcox.test(
        valid,
        alternative = "greater"
      )
    )


  ci <-
    quantile(
      valid,
      probs =
        c(
          0.025,
          0.5,
          0.975
        ),
      na.rm = TRUE
    )


  tibble(

    family =
      family,

    median_diff =
      unname(
        ci[2]
      ),

    ci_lo =
      unname(
        ci[1]
      ),

    ci_hi =
      unname(
        ci[3]
      ),

    wilcoxon_p =
      wt$p.value,

    n_reps =
      length(valid)
  )
}


# =============================================================================
# 15. MAIN ANALYSIS
# =============================================================================

message("")
message("============================================================")
message("TOPOLOGICAL DISCORDANCE ANALYSIS")
message("============================================================")
message("")


manifest <-
  classify_files(
    CONFIG
  )


audit_all <-
  list()


observed_all <-
  list()


null_all <-
  list()


boot_all <-
  list()


pruned_trees_all <-
  list()


for (
  fam in CONFIG$family_levels
) {

  message("")
  message("------------------------------------------------------------")
  message(
    "Processing family: ",
    fam
  )
  message("------------------------------------------------------------")


  trs <-
    load_family_trees(
      manifest,
      fam,
      CONFIG
    )


  available_regions <-
    names(trs)[
      purrr::map_lgl(
        trs,
        function(x) {

          !is.null(x)
        }
      )
    ]


  message(
    "Available regions: ",
    if (
      length(available_regions) > 0
    )
      paste(
        available_regions,
        collapse = ", "
      )
    else
      "NONE"
  )


  # ---------------------------------------------------------------------------
  # Require all three regions
  # ---------------------------------------------------------------------------

  if (
    length(available_regions) < 3
  ) {

    warning(
      fam,
      ": fewer than three regions available. Skipping."
    )


    next
  }


  ap <-
    audit_and_prune(
      trs,
      fam
    )


  audit_all[[fam]] <-
    ap$audit


  message(
    "Common taxa across all three regions: ",
    ap$n_common
  )


  if (
    ap$n_common < 8
  ) {

    warning(
      fam,
      ": fewer than 8 common taxa. Skipping."
    )


    next
  }


  pruned_trees_all[[fam]] <-
    ap$pruned


  # ---------------------------------------------------------------------------
  # OBSERVED
  # ---------------------------------------------------------------------------

  obs <-
    all_pairwise_distances(
      ap$pruned,
      CONFIG
    ) %>%
    mutate(

      family =
        fam,

      n_shared =
        ap$n_common
    )


  observed_all[[fam]] <-
    obs


  message("")
  message(
    "Observed distances for ",
    fam,
    ":"
  )


  print(
    obs
  )


  # ---------------------------------------------------------------------------
  # PERMUTATION NULL
  # ---------------------------------------------------------------------------

  combos <-
    list(

      c(
        "upstream",
        "cds"
      ),

      c(
        "cds",
        "downstream"
      ),

      c(
        "upstream",
        "downstream"
      )
    )


  null_rows <-
    purrr::map_dfr(

      combos,

      function(cmb) {

        r1 <-
          cmb[1]

        r2 <-
          cmb[2]


        if (
          is.null(
            ap$pruned[[r1]]
          ) ||
          is.null(
            ap$pruned[[r2]]
          )
        )
          return(NULL)


        message(
          "Permutation: ",
          r1,
          " vs ",
          r2,
          " (",
          CONFIG$n_null_perm,
          ")"
        )


        res <-
          permutation_null(
            ap$pruned[[r1]],
            ap$pruned[[r2]],
            CONFIG$n_null_perm,
            grf_distance
          )


        tibble(

          family =
            fam,

          comparison =
            paste0(
              r1,
              "_vs_",
              r2
            ),

          observed =
            res$observed,

          null_mean =
            res$null_mean,

          null_sd =
            res$null_sd,

          z =
            res$z,

          percentile =
            res$percentile
        )
      }
    )


  null_all[[fam]] <-
    null_rows


  # ---------------------------------------------------------------------------
  # BOOTSTRAP
  # ---------------------------------------------------------------------------

  message("")
  message(
    "Running ",
    CONFIG$n_bootstrap,
    " matched-taxon bootstrap replicates..."
  )


  boot <-
    matched_bootstrap(
      ap$pruned,
      CONFIG
    )


  if (
    nrow(boot) > 0
  ) {

    boot$family <-
      fam


    boot_all[[fam]] <-
      boot
  }


  message(
    "Finished family: ",
    fam
  )
}


# =============================================================================
# 16. COMBINE RESULTS
# =============================================================================

taxon_audit <-
  bind_rows(
    audit_all
  )


observed_distances <-
  bind_rows(
    observed_all
  )


null_zscores <-
  bind_rows(
    null_all
  )


bootstrap_results <-
  bind_rows(
    boot_all
  )


# =============================================================================
# 17. CENTRAL HYPOTHESIS RESULTS
# =============================================================================

central_hypothesis <-
  purrr::map_dfr(

    CONFIG$family_levels,

    function(fam) {

      if (
        is.null(
          boot_all[[fam]]
        )
      )
        return(NULL)


      test_central_hypothesis(
        boot_all[[fam]],
        fam
      )
    }
  )


if (
  nrow(central_hypothesis) > 0
) {

  central_hypothesis$p_adj_BH <-
    p.adjust(
      central_hypothesis$wilcoxon_p,
      method = "BH"
    )
}


# =============================================================================
# 18. WRITE CSV OUTPUTS
# =============================================================================

write.csv(
  taxon_audit,

  file.path(
    CONFIG$out_dir,
    "taxon_audit.csv"
  ),

  row.names = FALSE
)


write.csv(
  observed_distances,

  file.path(
    CONFIG$out_dir,
    "observed_distances.csv"
  ),

  row.names = FALSE
)


write.csv(
  null_zscores,

  file.path(
    CONFIG$out_dir,
    "standardized_discordance_zscores.csv"
  ),

  row.names = FALSE
)


write.csv(
  bootstrap_results,

  file.path(
    CONFIG$out_dir,
    "matched_bootstrap_replicates.csv"
  ),

  row.names = FALSE
)


write.csv(
  central_hypothesis,

  file.path(
    CONFIG$out_dir,
    "central_hypothesis_test.csv"
  ),

  row.names = FALSE
)


# =============================================================================
# 19. REPRESENTATIVE TANGLEGRAM
# =============================================================================

make_tanglegram <- function(
  pruned,
  family,
  cfg
) {

  if (
    is.null(
      pruned$upstream
    ) ||
    is.null(
      pruned$downstream
    ) ||
    is.null(
      pruned$cds
    )
  ) {

    return(NULL)
  }


  n <-
    Ntip(
      pruned$cds
    )


  if (
    n < 8
  )
    return(NULL)


  n_show <-
    min(
      cfg$tanglegram_max_tips,
      n
    )


  tips <-
    stratified_taxon_sample(
      pruned$cds,
      n_show
    )


  up <-
    ape::keep.tip(
      pruned$upstream,
      tips
    )


  down <-
    ape::keep.tip(
      pruned$downstream,
      tips
    )


  assoc <-
    cbind(
      up$tip.label,
      up$tip.label
    )


  outfile <-
    file.path(
      cfg$out_dir,
      paste0(
        "tanglegram_",
        family,
        "_upstream_vs_downstream.pdf"
      )
    )


  pdf(
    outfile,
    width = 10,
    height = max(
      6,
      n_show * 0.22
    )
  )


  ape::cophyloplot(
    up,
    down,

    assoc = assoc,

    space = 60,

    gap = 3,

    show.tip.label = TRUE,

    col = "grey40",

    lwd = 1,

    main =
      paste0(
        family,
        ": upstream vs downstream"
      )
  )


  dev.off()


  outfile
}


purrr::walk(

  names(
    pruned_trees_all
  ),

  function(fam) {

    make_tanglegram(
      pruned_trees_all[[fam]],
      fam,
      CONFIG
    )
  }
)


# =============================================================================
# 20. MAIN PLOT
# =============================================================================

if (
  nrow(
    bootstrap_results
  ) > 0
) {

  boot_summary <-
    bootstrap_results %>%
    filter(
      is.finite(grf)
    ) %>%
    group_by(
      family,
      comparison
    ) %>%
    summarise(

      median_grf =
        median(
          grf,
          na.rm = TRUE
        ),

      lo =
        quantile(
          grf,
          0.025,
          na.rm = TRUE
        ),

      hi =
        quantile(
          grf,
          0.975,
          na.rm = TRUE
        ),

      n_reps =
        sum(
          is.finite(grf)
        ),

      .groups = "drop"
    )


  boot_summary$family <-
    factor(
      boot_summary$family,
      levels =
        CONFIG$family_levels
    )


  comparison_levels <-
    c(
      "upstream_vs_cds",
      "cds_vs_downstream",
      "upstream_vs_downstream"
    )


  comparison_labels <-
    c(
      "Upstream-CDS",
      "CDS-Downstream",
      "Upstream-Downstream"
    )


  boot_summary$comparison <-
    factor(
      boot_summary$comparison,
      levels =
        comparison_levels,
      labels =
        comparison_labels
    )


  plot_data <-
    bootstrap_results


  plot_data$family <-
    factor(
      plot_data$family,
      levels =
        CONFIG$family_levels
    )


  plot_data$comparison <-
    factor(
      plot_data$comparison,
      levels =
        comparison_levels,
      labels =
        comparison_labels
    )


  p <-
    ggplot(
      boot_summary,

      aes(
        x = comparison,
        y = median_grf,
        colour = family
      )
    ) +

    geom_point(
      position =
        position_dodge(
          width = 0.6
        ),

      size = 2.5
    ) +

    geom_errorbar(
      aes(
        ymin = lo,
        ymax = hi
      ),

      position =
        position_dodge(
          width = 0.6
        ),

      width = 0,

      linewidth = 0.6
    ) +

    geom_jitter(
      data =
        plot_data,

      aes(
        x = comparison,
        y = grf,
        colour = family
      ),

      alpha = 0.06,

      size = 0.4,

      position =
        position_jitterdodge(
          dodge.width = 0.6,
          jitter.width = 0.15
        ),

      inherit.aes = FALSE
    ) +

    scale_colour_manual(
      values =
        CONFIG$palette
    ) +

    labs(

      x = NULL,

      y =
        "Normalized Clustering Information Distance",

      colour = NULL,

      title =
        "Topological discordance among upstream, CDS and downstream",

      subtitle =
        paste0(
          "Matched-taxon bootstrap; n = ",
          CONFIG$matched_n_cap,
          ", B = ",
          CONFIG$n_bootstrap
        )
    ) +

    theme_minimal(
      base_size = 11
    ) +

    theme(

      panel.grid.minor =
        element_blank(),

      panel.grid.major.x =
        element_blank(),

      legend.position =
        "top",

      plot.title =
        element_text(
          face = "bold"
        )
    )


  ggsave(
    file.path(
      CONFIG$out_dir,
      "main_discordance_plot.png"
    ),

    p,

    width = 8,

    height = 5,

    dpi = 400
  )


  ggsave(
    file.path(
      CONFIG$out_dir,
      "main_discordance_plot.pdf"
    ),

    p,

    width = 8,

    height = 5
  )
}


# =============================================================================
# 21. FINAL SUMMARY
# =============================================================================

message("")
message("============================================================")
message("PIPELINE COMPLETE")
message("============================================================")
message("")
message(
  "Results written to:"
)
message(
  CONFIG$out_dir
)
message("")
message("Output files:")
message("  tree_manifest.csv")
message("  taxon_audit.csv")
message("  observed_distances.csv")
message("  standardized_discordance_zscores.csv")
message("  matched_bootstrap_replicates.csv")
message("  central_hypothesis_test.csv")
message("  main_discordance_plot.png")
message("  main_discordance_plot.pdf")
message("")