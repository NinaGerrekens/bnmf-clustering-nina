#!/usr/bin/env Rscript

# Run-to-run stability analysis for a checkpointed bNMF analysis.
#
# Unlike subject-clustering consensus scripts, this analysis treats W rows as
# genetic variants and H columns as traits/features. It provides:
#   1. The post-shrinkage K distribution, including explicit modal ties.
#   2. A label-invariant variant co-assignment consensus matrix.
#   3. H- and W-side component alignment to the final reference solution.
#   4. Per-cluster recovery, split/merge, top-member, and consensus metrics.
#
# Usage:
#   Rscript scripts/assess_bNMF_stability.R RESULTS_DIR [OUTPUT_DIR]


TOP_N <- 20L
RECOVERY_THRESHOLD <- 0.80
STRONG_RECOVERY_THRESHOLD <- 0.90
PAC_LOWER <- 0.10
PAC_UPPER <- 0.90

## ===========================================================================


args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[1L] %in% c("-h", "--help")) {
  cat("Usage: Rscript scripts/assess_bNMF_stability.R RESULTS_DIR [OUTPUT_DIR]\n")
  quit(status = if (length(args)) 0L else 1L)
}
results_dir <- args[1L]
if (!dir.exists(results_dir)) stop("Results directory not found: ", results_dir)
results_dir <- normalizePath(results_dir)

output_dir <- if (length(args) >= 2L) {
  args[2L]
} else {
  file.path(results_dir, "stability_analysis")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)

required_packages <- c("ggplot2", "pheatmap", "clue")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

bnmf_path <- file.path(results_dir, "bnmf_out.rds")
cph_path <- file.path(results_dir, "cluster_post_hoc.rds")
if (!file.exists(bnmf_path)) stop("Missing bnmf_out.rds: ", bnmf_path)

message("Loading repetitions and final reference solution...")
bnmf_reps <- readRDS(bnmf_path)
cph <- if (file.exists(cph_path)) readRDS(cph_path) else NULL
if (!is.list(bnmf_reps) || !length(bnmf_reps)) {
  stop("bnmf_out.rds does not contain a non-empty repetition list")
}


## ---- helpers ---------------------------------------------------------------

last_list_value <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  value <- x[[length(x)]]
  if (!length(value)) NA_real_ else value
}

active_components <- function(rep_result) {
  if (!is.null(rep_result$active_components)) {
    active <- as.integer(rep_result$active_components)
  } else if (!is.null(rep_result$lambda_cut)) {
    lambda <- last_list_value(rep_result$n.lambda)
    if (length(lambda) <= 1L || any(!is.finite(lambda))) {
      stop("A repetition has missing or invalid final lambda values")
    }
    active <- which(lambda >= rep_result$lambda_cut)
  } else {
    active <- which(
      colSums(rep_result$W > 1e-10) > 0 &
        rowSums(rep_result$H > 1e-10) > 0
    )
  }

  active <- active[active >= 1L & active <= ncol(rep_result$W) &
                     active <= nrow(rep_result$H)]
  if (!length(active)) stop("A repetition has no valid active components")
  unique(active)
}

normalize_columns <- function(x) {
  norms <- sqrt(colSums(x^2))
  out <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  keep <- is.finite(norms) & norms > 0
  out[, keep] <- sweep(x[, keep, drop = FALSE], 2L, norms[keep], "/")
  out
}

normalize_rows <- function(x) {
  norms <- sqrt(rowSums(x^2))
  out <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  keep <- is.finite(norms) & norms > 0
  out[keep, ] <- sweep(x[keep, , drop = FALSE], 1L, norms[keep], "/")
  out
}

jaccard <- function(a, b) {
  union_size <- length(union(a, b))
  if (!union_size) return(NA_real_)
  length(intersect(a, b)) / union_size
}

top_names <- function(values, names_vector, n = TOP_N) {
  n <- min(as.integer(n), length(values))
  if (n <= 0L) return(character(0))
  names_vector[head(order(values, decreasing = TRUE), n)]
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x)) {
    unname(stats::quantile(x, probability, names = FALSE, type = 8))
  } else {
    NA_real_
  }
}


## ---- reference W and H -----------------------------------------------------

if (!is.null(cph)) {
  w_columns <- grep("^X[0-9]+$", colnames(cph$w), value = TRUE)
  if (!length(w_columns)) stop("No X<number> cluster columns found in cph$w")
  reference_W <- as.matrix(cph$w[, w_columns, drop = FALSE])
  rownames(reference_W) <- as.character(cph$w$variant)
  reference_H <- as.matrix(cph$h)
} else {
  message("cluster_post_hoc.rds is absent; selecting the minimum-evidence run at the modal K as the reference.")
  repetition_K <- vapply(bnmf_reps, function(x) length(active_components(x)), integer(1L))
  K_table <- table(repetition_K)
  modal_K <- as.integer(names(K_table)[which.max(K_table)])
  eligible_runs <- which(repetition_K == modal_K)
  eligible_evidence <- vapply(
    bnmf_reps[eligible_runs],
    function(x) as.numeric(last_list_value(x$n.evid)),
    numeric(1L)
  )
  reference_run <- eligible_runs[which.min(eligible_evidence)]
  reference_result <- bnmf_reps[[reference_run]]
  reference_active <- active_components(reference_result)
  reference_W <- reference_result$W[, reference_active, drop = FALSE]
  reference_H <- reference_result$H[reference_active, , drop = FALSE]
  colnames(reference_H) <- make.names(colnames(reference_H), unique = TRUE)
}
reference_K <- ncol(reference_W)
if (nrow(reference_H) != reference_K) {
  stop("Reference W/H disagree on K: ", reference_K, " versus ", nrow(reference_H))
}
rownames(reference_H) <- paste0("C", seq_len(reference_K))
colnames(reference_W) <- paste0("C", seq_len(reference_K))

variant_ids <- rownames(reference_W)
feature_ids <- colnames(reference_H)
if (anyNA(variant_ids) || anyDuplicated(variant_ids)) {
  stop("Reference W requires unique, non-missing variant identifiers")
}
if (is.null(feature_ids) || anyNA(feature_ids) || anyDuplicated(feature_ids)) {
  stop("Reference H requires unique, non-missing feature identifiers")
}

for (i in seq_along(bnmf_reps)) {
  rep_result <- bnmf_reps[[i]]
  if (!all(variant_ids %in% rownames(rep_result$W))) {
    stop("Repetition ", i, " is missing reference variants")
  }
  # cluster_post_hoc is reconstructed from a text table, so R sanitizes spaces,
  # hyphens, and the :: duplicate-name separator with make.names().
  rep_feature_ids <- make.names(colnames(rep_result$H), unique = TRUE)
  if (!all(feature_ids %in% rep_feature_ids)) {
    stop("Repetition ", i, " is missing reference features")
  }
}


## ---- 1. K and run-level diagnostics ---------------------------------------

message("Calculating K distribution and run diagnostics...")
run_diagnostics <- do.call(rbind, lapply(seq_along(bnmf_reps), function(i) {
  rep_result <- bnmf_reps[[i]]
  active <- active_components(rep_result)
  final_error <- last_list_value(rep_result$n.error)
  final_evidence <- last_list_value(rep_result$n.evid)
  data.frame(
    run = i,
    K = length(active),
    final_error = as.numeric(final_error),
    final_evidence = as.numeric(final_evidence),
    iterations = if (!is.null(rep_result$iterations)) rep_result$iterations else NA_integer_,
    elapsed_minutes = if (!is.null(rep_result$elapsed_seconds)) {
      rep_result$elapsed_seconds / 60
    } else {
      NA_real_
    },
    seed = if (!is.null(rep_result$seed)) rep_result$seed else NA_integer_,
    stringsAsFactors = FALSE
  )
}))
write.csv(run_diagnostics, file.path(output_dir, "stability_run_diagnostics.csv"),
          row.names = FALSE)

k_distribution <- as.data.frame(table(run_diagnostics$K), stringsAsFactors = FALSE)
colnames(k_distribution) <- c("K", "n_reps")
k_distribution$K <- as.integer(as.character(k_distribution$K))
k_distribution$fraction <- k_distribution$n_reps / sum(k_distribution$n_reps)
max_count <- max(k_distribution$n_reps)
modal_K_values <- k_distribution$K[k_distribution$n_reps == max_count]
k_distribution$is_modal <- k_distribution$K %in% modal_K_values
k_distribution$is_reference_K <- k_distribution$K == reference_K
write.csv(k_distribution, file.path(output_dir, "stability_k_distribution.csv"),
          row.names = FALSE)

k_plot <- ggplot2::ggplot(
  k_distribution,
  ggplot2::aes(x = factor(K), y = n_reps, fill = is_modal)
) +
  ggplot2::geom_col(alpha = 0.9) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.0f%%", 100 * fraction)),
    vjust = -0.3, size = 3.5
  ) +
  ggplot2::scale_fill_manual(
    values = c("TRUE" = "coral", "FALSE" = "steelblue"), guide = "none"
  ) +
  ggplot2::labs(
    title = sprintf("Post-shrinkage K across %d bNMF repetitions", length(bnmf_reps)),
    subtitle = sprintf(
      "Modal K: %s; final reference K: %d",
      paste(modal_K_values, collapse = " and "), reference_K
    ),
    x = "Number of active components", y = "Repetitions"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::expand_limits(y = max(k_distribution$n_reps) * 1.15)
ggplot2::ggsave(
  file.path(output_dir, "stability_k_distribution.png"), k_plot,
  width = 8, height = 5, dpi = 200
)


## ---- 2. Label-invariant variant consensus ---------------------------------

message("Building full variant co-assignment consensus matrix...")
rep_variant_labels <- lapply(bnmf_reps, function(rep_result) {
  active <- active_components(rep_result)
  rep_W <- rep_result$W[match(variant_ids, rownames(rep_result$W)), active,
                        drop = FALSE]
  max.col(rep_W, ties.method = "first")
})

consensus_counts <- matrix(0L, length(variant_ids), length(variant_ids))
for (labels in rep_variant_labels) {
  consensus_counts <- consensus_counts + outer(labels, labels, "==")
}
variant_consensus <- consensus_counts / length(rep_variant_labels)
rownames(variant_consensus) <- colnames(variant_consensus) <- variant_ids
saveRDS(
  variant_consensus,
  file.path(output_dir, "stability_variant_consensus_matrix.rds"),
  compress = FALSE
)

reference_variant_cluster <- max.col(reference_W, ties.method = "first")
reference_variant_max <- reference_W[
  cbind(seq_len(nrow(reference_W)), reference_variant_cluster)
]
reference_W_rowsum <- rowSums(reference_W)
reference_variant_fraction <- ifelse(
  reference_W_rowsum > 0, reference_variant_max / reference_W_rowsum, NA_real_
)
reference_variant_second <- apply(reference_W, 1L, function(x) {
  sorted <- sort(x, decreasing = TRUE)
  if (length(sorted) >= 2L) sorted[2L] else NA_real_
})
variant_assignments <- data.frame(
  variant = variant_ids,
  reference_cluster = reference_variant_cluster,
  max_W = reference_variant_max,
  second_W = reference_variant_second,
  max_to_second_ratio = reference_variant_max / pmax(reference_variant_second, 1e-300),
  max_fraction_of_row_W = reference_variant_fraction,
  stringsAsFactors = FALSE
)
write.csv(
  variant_assignments,
  file.path(output_dir, "stability_reference_variant_assignments.csv"),
  row.names = FALSE
)

upper <- upper.tri(variant_consensus)
same_reference_cluster <- outer(
  reference_variant_cluster, reference_variant_cluster, "=="
)
ambiguous_pair <- variant_consensus > PAC_LOWER & variant_consensus < PAC_UPPER

overall_PAC <- mean(ambiguous_pair[upper])
within_PAC <- mean(ambiguous_pair[upper & same_reference_cluster])
between_PAC <- mean(ambiguous_pair[upper & !same_reference_cluster])

consensus_by_cluster <- do.call(rbind, lapply(seq_len(reference_K), function(k) {
  members <- which(reference_variant_cluster == k)
  nonmembers <- which(reference_variant_cluster != k)
  within_values <- if (length(members) >= 2L) {
    block <- variant_consensus[members, members, drop = FALSE]
    block[upper.tri(block)]
  } else {
    numeric(0)
  }
  between_values <- if (length(members) && length(nonmembers)) {
    as.numeric(variant_consensus[members, nonmembers, drop = FALSE])
  } else {
    numeric(0)
  }
  data.frame(
    cluster = k,
    n_hard_assigned_variants = length(members),
    mean_within_consensus = safe_mean(within_values),
    median_within_consensus = safe_median(within_values),
    p05_within_consensus = safe_quantile(within_values, 0.05),
    mean_between_consensus = safe_mean(between_values),
    consensus_separation = safe_mean(within_values) - safe_mean(between_values),
    within_PAC = if (length(within_values)) {
      mean(within_values > PAC_LOWER & within_values < PAC_UPPER)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}))
write.csv(
  consensus_by_cluster,
  file.path(output_dir, "stability_variant_consensus_by_cluster.csv"),
  row.names = FALSE
)

consensus_order <- order(
  reference_variant_cluster,
  -reference_variant_fraction,
  -reference_variant_max
)
ordered_consensus <- variant_consensus[consensus_order, consensus_order]
plot_ids <- paste0("v", seq_along(consensus_order))
rownames(ordered_consensus) <- colnames(ordered_consensus) <- plot_ids
annotation <- data.frame(
  cluster = factor(reference_variant_cluster[consensus_order]),
  row.names = plot_ids
)
annotation_levels <- levels(annotation$cluster)
annotation_colors <- list(cluster = stats::setNames(
  grDevices::hcl.colors(length(annotation_levels), palette = "Dynamic"),
  annotation_levels
))
cluster_block_sizes <- as.integer(table(annotation$cluster))
cluster_gaps <- cumsum(cluster_block_sizes)
cluster_gaps <- cluster_gaps[-length(cluster_gaps)]
grDevices::png(
  file.path(output_dir, "stability_variant_consensus_heatmap.png"),
  width = 9, height = 8, units = "in", res = 200
)
pheatmap::pheatmap(
  ordered_consensus,
  color = grDevices::colorRampPalette(c("white", "steelblue", "navy"))(100),
  breaks = seq(0, 1, length.out = 101),
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  annotation_row = annotation, annotation_col = annotation,
  annotation_colors = annotation_colors,
  annotation_legend = FALSE,
  gaps_row = cluster_gaps, gaps_col = cluster_gaps,
  main = sprintf(
    "Variant co-assignment consensus across %d repetitions\nordered by final K=%d solution",
    length(bnmf_reps), reference_K
  ),
  border_color = NA
)
grDevices::dev.off()


## ---- 3. Align every repetition to final reference components --------------

message("Aligning repetition components to final W/H reference...")
reference_W_normalized <- normalize_columns(reference_W)
reference_H_normalized <- normalize_rows(reference_H)

reference_top_variants <- lapply(seq_len(reference_K), function(k) {
  top_names(reference_W[, k], variant_ids, TOP_N)
})
reference_top_features <- lapply(seq_len(reference_K), function(k) {
  top_names(reference_H[k, ], feature_ids, TOP_N)
})

match_rows <- vector("list", length(bnmf_reps))
best_rows <- vector("list", length(bnmf_reps))

for (run_index in seq_along(bnmf_reps)) {
  rep_result <- bnmf_reps[[run_index]]
  active <- active_components(rep_result)
  rep_K <- length(active)
  rep_W <- rep_result$W[
    match(variant_ids, rownames(rep_result$W)), active, drop = FALSE
  ]
  rep_H <- rep_result$H[
    active,
    match(feature_ids, make.names(colnames(rep_result$H), unique = TRUE)),
    drop = FALSE
  ]

  W_similarity <- crossprod(reference_W_normalized, normalize_columns(rep_W))
  H_similarity <- reference_H_normalized %*% t(normalize_rows(rep_H))
  joint_similarity <- (W_similarity + H_similarity) / 2

  # Pad to a square matrix so Hungarian matching works for K above or below the
  # reference K. Assignments to padded columns are recorded as unmatched.
  padded_K <- max(reference_K, rep_K)
  padded_similarity <- matrix(0, padded_K, padded_K)
  padded_similarity[seq_len(reference_K), seq_len(rep_K)] <- joint_similarity
  hungarian_assignment <- as.integer(clue::solve_LSAP(
    padded_similarity, maximum = TRUE
  ))

  unrestricted_best_rep <- max.col(joint_similarity, ties.method = "first")
  unrestricted_best_joint <- joint_similarity[
    cbind(seq_len(reference_K), unrestricted_best_rep)
  ]
  rep_best_reference <- max.col(t(joint_similarity), ties.method = "first")
  rep_best_joint <- joint_similarity[cbind(rep_best_reference, seq_len(rep_K))]

  match_rows[[run_index]] <- do.call(rbind, lapply(seq_len(reference_K), function(k) {
    matched_component <- hungarian_assignment[k]
    if (matched_component > rep_K) matched_component <- NA_integer_

    if (is.na(matched_component)) {
      h_cos <- w_cos <- joint_cos <- top_variant_jaccard <-
        top_feature_jaccard <- NA_real_
    } else {
      h_cos <- H_similarity[k, matched_component]
      w_cos <- W_similarity[k, matched_component]
      joint_cos <- joint_similarity[k, matched_component]
      top_variant_jaccard <- jaccard(
        reference_top_variants[[k]],
        top_names(rep_W[, matched_component], variant_ids, TOP_N)
      )
      top_feature_jaccard <- jaccard(
        reference_top_features[[k]],
        top_names(rep_H[matched_component, ], feature_ids, TOP_N)
      )
    }

    best_component <- unrestricted_best_rep[k]
    merge_multiplicity <- sum(
      unrestricted_best_rep == best_component &
        unrestricted_best_joint >= RECOVERY_THRESHOLD
    )
    split_count <- sum(
      rep_best_reference == k & rep_best_joint >= RECOVERY_THRESHOLD
    )

    data.frame(
      run = run_index,
      run_K = rep_K,
      reference_cluster = k,
      matched_component = matched_component,
      H_cosine = h_cos,
      W_cosine = w_cos,
      joint_cosine = joint_cos,
      top_variant_jaccard = top_variant_jaccard,
      top_feature_jaccard = top_feature_jaccard,
      unrestricted_best_component = best_component,
      unrestricted_best_joint_cosine = unrestricted_best_joint[k],
      split_component_count_ge_threshold = split_count,
      merge_reference_count_ge_threshold = merge_multiplicity,
      recovered_ge_0_80 = is.finite(joint_cos) &&
        joint_cos >= RECOVERY_THRESHOLD,
      recovered_ge_0_90 = is.finite(joint_cos) &&
        joint_cos >= STRONG_RECOVERY_THRESHOLD,
      stringsAsFactors = FALSE
    )
  }))

  best_rows[[run_index]] <- do.call(rbind, lapply(seq_len(rep_K), function(j) {
    best_reference <- rep_best_reference[j]
    data.frame(
      run = run_index,
      run_K = rep_K,
      repetition_component = j,
      best_reference_cluster = best_reference,
      H_cosine = H_similarity[best_reference, j],
      W_cosine = W_similarity[best_reference, j],
      joint_cosine = joint_similarity[best_reference, j],
      stringsAsFactors = FALSE
    )
  }))
}

component_matches <- do.call(rbind, match_rows)
repetition_component_best_matches <- do.call(rbind, best_rows)
write.csv(
  component_matches,
  file.path(output_dir, "stability_component_matches_long.csv"),
  row.names = FALSE
)
write.csv(
  repetition_component_best_matches,
  file.path(output_dir, "stability_repetition_components_best_match.csv"),
  row.names = FALSE
)


## ---- 4. Per-reference-cluster component stability -------------------------

message("Summarizing per-cluster recovery and split/merge behavior...")
component_stability <- do.call(rbind, lapply(seq_len(reference_K), function(k) {
  d <- component_matches[component_matches$reference_cluster == k, , drop = FALSE]
  data.frame(
    cluster = k,
    n_repetitions = nrow(d),
    recovery_rate_joint_ge_0_80 = mean(d$recovered_ge_0_80),
    recovery_rate_joint_ge_0_90 = mean(d$recovered_ge_0_90),
    mean_H_cosine = safe_mean(d$H_cosine),
    p05_H_cosine = safe_quantile(d$H_cosine, 0.05),
    median_H_cosine = safe_median(d$H_cosine),
    mean_W_cosine = safe_mean(d$W_cosine),
    p05_W_cosine = safe_quantile(d$W_cosine, 0.05),
    median_W_cosine = safe_median(d$W_cosine),
    mean_joint_cosine = safe_mean(d$joint_cosine),
    p05_joint_cosine = safe_quantile(d$joint_cosine, 0.05),
    median_joint_cosine = safe_median(d$joint_cosine),
    mean_top_variant_jaccard = safe_mean(d$top_variant_jaccard),
    mean_top_feature_jaccard = safe_mean(d$top_feature_jaccard),
    fraction_runs_split_ge_2 = mean(
      d$split_component_count_ge_threshold >= 2, na.rm = TRUE
    ),
    fraction_runs_merged_ge_2 = mean(
      d$merge_reference_count_ge_threshold >= 2, na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}))

component_stability <- merge(
  component_stability, consensus_by_cluster,
  by = "cluster", all.x = TRUE, sort = TRUE
)
component_stability$stability_flag <- ifelse(
  component_stability$recovery_rate_joint_ge_0_80 >= 0.80 &
    component_stability$mean_joint_cosine >= 0.80,
  "stable",
  ifelse(
    component_stability$recovery_rate_joint_ge_0_80 >= 0.50,
    "intermediate",
    "fragile"
  )
)
write.csv(
  component_stability,
  file.path(output_dir, "stability_per_cluster.csv"),
  row.names = FALSE
)

recovery_plot <- ggplot2::ggplot(
  component_stability,
  ggplot2::aes(x = factor(cluster), y = mean_joint_cosine,
               fill = stability_flag)
) +
  ggplot2::geom_col(na.rm = TRUE) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = p05_joint_cosine, ymax = mean_joint_cosine),
    width = 0.25
  ) +
  ggplot2::geom_hline(
    yintercept = RECOVERY_THRESHOLD, linetype = 2, color = "firebrick"
  ) +
  ggplot2::scale_fill_manual(values = c(
    stable = "#2E8B57", intermediate = "#E6A23C", fragile = "#C94C4C"
  )) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Final-cluster recovery across bNMF repetitions",
    subtitle = "Bars: mean joint H/W cosine; error bars: 5th percentile to mean",
    x = "Final reference cluster", y = "Joint cosine similarity", fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
ggplot2::ggsave(
  file.path(output_dir, "stability_component_recovery.png"), recovery_plot,
  width = 11, height = 6, dpi = 200
)

consensus_plot <- ggplot2::ggplot(
  component_stability,
  ggplot2::aes(x = factor(cluster), y = mean_within_consensus,
               fill = stability_flag)
) +
  ggplot2::geom_col(na.rm = TRUE) +
  ggplot2::geom_point(
    ggplot2::aes(y = mean_between_consensus), color = "black", size = 1.8
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::scale_fill_manual(values = c(
    stable = "#2E8B57", intermediate = "#E6A23C", fragile = "#C94C4C"
  )) +
  ggplot2::labs(
    title = "Variant consensus by final cluster",
    subtitle = "Bars: mean within-cluster consensus; black points: mean between-cluster consensus",
    x = "Final reference cluster", y = "Co-assignment consensus", fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
ggplot2::ggsave(
  file.path(output_dir, "stability_variant_consensus_by_cluster.png"),
  consensus_plot, width = 11, height = 6, dpi = 200
)


## ---- human-readable summary ------------------------------------------------

modal_description <- if (length(modal_K_values) == 1L) {
  sprintf("Unique modal K = %d (%d/%d repetitions).",
          modal_K_values, max_count, length(bnmf_reps))
} else {
  sprintf(
    "Modal tie: K = %s each occurred in %d/%d repetitions; reference K = %d.",
    paste(modal_K_values, collapse = " and "), max_count,
    length(bnmf_reps), reference_K
  )
}

flag_counts <- table(component_stability$stability_flag)
flag_description <- paste(
  paste(names(flag_counts), as.integer(flag_counts), sep = "="),
  collapse = ", "
)

summary_lines <- c(
  "bNMF RUN-TO-RUN STABILITY SUMMARY",
  "=================================",
  sprintf("Results directory: %s", results_dir),
  sprintf("Repetitions: %d", length(bnmf_reps)),
  sprintf("Final reference K: %d", reference_K),
  modal_description,
  sprintf(
    "K range: %d-%d; K within reference +/-2: %.1f%%",
    min(run_diagnostics$K), max(run_diagnostics$K),
    100 * mean(abs(run_diagnostics$K - reference_K) <= 2)
  ),
  sprintf("Overall variant PAC (%.2f-%.2f): %.3f",
          PAC_LOWER, PAC_UPPER, overall_PAC),
  sprintf("Within-final-cluster PAC: %.3f", within_PAC),
  sprintf("Between-final-cluster PAC: %.3f", between_PAC),
  sprintf("Per-cluster flags: %s", flag_description),
  sprintf("Recovery threshold: joint mean of H/W cosine >= %.2f", RECOVERY_THRESHOLD),
  sprintf("Top-member Jaccard uses top %d variants/features.", TOP_N),
  "",
  "Notes:",
  "- Variant consensus is label-invariant and uses hard max-W assignments.",
  "- Component recovery retains the soft W/H loadings and aligns labels with",
  "  maximum-weight Hungarian matching.",
  "- A K tie is reported explicitly; no majority is claimed when none exists.",
  "- Feature stability is represented by H cosine and top-feature Jaccard rather",
  "  than a 4,872 x 4,872 hard feature-consensus matrix.",
  "",
  sprintf("Outputs: %s", output_dir)
)
writeLines(summary_lines, file.path(output_dir, "stability_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
