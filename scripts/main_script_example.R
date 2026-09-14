#!/usr/bin/env Rscript

# =============================================================================
# Generic bNMF clustering example
# =============================================================================
#
# This script runs the bNMF pipeline using the toy summary statistics in
# example_data/. It is deliberately limited to variant selection, matrix
# preparation, bNMF, and summarize_bNMF(); project-specific post-processing
# belongs in a separate analysis script.
#
# By default, the example uses position-based clumping and keeps the resulting
# original variants. LDlink pruning and proxy replacement are optional because
# they require an LDlink token and, for proxies, a current per-chromosome rsID
# map. The supplied toy data can therefore run without either external input.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(furrr)
  library(magrittr)
  library(readr)
  library(readxl)
  library(softImpute)
  library(tidyr)
})

# =============================================================================
# 0. Configuration
# =============================================================================

# Find the repository from this script's location. This lets the example work
# whether it is run from the repository root or as:
#   Rscript scripts/main_script_example.R
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_file <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
} else if (!is.null(sys.frames()[[1L]]$ofile)) {
  script_file <- normalizePath(sys.frames()[[1L]]$ofile, mustWork = TRUE)
} else {
  script_file <- normalizePath(file.path("scripts", "main_script_example.R"),
                               mustWork = TRUE)
}

repo_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
example_dir <- file.path(repo_dir, "example_data")
scripts_dir <- file.path(repo_dir, "scripts")

version <- "toy_bnmf_example"
main_dir <- file.path(repo_dir, paste0(version, "_results"))
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)

# Variant-selection settings.
PVCUTOFF <- 5e-8
PVCUTOFF_PROXY <- 5e-6
PROXY_WINDOW_KB <- 500
CLUMP_WINDOW_BP <- 100e3

# External LD operations are disabled for the self-contained toy run.
RUN_LDLINK_PRUNING <- FALSE
RUN_PROXY_SEARCH <- FALSE
LD_POPS <- c("EUR", "EAS", "AFR", "AMR", "SAS")
LD_R2 <- 0.05
LD_MAF <- 0.001
LDLINK_TOKEN <- Sys.getenv("LDLINK_TOKEN", unset = "")

# Proxy search additionally requires files named chrN.txt with four columns:
# hg19 position (chrN:position), rsID, reference allele, alternate allele.
RSID_MAP_DIR <- file.path(example_dir, "rsid_maps_by_chr")
PROXY_MIN_NONMISSING <- 0.8
PROXY_MIN_R2 <- 0.8

# Trait and bNMF settings.
MIN_MEDIAN_N <- 5000
MAX_TRAIT_MISSING <- 0.30
CORRELATION_CUTOFF <- 0.8
IMPUTATION_HOLDOUT <- 0.10
RANDOM_SEED <- 123

BNMF_REPS <- 10
BNMF_K_INITIAL <- 15
BNMF_K0 <- 10
BNMF_TOLERANCE <- 1e-6
BNMF_PHI <- 1
BNMF_WORKERS <- max(1L, min(2L, parallelly::availableCores()))

if (RUN_PROXY_SEARCH && !RUN_LDLINK_PRUNING) {
  stop("RUN_PROXY_SEARCH requires RUN_LDLINK_PRUNING so original rsIDs are available.")
}
if ((RUN_LDLINK_PRUNING || RUN_PROXY_SEARCH) && !nzchar(LDLINK_TOKEN)) {
  stop("Set the LDLINK_TOKEN environment variable before enabling LDlink steps.")
}
if (RUN_PROXY_SEARCH && !dir.exists(RSID_MAP_DIR)) {
  stop("Proxy search requires a current per-chromosome rsID map directory: ",
       RSID_MAP_DIR)
}

# =============================================================================
# 1. Load helper functions and toy manifest
# =============================================================================

source(file.path(scripts_dir, "choose_variants_2025.R"))
source(file.path(scripts_dir, "prep_bNMF_2025.R"))
source(file.path(scripts_dir, "run_bNMF_2025.R"))
future::plan(future::sequential)

manifest_file <- file.path(example_dir, "clustering_data_source_example.xlsx")
gwas <- read_excel(manifest_file, sheet = "main_gwas") %>%
  mutate(
    ID = paste(study, trait, population, sep = "_"),
    # Resolve the example paths from the repository rather than trusting the
    # working directory encoded in the workbook.
    full_path = file.path(example_dir, "my_GWAS", basename(file))
  )

gwas_traits <- read_excel(manifest_file, sheet = "trait_gwas") %>%
  mutate(full_path = file.path(example_dir, "my_GWAS", basename(file)))

input_files <- c(gwas$full_path, gwas_traits$full_path)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop("Missing toy summary-statistics file(s): ",
       paste(missing_files, collapse = ", "))
}

main_rows <- which(toupper(gwas$largest) == "YES")
if (length(main_rows) != 1L) {
  stop("The main_gwas sheet must mark exactly one row as largest = 'Yes'.")
}
main_ss_filepath <- gwas$full_path[main_rows]

trait_ss_files <- setNames(gwas_traits$full_path, gwas_traits$trait)
trait_ss_size <- setNames(as.numeric(gwas_traits$sample_size), gwas_traits$trait)

# =============================================================================
# 2. Select and position-clump sentinel variants
# =============================================================================

# The broad set supplies possible proxy candidates. Sentinel eligibility is
# still determined using each discovery GWAS's p-value and PVCUTOFF below.
vars_sig <- get_sig_snps(
  gwas = gwas,
  rename_cols = NULL,
  PVCUTOFF = PVCUTOFF_PROXY
) %>%
  mutate(PVALUE = if_else(PVALUE == 0, 1e-300, PVALUE))

# Normalize the designated primary GWAS in memory. This works for both the
# uncompressed toy file and real data supplied as a data frame to
# fetch_summary_stats().
gwas_primary <- fread(main_ss_filepath, data.table = FALSE)
if (!"BETA" %in% names(gwas_primary)) {
  gwas_primary <- gwas_primary %>%
    mutate(BETA = log(as.numeric(ODDS_RATIO)))
}
if (!"SE" %in% names(gwas_primary)) {
  gwas_primary <- gwas_primary %>%
    mutate(SE = abs(BETA / qnorm(pmax(as.numeric(P_VALUE),
                                     .Machine$double.xmin) / 2)))
}
gwas_primary <- gwas_primary %>%
  mutate(P_VALUE = as.numeric(P_VALUE), BETA = as.numeric(BETA), SE = as.numeric(SE)) %>%
  separate(VAR_ID, into = c("CHR", "POS", "REF", "ALT"),
           sep = "_", remove = FALSE, convert = TRUE) %>%
  mutate(SNP = paste(CHR, POS, sep = ":"))

primary_pvalues <- gwas_primary %>%
  transmute(VAR_ID, PVALUE = P_VALUE)

# Mirror get_biggest_gwas(): retain the strongest discovery record for each
# variant, then require the variant to be represented in the primary GWAS.
vars_main <- vars_sig %>%
  arrange(PVALUE) %>%
  distinct(VAR_ID, .keep_all = TRUE) %>%
  rename(PVALUE.Pop = PVALUE) %>%
  inner_join(primary_pvalues, by = "VAR_ID") %>%
  separate(VAR_ID, into = c("CHR", "POS", "REF.primary", "ALT.primary"),
           sep = "_", remove = FALSE, convert = TRUE) %>%
  mutate(
    REF = REF.primary,
    ALT = ALT.primary,
    ChrPos = paste(CHR, POS, sep = ":"),
    ChrPos_LDlink = paste0("chr", ChrPos)
  ) %>%
  select(-REF.primary, -ALT.primary)

vars_no_hla <- vars_main %>%
  filter(!(CHR == 6 & between(POS, 28477797, 33448354)))

sentinel_candidates <- vars_no_hla %>%
  filter(PVALUE.Pop < PVCUTOFF)
if (nrow(sentinel_candidates) == 0L) {
  stop("No sentinel variants passed PVCUTOFF.")
}

clumped_ids <- snp_clump(
  sentinel_candidates,
  id = "VAR_ID",
  window = CLUMP_WINDOW_BP
)
vars_clumped <- sentinel_candidates %>%
  filter(VAR_ID %in% clumped_ids)

message(sprintf(
  "Variant selection: %d broad candidates, %d sentinels, %d after position clumping.",
  nrow(vars_no_hla), nrow(sentinel_candidates), nrow(vars_clumped)
))

# =============================================================================
# 3. Optional LDlink pruning
# =============================================================================

if (RUN_LDLINK_PRUNING) {
  ld_result_dir <- file.path(main_dir, "ld_pruning")
  dir.create(ld_result_dir, recursive = TRUE, showWarnings = FALSE)

  for (ld_pop in LD_POPS) {
    ld_pruning_SNP.clip(
      df_snps = vars_clumped,
      pop = ld_pop,
      output_dir = ld_result_dir,
      r2 = LD_R2,
      maf = LD_MAF,
      chr = 1:22,
      token = LDLINK_TOKEN
    )
  }

  # Read the explicitly named outputs. This also makes a missing or failed
  # population a hard error instead of silently treating its variants as absent.
  ld_by_pop <- lapply(LD_POPS, function(ld_pop) {
    files <- list.files(
      ld_result_dir,
      pattern = paste0("^snpClip_results_", ld_pop, "_chr[0-9]+\\.txt$"),
      full.names = TRUE
    )
    if (length(files) == 0L) {
      stop("No successful SNPclip outputs found for population ", ld_pop, ".")
    }
    bind_rows(lapply(files, fread)) %>%
      rename(RS_Number = any_of("RS Number")) %>%
      mutate(LD_population = ld_pop)
  })

  ld_kept <- bind_rows(ld_by_pop) %>%
    filter(Details == "Variant kept.")

  kept_in_every_panel <- ld_kept %>%
    distinct(Position, RS_Number, LD_population) %>%
    count(Position, RS_Number, name = "n_pop") %>%
    filter(n_pop == length(LD_POPS))

  pruned_vars <- vars_clumped %>%
    inner_join(kept_in_every_panel,
               by = c("ChrPos_LDlink" = "Position")) %>%
    select(-n_pop)
} else {
  message("LDlink pruning is disabled; using position-clumped toy sentinels.")
  pruned_vars <- vars_clumped %>%
    mutate(RS_Number = NA_character_)
}

if (nrow(pruned_vars) == 0L) {
  stop("No variants remain after pruning.")
}

# =============================================================================
# 4. Fetch and harmonize the multi-trait summary statistics
# =============================================================================

fetch_input <- window_to_sentinels(
  candidates = vars_no_hla,
  sentinels = pruned_vars,
  window_kb = PROXY_WINDOW_KB
) %>%
  mutate(SNP = ChrPos) %>%
  arrange(PVALUE) %>%
  distinct(SNP, .keep_all = TRUE)

pval_bonf_sentinel <- 0.05 / nrow(pruned_vars)
trait_checkpoint_dir <- file.path(main_dir, "trait_checkpoints")

z_n_mats <- fetch_summary_stats(
  df_input = fetch_input,
  gwas_ss_file = gwas_primary,
  trait_ss_files = trait_ss_files,
  trait_ss_size = trait_ss_size,
  pval_cutoff = 0.05,
  pval_bonf = pval_bonf_sentinel,
  checkpoint_dir = trait_checkpoint_dir
)

zmat_fullset <- as.matrix(z_n_mats$df_z)
Nmat_fullset <- as.matrix(z_n_mats$df_N)

missing_sentinels <- setdiff(pruned_vars$ChrPos, rownames(zmat_fullset))
if (length(missing_sentinels) > 0L) {
  stop("No z-matrix row was created for sentinel(s): ",
       paste(missing_sentinels, collapse = ", "))
}

zmat_pruned <- zmat_fullset[pruned_vars$ChrPos, , drop = FALSE]
Nmat_pruned <- Nmat_fullset[pruned_vars$ChrPos, , drop = FALSE]

median_N <- apply(Nmat_pruned, 2, median, na.rm = TRUE)
traits_high_N <- names(median_N)[is.finite(median_N) & median_N >= MIN_MEDIAN_N]
trait_missingness <- colMeans(is.na(zmat_pruned[, traits_high_N, drop = FALSE]))
traits_final <- names(trait_missingness)[trait_missingness <= MAX_TRAIT_MISSING]

if (length(traits_final) < 2L) {
  stop("Fewer than two traits remain after sample-size and missingness filtering.")
}

traits_removed_low_N <- setdiff(colnames(zmat_pruned), traits_high_N)
traits_removed_missing <- setdiff(traits_high_N, traits_final)

variant_nonmissingness <- rowMeans(!is.na(zmat_pruned[, traits_final, drop = FALSE]))
variant_nonmissingness <- setNames(variant_nonmissingness, pruned_vars$VAR_ID)

# =============================================================================
# 5. Optional proxy search and final variant assembly
# =============================================================================

if (RUN_PROXY_SEARCH) {
  proxies_needed <- find_variants_needing_proxies(
    gwas_variant_df = pruned_vars,
    var_nonmissingness = variant_nonmissingness
  ) %>%
    inner_join(pruned_vars %>% select(VAR_ID, Population), by = "VAR_ID") %>%
    mutate(
      search_population = case_when(
        Population %in% c("EUR", "TA", "MA") ~ "EUR",
        Population == "SA" ~ "SAS",
        TRUE ~ Population
      )
    )

  pruned_by_search_pop <- pruned_vars %>%
    mutate(
      search_population = case_when(
        Population %in% c("EUR", "TA", "MA") ~ "EUR",
        Population == "SA" ~ "SAS",
        TRUE ~ Population
      )
    )

  supported_pops <- c("EUR", "EAS", "AFR", "AMR", "SAS")
  proxy_pops <- unique(pruned_by_search_pop$search_population)
  unsupported_pops <- setdiff(proxy_pops, supported_pops)
  if (length(unsupported_pops) > 0L) {
    stop("Unsupported LDlink proxy population(s): ",
         paste(unsupported_pops, collapse = ", "))
  }

  proxy_sets <- lapply(proxy_pops, function(pop) {
    needed_pop <- proxies_needed %>%
      filter(search_population == pop) %>%
      select(-search_population)
    pruned_pop <- pruned_by_search_pop %>%
      filter(search_population == pop) %>%
      select(-search_population)

    if (nrow(needed_pop) == 0L) {
      return(list(pruned_pop$VAR_ID, NULL))
    }

    choose_proxies(
      need_proxies = needed_pop,
      rsid_map_dir = RSID_MAP_DIR,
      pruned_variants = pruned_pop,
      zmat_fullset = zmat_fullset[, traits_final, drop = FALSE],
      token = LDLINK_TOKEN,
      population = pop,
      frac_nonmissing_num = PROXY_MIN_NONMISSING,
      r2_num = PROXY_MIN_R2,
      variant_metadata = fetch_input,
      output_dir = file.path(main_dir, "proxy_search", pop)
    )
  })
  names(proxy_sets) <- proxy_pops

  original_ids <- unique(unlist(lapply(proxy_sets, `[[`, 1L)))
  proxy_rows <- bind_rows(lapply(proxy_sets, `[[`, 2L))

  final_originals <- pruned_vars %>%
    filter(VAR_ID %in% original_ids) %>%
    transmute(
      VAR_ID, ChrPos, rsID = RS_Number, Variant_Type = "original",
      original_SNP = VAR_ID, REF, ALT, PVALUE, Risk_Allele, GWAS, Population
    )

  if (nrow(proxy_rows) > 0L) {
    proxy_metadata <- fetch_input %>%
      select(VAR_ID, SNP, REF, ALT, PVALUE, Risk_Allele, GWAS, Population) %>%
      distinct(SNP, .keep_all = TRUE)

    final_proxies <- proxy_rows %>%
      transmute(
        original_SNP = VAR_ID,
        ChrPos = proxy_ChrPos,
        rsID = proxy_rsID,
        Variant_Type = "proxy"
      ) %>%
      inner_join(proxy_metadata, by = c("ChrPos" = "SNP")) %>%
      select(VAR_ID, ChrPos, rsID, Variant_Type, original_SNP,
             REF, ALT, PVALUE, Risk_Allele, GWAS, Population)
  } else {
    final_proxies <- final_originals[0, ]
  }

  final_snps <- bind_rows(final_originals, final_proxies)
} else {
  message("Proxy search is disabled; retaining all position/LD-pruned originals.")
  final_snps <- pruned_vars %>%
    transmute(
      VAR_ID, ChrPos, rsID = RS_Number, Variant_Type = "original",
      original_SNP = VAR_ID, REF, ALT, PVALUE, Risk_Allele, GWAS, Population
    )
}

if (anyDuplicated(final_snps$original_SNP)) {
  stop("Final assembly produced more than one retained variant for an original sentinel.")
}
if (!setequal(final_snps$original_SNP, pruned_vars$VAR_ID)) {
  stop("Final assembly did not account for every pruned sentinel exactly once.")
}

write_tsv(
  final_snps %>% select(VAR_ID, rsID, Variant_Type, original_SNP),
  file.path(main_dir, "rsID_map.txt")
)

# =============================================================================
# 6. Primary-GWAS alignment and final matrices
# =============================================================================

gwas_alignment_source <- z_n_mats$df_gwas %>%
  transmute(
    ChrPos = SNP,
    primary_REF = REF,
    primary_ALT = ALT,
    primary_Risk_Allele = Risk_Allele,
    P_VALUE,
    BETA,
    SE
  )

gwas_final <- final_snps %>%
  inner_join(gwas_alignment_source, by = "ChrPos") %>%
  mutate(
    BETA_aligned = case_when(
      primary_ALT == primary_Risk_Allele ~ BETA,
      primary_REF == primary_Risk_Allele ~ -BETA,
      TRUE ~ NA_real_
    )
  )

if (nrow(gwas_final) != nrow(final_snps) || anyNA(gwas_final$BETA_aligned)) {
  stop("One or more final variants failed the primary-GWAS alignment audit.")
}

write_csv(gwas_final, file.path(main_dir, "alignment_GWAS_summStats.csv"))

zmat_pre_imputed <- zmat_fullset[gwas_final$ChrPos, traits_final, drop = FALSE]
Nmat_pre_imputed <- Nmat_fullset[gwas_final$ChrPos, traits_final, drop = FALSE]

# =============================================================================
# 7. Impute missing z-scores and sample sizes
# =============================================================================

set.seed(RANDOM_SEED)
observed_indices <- which(!is.na(zmat_pre_imputed), arr.ind = TRUE)
holdout_n <- floor(IMPUTATION_HOLDOUT * nrow(observed_indices))
if (holdout_n < 1L) {
  stop("Too few observed z-scores to construct an imputation holdout set.")
}
cv_indices <- observed_indices[
  sample(seq_len(nrow(observed_indices)), size = holdout_n),
  , drop = FALSE
]

zmat_cv <- zmat_pre_imputed
zmat_cv[cv_indices] <- NA_real_

rank_max <- max(1L, min(50L, nrow(zmat_pre_imputed) - 1L,
                        ncol(zmat_pre_imputed) - 1L))
lambda_max <- softImpute::lambda0(zmat_pre_imputed)
lambda_grid <- seq(0.1 * lambda_max, lambda_max, length.out = 10L)

cv_rmse <- vapply(lambda_grid, function(lambda_value) {
  fit <- softImpute::softImpute(
    zmat_cv,
    rank.max = rank_max,
    lambda = lambda_value,
    type = "svd"
  )
  imputed_cv <- softImpute::complete(zmat_cv, fit)
  sqrt(mean((zmat_pre_imputed[cv_indices] - imputed_cv[cv_indices])^2))
}, numeric(1L))

best_lambda <- lambda_grid[which.min(cv_rmse)]
imputation_fit <- softImpute::softImpute(
  zmat_pre_imputed,
  rank.max = rank_max,
  lambda = best_lambda,
  type = "svd"
)
zmat_imputed <- softImpute::complete(zmat_pre_imputed, imputation_fit)

Nmat_imputed <- apply(Nmat_pre_imputed, 2, function(values) {
  replace(values, is.na(values), median(values, na.rm = TRUE))
})
Nmat_imputed <- as.matrix(Nmat_imputed)
rownames(Nmat_imputed) <- rownames(Nmat_pre_imputed)

write_tsv(
  tibble(lambda = lambda_grid, holdout_RMSE = cv_rmse),
  file.path(main_dir, "softImpute_cross_validation.txt")
)
saveRDS(zmat_pre_imputed, file.path(main_dir, "zmat_preImputed.rds"))
saveRDS(zmat_imputed, file.path(main_dir, "zmat_imputed.rds"))

# =============================================================================
# 8. Prepare and run bNMF
# =============================================================================

# prep_z_matrix() writes trait_cor_mat.txt in the working directory. Run it in
# the results directory so the example never leaves generated files in scripts/.
prep_in_results_dir <- function(...) {
  previous_dir <- getwd()
  on.exit(setwd(previous_dir), add = TRUE)
  setwd(main_dir)
  prep_z_matrix(...)
}

prep_z_output <- prep_in_results_dir(
  z_mat = zmat_imputed,
  N_mat = Nmat_imputed,
  corr_cutoff = CORRELATION_CUTOFF
)

final_zscore_matrix <- as.matrix(prep_z_output$final_z_mat)
trait_log <- prep_z_output$df_traits
if (length(traits_removed_low_N) > 0L) {
  trait_log <- bind_rows(
    trait_log,
    tibble(trait = traits_removed_low_N,
           result = "removed (low median N)", note = NA_character_)
  )
}
if (length(traits_removed_missing) > 0L) {
  trait_log <- bind_rows(
    trait_log,
    tibble(trait = traits_removed_missing,
           result = "removed (high missingness)", note = NA_character_)
  )
}
write_csv(trait_log, file.path(main_dir, "df_traits.csv"))
saveRDS(final_zscore_matrix, file.path(main_dir, "z_score_mat.rds"))

bnmf_settings <- list(
  n_reps = BNMF_REPS,
  K = BNMF_K_INITIAL,
  K0 = BNMF_K0,
  tolerance = BNMF_TOLERANCE,
  phi = BNMF_PHI,
  random_seed = RANDOM_SEED,
  workers = BNMF_WORKERS
)
saveRDS(bnmf_settings, file.path(main_dir, "bnmf_settings.rds"))

future::plan(future::multisession, workers = BNMF_WORKERS)
bnmf_out <- run_bNMF_parallel(
  final_zscore_matrix,
  n_reps = BNMF_REPS,
  K = BNMF_K_INITIAL,
  K0 = BNMF_K0,
  tolerance = BNMF_TOLERANCE,
  phi = BNMF_PHI,
  random_seed = RANDOM_SEED
)
future::plan(future::sequential)
saveRDS(bnmf_out, file.path(main_dir, "bnmf_out.rds"))

summarize_bNMF(bnmf_out, dir_save = main_dir)
