# Study-local optimized implementation of the L2/Euclidean Bayesian NMF used by
# the T2D clustering pipeline.
#
# The statistical updates are unchanged.  The expensive denominator products
# are evaluated through K x K Gram matrices:
#   t(W) %*% (W %*% H) == crossprod(W) %*% H
#   (W %*% H) %*% t(H) == W %*% tcrossprod(H)
# Reconstruction error is calculated with a Frobenius-norm identity, avoiding
# construction of the full W %*% H matrix on every iteration.


BayesNMF.L2EU_optimized <- function(
    V0, n.iter = 10000, a0 = 10, tol = 1e-7, K = 15, K0 = 15,
    phi = 1.0, window_size = 25, min_iter = 100, verbose = FALSE) {

  eps <- 1e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes, drop = FALSE]
  V <- V0 - min(V0)
  Vmax <- max(V)
  N <- nrow(V)
  M <- ncol(V)

  W <- matrix(runif(N * K) * Vmax, ncol = K)
  H <- matrix(runif(M * K) * Vmax, ncol = M)

  phi <- sd(V)^2 * phi
  C <- (N + M) / 2 + a0 + 1
  b0 <- 3.14 * (a0 - 1) * mean(V) / (2 * K0)
  lambda.bound <- b0 / C
  lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
  lambda.cut <- lambda.bound * 1.5
  sum_V_squared <- sum(V^2)

  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.active <- list()
  n.lambda[[1L]] <- lambda

  check_convergence <- function(errors, window_size) {
    if (length(errors) < window_size) return(FALSE)
    recent <- tail(unlist(errors), window_size)
    rel_change <- abs(mean(diff(recent)) / mean(recent))
    is.finite(rel_change) && rel_change < tol
  }

  iter <- 2L
  while ((iter <= min_iter || del >= tol) && iter < n.iter) {
    # H update. crossprod(W) %*% H is algebraically identical to
    # t(W) %*% (W %*% H), but avoids an N x M reconstruction.
    WtV <- crossprod(W, V)
    WtW <- crossprod(W)
    H <- H * WtV /
      (WtW %*% H + phi * sweep(H, 1L, lambda, "/") + eps)

    # W update. W %*% tcrossprod(H) is algebraically identical to
    # (W %*% H) %*% t(H).
    VHt <- V %*% t(H)
    HHt <- tcrossprod(H)
    W <- W * VHt /
      (W %*% HHt + phi * sweep(W, 2L, lambda, "/") + eps)

    lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
    lambda[!is.finite(lambda)] <- 1e-6

    previous_lambda <- n.lambda[[iter - 1L]]
    del <- if (!is.null(previous_lambda) && all(is.finite(previous_lambda))) {
      max(abs(lambda - previous_lambda) / (previous_lambda + 1e-10))
    } else {
      Inf
    }

    # ||V-WH||_F^2 = ||V||_F^2 - 2<tr(W'V H')> +
    #                  tr(W'W H H'). VHt and HHt are already available.
    WtW_updated <- crossprod(W)
    error <- sum_V_squared - 2 * sum(W * VHt) + sum(WtW_updated * HHt)
    # Protect against a tiny negative value from floating-point cancellation.
    error <- max(0, error)
    like <- error / 2

    n.like[[iter]] <- like
    n.evid[[iter]] <- like + phi * sum(
      (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / lambda +
        C * log(lambda)
    )
    n.error[[iter]] <- error
    n.lambda[[iter]] <- lambda
    n.active[[iter]] <- sum(lambda >= lambda.cut)

    if (verbose && iter %% 100L == 0L) {
      cat(sprintf(
        "Iter %d: Error=%g, Delta=%g, Active=%d, Evidence=%g\n",
        iter, error, del, n.active[[iter]], n.evid[[iter]]
      ))
    }

    if (iter > min_iter && iter %% 20L == 0L &&
        check_convergence(n.error, window_size)) {
      if (verbose) cat("Converged based on error stability\n")
      break
    }
    iter <- iter + 1L
  }

  final_active_components <- which(lambda >= lambda.cut)
  if (length(final_active_components) == 0L) {
    stop("bNMF ended without any active components.")
  }

  list(
    W = W,
    H = H,
    n.like = n.like,
    n.evid = n.evid,
    n.lambda = n.lambda,
    n.error = n.error,
    n.active = n.active,
    lambda_cut = lambda.cut,
    active_components = final_active_components,
    iterations = iter - 1L,
    converged = del < tol,
    implementation = "gram-matrix optimized L2EU v2"
  )
}


run_bNMF_parallel_checkpointed <- function(
    z_mat, n_reps = 10, random_seed = 1, K = 20, K0 = 10,
    tolerance = 1e-7, phi = 1, checkpoint_dir,
    n_workers = future::availableCores()) {

  if (missing(checkpoint_dir) || !nzchar(checkpoint_dir)) {
    stop("checkpoint_dir is required")
  }
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- normalizePath(checkpoint_dir)
  # The production caller validates this against UGER's NSLOTS allocation.
  # Do not clamp again with future::availableCores(), which may under-report
  # interactive SGE allocations because of unrelated local heuristics.
  n_workers <- max(1L, as.integer(n_workers))

  object_md5 <- function(object) {
    hash_file <- tempfile("bnmf_signature_", fileext = ".bin")
    on.exit(unlink(hash_file), add = TRUE)
    writeBin(serialize(object, NULL, version = 3), hash_file)
    unname(tools::md5sum(hash_file))
  }
  run_signature <- object_md5(list(
    implementation = "gram-matrix optimized L2EU v2",
    z_mat = z_mat,
    n_reps = n_reps,
    random_seed = random_seed,
    K = K,
    K0 = K0,
    tolerance = tolerance,
    phi = phi
  ))

  set.seed(random_seed)
  repetition_seeds <- sample.int(.Machine$integer.max, n_reps)
  checkpoint_path <- function(rep) {
    file.path(checkpoint_dir, sprintf("optimized_rep_%03d.rds", rep))
  }
  checkpoint_is_valid <- function(rep) {
    path <- checkpoint_path(rep)
    if (!file.exists(path)) return(FALSE)
    tryCatch({
      result <- readRDS(path)
      identical(result$rep, rep) &&
        identical(result$seed, repetition_seeds[rep]) &&
        identical(result$implementation, "gram-matrix optimized L2EU v2") &&
        identical(result$run_signature, run_signature) &&
        is.matrix(result$W) && is.matrix(result$H)
    }, error = function(e) FALSE)
  }
  completed <- which(vapply(seq_len(n_reps), function(rep) {
    checkpoint_is_valid(rep)
  }, logical(1)))
  pending <- setdiff(seq_len(n_reps), completed)

  message(sprintf(
    "Optimized bNMF: %d repetitions total, %d complete, %d pending, %d worker(s).",
    n_reps, length(completed), length(pending), n_workers
  ))
  if (length(completed)) {
    message("Resuming from checkpoints in: ", checkpoint_dir)
  }

  run_one <- function(rep) {
    set.seed(repetition_seeds[rep])
    started <- Sys.time()
    result <- BayesNMF.L2EU_optimized(
      V0 = z_mat,
      K = K,
      K0 = K0,
      tol = tolerance,
      phi = phi
    )
    result$rep <- rep
    result$seed <- repetition_seeds[rep]
    result$run_signature <- run_signature
    result$elapsed_seconds <- as.numeric(difftime(Sys.time(), started,
                                                  units = "secs"))
    final_path <- checkpoint_path(rep)
    temporary_path <- paste0(final_path, ".tmp-", Sys.getpid())
    saveRDS(result, temporary_path)
    if (!file.rename(temporary_path, final_path)) {
      stop("Could not atomically finalize checkpoint: ", final_path)
    }
    list(rep = rep, elapsed_seconds = result$elapsed_seconds)
  }

  if (length(pending)) {
    if (!requireNamespace("progressr", quietly = TRUE) ||
        !requireNamespace("progress", quietly = TRUE)) {
      stop("The local progress display requires the progressr and progress packages")
    }
    progress_handler <- progressr::handler_progress(
      format = "[:bar] :percent (:current/:total) :message",
      show_after = 0,
      clear = FALSE
    )

    progressr::with_progress(
      handlers = progress_handler,
      enable = TRUE,
      {
        progress <- progressr::progressor(steps = n_reps)
        if (length(completed)) {
          progress(
            amount = length(completed),
            message = sprintf("resumed %d completed repetition(s)",
                              length(completed))
          )
        }

        run_one_with_progress <- function(rep) {
          info <- run_one(rep)
          progress(message = sprintf(
            "finished rep %d/%d in %.1f min",
            info$rep, n_reps, info$elapsed_seconds / 60
          ))
          info$rep
        }

        if (n_workers > 1L) {
          future::plan(future::multisession, workers = n_workers)
          furrr::future_map(
            pending,
            run_one_with_progress,
            # Submit one future per repetition instead of one ten-repetition
            # chunk per worker. This lets progress conditions reach the main
            # process as each checkpoint completes and improves load balancing.
            .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
          )
          future::plan(future::sequential)
        } else {
          lapply(pending, run_one_with_progress)
        }
      }
    )
  }

  paths <- vapply(seq_len(n_reps), checkpoint_path, character(1))
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths)) {
    stop("Missing bNMF checkpoint(s): ", paste(missing_paths, collapse = ", "))
  }
  lapply(paths, readRDS)
}
