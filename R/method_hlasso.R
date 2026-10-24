# ============================================================
# method_hlasso.R
#
# Hierarchical / adaptive lasso for subgroup-specific effects.
#
# Raw-data / pseudo-outcome analogue of the hierarchical lasso
# penalty in He et al. (2015), Section 2.2 and 2.4.
#
# Target penalty:
#
#   lambda * sum_j sqrt( sum_k W_jk * |beta_kj| )
#
# Variational representation:
#
#   sqrt(S_j) = min_{gamma_j > 0} 0.5 * (S_j / gamma_j + gamma_j)
#
# Hence:
#
#   gamma_j = sqrt(S_j)
#
# and the beta-step is a weighted lasso with penalty factor:
#
#   0.5 * W_jk / gamma_j
#
# Acceleration:
#   If beta_init[, j] is zero for all subgroups, variable j is
#   removed from the optimization and restored as zero in output.
# ============================================================

# ------------------------------------------------------------
# Safe coefficient extractor for glmnet without intercept
# ------------------------------------------------------------
get_glmnet_coef_no_intercept_hlasso <- function(fit, pZ) {

  # Preferred route: glmnet stores beta as pZ x nlambda sparse matrix.
  # In hLASSO we fit one lambda at a time, so usually ncol(beta) = 1.
  if (!is.null(fit$beta)) {
    Bmat <- as.matrix(fit$beta)

    if (nrow(Bmat) == pZ) {
      return(as.numeric(Bmat[, 1]))
    }

    if (ncol(Bmat) == pZ) {
      return(as.numeric(Bmat[1, ]))
    }
  }

  # Fallback route: coef() may include an intercept entry.
  b_tmp <- as.numeric(stats::coef(fit))

  if (length(b_tmp) == pZ + 1) {
    return(b_tmp[-1])
  }

  if (length(b_tmp) == pZ) {
    return(b_tmp)
  }

  stop(
    "Unexpected coefficient length from glmnet in hLASSO. Expected ",
    pZ, " or ", pZ + 1,
    ", got ", length(b_tmp), ".",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Compute glmnet-style lambda_max for gaussian lasso
#
# glmnet gaussian objective is approximately:
#
#   1/(2n) * ||y - X beta||^2
#   + lambda * sum_j penalty.factor_j * |beta_j|
#
# glmnet also rescales penalty.factor internally so that
# sum(penalty.factor) = number of variables.
# ------------------------------------------------------------
glmnet_lambda_max <- function(X, y, penalty.factor) {
  n <- nrow(X)
  pZ <- ncol(X)

  pf <- as.numeric(penalty.factor)

  if (length(pf) != pZ) {
    stop(
      "Length of penalty.factor must equal ncol(X). Expected ",
      pZ, ", got ", length(pf), ".",
      call. = FALSE
    )
  }

  if (any(!is.finite(pf)) || any(pf < 0)) {
    stop("penalty.factor contains non-finite or negative values.", call. = FALSE)
  }

  if (any(pf == 0)) {
    stop(
      "This helper assumes all coefficients are penalized. ",
      "Found zero penalty.factor.",
      call. = FALSE
    )
  }

  # glmnet rescales penalty.factor to sum to pZ
  pf_adj <- pf * pZ / sum(pf)

  score <- abs(drop(crossprod(X, y))) / n

  lambda_max <- max(score / pf_adj)

  if (!is.finite(lambda_max) || lambda_max <= 0) {
    stop("Computed lambda_max is invalid.", call. = FALSE)
  }

  lambda_max
}

# ------------------------------------------------------------
# Main hLASSO method
# ------------------------------------------------------------
method_hlasso <- function(mydat,
                          beta_init = NULL,
                          wloss = NULL,
                          gamma_adapt = 1,
                          nlambda = 50,
                          lambda.min.ratio = 1e-4,
                          lambda.max.multiplier = 10,
                          lambda_seq = NULL,
                          maxit_outer = 50,
                          tol_outer = 1e-4,
                          eps_gamma = 1e-8,
                          eps_weight = 1e-4,
                          use_half_factor = TRUE,
                          verbose = FALSE,
                          zero_init_tol = 0) {

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop(
      "Package 'glmnet' is required. Please install it first:\n",
      "install.packages('glmnet')",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Extract full data
  # ----------------------------------------------------------
  Yt_full <- mydat$Ytilde
  At_full <- mydat$Atilde
  X_full  <- mydat$X
  G_full  <- mydat$G
  K <- mydat$K
  p_full <- mydat$p
  n <- mydat$n

  if (is.null(wloss)) {
    wloss <- rep(1, n)
  }

  if (length(wloss) != n) {
    stop(
      "Length of wloss must equal mydat$n. Expected ",
      n, ", got ", length(wloss), ".",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Initial beta
  # ----------------------------------------------------------
  if (is.null(beta_init)) {
    beta_init <- method_glasso(mydat, wloss = wloss)$beta
  }

  beta_init <- as.matrix(beta_init)

  if (!(nrow(beta_init) == K && ncol(beta_init) == p_full)) {
    if (nrow(beta_init) == p_full && ncol(beta_init) == K) {
      beta_init <- t(beta_init)
    } else {
      stop(
        "Dimension mismatch for beta_init in method_hlasso(). ",
        "Expected K x p = ", K, " x ", p_full,
        " or p x K = ", p_full, " x ", K,
        ", got ", nrow(beta_init), " x ", ncol(beta_init), ".",
        call. = FALSE
      )
    }
  }

  # ----------------------------------------------------------
  # Screen variables whose initial estimates are zero in all groups
  # ----------------------------------------------------------
  screen <- screen_zero_init_vars(
    mydat = mydat,
    beta_init = beta_init,
    zero_init_tol = zero_init_tol,
    caller = "method_hlasso()"
  )

  active_j <- screen$active_j
  dropped_j <- screen$dropped_j
  p_active <- screen$p_active

  if (p_active == 0) {
    beta_zero <- matrix(0, nrow = K, ncol = p_full)

    return(list(
      beta = beta_zero,
      beta_active = matrix(0, nrow = K, ncol = 0),
      lambda = NA_real_,
      bic = bic_zero_fit(mydat),
      df = 0L,
      beta_path = list(beta_zero),
      beta_active_path = list(matrix(0, nrow = K, ncol = 0)),
      bic_path = bic_zero_fit(mydat),
      df_path = 0L,
      lambda_path = NA_real_,
      iter_path = 0L,
      Wjk = matrix(0, nrow = K, ncol = 0),
      Wjk_full = matrix(Inf, nrow = K, ncol = p_full),
      pf0 = numeric(0),
      lambda.max.multiplier = lambda.max.multiplier,
      lambda.min.ratio = lambda.min.ratio,
      use_half_factor = use_half_factor,
      active_j = active_j,
      dropped_j = dropped_j,
      p_active = 0L,
      p_full = p_full,
      zero_init_tol = zero_init_tol
    ))
  }

  # Use screened data internally.
  mydat_active <- screen$mydat
  beta_init_active <- screen$beta_init

  Yt <- mydat_active$Ytilde
  At <- mydat_active$Atilde
  X  <- mydat_active$X
  G  <- mydat_active$G
  p  <- mydat_active$p

  # ----------------------------------------------------------
  # Build design matrix only for active variables
  # ----------------------------------------------------------
  Z <- build_block_design(X, At, G, K)

  wy <- weighted_xy(Z, Yt, wloss)
  Zw <- wy$Z
  yw <- wy$y

  pZ <- ncol(Z)

  if (pZ != K * p) {
    stop(
      "Unexpected design dimension in method_hlasso(). ",
      "Expected ncol(Z) = K * p = ", K * p,
      ", got ", pZ, ".",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Adaptive weights
  # ----------------------------------------------------------
  Wjk <- 1 / (abs(beta_init_active)^gamma_adapt + eps_weight)

  if (any(!is.finite(Wjk))) {
    stop("Non-finite adaptive weights Wjk detected.", call. = FALSE)
  }

  # ----------------------------------------------------------
  # Initial hLASSO penalty factor for lambda_max calculation
  #
  # S_j = sum_k W_jk |beta_kj|
  # gamma_j = sqrt(S_j)
  # PF_kj = 0.5 * W_jk / gamma_j
  # ----------------------------------------------------------
  S0 <- rep(0, p)

  for (j in seq_len(p)) {
    S0[j] <- sum(Wjk[, j] * abs(beta_init_active[, j]))
  }

  gamma0 <- sqrt(pmax(S0, eps_gamma))
  gamma0 <- pmax(gamma0, sqrt(eps_gamma))

  PF0 <- matrix(0, nrow = K, ncol = p)

  for (j in seq_len(p)) {
    PF0[, j] <- Wjk[, j] / gamma0[j]
  }

  if (isTRUE(use_half_factor)) {
    PF0 <- 0.5 * PF0
  }

  pf0 <- as.vector(t(PF0))

  if (length(pf0) != pZ) {
    stop(
      "Length mismatch for initial hLASSO penalty.factor. Expected ",
      pZ, ", got ", length(pf0), ".",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Lambda sequence
  # ----------------------------------------------------------
  if (is.null(lambda_seq)) {
    lambda_max0 <- glmnet_lambda_max(
      X = Zw,
      y = yw,
      penalty.factor = pf0
    )

    lambda_max <- lambda.max.multiplier * lambda_max0
    lambda_min <- lambda_max * lambda.min.ratio

    lambda_seq <- exp(seq(
      log(lambda_max),
      log(lambda_min),
      length.out = nlambda
    ))
  } else {
    lambda_seq <- as.numeric(lambda_seq)
    lambda_seq <- lambda_seq[is.finite(lambda_seq) & lambda_seq > 0]

    if (length(lambda_seq) == 0) {
      stop("User-supplied lambda_seq is empty or invalid.", call. = FALSE)
    }
  }

  lambda_seq <- sort(unique(lambda_seq), decreasing = TRUE)

  beta_store <- vector("list", length(lambda_seq))
  bic_store  <- rep(NA_real_, length(lambda_seq))
  df_store   <- rep(NA_integer_, length(lambda_seq))
  iter_store <- rep(NA_integer_, length(lambda_seq))

  # ----------------------------------------------------------
  # Main lambda loop
  # ----------------------------------------------------------
  for (ll in seq_along(lambda_seq)) {

    lam <- lambda_seq[ll]

    # Important:
    # Do NOT initialize from previous lambda solution.
    # For nonconvex hLASSO, if the previous larger-lambda solution
    # is exactly zero, warm-starting from it can trap the path at zero.
    B_old <- beta_init_active

    it_used <- maxit_outer

    for (it in seq_len(maxit_outer)) {

      # ------------------------------------------------------
      # Update gamma
      # ------------------------------------------------------
      S_j <- rep(0, p)

      for (j in seq_len(p)) {
        S_j[j] <- sum(Wjk[, j] * abs(B_old[, j]))
      }

      gamma_j <- sqrt(pmax(S_j, eps_gamma))
      gamma_j <- pmax(gamma_j, sqrt(eps_gamma))

      # ------------------------------------------------------
      # Dynamic penalty.factor
      # ------------------------------------------------------
      PF <- matrix(0, nrow = K, ncol = p)

      for (j in seq_len(p)) {
        PF[, j] <- Wjk[, j] / gamma_j[j]
      }

      if (isTRUE(use_half_factor)) {
        PF <- 0.5 * PF
      }

      pf_vec <- as.vector(t(PF))

      if (length(pf_vec) != pZ) {
        stop(
          "Length mismatch for penalty.factor in method_hlasso(). ",
          "Expected ", pZ, ", got ", length(pf_vec), ".",
          call. = FALSE
        )
      }

      if (any(!is.finite(pf_vec))) {
        stop("Non-finite penalty.factor detected.", call. = FALSE)
      }

      # ------------------------------------------------------
      # Weighted lasso step
      # ------------------------------------------------------
      fit <- glmnet::glmnet(
        x = Zw,
        y = yw,
        alpha = 1,
        lambda = lam,
        intercept = FALSE,
        standardize = FALSE,
        penalty.factor = pf_vec
      )

      b_new <- get_glmnet_coef_no_intercept_hlasso(
        fit = fit,
        pZ = pZ
      )

      B_new <- mat_beta(b_new, K, p)

      if (!(nrow(B_new) == K && ncol(B_new) == p)) {
        stop(
          "mat_beta() returned unexpected dimension. Expected ",
          K, " x ", p,
          ", got ", nrow(B_new), " x ", ncol(B_new), ".",
          call. = FALSE
        )
      }

      diff <- sqrt(sum((B_new - B_old)^2)) /
        (sqrt(sum(B_old^2)) + 1e-8)

      B_old <- B_new

      if (diff < tol_outer) {
        it_used <- it
        break
      }
    }

    # --------------------------------------------------------
    # Evaluation
    # --------------------------------------------------------
    b_final <- vec_beta(B_old)

    if (length(b_final) != pZ) {
      stop(
        "vec_beta() returned unexpected length. Expected ",
        pZ, ", got ", length(b_final), ".",
        call. = FALSE
      )
    }

    eta <- drop(Z %*% b_final)

    df <- sum(abs(B_old) > 1e-6)

    bic_store[ll] <- bic_lm_loss(Yt, eta, df, n)

    beta_store[[ll]] <- B_old
    df_store[ll] <- df
    iter_store[ll] <- it_used

    if (isTRUE(verbose)) {
      message(
        "lambda[", ll, "/", length(lambda_seq), "] = ",
        signif(lam, 4),
        ", df = ", df,
        ", bic = ", signif(bic_store[ll], 5),
        ", iter = ", it_used
      )
    }
  }

  id <- which.min(bic_store)

  beta_active_path <- beta_store
  beta_path_full <- lapply(
    beta_active_path,
    restore_beta_screen,
    active_j = active_j,
    K = K
  )

  beta_active <- beta_active_path[[id]]
  beta_full <- beta_path_full[[id]]

  Wjk_full <- matrix(Inf, nrow = K, ncol = p_full)
  Wjk_full[, active_j] <- Wjk

  list(
    beta = beta_full,
    beta_active = beta_active,
    lambda = lambda_seq[id],
    bic = bic_store[id],
    df = df_store[id],
    beta_path = beta_path_full,
    beta_active_path = beta_active_path,
    bic_path = bic_store,
    df_path = df_store,
    lambda_path = lambda_seq,
    iter_path = iter_store,
    Wjk = Wjk,
    Wjk_full = Wjk_full,
    pf0 = pf0,
    lambda.max.multiplier = lambda.max.multiplier,
    lambda.min.ratio = lambda.min.ratio,
    use_half_factor = use_half_factor,
    active_j = active_j,
    dropped_j = dropped_j,
    p_active = p_active,
    p_full = p_full,
    zero_init_tol = zero_init_tol
  )
}