# ============================================================
# Pooled adaptive Lasso estimator
#
# Initial estimator:
#   BIC-selected ordinary Lasso + weighted least-squares refit
#
# The refitted estimate is used to construct adaptive-Lasso
# penalty weights. The final adaptive-Lasso tuning parameter is
# selected by BIC.
#
# The pooled coefficient is replicated across all K groups in
# the returned K x p beta matrix.
# ============================================================

method_pooled_lasso <- function(mydat,
                                wloss = NULL,
                                gamma = 1,
                                nlambda = 50,
                                lambda.min.ratio = 0.02,
                                adaptive_eps = 1e-4,
                                zero_tol = 1e-6) {

  # ----------------------------------------------------------
  # 1. Extract data
  # ----------------------------------------------------------
  Yt <- as.numeric(mydat$Ytilde)
  At <- as.numeric(mydat$Atilde)
  X  <- as.matrix(mydat$X)

  K <- as.integer(mydat$K)
  p <- as.integer(mydat$p)
  n <- as.integer(mydat$n)

  storage.mode(X) <- "double"

  if (nrow(X) != n || ncol(X) != p) {
    stop(
      "Dimension mismatch in mydat$X. Expected ",
      n, " x ", p, ", got ",
      nrow(X), " x ", ncol(X), ".",
      call. = FALSE
    )
  }

  if (length(Yt) != n || length(At) != n) {
    stop(
      "Lengths of Ytilde and Atilde must equal mydat$n.",
      call. = FALSE
    )
  }

  if (any(!is.finite(Yt)) ||
      any(!is.finite(At)) ||
      any(!is.finite(X))) {
    stop(
      "Non-finite values were found in Ytilde, Atilde, or X.",
      call. = FALSE
    )
  }

  if (is.null(wloss)) {
    wloss <- rep(1, n)
  }

  wloss <- as.numeric(wloss)

  if (length(wloss) != n) {
    stop(
      "Length of wloss must equal mydat$n.",
      call. = FALSE
    )
  }

  if (any(!is.finite(wloss)) || any(wloss < 0)) {
    stop(
      "wloss must contain finite nonnegative values.",
      call. = FALSE
    )
  }

  if (all(wloss == 0)) {
    stop(
      "At least one element of wloss must be positive.",
      call. = FALSE
    )
  }

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop(
      "Package 'glmnet' is required.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # 2. Pooled R-learning design
  #
  # The pooled design has p columns and does not contain
  # group-specific blocks.
  # ----------------------------------------------------------
  Z <- X * At

  sw <- sqrt(wloss)
  Zw <- Z * sw
  yw <- Yt * sw

  # ----------------------------------------------------------
  # 3. Initial estimator:
  #    BIC-selected ordinary Lasso + refit
  # ----------------------------------------------------------

  # Ordinary Lasso path used to obtain the initial support
  fit_init_lasso <- glmnet::glmnet(
    x = Zw,
    y = yw,
    alpha = 1,
    intercept = FALSE,
    standardize = FALSE,
    nlambda = nlambda,
    lambda.min.ratio = lambda.min.ratio
  )

  nlam_init <- length(fit_init_lasso$lambda)
  BIC_init <- rep(NA_real_, nlam_init)

  if (!is.null(fit_init_lasso$beta)) {

    Bmat_init <- as.matrix(fit_init_lasso$beta)

    if (nrow(Bmat_init) != p && ncol(Bmat_init) == p) {
      Bmat_init <- t(Bmat_init)
    }

    if (nrow(Bmat_init) != p) {
      stop(
        "Unexpected dimension of initial Lasso coefficients.",
        call. = FALSE
      )
    }

    Eta_init <- Z %*% Bmat_init
    DF_init <- colSums(abs(Bmat_init) > zero_tol)

    for (ll in seq_len(nlam_init)) {
      BIC_init[ll] <- bic_lm_loss(
        y = Yt,
        eta = Eta_init[, ll],
        df = DF_init[ll],
        n = n
      )
    }

    id_init <- which.min(BIC_init)
    beta_lasso_init <- as.numeric(Bmat_init[, id_init])

  } else {

    for (ll in seq_len(nlam_init)) {
      beta_ll <- as.numeric(
        stats::coef(
          fit_init_lasso,
          s = fit_init_lasso$lambda[ll]
        )
      )[-1]

      eta_ll <- drop(Z %*% beta_ll)
      df_ll <- sum(abs(beta_ll) > zero_tol)

      BIC_init[ll] <- bic_lm_loss(
        y = Yt,
        eta = eta_ll,
        df = df_ll,
        n = n
      )
    }

    id_init <- which.min(BIC_init)

    beta_lasso_init <- as.numeric(
      stats::coef(
        fit_init_lasso,
        s = fit_init_lasso$lambda[id_init]
      )
    )[-1]
  }

  beta_lasso_init[
    abs(beta_lasso_init) <= zero_tol
  ] <- 0

  # ----------------------------------------------------------
  # 4. Weighted least-squares refit on the BIC-Lasso support
  # ----------------------------------------------------------

  support_init <- which(abs(beta_lasso_init) > zero_tol)

  beta_init <- numeric(p)

  if (length(support_init) > 0) {

    Z_support_w <- Zw[, support_init, drop = FALSE]

    # QR-based least-squares refit.
    # qr.solve may fail if the selected design is rank deficient,
    # in which case a generalized inverse is used.
    beta_support <- tryCatch(
      {
        as.numeric(
          qr.solve(
            crossprod(Z_support_w),
            crossprod(Z_support_w, yw)
          )
        )
      },
      error = function(e) {
        if (!requireNamespace("MASS", quietly = TRUE)) {
          stop(
            "The selected Lasso support is rank deficient and ",
            "package 'MASS' is required for generalized-inverse refitting.",
            call. = FALSE
          )
        }

        as.numeric(
          MASS::ginv(crossprod(Z_support_w)) %*%
            crossprod(Z_support_w, yw)
        )
      }
    )

    beta_init[support_init] <- beta_support
  }

  beta_init[
    abs(beta_init) <= zero_tol
  ] <- 0

  # ----------------------------------------------------------
  # 5. Construct adaptive-Lasso penalty factors
  # ----------------------------------------------------------
  penalty_factor <- 1 / (abs(beta_init)^gamma + adaptive_eps)

  # ----------------------------------------------------------
  # 6. Final adaptive Lasso
  # ----------------------------------------------------------
  fit <- glmnet::glmnet(
    x = Zw,
    y = yw,
    alpha = 1,
    intercept = FALSE,
    standardize = FALSE,
    penalty.factor = penalty_factor,
    nlambda = nlambda,
    lambda.min.ratio = lambda.min.ratio
  )

  nlam <- length(fit$lambda)
  BIC <- rep(NA_real_, nlam)

  # ----------------------------------------------------------
  # 7. Select final adaptive-Lasso lambda using BIC
  #
  # BIC is evaluated on the original residualized scale.
  # ----------------------------------------------------------
  if (!is.null(fit$beta)) {

    Bmat <- as.matrix(fit$beta)

    if (nrow(Bmat) != p && ncol(Bmat) == p) {
      Bmat <- t(Bmat)
    }

    if (nrow(Bmat) != p) {
      stop(
        "Unexpected dimension of fit$beta from glmnet.",
        call. = FALSE
      )
    }

    Eta <- Z %*% Bmat
    DF <- colSums(abs(Bmat) > zero_tol)

    for (ll in seq_len(nlam)) {
      BIC[ll] <- bic_lm_loss(
        y = Yt,
        eta = Eta[, ll],
        df = DF[ll],
        n = n
      )
    }

    id <- which.min(BIC)
    beta_pooled <- as.numeric(Bmat[, id])

  } else {

    for (ll in seq_len(nlam)) {

      beta_ll <- as.numeric(
        stats::coef(
          fit,
          s = fit$lambda[ll]
        )
      )[-1]

      eta_ll <- drop(Z %*% beta_ll)
      df_ll <- sum(abs(beta_ll) > zero_tol)

      BIC[ll] <- bic_lm_loss(
        y = Yt,
        eta = eta_ll,
        df = df_ll,
        n = n
      )
    }

    id <- which.min(BIC)

    beta_pooled <- as.numeric(
      stats::coef(
        fit,
        s = fit$lambda[id]
      )
    )[-1]
  }

  # Numerical cleanup
  beta_pooled[
    abs(beta_pooled) <= zero_tol
  ] <- 0

  # ----------------------------------------------------------
  # 8. Replicate pooled estimate across groups
  # ----------------------------------------------------------
  beta_hat <- matrix(
    rep(beta_pooled, times = K),
    nrow = K,
    ncol = p,
    byrow = TRUE
  )

  # ----------------------------------------------------------
  # 9. Return
  # ----------------------------------------------------------
  list(
    beta = beta_hat,
    beta_pooled = beta_pooled,

    # BIC-Lasso + refit initial estimate
    beta_init = beta_init,
    beta_lasso_init = beta_lasso_init,
    support_init = support_init,
    lambda_init = fit_init_lasso$lambda[id_init],
    lambda_init_index = id_init,
    bic_init = BIC_init[id_init],
    bic_init_path = BIC_init,
    lambda_init_path = fit_init_lasso$lambda,

    # Adaptive-Lasso information
    penalty_factor = penalty_factor,
    lambda = fit$lambda[id],
    lambda_index = id,
    bic = BIC[id],
    bic_path = BIC,
    lambda_path = fit$lambda,
    df = sum(abs(beta_pooled) > zero_tol),

    sigma2 = estimate_group_variance(mydat, beta_hat),

    # Model objects
    fit_init_lasso = fit_init_lasso,
    fit = fit
  )
}
