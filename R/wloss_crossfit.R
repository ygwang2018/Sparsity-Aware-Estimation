# ============================================================
# wloss_crossfit.R
#
# Cross-fitted subgroup inverse-variance loss weights.
#
# For each subgroup k:
#   1. Split observations with G = k into L folds.
#   2. For each fold, fit a local treatment-effect model on the
#      remaining L - 1 folds within subgroup k.
#   3. Estimate residual variance on the training part.
#   4. Assign wloss = 1 / sigma2_train to the held-out fold.
#
# This avoids using each observation's own residual to construct
# its own loss weight.
# ============================================================

# ------------------------------------------------------------
# Local adaptive lasso fit for one subgroup
#
# Model:
#   Ytilde_i = Atilde_i * X_i^T beta_k + error
#
# No intercept.
# ------------------------------------------------------------
fit_local_lasso_one_group <- function(Yt,
                                      At,
                                      X,
                                      gamma = 1,
                                      nlambda = 50,
                                      lambda.min.ratio = 0.02,
                                      eps_weight = 1e-4) {

  Yt <- as.numeric(Yt)
  At <- as.numeric(At)
  X <- as.matrix(X)

  n <- nrow(X)
  p <- ncol(X)

  if (n <= 1) {
    return(rep(0, p))
  }

  Z <- X * At

  # If design is degenerate, return zero beta.
  if (all(abs(Z) < 1e-12)) {
    return(rep(0, p))
  }

  fit0 <- tryCatch(
    {
      glmnet::glmnet(
        x = Z,
        y = Yt,
        alpha = 0,
        intercept = FALSE,
        standardize = FALSE
      )
    },
    error = function(e) NULL
  )

  if (is.null(fit0)) {
    return(rep(0, p))
  }

  if (!is.null(fit0$beta)) {
    b0 <- as.numeric(fit0$beta[, ncol(fit0$beta)])
  } else {
    b0 <- as.numeric(stats::coef(
      fit0,
      s = fit0$lambda[length(fit0$lambda)]
    ))[-1]
  }

  if (length(b0) != p || any(!is.finite(b0))) {
    b0 <- rep(0, p)
  }

  pf <- 1 / (abs(b0)^gamma + eps_weight)

  fit <- tryCatch(
    {
      glmnet::glmnet(
        x = Z,
        y = Yt,
        alpha = 1,
        intercept = FALSE,
        standardize = FALSE,
        penalty.factor = pf,
        nlambda = nlambda,
        lambda.min.ratio = lambda.min.ratio
      )
    },
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(b0)
  }

  nlam <- length(fit$lambda)
  BIC <- rep(NA_real_, nlam)

  if (!is.null(fit$beta)) {
    Bmat <- as.matrix(fit$beta)

    if (nrow(Bmat) != p && ncol(Bmat) == p) {
      Bmat <- t(Bmat)
    }

    if (nrow(Bmat) != p) {
      return(b0)
    }

    Eta <- Z %*% Bmat
    DF <- colSums(abs(Bmat) > 1e-6)

    for (ll in seq_len(nlam)) {
      BIC[ll] <- bic_lm_loss(Yt, Eta[, ll], DF[ll], n)
    }

    id <- which.min(BIC)
    bhat <- as.numeric(Bmat[, id])
  } else {
    for (ll in seq_len(nlam)) {
      b <- as.numeric(stats::coef(fit, s = fit$lambda[ll]))[-1]
      eta <- drop(Z %*% b)
      df <- sum(abs(b) > 1e-6)
      BIC[ll] <- bic_lm_loss(Yt, eta, df, n)
    }

    id <- which.min(BIC)
    bhat <- as.numeric(stats::coef(fit, s = fit$lambda[id]))[-1]
  }

  if (length(bhat) != p || any(!is.finite(bhat))) {
    bhat <- rep(0, p)
  }

  bhat
}

# ------------------------------------------------------------
# Main function: cross-fitted wloss from local subgroup fits
# ------------------------------------------------------------
make_wloss_cf_local <- function(mydat,
                                nfold = 5,
                                seed = 1,
                                gamma = 1,
                                nlambda = 50,
                                lambda.min.ratio = 0.02,
                                eps_weight = 1e-4,
                                sigma_floor = 1e-6,
                                sigma_cap = Inf,
                                w_floor = 1e-3,
                                w_cap = 1e3,
                                denom = c("n", "n-df"),
                                normalize = TRUE,
                                verbose = FALSE) {

  denom <- match.arg(denom)

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for make_wloss_cf_local().", call. = FALSE)
  }

  Yt <- as.numeric(mydat$Ytilde)
  At <- as.numeric(mydat$Atilde)
  X  <- as.matrix(mydat$X)
  G  <- as.integer(mydat$G)
  K  <- as.integer(mydat$K)
  p  <- as.integer(mydat$p)
  n  <- as.integer(mydat$n)

  if (length(Yt) != n || length(At) != n || length(G) != n || nrow(X) != n) {
    stop("Dimension mismatch in make_wloss_cf_local().", call. = FALSE)
  }

  if (nfold < 2) {
    stop("nfold must be at least 2.", call. = FALSE)
  }

  if (any(G < 1 | G > K)) {
    stop("G must contain subgroup labels from 1 to K.", call. = FALSE)
  }

  set.seed(seed)

  wloss <- rep(NA_real_, n)
  sigma2_oof <- rep(NA_real_, n)
  foldid_group <- rep(NA_integer_, n)

  beta_fold <- vector("list", K)
  sigma2_fold <- vector("list", K)

  for (k in seq_len(K)) {
    idx_k <- which(G == k)
    nk <- length(idx_k)

    beta_fold[[k]] <- list()
    sigma2_fold[[k]] <- numeric(0)

    if (nk == 0) {
      next
    }

    Lk <- min(nfold, nk)

    if (Lk < 2) {
      sigma2_k <- mean(Yt[idx_k]^2)
      sigma2_k <- min(max(sigma2_k, sigma_floor), sigma_cap)

      sigma2_oof[idx_k] <- sigma2_k
      wloss[idx_k] <- 1 / sigma2_k
      foldid_group[idx_k] <- 1L

      next
    }

    fold_k <- sample(rep(seq_len(Lk), length.out = nk))
    foldid_group[idx_k] <- fold_k

    if (isTRUE(verbose)) {
      message("Group ", k, ": nk = ", nk, ", folds = ", Lk)
    }

    for (ff in seq_len(Lk)) {
      idx_test <- idx_k[fold_k == ff]
      idx_train <- idx_k[fold_k != ff]

      if (length(idx_train) <= 1) {
        sigma2_train <- mean(Yt[idx_train]^2)
        beta_k <- rep(0, p)
      } else {
        beta_k <- fit_local_lasso_one_group(
          Yt = Yt[idx_train],
          At = At[idx_train],
          X = X[idx_train, , drop = FALSE],
          gamma = gamma,
          nlambda = nlambda,
          lambda.min.ratio = lambda.min.ratio,
          eps_weight = eps_weight
        )

        eta_train <- At[idx_train] *
          drop(X[idx_train, , drop = FALSE] %*% beta_k)

        res_train <- Yt[idx_train] - eta_train
        rss_train <- sum(res_train^2)

        df_k <- sum(abs(beta_k) > 1e-6)

        if (denom == "n-df") {
          den_k <- max(length(idx_train) - df_k, 1L)
        } else {
          den_k <- length(idx_train)
        }

        sigma2_train <- rss_train / den_k
      }

      if (!is.finite(sigma2_train)) {
        sigma2_train <- mean(Yt[idx_train]^2)
      }

      sigma2_train <- min(max(sigma2_train, sigma_floor), sigma_cap)

      sigma2_oof[idx_test] <- sigma2_train
      wloss[idx_test] <- 1 / sigma2_train

      beta_fold[[k]][[ff]] <- beta_k
      sigma2_fold[[k]][ff] <- sigma2_train
    }
  }

  # Fallback for any missing value.
  if (any(!is.finite(wloss))) {
    finite_w <- wloss[is.finite(wloss)]
    fallback_w <- ifelse(length(finite_w) > 0, stats::median(finite_w), 1)
    wloss[!is.finite(wloss)] <- fallback_w
  }

  wloss <- pmin(pmax(wloss, w_floor), w_cap)

  if (isTRUE(normalize)) {
    wloss <- wloss / mean(wloss)
  }

  # Group-level summary only for diagnostics.
  sigma2_group_mean <- rep(NA_real_, K)
  wloss_group_mean <- rep(NA_real_, K)

  for (k in seq_len(K)) {
    idx_k <- which(G == k)

    if (length(idx_k) > 0) {
      sigma2_group_mean[k] <- mean(sigma2_oof[idx_k], na.rm = TRUE)
      wloss_group_mean[k] <- mean(wloss[idx_k], na.rm = TRUE)
    }
  }

  list(
    wloss = as.numeric(wloss),
    sigma2_oof = as.numeric(sigma2_oof),
    sigma2_group_mean = as.numeric(sigma2_group_mean),
    wloss_group_mean = as.numeric(wloss_group_mean),
    foldid_group = foldid_group,
    beta_fold = beta_fold,
    sigma2_fold = sigma2_fold,
    nfold = nfold,
    denom = denom,
    normalize = normalize,
    sigma_floor = sigma_floor,
    sigma_cap = sigma_cap,
    w_floor = w_floor,
    w_cap = w_cap
  )
}