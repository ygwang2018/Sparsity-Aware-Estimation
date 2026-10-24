soft <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

group_soft <- function(x, lambda) {
  nx <- sqrt(sum(x^2))
  if (nx == 0) return(x * 0)
  pmax(1 - lambda / nx, 0) * x
}

make_AR_cov <- function(p, rho) {
  M <- matrix(rep(1:p, p), ncol = p, byrow = FALSE)
  rho ^ abs(M - t(M))
}

vec_beta <- function(B) {
  # B: K x p, stacked by subgroup: beta^(1),...,beta^(K)
  as.vector(t(B))
}

mat_beta <- function(b, K, p) {
  matrix(b, nrow = K, ncol = p, byrow = TRUE)
}

center_beta_j <- function(bj) {
  bj - mean(bj)
}

get_active_set <- function(beta, tol = 1e-6) {
  which(abs(beta) > tol, arr.ind = TRUE)
}

get_homo_set <- function(beta, tol = 1e-6) {
  # beta: K x p
  K <- nrow(beta)
  p <- ncol(beta)
  out <- rep(FALSE, p)
  for (j in 1:p) {
    out[j] <- sqrt(sum((beta[, j] - mean(beta[, j]))^2)) <= tol
  }
  out
}

eval_methods <- function(beta_hat_list, beta0, tol = 1e-6) {
  K <- nrow(beta0)
  p <- ncol(beta0)

  true_active <- abs(beta0) > tol
  true_active_var <- apply(true_active, 2, any)

  true_homo <- get_homo_set(beta0, tol)

  # Nonzero homogeneous effects only.
  true_homo_active <- true_homo & true_active_var

  out <- data.frame()

  for (nm in names(beta_hat_list)) {
    bh <- as.matrix(beta_hat_list[[nm]])

    if (!(nrow(bh) == K && ncol(bh) == p)) {
      stop(
        "Dimension mismatch for method ", nm,
        ". Expected ", K, " x ", p,
        ", got ", nrow(bh), " x ", ncol(bh), ".",
        call. = FALSE
      )
    }

    est_active <- abs(bh) > tol
    est_active_var <- apply(est_active, 2, any)

    est_homo <- get_homo_set(bh, tol)
    est_homo_active <- est_homo & est_active_var

    mse <- mean((bh - beta0)^2)

    # Variable selection at coefficient level.
    var_tpr <- ifelse(
      sum(true_active) == 0,
      NA_real_,
      sum(est_active & true_active) / sum(true_active)
    )

    var_fdr <- ifelse(
      sum(est_active) == 0,
      0,
      sum(est_active & !true_active) / sum(est_active)
    )

    # Homogeneous effect selection at variable level.
    # Positive discovery: estimated active and homogeneous.
    homo_tpr <- ifelse(
      sum(true_homo_active) == 0,
      NA_real_,
      sum(est_homo_active & true_homo_active) / sum(true_homo_active)
    )

    homo_fdr <- ifelse(
      sum(est_homo_active) == 0,
      0,
      sum(est_homo_active & !true_homo_active) / sum(est_homo_active)
    )

    out <- rbind(out, data.frame(
      method = nm,
      MSE = mse*p,
      Var_TPR = var_tpr,
      Var_FDR = var_fdr,
      Homo_TPR = homo_tpr,
      Homo_FDR = homo_fdr
    ))
  }

  out
}

bic_lm_loss <- function(y, eta, df, n) {
  rss <- sum((y - eta)^2)
  log(rss / n + 1e-12) + log(n) * df / n
}

build_block_design <- function(X, Atilde, G, K) {
  # Output Z: n x (Kp), ordered by subgroup: beta1, beta2,...
  n <- nrow(X)
  p <- ncol(X)
  Z <- matrix(0, n, K * p)

  XA <- X * as.numeric(Atilde)

  for (k in 1:K) {
    idx <- which(G == k)
    Z[idx, ((k - 1) * p + 1):(k * p)] <- XA[idx, , drop = FALSE]
  }

  Z
}

weighted_xy <- function(Z, y, wloss) {
  sw <- sqrt(wloss)
  list(
    Z = Z * sw,
    y = y * sw
  )
}

make_lambda_grid <- function(Z, y, penalty.factor = NULL, nlambda = 40, lambda.min.ratio = 0.02) {
  n <- nrow(Z)
  if (is.null(penalty.factor)) penalty.factor <- rep(1, ncol(Z))
  score <- abs(drop(t(Z) %*% y)) / n
  score <- score / pmax(penalty.factor, 1e-8)
  lam.max <- max(score, na.rm = TRUE)
  exp(seq(log(lam.max), log(lam.max * lambda.min.ratio), length.out = nlambda))
}


# ============================================================
# Convert beta_hat_list to wide format
#
# Input:
#   beta_hat_list: named list of beta matrices
#                  each beta usually has dimension K x p
#                  rows = subgroup, columns = covariates
#
# Output:
#   one row per method
#   columns:
#     method, simu, seed, beta_g1_x1, beta_g1_x2, ...
# ============================================================

beta_hat_to_wide <- function(beta_hat_list, simu = NA_integer_, seed = NA_integer_) {
  if (is.null(beta_hat_list) || !is.list(beta_hat_list)) {
    stop("beta_hat_list must be a non-null list.", call. = FALSE)
  }

  if (is.null(names(beta_hat_list))) {
    names(beta_hat_list) <- paste0("method", seq_along(beta_hat_list))
  }

  OUT <- NULL

  for (mm in seq_along(beta_hat_list)) {
    method_name <- names(beta_hat_list)[mm]

    beta <- beta_hat_list[[mm]]

    if (is.null(beta)) {
      warning("beta_hat for method ", method_name, " is NULL. Skip.")
      next
    }

    beta <- as.matrix(beta)

    if (any(!is.finite(beta))) {
      warning("beta_hat for method ", method_name, " contains non-finite values.")
    }

    K <- nrow(beta)
    p <- ncol(beta)

    # column names: beta_g1_x1, beta_g1_x2, ..., beta_gK_xp
    beta_names <- as.vector(
      outer(
        paste0("g", seq_len(K)),
        paste0("x", seq_len(p)),
        FUN = function(g, x) paste0("beta_", g, "_", x)
      )
    )

    # as.vector(t(beta)) gives row-wise order:
    # beta[1,1], beta[1,2], ..., beta[1,p], beta[2,1], ...
    beta_vec <- as.vector(t(beta))

    one_row <- data.frame(
      simu = simu,
      seed = seed,
      method = method_name,
      stringsAsFactors = FALSE
    )

    beta_df <- as.data.frame(
      as.list(beta_vec),
      stringsAsFactors = FALSE
    )

    colnames(beta_df) <- beta_names

    one_row <- cbind(one_row, beta_df)

    OUT <- rbind(OUT, one_row)
  }

  rownames(OUT) <- NULL

  OUT
}





# ------------------------------------------------------------
# Post-selection refit after local lasso, then estimate
# subgroup-specific variance and construct wloss.
#
# Input:
#   mydat      : prepared data object
#   beta_local : K x p local lasso beta estimate
#
# Main idea:
#   For each subgroup k:
#     active_k = {j: |beta_local[k,j]| > select_tol}
#     Refit Ytilde ~ Atilde * X_active_k without intercept
#     sigma2_k = RSS_k / denominator
#     w_k = 1 / sigma2_k
#
# Output:
#   list(
#     wloss       : individual-level weights, length n, mean normalized to 1
#     wloss_group : group-level weights, length K
#     sigma2      : group-level variance estimates, length K
#     beta_refit  : K x p refitted beta matrix
#     active      : list of selected variables per group
#     df          : refit rank per group
#   )
# ------------------------------------------------------------
make_wloss_from_local_refit <- function(mydat,
                                        beta_local,
                                        select_tol = 1e-6,
                                        sigma_floor = 1e-6,
                                        sigma_cap = Inf,
                                        denom = c("n-df", "n"),
                                        normalize = TRUE,
                                        selection = c("group", "union"),
                                        verbose = FALSE) {

  denom <- match.arg(denom)
  selection <- match.arg(selection)

  Yt <- as.numeric(mydat$Ytilde)
  At <- as.numeric(mydat$Atilde)
  X  <- as.matrix(mydat$X)
  G  <- as.integer(mydat$G)
  K  <- mydat$K
  p  <- mydat$p
  n  <- mydat$n

  beta_local <- as.matrix(beta_local)

  # Defensive dimension check
  if (!(nrow(beta_local) == K && ncol(beta_local) == p)) {
    if (nrow(beta_local) == p && ncol(beta_local) == K) {
      beta_local <- t(beta_local)
    } else {
      stop(
        "Dimension mismatch for beta_local. Expected K x p = ",
        K, " x ", p, " or p x K = ", p, " x ", K,
        ", got ", nrow(beta_local), " x ", ncol(beta_local), ".",
        call. = FALSE
      )
    }
  }

  beta_refit <- matrix(0, nrow = K, ncol = p)
  sigma2 <- rep(NA_real_, K)
  df <- rep(0L, K)
  active <- vector("list", K)

  # If selection = "union", use variables selected in any subgroup
  if (selection == "union") {
    active_union <- which(apply(abs(beta_local) > select_tol, 2, any))
  }

  for (k in seq_len(K)) {

    idx <- which(G == k)
    nk <- length(idx)

    if (nk == 0) {
      sigma2[k] <- NA_real_
      active[[k]] <- integer(0)
      next
    }

    if (selection == "group") {
      act_k <- which(abs(beta_local[k, ]) > select_tol)
    } else {
      act_k <- active_union
    }

    active[[k]] <- act_k

    # --------------------------------------------------------
    # Case 1: no selected variables in group k
    # Model fitted value is zero because Ytilde is residualized.
    # --------------------------------------------------------
    if (length(act_k) == 0) {
      res <- Yt[idx]
      rss <- sum(res^2)

      df[k] <- 0L

      if (denom == "n-df") {
        den_k <- max(nk - df[k], 1L)
      } else {
        den_k <- nk
      }

      sigma2[k] <- rss / den_k

      if (isTRUE(verbose)) {
        message("Group ", k, ": no active variable; sigma2 = ", signif(sigma2[k], 4))
      }

      next
    }

    # --------------------------------------------------------
    # Case 2: refit selected variables by least squares
    #
    # Local subgroup model:
    #   Ytilde_i = Atilde_i * X_i,act^T beta_k,act + error
    #
    # No intercept.
    # --------------------------------------------------------
    Zk <- X[idx, act_k, drop = FALSE] * At[idx]

    fit_k <- stats::lm.fit(x = Zk, y = Yt[idx])

    coef_k <- as.numeric(fit_k$coefficients)
    coef_k[is.na(coef_k)] <- 0

    beta_refit[k, act_k] <- coef_k

    eta_k <- drop(Zk %*% coef_k)
    res_k <- Yt[idx] - eta_k
    rss_k <- sum(res_k^2)

    # lm.fit rank handles collinearity
    rank_k <- fit_k$rank
    df[k] <- rank_k

    if (denom == "n-df") {
      den_k <- max(nk - rank_k, 1L)
    } else {
      den_k <- nk
    }

    sigma2[k] <- rss_k / den_k

    if (isTRUE(verbose)) {
      message(
        "Group ", k,
        ": active = ", length(act_k),
        ", rank = ", rank_k,
        ", sigma2 = ", signif(sigma2[k], 4)
      )
    }
  }

  # ----------------------------------------------------------
  # Stabilize sigma2 and convert to weights
  # ----------------------------------------------------------
  sigma2 <- pmax(sigma2, sigma_floor)
  sigma2 <- pmin(sigma2, sigma_cap)

  # If any NA due to empty group, replace by median finite sigma2
  if (any(!is.finite(sigma2))) {
    med_sig2 <- stats::median(sigma2[is.finite(sigma2)], na.rm = TRUE)
    if (!is.finite(med_sig2)) med_sig2 <- 1
    sigma2[!is.finite(sigma2)] <- med_sig2
  }

  wloss_group <- 1 / sigma2

  # Individual-level weights
  wloss <- wloss_group[G]

  # Normalize to mean 1, so lambda scale is comparable
  if (isTRUE(normalize)) {
    wloss <- wloss / mean(wloss)
    wloss_group <- wloss_group / mean(wloss_group[G])
  }

  list(
    wloss = as.numeric(wloss),
    wloss_group = as.numeric(wloss_group),
    sigma2 = as.numeric(sigma2),
    beta_refit = beta_refit,
    active = active,
    df = df,
    selection = selection,
    denom = denom
  )
}



# ------------------------------------------------------------
# Screen variables whose initial estimates are zero in all groups
# ------------------------------------------------------------
screen_zero_init_vars <- function(mydat, beta_init, zero_init_tol = 0,
                                  caller = "method") {
  K <- mydat$K
  p <- mydat$p

  beta_init <- as.matrix(beta_init)

  if (!(nrow(beta_init) == K && ncol(beta_init) == p)) {
    if (nrow(beta_init) == p && ncol(beta_init) == K) {
      beta_init <- t(beta_init)
    } else {
      stop(
        "Dimension mismatch for beta_init in ", caller, ". ",
        "Expected K x p = ", K, " x ", p,
        " or p x K = ", p, " x ", K,
        ", got ", nrow(beta_init), " x ", ncol(beta_init), ".",
        call. = FALSE
      )
    }
  }

  active_j <- colSums(abs(beta_init) > zero_init_tol) > 0

  mydat_active <- mydat
  mydat_active$X <- mydat$X[, active_j, drop = FALSE]
  mydat_active$p <- sum(active_j)

  list(
    mydat = mydat_active,
    beta_init = beta_init[, active_j, drop = FALSE],
    active_j = active_j,
    dropped_j = which(!active_j),
    original_p = p,
    p_active = sum(active_j)
  )
}

restore_beta_screen <- function(beta_active, active_j, K) {
  beta_active <- as.matrix(beta_active)
  beta_full <- matrix(0, nrow = K, ncol = length(active_j))

  if (sum(active_j) > 0) {
    beta_full[, active_j] <- beta_active
  }

  beta_full
}

restore_vector_screen <- function(x_active, active_j) {
  x_full <- rep(0, length(active_j))

  if (sum(active_j) > 0) {
    x_full[active_j] <- as.numeric(x_active)
  }

  x_full
}

bic_zero_fit <- function(mydat) {
  bic_lm_loss(
    y = mydat$Ytilde,
    eta = rep(0, mydat$n),
    df = 0,
    n = mydat$n
  )
}