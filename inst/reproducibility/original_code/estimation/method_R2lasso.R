# ============================================================
# method_R2lasso.R
#
# R2Lasso:
# Separately penalized reparameterized homogeneous and
# heterogeneous parameters.
#
# Model:
#
#   beta^(k) = mu + alpha^(k),   k = 1,...,K
#
# subject to:
#
#   sum_{k=1}^K alpha^(k) = 0
#
# Objective:
#
#   1/(2n) * sum_i w_i *
#     (Ytilde_i - Atilde_i X_i^T beta_{G_i})^2
#
#   + lambda1 * sum_j w_j |mu_j|
#   + lambda2 * sum_j v_j ||alpha_j||_2
#
# Acceleration:
#   If beta_init[, j] is zero for all subgroups, variable j is
#   removed from the ADMM optimization and restored as zero in output.
#
# This file assumes utils.R has already been sourced and provides:
#
#   build_block_design()
#   weighted_xy()
#   bic_lm_loss()
#   vec_beta()
#   mat_beta()
# ============================================================

# ------------------------------------------------------------
# Fallback helpers for zero-initialization screening
# ------------------------------------------------------------
if (!exists("screen_zero_init_vars", mode = "function")) {
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
}

if (!exists("restore_beta_screen", mode = "function")) {
  restore_beta_screen <- function(beta_active, active_j, K) {
    beta_active <- as.matrix(beta_active)
    beta_full <- matrix(0, nrow = K, ncol = length(active_j))

    if (sum(active_j) > 0) {
      beta_full[, active_j] <- beta_active
    }

    beta_full
  }
}

if (!exists("bic_zero_fit", mode = "function")) {
  bic_zero_fit <- function(mydat) {
    bic_lm_loss(
      y = mydat$Ytilde,
      eta = rep(0, mydat$n),
      df = 0,
      n = mydat$n
    )
  }
}

# ------------------------------------------------------------
# Soft-thresholding
# ------------------------------------------------------------
if (!exists("soft", mode = "function")) {
  soft <- function(x, lambda) {
    sign(x) * pmax(abs(x) - lambda, 0)
  }
}

# ------------------------------------------------------------
# Group soft-thresholding
# ------------------------------------------------------------
if (!exists("group_soft", mode = "function")) {
  group_soft <- function(x, lambda) {
    nx <- sqrt(sum(x^2))

    if (!is.finite(nx) || nx <= lambda) {
      return(rep(0, length(x)))
    }

    (1 - lambda / nx) * x
  }
}

# ------------------------------------------------------------
# Make sum-to-zero constraint matrix for alpha
# ------------------------------------------------------------
make_alpha_sum_constraint_R2 <- function(K, p) {
  S <- matrix(0, nrow = p, ncol = K * p)

  for (j in seq_len(p)) {
    idx <- j + (0:(K - 1)) * p
    S[j, idx] <- 1
  }

  S
}

# ------------------------------------------------------------
# Convert mu and alpha matrix to beta matrix
# ------------------------------------------------------------
recover_beta_R2 <- function(mu, Alpha) {
  K <- nrow(Alpha)
  p <- ncol(Alpha)

  matrix(rep(mu, each = K), nrow = K, ncol = p) + Alpha
}

# ------------------------------------------------------------
# Degrees of freedom for R2Lasso
# ------------------------------------------------------------
df_R2lasso <- function(mu, Alpha, tol = 1e-6) {
  K <- nrow(Alpha)
  p <- ncol(Alpha)

  df_mu <- sum(abs(mu) > tol)

  alpha_norm <- rep(0, p)

  for (j in seq_len(p)) {
    alpha_norm[j] <- sqrt(sum(Alpha[, j]^2))
  }

  df_alpha <- sum(alpha_norm > tol) * (K - 1)

  df_mu + df_alpha
}

# ------------------------------------------------------------
# Lambda1 maximum for mu penalty
# ------------------------------------------------------------
lambda1_max_R2lasso <- function(g_mu, wj, eps = 1e-8) {
  if (length(g_mu) != length(wj)) {
    stop(
      "Length mismatch in lambda1_max_R2lasso(). length(g_mu) = ",
      length(g_mu), ", length(wj) = ", length(wj), ".",
      call. = FALSE
    )
  }

  lam <- max(abs(g_mu) / pmax(wj, eps), na.rm = TRUE)

  if (!is.finite(lam) || lam <= 0) {
    lam <- 1
  }

  lam
}

# ------------------------------------------------------------
# Lambda2 maximum for alpha group penalty
# ------------------------------------------------------------
lambda2_max_R2lasso <- function(g_alpha, vj, K, p, eps = 1e-8) {
  if (length(g_alpha) != K * p) {
    stop(
      "Length mismatch in lambda2_max_R2lasso(). Expected ",
      K * p, ", got ", length(g_alpha), ".",
      call. = FALSE
    )
  }

  if (length(vj) != p) {
    stop(
      "Length of vj must equal p. Expected ",
      p, ", got ", length(vj), ".",
      call. = FALSE
    )
  }

  score <- rep(0, p)

  for (j in seq_len(p)) {
    idx <- j + (0:(K - 1)) * p

    gj <- g_alpha[idx]
    gj <- gj - mean(gj)

    score[j] <- sqrt(sum(gj^2)) / pmax(vj[j], eps)
  }

  lam <- max(score, na.rm = TRUE)

  if (!is.finite(lam) || lam <= 0) {
    lam <- 1
  }

  lam
}

# ------------------------------------------------------------
# Single R2Lasso ADMM fit
# ------------------------------------------------------------
method_R2lasso_single <- function(mydat,
                                  lambda1,
                                  lambda2,
                                  wj,
                                  vj,
                                  wloss = NULL,
                                  beta_init = NULL,
                                  rho_mu = 1,
                                  rho_alpha = 1,
                                  rho_con = 1,
                                  maxit = 1000,
                                  eps_pri = 1e-4,
                                  eps_dual = 1e-4,
                                  df_tol = 1e-6,
                                  verbose = FALSE) {

  Yt <- mydat$Ytilde
  At <- mydat$Atilde
  X  <- mydat$X
  G  <- mydat$G
  K  <- mydat$K
  p  <- mydat$p
  n  <- mydat$n

  q_mu <- p
  q_alpha <- K * p
  q <- q_mu + q_alpha

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

  if (length(wj) != p) {
    stop(
      "Length of wj must equal p. Expected ",
      p, ", got ", length(wj), ".",
      call. = FALSE
    )
  }

  if (length(vj) != p) {
    stop(
      "Length of vj must equal p. Expected ",
      p, ", got ", length(vj), ".",
      call. = FALSE
    )
  }

  if (any(!is.finite(wj)) || any(wj <= 0)) {
    stop("wj contains non-finite or non-positive values.", call. = FALSE)
  }

  if (any(!is.finite(vj)) || any(vj <= 0)) {
    stop("vj contains non-finite or non-positive values.", call. = FALSE)
  }

  Z_mu <- X * as.numeric(At)

  Z_alpha <- build_block_design(
    X = X,
    Atilde = At,
    G = G,
    K = K
  )

  Z_full <- cbind(Z_mu, Z_alpha)

  if (ncol(Z_full) != q) {
    stop(
      "Unexpected R2Lasso design dimension. Expected ",
      q, ", got ", ncol(Z_full), ".",
      call. = FALSE
    )
  }

  WZ <- Z_full * as.numeric(wloss)

  H <- crossprod(Z_full, WZ) / n
  g <- drop(crossprod(Z_full, as.numeric(wloss) * Yt)) / n

  S <- make_alpha_sum_constraint_R2(K, p)

  if (!(nrow(S) == p && ncol(S) == q_alpha)) {
    stop("Unexpected dimension of alpha sum constraint matrix.", call. = FALSE)
  }

  if (is.null(beta_init)) {
    mu <- rep(0, p)
    Alpha <- matrix(0, nrow = K, ncol = p)
  } else {
    beta_init <- as.matrix(beta_init)

    if (!(nrow(beta_init) == K && ncol(beta_init) == p)) {
      if (nrow(beta_init) == p && ncol(beta_init) == K) {
        beta_init <- t(beta_init)
      } else {
        stop(
          "Dimension mismatch for beta_init in method_R2lasso_single(). ",
          "Expected K x p = ", K, " x ", p,
          " or p x K = ", p, " x ", K,
          ", got ", nrow(beta_init), " x ", ncol(beta_init), ".",
          call. = FALSE
        )
      }
    }

    mu <- colMeans(beta_init)

    Alpha <- beta_init -
      matrix(rep(mu, each = K), nrow = K, ncol = p)
  }

  alpha_vec <- vec_beta(Alpha)
  theta <- c(mu, alpha_vec)

  if (length(theta) != q) {
    stop(
      "Length of initial theta must be ",
      q, ", got ", length(theta), ".",
      call. = FALSE
    )
  }

  u <- mu
  z <- alpha_vec

  d_mu <- rep(0, p)
  d_alpha <- rep(0, q_alpha)
  d_con <- rep(0, p)

  M <- H

  M[seq_len(p), seq_len(p)] <-
    M[seq_len(p), seq_len(p)] + rho_mu * diag(p)

  alpha_idx <- (p + 1):q

  M[alpha_idx, alpha_idx] <-
    M[alpha_idx, alpha_idx] +
    rho_alpha * diag(q_alpha) +
    rho_con * crossprod(S)

  Rchol <- tryCatch(
    chol(M),
    error = function(e) {
      chol(M + 1e-8 * diag(q))
    }
  )

  it_used <- maxit
  r_pri <- NA_real_
  r_dual <- NA_real_

  for (it in seq_len(maxit)) {

    u_old <- u
    z_old <- z

    rhs <- g

    rhs[seq_len(p)] <- rhs[seq_len(p)] +
      rho_mu * (u - d_mu)

    rhs[alpha_idx] <- rhs[alpha_idx] +
      rho_alpha * (z - d_alpha)

    rhs[alpha_idx] <- rhs[alpha_idx] -
      rho_con * drop(t(S) %*% d_con)

    theta <- backsolve(
      Rchol,
      forwardsolve(t(Rchol), rhs)
    )

    mu <- theta[seq_len(p)]
    alpha_vec <- theta[alpha_idx]

    u <- soft(
      mu + d_mu,
      lambda1 * wj / rho_mu
    )

    z_new <- rep(0, q_alpha)
    tmp_alpha <- alpha_vec + d_alpha

    for (j in seq_len(p)) {
      idx <- j + (0:(K - 1)) * p

      z_new[idx] <- group_soft(
        tmp_alpha[idx],
        lambda2 * vj[j] / rho_alpha
      )
    }

    z <- z_new

    d_mu <- d_mu + mu - u
    d_alpha <- d_alpha + alpha_vec - z
    d_con <- d_con + drop(S %*% alpha_vec)

    r_pri <- sqrt(
      sum((mu - u)^2) +
        sum((alpha_vec - z)^2) +
        sum(drop(S %*% alpha_vec)^2)
    )

    r_dual <- sqrt(
      rho_mu^2 * sum((u - u_old)^2) +
        rho_alpha^2 * sum((z - z_old)^2)
    )

    if (isTRUE(verbose) && it %% 50 == 0) {
      message(
        "R2Lasso ADMM iter = ", it,
        ", r_pri = ", signif(r_pri, 4),
        ", r_dual = ", signif(r_dual, 4)
      )
    }

    if (r_pri < eps_pri && r_dual < eps_dual) {
      it_used <- it
      break
    }
  }

  mu_hat <- u

  Alpha_hat <- mat_beta(z, K, p)

  for (j in seq_len(p)) {
    Alpha_hat[, j] <- Alpha_hat[, j] - mean(Alpha_hat[, j])
  }

  beta_hat <- recover_beta_R2(
    mu = mu_hat,
    Alpha = Alpha_hat
  )

  eta <- drop(
    Z_mu %*% mu_hat +
      Z_alpha %*% vec_beta(Alpha_hat)
  )

  df <- df_R2lasso(
    mu = mu_hat,
    Alpha = Alpha_hat,
    tol = df_tol
  )

  bic <- bic_lm_loss(Yt, eta, df, n)

  list(
    beta = beta_hat,
    mu = mu_hat,
    alpha = Alpha_hat,
    lambda1 = lambda1,
    lambda2 = lambda2,
    bic = bic,
    df = df,
    iter = it_used,
    r_pri = r_pri,
    r_dual = r_dual
  )
}

# ------------------------------------------------------------
# Main R2Lasso tuning wrapper
# ------------------------------------------------------------
method_R2lasso <- function(mydat,
                           beta_init = NULL,
                           wloss = NULL,
                           gamma = 1,
                           nlambda1 = 10,
                           nlambda2 = 10,
                           lambda.min.ratio = 1e-4,
                           lambda1.max.multiplier = 10,
                           lambda2.max.multiplier = 10,
                           lambda1_seq = NULL,
                           lambda2_seq = NULL,
                           eps_weight = 1e-4,
                           rho_mu = 1,
                           rho_alpha = 1,
                           rho_con = 1,
                           maxit = 1000,
                           eps_pri = 1e-4,
                           eps_dual = 1e-4,
                           df_tol = 1e-6,
                           verbose = FALSE,
                           zero_init_tol = 0) {

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

  if (is.null(beta_init)) {
    beta_init <- method_glasso(mydat, wloss = wloss)$beta
  }

  beta_init <- as.matrix(beta_init)

  if (!(nrow(beta_init) == K && ncol(beta_init) == p_full)) {
    if (nrow(beta_init) == p_full && ncol(beta_init) == K) {
      beta_init <- t(beta_init)
    } else {
      stop(
        "Dimension mismatch for beta_init in method_R2lasso(). ",
        "Expected K x p = ", K, " x ", p_full,
        " or p x K = ", p_full, " x ", K,
        ", got ", nrow(beta_init), " x ", ncol(beta_init), ".",
        call. = FALSE
      )
    }
  }

  screen <- screen_zero_init_vars(
    mydat = mydat,
    beta_init = beta_init,
    zero_init_tol = zero_init_tol,
    caller = "method_R2lasso()"
  )

  active_j <- screen$active_j
  dropped_j <- screen$dropped_j
  p_active <- screen$p_active

  if (p_active == 0) {
    beta_zero <- matrix(0, nrow = K, ncol = p_full)
    mu_zero <- rep(0, p_full)
    alpha_zero <- matrix(0, nrow = K, ncol = p_full)
    bic0 <- bic_zero_fit(mydat)

    return(list(
      beta = beta_zero,
      beta_active = matrix(0, nrow = K, ncol = 0),
      mu = mu_zero,
      mu_active = numeric(0),
      alpha = alpha_zero,
      alpha_active = matrix(0, nrow = K, ncol = 0),
      lambda1 = NA_real_,
      lambda2 = NA_real_,
      bic = bic0,
      df = 0L,
      iter = 0L,
      r_pri = NA_real_,
      r_dual = NA_real_,
      wj = numeric(0),
      vj = numeric(0),
      wj_full = rep(Inf, p_full),
      vj_full = rep(Inf, p_full),
      lambda1_seq = NA_real_,
      lambda2_seq = NA_real_,
      bic_mat = matrix(bic0, nrow = 1, ncol = 1),
      df_mat = matrix(0L, nrow = 1, ncol = 1),
      iter_mat = matrix(0L, nrow = 1, ncol = 1),
      path = list(list(
        beta = beta_zero,
        beta_active = matrix(0, nrow = K, ncol = 0),
        mu = mu_zero,
        mu_active = numeric(0),
        alpha = alpha_zero,
        alpha_active = matrix(0, nrow = K, ncol = 0),
        lambda1 = NA_real_,
        lambda2 = NA_real_,
        bic = bic0,
        df = 0L,
        iter = 0L,
        r_pri = NA_real_,
        r_dual = NA_real_
      )),
      gamma = gamma,
      lambda.min.ratio = lambda.min.ratio,
      lambda1.max.multiplier = lambda1.max.multiplier,
      lambda2.max.multiplier = lambda2.max.multiplier,
      df_tol = df_tol,
      active_j = active_j,
      dropped_j = dropped_j,
      p_active = 0L,
      p_full = p_full,
      zero_init_tol = zero_init_tol
    ))
  }

  mydat_active <- screen$mydat
  beta_init_active <- screen$beta_init

  Yt <- mydat_active$Ytilde
  At <- mydat_active$Atilde
  X  <- mydat_active$X
  G  <- mydat_active$G
  p  <- mydat_active$p

  # ----------------------------------------------------------
  # Adaptive weights on active variables only
  # ----------------------------------------------------------
  mu_init <- colMeans(beta_init_active)

  Alpha_init <- beta_init_active -
    matrix(rep(mu_init, each = K), nrow = K, ncol = p)

  wj <- 1 / (abs(mu_init)^gamma + eps_weight)

  alpha_norm <- rep(0, p)

  for (j in seq_len(p)) {
    alpha_norm[j] <- sqrt(sum(Alpha_init[, j]^2))
  }

  vj <- 1 / (alpha_norm^gamma + eps_weight)

  if (any(!is.finite(wj)) || any(wj <= 0)) {
    stop("Non-finite or non-positive wj detected.", call. = FALSE)
  }

  if (any(!is.finite(vj)) || any(vj <= 0)) {
    stop("Non-finite or non-positive vj detected.", call. = FALSE)
  }

  # ----------------------------------------------------------
  # Build design and gradient scale for lambda grids
  # ----------------------------------------------------------
  Z_mu <- X * as.numeric(At)

  Z_alpha <- build_block_design(
    X = X,
    Atilde = At,
    G = G,
    K = K
  )

  Z_full <- cbind(Z_mu, Z_alpha)

  wy <- weighted_xy(Z_full, Yt, wloss)

  Zw <- wy$Z
  yw <- wy$y

  g <- drop(crossprod(Zw, yw)) / n

  g_mu <- g[seq_len(p)]
  g_alpha <- g[(p + 1):(p + K * p)]

  # ----------------------------------------------------------
  # Lambda1 grid
  # ----------------------------------------------------------
  if (is.null(lambda1_seq)) {

    lambda1_max0 <- lambda1_max_R2lasso(
      g_mu = g_mu,
      wj = wj
    )

    lambda1_max <- lambda1.max.multiplier * lambda1_max0
    lambda1_min <- lambda1_max * lambda.min.ratio

    lambda1_seq <- exp(seq(
      log(lambda1_max),
      log(lambda1_min),
      length.out = nlambda1
    ))

  } else {

    lambda1_seq <- as.numeric(lambda1_seq)
    lambda1_seq <- lambda1_seq[is.finite(lambda1_seq) & lambda1_seq >= 0]

    if (length(lambda1_seq) == 0) {
      stop("User-supplied lambda1_seq is empty or invalid.", call. = FALSE)
    }
  }

  # ----------------------------------------------------------
  # Lambda2 grid
  # ----------------------------------------------------------
  if (is.null(lambda2_seq)) {

    lambda2_max0 <- lambda2_max_R2lasso(
      g_alpha = g_alpha,
      vj = vj,
      K = K,
      p = p
    )

    lambda2_max <- lambda2.max.multiplier * lambda2_max0
    lambda2_min <- lambda2_max * lambda.min.ratio

    lambda2_seq <- exp(seq(
      log(lambda2_max),
      log(lambda2_min),
      length.out = nlambda2
    ))

  } else {

    lambda2_seq <- as.numeric(lambda2_seq)
    lambda2_seq <- lambda2_seq[is.finite(lambda2_seq) & lambda2_seq >= 0]

    if (length(lambda2_seq) == 0) {
      stop("User-supplied lambda2_seq is empty or invalid.", call. = FALSE)
    }
  }

  lambda1_seq <- sort(unique(lambda1_seq), decreasing = TRUE)
  lambda2_seq <- sort(unique(lambda2_seq), decreasing = TRUE)

  # ----------------------------------------------------------
  # 2D grid search by BIC
  # ----------------------------------------------------------
  nfit <- length(lambda1_seq) * length(lambda2_seq)

  path <- vector("list", nfit)

  best <- NULL
  best_bic <- Inf

  counter <- 0L

  for (l1 in lambda1_seq) {
    for (l2 in lambda2_seq) {

      counter <- counter + 1L

      fit <- method_R2lasso_single(
        mydat = mydat_active,
        lambda1 = l1,
        lambda2 = l2,
        wj = wj,
        vj = vj,
        wloss = wloss,
        beta_init = beta_init_active,
        rho_mu = rho_mu,
        rho_alpha = rho_alpha,
        rho_con = rho_con,
        maxit = maxit,
        eps_pri = eps_pri,
        eps_dual = eps_dual,
        df_tol = df_tol,
        verbose = FALSE
      )

      fit$beta_active <- fit$beta
      fit$mu_active <- fit$mu
      fit$alpha_active <- fit$alpha

      fit$beta <- restore_beta_screen(
        beta_active = fit$beta_active,
        active_j = active_j,
        K = K
      )

      mu_full <- rep(0, p_full)
      mu_full[active_j] <- fit$mu_active
      fit$mu <- mu_full

      fit$alpha <- restore_beta_screen(
        beta_active = fit$alpha_active,
        active_j = active_j,
        K = K
      )

      fit$active_j <- active_j
      fit$dropped_j <- dropped_j
      fit$p_active <- p_active
      fit$p_full <- p_full

      path[[counter]] <- fit

      if (isTRUE(verbose)) {
        message(
          "[", counter, "/", nfit, "] ",
          "lambda1 = ", signif(l1, 4),
          ", lambda2 = ", signif(l2, 4),
          ", df = ", fit$df,
          ", bic = ", signif(fit$bic, 5),
          ", iter = ", fit$iter,
          ", r_pri = ", signif(fit$r_pri, 3),
          ", r_dual = ", signif(fit$r_dual, 3)
        )
      }

      if (fit$bic < best_bic) {
        best_bic <- fit$bic
        best <- fit
      }
    }
  }

  bic_mat <- matrix(
    sapply(path, function(x) x$bic),
    nrow = length(lambda1_seq),
    ncol = length(lambda2_seq),
    byrow = TRUE
  )

  df_mat <- matrix(
    sapply(path, function(x) x$df),
    nrow = length(lambda1_seq),
    ncol = length(lambda2_seq),
    byrow = TRUE
  )

  iter_mat <- matrix(
    sapply(path, function(x) x$iter),
    nrow = length(lambda1_seq),
    ncol = length(lambda2_seq),
    byrow = TRUE
  )

  wj_full <- rep(Inf, p_full)
  wj_full[active_j] <- wj

  vj_full <- rep(Inf, p_full)
  vj_full[active_j] <- vj

  best$wj <- wj
  best$vj <- vj
  best$wj_full <- wj_full
  best$vj_full <- vj_full

  best$lambda1_seq <- lambda1_seq
  best$lambda2_seq <- lambda2_seq

  best$bic_mat <- bic_mat
  best$df_mat <- df_mat
  best$iter_mat <- iter_mat

  best$path <- path

  best$gamma <- gamma
  best$lambda.min.ratio <- lambda.min.ratio
  best$lambda1.max.multiplier <- lambda1.max.multiplier
  best$lambda2.max.multiplier <- lambda2.max.multiplier
  best$df_tol <- df_tol

  best$active_j <- active_j
  best$dropped_j <- dropped_j
  best$p_active <- p_active
  best$p_full <- p_full
  best$zero_init_tol <- zero_init_tol

  best
}