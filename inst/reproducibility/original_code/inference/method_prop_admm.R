# ============================================================
# method_prop_admm.R
#
# Proposed ADMM method for subgroup-specific treatment effects.
#
# Objective:
#
#   1/(2n) * sum_i w_i (Ytilde_i - Z_i beta)^2
#   + lambda1 * sum_{k,j} W_kj |beta_kj|
#   + lambda2 * sum_j v_j || C beta_j ||_2
#
# where:
#
#   C = I_K - 1_K 1_K^T / K
#
# The second penalty encourages homogeneity across subgroups.
#
# BIC degrees of freedom:
#
#   df = number of unique nonzero values in beta matrix
#
# Acceleration:
#   If beta_init[, j] is zero for all subgroups, variable j is
#   removed from the ADMM optimization and restored as zero in output.
#
# This file assumes the same coefficient ordering as:
#
#   vec_beta()
#   mat_beta()
#   build_block_design()
#
# Specifically, it assumes:
#
#   vec_beta(B) = as.vector(t(B))
#
# for B being K x p.
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
# Construct Cstar
# ------------------------------------------------------------
make_Cstar <- function(K, p) {
  q <- K * p

  C <- diag(K) - matrix(1 / K, K, K)

  Cstar <- matrix(0, q, q)

  for (j in seq_len(p)) {
    idx <- j + (0:(K - 1)) * p
    Cstar[idx, idx] <- C
  }

  Cstar
}

# ------------------------------------------------------------
# Extract unique nonzero beta values
# ------------------------------------------------------------
get_unique_nonzero_beta_values <- function(Bhat, tol = 1e-6) {
  b <- as.numeric(Bhat)

  b_nz <- b[abs(b) > tol]

  if (length(b_nz) == 0) {
    return(numeric(0))
  }

  b_nz <- sort(b_nz)

  unique_vals <- b_nz[1]

  if (length(b_nz) >= 2) {
    for (ii in 2:length(b_nz)) {
      if (min(abs(b_nz[ii] - unique_vals)) > tol) {
        unique_vals <- c(unique_vals, b_nz[ii])
      }
    }
  }

  unique_vals
}

# ------------------------------------------------------------
# Degrees of freedom for BIC
# ------------------------------------------------------------
df_prop_admm <- function(Bhat, tol = 1e-6) {
  length(get_unique_nonzero_beta_values(Bhat, tol = tol))
}

# ------------------------------------------------------------
# Lambda1 maximum
# ------------------------------------------------------------
lambda1_max_prop_admm <- function(g, Wjk, eps = 1e-8) {
  w_vec <- vec_beta(Wjk)

  if (length(g) != length(w_vec)) {
    stop(
      "Length mismatch in lambda1_max_prop_admm(). length(g) = ",
      length(g), ", length(vec_beta(Wjk)) = ", length(w_vec), ".",
      call. = FALSE
    )
  }

  lam <- max(abs(g) / pmax(w_vec, eps))

  if (!is.finite(lam) || lam <= 0) {
    lam <- 1
  }

  lam
}

# ------------------------------------------------------------
# Lambda2 maximum
# ------------------------------------------------------------
lambda2_max_prop_admm <- function(g, Cstar, vj, K, p, eps = 1e-8) {
  if (length(vj) != p) {
    stop(
      "Length of vj must equal p. Expected ",
      p, ", got ", length(vj), ".",
      call. = FALSE
    )
  }

  Cg <- drop(Cstar %*% g)

  scores <- rep(0, p)

  for (j in seq_len(p)) {
    idx <- j + (0:(K - 1)) * p
    scores[j] <- sqrt(sum(Cg[idx]^2)) / pmax(vj[j], eps)
  }

  lam <- max(scores)

  if (!is.finite(lam) || lam <= 0) {
    lam <- 1
  }

  lam
}

# ------------------------------------------------------------
# Single ADMM fit
# ------------------------------------------------------------
method_prop_admm_single <- function(mydat,
                                    lambda1,
                                    lambda2,
                                    Wjk,
                                    vj,
                                    wloss = NULL,
                                    beta_init = NULL,
                                    rho1 = 1,
                                    rho2 = 1,
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

  q <- K * p

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

  if (!(nrow(Wjk) == K && ncol(Wjk) == p)) {
    stop(
      "Wjk must be K x p. Expected ",
      K, " x ", p, ", got ",
      nrow(Wjk), " x ", ncol(Wjk), ".",
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

  if (any(!is.finite(Wjk)) || any(Wjk <= 0)) {
    stop("Wjk contains non-finite or non-positive values.", call. = FALSE)
  }

  if (any(!is.finite(vj)) || any(vj <= 0)) {
    stop("vj contains non-finite or non-positive values.", call. = FALSE)
  }

  Z <- build_block_design(X, At, G, K)

  if (ncol(Z) != q) {
    stop(
      "Unexpected design dimension. Expected ncol(Z) = K * p = ",
      q, ", got ", ncol(Z), ".",
      call. = FALSE
    )
  }

  WZ <- Z * as.numeric(wloss)
  H <- crossprod(Z, WZ) / n
  g <- drop(crossprod(Z, as.numeric(wloss) * Yt)) / n

  Cstar <- make_Cstar(K, p)

  M <- H + rho1 * diag(q) + rho2 * crossprod(Cstar)

  Rchol <- tryCatch(
    chol(M),
    error = function(e) {
      chol(M + 1e-8 * diag(q))
    }
  )

  if (is.null(beta_init)) {
    beta <- rep(0, q)
  } else {
    beta_init <- as.matrix(beta_init)

    if (!(nrow(beta_init) == K && ncol(beta_init) == p)) {
      if (nrow(beta_init) == p && ncol(beta_init) == K) {
        beta_init <- t(beta_init)
      } else {
        stop(
          "Dimension mismatch for beta_init. Expected ",
          K, " x ", p, " or ", p, " x ", K,
          ", got ", nrow(beta_init), " x ", ncol(beta_init), ".",
          call. = FALSE
        )
      }
    }

    beta <- vec_beta(beta_init)
  }

  if (length(beta) != q) {
    stop(
      "Length of initial beta vector must be K * p = ",
      q, ", got ", length(beta), ".",
      call. = FALSE
    )
  }

  u <- beta
  z <- drop(Cstar %*% beta)

  d  <- rep(0, q)
  qq <- rep(0, q)

  w_vec <- vec_beta(Wjk)

  if (length(w_vec) != q) {
    stop(
      "Length mismatch for vec_beta(Wjk). Expected ",
      q, ", got ", length(w_vec), ".",
      call. = FALSE
    )
  }

  for (it in seq_len(maxit)) {
    u_old <- u
    z_old <- z

    a <- u - d
    b <- z - qq

    rhs <- g + rho1 * a + rho2 * drop(t(Cstar) %*% b)

    beta <- backsolve(
      Rchol,
      forwardsolve(t(Rchol), rhs)
    )

    u <- soft(beta + d, lambda1 * w_vec / rho1)

    cb <- drop(Cstar %*% beta)

    z_new <- rep(0, q)

    for (j in seq_len(p)) {
      idx <- j + (0:(K - 1)) * p
      tmp <- cb[idx] + qq[idx]
      z_new[idx] <- group_soft(tmp, lambda2 * vj[j] / rho2)
    }

    z <- z_new

    d  <- d  + beta - u
    qq <- qq + cb - z

    r_pri <- sqrt(
      sum((beta - u)^2) +
        sum((cb - z)^2)
    )

    r_dual <- sqrt(
      rho1^2 * sum((u - u_old)^2) +
        rho2^2 * sum(drop(t(Cstar) %*% (z - z_old))^2)
    )

    if (isTRUE(verbose) && it %% 50 == 0) {
      message(
        "ADMM iter = ", it,
        ", r_pri = ", signif(r_pri, 4),
        ", r_dual = ", signif(r_dual, 4)
      )
    }

    if (r_pri < eps_pri && r_dual < eps_dual) {
      break
    }
  }

  Bhat <- mat_beta(u, K, p)

  if (!(nrow(Bhat) == K && ncol(Bhat) == p)) {
    stop(
      "mat_beta() returned unexpected dimension. Expected ",
      K, " x ", p, ", got ",
      nrow(Bhat), " x ", ncol(Bhat), ".",
      call. = FALSE
    )
  }

  eta <- drop(Z %*% vec_beta(Bhat))

  unique_beta_values <- get_unique_nonzero_beta_values(
    Bhat,
    tol = df_tol
  )

  df <- length(unique_beta_values)

  bic <- bic_lm_loss(Yt, eta, df, n)

  list(
    beta = Bhat,
    iter = it,
    bic = bic,
    df = df,
    unique_beta_values = unique_beta_values,
    lambda1 = lambda1,
    lambda2 = lambda2,
    r_pri = r_pri,
    r_dual = r_dual
  )
}

# ------------------------------------------------------------
# Main tuning wrapper
# ------------------------------------------------------------
method_prop_admm <- function(mydat,
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
                             rho1 = 1,
                             rho2 = 1,
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
        "Dimension mismatch for beta_init in method_prop_admm(). ",
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
    caller = "method_prop_admm()"
  )

  active_j <- screen$active_j
  dropped_j <- screen$dropped_j
  p_active <- screen$p_active

  if (p_active == 0) {
    beta_zero <- matrix(0, nrow = K, ncol = p_full)
    bic0 <- bic_zero_fit(mydat)

    return(list(
      beta = beta_zero,
      beta_active = matrix(0, nrow = K, ncol = 0),
      iter = 0L,
      bic = bic0,
      df = 0L,
      unique_beta_values = numeric(0),
      lambda1 = NA_real_,
      lambda2 = NA_real_,
      r_pri = NA_real_,
      r_dual = NA_real_,
      Wjk = matrix(0, nrow = K, ncol = 0),
      vj = numeric(0),
      Wjk_full = matrix(Inf, nrow = K, ncol = p_full),
      vj_full = rep(Inf, p_full),
      lambda1_seq = NA_real_,
      lambda2_seq = NA_real_,
      bic_mat = matrix(bic0, nrow = 1, ncol = 1),
      df_mat = matrix(0L, nrow = 1, ncol = 1),
      iter_mat = matrix(0L, nrow = 1, ncol = 1),
      path = list(list(
        beta = beta_zero,
        beta_active = matrix(0, nrow = K, ncol = 0),
        iter = 0L,
        bic = bic0,
        df = 0L,
        unique_beta_values = numeric(0),
        lambda1 = NA_real_,
        lambda2 = NA_real_,
        r_pri = NA_real_,
        r_dual = NA_real_
      )),
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

  q <- K * p

  # ----------------------------------------------------------
  # Adaptive weights on active variables only
  # ----------------------------------------------------------
  Wjk <- 1 / (abs(beta_init_active)^gamma + eps_weight)

  Cbeta <- apply(
    beta_init_active,
    2,
    function(x) x - mean(x)
  )

  if (p == 1) {
    Cbeta <- matrix(Cbeta, nrow = K, ncol = 1)
  }

  Cnorm <- sqrt(colSums(Cbeta^2))

  vj <- 1 / (Cnorm^gamma + eps_weight)

  if (any(!is.finite(Wjk)) || any(Wjk <= 0)) {
    stop("Wjk contains non-finite or non-positive values.", call. = FALSE)
  }

  if (any(!is.finite(vj)) || any(vj <= 0)) {
    stop("vj contains non-finite or non-positive values.", call. = FALSE)
  }

  # ----------------------------------------------------------
  # Build design and compute gradient scale on active variables
  # ----------------------------------------------------------
  Z <- build_block_design(X, At, G, K)

  if (ncol(Z) != q) {
    stop(
      "Unexpected design dimension. Expected ncol(Z) = K * p = ",
      q, ", got ", ncol(Z), ".",
      call. = FALSE
    )
  }

  wy <- weighted_xy(Z, Yt, wloss)
  Zw <- wy$Z
  yw <- wy$y

  g <- drop(crossprod(Zw, yw)) / n

  Cstar <- make_Cstar(K, p)

  # ----------------------------------------------------------
  # Lambda1 grid
  # ----------------------------------------------------------
  if (is.null(lambda1_seq)) {
    lam1.max0 <- lambda1_max_prop_admm(
      g = g,
      Wjk = Wjk
    )

    lam1.max <- lambda1.max.multiplier * lam1.max0
    lam1.min <- lam1.max * lambda.min.ratio

    lambda1_seq <- exp(seq(
      log(lam1.max),
      log(lam1.min),
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
    lam2.max0 <- lambda2_max_prop_admm(
      g = g,
      Cstar = Cstar,
      vj = vj,
      K = K,
      p = p
    )

    lam2.max <- lambda2.max.multiplier * lam2.max0
    lam2.min <- lam2.max * lambda.min.ratio

    lambda2_seq <- exp(seq(
      log(lam2.max),
      log(lam2.min),
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
  # Grid search by BIC
  # ----------------------------------------------------------
  best <- NULL
  best_bic <- Inf

  nfit <- length(lambda1_seq) * length(lambda2_seq)

  path <- vector("list", nfit)

  counter <- 0L

  for (l1 in lambda1_seq) {
    for (l2 in lambda2_seq) {
      counter <- counter + 1L

      fit <- method_prop_admm_single(
        mydat = mydat_active,
        lambda1 = l1,
        lambda2 = l2,
        Wjk = Wjk,
        vj = vj,
        wloss = wloss,
        beta_init = beta_init_active,
        rho1 = rho1,
        rho2 = rho2,
        maxit = maxit,
        eps_pri = eps_pri,
        eps_dual = eps_dual,
        df_tol = df_tol,
        verbose = FALSE
      )

      fit$beta_active <- fit$beta
      fit$beta <- restore_beta_screen(
        beta_active = fit$beta_active,
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
          ", df_unique_beta = ", fit$df,
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

  # ----------------------------------------------------------
  # Store BIC and df matrices
  #
  # Rows correspond to lambda1_seq.
  # Columns correspond to lambda2_seq.
  # ----------------------------------------------------------
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

  Wjk_full <- matrix(Inf, nrow = K, ncol = p_full)
  Wjk_full[, active_j] <- Wjk

  vj_full <- rep(Inf, p_full)
  vj_full[active_j] <- vj

  best$Wjk <- Wjk
  best$vj <- vj
  best$Wjk_full <- Wjk_full
  best$vj_full <- vj_full

  best$lambda1_seq <- lambda1_seq
  best$lambda2_seq <- lambda2_seq

  best$bic_mat <- bic_mat
  best$df_mat <- df_mat
  best$iter_mat <- iter_mat

  best$path <- path

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