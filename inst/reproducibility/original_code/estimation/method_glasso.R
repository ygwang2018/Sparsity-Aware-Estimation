# ============================================================
# method_glasso.R
#
# Group lasso for subgroup-specific treatment effects.
#
# Model:
#   Ytilde_i approx Atilde_i * X_i^T beta_{G_i}
#
# Parameter order in design matrix Z:
#   beta_1^(1), beta_1^(2), ..., beta_1^(K),
#   beta_2^(1), beta_2^(2), ..., beta_2^(K),
#   ...
#   beta_p^(1), beta_p^(2), ..., beta_p^(K)
#
# No intercept is fitted.
#
# Acceleration:
#   If beta_init[, j] is zero for all subgroups, variable j is
#   removed from the optimization and restored as zero in output.
# ============================================================

# ------------------------------------------------------------
# Build design matrix with variable-major order
#
# For each subject i in subgroup G_i = k:
#   columns corresponding to variable j and subgroup k receive
#   Atilde_i * X_ij
#
# Z has dimension:
#   n x (p * K)
# ------------------------------------------------------------
build_var_order_design <- function(X, Atilde, G, K) {
  n <- nrow(X)
  p <- ncol(X)

  G <- as.integer(G)
  Atilde <- as.numeric(Atilde)

  if (length(G) != n) {
    stop("Length of G must equal nrow(X).", call. = FALSE)
  }

  if (length(Atilde) != n) {
    stop("Length of Atilde must equal nrow(X).", call. = FALSE)
  }

  if (any(G < 1 | G > K)) {
    stop("G must contain subgroup labels from 1 to K.", call. = FALSE)
  }

  XA <- X * Atilde

  Z <- matrix(0, nrow = n, ncol = p * K)

  # Column index for variable j and subgroup k:
  #   colid = (j - 1) * K + k
  for (j in seq_len(p)) {
    for (k in seq_len(K)) {
      idx <- which(G == k)
      colid <- (j - 1) * K + k
      Z[idx, colid] <- XA[idx, j]
    }
  }

  Z
}

# ------------------------------------------------------------
# Convert vector in variable-major order to K x p beta matrix
#
# Input b order:
#   beta_1^(1), ..., beta_1^(K),
#   beta_2^(1), ..., beta_2^(K),
#   ...
#   beta_p^(1), ..., beta_p^(K)
#
# Output:
#   beta[k, j] = coefficient of variable j in subgroup k
# ------------------------------------------------------------
var_order_to_beta <- function(b, K, p) {
  b <- as.numeric(b)

  if (length(b) != K * p) {
    stop(
      "Length mismatch in var_order_to_beta(). Expected ",
      K * p, ", got ", length(b), ".",
      call. = FALSE
    )
  }

  B <- matrix(0, nrow = K, ncol = p)

  for (j in seq_len(p)) {
    B[, j] <- b[((j - 1) * K + 1):(j * K)]
  }

  B
}

# ------------------------------------------------------------
# Convert K x p beta matrix to variable-major vector
# ------------------------------------------------------------
beta_to_var_order <- function(B) {
  B <- as.matrix(B)

  K <- nrow(B)
  p <- ncol(B)

  b <- rep(0, K * p)

  for (j in seq_len(p)) {
    b[((j - 1) * K + 1):(j * K)] <- B[, j]
  }

  b
}

# ------------------------------------------------------------
# Safe coefficient extractor for grpreg without intercept
#
# Since we fit with intercept = FALSE, expected length is pZ.
#
# However, some versions or methods may still return an intercept
# from coef(). If length is pZ + 1, we remove the first element.
# ------------------------------------------------------------
get_grpreg_coef_no_intercept <- function(fit, lambda, pZ) {
  b <- as.numeric(coef(fit, lambda = lambda))

  if (length(b) == pZ) {
    return(b)
  }

  if (length(b) == pZ + 1) {
    return(b[-1])
  }

  stop(
    "Unexpected coefficient length from grpreg. Expected ",
    pZ, " or ", pZ + 1,
    ", got ", length(b), ".",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Main group lasso method
# ------------------------------------------------------------
method_glasso <- function(mydat,
                          beta_init = NULL,
                          wloss = NULL,
                          gamma = 1,
                          nlambda = 50,
                          lambda.min.ratio = 0.02,
                          zero_init_tol = 0) {

  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop(
      "Package 'grpreg' is required. Please install it first:\n",
      "install.packages('grpreg')",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Extract full data
  # ----------------------------------------------------------
  Yt_full <- mydat$Ytilde
  At_full <- mydat$Atilde
  X_full <- mydat$X
  G_full <- mydat$G
  K <- mydat$K
  p_full <- mydat$p
  n <- mydat$n

  if (is.null(wloss)) {
    wloss <- rep(1, n)
  }

  if (length(wloss) != n) {
    stop("Length of wloss must equal mydat$n.", call. = FALSE)
  }

  # ----------------------------------------------------------
  # Initial beta for adaptive group weights
  # ----------------------------------------------------------
  if (is.null(beta_init)) {
    beta_init <- method_loc_lasso(mydat, wloss = wloss)$beta
  }

  beta_init <- as.matrix(beta_init)

  if (!(nrow(beta_init) == K && ncol(beta_init) == p_full)) {
    if (nrow(beta_init) == p_full && ncol(beta_init) == K) {
      beta_init <- t(beta_init)
    } else {
      stop(
        "Dimension mismatch for beta_init in method_glasso(). ",
        "Expected K x p = ", K, " x ", p_full,
        " or p x K = ", p_full, " x ", K,
        ", but got ",
        nrow(beta_init), " x ", ncol(beta_init), ".",
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
    caller = "method_glasso()"
  )

  active_j <- screen$active_j
  dropped_j <- screen$dropped_j
  p_active <- screen$p_active

  # If all variables are screened out, return a full zero estimate.
  if (p_active == 0) {
    beta_zero <- matrix(0, nrow = K, ncol = p_full)

    return(list(
      beta = beta_zero,
      beta_active = matrix(0, nrow = K, ncol = 0),
      lambda = NA_real_,
      bic = bic_zero_fit(mydat),
      wj = numeric(0),
      wj_full = rep(Inf, p_full),
      fit = NULL,
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
  X <- mydat_active$X
  G <- mydat_active$G
  p <- mydat_active$p

  # ----------------------------------------------------------
  # Build design matrix only for active variables
  # ----------------------------------------------------------
  Z <- build_var_order_design(
    X = X,
    Atilde = At,
    G = G,
    K = K
  )

  # Weighted design and response
  wy <- weighted_xy(Z, Yt, wloss)
  Zw <- wy$Z
  yw <- wy$y

  # group structure:
  # each original active variable j forms one group of size K
  group <- rep(seq_len(p), each = K)

  # ----------------------------------------------------------
  # Adaptive group weights
  # ----------------------------------------------------------
  wj <- rep(0, p)

  for (j in seq_len(p)) {
    wj[j] <- 1 / (sqrt(sum(beta_init_active[, j]^2))^gamma + 1e-4)
  }

  group.multiplier <- wj

  # ----------------------------------------------------------
  # Fit group lasso WITHOUT intercept
  # ----------------------------------------------------------
  fit <- grpreg::grpreg(
    X = Zw,
    y = yw,
    group = group,
    penalty = "grLasso",
    group.multiplier = group.multiplier,
    nlambda = nlambda,
    lambda.min = lambda.min.ratio,
    intercept = FALSE
  )

  pZ <- ncol(Z)

  # ----------------------------------------------------------
  # Select lambda by BIC
  # ----------------------------------------------------------
  BIC <- rep(NA_real_, length(fit$lambda))

  for (ll in seq_along(fit$lambda)) {

    b <- get_grpreg_coef_no_intercept(
      fit = fit,
      lambda = fit$lambda[ll],
      pZ = pZ
    )

    eta <- drop(Z %*% b)

    bh <- var_order_to_beta(
      b = b,
      K = K,
      p = p
    )

    # df: selected variables times K
    # A variable is selected if any subgroup coefficient is nonzero.
    selected_j <- apply(abs(bh) > 1e-6, 2, any)
    df <- sum(selected_j) * K

    BIC[ll] <- bic_lm_loss(Yt, eta, df, n)
  }

  id <- which.min(BIC)

  # ----------------------------------------------------------
  # Final coefficient at selected lambda
  # ----------------------------------------------------------
  bhat <- get_grpreg_coef_no_intercept(
    fit = fit,
    lambda = fit$lambda[id],
    pZ = pZ
  )

  beta_active <- var_order_to_beta(
    b = bhat,
    K = K,
    p = p
  )

  beta_hat <- restore_beta_screen(
    beta_active = beta_active,
    active_j = active_j,
    K = K
  )

  wj_full <- rep(Inf, p_full)
  wj_full[active_j] <- wj

  list(
    beta = beta_hat,
    beta_active = beta_active,
    lambda = fit$lambda[id],
    bic = BIC[id],
    wj = wj,
    wj_full = wj_full,
    fit = fit,
    active_j = active_j,
    dropped_j = dropped_j,
    p_active = p_active,
    p_full = p_full,
    zero_init_tol = zero_init_tol
  )
}