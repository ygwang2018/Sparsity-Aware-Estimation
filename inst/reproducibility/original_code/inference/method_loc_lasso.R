method_loc_lasso <- function(mydat, wloss = NULL,
                             gamma = 1,
                             nlambda = 50,
                             lambda.min.ratio = 0.02) {

  Yt <- mydat$Ytilde
  At <- mydat$Atilde
  X <- mydat$X
  G <- mydat$G
  K <- mydat$K
  p <- mydat$p
  n <- mydat$n

  if (is.null(wloss)) {
    wloss <- rep(1, n)
  }

  Z <- build_block_design(X, At, G, K)

  wy <- weighted_xy(Z, Yt, wloss)
  Zw <- wy$Z
  yw <- wy$y

  # Initial ridge for adaptive weights
  fit0 <- glmnet(
    Zw, yw,
    alpha = 0,
    intercept = FALSE,
    standardize = FALSE
  )

  # 原写法：
  # b0 <- as.numeric(coef(fit0, s = fit0$lambda[length(fit0$lambda)]))[-1]
  #
  # glmnet 对象中 fit0$beta 通常是 sparse matrix:
  # n_variables x n_lambda
  # 最后一列对应最小 lambda
  if (!is.null(fit0$beta)) {
    b0 <- as.numeric(fit0$beta[, ncol(fit0$beta)])
  } else {
    b0 <- as.numeric(coef(fit0, s = fit0$lambda[length(fit0$lambda)]))[-1]
  }

  pf <- 1 / (abs(b0)^gamma + 1e-4)

  fit <- glmnet(
    Zw, yw,
    alpha = 1,
    intercept = FALSE,
    standardize = FALSE,
    penalty.factor = pf,
    nlambda = nlambda,
    lambda.min.ratio = lambda.min.ratio
  )

  nlam <- length(fit$lambda)
  BIC <- rep(NA_real_, nlam)

  # 优先直接使用 fit$beta，避免循环中反复 coef()
  if (!is.null(fit$beta)) {

    Bmat <- as.matrix(fit$beta)

    # 防御性检查：glmnet 的 beta 应该是 ncol(Z) x nlambda
    if (nrow(Bmat) != ncol(Z) && ncol(Bmat) == ncol(Z)) {
      Bmat <- t(Bmat)
    }

    if (nrow(Bmat) != ncol(Z)) {
      stop("Unexpected dimension of fit$beta from glmnet.")
    }

    # 一次性计算所有 lambda 下的线性预测
    # Eta: n x nlambda
    Eta <- Z %*% Bmat

    # 一次性计算 df
    DF <- colSums(abs(Bmat) > 1e-6)

    for (ll in seq_len(nlam)) {
      BIC[ll] <- bic_lm_loss(Yt, Eta[, ll], DF[ll], n)
    }

    id <- which.min(BIC)
    bhat <- Bmat[, id]

  } else {

    # 备用方案：如果 glmnet 对象没有 fit$beta，则回退到 coef()
    for (ll in seq_len(nlam)) {
      b <- as.numeric(coef(fit, s = fit$lambda[ll]))[-1]

      eta <- drop(Z %*% b)
      df <- sum(abs(b) > 1e-6)

      BIC[ll] <- bic_lm_loss(Yt, eta, df, n)
    }

    id <- which.min(BIC)
    bhat <- as.numeric(coef(fit, s = fit$lambda[id]))[-1]
  }

  beta_hat <- mat_beta(bhat, K, p)

  list(
    beta = beta_hat,
    lambda = fit$lambda[id],
    bic = BIC[id],
    sigma2 = estimate_group_variance(mydat, beta_hat)
  )
}

estimate_group_variance <- function(mydat, beta_hat) {
  Yt <- mydat$Ytilde
  At <- mydat$Atilde
  X <- mydat$X
  G <- mydat$G
  K <- mydat$K
  n <- mydat$n

  # 原写法是对每个 k 循环：
  # for (k in 1:K) {
  #   idx <- which(G == k)
  #   eta <- At[idx] * drop(X[idx, , drop = FALSE] %*% beta_hat[k, ])
  #   res <- Yt[idx] - eta
  #   sig2[k] <- mean(res^2)
  # }
  #
  # 向量化写法：
  # 每个样本 i 使用 beta_hat[G[i], ] 作为对应 subgroup 的系数
  XB <- rowSums(X * beta_hat[G, , drop = FALSE])
  eta <- as.numeric(At) * XB

  res2 <- as.numeric((Yt - eta)^2)

  # 按 G 分组求平方残差均值
  rss_by_group <- rowsum(res2, group = G, reorder = FALSE)
  n_by_group <- tabulate(G, nbins = K)

  sig2 <- rep(NA_real_, K)

  groups_present <- as.integer(rownames(rss_by_group))
  valid_groups <- groups_present[
    groups_present >= 1L &
      groups_present <= K &
      n_by_group[groups_present] > 0L
  ]

  sig2[valid_groups] <- as.numeric(rss_by_group[as.character(valid_groups), 1]) /
    n_by_group[valid_groups]

  sig2
}


method_loc_lasso_ols_refit <- function(mydat,
                                       wloss = NULL,
                                       zero_tol = 0,
                                       df_tol = 1e-6,
                                       ...) {
  fit_loc <- method_loc_lasso(mydat, wloss = wloss, ...)

  beta_lasso <- as.matrix(fit_loc$beta)

  K <- mydat$K
  p <- mydat$p
  n <- mydat$n

  if (!(nrow(beta_lasso) == K && ncol(beta_lasso) == p)) {
    if (nrow(beta_lasso) == p && ncol(beta_lasso) == K) {
      beta_lasso <- t(beta_lasso)
    } else {
      stop(
        "Dimension mismatch for fit_loc$beta. Expected ",
        K, " x ", p, " or ", p, " x ", K,
        ", got ", nrow(beta_lasso), " x ", ncol(beta_lasso), ".",
        call. = FALSE
      )
    }
  }

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

  active_j <- colSums(abs(beta_lasso) > zero_tol) > 0
  dropped_j <- which(!active_j)
  p_active <- sum(active_j)

  beta_ols_full <- matrix(0, nrow = K, ncol = p)

  if (p_active == 0) {
    fit_loc$beta_lasso <- beta_lasso
    fit_loc$beta <- beta_ols_full
    fit_loc$beta_ols_active <- matrix(0, nrow = K, ncol = 0)
    fit_loc$active_j <- active_j
    fit_loc$dropped_j <- dropped_j
    fit_loc$p_active <- 0L
    fit_loc$p_full <- p
    fit_loc$refit <- "ols_after_loc_lasso"

    if (exists("bic_lm_loss", mode = "function")) {
      fit_loc$bic_ols <- bic_lm_loss(
        y = mydat$Ytilde,
        eta = rep(0, n),
        df = 0,
        n = n
      )
    }

    fit_loc$df_ols <- 0L

    return(fit_loc)
  }

  mydat_active <- mydat
  mydat_active$X <- mydat$X[, active_j, drop = FALSE]
  mydat_active$p <- p_active

  Z_active <- build_block_design(
    X = mydat_active$X,
    Atilde = mydat_active$Atilde,
    G = mydat_active$G,
    K = K
  )

  y <- mydat_active$Ytilde

  ols_fit <- stats::lm.wfit(
    x = Z_active,
    y = y,
    w = wloss
  )

  b_active <- as.numeric(ols_fit$coefficients)

  # lm.wfit returns NA for aliased columns if the active design is rank-deficient.
  # Setting them to zero keeps the output finite and preserves the selected format.
  b_active[!is.finite(b_active)] <- 0

  beta_ols_active <- mat_beta(b_active, K, p_active)

  beta_ols_full[, active_j] <- beta_ols_active

  eta_ols <- drop(build_block_design(
    X = mydat$X,
    Atilde = mydat$Atilde,
    G = mydat$G,
    K = K
  ) %*% vec_beta(beta_ols_full))

  df_ols <- sum(abs(beta_ols_full) > df_tol)

  fit_loc$beta_lasso <- beta_lasso
  fit_loc$beta <- beta_ols_full
  fit_loc$beta_ols_active <- beta_ols_active
  fit_loc$active_j <- active_j
  fit_loc$dropped_j <- dropped_j
  fit_loc$p_active <- p_active
  fit_loc$p_full <- p
  fit_loc$eta_ols <- eta_ols
  fit_loc$df_ols <- df_ols
  fit_loc$refit <- "ols_after_loc_lasso"

  if (exists("bic_lm_loss", mode = "function")) {
    fit_loc$bic_ols <- bic_lm_loss(
      y = mydat$Ytilde,
      eta = eta_ols,
      df = df_ols,
      n = n
    )
  }

  fit_loc
}
