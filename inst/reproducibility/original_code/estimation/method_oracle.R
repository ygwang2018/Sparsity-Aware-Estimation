method_oracle <- function(mydat, beta0, wloss = NULL) {
  Yt <- mydat$Ytilde
  At <- mydat$Atilde
  X <- mydat$X
  G <- mydat$G
  K <- mydat$K
  p <- mydat$p
  n <- mydat$n

  if (is.null(wloss)) wloss <- rep(1, n)

  true_active <- abs(beta0) > 1e-12
  true_homo <- get_homo_set(beta0, tol = 1e-12)

  cols <- list()
  names_cols <- c()

  # Heterogeneous active coefficients: each subgroup has its own coefficient
  # Homogeneous nonzero block: one common coefficient
  Zlist <- list()

  for (j in 1:p) {
    if (all(!true_active[, j])) next

    if (true_homo[j]) {
      z <- At * X[, j]
      Zlist[[length(Zlist) + 1]] <- z
      names_cols <- c(names_cols, paste0("homo_", j))
    } else {
      for (k in 1:K) {
        if (true_active[k, j]) {
          z <- At * X[, j] * as.numeric(G == k)
          Zlist[[length(Zlist) + 1]] <- z
          names_cols <- c(names_cols, paste0("het_", j, "_", k))
        }
      }
    }
  }

  if (length(Zlist) == 0) {
    return(list(beta = matrix(0, K, p)))
  }

  Zor <- do.call(cbind, Zlist)
  sw <- sqrt(wloss)

  fit <- lm.fit(Zor * sw, Yt * sw)
  theta <- coef(fit)
  theta[is.na(theta)] <- 0

  beta_hat <- matrix(0, K, p)

  id <- 1
  for (j in 1:p) {
    if (all(!true_active[, j])) next

    if (true_homo[j]) {
      beta_hat[, j] <- theta[id]
      id <- id + 1
    } else {
      for (k in 1:K) {
        if (true_active[k, j]) {
          beta_hat[k, j] <- theta[id]
          id <- id + 1
        }
      }
    }
  }

  list(beta = beta_hat)
}