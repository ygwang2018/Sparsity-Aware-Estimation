gdat <- function(n, K, p, rho, beta0, treat_prob = 0.3,
                 balance = TRUE, hetero_error = TRUE) {

  beta0 <- t(beta0)
  Sigma <- make_AR_cov(p, rho)

  if (balance) {
    nk <- rep(floor(n / K), K)
    nk[K] <- n - sum(nk[-K])
  } else {
    probG <- rep(c(2, 2, 3, 3),length.out=K)
    probG <- probG / sum(probG)
    nk <- as.vector(rmultinom(1, n, probG))
  }

  G <- rep(seq_len(K), times = nk)
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

  # subgroup-specific propensity score
f_e <- 0.20 * X[, 1] - 0.20 * (X[, 3]>0) +
       0.20 * (G - mean(seq_len(K)))

  lin_e <- uniroot(
    function(a) mean(plogis(a + f_e)) - treat_prob,
    interval = c(-20, 20)
  )$root + f_e

  e0 <- 1 / (1 + exp(-lin_e))
  e0 <- pmin(pmax(e0, 0.05), 0.95)

  A <- rbinom(n, 1, e0)

  # baseline mean m0(X,G) = E[Y | X,G,A=0]
  m0 <- -1 + 0.20 * (X[, 1] > 0) + 0.20 * X[, 2]^2 + 0.20 * G

  tau <- rep(0, n)
  for (k in seq_len(K)) {
    idx <- which(G == k)
    tau[idx] <- X[idx, , drop = FALSE] %*% beta0[k, ]
  }

  # true nuisance for residualization: E[Y | X, G]
  m_true <- m0 + e0 * tau

  if (hetero_error) {
    sigma_k <- sort(rep(c(0.8, 1.2), length.out = K))
  } else {
    sigma_k <- rep(1, K)
  }

  eps <- rnorm(n, sd = sigma_k[G])

  Y <- m0 + A * tau + eps

  out <- data.frame(
    Y = Y,
    A = A,
    G = G,
    X
  )

  attr(out, "m0") <- as.numeric(m0)
  attr(out, "m_true") <- as.numeric(m_true)
  attr(out, "e0") <- as.numeric(e0)
  attr(out, "tau") <- as.numeric(tau)
  attr(out, "sigma_k") <- as.numeric(sigma_k)

  out
}