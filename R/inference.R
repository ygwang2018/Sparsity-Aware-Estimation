safe_inv <- function(M, ridge = 1e-8) {
  M <- as.matrix(M)

  out <- tryCatch(
    solve(M),
    error = function(e) NULL
  )

  if (!is.null(out) && all(is.finite(out))) {
    return(out)
  }

  out <- tryCatch(
    solve(M + ridge * diag(nrow(M))),
    error = function(e) NULL
  )

  if (!is.null(out) && all(is.finite(out))) {
    return(out)
  }

  MASS::ginv(M)
}

normalize_beta_Kp <- function(beta, K, p, name = "beta") {
  beta <- as.matrix(beta)

  if (nrow(beta) == K && ncol(beta) == p) {
    return(beta)
  }

  if (nrow(beta) == p && ncol(beta) == K) {
    return(t(beta))
  }

  stop(
    "Dimension mismatch for ", name, ". Expected K x p = ",
    K, " x ", p, " or p x K = ", p, " x ", K,
    ", got ", nrow(beta), " x ", ncol(beta), ".",
    call. = FALSE
  )
}

normalize_contrast_a <- function(a, K, p) {
  if (is.null(a)) {
    stop("Contrast a cannot be NULL.", call. = FALSE)
  }

  if (is.matrix(a) || is.data.frame(a)) {
    a <- as.matrix(a)

    if (nrow(a) == K && ncol(a) == p) {
      return(a)
    }

    if (nrow(a) == p && ncol(a) == K) {
      return(t(a))
    }

    stop(
      "Dimension mismatch for contrast a. Expected K x p = ",
      K, " x ", p, " or p x K = ", p, " x ", K,
      ", got ", nrow(a), " x ", ncol(a), ".",
      call. = FALSE
    )
  }

  a <- as.numeric(a)

  if (length(a) != K * p) {
    stop(
      "Length mismatch for contrast a. Expected length K*p = ",
      K * p, ", got ", length(a), ".",
      call. = FALSE
    )
  }

  mat_beta(a, K, p)
}

make_beta1_contrast <- function(j, K, p) {
  if (length(j) != 1 || !is.finite(j) || j < 1 || j > p) {
    stop("j must be an integer in 1:p.", call. = FALSE)
  }

  a <- matrix(0, nrow = K, ncol = p)
  a[1, as.integer(j)] <- 1
  a
}

get_homo_nonzero_j_from_beta0 <- function(beta0,
                                          K,
                                          p,
                                          active_tol = 1e-12,
                                          homo_tol = 1e-12) {
  beta0 <- normalize_beta_Kp(beta0, K, p, name = "beta0")

  active <- abs(beta0) > active_tol
  nonzero_j <- apply(active, 2, any)
  homo_j <- get_homo_set(beta0, tol = homo_tol)

  which(nonzero_j & homo_j)
}

build_T_from_structure <- function(active, homo, K, p) {
  active <- as.matrix(active)
  homo <- as.logical(homo)

  if (!(nrow(active) == K && ncol(active) == p)) {
    stop("active must be K x p.", call. = FALSE)
  }

  if (length(homo) != p) {
    stop("homo must have length p.", call. = FALSE)
  }

  q <- K * p

  Tlist <- list()
  col_info <- data.frame(
    type = character(0),
    j = integer(0),
    k = integer(0),
    stringsAsFactors = FALSE
  )

  for (j in seq_len(p)) {
    if (all(!active[, j])) {
      next
    }

    if (isTRUE(homo[j])) {
      tcol <- rep(0, q)

      for (k in seq_len(K)) {
        idx <- (k - 1) * p + j
        tcol[idx] <- 1
      }

      Tlist[[length(Tlist) + 1]] <- tcol

      col_info <- rbind(
        col_info,
        data.frame(
          type = "homo",
          j = j,
          k = NA_integer_,
          stringsAsFactors = FALSE
        )
      )

    } else {
      for (k in seq_len(K)) {
        if (active[k, j]) {
          tcol <- rep(0, q)
          idx <- (k - 1) * p + j
          tcol[idx] <- 1

          Tlist[[length(Tlist) + 1]] <- tcol

          col_info <- rbind(
            col_info,
            data.frame(
              type = "het",
              j = j,
              k = k,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }

  if (length(Tlist) == 0) {
    T <- matrix(0, nrow = q, ncol = 0)
  } else {
    T <- do.call(cbind, Tlist)
  }

  list(
    T = T,
    col_info = col_info
  )
}

get_structure_from_beta <- function(beta_hat,
                                    active_tol = 1e-6,
                                    homo_tol = 1e-6) {
  beta_hat <- as.matrix(beta_hat)

  active <- abs(beta_hat) > active_tol
  homo <- get_homo_set(beta_hat, tol = homo_tol)

  list(
    active = active,
    homo = homo
  )
}

get_oracle_structure <- function(beta0,
                                 active_tol = 1e-12,
                                 homo_tol = 1e-12) {
  beta0 <- as.matrix(beta0)

  active <- abs(beta0) > active_tol
  homo <- get_homo_set(beta0, tol = homo_tol)

  list(
    active = active,
    homo = homo
  )
}

infer_contrast_structured <- function(mydat,
                                      a,
                                      structure,
                                      beta_input = NULL,
                                      beta0 = NULL,
                                      wloss = NULL,
                                      refit = TRUE,
                                      ridge = 1e-8,
                                      level = 0.95) {

  Yt <- as.numeric(mydat$Ytilde)
  At <- as.numeric(mydat$Atilde)
  X  <- as.matrix(mydat$X)
  G  <- as.integer(mydat$G)
  K  <- mydat$K
  p  <- mydat$p
  n  <- mydat$n

  if (is.null(wloss)) {
    wloss <- rep(1, n)
  }

  if (length(wloss) != n) {
    stop("Length of wloss must equal mydat$n.", call. = FALSE)
  }

  a_mat <- normalize_contrast_a(a, K, p)
  a_vec <- vec_beta(a_mat)

  Z <- build_block_design(
    X = X,
    Atilde = At,
    G = G,
    K = K
  )

  active <- structure$active
  homo <- structure$homo

  TT <- build_T_from_structure(
    active = active,
    homo = homo,
    K = K,
    p = p
  )

  Tmat <- TT$T
  d <- ncol(Tmat)

  if (d == 0) {
    beta_hat <- matrix(0, K, p)
    psi_hat <- sum(a_mat * beta_hat)

    psi0 <- NA_real_
    bias <- NA_real_
    cover <- NA

    if (!is.null(beta0)) {
      beta0 <- normalize_beta_Kp(beta0, K, p, name = "beta0")
      psi0 <- sum(a_mat * beta0)
      bias <- psi_hat - psi0
      cover <- as.numeric(psi0 >= psi_hat && psi0 <= psi_hat)
    }

    return(list(
      psi_hat = psi_hat,
      se = 0,
      var = 0,
      ci_low = psi_hat,
      ci_high = psi_hat,
      psi0 = psi0,
      bias = bias,
      cover = cover,
      beta_infer = beta_hat,
      theta_hat = numeric(0),
      d = 0,
      structure = structure,
      col_info = TT$col_info
    ))
  }

  Zred <- Z %*% Tmat

  if (isTRUE(refit)) {
    sw <- sqrt(wloss)

    fit <- stats::lm.fit(
      x = Zred * sw,
      y = Yt * sw
    )

    theta_hat <- as.numeric(fit$coefficients)
    theta_hat[is.na(theta_hat)] <- 0

  } else {
    if (is.null(beta_input)) {
      stop("beta_input is required when refit = FALSE.", call. = FALSE)
    }

    beta_input <- normalize_beta_Kp(beta_input, K, p, name = "beta_input")
    b_input <- vec_beta(beta_input)

    theta_hat <- safe_inv(crossprod(Tmat) + ridge * diag(d)) %*%
      crossprod(Tmat, b_input)

    theta_hat <- as.numeric(theta_hat)
  }

  b_hat_vec <- as.numeric(Tmat %*% theta_hat)
  beta_hat <- mat_beta(b_hat_vec, K, p)

  eta <- drop(Zred %*% theta_hat)
  eps_hat <- Yt - eta

  WZred <- Zred * as.numeric(wloss)

  Gamma_hat <- crossprod(Zred, WZred) / n

  Zomega <- Zred * as.numeric(wloss * eps_hat)

  Omega_hat <- crossprod(Zomega) / n

  Gamma_inv <- safe_inv(Gamma_hat, ridge = ridge)

  Sigma_theta_hat <- Gamma_inv %*% Omega_hat %*% Gamma_inv

  c_vec <- drop(t(Tmat) %*% a_vec)

  var_psi <- as.numeric(
    drop(t(c_vec) %*% Sigma_theta_hat %*% c_vec) / n
  )

  var_psi <- max(var_psi, 0)
  se_psi <- sqrt(var_psi)

  psi_hat <- sum(a_mat * beta_hat)

  alpha <- 1 - level
  zcrit <- stats::qnorm(1 - alpha / 2)

  ci_low <- psi_hat - zcrit * se_psi
  ci_high <- psi_hat + zcrit * se_psi

  psi0 <- NA_real_
  bias <- NA_real_
  cover <- NA

  if (!is.null(beta0)) {
    beta0 <- normalize_beta_Kp(beta0, K, p, name = "beta0")

    psi0 <- sum(a_mat * beta0)
    bias <- psi_hat - psi0
    cover <- as.numeric(ci_low <= psi0 && psi0 <= ci_high)
  }

  list(
    psi_hat = psi_hat,
    se = se_psi,
    var = var_psi,
    ci_low = ci_low,
    ci_high = ci_high,
    psi0 = psi0,
    bias = bias,
    cover = cover,
    beta_infer = beta_hat,
    theta_hat = theta_hat,
    Sigma_theta = Sigma_theta_hat,
    Gamma = Gamma_hat,
    Omega = Omega_hat,
    eps_hat = eps_hat,
    d = d,
    structure = structure,
    col_info = TT$col_info,
    refit = refit
  )
}

inference_oracle_prop <- function(mydat,
                                  fit_prop,
                                  beta0,
                                  a,
                                  wloss = NULL,
                                  fit_oracle = NULL,
                                  beta_oracle = NULL,
                                  prop_name = "Int-Prop",
                                  oracle_name = "Int-Oracle",
                                  refit_prop = TRUE,
                                  refit_oracle = TRUE,
                                  active_tol = 1e-6,
                                  homo_tol = 1e-6,
                                  level = 0.95,
                                  target = "custom_a",
                                  estimand = "a_beta",
                                  j = NA_integer_,
                                  k = NA_integer_) {

  K <- mydat$K
  p <- mydat$p

  beta0 <- normalize_beta_Kp(beta0, K, p, name = "beta0")

  str_oracle <- get_oracle_structure(
    beta0 = beta0,
    active_tol = 1e-12,
    homo_tol = 1e-12
  )

  beta_input_oracle <- NULL

  if (!isTRUE(refit_oracle)) {
    if (!is.null(beta_oracle)) {
      beta_input_oracle <- beta_oracle
    } else if (!is.null(fit_oracle) && !is.null(fit_oracle$beta)) {
      beta_input_oracle <- fit_oracle$beta
    } else {
      warning(
        "refit_oracle = FALSE but neither beta_oracle nor fit_oracle$beta was supplied. ",
        "Use beta0 as beta_input_oracle fallback. ",
        "This is usually not desired for evaluating estimator performance."
      )
      beta_input_oracle <- beta0
    }
  }

  inf_oracle <- infer_contrast_structured(
    mydat = mydat,
    a = a,
    structure = str_oracle,
    beta_input = beta_input_oracle,
    beta0 = beta0,
    wloss = wloss,
    refit = refit_oracle,
    level = level
  )

  if (is.null(fit_prop) || is.null(fit_prop$beta)) {
    stop("fit_prop$beta is required.", call. = FALSE)
  }

  str_prop <- get_structure_from_beta(
    beta_hat = fit_prop$beta,
    active_tol = active_tol,
    homo_tol = homo_tol
  )

  inf_prop <- infer_contrast_structured(
    mydat = mydat,
    a = a,
    structure = str_prop,
    beta_input = fit_prop$beta,
    beta0 = beta0,
    wloss = wloss,
    refit = refit_prop,
    level = level
  )

  out <- rbind(
    data.frame(
      method = oracle_name,
      target = target,
      estimand = estimand,
      j = j,
      k = k,
      psi_hat = inf_oracle$psi_hat,
      se = inf_oracle$se,
      var = inf_oracle$var,
      ci_low = inf_oracle$ci_low,
      ci_high = inf_oracle$ci_high,
      psi0 = inf_oracle$psi0,
      bias = inf_oracle$bias,
      cover = inf_oracle$cover,
      d = inf_oracle$d,
      refit = refit_oracle,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = prop_name,
      target = target,
      estimand = estimand,
      j = j,
      k = k,
      psi_hat = inf_prop$psi_hat,
      se = inf_prop$se,
      var = inf_prop$var,
      ci_low = inf_prop$ci_low,
      ci_high = inf_prop$ci_high,
      psi0 = inf_prop$psi0,
      bias = inf_prop$bias,
      cover = inf_prop$cover,
      d = inf_prop$d,
      refit = refit_prop,
      stringsAsFactors = FALSE
    )
  )

  attr(out, "inf_oracle") <- inf_oracle
  attr(out, "inf_prop") <- inf_prop

  out
}

inference_oracle_prop_with_homo <- function(mydat,
                                            fit_prop,
                                            beta0,
                                            a = NULL,
                                            wloss = NULL,
                                            fit_oracle = NULL,
                                            beta_oracle = NULL,
                                            prop_name = "Int-Prop",
                                            oracle_name = "Int-Oracle",
                                            refit_prop = TRUE,
                                            refit_oracle = TRUE,
                                            active_tol = 1e-6,
                                            homo_tol = 1e-6,
                                            homo_active_tol_beta0 = 1e-12,
                                            homo_tol_beta0 = 1e-12,
                                            level = 0.95,
                                            include_custom_a = TRUE,
                                            include_homo_beta1 = TRUE) {
  K <- mydat$K
  p <- mydat$p

  beta0 <- normalize_beta_Kp(beta0, K, p, name = "beta0")

  out <- NULL

  if (isTRUE(include_custom_a) && !is.null(a)) {
    out_a <- inference_oracle_prop(
      mydat = mydat,
      fit_prop = fit_prop,
      beta0 = beta0,
      a = a,
      wloss = wloss,
      fit_oracle = fit_oracle,
      beta_oracle = beta_oracle,
      prop_name = prop_name,
      oracle_name = oracle_name,
      refit_prop = refit_prop,
      refit_oracle = refit_oracle,
      active_tol = active_tol,
      homo_tol = homo_tol,
      level = level,
      target = "custom_a",
      estimand = "a_beta",
      j = NA_integer_,
      k = NA_integer_
    )

    out <- rbind(out, out_a)
  }

  if (isTRUE(include_homo_beta1)) {
    homo_j <- get_homo_nonzero_j_from_beta0(
      beta0 = beta0,
      K = K,
      p = p,
      active_tol = homo_active_tol_beta0,
      homo_tol = homo_tol_beta0
    )

    if (length(homo_j) > 0) {
      for (jj in homo_j) {
        aj <- make_beta1_contrast(j = jj, K = K, p = p)

        out_j <- inference_oracle_prop(
          mydat = mydat,
          fit_prop = fit_prop,
          beta0 = beta0,
          a = aj,
          wloss = wloss,
          fit_oracle = fit_oracle,
          beta_oracle = beta_oracle,
          prop_name = prop_name,
          oracle_name = oracle_name,
          refit_prop = refit_prop,
          refit_oracle = refit_oracle,
          active_tol = active_tol,
          homo_tol = homo_tol,
          level = level,
          target = paste0("beta_", jj, "^(1)"),
          estimand = "homo_beta1",
          j = jj,
          k = 1L
        )

        out <- rbind(out, out_j)
      }
    }
  }

  out
}