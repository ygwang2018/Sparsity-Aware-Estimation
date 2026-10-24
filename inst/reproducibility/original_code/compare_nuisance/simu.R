simu <- function(n, K, p, rho, beta0, treat_prob, balance, nui_app,
                 seedstart, Nsimu, taskID) {

  set.seed(seedstart)
  seednum <- sample(1:10^6, Nsimu)

  ALL_EVAL <- NULL
  ALL_BHAT_WIDE <- NULL

  for (b in 1:Nsimu) {
    cat("Simulation:", b, "/", Nsimu, "\n")

    one <- once_simu(
      n = n,
      K = K,
      p = p,
      rho = rho,
      beta0 = beta0,
      treat_prob = treat_prob,
      balance = balance,
      seed = seednum[b],nui_app
    )

    # -------------------------------------------------------
    # 1. 保存 evaluation results
    # -------------------------------------------------------
    ev <- one$eval
    ev$simu <- b
    ev$seed <- seednum[b]

    ALL_EVAL <- rbind(ALL_EVAL, ev)

    # -------------------------------------------------------
    # 2. 保存 beta_hat wide format
    # -------------------------------------------------------
    bhat_wide <- beta_hat_to_wide(
      beta_hat_list = one$beta_hat,
      simu = b,
      seed = seednum[b]
    )

    ALL_BHAT_WIDE <- rbind(ALL_BHAT_WIDE, bhat_wide)

    # -------------------------------------------------------
    # 3. 每次 simulation 后都写出，防止中途断掉
    # -------------------------------------------------------

    # raw evaluation
    write.csv(
      ALL_EVAL,
      paste0(taskID, "_AllEval.csv"),
      row.names = FALSE
    )

    # beta_hat wide
    write.csv(
      ALL_BHAT_WIDE,
      paste0(taskID, "_AllBhat_wide.csv"),
      row.names = FALSE
    )

    # summary
    method_order <- c(
  "Int-Oracle-KN",
  "Int-Oracle-KN-w",
  "Int-Oracle",
  "Int-Oracle-w",
  "Int-Prop-KN",
  "Int-Prop-KN-w",
  "Int-Prop",
  "Int-Prop-w",
  "Int-R2Lasso",
  "Int-R2Lasso-w",
  "Int-GLasso",
  "Int-GLasso-w",
  "Int-HLasso",
  "Int-HLasso-w",
  "Loc-Lasso"
)

ALL_EVAL$method <- factor(
  ALL_EVAL$method,
  levels = method_order
)

SUM <- aggregate(
  cbind(MSE, Var_TPR, Var_FDR, Homo_TPR, Homo_FDR) ~ method,
  data = ALL_EVAL,
  FUN = mean,
  na.rm = TRUE
)

SUM <- SUM[match(method_order, as.character(SUM$method)), ]
SUM <- SUM[!is.na(SUM$method), ]
SUM$method <- as.character(SUM$method)

    write.csv(
      SUM,
      paste0(taskID, "_Summary.csv"),
      row.names = FALSE
    )

    # RDS 总结果，也保存一份，方便之后 R 里直接读取
    saveRDS(
      list(
        ALL_EVAL = ALL_EVAL,
        ALL_BHAT_WIDE = ALL_BHAT_WIDE,
        seednum = seednum,
        settings = list(
          n = n,
          K = K,
          p = p,
          rho = rho,
          beta0 = beta0,
          treat_prob = treat_prob,
          balance = balance,
          seedstart = seedstart,
          Nsimu = Nsimu,
          taskID = taskID
        )
      ),
      paste0(taskID, "_AllResults.rds")
    )

    print(SUM)
  }

  invisible(list(
    ALL_EVAL = ALL_EVAL,
    ALL_BHAT_WIDE = ALL_BHAT_WIDE,
    seednum = seednum
  ))
} 




# ------------------------------------------------------------
# Prepare data when nuisance functions are known.
#
# Uses:
#   m_true(X,G) = E[Y | X,G]
#   e0(X,G)    = P(A = 1 | X,G)
#
# Output has the same structure as prepare_dat().
# ------------------------------------------------------------
prepare_dat_known_nuisance <- function(dat, K, trim_e = 1e-4) {

  if (!all(c("Y", "A", "G") %in% colnames(dat))) {
    stop("dat must contain columns Y, A, and G.", call. = FALSE)
  }

  Y <- as.numeric(dat$Y)
  A <- as.numeric(dat$A)
  G <- as.integer(dat$G)

  X <- get_X_from_dat(dat)
  X <- clean_matrix(X)

  n <- nrow(X)
  p <- ncol(X)

  m_true <- attr(dat, "m_true")
  e0 <- attr(dat, "e0")

  if (is.null(m_true)) {
    stop("dat does not contain attr(dat, 'm_true'). Please update gdat.R.", call. = FALSE)
  }

  if (is.null(e0)) {
    stop("dat does not contain attr(dat, 'e0'). Please update gdat.R.", call. = FALSE)
  }

  m_true <- as.numeric(m_true)
  e0 <- as.numeric(e0)

  if (length(m_true) != n) {
    stop("Length of m_true does not match sample size.", call. = FALSE)
  }

  if (length(e0) != n) {
    stop("Length of e0 does not match sample size.", call. = FALSE)
  }

  e0 <- pmin(pmax(e0, trim_e), 1 - trim_e)

  Ytilde <- Y - m_true
  Atilde <- A - e0

  if (any(!is.finite(Ytilde))) {
    stop("Non-finite Ytilde in prepare_dat_known_nuisance().", call. = FALSE)
  }

  if (any(!is.finite(Atilde))) {
    stop("Non-finite Atilde in prepare_dat_known_nuisance().", call. = FALSE)
  }

  list(
    Y = Y,
    A = A,
    G = G,
    X = X,
    K = as.integer(K),
    p = p,
    n = n,

    mhat = m_true,
    ehat = e0,
    Ytilde = Ytilde,
    Atilde = Atilde,

    foldid = rep(NA_integer_, n),
    nfold = NA_integer_,

    outcome_method = "known",
    propensity_method = "known",
    outcome_args = NULL,
    propensity_args = NULL,

    use_group_for_nuisance = TRUE,
    add_group_interactions = FALSE,
    calibrate_propensity = FALSE,
    trim_e = trim_e
  )
}







once_simu <- function(n, K, p, rho, beta0, treat_prob, balance, seed,
                      nui_app = "bagging") {

  set.seed(seed)

  dat <- gdat(
    n = n,
    K = K,
    p = p,
    rho = rho,
    beta0 = beta0,
    treat_prob = treat_prob,
    balance = balance
  )

  # ---------------------------------------------------------
  # A. Estimated nuisance version
  # ---------------------------------------------------------
   if (nui_app == "rf") {
    mydat <- prepare_dat(
      dat,
      K = K,
      nfold = 10,
      outcome_method = "rf",
      propensity_method = "rf",
      outcome_args = list(
        ntree = 150
      ),
      propensity_args = list(
        ntree = 150
      ),
      seed = seed + 1000,
      trim_e = 1e-4,
      use_group_for_nuisance = TRUE,
      add_group_interactions = FALSE,
      calibrate_propensity = TRUE
    )

  } else {     
    mydat <- prepare_dat(
      dat,
      K = K,
      nfold = 10,
      outcome_method = "bagging",
      propensity_method = "bagging",
      use_dc_screening = FALSE, 
      outcome_args = list(nbagg = 50),
      propensity_args = list(nbagg = 50),
      seed = seed + 1000,
      trim_e = 1e-4,
      use_group_for_nuisance = TRUE,
      add_group_interactions = FALSE,
      calibrate_propensity = TRUE
    )
    }
  }

  # ---------------------------------------------------------
  # B. Known nuisance version
  # ---------------------------------------------------------
  mydat_kn <- prepare_dat_known_nuisance(
    dat = dat,
    K = K,
    trim_e = 1e-4
  )

  # ---------------------------------------------------------
  # 1. Loc-Lasso
  # ---------------------------------------------------------
  fit_loc <- method_loc_lasso(mydat)
  
  fit_loc_ols <- method_loc_lasso_ols_refit(mydat)
  betahat_ini <- fit_loc_ols$beta
  # ---------------------------------------------------------
  # 2. Cross-fitted loss weights for estimated nuisance version
  # ---------------------------------------------------------
  fit_wloss <- make_wloss_cf_local(
    mydat = mydat,
    nfold = 10,
    seed = seed + 1000,
    gamma = 1,
    nlambda = 50,
    lambda.min.ratio = 0.02,
    eps_weight = 1e-4,
    sigma_floor = 1e-6,
    sigma_cap = Inf,
    w_floor = 1e-3,
    w_cap = 1e3,
    denom = "n",
    normalize = TRUE,
    verbose = FALSE
  )

  wloss <- fit_wloss$wloss

  # ---------------------------------------------------------
  # 3. Cross-fitted loss weights for known nuisance version
  # ---------------------------------------------------------
  fit_wloss_kn <- make_wloss_cf_local(
    mydat = mydat_kn,
    nfold = 10,
    seed = seed + 2000,
    gamma = 1,
    nlambda = 50,
    lambda.min.ratio = 0.02,
    eps_weight = 1e-4,
    sigma_floor = 1e-6,
    sigma_cap = Inf,
    w_floor = 1e-3,
    w_cap = 1e3,
    denom = "n",
    normalize = TRUE,
    verbose = FALSE
  )

  wloss_kn <- fit_wloss_kn$wloss

  # ---------------------------------------------------------
  # 4. Int-GLasso
  # ---------------------------------------------------------
  fit_glasso <- method_glasso(
    mydat,
    beta_init = betahat_ini,
    wloss = NULL,
    lambda.min.ratio = 1e-03
  )

  fit_glasso_w <- method_glasso(
    mydat,
    beta_init = betahat_ini,
    wloss = wloss,
    lambda.min.ratio = 1e-03
  )

  # ---------------------------------------------------------
  # 5. Int-HLasso
  # ---------------------------------------------------------
  fit_hlasso <- method_hlasso(
    mydat,
    beta_init = betahat_ini,
    wloss = NULL, nlambda = 100,
    lambda.min.ratio = 1e-03
  )

  fit_hlasso_w <- method_hlasso(
    mydat,
    beta_init = betahat_ini,
    wloss = wloss, nlambda = 60,
    lambda.min.ratio = 1e-02
  )

  # ---------------------------------------------------------
  # 6. Int-R2Lasso
  # ---------------------------------------------------------
  fit_R2lasso <- method_R2lasso(
    mydat,
    beta_init = betahat_ini,
    wloss = NULL,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-4
  )

  fit_R2lasso_w <- method_R2lasso(
    mydat,
    beta_init = betahat_ini,
    wloss = wloss,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-4
  )

  # ---------------------------------------------------------
  # 7. Int-Prop with estimated nuisance
  # ---------------------------------------------------------
  fit_prop <- method_prop_admm(
    mydat,
    beta_init = betahat_ini,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-4,
    wloss = NULL
  )

  fit_prop_w <- method_prop_admm(
    mydat,
    beta_init = betahat_ini,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-4,
    wloss = wloss
  )

  # ---------------------------------------------------------
  # 8. Int-Prop with known nuisance
  # ---------------------------------------------------------
  fit_prop_kn <- method_prop_admm(
    mydat_kn,
    beta_init = betahat_ini,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-7,
    wloss = NULL
  )

  fit_prop_kn_w <- method_prop_admm(
    mydat_kn,
    beta_init = betahat_ini,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-7,
    wloss = wloss_kn
  )

  # ---------------------------------------------------------
  # 9. Oracle with estimated nuisance
  # ---------------------------------------------------------
  fit_oracle <- method_oracle(
    mydat,
    beta0 = t(beta0),
    wloss = NULL
  )

  fit_oracle_w <- method_oracle(
    mydat,
    beta0 = t(beta0),
    wloss = wloss
  )

  # ---------------------------------------------------------
  # 10. Oracle with known nuisance
  # ---------------------------------------------------------
  fit_oracle_kn <- method_oracle(
    mydat_kn,
    beta0 = t(beta0),
    wloss = NULL
  )

  fit_oracle_kn_w <- method_oracle(
    mydat_kn,
    beta0 = t(beta0),
    wloss = wloss_kn
  )

  # ---------------------------------------------------------
  # 11. Collect beta estimates
  # ---------------------------------------------------------
  beta_hat_list <- list(
  "Int-Oracle-KN" = fit_oracle_kn$beta * (abs(fit_oracle_kn$beta) > 0.01),
  "Int-Oracle-KN-w" = fit_oracle_kn_w$beta * (abs(fit_oracle_kn_w$beta) > 0.01),

  "Int-Oracle" = fit_oracle$beta * (abs(fit_oracle$beta) > 0.01),
  "Int-Oracle-w" = fit_oracle_w$beta * (abs(fit_oracle_w$beta) > 0.01),

  "Int-Prop-KN" = fit_prop_kn$beta * (abs(fit_prop_kn$beta) > 0.01),
  "Int-Prop-KN-w" = fit_prop_kn_w$beta * (abs(fit_prop_kn_w$beta) > 0.01),

  "Int-Prop" = fit_prop$beta * (abs(fit_prop$beta) > 0.01),
  "Int-Prop-w" = fit_prop_w$beta * (abs(fit_prop_w$beta) > 0.01),

  "Int-R2Lasso" = fit_R2lasso$beta * (abs(fit_R2lasso$beta) > 0.01),
  "Int-R2Lasso-w" = fit_R2lasso_w$beta * (abs(fit_R2lasso_w$beta) > 0.01),

  "Int-GLasso" = fit_glasso$beta * (abs(fit_glasso$beta) > 0.01),
  "Int-GLasso-w" = fit_glasso_w$beta * (abs(fit_glasso_w$beta) > 0.01),

  "Int-HLasso" = fit_hlasso$beta * (abs(fit_hlasso$beta) > 0.01),
  "Int-HLasso-w" = fit_hlasso_w$beta * (abs(fit_hlasso_w$beta) > 0.01),

  "Loc-Lasso" = fit_loc$beta * (abs(fit_loc$beta) > 0.01)
)

  ev <- eval_methods(beta_hat_list, t(beta0))

  list(
    beta_hat = beta_hat_list,
    eval = ev,

    mydat = mydat,
    mydat_kn = mydat_kn,

    wloss = wloss,
    wloss_kn = wloss_kn,

    wloss_fit = fit_wloss,
    wloss_fit_kn = fit_wloss_kn
  )
}