simu <- function(n, K, p, rho, beta0, treat_prob, balance, nui_app,
                 seedstart, Nsimu, taskID, a = NULL,
                 save_raw_each = TRUE,
                 raw_dir = NULL) {

  set.seed(seedstart)
  seednum <- sample(1:10^6, Nsimu)

  ALL_EVAL <- NULL
  ALL_BHAT_WIDE <- NULL
  ALL_INFER <- NULL
  ALL_RAW_INDEX <- NULL

  if (is.null(raw_dir)) {
    raw_dir <- paste0(taskID, "_RawEachSim")
  }

  if (isTRUE(save_raw_each) && !dir.exists(raw_dir)) {
    dir.create(raw_dir, recursive = TRUE)
  }

  for (b in seq_len(Nsimu)) {
    cat("Simulation:", b, "/", Nsimu, "\n")

    one <- once_simu(
      n = n,
      K = K,
      p = p,
      rho = rho,
      beta0 = beta0,
      treat_prob = treat_prob,
      balance = balance,
      seed = seednum[b],
      nui_app = nui_app,
      a = a
    )

    raw_file <- NA_character_

    if (isTRUE(save_raw_each)) {
      raw_file <- file.path(
        raw_dir,
        paste0(taskID, "_sim", b, "_seed", seednum[b], "_raw.rds")
      )

      saveRDS(
        one,
        raw_file
      )
    }

    raw_index_one <- data.frame(
      simu = b,
      seed = seednum[b],
      raw_file = raw_file,
      stringsAsFactors = FALSE
    )

    ALL_RAW_INDEX <- rbind(ALL_RAW_INDEX, raw_index_one)

    write.csv(
      ALL_RAW_INDEX,
      paste0(taskID, "_RawIndex.csv"),
      row.names = FALSE
    )

    ev <- one$eval
    ev$simu <- b
    ev$seed <- seednum[b]

    ALL_EVAL <- rbind(ALL_EVAL, ev)

    write.csv(
      ALL_EVAL,
      paste0(taskID, "_AllEval.csv"),
      row.names = FALSE
    )

    bhat_wide <- beta_hat_to_wide(
      beta_hat_list = one$beta_hat,
      simu = b,
      seed = seednum[b]
    )

    ALL_BHAT_WIDE <- rbind(ALL_BHAT_WIDE, bhat_wide)

    write.csv(
      ALL_BHAT_WIDE,
      paste0(taskID, "_AllBhat_wide.csv"),
      row.names = FALSE
    )

    if (!is.null(one$infer) && nrow(one$infer) > 0) {
      inf <- one$infer
      inf$simu <- b
      inf$seed <- seednum[b]

      ALL_INFER <- rbind(ALL_INFER, inf)

      write.csv(
        ALL_INFER,
        paste0(taskID, "_AllInfer.csv"),
        row.names = FALSE
      )
    }

    SUM <- aggregate(
      cbind(MSE, Var_SP, Var_SN, Var_FDR, Homo_SP, Homo_SN) ~ method,
      data = ALL_EVAL,
      FUN = mean,
      na.rm = TRUE
    )

    write.csv(
      SUM,
      paste0(taskID, "_Summary.csv"),
      row.names = FALSE
    )

    INFER_SUM <- NULL

    if (!is.null(ALL_INFER) && nrow(ALL_INFER) > 0) {
      infer_mean <- aggregate(
        cbind(psi_hat, se, var, bias, cover, d) ~ method + target + estimand,
        data = ALL_INFER,
        FUN = mean,
        na.rm = TRUE
      )

      colnames(infer_mean)[colnames(infer_mean) == "psi_hat"] <- "mean_psi_hat"
      colnames(infer_mean)[colnames(infer_mean) == "se"] <- "mean_se"
      colnames(infer_mean)[colnames(infer_mean) == "var"] <- "mean_var"
      colnames(infer_mean)[colnames(infer_mean) == "bias"] <- "mean_bias"
      colnames(infer_mean)[colnames(infer_mean) == "cover"] <- "coverage"
      colnames(infer_mean)[colnames(infer_mean) == "d"] <- "mean_d"

      infer_sd <- aggregate(
        psi_hat ~ method + target + estimand,
        data = ALL_INFER,
        FUN = sd,
        na.rm = TRUE
      )

      colnames(infer_sd)[colnames(infer_sd) == "psi_hat"] <- "sd_psi_hat"

      infer_n <- aggregate(
        psi_hat ~ method + target + estimand,
        data = ALL_INFER,
        FUN = function(x) sum(is.finite(x))
      )

      colnames(infer_n)[colnames(infer_n) == "psi_hat"] <- "nsimu_used"

      psi0_tab <- aggregate(
        psi0 ~ method + target + estimand,
        data = ALL_INFER,
        FUN = function(x) {
          ux <- unique(x[is.finite(x)])
          if (length(ux) == 0) return(NA_real_)
          ux[1]
        }
      )

      j_tab <- aggregate(
        j ~ method + target + estimand,
        data = ALL_INFER,
        FUN = function(x) {
          ux <- unique(x[is.finite(x)])
          if (length(ux) == 0) return(NA_integer_)
          ux[1]
        }
      )

      k_tab <- aggregate(
        k ~ method + target + estimand,
        data = ALL_INFER,
        FUN = function(x) {
          ux <- unique(x[is.finite(x)])
          if (length(ux) == 0) return(NA_integer_)
          ux[1]
        }
      )

      INFER_SUM <- merge(
        infer_mean,
        infer_sd,
        by = c("method", "target", "estimand"),
        all.x = TRUE
      )

      INFER_SUM <- merge(
        INFER_SUM,
        infer_n,
        by = c("method", "target", "estimand"),
        all.x = TRUE
      )

      INFER_SUM <- merge(
        INFER_SUM,
        psi0_tab,
        by = c("method", "target", "estimand"),
        all.x = TRUE
      )

      INFER_SUM <- merge(
        INFER_SUM,
        j_tab,
        by = c("method", "target", "estimand"),
        all.x = TRUE
      )

      INFER_SUM <- merge(
        INFER_SUM,
        k_tab,
        by = c("method", "target", "estimand"),
        all.x = TRUE
      )

      INFER_SUM$se_over_sd <- INFER_SUM$mean_se / INFER_SUM$sd_psi_hat

      write.csv(
        INFER_SUM,
        paste0(taskID, "_InferSummary.csv"),
        row.names = FALSE
      )
    }

    saveRDS(
      list(
        ALL_EVAL = ALL_EVAL,
        ALL_BHAT_WIDE = ALL_BHAT_WIDE,
        ALL_INFER = ALL_INFER,
        ALL_RAW_INDEX = ALL_RAW_INDEX,
        seednum = seednum,
        settings = list(
          n = n,
          K = K,
          p = p,
          rho = rho,
          beta0 = beta0,
          treat_prob = treat_prob,
          balance = balance,
          nui_app = nui_app,
          seedstart = seedstart,
          Nsimu = Nsimu,
          taskID = taskID,
          a = a,
          save_raw_each = save_raw_each,
          raw_dir = raw_dir
        )
      ),
      paste0(taskID, "_AllResults.rds")
    )

    print(SUM)

    if (!is.null(INFER_SUM)) {
      print(INFER_SUM)
    }
  }

  invisible(list(
    ALL_EVAL = ALL_EVAL,
    ALL_BHAT_WIDE = ALL_BHAT_WIDE,
    ALL_INFER = ALL_INFER,
    ALL_RAW_INDEX = ALL_RAW_INDEX,
    seednum = seednum
  ))
}





# ------------------------------------------------------------
# Prepare data when nuisance functions are known.
# Uses true m(X,G)=E[Y|X,G] and e(X,G)=P(A=1|X,G).
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

    foldid = rep(1L, n),
    nfold = 1L,

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
                      nui_app = "bagging", a = NULL) {

  set.seed(seed)

  beta0_Kp <- t(beta0)

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
  # Estimated nuisance version
  # ---------------------------------------------------------
  if (nui_app == "rf") {
    mydat <- prepare_dat(
      dat,
      K = K,
      nfold = 10,
      outcome_method = "rf",
      propensity_method = "rf",
      outcome_args = list(
        ntree = 500, nodesize = 30, mtry = 8
      ),
      propensity_args = list(
        ntree = 500, nodesize = 80, mtry = 8
      ),
      seed = seed,
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
         outcome_args = list(nbagg = 200, mtry = 8),
         propensity_args = list(nbagg = 200, mtry = 8),
         seed = seed,
         trim_e = 1e-4,
         use_group_for_nuisance = TRUE,
         add_group_interactions = !TRUE,
         calibrate_propensity = TRUE )
  }

  # ---------------------------------------------------------
  # Known nuisance version
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
  # 2. Loss weights for estimated nuisance
  # ---------------------------------------------------------
  fit_wloss <- make_wloss_cf_local(
    mydat = mydat,
    nfold = 10,
    seed = seed + 10000,
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
  # 3. Loss weights for known nuisance
  # ---------------------------------------------------------
  fit_wloss_kn <- make_wloss_cf_local(
    mydat = mydat_kn,
    nfold = 10,
    seed = seed + 20000,
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
  # 6. Int-Prop with estimated nuisance
  # ---------------------------------------------------------
  fit_prop <- method_prop_admm(
    mydat,
    beta_init = betahat_ini ,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-4,
    wloss = NULL
  )

  fit_prop_w <- method_prop_admm(
    mydat,
    beta_init = betahat_ini ,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-4,
    wloss = wloss
  )

  # ---------------------------------------------------------
  # 7. Int-Prop with known nuisance
  # ---------------------------------------------------------
  fit_prop_kn <- method_prop_admm(
    mydat_kn,
    beta_init = betahat_ini,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-4,
    wloss = NULL
  )

  fit_prop_kn_w <- method_prop_admm(
    mydat_kn,
    beta_init = betahat_ini,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-4,
    wloss = wloss_kn
  )

  # ---------------------------------------------------------
  # 8. Oracle with estimated nuisance
  # ---------------------------------------------------------
  fit_oracle <- method_oracle(
    mydat,
    beta0 = beta0_Kp,
    wloss = NULL
  )

  fit_oracle_w <- method_oracle(
    mydat,
    beta0 = beta0_Kp,
    wloss = wloss
  )

  # ---------------------------------------------------------
  # 9. Oracle with known nuisance
  # ---------------------------------------------------------
  fit_oracle_kn <- method_oracle(
    mydat_kn,
    beta0 = beta0_Kp,
    wloss = NULL
  )

  fit_oracle_kn_w <- method_oracle(
    mydat_kn,
    beta0 = beta0_Kp,
    wloss = wloss_kn
  )

  # ---------------------------------------------------------
  # 10. Inference
  # ---------------------------------------------------------
  infer_unw <- inference_oracle_prop_with_homo(
    mydat = mydat,
    fit_prop = fit_prop,
    fit_oracle = fit_oracle,
    beta0 = beta0_Kp,
    a = a,
    wloss = NULL,
    prop_name = "Int-Prop",
    oracle_name = "Int-Oracle",
    refit_prop = FALSE,
    refit_oracle = FALSE,
    active_tol = 1e-6,
    homo_tol = 1e-4,
    level = 0.95,
    include_custom_a = !is.null(a),
    include_homo_beta1 = TRUE
  )

  infer_w <- inference_oracle_prop_with_homo(
    mydat = mydat,
    fit_prop = fit_prop_w,
    fit_oracle = fit_oracle_w,
    beta0 = beta0_Kp,
    a = a,
    wloss = wloss,
    prop_name = "Int-Prop-w",
    oracle_name = "Int-Oracle-w",
    refit_prop = FALSE,
    refit_oracle = FALSE,
    active_tol = 1e-6,
    homo_tol = 1e-4,
    level = 0.95,
    include_custom_a = !is.null(a),
    include_homo_beta1 = TRUE
  )

  infer_unw_kn <- inference_oracle_prop_with_homo(
    mydat = mydat_kn,
    fit_prop = fit_prop_kn,
    fit_oracle = fit_oracle_kn,
    beta0 = beta0_Kp,
    a = a,
    wloss = NULL,
    prop_name = "Int-Prop-KN",
    oracle_name = "Int-Oracle-KN",
    refit_prop = FALSE,
    refit_oracle = FALSE,
    active_tol = 1e-6,
    homo_tol = 1e-4,
    level = 0.95,
    include_custom_a = !is.null(a),
    include_homo_beta1 = TRUE
  )

  infer_w_kn <- inference_oracle_prop_with_homo(
    mydat = mydat_kn,
    fit_prop = fit_prop_kn_w,
    fit_oracle = fit_oracle_kn_w,
    beta0 = beta0_Kp,
    a = a,
    wloss = wloss_kn,
    prop_name = "Int-Prop-KN-w",
    oracle_name = "Int-Oracle-KN-w",
    refit_prop = FALSE,
    refit_oracle = FALSE,
    active_tol = 1e-6,
    homo_tol = 1e-4,
    level = 0.95,
    include_custom_a = !is.null(a),
    include_homo_beta1 = TRUE
  )

  infer <- rbind(infer_unw, infer_w, infer_unw_kn, infer_w_kn)

  # ---------------------------------------------------------
  # 11. Collect beta estimates
  # ---------------------------------------------------------
  beta_hat_list <- list(
    "Loc-Lasso" = fit_loc$beta,
    "Int-Prop" = fit_prop$beta,
    "Int-Prop-w" = fit_prop_w$beta,
    "Int-Prop-KN" = fit_prop_kn$beta,
    "Int-Prop-KN-w" = fit_prop_kn_w$beta,
    "Int-Oracle" = fit_oracle$beta,
    "Int-Oracle-w" = fit_oracle_w$beta,
    "Int-Oracle-KN" = fit_oracle_kn$beta,
    "Int-Oracle-KN-w" = fit_oracle_kn_w$beta
  )

  ev <- eval_methods(beta_hat_list, beta0_Kp)

  raw_fits <- list(
    fit_loc = fit_loc,
    fit_wloss = fit_wloss,
    fit_wloss_kn = fit_wloss_kn,
    fit_prop = fit_prop,
    fit_prop_w = fit_prop_w,
    fit_prop_kn = fit_prop_kn,
    fit_prop_kn_w = fit_prop_kn_w,
    fit_oracle = fit_oracle,
    fit_oracle_w = fit_oracle_w,
    fit_oracle_kn = fit_oracle_kn,
    fit_oracle_kn_w = fit_oracle_kn_w
  )

  list(
    beta_hat = beta_hat_list,
    eval = ev,
    infer = infer,
    dat = dat,
    mydat = mydat,
    mydat_kn = mydat_kn,
    wloss = wloss,
    wloss_kn = wloss_kn,
    raw_fits = raw_fits,
    beta0_Kp = beta0_Kp,
    seed = seed,
    settings = list(
      n = n,
      K = K,
      p = p,
      rho = rho,
      beta0 = beta0,
      treat_prob = treat_prob,
      balance = balance,
      nui_app = nui_app,
      a = a
    )
  )
}