# ============================================================
# Available simulation methods
# ============================================================

get_available_simu_methods <- function() {
  c(
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
    "Loc-Lasso",
    "Pooled-Lasso"
  )
}


# ============================================================
# Validate the user-specified method vector
# ============================================================

validate_simu_methods <- function(methods = NULL) {

  available <- get_available_simu_methods()

  if (is.null(methods)) {
    return(available)
  }

  methods <- as.character(methods)
  methods <- unique(methods)

  if (length(methods) == 0) {
    stop(
      "methods must contain at least one method.",
      call. = FALSE
    )
  }

  unknown <- setdiff(methods, available)

  if (length(unknown) > 0) {
    stop(
      "Unknown methods: ",
      paste(unknown, collapse = ", "),
      ". Available methods are: ",
      paste(available, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  methods
}


# ============================================================
# Row-bind data frames while allowing the first object to be NULL
# ============================================================

append_df <- function(old, new) {
  if (is.null(old)) {
    rownames(new) <- NULL
    return(new)
  }

  out <- rbind(old, new)
  rownames(out) <- NULL
  out
}


# ============================================================
# Summarize evaluation results in the requested method order
# ============================================================

summarize_evaluation <- function(all_eval, method_order) {

  metric_names <- c(
    "MSE",
    "Var_TPR",
    "Var_FDR",
    "Homo_TPR",
    "Homo_FDR"
  )

  out_list <- lapply(method_order, function(mm) {

    dat_mm <- all_eval[
      as.character(all_eval$method) == mm,
      ,
      drop = FALSE
    ]

    if (nrow(dat_mm) == 0) {
      return(NULL)
    }

    ans <- data.frame(
      method = mm,
      stringsAsFactors = FALSE
    )

    for (metric in metric_names) {
      x <- dat_mm[[metric]]

      ans[[metric]] <- if (all(is.na(x))) {
        NA_real_
      } else {
        mean(x, na.rm = TRUE)
      }
    }

    ans
  })

  out_list <- Filter(Negate(is.null), out_list)

  if (length(out_list) == 0) {
    return(data.frame())
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}


# ============================================================
# Summarize dimension-specific homogeneity decisions
# ============================================================

summarize_homo_dims <- function(all_homo_dims, method_order) {

  if (is.null(all_homo_dims) || nrow(all_homo_dims) == 0) {
    return(data.frame())
  }

  homo_cols <- grep(
    "^Homo_x",
    colnames(all_homo_dims),
    value = TRUE
  )

  if (length(homo_cols) == 0) {
    return(data.frame())
  }

  out_list <- lapply(method_order, function(mm) {

    dat_mm <- all_homo_dims[
      as.character(all_homo_dims$method) == mm,
      ,
      drop = FALSE
    ]

    if (nrow(dat_mm) == 0) {
      return(NULL)
    }

    ans <- data.frame(
      method = mm,
      stringsAsFactors = FALSE
    )

    for (cc in homo_cols) {
      x <- dat_mm[[cc]]

      ans[[cc]] <- if (all(is.na(x))) {
        NA_real_
      } else {
        mean(x, na.rm = TRUE)
      }
    }

    ans
  })

  out_list <- Filter(Negate(is.null), out_list)

  if (length(out_list) == 0) {
    return(data.frame())
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}


# ============================================================
# Prepare data when nuisance functions are known
#
# Uses:
#   m_true(X,G) = E[Y | X,G]
#   e0(X,G)     = P(A = 1 | X,G)
#
# Output has the same structure as prepare_dat().
# ============================================================

prepare_dat_known_nuisance <- function(dat, K, trim_e = 1e-4) {

  if (!all(c("Y", "A", "G") %in% colnames(dat))) {
    stop(
      "dat must contain columns Y, A, and G.",
      call. = FALSE
    )
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
    stop(
      "dat does not contain attr(dat, 'm_true'). ",
      "Please update gdat.R.",
      call. = FALSE
    )
  }

  if (is.null(e0)) {
    stop(
      "dat does not contain attr(dat, 'e0'). ",
      "Please update gdat.R.",
      call. = FALSE
    )
  }

  m_true <- as.numeric(m_true)
  e0 <- as.numeric(e0)

  if (length(m_true) != n) {
    stop(
      "Length of m_true does not match sample size.",
      call. = FALSE
    )
  }

  if (length(e0) != n) {
    stop(
      "Length of e0 does not match sample size.",
      call. = FALSE
    )
  }

  e0 <- pmin(pmax(e0, trim_e), 1 - trim_e)

  Ytilde <- Y - m_true
  Atilde <- A - e0

  if (any(!is.finite(Ytilde))) {
    stop(
      "Non-finite Ytilde in prepare_dat_known_nuisance().",
      call. = FALSE
    )
  }

  if (any(!is.finite(Atilde))) {
    stop(
      "Non-finite Atilde in prepare_dat_known_nuisance().",
      call. = FALSE
    )
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


# ============================================================
# Run one simulation
#
# methods:
#   Character vector containing any subset of
#   get_available_simu_methods().
# ============================================================

once_simu <- function(n, K, p, rho, beta0,
                      treat_prob, balance, seed,
                      nui_app = "bagging",
                      methods = NULL) {

  methods <- validate_simu_methods(methods)

  set.seed(seed)

  beta0 <- as.matrix(beta0)

  if (!(nrow(beta0) == p && ncol(beta0) == K)) {
    stop(
      "beta0 must be a p x K matrix. Expected ",
      p, " x ", K, ", got ",
      nrow(beta0), " x ", ncol(beta0), ".",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # 1. Generate data
  # ----------------------------------------------------------
  dat <- gdat(
    n = n,
    K = K,
    p = p,
    rho = rho,
    beta0 = beta0,
    treat_prob = treat_prob,
    balance = balance
  )

  # ----------------------------------------------------------
  # 2. Determine which intermediate objects are needed
  # ----------------------------------------------------------
  known_methods <- c(
    "Int-Oracle-KN",
    "Int-Oracle-KN-w",
    "Int-Prop-KN",
    "Int-Prop-KN-w"
  )

  estimated_methods <- setdiff(
    get_available_simu_methods(),
    known_methods
  )

  init_methods <- c(
    "Int-Prop-KN",
    "Int-Prop-KN-w",
    "Int-Prop",
    "Int-Prop-w",
    "Int-R2Lasso",
    "Int-R2Lasso-w",
    "Int-GLasso",
    "Int-GLasso-w",
    "Int-HLasso",
    "Int-HLasso-w"
  )

  weighted_estimated_methods <- c(
    "Int-Oracle-w",
    "Int-Prop-w",
    "Int-R2Lasso-w",
    "Int-GLasso-w",
    "Int-HLasso-w"
  )

  weighted_known_methods <- c(
    "Int-Oracle-KN-w",
    "Int-Prop-KN-w"
  )

  need_known <- any(methods %in% known_methods)

  # Estimated nuisance data are also used to construct the common
  # initial estimator for the integrated methods. This preserves
  # the behavior of the original simulation code and ensures that
  # a method does not change when other methods are omitted.
  need_estimated <- any(methods %in% estimated_methods) ||
    any(methods %in% init_methods)

  need_beta_init <- any(methods %in% init_methods)

  need_wloss <- any(methods %in% weighted_estimated_methods)
  need_wloss_kn <- any(methods %in% weighted_known_methods)

  mydat <- NULL
  mydat_kn <- NULL

  fit_wloss <- NULL
  fit_wloss_kn <- NULL

  wloss <- NULL
  wloss_kn <- NULL

  fit_loc_init <- NULL
  betahat_ini <- NULL

  # ----------------------------------------------------------
  # 3. Estimated nuisance version
  # ----------------------------------------------------------
  if (need_estimated) {

    if (identical(nui_app, "rf")) {

      mydat <- prepare_dat(
        dat,
        K = K,
        nfold = 10,
        outcome_method = "rf",
        propensity_method = "rf",
        outcome_args = list(ntree = 150),
        propensity_args = list(ntree = 150),
        seed = seed + 1000,
        trim_e = 1e-4,
        use_group_for_nuisance = TRUE,
        add_group_interactions = FALSE,
        calibrate_propensity = TRUE
      )

    } else if (identical(nui_app, "bagging")) {

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

    } else {
      stop(
        "nui_app must be either 'bagging' or 'rf'.",
        call. = FALSE
      )
    }
  }

  # ----------------------------------------------------------
  # 4. Known nuisance version
  # ----------------------------------------------------------
  if (need_known) {
    mydat_kn <- prepare_dat_known_nuisance(
      dat = dat,
      K = K,
      trim_e = 1e-4
    )
  }

  # ----------------------------------------------------------
  # 5. Initial estimator required by integrated methods
  # ----------------------------------------------------------
  if (need_beta_init) {

    if (is.null(mydat)) {
      stop(
        "Estimated-nuisance data are required to construct ",
        "the initial estimator.",
        call. = FALSE
      )
    }

    fit_loc_init <- method_loc_lasso_ols_refit(mydat)
    betahat_ini <- fit_loc_init$beta
  }

  # ----------------------------------------------------------
  # 6. Cross-fitted loss weights: estimated nuisance
  # ----------------------------------------------------------
  if (need_wloss) {

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
  }

  # ----------------------------------------------------------
  # 7. Cross-fitted loss weights: known nuisance
  # ----------------------------------------------------------
  if (need_wloss_kn) {

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
  }

  # ----------------------------------------------------------
  # 8. Fit only the requested methods
  # ----------------------------------------------------------
  beta_hat_list <- list()
  fit_list <- list()

  add_result <- function(method_name, fit_object, threshold = 0.01) {

    beta_hat <- as.matrix(fit_object$beta)

    beta_hat <- beta_hat *
      (abs(beta_hat) > threshold)

    beta_hat_list[[method_name]] <<- beta_hat
    fit_list[[method_name]] <<- fit_object

    invisible(NULL)
  }

  # ----------------------------------------------------------
  # 8.1 Known-nuisance Oracle
  # ----------------------------------------------------------
  if ("Int-Oracle-KN" %in% methods) {
    fit <- method_oracle(
      mydat_kn,
      beta0 = t(beta0),
      wloss = NULL
    )
    add_result("Int-Oracle-KN", fit)
  }

  if ("Int-Oracle-KN-w" %in% methods) {
    fit <- method_oracle(
      mydat_kn,
      beta0 = t(beta0),
      wloss = wloss_kn
    )
    add_result("Int-Oracle-KN-w", fit)
  }

  # ----------------------------------------------------------
  # 8.2 Estimated-nuisance Oracle
  # ----------------------------------------------------------
  if ("Int-Oracle" %in% methods) {
    fit <- method_oracle(
      mydat,
      beta0 = t(beta0),
      wloss = NULL
    )
    add_result("Int-Oracle", fit)
  }

  if ("Int-Oracle-w" %in% methods) {
    fit <- method_oracle(
      mydat,
      beta0 = t(beta0),
      wloss = wloss
    )
    add_result("Int-Oracle-w", fit)
  }

  # ----------------------------------------------------------
  # 8.3 Known-nuisance proposed estimator
  # ----------------------------------------------------------
  if ("Int-Prop-KN" %in% methods) {
    fit <- method_prop_admm(
      mydat_kn,
      beta_init = betahat_ini,
      nlambda1 = 25,
      nlambda2 = 15,
      lambda.min.ratio = 1e-4,
      wloss = NULL
    )
    add_result("Int-Prop-KN", fit)
  }

  if ("Int-Prop-KN-w" %in% methods) {
    fit <- method_prop_admm(
      mydat_kn,
      beta_init = betahat_ini,
      nlambda1 = 25,
      nlambda2 = 15,
      lambda.min.ratio = 1e-4,
      wloss = wloss_kn
    )
    add_result("Int-Prop-KN-w", fit)
  }

  # ----------------------------------------------------------
  # 8.4 Estimated-nuisance proposed estimator
  # ----------------------------------------------------------
  if ("Int-Prop" %in% methods) {
    fit <- method_prop_admm(
      mydat,
      beta_init = betahat_ini,
      nlambda1 = 25,
      nlambda2 = 15,
      lambda.min.ratio = 1e-4,
      wloss = NULL
    )
    add_result("Int-Prop", fit)
  }

  if ("Int-Prop-w" %in% methods) {
    fit <- method_prop_admm(
      mydat,
      beta_init = betahat_ini,
      nlambda1 = 25,
      nlambda2 = 15,
      lambda.min.ratio = 1e-4,
      wloss = wloss
    )
    add_result("Int-Prop-w", fit)
  }

  # ----------------------------------------------------------
  # 8.5 R2Lasso
  # ----------------------------------------------------------
  if ("Int-R2Lasso" %in% methods) {
    fit <- method_R2lasso(
      mydat,
      beta_init = betahat_ini,
      wloss = NULL,
      nlambda1 = 25,
      nlambda2 = 15,
      lambda.min.ratio = 1e-4
    )
    add_result("Int-R2Lasso", fit)
  }

  if ("Int-R2Lasso-w" %in% methods) {
    fit <- method_R2lasso(
      mydat,
      beta_init = betahat_ini,
      wloss = wloss,
      nlambda1 = 25,
      nlambda2 = 15,
      lambda.min.ratio = 1e-4
    )
    add_result("Int-R2Lasso-w", fit)
  }

  # ----------------------------------------------------------
  # 8.6 Group Lasso
  # ----------------------------------------------------------
  if ("Int-GLasso" %in% methods) {
    fit <- method_glasso(
      mydat,
      beta_init = betahat_ini,
      wloss = NULL,
      lambda.min.ratio = 1e-4
    )
    add_result("Int-GLasso", fit)
  }

  if ("Int-GLasso-w" %in% methods) {
    fit <- method_glasso(
      mydat,
      beta_init = betahat_ini,
      wloss = wloss,
      lambda.min.ratio = 1e-4
    )
    add_result("Int-GLasso-w", fit)
  }

  # ----------------------------------------------------------
  # 8.7 Homogeneous Lasso
  # ----------------------------------------------------------
  if ("Int-HLasso" %in% methods) {
    fit <- method_hlasso(
      mydat,
      beta_init = betahat_ini,
      wloss = NULL,
      nlambda = 100,
      lambda.min.ratio = 1e-4
    )
    add_result("Int-HLasso", fit)
  }

  if ("Int-HLasso-w" %in% methods) {
    fit <- method_hlasso(
      mydat,
      beta_init = betahat_ini,
      wloss = wloss,
      nlambda = 60,
      lambda.min.ratio = 1e-4
    )
    add_result("Int-HLasso-w", fit)
  }

  # ----------------------------------------------------------
  # 8.8 Separate local adaptive Lasso
  # ----------------------------------------------------------
  if ("Loc-Lasso" %in% methods) {

    # If the initial refit has already been computed, its stored
    # beta_lasso is exactly the corresponding local Lasso fit.
    if (!is.null(fit_loc_init) &&
        !is.null(fit_loc_init$beta_lasso)) {

      fit <- fit_loc_init
      fit$beta <- fit_loc_init$beta_lasso

    } else {
      fit <- method_loc_lasso(mydat)
    }

    add_result("Loc-Lasso", fit)
  }

  # ----------------------------------------------------------
  # 8.9 Pooled adaptive Lasso
  # ----------------------------------------------------------
  if ("Pooled-Lasso" %in% methods) {
    fit <- method_pooled_lasso(
      mydat = mydat,
      wloss = NULL,
      gamma = 1,
      nlambda = 50,
      lambda.min.ratio = 1e-4
    )
    add_result("Pooled-Lasso", fit)
  }

  # Restore exactly the order specified by the user.
  beta_hat_list <- beta_hat_list[methods]
  fit_list <- fit_list[methods]

  # ----------------------------------------------------------
  # 9. Evaluate
  # ----------------------------------------------------------
  ev <- eval_methods(
    beta_hat_list = beta_hat_list,
    beta0 = t(beta0)
  )

  list(
    beta_hat = beta_hat_list,
    eval = ev,
    fits = fit_list,

    methods = methods,
    dat = dat,
    mydat = mydat,
    mydat_kn = mydat_kn,

    beta_init = betahat_ini,
    fit_loc_init = fit_loc_init,

    wloss = wloss,
    wloss_kn = wloss_kn,

    wloss_fit = fit_wloss,
    wloss_fit_kn = fit_wloss_kn
  )
}


# ============================================================
# Main simulation function
# ============================================================

simu <- function(n, K, p, rho, beta0,
                 treat_prob, balance, nui_app,
                 seedstart, Nsimu, taskID,
                 methods = NULL) {

  methods <- validate_simu_methods(methods)
  method_order <- methods

  if (length(Nsimu) != 1 ||
      !is.finite(Nsimu) ||
      Nsimu < 1) {
    stop(
      "Nsimu must be a positive integer.",
      call. = FALSE
    )
  }

  Nsimu <- as.integer(Nsimu)

  set.seed(seedstart)
  seednum <- sample.int(10^6, Nsimu)

  ALL_EVAL <- NULL
  ALL_BHAT_WIDE <- NULL
  ALL_HOMO_DIMS <- NULL

  SUM <- NULL
  HOMO_DIM_SUM <- NULL

  # Number of active covariate dimensions. This follows the
  # original simulation design, in which active variables occur
  # in the first supp_dim positions.
  supp_dim <- sum(rowSums(as.matrix(beta0)^2) > 0)

  if (supp_dim > 0) {
    homo_dims_to_save <- seq_len(supp_dim)
  } else {
    homo_dims_to_save <- integer(0)
  }

  for (b in seq_len(Nsimu)) {

    cat(
      "\nSimulation:", b, "/", Nsimu,
      "  Seed:", seednum[b], "\n"
    )

    cat(
      "Methods:",
      paste(methods, collapse = ", "),
      "\n"
    )

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
      methods = methods
    )

    # --------------------------------------------------------
    # 1. Save evaluation results
    # --------------------------------------------------------
    ev <- one$eval
    ev$simu <- b
    ev$seed <- seednum[b]

    # Put identifiers first.
    ev <- ev[
      ,
      c(
        "simu",
        "seed",
        setdiff(colnames(ev), c("simu", "seed"))
      ),
      drop = FALSE
    ]

    ALL_EVAL <- append_df(ALL_EVAL, ev)

    # --------------------------------------------------------
    # 2. Save beta estimates in wide format
    # --------------------------------------------------------
    bhat_wide <- beta_hat_to_wide(
      beta_hat_list = one$beta_hat,
      simu = b,
      seed = seednum[b]
    )

    ALL_BHAT_WIDE <- append_df(
      ALL_BHAT_WIDE,
      bhat_wide
    )

    # --------------------------------------------------------
    # 3. Save dimension-specific homogeneity decisions
    # --------------------------------------------------------
    if (length(homo_dims_to_save) > 0) {

      homo_dims <- eval_homo_first_dims(
        beta_hat_list = one$beta_hat,
        dims = homo_dims_to_save,
        simu = b,
        seed = seednum[b],
        tol = 1e-6
      )

      ALL_HOMO_DIMS <- append_df(
        ALL_HOMO_DIMS,
        homo_dims
      )
    }

    # --------------------------------------------------------
    # 4. Write raw results after every replication
    # --------------------------------------------------------
    write.csv(
      ALL_EVAL,
      paste0(taskID, "_AllEval.csv"),
      row.names = FALSE
    )

    write.csv(
      ALL_BHAT_WIDE,
      paste0(taskID, "_AllBhat_wide.csv"),
      row.names = FALSE
    )

    if (!is.null(ALL_HOMO_DIMS)) {
      write.csv(
        ALL_HOMO_DIMS,
        paste0(taskID, "_AllHomoDims.csv"),
        row.names = FALSE
      )
    }

    # --------------------------------------------------------
    # 5. Evaluation summary
    # --------------------------------------------------------
    SUM <- summarize_evaluation(
      all_eval = ALL_EVAL,
      method_order = method_order
    )

    # --------------------------------------------------------
    # 6. Homogeneity summary
    # --------------------------------------------------------
    HOMO_DIM_SUM <- summarize_homo_dims(
      all_homo_dims = ALL_HOMO_DIMS,
      method_order = method_order
    )

    if (nrow(HOMO_DIM_SUM) > 0) {

      write.csv(
        HOMO_DIM_SUM,
        paste0(taskID, "_HomoDimsSummary.csv"),
        row.names = FALSE
      )

      SUM <- merge(
        SUM,
        HOMO_DIM_SUM,
        by = "method",
        all.x = TRUE,
        sort = FALSE
      )

      SUM <- SUM[
        match(method_order, SUM$method),
        ,
        drop = FALSE
      ]

      SUM <- SUM[
        !is.na(SUM$method),
        ,
        drop = FALSE
      ]
    }

    # --------------------------------------------------------
    # 7. Save summary
    # --------------------------------------------------------
    write.csv(
      SUM,
      paste0(taskID, "_Summary.csv"),
      row.names = FALSE
    )

    # --------------------------------------------------------
    # 8. Save RDS
    # --------------------------------------------------------
    saveRDS(
      list(
        ALL_EVAL = ALL_EVAL,
        ALL_BHAT_WIDE = ALL_BHAT_WIDE,
        ALL_HOMO_DIMS = ALL_HOMO_DIMS,
        HOMO_DIM_SUM = HOMO_DIM_SUM,
        SUMMARY = SUM,
        seednum = seednum,
        completed_simulations = b,

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
          methods = methods
        )
      ),
      paste0(taskID, "_AllResults.rds")
    )

    print(SUM)
  }

  invisible(
    list(
      ALL_EVAL = ALL_EVAL,
      ALL_BHAT_WIDE = ALL_BHAT_WIDE,
      ALL_HOMO_DIMS = ALL_HOMO_DIMS,
      HOMO_DIM_SUM = HOMO_DIM_SUM,
      SUMMARY = SUM,
      seednum = seednum,
      methods = methods
    )
  )
}
