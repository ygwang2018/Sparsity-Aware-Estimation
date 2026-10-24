# ============================================================
# simu_nuisance_compare.R
#
# Compare different nuisance estimators for proposed method.
# Only output:
#   Prop and Oracle, with and without loss weights.
#
# Nuisance estimators:
#   known
#   bagging
#   rf
#   glm_lasso
#
# Modeling modes:
#   pooled: use G as nuisance covariate
#   byG:    fit nuisance models separately within each G
# ============================================================

# ------------------------------------------------------------
# Known nuisance
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
    stop("dat does not contain attr(dat, 'm_true').", call. = FALSE)
  }

  if (is.null(e0)) {
    stop("dat does not contain attr(dat, 'e0').", call. = FALSE)
  }

  m_true <- as.numeric(m_true)
  e0 <- as.numeric(e0)

  e0 <- pmin(pmax(e0, trim_e), 1 - trim_e)

  Ytilde <- Y - m_true
  Atilde <- A - e0

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
    nuisance_mode = "known",

    outcome_args = NULL,
    propensity_args = NULL,

    use_group_for_nuisance = TRUE,
    add_group_interactions = FALSE,
    calibrate_propensity = FALSE,
    trim_e = trim_e
  )
}

# ------------------------------------------------------------
# RF regression for A: E[A | X, G]
# Not classification probability.
# ------------------------------------------------------------
fit_predict_rf_reg_A <- function(Xtrain, atrain, Xtest, args = NULL) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    warning("Package 'randomForest' is not available. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)
  colnames(Xtest) <- colnames(Xtrain)

  atrain <- as.numeric(atrain)

  if (length(unique(atrain)) < 2) {
    return(rep(mean(atrain, na.rm = TRUE), nrow(Xtest)))
  }

  if (is.null(args)) args <- list()

  fit <- tryCatch(
    do.call(
      randomForest::randomForest,
      c(
        list(
          x = as.data.frame(Xtrain, check.names = FALSE),
          y = atrain
        ),
        args
      )
    ),
    error = function(e) {
      warning("RF regression for A failed. Fallback to glm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred <- tryCatch(
    stats::predict(fit, newdata = as.data.frame(Xtest, check.names = FALSE)),
    error = function(e) {
      warning("RF regression prediction for A failed. Fallback to glm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(pred)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred <- as.numeric(pred)

  if (length(pred) != nrow(Xtest) || any(!is.finite(pred))) {
    warning("RF regression prediction for A invalid. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred
}

# ------------------------------------------------------------
# Bagging regression for A: E[A | X, G]
# If you already source prepare_reg_bag.R, this is consistent.
# ------------------------------------------------------------
fit_predict_bag_reg_A <- function(Xtrain, atrain, Xtest, args = NULL) {
  if (!requireNamespace("ipred", quietly = TRUE)) {
    warning("Package 'ipred' is not available. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)
  colnames(Xtest) <- colnames(Xtrain)

  atrain <- as.numeric(atrain)

  if (length(unique(atrain)) < 2) {
    return(rep(mean(atrain, na.rm = TRUE), nrow(Xtest)))
  }

  if (is.null(args)) args <- list()

  df_train <- data.frame(
    a = atrain,
    as.data.frame(Xtrain, check.names = FALSE)
  )

  df_test <- as.data.frame(Xtest, check.names = FALSE)

  fit <- tryCatch(
    do.call(
      ipred::bagging,
      c(
        list(
          formula = a ~ . - 1,
          data = df_train
        ),
        args
      )
    ),
    error = function(e) {
      warning("Bagging regression for A failed. Fallback to glm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred <- tryCatch(
    stats::predict(fit, newdata = df_test),
    error = function(e) {
      warning("Bagging regression prediction for A failed. Fallback to glm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(pred)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred <- as.numeric(pred)

  if (length(pred) != nrow(Xtest) || any(!is.finite(pred))) {
    warning("Bagging regression prediction for A invalid. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred
}

# ------------------------------------------------------------
# glm + lasso for continuous outcome/regression target
# Used for both Y and numeric A.
# ------------------------------------------------------------
fit_predict_glm_lasso_reg <- function(Xtrain, ytrain, Xtest, args = NULL) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    warning("Package 'glmnet' is not available. Fallback to lm.")
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)
  colnames(Xtest) <- colnames(Xtrain)

  ytrain <- as.numeric(ytrain)

  if (is.null(args)) args <- list()

  alpha <- ifelse(is.null(args$alpha), 1, args$alpha)
  nfolds <- ifelse(is.null(args$nfolds), 5, args$nfolds)

  fit <- tryCatch(
    glmnet::cv.glmnet(
      x = Xtrain,
      y = ytrain,
      family = "gaussian",
      alpha = alpha,
      nfolds = min(nfolds, nrow(Xtrain)),
      intercept = TRUE,
      standardize = TRUE
    ),
    error = function(e) {
      warning("glmnet lasso regression failed. Fallback to lm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred <- tryCatch(
    as.numeric(stats::predict(fit, newx = Xtest, s = "lambda.min")),
    error = function(e) {
      warning("glmnet lasso prediction failed. Fallback to lm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(pred) || length(pred) != nrow(Xtest) || any(!is.finite(pred))) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred
}

# ------------------------------------------------------------
# One unified prediction interface
# ------------------------------------------------------------
fit_predict_nuisance_reg <- function(method, Xtrain, ytrain, Xtest, args = NULL,
                                     target = c("Y", "A")) {
  target <- match.arg(target)

  if (method %in% c("lm", "linear")) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args))
  }

  if (method %in% c("glm_lasso", "lasso", "glmnet")) {
    return(fit_predict_glm_lasso_reg(Xtrain, ytrain, Xtest, args))
  }

  if (method %in% c("rf", "randomForest", "randomforest")) {
    if (target == "A") {
      return(fit_predict_rf_reg_A(Xtrain, ytrain, Xtest, args))
    } else {
      return(fit_predict_rf_reg(Xtrain, ytrain, Xtest, args))
    }
  }

  if (method %in% c("bagging", "bag")) {
    if (target == "A") {
      return(fit_predict_bag_reg_A(Xtrain, ytrain, Xtest, args))
    } else {
      return(fit_predict_bag_reg(Xtrain, ytrain, Xtest, args))
    }
  }

  stop("Unknown nuisance method: ", method, call. = FALSE)
}

# ------------------------------------------------------------
# Prepare data for nuisance comparison
#
# nuisance_mode:
#   pooled: include G in nuisance covariates
#   byG:    fit nuisance models separately within each G
# ------------------------------------------------------------
prepare_dat_nuisance_compare <- function(dat,
                                         K,
                                         nuisance_method = c("known", "bagging", "rf", "glm_lasso"),
                                         nuisance_mode = c("pooled", "byG"),
                                         nfold = 5,
                                         outcome_args = NULL,
                                         propensity_args = NULL,
                                         seed = 1,
                                         trim_e = 1e-4,
                                         add_group_interactions = FALSE,
                                         calibrate_propensity = TRUE) {

  nuisance_method <- match.arg(nuisance_method)
  nuisance_mode <- match.arg(nuisance_mode)

  if (nuisance_method == "known") {
    out <- prepare_dat_known_nuisance(dat, K = K, trim_e = trim_e)
    out$nuisance_method <- "known"
    out$nuisance_mode <- "known"
    return(out)
  }

  set.seed(seed)

  Y <- as.numeric(dat$Y)
  A <- as.numeric(dat$A)
  G <- as.integer(dat$G)

  X <- get_X_from_dat(dat)
  X <- clean_matrix(X)

  n <- nrow(X)
  p <- ncol(X)

  if (nfold < 2) stop("nfold must be at least 2.", call. = FALSE)
  if (nfold > n) stop("nfold cannot exceed n.", call. = FALSE)

  foldid <- sample(rep(seq_len(nfold), length.out = n))

  mhat <- rep(NA_real_, n)
  ehat <- rep(NA_real_, n)

  # ----------------------------------------------------------
  # Pooled nuisance: G enters as predictor.
  # ----------------------------------------------------------
  if (nuisance_mode == "pooled") {

    Xnuis <- build_nuisance_X(
      X = X,
      G = G,
      K = K,
      add_group_interactions = add_group_interactions
    )

    Xnuis <- clean_matrix(Xnuis)

    for (ff in seq_len(nfold)) {
      idx_test <- which(foldid == ff)
      idx_train <- setdiff(seq_len(n), idx_test)

      Xtr <- Xnuis[idx_train, , drop = FALSE]
      Xte <- Xnuis[idx_test, , drop = FALSE]

      Ytr <- Y[idx_train]
      Atr <- A[idx_train]
      Gtr <- G[idx_train]
      Gte <- G[idx_test]

      mhat[idx_test] <- fit_predict_nuisance_reg(
        method = nuisance_method,
        Xtrain = Xtr,
        ytrain = Ytr,
        Xtest = Xte,
        args = outcome_args,
        target = "Y"
      )

      ehat_ff <- fit_predict_nuisance_reg(
        method = nuisance_method,
        Xtrain = Xtr,
        ytrain = Atr,
        Xtest = Xte,
        args = propensity_args,
        target = "A"
      )

      ehat_ff <- pmin(pmax(ehat_ff, trim_e), 1 - trim_e)

      if (isTRUE(calibrate_propensity)) {
        ehat_ff <- calibrate_ehat_by_train_group(
          ehat_test = ehat_ff,
          Atrain = Atr,
          Gtrain = Gtr,
          Gtest = Gte,
          K = K,
          trim_e = trim_e
        )
      }

      ehat[idx_test] <- ehat_ff
    }
  }

  # ----------------------------------------------------------
  # byG nuisance: fit a separate nuisance model within each G.
  # G is not used as a covariate.
  # ----------------------------------------------------------
  if (nuisance_mode == "byG") {

    for (k in seq_len(K)) {
      idx_k <- which(G == k)
      nk <- length(idx_k)

      if (nk == 0) next

      Lk <- min(nfold, nk)

      if (Lk < 2) {
        mhat[idx_k] <- mean(Y[idx_k])
        ehat[idx_k] <- mean(A[idx_k])
        next
      }

      fold_k <- sample(rep(seq_len(Lk), length.out = nk))
      foldid[idx_k] <- fold_k

      for (ff in seq_len(Lk)) {
        idx_test <- idx_k[fold_k == ff]
        idx_train <- idx_k[fold_k != ff]

        Xtr <- X[idx_train, , drop = FALSE]
        Xte <- X[idx_test, , drop = FALSE]

        Ytr <- Y[idx_train]
        Atr <- A[idx_train]

        mhat[idx_test] <- fit_predict_nuisance_reg(
          method = nuisance_method,
          Xtrain = Xtr,
          ytrain = Ytr,
          Xtest = Xte,
          args = outcome_args,
          target = "Y"
        )

        ehat_ff <- fit_predict_nuisance_reg(
          method = nuisance_method,
          Xtrain = Xtr,
          ytrain = Atr,
          Xtest = Xte,
          args = propensity_args,
          target = "A"
        )

        ehat_ff <- pmin(pmax(ehat_ff, trim_e), 1 - trim_e)

        # byG 中每个 G 单独建模，这里可用该 G training fold 的 A 均值做校准
        if (isTRUE(calibrate_propensity)) {
          target <- mean(Atr)
          target <- pmin(pmax(target, trim_e), 1 - trim_e)

          current <- mean(ehat_ff)
          current <- pmin(pmax(current, trim_e), 1 - trim_e)

          delta <- stats::qlogis(target) - stats::qlogis(current)
          ehat_ff <- stats::plogis(stats::qlogis(ehat_ff) + delta)
          ehat_ff <- pmin(pmax(ehat_ff, trim_e), 1 - trim_e)
        }

        ehat[idx_test] <- ehat_ff
      }
    }
  }

  if (any(!is.finite(mhat))) {
    stop("Non-finite mhat in prepare_dat_nuisance_compare().", call. = FALSE)
  }

  if (any(!is.finite(ehat))) {
    stop("Non-finite ehat in prepare_dat_nuisance_compare().", call. = FALSE)
  }

  ehat <- pmin(pmax(ehat, trim_e), 1 - trim_e)

  Ytilde <- Y - mhat
  Atilde <- A - ehat

  list(
    Y = Y,
    A = A,
    G = G,
    X = X,
    K = as.integer(K),
    p = p,
    n = n,

    mhat = mhat,
    ehat = ehat,
    Ytilde = Ytilde,
    Atilde = Atilde,

    foldid = foldid,
    nfold = nfold,

    outcome_method = nuisance_method,
    propensity_method = nuisance_method,
    nuisance_method = nuisance_method,
    nuisance_mode = nuisance_mode,

    outcome_args = outcome_args,
    propensity_args = propensity_args,

    use_group_for_nuisance = nuisance_mode == "pooled",
    add_group_interactions = add_group_interactions,
    calibrate_propensity = calibrate_propensity,
    trim_e = trim_e
  )
}

# ------------------------------------------------------------
# Fit Prop and Oracle for one prepared mydat
# ------------------------------------------------------------
fit_prop_oracle_one_nuisance <- function(mydat,
                                         beta0_Kp,
                                         method_prefix,
                                         seed = 1,
                                         use_weighted = TRUE,
                                         beta_init = NULL,
                                         threshold = 0.01) {

  if (is.null(beta_init)) {
    fit_loc_ols <- method_loc_lasso_ols_refit(mydat)
    beta_init <- fit_loc_ols$beta
  }

  beta_hat_list <- list()
  wloss_fit <- NULL
  wloss <- NULL

  if (isTRUE(use_weighted)) {
    wloss_fit <- make_wloss_cf_local(
      mydat = mydat,
      nfold = 10,
      seed = seed + 5000,
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

    wloss <- wloss_fit$wloss
  }

  # Prop, unweighted
  fit_prop <- method_prop_admm(
    mydat,
    beta_init = beta_init,
    nlambda1 = 25,
    nlambda2 = 15,
    lambda.min.ratio = 1e-4,
    wloss = NULL
  )

  # Oracle, unweighted
  fit_oracle <- method_oracle(
    mydat,
    beta0 = beta0_Kp,
    wloss = NULL
  )

  beta_hat_list[[paste0("Prop-", method_prefix)]] <-
    fit_prop$beta * (abs(fit_prop$beta) > threshold)

  beta_hat_list[[paste0("Oracle-", method_prefix)]] <-
    fit_oracle$beta * (abs(fit_oracle$beta) > threshold)

  # Weighted versions
  if (isTRUE(use_weighted)) {
    fit_prop_w <- method_prop_admm(
      mydat,
      beta_init = beta_init,
      nlambda1 = 25,
      nlambda2 = 15,
      lambda.min.ratio = 1e-4,
      wloss = wloss
    )

    fit_oracle_w <- method_oracle(
      mydat,
      beta0 = beta0_Kp,
      wloss = wloss
    )

    beta_hat_list[[paste0("Prop-", method_prefix, "-w")]] <-
      fit_prop_w$beta * (abs(fit_prop_w$beta) > threshold)

    beta_hat_list[[paste0("Oracle-", method_prefix, "-w")]] <-
      fit_oracle_w$beta * (abs(fit_oracle_w$beta) > threshold)
  }

  list(
    beta_hat = beta_hat_list,
    wloss = wloss,
    wloss_fit = wloss_fit
  )
}

# ------------------------------------------------------------
# One simulation for nuisance comparison
# ------------------------------------------------------------
once_simu_nuisance_compare <- function(n,
                                       K,
                                       p,
                                       rho,
                                       beta0,
                                       treat_prob,
                                       balance,
                                       seed,
                                       nuisance_methods = c("known", "bagging", "rf", "glm_lasso"),
                                       nuisance_modes = c("pooled", "byG"),
                                       nfold = 10,
                                       use_weighted = TRUE) {

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

  beta0_Kp <- t(beta0)

  beta_hat_all <- list()
  mydat_all <- list()
  wloss_all <- list()

  for (mm in nuisance_methods) {

    if (mm == "known") {
      mode_vec <- "known"
    } else {
      mode_vec <- nuisance_modes
    }

    for (mode in mode_vec) {

      cat("  nuisance:", mm, ", mode:", mode, "\n")

      # method label
      prefix <- if (mm == "known") {
        "known"
      } else {
        paste0(mm, "-", mode)
      }

      # arguments
      if (mm == "bagging") {
        outcome_args <- list(nbagg = 50)
        propensity_args <- list(nbagg = 50)
      } else if (mm == "rf") {
        outcome_args <- list(ntree = 150)
        propensity_args <- list(ntree = 150)
      } else if (mm == "glm_lasso") {
        outcome_args <- list(alpha = 1, nfolds = 5)
        propensity_args <- list(alpha = 1, nfolds = 5)
      } else {
        outcome_args <- NULL
        propensity_args <- NULL
      }

      mydat <- prepare_dat_nuisance_compare(
        dat = dat,
        K = K,
        nuisance_method = mm,
        nuisance_mode = ifelse(mm == "known", "pooled", mode),
        nfold = nfold,
        outcome_args = outcome_args,
        propensity_args = propensity_args,
        seed = seed + 1000,
        trim_e = 1e-4,
        add_group_interactions = FALSE,
        calibrate_propensity = TRUE
      )

      fit_one <- fit_prop_oracle_one_nuisance(
        mydat = mydat,
        beta0_Kp = beta0_Kp,
        method_prefix = prefix,
        seed = seed,
        use_weighted = use_weighted,
        beta_init = NULL,
        threshold = 0.01
      )

      beta_hat_all <- c(beta_hat_all, fit_one$beta_hat)
      mydat_all[[prefix]] <- mydat
      wloss_all[[prefix]] <- fit_one$wloss
    }
  }

  ev <- eval_methods(beta_hat_all, beta0_Kp)

  list(
    beta_hat = beta_hat_all,
    eval = ev,
    mydat = mydat_all,
    wloss = wloss_all
  )
}

# ------------------------------------------------------------
# Main simulation wrapper for nuisance comparison
# ------------------------------------------------------------
simu_nuisance_compare <- function(n,
                                  K,
                                  p,
                                  rho,
                                  beta0,
                                  treat_prob,
                                  balance,
                                  seedstart,
                                  Nsimu,
                                  taskID,
                                  nuisance_methods = c("known", "bagging", "rf", "glm_lasso"),
                                  nuisance_modes = c("pooled", "byG"),
                                  nfold = 10,
                                  use_weighted = TRUE) {

  set.seed(seedstart)
  seednum <- sample(1:10^6, Nsimu)

  ALL_EVAL <- NULL
  ALL_BHAT_WIDE <- NULL

  method_order <- c(
    "Oracle-known",
    "Oracle-known-w",
    "Prop-known",
    "Prop-known-w",

    "Oracle-bagging-pooled",
    "Oracle-bagging-pooled-w",
    "Prop-bagging-pooled",
    "Prop-bagging-pooled-w",

    "Oracle-bagging-byG",
    "Oracle-bagging-byG-w",
    "Prop-bagging-byG",
    "Prop-bagging-byG-w",

    "Oracle-rf-pooled",
    "Oracle-rf-pooled-w",
    "Prop-rf-pooled",
    "Prop-rf-pooled-w",

    "Oracle-rf-byG",
    "Oracle-rf-byG-w",
    "Prop-rf-byG",
    "Prop-rf-byG-w",

    "Oracle-glm_lasso-pooled",
    "Oracle-glm_lasso-pooled-w",
    "Prop-glm_lasso-pooled",
    "Prop-glm_lasso-pooled-w",

    "Oracle-glm_lasso-byG",
    "Oracle-glm_lasso-byG-w",
    "Prop-glm_lasso-byG",
    "Prop-glm_lasso-byG-w"
  )

  if (!isTRUE(use_weighted)) {
    method_order <- method_order[!grepl("-w$", method_order)]
  }

  for (b in seq_len(Nsimu)) {
    cat("Simulation:", b, "/", Nsimu, "\n")

    one <- once_simu_nuisance_compare(
      n = n,
      K = K,
      p = p,
      rho = rho,
      beta0 = beta0,
      treat_prob = treat_prob,
      balance = balance,
      seed = seednum[b],
      nuisance_methods = nuisance_methods,
      nuisance_modes = nuisance_modes,
      nfold = nfold,
      use_weighted = use_weighted
    )

    ev <- one$eval
    ev$simu <- b
    ev$seed <- seednum[b]

    ALL_EVAL <- rbind(ALL_EVAL, ev)

    bhat_wide <- beta_hat_to_wide(
      beta_hat_list = one$beta_hat,
      simu = b,
      seed = seednum[b]
    )

    ALL_BHAT_WIDE <- rbind(ALL_BHAT_WIDE, bhat_wide)

    write.csv(
      ALL_EVAL,
      paste0(taskID, "_NuisanceCompare_AllEval.csv"),
      row.names = FALSE
    )

    write.csv(
      ALL_BHAT_WIDE,
      paste0(taskID, "_NuisanceCompare_AllBhat_wide.csv"),
      row.names = FALSE
    )

    ALL_EVAL$method <- factor(ALL_EVAL$method, levels = method_order)

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
      paste0(taskID, "_NuisanceCompare_Summary.csv"),
      row.names = FALSE
    )

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
          taskID = taskID,
          nuisance_methods = nuisance_methods,
          nuisance_modes = nuisance_modes,
          nfold = nfold,
          use_weighted = use_weighted
        )
      ),
      paste0(taskID, "_NuisanceCompare_AllResults.rds")
    )

    print(SUM)
  }

  invisible(list(
    ALL_EVAL = ALL_EVAL,
    ALL_BHAT_WIDE = ALL_BHAT_WIDE,
    seednum = seednum
  ))
}