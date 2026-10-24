# ============================================================
# prepare.R
#
# Prepare residualized outcome and treatment for subgroup causal
# estimation.
#
# Input dat should contain:
#   Y, A, G, X1, X2, ..., Xp
#
# Output mydat contains:
#   Ytilde = Y - mhat(X, G)
#   Atilde = A - ehat(X, G)
#   X      = original covariate matrix
#   G, A, Y, K, p, n
#
# This version completely removes xgboost.
#
# Supported nuisance methods:
#
# Outcome model m(X,G):
#   "lm", "linear", "glm"
#   "rf", "randomForest", "randomforest"
#   "bagging", "bag"
#
# Propensity model e(X,G):
#   "glm", "logit", "logistic"
#   "rf", "randomForest", "randomforest"
#   "bagging", "bag"
#
# Important:
#   ehat is always probability P(A = 1 | X, G), not class labels.
# ============================================================

# ------------------------------------------------------------
# Extract original X matrix from dat
# ------------------------------------------------------------
get_X_from_dat <- function(dat) {
  if (!all(c("Y", "A", "G") %in% colnames(dat))) {
    stop("dat must contain columns Y, A, and G.", call. = FALSE)
  }

  xcols <- setdiff(colnames(dat), c("Y", "A", "G"))

  if (length(xcols) == 0) {
    stop("No covariate columns found in dat.", call. = FALSE)
  }

  X <- dat[, xcols, drop = FALSE]
  X <- as.data.frame(X)

  for (jj in seq_along(X)) {
    X[[jj]] <- as.numeric(X[[jj]])
  }

  X <- as.matrix(X)
  storage.mode(X) <- "double"

  colnames(X) <- make.names(xcols, unique = TRUE)

  X
}

# ------------------------------------------------------------
# Clean matrix
# ------------------------------------------------------------
clean_matrix <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }

  colnames(X) <- make.names(colnames(X), unique = TRUE)

  X
}

# ------------------------------------------------------------
# Build nuisance feature matrix
#
# Nuisance models use:
#   X
#   subgroup dummies G
#   optional subgroup interactions X:G
#
# Treatment-effect models still use original X only.
# ------------------------------------------------------------
build_nuisance_X <- function(X, G, K, add_group_interactions = TRUE) {
  X <- clean_matrix(X)

  G <- as.integer(G)

  if (any(G < 1 | G > K)) {
    stop("G must contain subgroup labels from 1 to K.", call. = FALSE)
  }

  Gfac <- factor(G, levels = seq_len(K))

  Gdum <- stats::model.matrix(~ Gfac - 1)
  colnames(Gdum) <- paste0("G", seq_len(K))

  out <- cbind(X, Gdum)

  if (isTRUE(add_group_interactions)) {
    XG <- do.call(
      cbind,
      lapply(seq_len(K), function(k) {
        X * as.numeric(G == k)
      })
    )

    colnames(XG) <- as.vector(
      outer(
        colnames(X),
        paste0("G", seq_len(K)),
        FUN = function(x, g) paste0(x, "_x_", g)
      )
    )

    out <- cbind(out, XG)
  }

  out <- as.matrix(out)
  storage.mode(out) <- "double"

  colnames(out) <- make.names(colnames(out), unique = TRUE)

  out
}

# ------------------------------------------------------------
# Fold-wise subgroup calibration for propensity
#
# This calibrates test-fold ehat using only training-fold A/G.
# It applies a group-specific intercept shift on logit(ehat).
# ------------------------------------------------------------
calibrate_ehat_by_train_group <- function(ehat_test,
                                          Atrain,
                                          Gtrain,
                                          Gtest,
                                          K,
                                          trim_e = 1e-4) {
  ehat_test <- pmin(pmax(as.numeric(ehat_test), trim_e), 1 - trim_e)

  Atrain <- as.numeric(Atrain)
  Gtrain <- as.integer(Gtrain)
  Gtest <- as.integer(Gtest)

  out <- ehat_test

  global_target <- mean(Atrain)
  global_target <- pmin(pmax(global_target, trim_e), 1 - trim_e)

  for (k in seq_len(K)) {
    idx_test_k <- which(Gtest == k)

    if (length(idx_test_k) == 0) {
      next
    }

    idx_train_k <- which(Gtrain == k)

    if (length(idx_train_k) > 0) {
      target <- mean(Atrain[idx_train_k])
    } else {
      target <- global_target
    }

    target <- pmin(pmax(target, trim_e), 1 - trim_e)

    current <- mean(ehat_test[idx_test_k])
    current <- pmin(pmax(current, trim_e), 1 - trim_e)

    delta <- stats::qlogis(target) - stats::qlogis(current)

    out[idx_test_k] <- stats::plogis(stats::qlogis(ehat_test[idx_test_k]) + delta)
  }

  pmin(pmax(out, trim_e), 1 - trim_e)
}

# ------------------------------------------------------------
# Linear regression nuisance model for outcome
#
# Returns:
#   mhat = E[Y | X, G]
# ------------------------------------------------------------
fit_predict_lm_reg <- function(Xtrain, ytrain, Xtest, args = NULL) {
  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)

  colnames(Xtest) <- colnames(Xtrain)

  ytrain <- as.numeric(ytrain)

  if (length(ytrain) != nrow(Xtrain)) {
    stop("Length of ytrain must equal nrow(Xtrain).", call. = FALSE)
  }

  df_train <- data.frame(
    y = ytrain,
    as.data.frame(Xtrain, check.names = FALSE)
  )

  df_test <- as.data.frame(Xtest, check.names = FALSE)

  fit <- tryCatch(
    {
      stats::lm(y ~ . - 1, data = df_train)
    },
    error = function(e) {
      warning(
        "lm outcome model failed. Use mean outcome. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(rep(mean(ytrain, na.rm = TRUE), nrow(Xtest)))
  }

  pred <- tryCatch(
    {
      stats::predict(fit, newdata = df_test)
    },
    error = function(e) {
      warning(
        "lm outcome prediction failed. Use mean outcome. Original error: ",
        conditionMessage(e)
      )
      rep(mean(ytrain, na.rm = TRUE), nrow(Xtest))
    }
  )

  pred <- as.numeric(pred)

  if (length(pred) != nrow(Xtest)) {
    warning("lm prediction length mismatch. Use mean outcome.")
    pred <- rep(mean(ytrain, na.rm = TRUE), nrow(Xtest))
  }

  if (any(!is.finite(pred))) {
    pred[!is.finite(pred)] <- mean(ytrain, na.rm = TRUE)
  }

  pred
}

# ------------------------------------------------------------
# Logistic regression nuisance model for propensity
#
# Returns:
#   ehat = P[A = 1 | X, G]
# ------------------------------------------------------------
fit_predict_glm_bin <- function(Xtrain, atrain, Xtest, args = NULL) {
  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)

  colnames(Xtest) <- colnames(Xtrain)

  atrain <- as.numeric(atrain)

  if (length(atrain) != nrow(Xtrain)) {
    stop("Length of atrain must equal nrow(Xtrain).", call. = FALSE)
  }

  if (any(!atrain %in% c(0, 1))) {
    stop("A must be binary 0/1 for propensity model.", call. = FALSE)
  }

  if (length(unique(atrain)) < 2) {
    pred <- rep(mean(atrain, na.rm = TRUE), nrow(Xtest))
    return(pred)
  }

  df_train <- data.frame(
    a = atrain,
    as.data.frame(Xtrain, check.names = FALSE)
  )

  df_test <- as.data.frame(Xtest, check.names = FALSE)

  fit <- tryCatch(
    {
      suppressWarnings(
        stats::glm(
          a ~ . - 1,
          data = df_train,
          family = stats::binomial()
        )
      )
    },
    error = function(e) {
      warning(
        "glm propensity model failed. Use mean treatment probability. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(fit)) {
    pred <- rep(mean(atrain, na.rm = TRUE), nrow(Xtest))
  } else {
    pred <- tryCatch(
      {
        stats::predict(
          fit,
          newdata = df_test,
          type = "response"
        )
      },
      error = function(e) {
        warning(
          "glm propensity prediction failed. Use mean treatment probability. Original error: ",
          conditionMessage(e)
        )
        rep(mean(atrain, na.rm = TRUE), nrow(Xtest))
      }
    )
  }

  pred <- as.numeric(pred)

  if (length(pred) != nrow(Xtest)) {
    warning("glm propensity prediction length mismatch. Use mean treatment probability.")
    pred <- rep(mean(atrain, na.rm = TRUE), nrow(Xtest))
  }

  if (any(!is.finite(pred))) {
    pred[!is.finite(pred)] <- mean(atrain, na.rm = TRUE)
  }

  pred
}

# ------------------------------------------------------------
# Random forest regression nuisance model for outcome
#
# Returns:
#   mhat = E[Y | X, G]
# ------------------------------------------------------------
fit_predict_rf_reg <- function(Xtrain, ytrain, Xtest, args = NULL) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    warning("Package 'randomForest' is not available. Fallback to lm.")
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)

  colnames(Xtest) <- colnames(Xtrain)

  ytrain <- as.numeric(ytrain)

  if (is.null(args)) {
    args <- list()
  }

  fit <- tryCatch(
    {
      do.call(
        randomForest::randomForest,
        c(
          list(
            x = as.data.frame(Xtrain, check.names = FALSE),
            y = ytrain
          ),
          args
        )
      )
    },
    error = function(e) {
      warning(
        "randomForest outcome model failed. Fallback to lm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred <- tryCatch(
    {
      stats::predict(
        fit,
        newdata = as.data.frame(Xtest, check.names = FALSE)
      )
    },
    error = function(e) {
      warning(
        "randomForest outcome prediction failed. Fallback to lm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(pred)) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred <- as.numeric(pred)

  if (length(pred) != nrow(Xtest)) {
    warning("randomForest outcome prediction length mismatch. Fallback to lm.")
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  if (any(!is.finite(pred))) {
    warning("randomForest outcome prediction contains non-finite values. Fallback to lm.")
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred
}

# ------------------------------------------------------------
# Random forest classification nuisance model for propensity
#
# Returns:
#   ehat = P[A = 1 | X, G]
# ------------------------------------------------------------
fit_predict_rf_bin <- function(Xtrain, atrain, Xtest, args = NULL) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    warning("Package 'randomForest' is not available. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)

  colnames(Xtest) <- colnames(Xtrain)

  atrain <- as.numeric(atrain)

  if (any(!atrain %in% c(0, 1))) {
    stop("A must be binary 0/1 for propensity model.", call. = FALSE)
  }

  if (length(unique(atrain)) < 2) {
    pred <- rep(mean(atrain, na.rm = TRUE), nrow(Xtest))
    return(pred)
  }

  if (is.null(args)) {
    args <- list()
  }

  fit <- tryCatch(
    {
      do.call(
        randomForest::randomForest,
        c(
          list(
            x = as.data.frame(Xtrain, check.names = FALSE),
            y = factor(atrain, levels = c(0, 1))
          ),
          args
        )
      )
    },
    error = function(e) {
      warning(
        "randomForest propensity model failed. Fallback to glm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pr <- tryCatch(
    {
      stats::predict(
        fit,
        newdata = as.data.frame(Xtest, check.names = FALSE),
        type = "prob"
      )
    },
    error = function(e) {
      warning(
        "randomForest propensity probability prediction failed. Fallback to glm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(pr)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pr <- as.data.frame(pr)

  if ("1" %in% colnames(pr)) {
    pred <- as.numeric(pr[, "1"])
  } else if (ncol(pr) >= 2) {
    pred <- as.numeric(pr[, 2])
  } else {
    warning("randomForest propensity probability matrix has unexpected columns. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  if (length(pred) != nrow(Xtest)) {
    warning("randomForest propensity prediction length mismatch. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  if (any(!is.finite(pred))) {
    warning("randomForest propensity prediction contains non-finite values. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred
}

# ------------------------------------------------------------
# Bagging regression nuisance model for outcome using ipred
#
# Returns:
#   mhat = E[Y | X, G]
# ------------------------------------------------------------
fit_predict_bag_reg <- function(Xtrain, ytrain, Xtest, args = NULL) {
  if (!requireNamespace("ipred", quietly = TRUE)) {
    warning("Package 'ipred' is not available. Fallback to lm.")
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)

  colnames(Xtest) <- colnames(Xtrain)

  ytrain <- as.numeric(ytrain)

  if (is.null(args)) {
    args <- list()
  }

  df_train <- data.frame(
    y = ytrain,
    as.data.frame(Xtrain, check.names = FALSE)
  )

  df_test <- as.data.frame(Xtest, check.names = FALSE)

  fit <- tryCatch(
    {
      do.call(
        ipred::bagging,
        c(
          list(
            formula = y ~ . - 1,
            data = df_train
          ),
          args
        )
      )
    },
    error = function(e) {
      warning(
        "bagging outcome model failed. Fallback to lm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred <- tryCatch(
    {
      stats::predict(
        fit,
        newdata = df_test
      )
    },
    error = function(e) {
      warning(
        "bagging outcome prediction failed. Fallback to lm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(pred)) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred <- as.numeric(pred)

  if (length(pred) != nrow(Xtest)) {
    warning("bagging outcome prediction length mismatch. Fallback to lm.")
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  if (any(!is.finite(pred))) {
    warning("bagging outcome prediction contains non-finite values. Fallback to lm.")
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred
}

# ------------------------------------------------------------
# Bagging classification nuisance model for propensity using ipred
#
# Returns:
#   ehat = P[A = 1 | X, G]
# ------------------------------------------------------------
fit_predict_bag_bin <- function(Xtrain, atrain, Xtest, args = NULL) {
  if (!requireNamespace("ipred", quietly = TRUE)) {
    warning("Package 'ipred' is not available. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)

  colnames(Xtest) <- colnames(Xtrain)

  atrain <- as.numeric(atrain)

  if (any(!atrain %in% c(0, 1))) {
    stop("A must be binary 0/1 for propensity model.", call. = FALSE)
  }

  if (length(unique(atrain)) < 2) {
    pred <- rep(mean(atrain, na.rm = TRUE), nrow(Xtest))
    return(pred)
  }

  if (is.null(args)) {
    args <- list()
  }

  df_train <- data.frame(
    a = factor(atrain, levels = c(0, 1)),
    as.data.frame(Xtrain, check.names = FALSE)
  )

  df_test <- as.data.frame(Xtest, check.names = FALSE)

  fit <- tryCatch(
    {
      do.call(
        ipred::bagging,
        c(
          list(
            formula = a ~ . - 1,
            data = df_train
          ),
          args
        )
      )
    },
    error = function(e) {
      warning(
        "bagging propensity model failed. Fallback to glm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pr <- tryCatch(
    {
      stats::predict(
        fit,
        newdata = df_test,
        type = "prob"
      )
    },
    error = function(e) {
      warning(
        "bagging propensity probability prediction failed. Fallback to glm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(pr)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  if (is.vector(pr)) {
    warning("bagging propensity prediction returned a vector, not probability matrix. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pr <- as.data.frame(pr)

  if ("1" %in% colnames(pr)) {
    pred <- as.numeric(pr[, "1"])
  } else if (ncol(pr) >= 2) {
    pred <- as.numeric(pr[, 2])
  } else {
    warning("bagging propensity probability matrix has unexpected columns. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  if (length(pred) != nrow(Xtest)) {
    warning("bagging propensity prediction length mismatch. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  if (any(!is.finite(pred))) {
    warning("bagging propensity prediction contains non-finite values. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred
}

# ------------------------------------------------------------
# Main prepare_dat function
# ------------------------------------------------------------
prepare_dat <- function(dat,
                        K,
                        nfold = 5,
                        outcome_method = "lm",
                        propensity_method = "glm",
                        outcome_args = NULL,
                        propensity_args = NULL,
                        seed = 1,
                        trim_e = 1e-4,
                        use_group_for_nuisance = TRUE,
                        add_group_interactions = TRUE,
                        calibrate_propensity = TRUE) {

  # ----------------------------------------------------------
  # Compatibility:
  # If old simu.R/main.R still passes xgboost/xgb,
  # do not call xgboost. Safely map to lm/glm.
  # ----------------------------------------------------------
  if (outcome_method %in% c("xgboost", "xgb")) {
    warning("xgboost has been removed from prepare.R. Use lm for outcome nuisance model.")
    outcome_method <- "lm"
    outcome_args <- NULL
  }

  if (propensity_method %in% c("xgboost", "xgb")) {
    warning("xgboost has been removed from prepare.R. Use glm for propensity nuisance model.")
    propensity_method <- "glm"
    propensity_args <- NULL
  }

  # ----------------------------------------------------------
  # Basic checks
  # ----------------------------------------------------------
  if (!all(c("Y", "A", "G") %in% colnames(dat))) {
    stop("dat must contain columns Y, A, and G.", call. = FALSE)
  }

  if (!is.numeric(K) || length(K) != 1 || K < 1) {
    stop("K must be a positive integer.", call. = FALSE)
  }

  K <- as.integer(K)

  if (!is.numeric(nfold) || length(nfold) != 1) {
    stop("nfold must be a numeric scalar.", call. = FALSE)
  }

  nfold <- as.integer(nfold)

  if (!is.numeric(trim_e) || length(trim_e) != 1 || !is.finite(trim_e)) {
    stop("trim_e must be a finite numeric scalar.", call. = FALSE)
  }

  if (trim_e <= 0 || trim_e >= 0.5) {
    stop("trim_e must be in (0, 0.5).", call. = FALSE)
  }

  set.seed(seed)

  Y <- as.numeric(dat$Y)
  A <- as.numeric(dat$A)
  G <- as.integer(dat$G)

  X <- get_X_from_dat(dat)
  X <- clean_matrix(X)

  n <- nrow(X)
  p <- ncol(X)

  if (length(Y) != n) {
    stop("Length of Y must equal nrow(X).", call. = FALSE)
  }

  if (length(A) != n) {
    stop("Length of A must equal nrow(X).", call. = FALSE)
  }

  if (length(G) != n) {
    stop("Length of G must equal nrow(X).", call. = FALSE)
  }

  if (any(!is.finite(Y))) {
    stop("Y contains non-finite values.", call. = FALSE)
  }

  if (any(!is.finite(A))) {
    stop("A contains non-finite values.", call. = FALSE)
  }

  if (any(!is.finite(G))) {
    stop("G contains non-finite values.", call. = FALSE)
  }

  if (any(G < 1 | G > K)) {
    stop("G must contain subgroup labels from 1 to K.", call. = FALSE)
  }

  if (any(!A %in% c(0, 1))) {
    stop("A must be binary 0/1.", call. = FALSE)
  }

  if (nfold < 2) {
    stop("nfold must be at least 2.", call. = FALSE)
  }

  if (nfold > n) {
    stop("nfold cannot exceed sample size n.", call. = FALSE)
  }

  # ----------------------------------------------------------
  # Nuisance feature matrix
  # ----------------------------------------------------------
  if (isTRUE(use_group_for_nuisance)) {
    Xnuis <- build_nuisance_X(
      X = X,
      G = G,
      K = K,
      add_group_interactions = add_group_interactions
    )
  } else {
    Xnuis <- X
  }

  Xnuis <- clean_matrix(Xnuis)

  # ----------------------------------------------------------
  # Cross-fitting fold assignment
  # ----------------------------------------------------------
  foldid <- sample(rep(seq_len(nfold), length.out = n))

  mhat <- rep(NA_real_, n)
  ehat <- rep(NA_real_, n)

  # ----------------------------------------------------------
  # Cross-fitting loop
  # ----------------------------------------------------------
  for (ff in seq_len(nfold)) {
    idx_test <- which(foldid == ff)
    idx_train <- setdiff(seq_len(n), idx_test)

    Xtr <- Xnuis[idx_train, , drop = FALSE]
    Xte <- Xnuis[idx_test, , drop = FALSE]

    Ytr <- Y[idx_train]
    Atr <- A[idx_train]
    Gtr <- G[idx_train]
    Gte <- G[idx_test]

    # --------------------------------------------------------
    # Outcome model m(X,G) = E[Y | X, G]
    # --------------------------------------------------------
    if (outcome_method %in% c("lm", "linear", "glm")) {

      mhat[idx_test] <- fit_predict_lm_reg(
        Xtrain = Xtr,
        ytrain = Ytr,
        Xtest = Xte,
        args = outcome_args
      )

    } else if (outcome_method %in% c("rf", "randomForest", "randomforest")) {

      mhat[idx_test] <- fit_predict_rf_reg(
        Xtrain = Xtr,
        ytrain = Ytr,
        Xtest = Xte,
        args = outcome_args
      )

    } else if (outcome_method %in% c("bagging", "bag")) {

      mhat[idx_test] <- fit_predict_bag_reg(
        Xtrain = Xtr,
        ytrain = Ytr,
        Xtest = Xte,
        args = outcome_args
      )

    } else {

      stop(
        "Unknown outcome_method: ", outcome_method, ". ",
        "Supported outcome methods are: lm, glm, linear, rf, randomForest, bagging, bag.",
        call. = FALSE
      )

    }

    # --------------------------------------------------------
    # Propensity model e(X,G) = P[A = 1 | X, G]
    # --------------------------------------------------------
    if (propensity_method %in% c("glm", "logit", "logistic")) {

      ehat_ff <- fit_predict_glm_bin(
        Xtrain = Xtr,
        atrain = Atr,
        Xtest = Xte,
        args = propensity_args
      )

    } else if (propensity_method %in% c("rf", "randomForest", "randomforest")) {

      ehat_ff <- fit_predict_rf_bin(
        Xtrain = Xtr,
        atrain = Atr,
        Xtest = Xte,
        args = propensity_args
      )

    } else if (propensity_method %in% c("bagging", "bag")) {

      ehat_ff <- fit_predict_bag_bin(
        Xtrain = Xtr,
        atrain = Atr,
        Xtest = Xte,
        args = propensity_args
      )

    } else {

      stop(
        "Unknown propensity_method: ", propensity_method, ". ",
        "Supported propensity methods are: glm, logit, logistic, rf, randomForest, bagging, bag.",
        call. = FALSE
      )

    }

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

  # ----------------------------------------------------------
  # Safety checks
  # ----------------------------------------------------------
  if (any(!is.finite(mhat))) {
    stop("Non-finite mhat detected in prepare_dat().", call. = FALSE)
  }

  if (any(!is.finite(ehat))) {
    stop("Non-finite ehat detected in prepare_dat().", call. = FALSE)
  }

  ehat <- pmin(pmax(ehat, trim_e), 1 - trim_e)

  Ytilde <- Y - mhat
  Atilde <- A - ehat

  if (any(!is.finite(Ytilde))) {
    stop("Non-finite Ytilde detected in prepare_dat().", call. = FALSE)
  }

  if (any(!is.finite(Atilde))) {
    stop("Non-finite Atilde detected in prepare_dat().", call. = FALSE)
  }

  # ----------------------------------------------------------
  # Return mydat object
  # ----------------------------------------------------------
  list(
    Y = Y,
    A = A,
    G = G,
    X = X,
    K = K,
    p = p,
    n = n,

    mhat = mhat,
    ehat = ehat,
    Ytilde = Ytilde,
    Atilde = Atilde,

    foldid = foldid,
    nfold = nfold,

    outcome_method = outcome_method,
    propensity_method = propensity_method,
    outcome_args = outcome_args,
    propensity_args = propensity_args,

    use_group_for_nuisance = use_group_for_nuisance,
    add_group_interactions = add_group_interactions,
    calibrate_propensity = calibrate_propensity,
    trim_e = trim_e
  )
}