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
# Nuisance models:
#   Outcome:    lm/glm/linear, rf/randomForest, bagging/bag
#   Propensity: glm/logit/logistic, rf/randomForest, bagging/bag
#
# Optional DC screening:
#   Before cross-fitting, sample at most dc_screen_n rows.
#   Compute distance-correlation scores for each original X column:
#     score_y = dc(Xj, Y)
#     score_a = dc(Xj, A)
#     score   = pmax(score_y, score_a)
#   Sort scores decreasingly.
#   Use the last dc_tail_n scores to define cutoff:
#     cutoff = mean(tail_scores) + dc_z * sd(tail_scores)
#   Select variables with score > cutoff.
#   Cross-fitting then uses the same selected variables in all folds.
# ============================================================

get_X_from_dat <- function(dat) {
  if (!all(c("Y", "A", "G") %in% colnames(dat))) {
    stop("dat must contain columns Y, A, and G.", call. = FALSE)
  }

  xcols <- setdiff(colnames(dat), c("Y", "A", "G"))
  if (length(xcols) == 0) {
    stop("No covariate columns found in dat.", call. = FALSE)
  }

  X <- as.data.frame(dat[, xcols, drop = FALSE])
  for (jj in seq_along(X)) {
    X[[jj]] <- as.numeric(X[[jj]])
  }

  X <- as.matrix(X)
  storage.mode(X) <- "double"
  colnames(X) <- make.names(xcols, unique = TRUE)
  X
}

clean_matrix <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }

  colnames(X) <- make.names(colnames(X), unique = TRUE)
  X
}

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
# Distance correlation screening
# ------------------------------------------------------------

.double_center_dist <- function(z) {
  z <- as.numeric(z)
  D <- abs(outer(z, z, "-"))
  D - rowMeans(D) - matrix(colMeans(D), nrow(D), ncol(D), byrow = TRUE) + mean(D)
}

.fast_dcor_1d <- function(x, y_centered_dist) {
  x <- as.numeric(x)

  if (length(unique(x)) <= 1) {
    return(0)
  }

  A <- .double_center_dist(x)

  dcov2 <- mean(A * y_centered_dist)
  dvarx <- mean(A * A)
  dvary <- mean(y_centered_dist * y_centered_dist)

  if (!is.finite(dcov2) || !is.finite(dvarx) || !is.finite(dvary)) {
    return(0)
  }

  if (dvarx <= 0 || dvary <= 0) {
    return(0)
  }

  out <- sqrt(max(dcov2, 0)) / sqrt(sqrt(dvarx * dvary))

  if (!is.finite(out)) {
    out <- 0
  }

  out
}
dc_screen_variables <- function(X,
                                target,
                                screen_n = 500,
                                tail_n = 30,
                                z = 3,
                                seed = 1,
                                min_selected = 1,
                                verbose = FALSE,
                                target_name = "target") {
  X <- clean_matrix(X)
  target <- as.numeric(target)

  n <- nrow(X)
  p <- ncol(X)

  if (length(target) != n) {
    stop("Length of target must equal nrow(X).", call. = FALSE)
  }

  if (!is.numeric(screen_n) || length(screen_n) != 1 || !is.finite(screen_n)) {
    stop("screen_n must be a finite numeric scalar.", call. = FALSE)
  }

  if (!is.numeric(tail_n) || length(tail_n) != 1 || !is.finite(tail_n)) {
    stop("tail_n must be a finite numeric scalar.", call. = FALSE)
  }

  if (!is.numeric(z) || length(z) != 1 || !is.finite(z)) {
    stop("z must be a finite numeric scalar.", call. = FALSE)
  }

  if (!is.numeric(min_selected) || length(min_selected) != 1 || !is.finite(min_selected)) {
    stop("min_selected must be a finite numeric scalar.", call. = FALSE)
  }

  screen_n <- as.integer(screen_n)
  tail_n <- as.integer(tail_n)
  min_selected <- as.integer(min_selected)

  if (screen_n < 2) {
    stop("screen_n must be at least 2.", call. = FALSE)
  }

  if (tail_n < 1) {
    stop("tail_n must be at least 1.", call. = FALSE)
  }

  if (min_selected < 1) {
    stop("min_selected must be at least 1.", call. = FALSE)
  }

  set.seed(seed)
  idx_screen <- sample(seq_len(n), size = min(screen_n, n), replace = FALSE)

  Xs <- X[idx_screen, , drop = FALSE]
  ts <- target[idx_screen]

  Bt <- .double_center_dist(ts)

  score <- numeric(p)

  for (j in seq_len(p)) {
    score[j] <- .fast_dcor_1d(Xs[, j], Bt)
  }

  score[!is.finite(score)] <- 0

  ord <- order(score, decreasing = TRUE)

  tail_m <- min(tail_n, p)
  tail_scores <- score[ord][seq.int(from = p - tail_m + 1, to = p)]

  tail_mean <- mean(tail_scores)
  tail_sd <- stats::sd(tail_scores)

  if (!is.finite(tail_sd)) {
    tail_sd <- 0
  }

  cutoff <- tail_mean + z * tail_sd

  selected_idx <- which(score > cutoff)

  if (length(selected_idx) < min_selected) {
    selected_idx <- ord[seq_len(min(min_selected, p))]
  } else {
    selected_idx <- selected_idx[order(score[selected_idx], decreasing = TRUE)]
  }

  selected_vars <- colnames(X)[selected_idx]

  if (isTRUE(verbose)) {
    message(
      "DC screening for ",
      target_name,
      " selected ",
      length(selected_idx),
      " / ",
      p,
      " variables. cutoff = ",
      signif(cutoff, 4)
    )
  }

  list(
    selected_idx = selected_idx,
    selected_vars = selected_vars,
    score = score,
    order = ord,
    cutoff = cutoff,
    tail_mean = tail_mean,
    tail_sd = tail_sd,
    tail_n = tail_m,
    screen_idx = idx_screen,
    screen_n = length(idx_screen),
    z = z,
    target_name = target_name
  )
}

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

fit_predict_lm_reg <- function(Xtrain, ytrain, Xtest, args = NULL) {
  Xtrain <- clean_matrix(Xtrain)
  Xtest <- clean_matrix(Xtest)
  colnames(Xtest) <- colnames(Xtrain)

  ytrain <- as.numeric(ytrain)

  if (length(ytrain) != nrow(Xtrain)) {
    stop("Length of ytrain must equal nrow(Xtrain).", call. = FALSE)
  }

  df_train <- data.frame(y = ytrain, as.data.frame(Xtrain, check.names = FALSE))
  df_test <- as.data.frame(Xtest, check.names = FALSE)

  fit <- tryCatch(
    stats::lm(y ~ . - 1, data = df_train),
    error = function(e) {
      warning("lm outcome model failed. Use mean outcome. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(rep(mean(ytrain, na.rm = TRUE), nrow(Xtest)))
  }

  pred <- tryCatch(
    stats::predict(fit, newdata = df_test),
    error = function(e) {
      warning("lm outcome prediction failed. Use mean outcome. Original error: ", conditionMessage(e))
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
    return(rep(mean(atrain, na.rm = TRUE), nrow(Xtest)))
  }

  df_train <- data.frame(a = atrain, as.data.frame(Xtrain, check.names = FALSE))
  df_test <- as.data.frame(Xtest, check.names = FALSE)

  fit <- tryCatch(
    suppressWarnings(stats::glm(a ~ . - 1, data = df_train, family = stats::binomial())),
    error = function(e) {
      warning("glm propensity model failed. Use mean treatment probability. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    pred <- rep(mean(atrain, na.rm = TRUE), nrow(Xtest))
  } else {
    pred <- tryCatch(
      stats::predict(fit, newdata = df_test, type = "response"),
      error = function(e) {
        warning("glm propensity prediction failed. Use mean treatment probability. Original error: ", conditionMessage(e))
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
    do.call(
      randomForest::randomForest,
      c(list(x = as.data.frame(Xtrain, check.names = FALSE), y = ytrain), args)
    ),
    error = function(e) {
      warning("randomForest outcome model failed. Fallback to lm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred <- tryCatch(
    stats::predict(fit, newdata = as.data.frame(Xtest, check.names = FALSE)),
    error = function(e) {
      warning("randomForest outcome prediction failed. Fallback to lm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(pred)) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred <- as.numeric(pred)

  if (length(pred) != nrow(Xtest) || any(!is.finite(pred))) {
    warning("randomForest outcome prediction invalid. Fallback to lm.")
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred
}

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
    return(rep(mean(atrain, na.rm = TRUE), nrow(Xtest)))
  }

  if (is.null(args)) {
    args <- list()
  }

  fit <- tryCatch(
    do.call(
      randomForest::randomForest,
      c(
        list(
          x = as.data.frame(Xtrain, check.names = FALSE),
          y = factor(atrain, levels = c(0, 1))
        ),
        args
      )
    ),
    error = function(e) {
      warning("randomForest propensity model failed. Fallback to glm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pr <- tryCatch(
    stats::predict(fit, newdata = as.data.frame(Xtest, check.names = FALSE), type = "prob"),
    error = function(e) {
      warning("randomForest propensity probability prediction failed. Fallback to glm. Original error: ", conditionMessage(e))
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

  if (length(pred) != nrow(Xtest) || any(!is.finite(pred))) {
    warning("randomForest propensity prediction invalid. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred
}

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

  df_train <- data.frame(y = ytrain, as.data.frame(Xtrain, check.names = FALSE))
  df_test <- as.data.frame(Xtest, check.names = FALSE)

  fit <- tryCatch(
    do.call(ipred::bagging, c(list(formula = y ~ . - 1, data = df_train), args)),
    error = function(e) {
      warning("bagging outcome model failed. Fallback to lm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred <- tryCatch(
    stats::predict(fit, newdata = df_test),
    error = function(e) {
      warning("bagging outcome prediction failed. Fallback to lm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(pred)) {
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred <- as.numeric(pred)

  if (length(pred) != nrow(Xtest) || any(!is.finite(pred))) {
    warning("bagging outcome prediction invalid. Fallback to lm.")
    return(fit_predict_lm_reg(Xtrain, ytrain, Xtest, args = NULL))
  }

  pred
}

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
    return(rep(mean(atrain, na.rm = TRUE), nrow(Xtest)))
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
    do.call(ipred::bagging, c(list(formula = a ~ . - 1, data = df_train), args)),
    error = function(e) {
      warning("bagging propensity model failed. Fallback to glm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pr <- tryCatch(
    stats::predict(fit, newdata = df_test, type = "prob"),
    error = function(e) {
      warning("bagging propensity probability prediction failed. Fallback to glm. Original error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(pr) || is.vector(pr)) {
    warning("bagging propensity prediction invalid. Fallback to glm.")
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

  if (length(pred) != nrow(Xtest) || any(!is.finite(pred))) {
    warning("bagging propensity prediction invalid. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred
}

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
                        calibrate_propensity = TRUE,
                        use_dc_screening = TRUE,
                        dc_screen_n = 500,
                        dc_tail_n = 30,
                        dc_z = 3,
                        dc_min_selected = 1,
                        dc_seed = seed,
                        dc_verbose = FALSE) {
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

  dc_y_info <- NULL
dc_a_info <- NULL

X_model_y <- X
X_model_a <- X

if (isTRUE(use_dc_screening)) {
  dc_y_info <- dc_screen_variables(
    X = X,
    target = Y,
    screen_n = dc_screen_n,
    tail_n = dc_tail_n,
    z = dc_z,
    seed = dc_seed,
    min_selected = dc_min_selected,
    verbose = dc_verbose,
    target_name = "Y"
  )

  dc_a_info <- dc_screen_variables(
    X = X,
    target = A,
    screen_n = dc_screen_n,
    tail_n = dc_tail_n,
    z = dc_z,
    seed = dc_seed,
    min_selected = dc_min_selected,
    verbose = dc_verbose,
    target_name = "A"
  )

  X_model_y <- X[, dc_y_info$selected_idx, drop = FALSE]
  X_model_a <- X[, dc_a_info$selected_idx, drop = FALSE]

  X_model_y <- clean_matrix(X_model_y)
  X_model_a <- clean_matrix(X_model_a)
}

p_model_y <- ncol(X_model_y)
p_model_a <- ncol(X_model_a)

if (isTRUE(use_group_for_nuisance)) {
  Xnuis_y <- build_nuisance_X(
    X = X_model_y,
    G = G,
    K = K,
    add_group_interactions = add_group_interactions
  )

  Xnuis_a <- build_nuisance_X(
    X = X_model_a,
    G = G,
    K = K,
    add_group_interactions = add_group_interactions
  )
} else {
  Xnuis_y <- X_model_y
  Xnuis_a <- X_model_a
}

Xnuis_y <- clean_matrix(Xnuis_y)
Xnuis_a <- clean_matrix(Xnuis_a)

  foldid <- sample(rep(seq_len(nfold), length.out = n))

  mhat <- rep(NA_real_, n)
  ehat <- rep(NA_real_, n)

  for (ff in seq_len(nfold)) {
    idx_test <- which(foldid == ff)
    idx_train <- setdiff(seq_len(n), idx_test)

    Xtr_y <- Xnuis_y[idx_train, , drop = FALSE]
Xte_y <- Xnuis_y[idx_test, , drop = FALSE]

Xtr_a <- Xnuis_a[idx_train, , drop = FALSE]
Xte_a <- Xnuis_a[idx_test, , drop = FALSE]

    Ytr <- Y[idx_train]
    Atr <- A[idx_train]
    Gtr <- G[idx_train]
    Gte <- G[idx_test]

    if (outcome_method %in% c("lm", "linear", "glm")) {
      mhat[idx_test] <- fit_predict_lm_reg(Xtr_y, Ytr, Xte_y, outcome_args)
    } else if (outcome_method %in% c("rf", "randomForest", "randomforest")) {
      mhat[idx_test] <- fit_predict_rf_reg(Xtr_y, Ytr, Xte_y, outcome_args)
    } else if (outcome_method %in% c("bagging", "bag")) {
      mhat[idx_test] <- fit_predict_bag_reg(Xtr_y, Ytr, Xte_y, outcome_args)
    } else {
      stop(
        "Unknown outcome_method: ", outcome_method, ". ",
        "Supported outcome methods are: lm, glm, linear, rf, randomForest, bagging, bag.",
        call. = FALSE
      )
    }

    if (propensity_method %in% c("glm", "logit", "logistic")) {
      ehat_ff <- fit_predict_glm_bin(Xtr_a, Atr, Xte_a, propensity_args)
    } else if (propensity_method %in% c("rf", "randomForest", "randomforest")) {
      ehat_ff <- fit_predict_rf_bin(Xtr_a, Atr, Xte_a, propensity_args)
    } else if (propensity_method %in% c("bagging", "bag")) {
      ehat_ff <- fit_predict_bag_bin(Xtr_a, Atr, Xte_a, propensity_args)
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

  out <- list(
    Y = Y,
    A = A,
    G = G,
    X = X,
    X_model_y = X_model_y,
X_model_a = X_model_a,
    K = K,
    p = p,
    p_model_y = p_model_y,
p_model_a = p_model_a,
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
    trim_e = trim_e,

    use_dc_screening = use_dc_screening
  )

 if (isTRUE(use_dc_screening)) {
  out$dc_y_score <- dc_y_info$score
  out$dc_y_order <- dc_y_info$order
  out$dc_y_selected_idx <- dc_y_info$selected_idx
  out$dc_y_selected_vars <- dc_y_info$selected_vars
  out$dc_y_cutoff <- dc_y_info$cutoff
  out$dc_y_tail_mean <- dc_y_info$tail_mean
  out$dc_y_tail_sd <- dc_y_info$tail_sd
  out$dc_y_tail_n <- dc_y_info$tail_n
  out$dc_y_screen_idx <- dc_y_info$screen_idx
  out$dc_y_screen_n <- dc_y_info$screen_n

  out$dc_a_score <- dc_a_info$score
  out$dc_a_order <- dc_a_info$order
  out$dc_a_selected_idx <- dc_a_info$selected_idx
  out$dc_a_selected_vars <- dc_a_info$selected_vars
  out$dc_a_cutoff <- dc_a_info$cutoff
  out$dc_a_tail_mean <- dc_a_info$tail_mean
  out$dc_a_tail_sd <- dc_a_info$tail_sd
  out$dc_a_tail_n <- dc_a_info$tail_n
  out$dc_a_screen_idx <- dc_a_info$screen_idx
  out$dc_a_screen_n <- dc_a_info$screen_n

  out$dc_z <- dc_z
} else {
  out$dc_y_score <- NULL
  out$dc_y_order <- NULL
  out$dc_y_selected_idx <- seq_len(p)
  out$dc_y_selected_vars <- colnames(X)
  out$dc_y_cutoff <- NULL
  out$dc_y_tail_mean <- NULL
  out$dc_y_tail_sd <- NULL
  out$dc_y_tail_n <- NULL
  out$dc_y_screen_idx <- NULL
  out$dc_y_screen_n <- NULL

  out$dc_a_score <- NULL
  out$dc_a_order <- NULL
  out$dc_a_selected_idx <- seq_len(p)
  out$dc_a_selected_vars <- colnames(X)
  out$dc_a_cutoff <- NULL
  out$dc_a_tail_mean <- NULL
  out$dc_a_tail_sd <- NULL
  out$dc_a_tail_n <- NULL
  out$dc_a_screen_idx <- NULL
  out$dc_a_screen_n <- NULL

  out$dc_z <- NULL
}

  out
}