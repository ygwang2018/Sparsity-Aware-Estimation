# ------------------------------------------------------------
# Bagging regression nuisance model for propensity using ipred
#
# This estimates:
#   ehat = E[A | X, G]
#
# Since A is binary 0/1:
#   E[A | X, G] = P[A = 1 | X, G]
#
# Important:
#   This function treats A as numeric, not as a classification label.
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

  if (length(atrain) != nrow(Xtrain)) {
    stop("Length of atrain must equal nrow(Xtrain).", call. = FALSE)
  }

  if (any(!atrain %in% c(0, 1))) {
    stop("A must be binary 0/1 for propensity model.", call. = FALSE)
  }

  # If training fold has only one treatment class, the conditional mean
  # is estimated as the observed treatment mean in that fold.
  if (length(unique(atrain)) < 2) {
    pred <- rep(mean(atrain, na.rm = TRUE), nrow(Xtest))
    return(pred)
  }

  if (is.null(args)) {
    args <- list()
  }

  df_train <- data.frame(
    a = atrain,
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
        "bagging propensity regression model failed. Fallback to glm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
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
        "bagging propensity regression prediction failed. Fallback to glm. Original error: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(pred)) {
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred <- as.numeric(pred)

  if (length(pred) != nrow(Xtest)) {
    warning("bagging propensity regression prediction length mismatch. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  if (any(!is.finite(pred))) {
    warning("bagging propensity regression prediction contains non-finite values. Fallback to glm.")
    return(fit_predict_glm_bin(Xtrain, atrain, Xtest, args = NULL))
  }

  pred
}