#' Generate Data
#' Creates covariate matrix and response for linear or logistic models.
#' @param type A character (linear, logistic) that specifies the type of model
#' @param n integer for number of observations, i.e., rows
#' @param p integer for number of predictors, i.e., columns, (not including intercept)
#' @param rho correlation between covariates if using multivariate normal matrix
#' @param mu means of predictors if using multivariate normal matrix
#' @param stdev standard deviation of predictors if using multivariate normal matrix
#' @param spars sparsity of covariate matrix by percentage
#' @param betas vector of true coefficients
#' @param train_size proportion of data to be used in training
#' @param sigma standard deviation of residuals for linear models
#' @param skew distribution of covariates. If null, multivariate normal is used.
#' @return Generated data
#' @export
gen_data <- function(type, n, p, rho, mu,
                    stdev, spars,
                    betas, train_size, sigma = NULL,
                    skew = NULL){
  if(is.null(skew)){
    X <- faux::rnorm_multi(n, vars = p, mu = mu, sd = stdev, r = rho)
    X <- as.matrix(X)
  }

  if(!is.null(skew)){
    if(skew == "exponential"){
      distribution <- function(x) pexp(x)
    }

    if(skew == "gamma"){
      distribution <- function(x) pgamma(x, shape = 2, scale = 2)
    }

    if(skew == "beta"){
      distribution <- function(x) pbeta(x, shape1 = 2, shape2 = 5)
    }

    corr_mat <- randcorr::randcorr(p)
    X <- generate_matrix_of_observations(corr_mat = corr_mat,
                                         cdf = distribution,
                                         num_rows = n)
  }

  BETAS_i <- betas
  BETAS_i[sample(1:p, round(p*spars))] <- 0

  if(type == "logistic"){
    etas <- c(data.matrix(X) %*% BETAS_i)
    P <- plogis(etas)
    Y <- rbinom(n, 1, P)
  }

  # Need error handling to make sure sigma is specified for linear models
  if(type == "linear"){
    Y <- X %*% BETAS_i + rnorm(n, mean = 0, sd = sigma)
  }

  n_train <- as.integer(length(Y)*train_size)
  train_idx <- sample.int(length(Y), size = n_train, replace = FALSE)
  X_train <- X[train_idx,]
  X_test <- X[-train_idx,]
  Y_train <- Y[train_idx]
  Y_test <- Y[-train_idx]

  list(
    xtrain = X_train,
    xtest = X_test,
    ytrain = Y_train,
    ytest = Y_test
  )
}

#' Fit and test regression model
#' @param type linear or logistic
#' @param xtrain training set covariate matrix
#' @param xtest test set covariate matrix
#' @param ytrain training set response
#' @param ytest test set response
#' @param nfolds number of folds to use in cross-validation
#' @return coefficients of fitted models and test set errors (RMSE or accuracy)
#' @export
fit_regressions <- function(type, xtrain, xtest, ytrain,
                           ytest, nfolds){
  if(type == "logistic"){
    family <- "binomial"
  }

  if(type == "linear"){
    family <- "gaussian"
  }

  lasso_fit <- glmnet::cv.glmnet(xtrain, ytrain, family = family,
                        alpha = 1, nfolds = nfolds)
  ridge_fit <- glmnet::cv.glmnet(xtrain, ytrain, family = family,
                        alpha = 0, nfolds = nfolds)

  lasso.1se <- lasso_fit$lambda.1se
  lasso.min <- lasso_fit$lambda.min
  ridge.1se <- ridge_fit$lambda.1se
  ridge.min <- ridge_fit$lambda.min

  if(type == "logistic"){
    resp_type <- "class"
  }

  if(type == "linear"){
    resp_type <- "response"
  }

  y_pred_lasso_1se <- as.numeric(predict(lasso_fit, s = lasso.1se,
                                        newx = xtest, type = resp_type))
  y_pred_lasso_min <- as.numeric(predict(lasso_fit, s = lasso.min,
                                        newx = xtest, type = resp_type))
  y_pred_ridge_1se <- as.numeric(predict(ridge_fit, s = ridge.1se,
                                    newx = xtest, type = resp_type))
  y_pred_ridge_min <- as.numeric(predict(ridge_fit, s = ridge.min,
                                        newx = xtest, type = resp_type))

  if(type == "logistic"){
    acc_lam_1se <- mean(y_pred_lasso_1se == ytest)
    acc_lam_min <- mean(y_pred_lasso_min == ytest)
    acc_ridge_1se <- mean(y_pred_ridge_1se == ytest)
    acc_ridge_min <- mean(y_pred_ridge_min == ytest)
  }

  if(type == "linear"){
    acc_lam_1se <- sqrt(mean((y_pred_lasso_1se - ytest)^2))
    acc_lam_min <- sqrt(mean((y_pred_lasso_min - ytest)^2))
    acc_ridge_1se <- sqrt(mean((y_pred_ridge_1se - ytest)^2))
    acc_ridge_min <- sqrt(mean((y_pred_ridge_min - ytest)^2))

  }

  list(
    lasso_min_coefs = as.vector(coef(lasso_fit, s = lasso.min)),
    lasso_1se_coefs = as.vector(coef(lasso_fit, s = lasso.1se)),
    ridge_min_coefs = as.vector(coef(ridge_fit, s = ridge.min)),
    ridge_1se_coefs = as.vector(coef(ridge_fit, s = ridge.1se)),
    acc.lasso.1se = acc_lam_1se,
    acc.lasso.min = acc_lam_min,
    acc.ridge.1se = acc_ridge_1se,
    acc.ridge.min = acc_ridge_min
  )
}

#' Runs simulation study to assess models
#' @param type linear or logistic model
#' @param ntrials number of trials in simulation
#' @param n number of observations
#' @param p number of predictors
#' @param rho correlation of predictors if using multivariate normal
#' @param mu means of predictors if using multivariate normal
#' @param stdev standard deviation of predictors if using multivariate normal
#' @param spars sparsity of covariate matrix as a proportion
#' @param betas true values of the coefficients
#' @param train_size proportion of data to use in training
#' @param nfolds number of folds for cross-validation
#' @param sigma standard deviation of the residuals for linear models
#' @param skew distribution of covariate matrix
#' @return matrix of fitted coefficients, mean test set accuracies
#' @export
run_simulations <- function(type, ntrials = 20, n = 1000, p = 10, rho = 0.5,
                           mu = runif(p, 1, 5),
                           stdev = runif(p, 0.2,20), spars = 0.5,
                           betas = runif(p, 0.01,0.2),
                           train_size = 0.75, nfolds = 5, sigma = NULL,
                           skew = NULL){
  lasso_1se_coefs <- matrix(0, nrow = ntrials, ncol = p+1)
  lasso_min_coefs <- matrix(0, nrow = ntrials, ncol = p+1)
  ridge_1se_coefs <- matrix(0, nrow = ntrials, ncol = p+1)
  ridge_min_coefs <- matrix(0, nrow = ntrials, ncol = p+1)

  lasso_1se_acc <- rep(NA, ntrials)
  lasso_min_acc <- rep(NA, ntrials)
  ridge_min_acc <- rep(NA, ntrials)
  ridge_1se_acc <- rep(NA, ntrials)

  for(i in 1:ntrials){
    data <- gen_data(type = type, n = n, p = p, rho = rho, mu = mu,
                    stdev = stdev, spars = spars,
                    betas = betas, train_size = train_size, sigma = sigma,
                    skew = skew)

    results <- fit_regressions(type = type, data$xtrain, data$xtest,
                              data$ytrain, data$ytest,
                              nfolds = nfolds)

    lasso_1se_coefs[i,] <- results$lasso_1se_coefs
    lasso_min_coefs[i,] <- results$lasso_min_coefs
    ridge_1se_coefs[i,] <- results$ridge_1se_coefs
    ridge_min_coefs[i,] <- results$ridge_min_coefs

    lasso_1se_acc[i] <- results$acc.lasso.1se
    lasso_min_acc[i] <- results$acc.lasso.min
    ridge_1se_acc[i] <- results$acc.ridge.1se
    ridge_min_acc[i] <-results$acc.ridge.min
  }

  list(
    lasso.1se.coefs = lasso_1se_coefs,
    lasso.min.coefs = lasso_min_coefs,
    ridge.1se.coefs = ridge_1se_coefs,
    ridge.min.coefs = ridge_min_coefs,
    mean.lasso.1se = mean(lasso_1se_acc),
    mean.lasso.min = mean(lasso_min_acc),
    mean.ridge.1se = mean(ridge_1se_acc),
    mean.ridge.min = mean(ridge_min_acc)
  )
}
