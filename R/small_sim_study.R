#' Simulate regressions under different conditions
#' @param type linear or logistic
#' @param n_folds vector of number of folds to vary
#' @param n_obs vector of number of observations to vary
#' @param p number of predictors
#' @param rho correlation of predictors for MVN
#' @param mu means of predictors for MVN
#' @param stddev standard deviation of predictors for MVN
#' @param spars sparsity of predictors as a proportion
#' @param betas true coefficients
#' @param sigma standard deviation of residuals for linear models
#' @param skew distribution of covariate matrix
#' @return tidy dataframe of mean number of non-zero coefficients and associated errors
#' @export
#' @examples
#' study_data_no_skew <- simulate_regressions(type = "linear", sigma=1)
#' study_data_exp <- simulate_regressions(type = "linear", sigma = 1, skew = "exponential")
#' study_data_beta <- simulate_regressions(type = "linear", sigma = 1, skew = "beta")
#' combined <- dplyr::bind_rows(study_data_beta, study_data_exp, study_data_no_skew)
#' combined %>%
#' dplyr::mutate(skew = dplyr::case_when(is.na(skew) ~ "normal",
#' skew == "exponential" ~ "exponential", skew == "beta" ~ "beta")) %>%
#' dplyr::group_by(model, skew) %>% dplyr::summarise(mean_rmse = mean(avgError),
#' mean_params = mean(params))
simulate_regressions <- function(type,
                                 n_folds = c(3, 5, 10),
                                 n_obs = c(100, 500, 1000),
                                 p = 10,
                                 rho = 0.5,
                                 mu = c(1, 5),
                                 stddev = c(0.2, 20),
                                 spars = 0.5,
                                 betas = c(0.01, 0.2),
                                 sigma = NULL,
                                 skew = NULL){

  mu_ <- runif(p, mu[1], mu[2])
  sigma_ <- runif(p, stddev[1], stddev[2])
  beta_ <- runif(p, betas[1], betas[2])

  obs_per_fold <- length(n_obs)
  final_df <- tibble::tibble()

  for(i in 1:length(n_folds)){

    df <- as.data.frame(matrix(NA, nrow = 4, ncol = 4),
                        row.names = c('lasso.1se', 'lasso.min',
                                      'ridge.1se', 'ridge.min'))

    colnames(df) <- c('avgError', 'params', 'nfolds', 'nobs')

    lasso_1se_mean_coeffs <- rep(NA, obs_per_fold)
    lasso_min_mean_coeffs <- rep(NA, obs_per_fold)
    lasso_1se_means <- rep(NA, obs_per_fold)
    lasso_min_means <- rep(NA, obs_per_fold)

    ridge_1se_mean_coeffs <- rep(NA, obs_per_fold)
    ridge_min_mean_coeffs <- rep(NA, obs_per_fold)
    ridge_1se_means <- rep(NA, obs_per_fold)
    ridge_min_means <- rep(NA, obs_per_fold)

    for(j in 1:length(n_obs)){

      result <- run_simulations(type = type, n = n_obs[j],
                                nfolds = n_folds[i],
                                p = p,
                                rho = rho,
                                mu = mu_,
                                stdev = sigma_,
                                spars = spars,
                                betas = beta_,
                                sigma = sigma,
                                skew = skew)

      lasso_1se_mean_coeffs[j] <- mean(colSums(t(result$lasso.1se.coefs) != 0))
      lasso_min_mean_coeffs[j] <- mean(colSums(t(result$lasso.min.coefs) != 0))
      lasso_1se_means[j] <- result$mean.lasso.1se
      lasso_min_means[j] <- result$mean.lasso.min

      ridge_1se_mean_coeffs[j] <- mean(colSums(t(result$ridge.1se.coefs) != 0))
      ridge_min_mean_coeffs[j] <- mean(colSums(t(result$ridge.min.coefs) != 0))
      ridge_1se_means[j] <- result$mean.ridge.1se
      ridge_min_means[j] <- result$mean.ridge.min

      df[, 'nfolds'] <- n_folds[i]
      df[, 'nobs'] <- n_obs[j]

      df['lasso.1se', 'avgError'] <- mean(lasso_1se_means, na.rm=TRUE)
      df['lasso.min', 'avgError'] <- mean(lasso_min_means, na.rm=TRUE)
      df['lasso.1se', 'params'] <- mean(lasso_1se_mean_coeffs, na.rm=TRUE)
      df['lasso.min', 'params'] <- mean(lasso_min_mean_coeffs, na.rm=TRUE)

      df['ridge.1se', 'avgError'] <- mean(ridge_1se_means, na.rm=TRUE)
      df['ridge.min', 'avgError'] <- mean(ridge_min_means, na.rm=TRUE)
      df['ridge.1se', 'params'] <- mean(ridge_1se_mean_coeffs, na.rm=TRUE)
      df['ridge.min', 'params'] <- mean(ridge_min_mean_coeffs, na.rm=TRUE)

      tidy_df <- tibble::rownames_to_column(df, var = "model")
      final_df <- dplyr::bind_rows(final_df, tidy_df)
    }
  }
  final_df <- final_df %>% dplyr::mutate(skew = skew)
}
