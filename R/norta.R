#' Determine x for NORTA
#' @param z value to evaluate cdf
#' @param x_cdf arbitrary cdf
#' @param tol numerical tolerance
#' @return lower
#' @export
determine_x <- function(z, x_cdf, tol = 10^-10){
  phi_z <- pnorm(z)

  lower <- -1

  while(x_cdf(lower) >= phi_z){
    lower <- lower * -2
  }

  upper <- 1

  while(x_cdf(upper) <= phi_z){
    upper <- upper * 2
  }

  while(upper - lower >= tol){
    middle <- lower + (upper - lower)/2.0
    x_cdf_middle <- x_cdf(middle)

    if(x_cdf_middle < phi_z){
      lower <- middle
    }

    else{
      upper <- middle
    }
  }
  lower
}

#' NORTA Observations
#' Generates covariate matrix using the normal-to-anything algorithm
#' @param corr_mat symmetric, positive definite correlation matrix
#' @param cdf cumulative distribution function of a random variable
#' @param num_rows number of observations, i.e., rows
#' @return (nxp) matrix of observations
#' @export
generate_matrix_of_observations <- function(corr_mat, cdf,
                                            num_rows = 1){
  lower_chol <- t(chol(corr_mat))
  k <- dim(lower_chol)[1]

  res <- matrix(0, num_rows, k)

  for(i in 1:num_rows){
    W <- rnorm(k)
    Z <- lower_chol %*% W

    for(j in 1:k){
      res[i,j] <- determine_x(Z[j], cdf)
    }
  }
  res
}
