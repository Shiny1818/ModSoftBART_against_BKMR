#' Create a list of hyperparameter values
#'
#' Creates a list which holds all the hyperparameters for use with the
#' model-fitting functions and with the \code{MakeForest} functionality.
#'
#' @param X A matrix of training data covariates.
#' @param Y A vector of training data responses.
#' @param group group indicators. For each column of \code{X}, \code{group} gives the associated group.
#' @param alpha Positive constant controlling the sparsity level among groups if group was specified.
#' @param alpha_comp Positive constant controlling the sparsity level across components if group exists.
#' @param beta Parameter penalizing tree depth in the branching process prior.
#' @param gamma Parameter penalizing new nodes in the branching process prior.
#' @param k Related to the signal-to-noise ratio, \code{sigma_mu = 0.5 / (sqrt(num_tree) * k)}. BART defaults to \code{k = 2} after applying the max/min normalization to the outcome.
#' @param sigma_hat A prior guess at the conditional variance of \code{Y} given \code{X}. If not provided, this is estimated empirically by linear regression.
#' @param shape Shape parameter for gating probabilities.
#' @param width Bandwidth of gating probabilities.
#' @param num_tree Number of trees in the ensemble.
#' @param alpha_scale Scale of the prior for \code{alpha}; if not provided, defaults to the number of groups.
#' @param alpha_comp_scale Scale of the prior for \code{alpha}; if not provided, defaults to the number of predictors in each group.
#' @param alpha_shape_1 Shape parameter for prior on \code{alpha}; if not provided, defaults to 0.5.
#' @param alpha_shape_2 Shape parameter for prior on \code{alpha}; if not provided, defaults to 1.0.
#' @param alpha_comp_shape_1 Shape parameter for prior on \code{alpha_comp}; if not provided, defaults to 0.5.
#' @param alpha_comp_shape_2 Shape parameter for prior on \code{alpha_comp}; if not provided, defaults to 1.0.
#' @param tau_rate Rate parameter for the bandwidths of the trees with an exponential prior; defaults to 10.
#' @param num_tree_prob Parameter for geometric prior on number of tree.
#' @param temperature The temperature applied to the posterior distribution; set to 1 unless you know what you are doing.
#' @param weights Only used by the function \code{softbart}, this is a vector of weights to be used in heteroskedastic regression models, with the variance of an observation given by \code{sigma_sq / weight}.
#' @param normalize_Y Do you want to compute \code{sigma_hat} after applying the standard BART max/min normalization to \eqn{(-0.5, 0.5)} for the outcome? If \code{FALSE}, no normalization is applied. This might be useful for fitting custom models where the outcome is normalized by hand.
#'
#' @return Returns a list containing the function arguments.
Hypers <- function(X,Y, group = NULL, alpha = 1, alpha_comp = 1, beta = 2, gamma = 0.95, k = 2,
                   sigma_hat = NULL, shape = 1, width = 0.1, num_tree = 20,
                   alpha_scale = NULL, alpha_comp_scale = NULL, alpha_shape_1 = 0.5,
                   alpha_shape_2 = 1, alpha_comp_shape_1 = .5, alpha_comp_shape_2 = 1, 
                   tau_rate = 10, num_tree_prob = NULL,
                   temperature = 1.0, weights = NULL, normalize_Y = TRUE) {

  # if(is.null(alpha_scale)) alpha_scale <- ncol(X)
  alpha_scale <- ifelse(is.null(alpha_scale), ifelse(is.null(group), ncol(X), max(group)), alpha_scale)
  if(is.null(group)){
    alpha_comp_scale <- as.vector(rep(ncol(X), ncol(X)))
  }else{
    if(length(alpha_comp_scale) != max(group)) alpha_comp_scale <- as.vector(table(group))
  }
  if(is.null(num_tree_prob)) num_tree_prob <- 2.0 / num_tree
  if(is.null(weights)) weights <- rep(1, length(Y))

  out                                  <- list()
  out$weights                          <- weights
  out$alpha                            <- alpha
  out$alpha_comp                       <- alpha_comp
  out$beta                             <- beta
  out$gamma                            <- gamma
  out$sigma_mu                         <- 0.5 / (k * sqrt(num_tree))
  out$k                                <- k
  out$num_tree                         <- num_tree
  out$shape                            <- shape
  out$width                            <- width
  if(is.null(group)) {
    out$group                          <- 1:ncol(X) - 1
  } else {
    out$group                          <- group - 1
  }
  

  if(normalize_Y) {
    Y                                  <- normalize_bart(Y)
  }
  if(is.null(sigma_hat))
    sigma_hat                          <- GetSigma(X,Y, weights = weights)

  out$sigma                            <- sigma_hat
  out$sigma_hat                        <- sigma_hat

  out$alpha_scale                      <- alpha_scale
  out$alpha_shape_1                    <- alpha_shape_1
  out$alpha_shape_2                    <- alpha_shape_2
  out$alpha_comp_scale                 <- alpha_comp_scale
  out$alpha_comp_shape_1               <- alpha_comp_shape_1
  out$alpha_comp_shape_2               <- alpha_comp_shape_2
  out$tau_rate                         <- tau_rate
  out$num_tree_prob                    <- num_tree_prob
  out$temperature                      <- temperature

  return(out)

}

#' MCMC options for SoftBart
#'
#' Creates a list that provides the parameters for running the Markov chain.
#'
#' @param num_burn Number of warmup iterations for the chain.
#' @param num_thin Thinning interval for the chain.
#' @param num_save The number of samples to collect; in total, \code{num_burn + num_save * num_thin} iterations are run.
#' @param num_print Interval for how often to print the chain's progress.
#' @param update_sigma_mu If \code{TRUE}, \code{sigma_mu} is  updated, with a half-Cauchy prior on \code{sigma_mu} centered at the initial guess.
#' @param update_sigma If \code{TRUE}, \code{sigma} is updated, with a half-Cauchy prior on \code{sigma} centered at the initial guess.
#' @param update_s If \code{TRUE}, \code{s} is updated using the Dirichlet prior \eqn{s \sim D(\alpha / P, \ldots, \alpha / P)} where \eqn{P} is the number of covariates.
#' @param update_alpha If \code{TRUE}, \code{alpha} and \code{alpha_comp} are updated using a scaled beta prime prior.
#' @param update_beta If \code{TRUE}, \code{beta} is updated using a normal prior with mean 0 and variance 4.
#' @param update_gamma If \code{TRUE}, gamma is updated using a Uniform(0.5, 1) prior.
#' @param update_tau If \code{TRUE}, the bandwidth \code{tau} is updated for each tree
#' @param update_tau_mean If \code{TRUE}, the mean of \code{tau} is updated
#' @param cache_trees If \code{TRUE}, we save the trees for each MCMC iteration when using the MakeForest interface
#'
#' @return Returns a list containing the function arguments.
Opts <- function(num_burn = 2500, num_thin = 1, num_save = 2500, num_print = 100,
                 update_sigma_mu = TRUE, update_s = TRUE, update_alpha = TRUE,
                 update_beta = FALSE, update_gamma = FALSE, update_tau = TRUE,
                 update_tau_mean = FALSE, update_sigma = TRUE,
                 cache_trees = TRUE) {
  out <- list()
  out$num_burn        <- num_burn
  out$num_thin        <- num_thin
  out$num_save        <- num_save
  out$num_print       <- num_print
  out$update_sigma_mu <- update_sigma_mu
  out$update_s         <- update_s
  out$update_alpha    <- update_alpha
  out$update_beta     <- update_beta
  out$update_gamma    <- update_gamma
  out$update_tau      <- update_tau
  out$update_tau_mean <- update_tau_mean
  # out$update_num_tree <- update_num_tree
  out$update_num_tree <- FALSE
  out$update_sigma    <- update_sigma
  out$cache_trees     <- cache_trees

  return(out)

}

normalize_bart <- function(y) {
  a <- min(y)
  b <- max(y)
  y <- (y - a) / (b - a) - 0.5
  return(y)
}

unnormalize_bart <- function(z, a, b) {
  y <- (b - a) * (z + 0.5) + a
  return(y)
}


GetSigma <- function(X,Y, weights = NULL) {
  
  if(is.null(weights)) weights <- rep(1, length(Y))
  stopifnot(is.matrix(X) | is.data.frame(X))

  if(is.data.frame(X)) {
    X <- model.matrix(~.-1, data = X)
  }


  fit <- cv.glmnet(x = X, y = Y, weights = weights)
  fitted <- predict_glmnet(fit, X)
  sigma_hat <- sqrt(mean((fitted - Y)^2))
  # sigma_hat <- 0
  # if(nrow(X) > 2 * ncol(X)) {
  #   fit <- lm(Y ~ X)
  #   sigma_hat <- summary(fit)$sigma
  # } else {
  #   sigma_hat <- sd(Y)
  # }

  return(sigma_hat)

}

predict_glmnet <- function (object, newx, s = c("lambda.1se", "lambda.min"), ...) {
  if (is.numeric(s)) 
    lambda = s
  else if (is.character(s)) {
    s = match.arg(s)
    lambda = object[[s]]
    names(lambda) = s
  }
  else stop("Invalid form for s")
  predict(object$glmnet.fit, newx, s = lambda, ...)
}


