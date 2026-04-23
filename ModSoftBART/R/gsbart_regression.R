#' General SoftBart Regression with group handling
#' 
#' Fits a semiparametric regression model with grouped predictors and linear covariates using a SoftBart model.
#'
#' @param Y A vector of outcome data.
#' @param Z A matrix of predictor variables to be included in the exposure-response function.
#' @param X A matrix of linear covariate included in the mean response.
#' @param Znew Optional matrix of new predictor values.
#' @param Xnew Optional matrix of new covariate values.
#' @param num_tree The number of trees in the ensemble to use.
#' @param k Determines the standard deviation of the leaf node parameters, which is given by \code{3 / k / sqrt(num_tree)}.
#' @param alpha_val Positive constant controlling the sparsity level among groups.
#' @param alpha_comp Positive constant controlling the sparsity level across group components.
#' @param alpha_scale Scale of the prior for \code{alpha}.
#' @param alpha_comp_scale Scale of the prior for \code{alpha_comp}.
#' @param num_burn Number of warmup iterations for the chain.
#' @param num_thin Thinning interval for the chain.
#' @param num_save The number of samples to collect; in total, \code{num_burn + num_save * num_thin} iterations are run.
#' @param hypers A list of hyperparameters constructed from the \code{Hypers()} function (\code{num_tree}, \code{k}, \code{alpha_val}, \code{alpha_comp}, \code{alpha_scale} and \code{alpha_comp_scale} are overridden by this function).
#' @param opts A list of options for running the chain constructed from the \code{Opts()} function (\code{num_burn}, \code{num_thin}, and \code{num_save} are overridden by this function).
#' @param group Optional vector of group indicators used for hierarchical variable selection. If not specified, variable selection is performed at the component level.
#' @param include_intercept If \code{FALSE} then any intercept term in the linear formula will be removed, with the overall location of the outcome captured by the nonparametric function.


#' @return Returns a list with the following components
#' \itemize{
#'   \item \code{h_train}: samples of the nonparametric function evaluated on the training set.
#'   \item \code{h_test}: samples of the nonparametric function evaluated on the test set.
#'   \item \code{eta_train}: samples of the linear predictor on the training set.
#'   \item \code{eta_test}: samples of the linear predictor on the test set.
#'   \item \code{beta_hat}: samples of the regression coefficients.
#'   \item \code{sigma}: samples of the error standard deviation.
#'   \item \code{sigma_mu}: samples of the standard deviation of the leaf node parameters.
#'   \item \code{var_counts}: a matrix with a column for each nonparametric predictor containing the number of times that predictor is used in the ensemble at each iteration.
#'   \item \code{tree_counts}: a list of counts at tree level. If \code{group} is \code{NULL}, each element is a matrix with one column per tree and one row per nonparametric predictor, giving the number of times each predictor is used in each tree. If \code{group} is specified, each element contains two matrices: one with group counts by tree and one with predictor counts by tree.
#'   \item \code{opts}: the options used when running the chain.
#'   \item \code{hypers}: the hypers used when running the chain.
#' }
#' @export
#'




gsbart_regression <- function(X, Y, Z, Znew = Z, Xnew = X, 
                              num_tree = 50, k = 2, alpha_val = 1.0, alpha_comp = 1.0,
                              alpha_scale = NULL, alpha_comp_scale = NULL,
                              num_burn = 5000, num_thin = 10, num_save = 500,
                              hypers = NULL, opts = NULL, group = NULL,
                              include_intercept = FALSE){
  
  
  # Define Bayesian Linear regression function
  blr <- function(X, Y, sigma_0, intercept = TRUE, update_sigma = TRUE){
    
    fit_lm <- lm(Y ~ X)
    if(!intercept){
      fit_lm <- lm(Y ~ X - 1)
    }
    
    beta_hat <- coef(fit_lm)
    sigma_hat <- sigma(fit_lm)
    V <- vcov(fit_lm) / sigma_hat^2 * sigma_0^2
    
    beta <- MASS::mvrnorm(n = 1, mu = beta_hat, Sigma = V)
    
    mu_Y <- as.numeric(model.matrix(fit_lm) %*% beta)
    
    if(update_sigma){
      RSS <- sum((Y-mu_Y)^2)
      prec <- rgamma(1, length(Y) / 2, RSS / 2)
      sigma_0 <- 1 / sqrt(prec)
    }
    
    return(list(beta = beta, mu_Y = mu_Y, sigma_0 = sigma_0))
    
  }
  
  if(include_intercept){
    Xnew <- model.matrix(~., data = as.data.frame(Xnew))
  }
  
  
  ## Normalize!
  make_01_norm <- function(x){
    a <- min(x)
    b <- max(x)
    return(function(y) (y - a) / (b - a))
  }
  
  ecdfs     <- list()
  Z_train   <- as.matrix(Z)
  Z_test    <- as.matrix(Znew)
  Z_train_q <- Z_train
  Z_test_q  <- Z_test
  
  for (i in 1:ncol(Z_train)) {
    ecdfs[[i]]     <- ecdf(Z_train[,i])
    if(length(unique(Z_train[,i]) == 1)) ecdfs[[i]] <- identity
    if(length(unique(Z_train[,i]) == 2)) ecdfs[[i]] <- make_01_norm(Z_train[,i])
  }
  
  for (i in 1:ncol(Z_train)) {
    Z_train_q[,i]  <- ecdfs[[i]](Z_train[,i])
    Z_test_q[,i]   <- ecdfs[[i]](Z_test[,i])
    
  }
  
  
  ## construct the forest
  mu_Ytrain <- mean(Y)
  sd_Ytrain <- sd(Y)
  y_scale <- (Y - mu_Ytrain) / sd_Ytrain
  
  # run Soft BART
  if(!is.null(hypers)){
    hypers <- hypers
    hypers$X <- Z_train_q
    hypers$Y <- y_scale
    hypers$normalize_Y <- FALSE
  }else{
    hypers <- Hypers(X = Z_train_q, 
                     Y = y_scale, 
                     group = group,
                     k = k,
                     alpha = alpha_val,
                     alpha_comp = alpha_comp,
                     alpha_scale = alpha_scale, 
                     alpha_comp_scale = alpha_comp_scale,
                     num_tree = num_tree, 
                     normalize_Y = FALSE,
                     sigma_hat = 1
    )
  }
  
  hypers$sigma_mu <- 3 / k / sqrt(num_tree)
  hypers$num_tree <- num_tree
  
  
  if(!is.null(opts)){
    opts <- opts
    
  }else{
    opts   <- Opts(update_s = TRUE,
                   num_burn = num_burn,
                   num_thin = num_thin,
                   num_save = num_save)
  }
  
  opts$num_print <- .Machine$integer.max
  num_iter <- opts$num_burn + opts$num_thin * opts$num_save
  
  mk_forest <- MakeForest(hypers = hypers, opts = opts, warn = FALSE)
  
  # initialize output
  beta_samp <- matrix(nrow = opts$num_save, ncol = ncol(Xnew))
  h_samp_train <- matrix(nrow = opts$num_save, ncol = nrow(Z_train))
  h_samp_test <- matrix(nrow = opts$num_save, ncol = nrow(Z_test))
  eta_samp_train <- matrix(nrow = opts$num_save, ncol = nrow(Z_train))
  eta_samp_test <- matrix(nrow = opts$num_save, ncol = nrow(Z_test))
  sigma_out <- numeric(opts$num_save)
  sigma_mu_out <- numeric(opts$num_save)
  
  varcounts <- matrix(nrow = opts$num_save, ncol = ncol(Z_train)) 
  treecounts  <- vector("list", length = opts$num_save)
  
  
  ## The chain ----
  r <- as.numeric(mk_forest$do_predict(Z_train_q))
  
  for (t in 1:num_iter) {
    
    c(beta, mu_y, sigma_Y) %<-%
      blr(X, Y = y_scale - r, sigma_0 = mk_forest$get_sigma(), 
          intercept = include_intercept, update_sigma = TRUE)
    
    # update the forest
    temp <- y_scale - mu_y
    r <- as.numeric(mk_forest$do_gibbs(Z_train_q, temp, Z_train_q, 1))
    
    h_hat <- as.numeric(mk_forest$do_predict(Z_test_q))
    
    # save posterior samples
    if(t > opts$num_burn){
      if((t-opts$num_burn) %% opts$num_thin == 0){
        tt <- (t - opts$num_burn) / opts$num_thin
        
        h_samp_train[tt, ] <- r * sd_Ytrain + mu_Ytrain
        h_samp_test[tt,] <- h_hat * sd_Ytrain + mu_Ytrain
        eta_samp_train[tt,] <- mu_y * sd_Ytrain
        eta_samp_test[tt,] <- as.numeric(Xnew %*% beta) * sd_Ytrain
        beta_samp[tt,] <- beta * sd_Ytrain
        sigma_out[tt] <- mk_forest$get_sigma() * sd_Ytrain
        sigma_mu_out[tt] <- mk_forest$get_sigma_mu() * sd_Ytrain
        if(is.null(group)){
          varcounts[tt, ] <- as.numeric(mk_forest$get_counts())
          treecounts[[tt]] <- mk_forest$get_tree_counts()
        }else{
          varcounts[tt, ] <- as.numeric(rowSums(mk_forest$get_tree_counts()$Variable_tree_counts))
          treecounts[[tt]] <- mk_forest$get_tree_counts()
        }
        
        
      }
      
    }
    
  }
  
  out <- list(h_train = h_samp_train, 
              h_test = h_samp_test, 
              eta_train = eta_samp_train,
              eta_test = eta_samp_test, 
              beta_hat = beta_samp, 
              sigma = sigma_out,
              sigma_mu = sigma_mu_out, 
              var_counts = varcounts, 
              tree_counts = treecounts,
              hypers = hypers, 
              opts = opts)
  
  return(out)
  
}

