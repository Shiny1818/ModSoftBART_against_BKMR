
################################################################
################### Simulation Example #########################
SimData_v2 <- function (n = 100, M = 15, sigsq.true = 0.5, beta.true = c(.5, .2), hfun = 1, 
                        Zgen = "norm", indx = c(1:2, 4:6), family = "gaussian") {
  # revised HFun to incorporate more variables
  
  sfunc <- function(x,a = 1, b=.3){10*(2/(1+exp(-a*x)) - b*x-1)}
  HFun1 <- function(z, ind = indx){
    # u <- (.5*(z[ind[1]] + z[ind[2]]) + .3*(z[ind[3]] + z[ind[4]] + z[ind[5]]))/.15
    u <- (.5*(z[ind[1]] + z[ind[2]]) + .3*(z[ind[3]] + z[ind[4]] + z[ind[5]]))/1
    return(sfunc(x=u))
  }
  
  
  HFun3 <- function(z, ind = indx){
    u <- 1/6 * (z[ind[1]] + z[ind[2]] + 2*z[ind[3]] + 2*z[ind[4]] + 3*z[ind[5]] +
                  1/2 * z[ind[1]] * z[ind[2]] + 1/3 * z[ind[3]] *z[ind[4]] *z[ind[5]])
    
    return(sfunc(x=u))
  }
  
  
  
  HFun2 <- function(z, ind = indx){
    1/6 * (z[ind[1]] + z[ind[2]] + 2*z[ind[3]] + 2*z[ind[4]] + 3*z[ind[5]] +
             1/2 * z[ind[1]] * z[ind[2]] + 1/3 * z[ind[3]] *z[ind[4]] *z[ind[5]])
  }
  
  HFun_linear <- function(z, ind = indx){
    1/6 * (z[ind[1]] + z[ind[2]] + 2*z[ind[3]] + 2*z[ind[4]] + 3*z[ind[5]])
  }
  
  HFun_linearInter <- function(z, ind = indx){
    1/6 * (z[ind[1]] + z[ind[2]] + 2*z[ind[3]] + 2*z[ind[4]] + 3*z[ind[5]] +
             1/2 * z[ind[1]] * z[ind[2]])
  }
  
  ## run data generation
  stopifnot(n > 0, M > 0, sigsq.true >= 0, family %in% c("gaussian", 
                                                         "binomial"))
  if (family == "binomial") {
    sigsq.true <- 1
  }
  
  
  if (hfun == 1) {
    HFun <- HFun1
  }else if (hfun == 2) {
    HFun <- HFun2
  }else if (hfun == 3) {
    HFun <- HFun3
  }else if (hfun == "linear") {
    HFun <- HFun_linear
  }else if (hfun == "linearInter"){
    HFun <- HFun_linearInter
  }else{
    stop("hfun must be an integer from 1 to 3")
  }
  if (Zgen == "unif") {
    Z <- matrix(runif(n * M, -2, 2), n, M)
  }else if (Zgen == "norm") {
    Z <- matrix(rnorm(n * M), n, M)
  }else if (Zgen == "corr") {
    if (M < 3) {
      stop("M must be an integer > 2 for Zgen = 'corr'")
    }
    Sigma <- diag(1, M, M)
    if(M >= 15){
      # group1: Z1, Z2, Z3
      Sigma[1, 3] <- Sigma[3, 1] <- Sigma[2, 3] <- Sigma[3, 2] <- Sigma[1, 2] <- Sigma[2, 1] <- 0.6
      # group2: Z4, Z5, Z6, Z7
      Sigma[4, 5] <- Sigma[5, 4] <- Sigma[4, 6] <- Sigma[6, 4] <- Sigma[4, 7] <- Sigma[7, 4] <- 0.3
      Sigma[5, 6] <- Sigma[6, 5] <- Sigma[5, 7] <- Sigma[7, 5] <- Sigma[6, 7] <- Sigma[7, 6] <- 0.3
      # group3: Z8 through Z12
      Sigma[8, 9] <- Sigma[9, 8] <- Sigma[8, 10] <- Sigma[10, 8] <- Sigma[8, 11] <- Sigma[11, 8] <- 0.75
      Sigma[8, 12] <- Sigma[12, 8] <- Sigma[9, 10] <- Sigma[10, 9] <- Sigma[9, 11] <- Sigma[11, 9] <- 0.75
      Sigma[9, 12] <- Sigma[12, 9] <- Sigma[10, 11] <- Sigma[11, 10] <- Sigma[10, 12] <- Sigma[12, 10] <- 0.75
      Sigma[11, 12] <- Sigma[12, 11] <- .75
      # group4: Z13 through Z15
      Sigma[13, 14] <- Sigma[14, 13] <- Sigma[13, 15] <- Sigma[15, 13] <- Sigma[14, 15] <- Sigma[15, 14] <- 0.5
      
    }else{
      Sigma[1, 3] <- Sigma[3, 1] <- 0.95
      Sigma[2, 3] <- Sigma[3, 2] <- 0.3
      Sigma[1, 2] <- Sigma[2, 1] <- 0.1
      
    }
    
    
    Z <- MASS::mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
  }else if (Zgen == "realistic") {
    VarRealistic <- structure(c(0.72, 0.65, 0.45, 0.48, 0.08, 
                                0.14, 0.16, 0.42, 0.2, 0.11, 0.35, 0.1, 0.11, 0.65, 
                                0.78, 0.48, 0.55, 0.06, 0.09, 0.17, 0.2, 0.16, 0.11, 
                                0.32, 0.12, 0.12, 0.45, 0.48, 0.56, 0.43, 0.11, 0.15, 
                                0.23, 0.25, 0.28, 0.16, 0.31, 0.15, 0.14, 0.48, 0.55, 
                                0.43, 0.71, 0.2, 0.23, 0.32, 0.22, 0.29, 0.14, 0.3, 
                                0.22, 0.18, 0.08, 0.06, 0.11, 0.2, 0.95, 0.7, 0.45, 
                                0.22, 0.29, 0.16, 0.24, 0.2, 0.13, 0.14, 0.09, 0.15, 
                                0.23, 0.7, 0.8, 0.36, 0.3, 0.35, 0.13, 0.23, 0.17, 
                                0.1, 0.16, 0.17, 0.23, 0.32, 0.45, 0.36, 0.83, 0.24, 
                                0.37, 0.2, 0.36, 0.34, 0.25, 0.42, 0.2, 0.25, 0.22, 
                                0.22, 0.3, 0.24, 1.03, 0.41, 0.13, 0.39, 0.1, 0.1, 
                                0.2, 0.16, 0.28, 0.29, 0.29, 0.35, 0.37, 0.41, 0.65, 
                                0.18, 0.3, 0.18, 0.16, 0.11, 0.11, 0.16, 0.14, 0.16, 
                                0.13, 0.2, 0.13, 0.18, 0.6, 0.18, 0.13, 0.08, 0.35, 
                                0.32, 0.31, 0.3, 0.24, 0.23, 0.36, 0.39, 0.3, 0.18, 
                                0.79, 0.42, 0.12, 0.1, 0.12, 0.15, 0.22, 0.2, 0.17, 
                                0.34, 0.1, 0.18, 0.13, 0.42, 1.27, 0.1, 0.11, 0.12, 
                                0.14, 0.18, 0.13, 0.1, 0.25, 0.1, 0.16, 0.08, 0.12, 
                                0.1, 0.67), .Dim = c(13L, 13L))
    if (M > ncol(VarRealistic)) {
      stop("Currently can only generate exposure data based on a realistic correlation structure with M = 13 or fewer. Please set M = 13 or use Zgen = c('unif','norm'")
    }else if (M <= 13) {
      Sigma <- VarRealistic[1:M, 1:M]
    }
    Z <- MASS::mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
  }
  colnames(Z) <- paste0("z", 1:M)
  # X <- cbind(3 * cos(Z[, 1]) + 2 * rnorm(n))
  X <- cbind(rnorm(n), rbinom(n,1,.5))
  eps <- rnorm(n, sd = sqrt(sigsq.true))
  h <- apply(Z, 1, HFun)
  
  mu <- X %*% beta.true + h
  y <- drop(mu + eps)
  if (family == "binomial") {
    ystar <- y
    y <- ifelse(ystar > 0, 1, 0)
  }
  dat <- list(n = n, M = M, sigsq.true = sigsq.true, beta.true = beta.true, 
              Z = Z, h = h, X = X, y = y, hfun = hfun, HFun = HFun, mu = mu,
              family = family)
  if (family == "binomial") {
    dat$ystar <- ystar
  }
  dat
}

#########################################################################








