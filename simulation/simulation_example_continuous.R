# Single-run simulation example for the paper.
# The full simulation study in the paper is based on 500 replicated runs.


library(bkmr)
library(ModSoftBart)
library(progress)
library(tictoc)
library(dplyr)
library(ggplot2)
library(stringr)
library(zeallot)
library(caret)
library(lubridate)

set.seed(19816)
SimData_v2 <- function (n = 100, M = 15, sigsq.true = 0.5, beta.true = c(.5, .2), hfun = 3, 
                        Zgen = "norm", indx = c(1:2, 4:6), family = "gaussian") {
  
  sfunc <- function(x,a = 1, b=.3){10*(2/(1+exp(-a*x)) - b*x-1)}
  HFun1 <- function(z, ind = indx){
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


nsize = 500; msize = 15; hfunc = 1; Zgenr = "corr"; index = c(1:2, 4:6)

dat_train <- SimData_v2(n = nsize, M = msize, hfun = hfunc, Zgen = Zgenr, indx = index)
Ytr   <- dat_train$y
Ztr   <- dat_train$Z
Xtr   <- dat_train$X

dat_test <- SimData_v2(n = nsize, M = msize, hfun = hfunc, Zgen = Zgenr, indx = index)
Ytest <- dat_test$y
Ztest <- dat_test$Z
Xtest <- dat_test$X

## BKMR
start_time1 <- now()
fitkm <- kmbayes(y=Ytr, Z=Ztr, X=Xtr, iter=10000, verbose=TRUE, varsel=TRUE, 
                 groups=c(rep(1,times=3), rep(2,times=4), rep(3,times=5),rep(4,times=3)))

h_est_train <- SamplePred(fitkm, Znew = Ztr, Xnew = matrix(rep(0, nrow(Ztr)*ncol(Xtr)), ncol = ncol(Xtr)))
h_train_bkmr <- colMeans(h_est_train)
h_est_val <- SamplePred(fitkm, Znew = Ztest, Xnew = matrix(rep(0, nrow(Ztest)*ncol(Xtest)), ncol = ncol(Xtest)))
h_test_bkmr <- colMeans(h_est_val)

time_taken1 <- as.duration(now() - start_time1)

## General soft BART (gsBART)
start_time2 <- now()

fit_gsbart <- gsbart_regression(X = Xtr, Y = Ytr, Z = Ztr, Znew = Ztest, Xnew = as.matrix(Xtest),
                                num_burn = 5000, num_save = 500, num_thin = 10, alpha_val = 1e-30, alpha_comp = 1e-10,
                                group = c(rep(1,times=3), rep(2,times=4), rep(3,times=5),rep(4,times=3)))

h_train_soft <- colMeans(fit_gsbart$h_train)
h_test_soft <- colMeans(fit_gsbart$h_test)

time_taken2 <- as.duration(now() - start_time2)

################################
## extract PIPs for both models
### BKMR
PIPs_bkmr <- ExtractPIPs(fitkm)

### gsBART
PIPs_gsbart <- calPIPs_gsbart(fit = fit_gsbart, group = TRUE, vars = Ztr)


## evaluation based on the regression of h_hat on h
h_true_train <- dat_train$h
h_true_test <- dat_test$h

hfit_bkmr_train <- lm(h_train_bkmr~h_true_train)

c(inter_bkmr_train, slope_bkmr_train) %<-% coef(hfit_bkmr_train)
rs_bkmr_train <- summary(hfit_bkmr_train)$r.squared
se_bkmr_train <- summary(hfit_bkmr_train)$sigma


hfit_bkmr_test <- lm(h_test_bkmr~h_true_test)

c(inter_bkmr_test, slope_bkmr_test) %<-% coef(hfit_bkmr_test)
rs_bkmr_test <- summary(hfit_bkmr_test)$r.squared
se_bkmr_test <- summary(hfit_bkmr_test)$sigma

performance_bkmr <- data.frame(BKMR_train = c(inter_bkmr_train, slope_bkmr_train,
                                              rs_bkmr_train, se_bkmr_train),
                               BKMR_test = c(inter_bkmr_test, slope_bkmr_test,
                                             rs_bkmr_test, se_bkmr_test))

row.names(performance_bkmr) <- c("Intercept", "Slope", "R-squared", "StdError")


hfit_bart_train <- lm(h_train_soft~h_true_train)

c(inter_bart_train, slope_bart_train) %<-% coef(hfit_bart_train)
rs_bart_train <- summary(hfit_bart_train)$r.squared
se_bart_train <- summary(hfit_bart_train)$sigma

hfit_test_bart <- lm(h_test_soft~h_true_test)

c(inter_bart_test, slope_bart_test) %<-% coef(hfit_test_bart)
rs_bart_test <- summary(hfit_test_bart)$r.squared
se_bart_test <- summary(hfit_test_bart)$sigma

performance_softbart <- data.frame(SoftBART_train = c(inter_bart_train, slope_bart_train,
                                                      rs_bart_train, se_bart_train),
                                   SoftBART_test = c(inter_bart_test, slope_bart_test,
                                                     rs_bart_test, se_bart_test))

row.names(performance_softbart) <- c("Intercept", "Slope", "R-squared", "StdError")


print(performance_bkmr)
print(performance_softbart)

###############################################################################
# GAM approximate the posterior samples.
library(mgcv)
library(dplyr)
library(ggplot2)

## BKMR posterior samples
h_est_val <- SamplePred(fitkm, Znew = Ztest, Xnew = matrix(rep(0, nrow(Ztest)*ncol(Xtest)), ncol = ncol(Xtest)))

htest_bkmr_sampl <- t(as.matrix(h_est_val))
htest_bkmr_mean <- as.matrix(colMeans(h_est_val))

GAMfit_bkmr <- gam(htest_bkmr_mean ~
                     s(z1) + 
                     s(z2) +
                     s(z3) + 
                     s(z4) + 
                     s(z5)+ 
                     s(z6) + 
                     s(z7) + 
                     s(z8)+
                     s(z9) + 
                     s(z10) + s(z11) + s(z12) + s(z13) + s(z14) + s(z15),
                   data = as.data.frame(data_test$Z))

# model matrix and covariance matrix
Xs_bkmr <- model.matrix(GAMfit_bkmr)
V_bkmr <- vcov(GAMfit_bkmr, dispersion = 1)

## Point summary and posterior draws of projection
q_bkmr <- crossprod(V_bkmr, crossprod(Xs_bkmr, htest_bkmr_mean))
Q_bkmr <- crossprod(V_bkmr, crossprod(Xs_bkmr, htest_bkmr_sampl))

## Posterior draws of the summary
gammaSamples_bkmr <- Xs_bkmr %*% Q_bkmr


## gsBART posterior samples
htest_gsbart_sampl <- t(as.matrix(fit_gsbart$h_test))
htest_gsbart_mean  <- as.matrix(colMeans(fit_gsbart$h_test))

GAMfit_gsbart <- gam(htest_gsbart_mean ~
                       s(z1) + 
                       s(z2) +
                       s(z3) + 
                       s(z4) + 
                       s(z5)+ 
                       s(z6) + 
                       s(z7) + 
                       s(z8)+
                       s(z9) + 
                       s(z10) + s(z11) + s(z12) + s(z13) + s(z14) + s(z15),
                     data = as.data.frame(data_test$Z))

Xs_gsbart <- model.matrix(GAMfit_gsbart)
V_gsbart <- vcov(GAMfit_gsbart, dispersion = 1)

q_gsbart <- crossprod(V_gsbart, crossprod(Xs_gsbart, htest_gsbart_mean))
Q_gsbart <- crossprod(V_gsbart, crossprod(Xs_gsbart, htest_gsbart_sampl))

gammaSamples_gsbart <- Xs_gsbart %*% Q_gsbart

# plot smooth items
gam.pred.bkmr <- predict(GAMfit_bkmr, newdata = as.data.frame(Ztest), type = "terms",
                         se.fit = TRUE)
gam.pred.gsbart <- predict(GAMfit_gsbart, newdata = as.data.frame(Ztest), type = "terms",
                           se.fit = TRUE)

