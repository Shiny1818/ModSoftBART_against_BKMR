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

set.seed(32517)
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

dat_train <- SimData_v2(n = nsize, M = msize, sigsq.true = 1, hfun = hfunc, 
                        Zgen = Zgenr, indx = index, family = "binomial")
Ytr   <- dat_train$y
Ztr   <- dat_train$Z
Xtr   <- dat_train$X

dat_test <- SimData_v2(n = nsize, M = msize, sigsq.true = 1, hfun = hfunc, 
                       Zgen = Zgenr, indx = index, family = "binomial")
Ytest <- dat_test$y
Ztest <- dat_test$Z
Xtest <- dat_test$X

h_true_train <- dat_train$h
h_true_test <- dat_test$h

start_time1 <- now()
fitkm.pb <- try(kmbayes(y = Ytr, Z = Ztr, X = Xtr, iter = 10000, verbose = FALSE, varsel = TRUE, 
                        groups = group, family = "binomial", 
                        control.params = list(r.jump2 = 0.5), est.h = TRUE), silent = TRUE)

if(!inherits(fitkm.pb, "try-error")){
  
  h_est_train <- SamplePred(fitkm.pb, Znew = Ztr, Xnew = matrix(rep(0, nrow(Ztr)*ncol(Xtr)), ncol = ncol(Xtr)))
  h_train_bkmr <- colMeans(h_est_train)
  h_est_val <- SamplePred(fitkm.pb, Znew = Ztest, Xnew = matrix(rep(0, nrow(Ztest)*ncol(Xtest)), ncol = ncol(Xtest)))
  h_test_bkmr <- colMeans(h_est_val)
  
  time_taken1 <- as.duration(now() - start_time1)
  
  ##############################################################################
  ## calculate PIPs
  PIPs_bkmr.pb <- ExtractPIPs(fitkm.pb)
  ##############################################################################
  
  # performance check for BKMR estimates
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
  
}else{
  error_msg <- attr(fitkm.pb, "condition")$message
  log_errors(nseed, error_msg, pid = pids, results_folder = results_folder)
  
  PIPs_bkmr <- NULL
  performance_bkmr <- NULL
}

start_time2 <- now()
fit_gsbart.pb <- pb_gsbart(X = Xtr, Y = Ytr, Z = Ztr, Znew = Ztest, Xnew = as.matrix(Xtest),
                           num_burn = num_burn, num_save = num_save, num_thin = num_thin, 
                           alpha_val = alpha_val, alpha_comp = alpha_comp, 
                           alpha_scale = alpha_scale, alpha_comp_scale = alpha_comp_scale,
                           group = group, num_tree = 20)


h_train_soft <- colMeans(fit_gsbart.pb$h_train)
h_test_soft <- colMeans(fit_gsbart.pb$h_test)

time_taken2 <- as.duration(now() - start_time2)

PIPs_gsbart.pb <- calPIPs_gsbart(fit = fit_gsbart.pb, group = group, vars = Ztr)

# performance check for soft BART estimates
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

##########################################
# ggplot GAM results.

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(stringr)


sfunc <- function(u,a=1, b = .3){
  10*(2/(1+exp(-a*u)) - b*u - 1)
}

####### For h1
Href <- function(z, index = c(1:2, 4:6), idx = 1){
  
  if(idx ==1){
    u = .5*(z[idx] + Zquant[index[2]]) + .3*(Zquant[index[3]] + Zquant[index[4]] + Zquant[index[5]])
    # 
  }else if(idx ==2){
    u = .5*(Zquant[index[1]] + z[idx]) + .3*(Zquant[index[3]] + Zquant[index[4]] + Zquant[index[5]])
    # 
  }else if(idx ==4){
    u = .5*(Zquant[index[1]] + Zquant[index[2]]) + .3*(z[idx] + Zquant[index[4]] + Zquant[index[5]])
    # 
  }else if(idx ==5){
    u = .5*(Zquant[index[1]] + Zquant[index[2]]) + .3*(Zquant[index[3]] + z[idx] + Zquant[index[5]])
    # 
  }else if(idx ==6){
    u = .5*(Zquant[index[1]] + Zquant[index[2]]) + .3*(Zquant[index[3]] + Zquant[index[4]] + z[idx])
    # 
  }else{
    stop("The true H function replies only on 5 exposures. Indice must belong to c(1:2, 4:6)")
  }
  
  
  return(sfunc(u))
  
}

####### For h3&h2
Href <- function(z, index = c(1:2, 4:6), idx = 1, h = 3){
  
  if(idx ==1){
    u = 1/6*(z[idx] + Zquant[index[2]] + 2*Zquant[index[3]] + 2*Zquant[index[4]] + 3*Zquant[index[5]] +
               1/2*z[idx]*Zquant[index[2]] + 1/3*Zquant[index[3]]*Zquant[index[4]]*Zquant[index[5]])   # H2 per se.
    
  }else if(idx ==2){
    u = 1/6*(Zquant[index[1]] + z[idx] + 2*Zquant[index[3]] + 2*Zquant[index[4]] + 3*Zquant[index[5]] +
               1/2*Zquant[index[1]]*z[idx] + 1/3*Zquant[index[3]]*Zquant[index[4]]*Zquant[index[5]])
    
  }else if(idx ==4){
    u = 1/6*(Zquant[index[1]] + Zquant[index[2]] + 2*z[idx] + 2*Zquant[index[4]] + 3*Zquant[index[5]] +
               1/2*Zquant[index[1]]*Zquant[index[2]] + 1/3*z[idx]*Zquant[index[4]]*Zquant[index[5]])
    
  }else if(idx ==5){
    u = 1/6*(Zquant[index[1]] + Zquant[index[2]] + 2*Zquant[index[3]] + 2*z[idx] + 3*Zquant[index[5]] +
               1/2*Zquant[index[1]]*Zquant[index[2]] + 1/3*Zquant[index[3]]*z[idx]*Zquant[index[5]])
    
  }else if(idx ==6){
    u = 1/6*(Zquant[index[1]] + Zquant[index[2]] + 2*Zquant[index[3]] + 2*Zquant[index[4]] + 3*z[idx] +
               1/2*Zquant[index[1]]*Zquant[index[2]] + 1/3*Zquant[index[3]]*Zquant[index[4]]*z[idx])
    
  }else{
    stop("The true H function replies only on 5 exposures. Indice must belong to c(1:2, 4:6)")
  }
  
  if(h == 2){
    return(u)
  }
  
  return(sfunc(u))
  
}

##### For H2
### reference line at median
Zquant <- apply(Ztest, 2, \(x) quantile(x, probs = .5))
Href1 <- apply(Ztest, 1, \(x) Href(x, idx = 1, h = 2))
Href2 <- apply(Ztest, 1, \(x) Href(x, idx = 2, h = 2))
Href4 <- apply(Ztest, 1, \(x) Href(x, idx = 4, h = 2))
Href5 <- apply(Ztest, 1, \(x) Href(x, idx = 5, h = 2))
Href6 <- apply(Ztest, 1, \(x) Href(x, idx = 6, h = 2))
### reference line at quant 25.
Zquant <- apply(Ztest, 2, \(x) quantile(x, probs = .25))
Href1_25 <- apply(Ztest, 1, \(x) Href(x, idx = 1, h = 2))
Href2_25 <- apply(Ztest, 1, \(x) Href(x, idx = 2, h = 2))
Href4_25 <- apply(Ztest, 1, \(x) Href(x, idx = 4, h = 2))
Href5_25 <- apply(Ztest, 1, \(x) Href(x, idx = 5, h = 2))
Href6_25 <- apply(Ztest, 1, \(x) Href(x, idx = 6, h = 2))
### reference line at quant 75
Zquant <- apply(Ztest, 2, \(x) quantile(x, probs = .75))
Href1_75 <- apply(Ztest, 1, \(x) Href(x, idx = 1, h = 2))
Href2_75 <- apply(Ztest, 1, \(x) Href(x, idx = 2, h = 2))
Href4_75 <- apply(Ztest, 1, \(x) Href(x, idx = 4, h = 2))
Href5_75 <- apply(Ztest, 1, \(x) Href(x, idx = 5, h = 2))
Href6_75 <- apply(Ztest, 1, \(x) Href(x, idx = 6, h = 2))


################################################################################
# For h3, h1
### reference line at median
Zquant <- apply(Ztest, 2, \(x) quantile(x, probs = .5))
Href1 <- apply(Ztest, 1, \(x) Href(x, idx = 1))
Href2 <- apply(Ztest, 1, \(x) Href(x, idx = 2))
Href4 <- apply(Ztest, 1, \(x) Href(x, idx = 4))
Href5 <- apply(Ztest, 1, \(x) Href(x, idx = 5))
Href6 <- apply(Ztest, 1, \(x) Href(x, idx = 6))
### reference line at quant 25.
Zquant <- apply(Ztest, 2, \(x) quantile(x, probs = .25))
Href1_25 <- apply(Ztest, 1, \(x) Href(x, idx = 1))
Href2_25 <- apply(Ztest, 1, \(x) Href(x, idx = 2))
Href4_25 <- apply(Ztest, 1, \(x) Href(x, idx = 4))
Href5_25 <- apply(Ztest, 1, \(x) Href(x, idx = 5))
Href6_25 <- apply(Ztest, 1, \(x) Href(x, idx = 6))
### reference line at quant 75
Zquant <- apply(Ztest, 2, \(x) quantile(x, probs = .75))
Href1_75 <- apply(Ztest, 1, \(x) Href(x, idx = 1))
Href2_75 <- apply(Ztest, 1, \(x) Href(x, idx = 2))
Href4_75 <- apply(Ztest, 1, \(x) Href(x, idx = 4))
Href5_75 <- apply(Ztest, 1, \(x) Href(x, idx = 5))
Href6_75 <- apply(Ztest, 1, \(x) Href(x, idx = 6))


# plots for relevant exposures z1, z2, z4, z5, z6.

df_z1 <- data.frame(z1 = Ztest[,"z1"], BKMR = gam.pred.smt.bkmr[,"s(z1)"], 
                    SoftBART = gam.pred.smt.gsbart_ntr50[,"s(z1)"],
                    Quant50 = Href1, Quant25 = Href1_25, Quant75 = Href1_75) %>% 
  pivot_longer(cols = c(BKMR, SoftBART, Quant50, Quant25, Quant75), names_to = "Label", values_to = "smooth") %>% 
  mutate(Type = ifelse(Label %in% c("BKMR", "SoftBART"), "Model", "Quantile"))


df_z2 <- data.frame(z2 = Ztest[,"z2"], BKMR = gam.pred.smt.bkmr[,"s(z2)"], 
                    SoftBART = gam.pred.smt.gsbart_ntr50[,"s(z2)"],
                    Quant50 = Href2, Quant25 = Href2_25, Quant75 = Href2_75) %>% 
  pivot_longer(cols = c(BKMR, SoftBART, Quant50, Quant25, Quant75), names_to = "Label", values_to = "smooth") %>% 
  mutate(Type = ifelse(Label %in% c("BKMR", "SoftBART"), "Model", "Quantile"))

df_z4 <- data.frame(z4 = Ztest[,"z4"], BKMR = gam.pred.smt.bkmr[,"s(z4)"], 
                    SoftBART = gam.pred.smt.gsbart_ntr50[,"s(z4)"],
                    Quant50 = Href4, Quant25 = Href4_25, Quant75 = Href4_75) %>% 
  pivot_longer(cols = c(BKMR, SoftBART, Quant50, Quant25, Quant75), names_to = "Label", values_to = "smooth") %>% 
  mutate(Type = ifelse(Label %in% c("BKMR", "SoftBART"), "Model", "Quantile"))

df_z5 <- data.frame(z5 = Ztest[,"z5"], BKMR = gam.pred.smt.bkmr[,"s(z5)"], 
                    SoftBART = gam.pred.smt.gsbart_ntr50[,"s(z5)"],
                    Quant50 = Href5, Quant25 = Href5_25, Quant75 = Href5_75) %>% 
  pivot_longer(cols = c(BKMR, SoftBART, Quant50, Quant25, Quant75), names_to = "Label", values_to = "smooth") %>% 
  mutate(Type = ifelse(Label %in% c("BKMR", "SoftBART"), "Model", "Quantile"))

df_z6 <- data.frame(z6 = Ztest[,"z6"], BKMR = gam.pred.smt.bkmr[,"s(z6)"], 
                    SoftBART = gam.pred.smt.gsbart_ntr50[,"s(z6)"],
                    Quant50 = Href6, Quant25 = Href6_25, Quant75 = Href6_75) %>% 
  pivot_longer(cols = c(BKMR, SoftBART, Quant50, Quant25, Quant75), names_to = "Label", values_to = "smooth") %>% 
  mutate(Type = ifelse(Label %in% c("BKMR", "SoftBART"), "Model", "Quantile"))



common_theme <- theme(
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 12),
  legend.title = element_blank(),
  legend.text = element_text(size = 10),
  legend.position = "none",
  plot.title = element_blank()
)

# for h1
label_map <- c("SoftBART" = '"GAM Approximation for modified BART"', 
               "BKMR" = '"GAM Approximation for BKMR"', 
               "Quant25" = "True~h[1]~'with Remaining at 1st Quartile'",
               "Quant50" = "True~h[1]~'with Remaining at Median'",
               "Quant75" = "True~h[1]~'with Remaining at 3rd Quartile'"
)

# for h2
label_map <- c("SoftBART" = '"GAM Approximation for modified BART"', 
               "BKMR" = '"GAM Approximation for BKMR"', 
               "Quant25" = "True~h[2]~'with Remaining at 1st Quartile'",
               "Quant50" = "True~h[2]~'with Remaining at Median'",
               "Quant75" = "True~h[2]~'with Remaining at 3rd Quartile'"
)


# for h3
label_map <- c("SoftBART" = '"GAM Approximation for modified BART"', 
               "BKMR" = '"GAM Approximation for BKMR"', 
               "Quant25" = "True~h[3]~'with Remaining at 1st Quartile'",
               "Quant50" = "True~h[3]~'with Remaining at Median'",
               "Quant75" = "True~h[3]~'with Remaining at 3rd Quartile'"
)


# prepare legend
colors_def <- scales::hue_pal()(length(names(label_map)))
color_map <- c("SoftBART" = colors_def[length(colors_def)], 
               "BKMR" = colors_def[1], 
               "Quant25" = colors_def[2],
               "Quant50" = colors_def[3],
               "Quant75" = colors_def[4]
)

prepare_legend <- function(df, label_map){
  model_level <- sort(unique(df$Label[df$Type=="Model"]), decreasing = TRUE)
  quant_level <- sort(unique(df$Label[df$Type=="Quantile"]))
  custom_order <- c(model_level, quant_level)
  
  df$Label <- factor(df$Label, levels = custom_order)
  display_labels <- label_map[levels(df$Label)]
  display_labels <- lapply(display_labels, function(lbl) parse(text = lbl)[[1]])
  
  list(data = df, labels = display_labels)
}


out <- prepare_legend(df_z1, label_map)
p1 <- ggplot(out$data, aes(z1, smooth, color = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = expression(Z[1])) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(df_z2, label_map)
p2 <- ggplot(out$data, aes(z2, smooth, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = expression(Z[2])) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(df_z4, label_map)
p3 <- ggplot(out$data, aes(z4, smooth, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = expression(Z[4])) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(df_z5, label_map)
p4 <- ggplot(out$data, aes(z5, smooth, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = expression(Z[5])) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)



# Use long description legend labels
common_theme <- theme(
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 12),
  legend.title = element_blank(),
  legend.text = element_text(size = 10),
  legend.position = "right",
  plot.title = element_blank()
)

out_legend <- prepare_legend(df_z1, label_map)
legend_only <- as_ggplot(get_legend(
  ggplot(out_legend$data, aes(z1, smooth, colour = Label)) +
    geom_line(data = out_legend$data %>% filter(Type == "Model"), linewidth = 1.2) + 
    geom_point(data = out_legend$data %>% filter(Type == "Model"), size = 3) + 
    geom_line(data = out_legend$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
    labs(y = "Marginal Effects") + common_theme + 
    scale_color_manual(values = color_map, labels = out$labels)
  
))



final_plot2 <- (p1 + p2) / (p3 + legend_only) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 16, face = "bold")
    ))
final_plot2 + plot_layout(guides = "collect")







