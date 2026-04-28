# Download the NHANES data from:
# https://github.com/lizzyagibson/SHARP.Mixtures.Workshop/tree/master/Supervised/BKMR
#
# Then follow the data manipulation steps below, copied from Gibson et al.,
# to create the `.RData` file used in this project.

# nhanes <- na.omit(read.csv("studypop.csv"))
# 
# ## center/scale continous covariates and create indicators for categorical covariates
# nhanes$age_z         <- scale(nhanes$age_cent)         ## center and scale age
# nhanes$agez_sq       <- nhanes$age_z^2                 ## square this age variable
# nhanes$bmicat2       <- as.numeric(nhanes$bmi_cat3==2) ## 25 <= BMI < 30
# nhanes$bmicat3       <- as.numeric(nhanes$bmi_cat3==3) ## BMI >= 30 (BMI < 25 is the reference)
# nhanes$educat1       <- as.numeric(nhanes$edu_cat==1)  ## no high school diploma
# nhanes$educat3       <- as.numeric(nhanes$edu_cat==3)  ## some college or AA degree
# nhanes$educat4       <- as.numeric(nhanes$edu_cat==4)  ## college grad or above (reference is high schol grad/GED or equivalent)
# nhanes$otherhispanic <- as.numeric(nhanes$race_cat==1) ## other Hispanic or other race - including multi-racial
# nhanes$mexamerican   <- as.numeric(nhanes$race_cat==2) ## Mexican American
# nhanes$black         <- as.numeric(nhanes$race_cat==3) ## non-Hispanic Black (non-Hispanic White as reference group)
# nhanes$wbcc_z        <- scale(nhanes$LBXWBCSI)
# nhanes$lymphocytes_z <- scale(nhanes$LBXLYPCT)
# nhanes$monocytes_z   <- scale(nhanes$LBXMOPCT)
# nhanes$neutrophils_z <- scale(nhanes$LBXNEPCT)
# nhanes$eosinophils_z <- scale(nhanes$LBXEOPCT)
# nhanes$basophils_z   <- scale(nhanes$LBXBAPCT)
# nhanes$lncotinine_z  <- scale(nhanes$ln_lbxcot)         ## to access smoking status, scaled ln cotinine levels
# 
# 
# ## our y variable - ln transformed and scaled mean telomere length
# lnLTL_z <- scale(log(nhanes$TELOMEAN))
# 
# ## our Z matrix
# mixture <- with(nhanes, cbind(LBX074LA, LBX099LA, LBX118LA, LBX138LA, LBX153LA, LBX170LA, LBX180LA, LBX187LA,
#                               LBX194LA, LBXHXCLA, LBXPCBLA,
#                               LBXD03LA, LBXD05LA, LBXD07LA,
#                               LBXF03LA, LBXF04LA, LBXF05LA, LBXF08LA))
# lnmixture   <- apply(mixture, 2, log)
# lnmixture_z <- scale(lnmixture)
# colnames(lnmixture_z) <- c(paste0("PCB",c(74, 99, 118, 138, 153, 170, 180, 187, 194, 169, 126)),
#                            paste0("Dioxin",1:3), paste0("Furan",1:4))
# 
# ## our X matrix
# covariates <- with(nhanes, cbind(age_z, agez_sq, male, bmicat2, bmicat3, educat1, educat3, educat4,
#                                  otherhispanic, mexamerican, black, wbcc_z, lymphocytes_z, monocytes_z,
#                                  neutrophils_z, eosinophils_z, basophils_z, lncotinine_z))
# 
# ### create knots matrix for Gaussian predictive process (to speed up BKMR with large datasets)
# set.seed(10)
# knots100     <- fields::cover.design(lnmixture_z, nd = 100)$design
# save(list = ls(), file="BKMR_NHANES_knots100.RData")


load("BKMR_NHANES_knots100.RData")


### fit general softBART regression
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


# for grouped regression.
## run BKMR.
start_bkmr = now()
fit_bkmr_knots100_grouped <-  kmbayes(y=lnLTL_z, Z=lnmixture_z, X=covariates, iter=100000, verbose=TRUE, varsel=TRUE,
                                      groups=c(rep(1,times=2), 2, rep(1,times=6), rep(3,times=2),rep(2,times=7)),
                                      knots=knots100)
time_bkmr_grouped <- as.duration(now() - start_bkmr)


## run gsBART
start_time <- now()
fit_gsbart_grouped <- gsbart_regression(X = as.matrix(covariates), Y = lnLTL_z, Z = as.matrix(lnmixture_z), num_tree = 20,
                                              num_burn = 50000, num_save = 1000, num_thin = 50, alpha_val = 1, alpha_comp = 1,
                                              group = c(rep(1,times=2), 2, rep(1,times=6), rep(3,times=2),rep(2,times=7)))
time_gsbart_grouped <- as.duration(now() - start_time)

## Calculate PIPs - grouped
PIPs_gsbart_grouped <- calPIPs_gsbart(fit = fit_gsbart_grouped, group = TRUE,
                                            vars = lnmixture_z)

PIPs_bkmr_grouped <- ExtractPIPs(fit_bkmr_knots100_grouped)



################################################################################
### plot the PIPs
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)

colnames(PIPs_bkmr_grouped) <- c("variable", "group", "groupPIP_bkmr", "condPIP_bkmr")
colnames(PIPs_gsbart_grouped) <- c("variable", "group", "groupPIP_gsbart", "condPIP_gsbart")

PIPs_merged <- merge(PIPs_bkmr_grouped, PIPs_gsbart_grouped, by = c("variable", "group"))
colnames(PIPs_merged)


PIPs_merged <- PIPs_merged %>% 
  arrange(group) %>%
  mutate(variable = case_when(
    variable == "Dioxin1" ~ "1,2,3,6,7,8-hxcdd" ,
    variable == "Dioxin2" ~ "1,2,3,4,6,7,8-hpcdd",
    variable == "Dioxin3" ~ "1,2,3,4,6,7,8,9-ocdd",
    variable == "Furan1" ~ "2,3,4,7,8-pncdf",
    variable == "Furan2" ~ "1,2,3,4,7,8-hxcdf",
    variable == "Furan3" ~ "1,2,3,6,7,8-hxcdf",
    variable == "Furan4" ~ "1,2,3,4,6,7,8-hxcdf",
    TRUE ~ variable
    
  ))

PIPs_merged <- PIPs_merged %>% 
  mutate(groupPIP_gsbart = as.numeric(groupPIP_gsbart),
         condPIP_gsbart = as.numeric(condPIP_gsbart))

PIPs_merged <- PIPs_merged %>% 
  mutate(group = paste0("Group", group))

y_axis_order <- c()
for (g in unique(PIPs_merged$group)) {
  y_axis_order <- c(y_axis_order, g, PIPs_merged$variable[PIPs_merged$group==g])
}


df_group <- PIPs_merged %>% 
  select(group, starts_with("groupPIP_")) %>% 
  pivot_longer(cols = starts_with("groupPIP_"),
               names_to = "PIP_Type",
               values_to = "PIP_Value") %>% 
  distinct() %>% 
  mutate(y_axis = group,
         Model = ifelse(grepl("bkmr", PIP_Type), "BKMR", "Soft BART"),
         type = "Group")


df_vars <- PIPs_merged %>% 
  select(variable, starts_with("condPIP_")) %>% 
  pivot_longer(cols = starts_with("condPIP_"),
               names_to = "PIP_Type",
               values_to = "PIP_Value") %>% 
  mutate(y_axis = variable,
         Model = ifelse(grepl("bkmr", PIP_Type), "BKMR", "Soft BART"),
         type = "Conditional")

# combine both datasets for plotting
df_long <- bind_rows(df_group, df_vars) %>% 
  mutate(y_axis = factor(y_axis, levels = rev(y_axis_order)))

bold_labels <- c("Group1", "Group2", "Group3")

ggplot(df_long, aes(x = PIP_Value, y = y_axis, color = type, shape = Model)) + 
  geom_point(size = 3) + 
  theme_minimal() + 
  labs(title = "Figure 3: Group and Conditional PIPs for \n Soft BART and BKMR Models",
       x = "PIP Value", y = "Conditional / Group", 
       color = "Type",
       shape = "Model") + 
  scale_y_discrete(labels = function(y) {
    ifelse(y %in% bold_labels, paste0("**", y, "**"), y)
  }) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_markdown(size = 12, color = "black", face = "plain"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12))


################################################################################
# GAM results for NHANES data results.
library(bkmr)
library(mgcv)
library(dplyr)
library(ggplot2)
library(grid)
library(tidyr)
library(patchwork)
library(cowplot)
library(ggpubr)



sel<-seq(50001,100000,by=50)

## with group BKMR GAM results
htest_bkmr_grp <- SamplePred(fit_bkmr_knots100_grouped, Znew = lnmixture_z, Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                             sel = sel)
htest_bkmr_grp_mean <- as.matrix(colMeans(htest_bkmr_grp))
htest_bkmr_grp_sampl <- t(as.matrix(htest_bkmr_grp))

GAMfit_bkmr_grp <- gam(htest_bkmr_grp_mean ~
                         s(PCB74) + s(PCB99) + s(PCB118) +
                         s(PCB138) + s(PCB153) + s(PCB170) +
                         s(PCB180) + s(PCB187) + s(PCB194) +
                         s(PCB169)+s(PCB126) + s(Dioxin1)+
                         s(Dioxin2)+s(Dioxin3)+s(Furan1)+
                         s(Furan2)+s(Furan3)+s(Furan4),
                       data = as.data.frame(lnmixture_z))

## with group gsBART GAM results
fit_gsbart <- fit_gsbart_grouped; htest_bkmr_mean <- htest_bkmr_grp_mean; htest_bkmr_sampl <- htest_bkmr_grp_sampl
htest_gsbart <- fit_gsbart$h_test
htest_gsbart_mean <- as.matrix(colMeans(htest_gsbart))
htest_gsbart_sampl <- t(as.matrix(htest_gsbart))

GAMfit_gsbart_grp <- gam(htest_gsbart_mean ~
                       s(PCB74) + s(PCB99) + s(PCB118) +
                       s(PCB138) + s(PCB153) + s(PCB170) +
                       s(PCB180) + s(PCB187) + s(PCB194) +
                       s(PCB169)+s(PCB126) + s(Dioxin1)+
                       s(Dioxin2)+s(Dioxin3)+s(Furan1)+
                       s(Furan2)+s(Furan3)+s(Furan4),
                     data = as.data.frame(lnmixture_z))

# predicted smooth items
gam.pred.bkmr <- predict(GAMfit_bkmr_grp, newdata = as.data.frame(lnmixture_z), type = "terms",
                         se.fit = TRUE)
gam.pred.gsbart <- predict(GAMfit_gsbart_grp, newdata = as.data.frame(lnmixture_z), type = "terms",
                           se.fit = TRUE)


gam.pred.smt.bkmr <- gam.pred.bkmr$fit
gam.pred.smt.gsbart <- gam.pred.gsbart$fit

fit_bkmr_knots100 <- fit_bkmr_knots100_grouped

# reference line plotting data
RLplotdata <- function(data, q.fixed = c(.25, .5, .75), keep_col){
  
  if(is.matrix(data)) data <- as.data.frame(data)
  
  data_fixed <- list()
  for (i in 1:length(q.fixed)) {
    data_fixed[[as.character(q.fixed[i])]] <- data
    for(col in colnames(data)){
      if(col != keep_col){
        data_fixed[[as.character(q.fixed[i])]][[col]] <- quantile(data[[col]], probs = q.fixed[i])
      }
    }
  }
  
  
  
  return(data_fixed)
}

newdata = as.data.frame(lnmixture_z)

##
data_fixed.PCB74 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB74")
hnew.bkmr.PCB74_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB74$`0.25`,
                                           Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                           sel = sel))
hnew.bkmr.PCB74_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB74$`0.5`,
                                           Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                           sel = sel))
hnew.bkmr.PCB74_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB74$`0.75`,
                                           Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                           sel = sel))

pred.bkmr.PCB74 <- data.frame(PCB74 = newdata[["PCB74"]], 
                              Pred = hnew.bkmr.PCB74_q50
                              ) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB74_q25 <- data.frame(PCB74 = newdata[["PCB74"]],
                                  Pred = hnew.bkmr.PCB74_q25
                                  ) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB74_q75 <- data.frame(PCB74 = newdata[["PCB74"]],
                                  Pred = hnew.bkmr.PCB74_q75
                                  ) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB74 <- bind_rows(pred.bkmr.PCB74, pred.bkmr.PCB74_q25, pred.bkmr.PCB74_q75)


##
data_fixed.PCB170 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB170")
hnew.bkmr.PCB170_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB170$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB170_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB170$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB170_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB170$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.PCB170 <- data.frame(PCB170 =  newdata[["PCB170"]],
                               Pred = hnew.bkmr.PCB170_q50
                               ) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB170_q25 <- data.frame(PCB170 = newdata[["PCB170"]],
                                   Pred = hnew.bkmr.PCB170_q25
                                   ) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB170_q75 <- data.frame(PCB170 = newdata[["PCB170"]],
                                   Pred = hnew.bkmr.PCB170_q75
                                   ) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB170 <- bind_rows(pred.bkmr.PCB170, pred.bkmr.PCB170_q25, pred.bkmr.PCB170_q75)


##
data_fixed.PCB126 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB126")

hnew.bkmr.PCB126_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB126$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB126_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB126$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB126_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB126$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.PCB126 <- data.frame(PCB126 = newdata[["PCB126"]],
                               Pred = hnew.bkmr.PCB126_q50
                               ) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB126_q25 <- data.frame(PCB126 = newdata[["PCB126"]],
                                   Pred = hnew.bkmr.PCB126_q25
                                   ) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB126_q75 <- data.frame(PCB126 = newdata[["PCB126"]],
                                   Pred = hnew.bkmr.PCB126_q75
                                   ) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB126 <- bind_rows(pred.bkmr.PCB126, pred.bkmr.PCB126_q25, pred.bkmr.PCB126_q75)


##
data_fixed.Furan2 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "Furan2")
hnew.bkmr.Furan2_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan2$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.Furan2_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan2$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.Furan2_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan2$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.Furan2 <- data.frame(Furan2 = newdata[["Furan2"]],
                               Pred = hnew.bkmr.Furan2_q50
                               ) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.Furan2_q25 <- data.frame(Furan2 = newdata[["Furan2"]],
                                   Pred = hnew.bkmr.Furan2_q25
                                   ) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.Furan2_q75 <- data.frame(Furan2 = newdata[["Furan2"]],
                                   Pred = hnew.bkmr.Furan2_q75
                                   ) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_Furan2 <- bind_rows(pred.bkmr.Furan2, pred.bkmr.Furan2_q25, pred.bkmr.Furan2_q75)


##
data_fixed.PCB99 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB99")

hnew.bkmr.PCB99_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB99$`0.25`,
                                           Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                           sel = sel))
hnew.bkmr.PCB99_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB99$`0.5`,
                                           Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                           sel = sel))
hnew.bkmr.PCB99_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB99$`0.75`,
                                           Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                           sel = sel))

pred.bkmr.PCB99 <- data.frame(PCB99 = newdata[["PCB99"]],
                              Pred = hnew.bkmr.PCB99_q50
                              ) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB99_q25 <- data.frame(PCB99 = newdata[["PCB99"]],
                                  Pred = hnew.bkmr.PCB99_q25
                                  ) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB99_q75 <- data.frame(PCB99 = newdata[["PCB99"]],
                                  Pred = hnew.bkmr.PCB99_q75
                                  ) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB99 <- bind_rows(pred.bkmr.PCB99, pred.bkmr.PCB99_q25, pred.bkmr.PCB99_q75)


##
data_fixed.PCB118 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB118")

hnew.bkmr.PCB118_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB118$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB118_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB118$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB118_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB118$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.PCB118 <- data.frame(PCB118 = newdata[["PCB118"]],
                               Pred = hnew.bkmr.PCB118_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB118_q25 <- data.frame(PCB118 = newdata[["PCB118"]],
                                   Pred = hnew.bkmr.PCB118_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB118_q75 <- data.frame(PCB118 = newdata[["PCB118"]],
                                   Pred = hnew.bkmr.PCB118_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB118 <- bind_rows(pred.bkmr.PCB118, pred.bkmr.PCB118_q25, pred.bkmr.PCB118_q75)


##
data_fixed.PCB138 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB138")

hnew.bkmr.PCB138_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB138$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB138_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB138$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB138_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB138$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.PCB138 <- data.frame(PCB138 = newdata[["PCB138"]],
                               Pred = hnew.bkmr.PCB138_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB138_q25 <- data.frame(PCB138 = newdata[["PCB138"]],
                                   Pred = hnew.bkmr.PCB138_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB138_q75 <- data.frame(PCB138 = newdata[["PCB138"]],
                                   Pred = hnew.bkmr.PCB138_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB138 <- bind_rows(pred.bkmr.PCB138, pred.bkmr.PCB138_q25, pred.bkmr.PCB138_q75)


##
data_fixed.PCB153 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB153")

hnew.bkmr.PCB153_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB153$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB153_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB153$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB153_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB153$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.PCB153 <- data.frame(PCB153 = newdata[["PCB153"]],
                               Pred = hnew.bkmr.PCB153_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB153_q25 <- data.frame(PCB153 = newdata[["PCB153"]],
                                   Pred = hnew.bkmr.PCB153_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB153_q75 <- data.frame(PCB153 = newdata[["PCB153"]],
                                   Pred = hnew.bkmr.PCB153_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB153 <- bind_rows(pred.bkmr.PCB153, pred.bkmr.PCB153_q25, pred.bkmr.PCB153_q75)


##
data_fixed.PCB180 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB180")

hnew.bkmr.PCB180_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB180$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB180_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB180$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB180_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB180$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.PCB180 <- data.frame(PCB180 = newdata[["PCB180"]],
                               Pred = hnew.bkmr.PCB180_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB180_q25 <- data.frame(PCB180 = newdata[["PCB180"]],
                                   Pred = hnew.bkmr.PCB180_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB180_q75 <- data.frame(PCB180 = newdata[["PCB180"]],
                                   Pred = hnew.bkmr.PCB180_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB180 <- bind_rows(pred.bkmr.PCB180, pred.bkmr.PCB180_q25, pred.bkmr.PCB180_q75)


##
data_fixed.PCB187 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB187")
hnew.bkmr.PCB187_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB187$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB187_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB187$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB187_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB187$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.PCB187 <- data.frame(PCB187 = newdata[["PCB187"]],
                               Pred = hnew.bkmr.PCB187_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB187_q25 <- data.frame(PCB187 = newdata[["PCB187"]],
                                   Pred = hnew.bkmr.PCB187_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB187_q75 <- data.frame(PCB187 = newdata[["PCB187"]],
                                   Pred = hnew.bkmr.PCB187_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB187 <- bind_rows(pred.bkmr.PCB187, pred.bkmr.PCB187_q25, pred.bkmr.PCB187_q75)


##
data_fixed.PCB194 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB194")
hnew.bkmr.PCB194_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB194$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB194_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB194$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB194_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB194$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.PCB194 <- data.frame(PCB194 = newdata[["PCB194"]],
                               Pred = hnew.bkmr.PCB194_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB194_q25 <- data.frame(PCB194 = newdata[["PCB194"]],
                                   Pred = hnew.bkmr.PCB194_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB194_q75 <- data.frame(PCB194 = newdata[["PCB194"]],
                                   Pred = hnew.bkmr.PCB194_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB194 <- bind_rows(pred.bkmr.PCB194, pred.bkmr.PCB194_q25, pred.bkmr.PCB194_q75)


##
data_fixed.PCB169 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "PCB169")
hnew.bkmr.PCB169_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB169$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB169_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB169$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.PCB169_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.PCB169$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.PCB169 <- data.frame(PCB169 = newdata[["PCB169"]],
                               Pred = hnew.bkmr.PCB169_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.PCB169_q25 <- data.frame(PCB169 = newdata[["PCB169"]],
                                   Pred = hnew.bkmr.PCB169_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.PCB169_q75 <- data.frame(PCB169 = newdata[["PCB169"]],
                                   Pred = hnew.bkmr.PCB169_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_PCB169 <- bind_rows(pred.bkmr.PCB169, pred.bkmr.PCB169_q25, pred.bkmr.PCB169_q75)


##
data_fixed.Dioxin1 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "Dioxin1")
hnew.bkmr.Dioxin1_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Dioxin1$`0.25`,
                                             Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                             sel = sel))
hnew.bkmr.Dioxin1_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Dioxin1$`0.5`,
                                             Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                             sel = sel))
hnew.bkmr.Dioxin1_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Dioxin1$`0.75`,
                                             Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                             sel = sel))

pred.bkmr.Dioxin1 <- data.frame(Dioxin1 = newdata[["Dioxin1"]],
                                Pred = hnew.bkmr.Dioxin1_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.Dioxin1_q25 <- data.frame(Dioxin1 = newdata[["Dioxin1"]],
                                    Pred = hnew.bkmr.Dioxin1_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.Dioxin1_q75 <- data.frame(Dioxin1 = newdata[["Dioxin1"]],
                                    Pred = hnew.bkmr.Dioxin1_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_Dioxin1 <- bind_rows(pred.bkmr.Dioxin1, pred.bkmr.Dioxin1_q25, pred.bkmr.Dioxin1_q75)


##
data_fixed.Dioxin2 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "Dioxin2")
hnew.bkmr.Dioxin2_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Dioxin2$`0.25`,
                                             Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                             sel = sel))
hnew.bkmr.Dioxin2_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Dioxin2$`0.5`,
                                             Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                             sel = sel))
hnew.bkmr.Dioxin2_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Dioxin2$`0.75`,
                                             Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                             sel = sel))

pred.bkmr.Dioxin2 <- data.frame(Dioxin2 = newdata[["Dioxin2"]],
                                Pred = hnew.bkmr.Dioxin2_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.Dioxin2_q25 <- data.frame(Dioxin2 = newdata[["Dioxin2"]],
                                    Pred = hnew.bkmr.Dioxin2_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.Dioxin2_q75 <- data.frame(Dioxin2 = newdata[["Dioxin2"]],
                                    Pred = hnew.bkmr.Dioxin2_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_Dioxin2 <- bind_rows(pred.bkmr.Dioxin2, pred.bkmr.Dioxin2_q25, pred.bkmr.Dioxin2_q75)


##
data_fixed.Dioxin3 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "Dioxin3")
hnew.bkmr.Dioxin3_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Dioxin3$`0.25`,
                                             Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                             sel = sel))
hnew.bkmr.Dioxin3_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Dioxin3$`0.5`,
                                             Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                             sel = sel))
hnew.bkmr.Dioxin3_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Dioxin3$`0.75`,
                                             Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                             sel = sel))

pred.bkmr.Dioxin3 <- data.frame(Dioxin3 = newdata[["Dioxin3"]],
                                Pred = hnew.bkmr.Dioxin3_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.Dioxin3_q25 <- data.frame(Dioxin3 = newdata[["Dioxin3"]],
                                    Pred = hnew.bkmr.Dioxin3_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.Dioxin3_q75 <- data.frame(Dioxin3 = newdata[["Dioxin3"]],
                                    Pred = hnew.bkmr.Dioxin3_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_Dioxin3 <- bind_rows(pred.bkmr.Dioxin3, pred.bkmr.Dioxin3_q25, pred.bkmr.Dioxin3_q75)


##
data_fixed.Furan1 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "Furan1")
hnew.bkmr.Furan1_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan1$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.Furan1_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan1$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.Furan1_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan1$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.Furan1 <- data.frame(Furan1 = newdata[["Furan1"]],
                               Pred = hnew.bkmr.Furan1_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.Furan1_q25 <- data.frame(Furan1 = newdata[["Furan1"]],
                                   Pred = hnew.bkmr.Furan1_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.Furan1_q75 <- data.frame(Furan1 = newdata[["Furan1"]],
                                   Pred = hnew.bkmr.Furan1_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_Furan1 <- bind_rows(pred.bkmr.Furan1, pred.bkmr.Furan1_q25, pred.bkmr.Furan1_q75)


##
data_fixed.Furan3 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "Furan3")
hnew.bkmr.Furan3_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan3$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.Furan3_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan3$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.Furan3_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan3$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.Furan3 <- data.frame(Furan3 = newdata[["Furan3"]],
                               Pred = hnew.bkmr.Furan3_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.Furan3_q25 <- data.frame(Furan3 = newdata[["Furan3"]],
                                   Pred = hnew.bkmr.Furan3_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.Furan3_q75 <- data.frame(Furan3 = newdata[["Furan3"]],
                                   Pred = hnew.bkmr.Furan3_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_Furan3 <- bind_rows(pred.bkmr.Furan3, pred.bkmr.Furan3_q25, pred.bkmr.Furan3_q75)


##
data_fixed.Furan4 <- RLplotdata(newdata, q.fixed = c(.25, .5, .75), keep_col = "Furan4")
hnew.bkmr.Furan4_q25 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan4$`0.25`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.Furan4_q50 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan4$`0.5`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))
hnew.bkmr.Furan4_q75 <- colMeans(SamplePred(fit_bkmr_knots100, Znew = data_fixed.Furan4$`0.75`,
                                            Xnew = matrix(rep(0, nrow(covariates)*ncol(covariates)), ncol = ncol(covariates)),
                                            sel = sel))

pred.bkmr.Furan4 <- data.frame(Furan4 = newdata[["Furan4"]],
                               Pred = hnew.bkmr.Furan4_q50
) %>% 
  mutate(Type = "Quantile", Label = "50th")

pred.bkmr.Furan4_q25 <- data.frame(Furan4 = newdata[["Furan4"]],
                                   Pred = hnew.bkmr.Furan4_q25
) %>% 
  mutate(Type = "Quantile", Label = "25th")

pred.bkmr.Furan4_q75 <- data.frame(Furan4 = newdata[["Furan4"]],
                                   Pred = hnew.bkmr.Furan4_q75
) %>% 
  mutate(Type = "Quantile", Label = "75th")

combined_df_Furan4 <- bind_rows(pred.bkmr.Furan4, pred.bkmr.Furan4_q25, pred.bkmr.Furan4_q75)





#############################
df_PCB74 <- data.frame(PCB74 = newdata[["PCB74"]], BKMR = gam.pred.smt.bkmr[,"s(PCB74)"], 
                       SoftBART = gam.pred.smt.gsbart[,"s(PCB74)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB74 <- bind_rows(df_PCB74, combined_df_PCB74)
combined_df_PCB74$Label <- as.factor(combined_df_PCB74$Label)


df_PCB170 <- data.frame(PCB170 = newdata[["PCB170"]], BKMR = gam.pred.smt.bkmr[,"s(PCB170)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(PCB170)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB170 <- bind_rows(df_PCB170, combined_df_PCB170)
combined_df_PCB170$Label <- as.factor(combined_df_PCB170$Label)


df_PCB126 <- data.frame(PCB126 = newdata[["PCB126"]], BKMR = gam.pred.smt.bkmr[,"s(PCB126)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(PCB126)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB126 <- bind_rows(df_PCB126, combined_df_PCB126)
combined_df_PCB126$Label <- as.factor(combined_df_PCB126$Label)


df_Furan2 <- data.frame(Furan2 = newdata[["Furan2"]], BKMR = gam.pred.smt.bkmr[,"s(Furan2)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(Furan2)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_Furan2 <- bind_rows(df_Furan2, combined_df_Furan2)
combined_df_Furan2$Label <- as.factor(combined_df_Furan2$Label)



df_PCB99 <- data.frame(PCB99 = newdata[["PCB99"]], BKMR = gam.pred.smt.bkmr[,"s(PCB99)"], 
                       SoftBART = gam.pred.smt.gsbart[,"s(PCB99)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB99 <- bind_rows(df_PCB99, combined_df_PCB99)
combined_df_PCB99$Label <- as.factor(combined_df_PCB99$Label)


df_PCB118 <- data.frame(PCB118 = newdata[["PCB118"]], BKMR = gam.pred.smt.bkmr[,"s(PCB118)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(PCB118)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB118 <- bind_rows(df_PCB118, combined_df_PCB118)
combined_df_PCB118$Label <- as.factor(combined_df_PCB118$Label)


df_PCB138 <- data.frame(PCB138 = newdata[["PCB138"]], BKMR = gam.pred.smt.bkmr[,"s(PCB138)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(PCB138)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB138 <- bind_rows(df_PCB138, combined_df_PCB138)
combined_df_PCB138$Label <- as.factor(combined_df_PCB138$Label)


df_PCB153 <- data.frame(PCB153 = newdata[["PCB153"]], BKMR = gam.pred.smt.bkmr[,"s(PCB153)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(PCB153)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB153 <- bind_rows(df_PCB153, combined_df_PCB153)
combined_df_PCB153$Label <- as.factor(combined_df_PCB153$Label)


df_PCB180 <- data.frame(PCB180 = newdata[["PCB180"]], BKMR = gam.pred.smt.bkmr[,"s(PCB180)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(PCB180)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB180 <- bind_rows(df_PCB180, combined_df_PCB180)
combined_df_PCB180$Label <- as.factor(combined_df_PCB180$Label)


df_PCB187 <- data.frame(PCB187 = newdata[["PCB187"]], BKMR = gam.pred.smt.bkmr[,"s(PCB187)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(PCB187)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB187 <- bind_rows(df_PCB187, combined_df_PCB187)
combined_df_PCB187$Label <- as.factor(combined_df_PCB187$Label)


df_PCB194 <- data.frame(PCB194 = newdata[["PCB194"]], BKMR = gam.pred.smt.bkmr[,"s(PCB194)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(PCB194)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB194 <- bind_rows(df_PCB194, combined_df_PCB194)
combined_df_PCB194$Label <- as.factor(combined_df_PCB194$Label)


df_PCB169 <- data.frame(PCB169 = newdata[["PCB169"]], BKMR = gam.pred.smt.bkmr[,"s(PCB169)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(PCB169)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_PCB169 <- bind_rows(df_PCB169, combined_df_PCB169)
combined_df_PCB169$Label <- as.factor(combined_df_PCB169$Label)


df_Dioxin1 <- data.frame(Dioxin1 = newdata[["Dioxin1"]], BKMR = gam.pred.smt.bkmr[,"s(Dioxin1)"], 
                         SoftBART = gam.pred.smt.gsbart[,"s(Dioxin1)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_Dioxin1 <- bind_rows(df_Dioxin1, combined_df_Dioxin1)
combined_df_Dioxin1$Label <- as.factor(combined_df_Dioxin1$Label)


df_Dioxin2 <- data.frame(Dioxin2 = newdata[["Dioxin2"]], BKMR = gam.pred.smt.bkmr[,"s(Dioxin2)"], 
                         SoftBART = gam.pred.smt.gsbart[,"s(Dioxin2)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_Dioxin2 <- bind_rows(df_Dioxin2, combined_df_Dioxin2)
combined_df_Dioxin2$Label <- as.factor(combined_df_Dioxin2$Label)


df_Dioxin3 <- data.frame(Dioxin3 = newdata[["Dioxin3"]], BKMR = gam.pred.smt.bkmr[,"s(Dioxin3)"], 
                         SoftBART = gam.pred.smt.gsbart[,"s(Dioxin3)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_Dioxin3 <- bind_rows(df_Dioxin3, combined_df_Dioxin3)
combined_df_Dioxin3$Label <- as.factor(combined_df_Dioxin3$Label)


df_Furan1 <- data.frame(Furan1 = newdata[["Furan1"]], BKMR = gam.pred.smt.bkmr[,"s(Furan1)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(Furan1)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_Furan1 <- bind_rows(df_Furan1, combined_df_Furan1)
combined_df_Furan1$Label <- as.factor(combined_df_Furan1$Label)


df_Furan3 <- data.frame(Furan3 = newdata[["Furan3"]], BKMR = gam.pred.smt.bkmr[,"s(Furan3)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(Furan3)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_Furan3 <- bind_rows(df_Furan3, combined_df_Furan3)
combined_df_Furan3$Label <- as.factor(combined_df_Furan3$Label)


df_Furan4 <- data.frame(Furan4 = newdata[["Furan4"]], BKMR = gam.pred.smt.bkmr[,"s(Furan4)"], 
                        SoftBART = gam.pred.smt.gsbart[,"s(Furan4)"]) %>% 
  pivot_longer(cols = c(BKMR, SoftBART), names_to = "Model", values_to = "Pred") %>% 
  mutate(Type = "Model", Label = Model)

combined_df_Furan4 <- bind_rows(df_Furan4, combined_df_Furan4)
combined_df_Furan4$Label <- as.factor(combined_df_Furan4$Label)



common_theme <- theme(
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 12),
  legend.title = element_blank(),
  legend.text = element_text(size = 10),
  legend.position = "none",
  plot.title = element_blank()
)


label_map = c("SoftBART" = "GAM Approximation for modified BART",
              "BKMR" = "GAM Approximation for BKMR",
              "25th" = "BKMR Summarization with Remaining at 1st Quartile",
              "50th" = "BKMR Summarization with Remaining at Median",
              "75th" = "BKMR Summarization with Remaining at 3rd Quartile")


colors_def <- scales::hue_pal()(length(names(label_map)))
color_map <- c("SoftBART" = colors_def[length(colors_def)], 
               "BKMR" = colors_def[1], 
               "25th" = colors_def[2],
               "50th" = colors_def[3],
               "75th" = colors_def[4]
)



prepare_legend <- function(df, label_map){
  model_level <- sort(unique(df$Label[df$Type=="Model"]), decreasing = TRUE)
  quant_level <- sort(unique(df$Label[df$Type=="Quantile"]))
  custom_order <- c(model_level, quant_level)
  
  df$Label <- factor(df$Label, levels = custom_order)
  display_labels <- label_map[levels(df$Label)]
  
  list(data = df, labels = display_labels)
}

out <- prepare_legend(combined_df_PCB74, label_map)
p1 <- ggplot(out$data, aes(PCB74, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB170, label_map)
p2 <- ggplot(out$data, aes(PCB170, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB126, label_map)
p3 <- ggplot(out$data, aes(PCB126, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)




chem_dic <- list(Dioxin1 = "1,2,3,6,7,8-hxcdd",
                 Dioxin2 = "1,2,3,4,6,7,8-hpcdd",
                 Dioxin3 = "1,2,3,4,6,7,8,9-ocdd",
                 Furan1 = "2,3,4,7,8-pncdf",
                 Furan2 = "1,2,3,4,7,8-hxcdf",
                 Furan3 = "1,2,3,6,7,8-hxcdf",
                 Furan4 = "1,2,3,4,6,7,8-hxcdf")


out <- prepare_legend(combined_df_Furan2, label_map)
p4 <- ggplot(out$data, aes(Furan2, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = chem_dic[["Furan2"]]) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB99, label_map)
p5 <- ggplot(out$data, aes(PCB99, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB118, label_map)
p6 <- ggplot(out$data, aes(PCB118, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB138, label_map)
p7 <- ggplot(out$data, aes(PCB138, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB153, label_map)
p8 <- ggplot(out$data, aes(PCB153, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB180, label_map)
p9 <- ggplot(out$data, aes(PCB180, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB187, label_map)
p10 <- ggplot(out$data, aes(PCB187, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB194, label_map)
p11 <- ggplot(out$data, aes(PCB194, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_PCB169, label_map)
p12 <- ggplot(out$data, aes(PCB169, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_Dioxin1, label_map)
p13 <- ggplot(out$data, aes(Dioxin1, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = chem_dic[["Dioxin1"]]) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_Dioxin2, label_map)
p14 <- ggplot(out$data, aes(Dioxin2, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = chem_dic[["Dioxin2"]]) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_Dioxin3, label_map)
p15 <- ggplot(out$data, aes(Dioxin3, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = chem_dic[["Dioxin3"]]) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_Furan1, label_map)
p16 <- ggplot(out$data, aes(Furan1, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = chem_dic[["Furan1"]]) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_Furan3, label_map)
p17 <- ggplot(out$data, aes(Furan3, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = chem_dic[["Furan3"]]) + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)

out <- prepare_legend(combined_df_Furan4, label_map)
p18 <- ggplot(out$data, aes(Furan4, Pred, colour = Label)) +
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = chem_dic[["Furan4"]]) + common_theme + 
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

out_legend <- prepare_legend(combined_df_PCB74, label_map)
legend_only <- as_ggplot(get_legend(
  ggplot(out_legend$data, aes(PCB74, Pred, colour = Label)) +
    geom_line(data = out_legend$data %>% filter(Type == "Model"), linewidth = 1.2) + 
    geom_point(data = out_legend$data %>% filter(Type == "Model"), size = 3) + 
    geom_line(data = out_legend$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
    labs(y = "Marginal Effects") + common_theme + 
    scale_color_manual(values = color_map, labels = out$labels)
))


blank_plot <- plot_spacer()

allplots <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, 
                 p11, p12, p13, p14, p15, p16, p17, p18)

wrap_plots(c(allplots, list(legend_only, blank_plot)), nrow = 5, ncol = 4)







