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
fit_gsbart_grouped_ntr20 <- gsbart_regression(X = as.matrix(covariates), Y = lnLTL_z, Z = as.matrix(lnmixture_z), num_tree = 20,
                                              num_burn = 50000, num_save = 1000, num_thin = 50, alpha_val = 1, alpha_comp = 1,
                                              group = c(rep(1,times=2), 2, rep(1,times=6), rep(3,times=2),rep(2,times=7)))
time_gsbart_grouped_ntr20 <- as.duration(now() - start_time)

## Calculate PIPs - grouped
PIPs_gsbart_grouped_ntr20 <- calPIPs_gsbart(fit = fit_gsbart_grouped_ntr20, group = c(rep(1,times=2), 2, rep(1,times=6), rep(3,times=2),rep(2,times=7)),
                                            vars = lnmixture_z)

PIPs_bkmr_grouped <- ExtractPIPs(fit_bkmr_knots100_grouped)



################################################################################
### plot out the PIPs
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
fit_gsbart <- fit_gsbart_grouped_ntr20; htest_bkmr_mean <- htest_bkmr_grp_mean; htest_bkmr_sampl <- htest_bkmr_grp_sampl
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

# plot smooth items
gam.pred.bkmr <- predict(GAMfit_bkmr_grp, newdata = as.data.frame(lnmixture_z), type = "terms",
                         se.fit = TRUE)
gam.pred.gsbart <- predict(GAMfit_gsbart_grp, newdata = as.data.frame(lnmixture_z), type = "terms",
                           se.fit = TRUE)


gam.pred.smt.bkmr <- gam.pred.bkmr$fit
gam.pred.smt.gsbart <- gam.pred.gsbart$fit

fit_bkmr_knots100 <- fit_bkmr_knots100_grouped

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

out <- prepare_legend(combined_df_Furan2, label_map)
p4 <- ggplot(out$data, aes(Furan2, Pred, colour = Label)) +
  
  geom_line(data = out$data %>% filter(Type == "Model"), linewidth = 1.2) + 
  geom_point(data = out$data %>% filter(Type == "Model"), size = 3) + 
  geom_line(data = out$data %>% filter(Type == "Quantile"), linewidth = 0.8) +
  labs(y = "Marginal Effects", x = "1,2,3,4,7,8-hxcdf") + common_theme + 
  scale_color_manual(values = color_map, labels = out$labels)


chem_dic <- list(Dioxin1 = "1,2,3,6,7,8-hxcdd",
                 Dioxin2 = "1,2,3,4,6,7,8-hpcdd",
                 Dioxin3 = "1,2,3,4,6,7,8,9-ocdd",
                 Furan1 = "2,3,4,7,8-pncdf",
                 Furan2 = "1,2,3,4,7,8-hxcdf",
                 Furan3 = "1,2,3,6,7,8-hxcdf",
                 Furan4 = "1,2,3,4,6,7,8-hxcdf")


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

wrap_plots(list(p3, legend_only, p16, blank_plot), nrow = 2, ncol = 2)







