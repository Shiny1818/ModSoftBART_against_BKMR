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



# for ungrouped regression.
## run BKMR.
start_bkmr = now()
fit_bkmr_knots100_ungrouped <-  kmbayes(y=lnLTL_z, Z=lnmixture_z, X=covariates, iter=100000, verbose=TRUE, varsel=TRUE, 
                                        groups=NULL, 
                                        knots=knots100)
time_bkmr_ungrouped <- as.duration(now() - start_bkmr)

## run gsBART
start_time <- now()
fit_gsbart_ungrouped_ntr20 <- gsbart_regression(X = as.matrix(covariates), Y = lnLTL_z, Z = as.matrix(lnmixture_z), num_tree = 20,
                                                num_burn = 50000, num_save = 1000, num_thin = 50, alpha_val = 1, alpha_comp = 1,
                                                group = NULL)
time_gsbart_ungrouped_ntr20 <- as.duration(now() - start_time)

## Calculate PIPs - ungrouped
PIPs_gsbart_ungrouped_ntr20 <- calPIPs_gsbart(fit = fit_gsbart_ungrouped_ntr20, group =  NULL, vars = lnmixture_z)

PIPs_gsbart_ungrouped_ntr50 <- calPIPs_gsbart(fit = fit_gsbart_ungrouped_ntr50, group =  NULL, vars = lnmixture_z)

PIPs_bkmr_ungrouped <- ExtractPIPs(fit_bkmr_knots100_ungrouped)


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
        # axis.title.y = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12))


################################################################################
# GAM results for NHANES data fitting results.
library(bkmr)
library(mgcv)
library(dplyr)
library(ggplot2)
library(grid)



