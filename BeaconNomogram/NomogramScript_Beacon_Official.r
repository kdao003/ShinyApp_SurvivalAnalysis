library(tidyverse)
library(ggplot2)
library(survival)

#### Analysis for c_JD samples ####

# Load data
df_beacon <- read.csv("BEACON_Master_Meta_Data_2846_relabeledSurvival.csv", header = TRUE, sep = ',')

#Remove ID Cols
df_beacon1 <- df_beacon %>%
  select(-contains("ID")) 


# Define your categorical predictors
cats <- c("tstage_c_JD", "nstage_c_JD", "mstage_c_JD")

# Find valid levels (levels with both 0 and 1 outcomes)
valid_levels <- df_beacon1 %>%
  filter(!is.na(clinical_outcome_JD_relabeled)) %>%
  group_by(across(all_of(cats)), clinical_outcome_JD_relabeled) %>%
  tally() %>%
  filter(n > 0) %>%
  tidyr::pivot_wider(
    names_from = clinical_outcome_JD_relabeled, 
    values_from = n, 
    values_fill = 0
  ) %>%
  filter(`0` > 0, `1` > 0) %>%   # keep only levels that have both outcomes
  select(all_of(cats)) 

#Filter the dataset to only these valid combinations
df_tmp <- df_beacon1 %>%
  filter(!is.na(tstage_c_JD),
         !is.na(nstage_c_JD),
         !is.na(mstage_c_JD),
         !is.na(clinical_outcome_JD_relabeled)) %>%
  semi_join(valid_levels, by = cats)

#Shows if the subcategorial varubales in these 3 predictor variables 
#have an associated 1 and 0 survival outcome
table(df_tmp$tstage_c_JD, df_tmp$clinical_outcome_JD_relabeled)
table(df_tmp$nstage_c_JD, df_tmp$clinical_outcome_JD_relabeled)
table(df_tmp$mstage_c_JD, df_tmp$clinical_outcome_JD_relabeled)

#Convert relevant columns to appropriate types
df_tmp$clinical_os_JD <- as.numeric(df_tmp$clinical_os_JD)
df_tmp$clinical_outcome_JD_relabeled <- as.numeric(df_tmp$clinical_outcome_JD_relabeled)
df_tmp$tstage_c_JD <- as.factor(df_tmp$tstage_c_JD)
df_tmp$nstage_c_JD <- as.factor(df_tmp$nstage_c_JD)
df_tmp$mstage_c_JD <- as.factor(df_tmp$mstage_c_JD)

#Remove rows that contain unknown in the predictor columns
df_tmp2 <- df_tmp %>%
  filter(
    !tstage_c_JD %in% c("unknown"),
    !nstage_c_JD %in% c("unknown"),
    !mstage_c_JD %in% c("unknown")
  ) %>%
  droplevels() %>%
  select(where(~ sum(!is.na(.)) >= 2))   # removes columns that has <2 non-missing values

#Creating survival function using cph (semiparametic survival model)
#Popular in medical literature (PMID: 28603726)
ddist <- datadist(df_tmp2)
options(datadist = 'ddist')
mod.cox <- cph(Surv(clinical_os_JD,clinical_outcome_JD_relabeled) ~  
                 tstage_c_JD + nstage_c_JD + mstage_c_JD, df_tmp2, 
               x = T, y = T, surv = T)
surv.cox <- Survival(mod.cox)

#creating/plotting nomogram
nom.cox1 <- nomogram(mod.cox, fun = list(function(x) surv.cox(12,x),
                                         function(x) surv.cox(24,x),
                                         function(x) surv.cox(36,x),
                                         function(x) surv.cox(48,x),
                                         function(x) surv.cox(60,x)),
                     funlabel = c("12-Month Survival Probability",
                                  "24-Month Survival Probability",
                                  "36-Month Survival Probability",
                                  "48-Month Survival Probability",
                                  "60-Month Survival Probability"),
                     lp = F)

plot(nom.cox1)






####Analysis 2#####

#Regroup subcategorical groups into larger groups
df_regrouped <- df_tmp2 %>%
  mutate(
    tstage_c_JD = case_when(
      tstage_c_JD %in% c("T2", "T2a", "T2b") ~ "T2",
      tstage_c_JD %in% c("T3", "T3a", "T3b") ~ "T3",
      tstage_c_JD %in% c("T4", "T4a", "T4b") ~ "T4",
      TRUE ~ tstage_c_JD
    )
  ) %>%
  mutate(tstage_c_JD = factor(tstage_c_JD))   # refresh factor levels
  
table(df_regrouped$tstage_c_JD, df_regrouped$clinical_outcome_JD_relabeled)
table(df_regrouped$nstage_c_JD, df_regrouped$clinical_outcome_JD_relabeled)

#plot survival plots of each subcategorical group
fit <- survfit(Surv(clinical_os_JD, clinical_outcome_JD_relabeled) ~ tstage_c_JD, 
               data = df_regrouped, type = "kaplan-meier")
summary(fit)
ggsurvplot(
  fit,
  data = df_regrouped,
  risk.table = TRUE,
  conf.int = FALSE,
  pval = TRUE,
  ggtheme = theme_minimal(),
  palette = "Dark2",
  linetype = 1,           # fixed, not mapped
  linewidth = 0.8,        # fixed width
  conf.int.fill = "gray80",  # use fill instead of mapping CI by group
  conf.int.alpha = 0.25,     # transparency for CI
  conf.int.style = "ribbon"  # keep CI bands simple
)
#Remove Tis and Ta variables


fit1 <- survfit(Surv(clinical_os_JD, clinical_outcome_JD_relabeled) ~ nstage_c_JD, 
               data = df_regrouped, type = "kaplan-meier")
summary(fit1)

ggsurvplot(
  fit1,
  data = df_regrouped,
  risk.table = TRUE,
  conf.int = FALSE,
  pval = TRUE,
  ggtheme = theme_minimal(),
  palette = "Dark2",
  linetype = 1,           # fixed, not mapped
  linewidth = 0.8,        # fixed width
  conf.int.fill = "gray80",  # use fill instead of mapping CI by group
  conf.int.alpha = 0.25,     # transparency for CI
  conf.int.style = "ribbon"  # keep CI bands simple
)
#Remove N+ variable

fit2 <- survfit(Surv(clinical_os_JD, clinical_outcome_JD_relabeled) ~ mstage_c_JD, 
               data = df_regrouped, type = "kaplan-meier")
summary(fit2)

ggsurvplot(
  fit2,
  data = df_regrouped,
  risk.table = TRUE,
  conf.int = FALSE,
  pval = TRUE,
  ggtheme = theme_minimal(),
  palette = "Dark2",
  linetype = 1,           # fixed, not mapped
  linewidth = 0.8,        # fixed width
  conf.int.fill = "gray80",  # use fill instead of mapping CI by group
  conf.int.alpha = 0.25,     # transparency for CI
  conf.int.style = "ribbon"  # keep CI bands simple
)
#All variables look fine, no need to remove any

#Omit N+ values, N+ is only 1 sample, survival plot drops off at beginning
#Omit Tis and Ta, too little sample to have valuable meaning
df_regrouped_omitted <- df_regrouped %>%
 filter(!nstage_c_JD %in% "N+",
  !tstage_c_JD %in% c("Tis", "Ta")) %>%
droplevels()

#Plotting new Nomogram after additional grouping and and removing nonmeaningful data
ddist <- datadist(df_regrouped_omitted)
options(datadist = 'ddist')
mod.cox2 <- cph(Surv(clinical_os_JD,clinical_outcome_JD_relabeled) ~  
                 tstage_c_JD + nstage_c_JD + mstage_c_JD, df_regrouped_omitted, 
               x = T, y = T, surv = T)
surv.cox <- Survival(mod.cox2)

#creating/plotting nomogram
nom.cox2 <- nomogram(mod.cox2, fun = list(function(x) surv.cox(12,x),
                                         function(x) surv.cox(24,x),
                                         function(x) surv.cox(36,x),
                                         function(x) surv.cox(48,x),
                                         function(x) surv.cox(60,x)),
                     funlabel = c("12-Month Survival Probability",
                                  "24-Month Survival Probability",
                                  "36-Month Survival Probability",
                                  "48-Month Survival Probability",
                                  "60-Month Survival Probability"),
                     lp = F)

plot(nom.cox2)

#Group by TURBT
df_TURBT <- df_regrouped_omitted %>%
  filter(turbt_cystectomy_JD %in% "TURBT")

#Plotting survival curve after TURBT grouping and omitting of subcategorical variables
TURBT_omit_plot1 <- survfit(Surv(clinical_os_JD, clinical_outcome_JD_relabeled) ~ mstage_c_JD, 
                     data = df_TURBT, type = "kaplan-meier")
summary(TURBT_omit_plot1)

ggsurvplot(
  TURBT_omit_plot1,
  data = df_TURBT,
  risk.table = TRUE,
  conf.int = FALSE,
  pval = TRUE,
  ggtheme = theme_minimal(),
  palette = "Dark2",
  linetype = 1,           # fixed, not mapped
  linewidth = 0.8,        # fixed width
  conf.int.fill = "gray80",  # use fill instead of mapping CI by group
  conf.int.alpha = 0.25,     # transparency for CI
  conf.int.style = "ribbon"  # keep CI bands simple
)

#Nomogram for TURBT
ddist <- datadist(df_TURBT)
options(datadist = 'ddist')
turbt.cox <- cph(Surv(clinical_os_JD,clinical_outcome_JD_relabeled) ~  
                  tstage_c_JD + nstage_c_JD + mstage_c_JD, df_TURBT, 
                x = T, y = T, surv = T)
turbt.surv.cox <- Survival(turbt.cox)

#creating/plotting nomogram
turbt.nom.cox <- nomogram(turbt.cox, fun = list(function(x) turbt.surv.cox(12,x),
                                          function(x) turbt.surv.cox(24,x),
                                          function(x) turbt.surv.cox(36,x),
                                          function(x) turbt.surv.cox(48,x),
                                          function(x) turbt.surv.cox(60,x)),
                     funlabel = c("12-Month Survival Probability",
                                  "24-Month Survival Probability",
                                  "36-Month Survival Probability",
                                  "48-Month Survival Probability",
                                  "60-Month Survival Probability"),
                     lp = F)

plot(turbt.nom.cox)


x


###Analysis for p_JD samples###

df_p_JD <- df_beacon %>%
  select(-contains("ID")) 

# Define your categorical predictors
cats_p_JD <- c("tstage_p_JD", "nstage_p_JD", "mstage_p_JD")

# Find valid levels (levels with both 0 and 1 outcomes)
valid_levels_p_JD <- df_p_JD %>%
  filter(!is.na(clinical_outcome_JD_relabeled)) %>%
  group_by(across(all_of(cats_p_JD)), clinical_outcome_JD_relabeled) %>%
  tally() %>%
  filter(n > 0) %>%
  tidyr::pivot_wider(
    names_from = clinical_outcome_JD_relabeled, 
    values_from = n, 
    values_fill = 0
  ) %>%
  filter(`0` > 0, `1` > 0) %>%   # keep only levels that have both outcomes
  select(all_of(cats_p_JD)) 

#Filter the dataset to only these valid combinations
df_p_JD1 <- df_p_JD  %>%
  filter(!is.na(tstage_p_JD),
         !is.na(nstage_p_JD),
         !is.na(mstage_p_JD),
         !is.na(clinical_outcome_JD_relabeled)) %>%
  semi_join(valid_levels_p_JD, by = cats_p_JD) %>%
  filter(
    !tstage_p_JD %in% c("unknown"), #Remove rows that contain unknown in the predictor columns
    !nstage_p_JD %in% c("unknown"),
    !mstage_p_JD %in% c("unknown")
  ) %>% #filter to only rows that have cystectomy
  droplevels() %>%
  select(where(~ sum(!is.na(.)) >= 2))  # removes columns that has <2 non-missing values

#Shows if the subcategorial varubales in these 3 predictor variables 
#have an associated 1 and 0 survival outcome
table(df_p_JD1$tstage_p_JD, df_p_JD1$clinical_outcome_JD_relabeled)
table(df_p_JD1$nstage_p_JD, df_p_JD1$clinical_outcome_JD_relabeled)
table(df_p_JD1$mstage_p_JD, df_p_JD1$clinical_outcome_JD_relabeled)

#Convert relevant columns to appropriate types
df_p_JD1$clinical_os_JD <- as.numeric(df_p_JD1$clinical_os_JD)
df_p_JD1$clinical_outcome_JD_relabeled <- as.numeric(df_p_JD1$clinical_outcome_JD_relabeled)
df_p_JD1$tstage_p_JD <- as.factor(df_p_JD1$tstage_p_JD)
df_p_JD1$nstage_p_JD <- as.factor(df_p_JD1$nstage_p_JD)
df_p_JD1$mstage_p_JD <- as.factor(df_p_JD1$mstage_p_JD)

#Grouping together T variables into larger categories
df_p_JD2 <- df_p_JD1 %>%
  mutate(
    tstage_p_JD = case_when(
      tstage_p_JD %in% c("T2", "T2b") ~ "T2",
      tstage_p_JD %in% c("T3", "T3a", "T3b") ~ "T3",
      tstage_p_JD %in% c("T4", "T4a") ~ "T4",
      TRUE ~ tstage_p_JD
    )
  ) %>%
  mutate(tstage_p_JD = factor(tstage_p_JD))

#Plot survival curves and summary for each predictor p_JD variable
fit_tpJD <- survfit(Surv(clinical_os_JD, clinical_outcome_JD_relabeled) ~ tstage_p_JD, 
                data = df_p_JD2, type = "kaplan-meier")
summary(fit_tpJD)

fit_npJD <- survfit(Surv(clinical_os_JD, clinical_outcome_JD_relabeled) ~ nstage_p_JD, 
                    data = df_p_JD2, type = "kaplan-meier")
summary(fit_npJD)

fit_mpJD <- survfit(Surv(clinical_os_JD, clinical_outcome_JD_relabeled) ~ mstage_p_JD, 
                    data = df_p_JD2, type = "kaplan-meier")
summary(fit_mpJD)

ggsurvplot(
  fit_mpJD,
  data = df_p_JD2,
  risk.table = TRUE,
  conf.int = FALSE,
  pval = TRUE,
  ggtheme = theme_minimal(),
  palette = "Dark2",
  linetype = 1,           # fixed, not mapped
  linewidth = 0.8,        # fixed width
  conf.int.fill = "gray80",  # use fill instead of mapping CI by group
  conf.int.alpha = 0.25,     # transparency for CI
  conf.int.style = "ribbon"  # keep CI bands simple
)

#Remove Ta subcategorical variable due to 1 sample
df_p_JD3 <- df_p_JD2 %>%
  filter(!tstage_p_JD %in% "Ta")%>%
  droplevels() 

#Omitted mstage_p_JD category due to having <2 category levels
#Creating survival function using cph (semiparametic survival model)
#Popular in medical literature (PMID: 28603726)
ddist <- datadist(df_p_JD3)
options(datadist = 'ddist')
p_JD.cox <- cph(Surv(clinical_os_JD,clinical_outcome_JD_relabeled) ~  
                 tstage_p_JD + nstage_p_JD, df_p_JD3, 
               x = T, y = T, surv = T)
p_JD.surv.cox <- Survival(p_JD.cox)

#creating/plotting nomogram
pJD_nom.cox <- nomogram(p_JD.cox, fun = list(function(x) p_JD.surv.cox(12,x),
                                         function(x) p_JD.surv.cox(24,x),
                                         function(x) p_JD.surv.cox(36,x),
                                         function(x) p_JD.surv.cox(48,x),
                                         function(x) p_JD.surv.cox(60,x)),
                     funlabel = c("12-Month Survival Probability",
                                  "24-Month Survival Probability",
                                  "36-Month Survival Probability",
                                  "48-Month Survival Probability",
                                  "60-Month Survival Probability"),
                     lp = F)

plot(pJD_nom.cox)


#####Cystectomy Grouping####
df_cys <- df_p_JD %>%
  filter(turbt_cystectomy_JD %in% "Cystectomy")

# Find valid levels (levels with both 0 and 1 outcomes)
valid_levels_p_JD <- df_p_JD %>%
  filter(!is.na(clinical_outcome_JD_relabeled)) %>%
  group_by(across(all_of(cats_p_JD)), clinical_outcome_JD_relabeled) %>%
  tally() %>%
  filter(n > 0) %>%
  tidyr::pivot_wider(
    names_from = clinical_outcome_JD_relabeled, 
    values_from = n, 
    values_fill = 0
  ) %>%
  filter(`0` > 0, `1` > 0) %>%   # keep only levels that have both outcomes
  select(all_of(cats_p_JD)) 

#keep unknowns in mstage so we keep data in the tstage and nstage
df_cys1 <- df_cys %>%
  semi_join(valid_levels_p_JD, by = cats_p_JD) %>%
  filter(!is.na(tstage_p_JD),
         !is.na(nstage_p_JD),
         !is.na(mstage_p_JD),
         !is.na(clinical_outcome_JD_relabeled)) %>%
  filter(
    !tstage_p_JD %in% c("unknown"), #Remove rows that contain unknown in the predictor columns
    !nstage_p_JD %in% c("unknown"),
    #!mstage_p_JD %in% c("unknown")
  ) %>% #filter to only rows that have cystectomy
  droplevels() %>%
  select(where(~ sum(!is.na(.)) >= 2)) # removes columns that has <2 non-missing values

#Converting to specific class types
df_cys1$clinical_os_JD <- as.numeric(df_cys1$clinical_os_JD)
df_cys1$clinical_outcome_JD_relabeled <- as.numeric(df_cys1$clinical_outcome_JD_relabeled)
df_cys1$tstage_p_JD <- as.factor(df_cys1$tstage_p_JD)
df_cys1$nstage_p_JD <- as.factor(df_cys1$nstage_p_JD)
df_cys1$mstage_p_JD <- as.factor(df_cys1$mstage_p_JD)


#Plotting survival plot for tstage and nstage predictors
fit_cys_tpJD <- survfit(Surv(clinical_os_JD, clinical_outcome_JD_relabeled) ~ tstage_p_JD, 
                    data = df_cys1, type = "kaplan-meier")

summary(fit_cys_tpJD)

fit_cys_npJD <- survfit(Surv(clinical_os_JD, clinical_outcome_JD_relabeled) ~ nstage_p_JD, 
                    data = df_cys1, type = "kaplan-meier")

summary(fit_cys_npJD)

ggsurvplot(
  fit_cys_npJD,
  data = df_cys1,
  risk.table = TRUE,
  conf.int = FALSE,
  pval = TRUE,
  ggtheme = theme_minimal(),
  palette = "Dark2",
  linetype = 1,           # fixed, not mapped
  linewidth = 0.8,        # fixed width
  conf.int.fill = "gray80",  # use fill instead of mapping CI by group
  conf.int.alpha = 0.25,     # transparency for CI
  conf.int.style = "ribbon"  # keep CI bands simple
)
  
#Nomogram for Cystectomy Patients
ddist <- datadist(df_cys1)
options(datadist = 'ddist')
cys_p_JD.cox <- cph(Surv(clinical_os_JD,clinical_outcome_JD_relabeled) ~  
                  tstage_p_JD + nstage_p_JD, df_cys1, 
                x = T, y = T, surv = T)
cys_p_JD.surv.cox <- Survival(cys_p_JD.cox)

#creating/plotting nomogram
cys_pJD_nom.cox <- nomogram(cys_p_JD.cox, fun = list(function(x) cys_p_JD.surv.cox(12,x),
                                             function(x) cys_p_JD.surv.cox(24,x),
                                             function(x) cys_p_JD.surv.cox(36,x),
                                             function(x) cys_p_JD.surv.cox(48,x),
                                             function(x) cys_p_JD.surv.cox(60,x)),
                        funlabel = c("12-Month Survival Probability",
                                     "24-Month Survival Probability",
                                     "36-Month Survival Probability",
                                     "48-Month Survival Probability",
                                     "60-Month Survival Probability"),
                        lp = F)

plot(cys_pJD_nom.cox)

table(df_cys1$tstage_p_JD, df_cys1$clinical_outcome_JD_relabeled)
table(df_cys1$nstage_p_JD, df_cys1$clinical_outcome_JD_relabeled)




