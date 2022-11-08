## R code for n-3PUFA Meta-analysis submitted to Advances in Nutrition.
# Code authors: Gabriel P Esteves, Gabriel H. C. Barreto.
# Contact: gabriel.perri.esteves@usp.br; gabriel.castanhob@usp.br.

## Please note that this code was executed using R version 4.2.0,
## and Rstudio version 2022.07.2. Both pipes (the base |> and {magrittr} %>% pipes)
## are used throughout the code, so make sure you have the necessary packages and updated 
## version of R in order to use both correctly.

## Loading Packages

library(readxl); library(metafor); library(tidyverse); library(ggplot2); 
library(Cairo); library(viridis); library(here); library(broom); library(writexl)
library(ggsci); library(ggforestplot); library(patchwork); library(glue)

options(scipen=999) # removing scientific notation

## Functions

# dppc2 and var, see Morris (2007)

prepostcon <- function(rho, meanpretreat, meanposttreat, meanprecon, meanpostcon, sdpretreat, sdprecon, ntreat, ncon){
  pooledsd =  sqrt ((((ntreat-1)*sdpretreat^2) + ((ncon-1)*sdprecon^2) ) / (ntreat + ncon - 2))
  
  correction <- 1 - (3 / ((4 * (ntreat + ncon - 2)) - 1) )
  
  ES <-  correction * (((meanposttreat - meanpretreat) - (meanpostcon - meanprecon)) / pooledsd)
  
  var <- (2 * (correction^2)) * (1 - rho) * ((ntreat + ncon) / (ntreat*ncon)) * ((ntreat + ncon - 2) / (ntreat + ncon - 4)) * (1 + (ES^2 / (2 * (1 - rho) * ((ntreat+ncon)/(ntreat*ncon))))) - ES^2
  
  se <- sqrt(var)
  lowCI <- ES - 1.96*se
  highCI <- ES + 1.96*se
  
  return(list(ES, var, se, lowCI, highCI))
}

# compute I^2 and CIs. See https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-August/000150.html

compI2 <- function(res) {
  sav <- confint.rma.mv(res)
  W <- diag(1/res$vi)
  X <- model.matrix(res)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2 <- tibble(100 * res$sigma2 / (res$sigma2 + (res$k-res$p)/sum(diag(P)))) |> 
    rename("I2" = everything()) |> 
    mutate(level=c("study", "outcome"))
  ci.study <- enframe(100 * sav[[1]]$random[1,2:3] / (sav[[1]]$random[1,2:3] + (res$k-res$p)/sum(diag(P)))) |> mutate(level="study")### CI for district-level I^2
  ci.out <- enframe(100 * sav[[2]]$random[1,2:3] / (sav[[2]]$random[1,2:3] + (res$k-res$p)/sum(diag(P)))) |> mutate(level="outcome")
  result <- full_join(ci.study, ci.out) |> pivot_wider() |> full_join(I2) |> select(I2, ci.lb, ci.ub, level)
  return(result)
}

## Data

n3data <- read_excel(here("Data/n3PUFA-data.xlsx")) |># place the spreadsheet into the location of the script or R project
  mutate(authoryear = paste(author, " ", "(", refnumber, ")", sep="")) 
n3data <- n3data[order(n3data$doseclassification, n3data$year),] # reordering for later when we make forest plots

# inverting values for function tests where less is better

n3data$meanpre_pla[n3data$units == "seconds"] <- -n3data$meanpre_pla[n3data$units == "seconds"]
n3data$meanpre_n3 [n3data$units == "seconds"] <- -n3data$meanpre_n3 [n3data$units == "seconds"]
n3data$meanpost_pla[n3data$units == "seconds"] <- -n3data$meanpost_pla[n3data$units == "seconds"]
n3data$meanpost_n3[n3data$units == "seconds"] <- -n3data$meanpost_n3[n3data$units == "seconds"]

# Estimating a correlation value via pilot data

# In order to calculate SMD for pre-post designs, we need a pre-post correlation value.
# Here we used a pilot data-set from an ongoing RCT from our research group to
# derive these values. Dataframes below are built with these exact values.

df_lean <- tibble(x=c(57.62,
                      65.431,
                      57.526,
                      72.228,
                      69.114,
                      61.509,
                      48.874,
                      57.825),
                  y=c(60.169,
                      65.953,
                      57.492,
                      71.079,
                      69.591,
                      63.455,
                      51.734,
                      57.958))

df_rm <- tibble(x=c(260,
                    270,
                    255,
                    450,
                    350,
                    280,
                    350,
                    315,
                    380),
                y=c(345,
                    390,
                    395,
                    520,
                    435,
                    470,
                    430,
                    510,
                    435))

lean_cor <- cor.test(df_lean$x, df_lean$y, method="pearson")
lean_cor$estimate # correlation for muscle mass outcome.

str_cor <- cor.test(df_rm$x, df_rm$y, method="pearson")
str_cor$estimate # correlation for strength outcome.

# We'll be using:
# correlation for lean mass: 0.90 (slightly more conservative than observed);
# correlation for strength: 0.65 (observed value from data-set);
# correlation for function: 0.65 (assuming a similar correlation to strength tests);

n3data <- n3data |> 
  mutate(ri = case_when(outcome_type == "body composition" ~ 0.9,
                        outcome_type == "strength" ~ 0.65,
                        outcome_type == "function" ~ 0.65))

## calculating SMDs
# We'll use our function made earlier to calculate each SMD and join them with the main dataframe

df2 <- data.frame(prepostcon(n3data$ri, 
                             n3data$meanpre_n3, 
                             n3data$meanpost_n3, 
                             n3data$meanpre_pla, 
                             n3data$meanpost_pla, 
                             n3data$sdpre_n3, 
                             n3data$sdpre_pla,
                             n3data$ntotal_n3, 
                             n3data$ntotal_pla)) %>%
  `colnames<-`(c("yi_2", "vi_2", "se", "lowCI", "highCI")) 

n3data <- bind_cols(n3data, df2)

## Meta-analysis
# We'll now run and save each meta-analysis model according to outcome type

# outcome: muscle mass (referred to in the code as 'lean mass' or 'muscle mass')

n3data_mass <- n3data %>% 
  filter(outcome_type=="body composition")

# We'll be using the rma.mv function in order to create a three-level meta-analysis
# This can be done through the 'random' argument, as shown below.
# For more information on three-level models with metafor, see the rma.mv documentation or 
# https://www.metafor-project.org/doku.php/analyses:konstantopoulos2011

# Here we use the structure authoryear/es.id, which informs the function that 
# our effect sizes are nested within studies.
# The same results could be obtained using the 'paperid' column instead of 'authoryear'.

meta_lean_mass <- rma.mv(yi = yi_2, 
                     V = vi_2, 
                     slab = authoryear,
                     data = n3data_mass,
                     random = ~ 1 | authoryear/es.id, 
                     test = "t", 
                     method = "REML")

# Creating a tidy summary of meta results for later!
lean <- tidy(meta_lean_mass, conf.int=T, include_studies=T) %>% mutate(outcome="lean_mass") 

# A quick peak into what results look like.
forest(meta_lean_mass)

# Counting number of studies in each outcome to see what moderator analyses will be done.
n3data_mass %>%
  group_by(doseclassification) %>% summarise(n())

n3data_mass %>%
  group_by(ageclass) %>% summarise(n())

n3data_mass %>%
  group_by(intervention_type2) %>% summarise(n()) 

# moderators analysis for intervention type only
# In order to accomplish this, we use the 'mods' argument within rma.mv
meta_lean_mass_training <- rma.mv(yi = yi_2, 
                         V = vi_2, 
                         mods = ~ factor(intervention_type2), 
                         slab = authoryear,
                         data = n3data_mass,
                         random = ~ 1 | authoryear/es.id, 
                         test = "t", 
                         method = "REML")

lean_training <- tidy(meta_lean_mass_training, conf.int=T, include_studies=T) %>% mutate(outcome="lean_mass_training")

# outcome: strength

n3data_str <- n3data %>%
  filter(outcome_type=="strength")

meta_str <- rma.mv(yi_2, 
                 vi_2,
                 slab=authoryear,
                 random= ~ 1 | authoryear/es.id,
                 method="REML",
                 test="t",
                 data=n3data_str)

forest(meta_str)

str <- tidy(meta_str, conf.int=T, include_studies=T) %>% mutate(outcome="str")

n3data_str %>%
  group_by(doseclassification) %>% summarise(n())

n3data_str %>%
  group_by(ageclass) %>% summarise(n())

n3data_str %>%
  group_by(intervention_type2) %>% summarise(n()) 

# moderator analyses for all factors
meta_str_dose <- rma.mv(yi_2, 
                        vi_2,
                        mods = ~ factor(doseclassification),
                        slab=authoryear,
                        random= ~ 1 | authoryear/es.id,
                        method="REML",
                        test="t",
                        data=n3data_str)

# quick sensitivity test by removing Cornish & Chilibeck (2009)
n3 <- n3data_str %>%
  filter(authoryear != "Cornish & Chilibeck (2009)")

meta_str_dose2 <- rma.mv(yi_2, 
                        vi_2,
                        mods = ~ factor(doseclassification),
                        slab=authoryear,
                        random= ~ 1 | authoryear/es.id,
                        method="REML",
                        test="t",
                        data=n3)

summary(meta_str_dose2)
summary(meta_str_dose)

# As expected, no real change in results, so we'll go with the initial model.
str_dose <- tidy(meta_str_dose, conf.int=T, include_studies=T) %>% mutate(outcome="str_dose")

# continuing with the other sub-analyses:
meta_str_age <- rma.mv(yi_2, 
                       vi_2,
                       mods = ~ factor(ageclass),
                       slab=authoryear,
                       random= ~ 1 | authoryear/es.id,
                       method="REML",
                       test="t",
                       data=n3data_str)

summary(meta_str_age)
str_age <- tidy(meta_str_age, conf.int=T, include_studies=T) %>% mutate(outcome="str_age")

meta_str_training <- rma.mv(yi_2, 
                            vi_2,
                            mods = ~ factor(intervention_type2),
                            slab=authoryear,
                            random= ~ 1 | authoryear/es.id,
                            method="REML",
                            test="t",
                            data=n3data_str)

summary(meta_str_training)
str_training <- tidy(meta_str_training, conf.int=T, include_studies=T) %>% mutate(outcome="str_training")

# outcome: function

n3data_func <- n3data %>%
  filter(outcome_type=="function")

meta_func <- rma.mv(yi_2, 
                     vi_2, 
                     slab = authoryear,
                     data = n3data_func,
                     random = ~ 1 | authoryear/es.id, 
                     test = "t", 
                     method = "REML")

forest(meta_func)

func <- tidy(meta_func, conf.int=T, include_studies=T) %>% mutate(outcome="func")

n3data_func %>%
  group_by(doseclassification) %>% summarise(n())

n3data_func %>%
  group_by(ageclass) %>% summarise(n())

n3data_func %>%
  group_by(intervention_type2) %>% summarise(n()) 

# moderator analysis for dose and intervention type
meta_func_dose <- rma.mv(yi_2, 
                          vi_2,
                          mods = ~ factor(doseclassification),
                          slab = authoryear,
                          data = n3data_func,
                          random = ~ 1 | authoryear/es.id, 
                          test = "t", 
                          method = "REML")

func_dose <- tidy(meta_func_dose, conf.int=T, include_studies=T) %>% mutate(outcome="func_dose")

meta_func_training <- rma.mv(yi_2, 
                         vi_2, 
                         mods = ~ factor(intervention_type2),
                         slab = authoryear,
                         data = n3data_func,
                         random = ~ 1 | authoryear/es.id, 
                         test = "t", 
                         method = "REML")

func_training <- tidy(meta_func_training, conf.int=T, include_studies=T) %>% mutate(outcome="func_training")

## Making a table summary of all models (Table 2 in the paper as of current version)

# get variances

var_meta_lean_mass <- enframe(as.data.frame(confint(meta_lean_mass))) |> mutate(type="lean_mass")
var_meta_lean_mass_training <- enframe(as.data.frame(confint(meta_lean_mass_training))) |> mutate(type="lean_mass_training")
var_meta_str <- enframe(as.data.frame(confint(meta_str))) |> mutate(type="str")
var_meta_str_dose <- enframe(as.data.frame(confint(meta_str_dose))) |> mutate(type="str_dose")
var_meta_str_training <- enframe(as.data.frame(confint(meta_str_training))) |> mutate(type="str_training")
var_meta_str_age <- enframe(as.data.frame(confint(meta_str_age))) |> mutate(type="str_age")
var_meta_func <- enframe(as.data.frame(confint(meta_func))) |> mutate(type="func")
var_meta_func_dose <- enframe(as.data.frame(confint(meta_func_dose))) |> mutate(type="func_dose")
var_meta_func_training <- enframe(as.data.frame(confint(meta_func_training))) |> mutate(type="func_training")

var_all <- reduce(list(var_meta_lean_mass,
                       var_meta_lean_mass_training,
                       var_meta_str,
                       var_meta_str_dose,
                       var_meta_str_training,
                       var_meta_str_age,
                       var_meta_func,
                       var_meta_func_dose,
                       var_meta_func_training), full_join) |> 
  mutate_if(is.numeric, round, 3) |> 
  mutate(value = paste(value[,"estimate"], " ","(",value[,"ci.lb"],";"," ",value[,"ci.ub"],")", sep=""))

var_outcome <- var_all |> filter(name=="sigma^2.2")
var_study <- var_all |> filter(name=="sigma^2.1")

# make the table

meta_summary <- reduce(list(lean, lean_training, str, str_dose, str_training, str_age,
                            func, func_dose, func_training), full_join) |> 
  filter(type == "summary") |> 
  mutate(estimate = round(estimate, 3)) |>  
  mutate(conf.low = round(conf.low, 3)) |> 
  mutate(conf.high = round(conf.high, 3)) |> 
  mutate(std.error = round(std.error, 3)) |> 
  mutate(result = paste(estimate, " ", "(", conf.low, ";", " ", conf.high, ")", sep="")) |> 
  select(outcome, term, result, std.error) |> 
  mutate(Fstat = c(NA_real_, 
                   paste(round(meta_lean_mass_training$QM, 2), meta_lean_mass_training$QMdf),
                   NA_real_,
                   paste(round(meta_str_dose$QM, 2), meta_str_dose$QMdf),
                   paste(round(meta_str_training$QM, 2), meta_str_training$QMdf),
                   paste(round(meta_str_age$QM, 2), meta_str_age$QMdf),
                   NA_real_,
                   paste(round(meta_func_dose$QM, 2), meta_func_dose$QMdf),
                   paste(round(meta_func_training$QM, 2), meta_func_training$QMdf)
                   ),
         outcome_var = c(var_outcome$value[1], 
                         var_outcome$value[2], 
                         NA_real_,
                         var_outcome$value[3],
                         var_outcome$value[4],
                         NA_real_,
                         var_outcome$value[5],
                         NA_real_,
                         var_outcome$value[6],
                         NA_real_,
                         var_outcome$value[7],
                         var_outcome$value[8],
                         NA_real_,
                         var_outcome$value[9],
                         NA_real_),
         study_var = c(var_study$value[1], 
                       var_study$value[2], 
                       NA_real_,
                       var_study$value[3],
                       var_study$value[4],
                       NA_real_,
                       var_study$value[5],
                       NA_real_,
                       var_study$value[6],
                       NA_real_,
                       var_study$value[7],
                       var_study$value[8],
                       NA_real_,
                       var_study$value[9],
                       NA_real_),
         QE = c(meta_lean_mass$QE,
         meta_lean_mass_training$QE,
         NA_real_,
         meta_str$QE,
         meta_str_dose$QE,
         NA_real_,
         meta_str_training$QE,
         NA_real_,
         meta_str_age$QE,
         NA_real_,
         meta_func$QE,
         meta_func_dose$QE,
         NA_real_,
         meta_func_training$QE,
         NA_real_))
  
# write_xlsx(meta_summary, here("Tables/table_summary.xlsx")) | code for saving the table, commented out!

## Plots

# ggplot forest plot

# get MA study weights

w_lean <- enframe(weights(meta_lean_mass)) %>% rename(term = name, weight = value) %>%
  mutate(outcome = "lean_mass")

w_str <- enframe(weights(meta_str)) %>% rename(term = name, weight = value) %>%
  mutate(outcome = "str")

w_func <- enframe(weights(meta_func)) %>% rename(term = name, weight = value) %>%
  mutate(outcome = "func")

w_all <- reduce(list(w_lean, w_str, w_func), full_join)

# make a dataframe

df_forest <- reduce(list(lean, str, func), full_join) %>%
  mutate(term = as.factor(case_when(term == "overall" ~ "Overall", 
                                     term != "overall" ~ term))) %>%
  mutate(term = as.factor(term)) %>%
  filter(!is.na(estimate)) %>%
  mutate(estimate = format(round(estimate, digits=2), nsmall = 2)) %>%
  mutate(conf.low = format(round(conf.low, digits=2), nsmall = 2)) %>%
  mutate(conf.high = format(round(conf.high, digits=2), nsmall = 2)) %>%
  mutate(result = paste(estimate, " ", "[", conf.low, ";", " ", conf.high, "]", sep="")) %>%
  mutate(estimate = as.numeric(estimate)) %>%
  mutate(conf.low = as.numeric(conf.low)) %>%
  mutate(conf.high = as.numeric(conf.high)) %>%
  full_join(w_all, by=c("term", "outcome")) %>%
  mutate(weight = case_when(term=="Overall" ~ 40, TRUE ~ weight)) %>%
  mutate(term = reorder(term, desc(term))) %>%
  mutate(term = relevel(term, "Overall")) %>%
  mutate(size = case_when(type == "study" ~ .25, type == "summary" ~ 10))

labels <- c("lean_mass" = "Muscle mass",
            "str" = "Strength",
            "func" = "Function")

(  
lean_fp<-df_forest %>%
  filter(!is.na(estimate) & outcome=="lean_mass") %>%
  ggplot(aes(x=estimate, xmin=conf.low, xmax=conf.high, 
             y=term, size=weight, shape=type)) +
  facet_wrap(~outcome, labeller = as_labeller(labels)) +
  geom_stripes(inherit.aes=F, aes(y=term, xmin=-Inf, xmax=Inf), 
               odd = "#ced3d6", even = "#00000000", alpha=.05) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey30") +
  geom_errorbar(size=.3, width=.25) +
  geom_point(aes(color=factor(type))) +  
  geom_text(aes(x=Inf, label=result), hjust = 1, size=2.95, color="black") +
  annotate("text", x=Inf, y=13, label="SMD [95% CI]", hjust=1, vjust=-4.95, size=3) +
  theme_classic() +
  scale_color_manual(values=c("#00A1D5", "black")) +
  scale_size_continuous(range=c(2, 4)) +
  scale_shape_manual(values=c(20, 18)) +
  scale_x_continuous(expand = expansion(mult = c(.05, .4)),
                      breaks = seq(-1, .5, .25)) +
  coord_cartesian(clip="off") +
  guides(color="none", size="none", shape="none") +
  labs(x=NULL, y=NULL) +
    theme(text = element_text(size = 11),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(color="black"),
          plot.margin = unit(c(18, 5.5, 5.5, 5.5), "pt"))
)

(
func_fp<-df_forest %>%
  filter(!is.na(estimate) & outcome=="func") %>%
  ggplot(aes(x=estimate, xmin=conf.low, xmax=conf.high, 
             y=term, size=weight, shape=type)) +
  facet_wrap(~outcome, labeller = as_labeller(labels)) +
  geom_stripes(inherit.aes=F, aes(y=term, xmin=-Inf, xmax=Inf), 
                 odd = "#ced3d6", even = "#00000000", alpha=.05) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey30") +
  geom_errorbar(size=.3, width=.25) +
  geom_point(aes(color=factor(type))) +
  geom_text(aes(x=Inf, label=result), hjust = 1, size=2.95, color="black") +
  theme_classic() +
  scale_color_manual(values=c("#DF8F44", "black")) +
  scale_size_continuous(range=c(2, 4)) +
  coord_cartesian(clip="off") +
  scale_shape_manual(values=c(20, 18)) +
  scale_x_continuous(expand = expansion(mult = c(.05, .375)),
                     breaks = seq(-1, 1, .5)) +
  guides(color="none", size="none", shape="none") +
  labs(x="Standardized Mean Difference (SMD)", y=NULL) +
  theme(text = element_text(size = 11),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(color="black"))
)

(
str_fp<-df_forest %>%
  filter(!is.na(estimate) & outcome=="str") %>%
  ggplot(aes(x=estimate, xmin=conf.low, xmax=conf.high, 
             y=term, size=weight, shape=type)) +
  facet_wrap(~outcome, labeller = as_labeller(labels)) +
  geom_stripes(inherit.aes=F, aes(y=term, xmin=-Inf, xmax=Inf), 
               odd = "#ced3d6", even = "#00000000", alpha=.05) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey30") +
  geom_errorbar(size=.3, width=.25) +
  geom_point(aes(color=factor(type))) +
  geom_text(aes(x=Inf, label=result), hjust = 1, size=2.95, color="black") +
  theme_classic() +
  scale_color_manual(values=c("#374E55", "black")) +
  scale_size_continuous(range=c(2, 4)) +
  coord_cartesian(clip="off") +
  scale_shape_manual(values=c(20, 18)) +
  scale_x_continuous(expand = expansion(mult = c(.05, .4)),
                     breaks = seq(-1, 1.5, .5)) +
  guides(color="none", size="none", shape="none") +
  labs(x=NULL, y=NULL) +
  theme(text = element_text(size = 11),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(color="black"))
  )

fp_all <- (lean_fp / str_fp / func_fp)

fp_all

# ggsave(here("Figures/Figure3_Forestplot.png"), units="in", width=6, height=11, dpi=600, type="cairo")

## Publication bias and small-study effects

# Eggers test

eggers <- df_forest %>% ## linear regression between the ESs normalized by SE
  filter(term != "overall") %>% ## against precision (reciprocal of SE)
  mutate(y=estimate/std.error, x=1/std.error) %>% ## Mind that this code will break if using |> instead of %>%
  lm(y ~ x, data = .) %>%
  tidy(conf.int=T)

est <- format(round(eggers$estimate[1], digits=2), nsmall=2)
est_lb <- format(round(eggers$conf.low[1], digits=2), nsmall=2)
est_ub <- format(round(eggers$conf.high[1], digits=2), nsmall=2)
eggers_p <- format(round(eggers$p.value[1], digits=2), nsmall=2)

eggers_text <- enframe(glue("Egger's test: {est} (95% CI: {est_lb}; {est_ub}), p-value = {eggers_p}"))

# Funnel plot of all outcomes
# Special thanks to John Sakaluk for the code, tutorial here: 
# https://sakaluk.wordpress.com/2016/02/16/7-make-it-pretty-plots-for-meta-analysis/

full.model <- rma.mv(yi = yi_2, 
                     V = vi_2, 
                     slab = authoryear,
                     data = n3data,
                     random = ~ 1 | authoryear/es.id, 
                     test = "t", 
                     method = "REML")

full.model.sum <- summary(full.model)
estimate <- full.model.sum$b[1]    
se <- full.model.sum$se[1]
variance <- as.data.frame(full.model$vi)
se_all <- sqrt(variance)
se.seq <- seq(0, max(se_all), 0.001)
ll95 = estimate-(1.96*se.seq)
ul95 = estimate+(1.96*se.seq)
meanll95 = estimate-(1.96*se)
meanul95 = estimate+(1.96*se)
funnel.full.data = data.frame(ll95, ul95, se.seq, estimate, meanll95, meanul95)

ggplot(aes(x = sqrt(vi_2), y=yi_2), data=n3data) +
  geom_segment(size = .7, aes(x = min(se.seq), xend = max(se.seq), 
                              y = estimate, yend = estimate), linetype = "33")+
  geom_segment(size = .7, aes(x = .034, xend = max(se.seq), 
                              y = 0, yend = 0), linetype = "33", color="red")+
  geom_point(alpha=0.75, size=2, aes(color=outcome_type)) +
  xlab('Standard Error') + ylab('Standardized Mean Difference (SMD)')+ labs(color="Outcome")+
  geom_line(aes(x = se.seq, y = ll95), linetype = "dashed", data = funnel.full.data) +
  geom_line(aes(x = se.seq, y = ul95), linetype = "dashed", data = funnel.full.data) +
  scale_x_reverse()+
  scale_y_continuous(breaks=seq(-1.5,2,0.25))+
  coord_flip()+
  theme_classic()+
  geom_text(inherit.aes=F, data=eggers_text, alpha=.65, size=2.75,
            aes(x=Inf, y=Inf, label=value), hjust=1.96, vjust=-0.4, fontface="italic") +
  scale_color_jama(labels=c("Muscle mass", "Function", "Strength"))

# ggsave(here("Figures/Figure2_Funnelplot.png"), units="in", width=7, height=5, dpi=600, type="cairo")

### Supplementary material

## Correlation sensitivity analysis using a 0.7 value

n3data <- n3data |> 
  mutate(ri_sens = 0.7)

ri_sens_df <- data.frame(prepostcon(n3data$ri_sens, 
                             n3data$meanpre_n3, 
                             n3data$meanpost_n3, 
                             n3data$meanpre_pla, 
                             n3data$meanpost_pla, 
                             n3data$sdpre_n3, 
                             n3data$sdpre_pla,
                             n3data$ntotal_n3, 
                             n3data$ntotal_pla)) %>%
  `colnames<-`(c("yi_sens", "vi_sens", "se_sens", "lowCI_sens", "highCI_sens")) 

n3data <- bind_cols(n3data, ri_sens_df)

n3data_mass <- n3data %>% 
  filter(outcome_type=="body composition")

n3data_str <- n3data %>%
  filter(outcome_type=="strength")

n3data_func <- n3data %>%
  filter(outcome_type=="function")

# MA models

meta_lean_mass_sens <- rma.mv(yi = yi_sens, 
                         V = vi_sens, 
                         slab = authoryear,
                         data = n3data_mass,
                         random = ~ 1 | authoryear/es.id, 
                         test = "t", 
                         method = "REML")

meta_str_sens <- rma.mv(yi_sens, 
                   vi_sens,
                   slab=authoryear,
                   random= ~ 1 | authoryear/es.id,
                   method="REML",
                   test="t",
                   data=n3data_str)

meta_func_sens <- rma.mv(yi_sens, 
                    vi_sens, 
                    slab = authoryear,
                    data = n3data_func,
                    random = ~ 1 | authoryear/es.id, 
                    test = "t", 
                    method = "REML")

sens_cor_df <- reduce(list(tidy(meta_lean_mass_sens, conf.int=TRUE),
            tidy(meta_str_sens, conf.int=TRUE),
            tidy(meta_func_sens, conf.int=TRUE)), full_join) |> 
  mutate(meta_type = rep(c("Muscle mass", "Strength", "Function")),
         meta_type = fct_relevel(meta_type, "Muscle mass", "Strength", "Function"))

# making a table with results 

sens_cor_tbl <- sens_cor_df |> select(meta_type, -term, -type, 
                      estimate, std.error, statistic, p.value, 
                      conf.low, conf.high) |> 
  mutate_if(is.numeric, round, digits=3) |> 
  arrange(meta_type)

# write_xlsx(sens_cor_tbl, here("Tables/sens_cor_table.xlsx"))

## Supplementary table for I^2 values

I2_tbl <- reduce(list(compI2(meta_lean_mass), # Takes a while to compute everything!
            compI2(meta_lean_mass_training),
            compI2(meta_str),
            compI2(meta_str_dose),
            compI2(meta_str_training),
            compI2(meta_str_age),
            compI2(meta_func),
            compI2(meta_func_dose),
            compI2(meta_func_training)), full_join) |> 
  mutate(type=c("lean_mass", "lean_mass", "lean_mass_training", "lean_mass_training",
                "str", "str", "str_dose", "str_dose", "str_training", "str_training", 
                "str_age", "str_age", "func", "func", 
                "func_dose", "func_dose", "func_training", "func_training")) |> 
  mutate_if(is.numeric, round, 2) |> 
  pivot_wider(names_from = level, values_from = c(I2, ci.lb, ci.ub)) |> 
  select(type, I2_outcome, ci.lb_outcome, ci.ub_outcome, I2_study, ci.lb_study, ci.ub_study)

# write_xlsx(I2_tbl, here("Tables/I2_table.xlsx"))

# ------------------------------------------------------------------------------------------#