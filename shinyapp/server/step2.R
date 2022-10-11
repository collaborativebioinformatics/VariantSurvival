de2
library(tidyverse)
library(janitor)
library(tidyquant)
library(patchwork)
library(survival)
library(survminer)
#survival analysis
de2
library(tidyverse)
library(janitor)
library(tidyquant)
library(patchwork)
library(survival)
library(survminer)
#survival analysis
df <- de2 %>%
  mutate(significant_variant= ifelse(significant_variant=="Yes", 1, 0))

df <- df %>% drop_na(Time_to_death_or_last_followup_days)  #why do not work ? # need to remove NAs


head(df[, c("significant_variant", "Time_to_death_or_last_followup_days", "Phenotype")])

survdiff(Surv(Time_to_death_or_last_followup_days, Phenotype) ~ significant_variant, data = df)
#
c <- coxph(Surv(Time_to_death_or_last_followup_days,Phenotype)~ significant_variant, data = df)

#plot
coxph(Surv(Time_to_death_or_last_followup_days,Phenotype) ~ significant_variant, data = df) %>% 
  tbl_regression(exp = TRUE) 

#
mv_fit <- coxph(Surv(Time_to_death_or_last_followup_days,Phenotype) ~ significant_variant, data = df)
cz <- cox.zph(mv_fit)
print(cz)
plot(cz)