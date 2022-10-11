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

#plot1

output$plot1 <- renderPlot({
  x <- coxph(Surv(Time_to_death_or_last_followup_days,Phenotype) ~ significant_variant, data = df) %>% 
    tbl_regression(exp = TRUE) 
  x
})
# plot2
output$plot2 <- renderPlot({
  mv_fit <- coxph(Surv(Time_to_death_or_last_followup_days,Phenotype) ~ significant_variant, data = df)
  cz <- cox.zph(mv_fit)
  print(cz)
  p<-plot(cz)
  p
})

s <- survfit(Surv(Time_to_death_or_last_followup_days, significant_variant)~ Phenotype, data = df)
g <- ggsurvplot(
  s,
  conf.int = TRUE,
  data = df,
  risk.table = TRUE
)

s <- survfit(Surv(Time_to_death_or_last_followup_days, significant_variant)~ Phenotype, data = df)

g3 <- ggsurvplot_facet(
  s,
  conf.int = TRUE,
  data = df,
  facet.by = "gene",
  nrow = 1,
)
g3

