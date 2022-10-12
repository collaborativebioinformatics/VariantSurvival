df <-de2

#survival analysis
df <- de2 %>%
  mutate(significant= ifelse(significant=="Yes", 1, 0))

#df <- df %>% drop_na(Time_to_death_or_last_followup_days)  #why do not work ? # need to remove NAs


#head(df[, c("significant_variant", "Time_to_death_or_last_followup_days", "Phenotype")])

#survdiff(Surv(Time_to_death_or_last_followup_days, Phenotype) ~ significant_variant, data = df)
#
output$plot4 <- renderPlot({
  s <- survfit(Surv(Time_to_death_or_last_followup_days,Phenotype)~ significant, data = df)
  s
  g <- ggsurvplot(
    s,
    conf.int = TRUE,
    data = df,
    risk.table = TRUE
  )
  g
})
# sv count
output$plot5 <- renderPlot({
  s2 <- survfit(Surv(Time_to_death_or_last_followup_days, significant)~ Phenotype, data = df)
  
  g2 <- ggsurvplot_facet(
    s2,
    conf.int = TRUE,
    data = df,
    facet.by = "SV_count",
    nrow = 1,
  )
  g2

})

#
output$plot6 <- renderTable({
  x <- coxph(Surv(Time_to_death_or_last_followup_days,Phenotype) ~ significant, data = df) %>% 
    tbl_regression(exp = TRUE) 
  x
})


