df <-de2

#survival analysis
df <- de2 %>%
  mutate(variant= ifelse(variant=="Yes", 1, 0))
#
output$plot4 <- renderPlot({
  s <- survfit(Surv(Time_to_death_or_last_followup_days, variant)~ Phenotype, data = df)
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
  s2 <- survfit(Surv(Time_to_death_or_last_followup_days, variant)~ Phenotype, data = df)
  
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
#output$plot6 <- renderTable({
  x <- coxph(Surv(Time_to_death_or_last_followup_days,Phenotype) ~ variant, data = df) %>% 
    tbl_regression(exp = TRUE) 
  t <- as_data_frame(x)
  t
  
#}) 
output$table <- DT::renderDataTable({
  t
}) 
#output$plot6 <- renderTable({
x2 <- coxph(Surv(Time_to_death_or_last_followup_days, variant)~ Phenotype, data = df) %>% 
  tbl_regression(exp = TRUE) 
t2 <- as_data_frame(x2)
t2

#}) 
output$table2 <- DT::renderDataTable({
  t2
}) 
