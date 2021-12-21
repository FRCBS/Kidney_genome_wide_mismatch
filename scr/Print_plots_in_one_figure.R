###############################################################################
# Printing the survival plots and arraging them in one figure
# First plotting them separately

# Plottng the kidney quartiles
kidney_quartiles_plot <- ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_kidney), 
  title = "A",
  xlab = "Time to event (months)", 
  ylab = "Rejection-free survival",
  font.x = c(face = "bold"),
  font.y = c(face = "bold"),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Q1: MM sum 509-584", "Q2: MM sum 585-604", "Q3: MM sum 605-626", "Q4: MM sum 627-757"),
  legend = c(0.20,0.50))
kidney_quartiles_plot

# The same plot with additional details written in the picture (HR and 95% CI)
kidney_quartiles_plot$plot <- kidney_quartiles_plot$plot+ 
  ggplot2::annotate("text", 
                    x = 40, y = 0.15, # x and y coordinates of the text
                    label = "HR 1.15, 95% CI 1.01-1.30, P-value 0.029",
                    size = 4)
kidney_quartiles_plot

# Print out the ggsurvplot
ggsave("/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/kidney_plot", plot = print(kidney_quartiles_plot),
       device = "pdf", dpi = 1000)

# Plotting the rs7542235 deletion
rs7542235_del_col <- ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs7542235_col, data = R_dos_pheno_dels_collision), 
  title = "B",
  xlab = "Time to event (months)", 
  ylab = "Rejection-free survival",
  font.x = c(face = "bold"),
  font.y = c(face = "bold"),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("no mismatch", "mismatch"),
  legend = c(0.20,0.35))
rs7542235_del_col

rs7542235_del_col$plot <- rs7542235_del_col$plot+ 
  ggplot2::annotate("text", 
                    x = 40, y = 0.15, # x and y coordinates of the text
                    label = "HR 3.10, 95% CI 1.53-6.29, P-value 0.002",
                    size = 4)
rs7542235_del_col

# Print out the ggsurvplot
ggsave("/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/rs7542235_plot", plot = print(rs7542235_del_col),
       device = "pdf", dpi = 1000)

###############################################################################
# Rearranging the plots in one figure

splots <- list()
splots[[1]] <- kidney_quartiles_plot
splots[[2]] <- rs7542235_del_col

# In two columns
plots <- arrange_ggsurvplots(splots, print = FALSE, ncol = 2, nrow = 1, risk.table.height = 0.25)
ggsave("/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/plots.pdf", plots,
       width = 13, height = 8, units = "in", dpi = 1000)

# In two rows
plots2 <- arrange_ggsurvplots(splots, print = FALSE, ncol = 1, nrow = 2, risk.table.height = 0.25)
ggsave("/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/plots2.pdf", plots2,
       width = 8, height = 13, units = "in", dpi = 1000)

