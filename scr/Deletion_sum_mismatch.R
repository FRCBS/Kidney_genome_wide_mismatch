###############################################################################
# 210525
# Deletion sum between donor and recipient
# (genomic collision)

# CONTINUING WITH THE COVARIATE FILES FROM PREVIOUS script:
# 'Survival_analysis_for_deletions_mismatch.R'
# where the mismatch has already been counted
###############################################################################

# Importing the covariate files 
R_dos_pheno_dels_collision <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/R_covariates_deletions_collision.txt")
D_dos_pheno_dels_collision <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/D_covariates_deletions_collision.txt")

# The collision value for mismatch is 1, other cases are 0
# recipient 0 + donor 0 = 0
# recipient 0 + donor 1 = 0
# recipient 1 + donor 0 = 1
# recipient 1 + donor 1 = 0

# Calculate the sum of deletions/mismatches (40 SNPs) per individual
R_dos_pheno_dels_collision$Del_sum_collision <- c(Del_sum_collision = rowSums(R_dos_pheno_dels_collision[207:246]))
D_dos_pheno_dels_collision$Del_sum_collision <- c(Del_sum_collision = rowSums(D_dos_pheno_dels_collision[207:246]))

# Write out the tables
write.table(R_dos_pheno_dels_collision, 
            file = "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/R_covariates_deletions_sum_collision.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(D_dos_pheno_dels_collision, 
            file = "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/D_covariates_deletions_sum_collision.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################
### Logistic regression analyses for deletion mismatch sum
### and acute rejection endpoint 

# Logistic regression for deletion mismatch, unadjusted
glm_deletionsumcol <- glm(Rejection ~ Del_sum_collision,
                          family = binomial(link = "logit"), data = R_dos_pheno_dels_collision)
summary(glm_deletionsumcol)

# Logistic regression for deletion mismatch, unadjusted
glm_delsumcol_adjusted <- glm(Rejection ~ Del_sum_collision + R_Gender + D_Gender + R_Age + D_Age + Cold_ischemia + 
                                Eplets_total_HLAI + Eplets_total_HLAII,
                              family = binomial(link = "logit"), data = R_dos_pheno_dels_collision)
summary(glm_delsumcol_adjusted)

write.table(tidy(glm_delsumcol_adjusted), 
            "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/glm_delsumcol_adjusted",
            sep = "\t", quote = F, row.names = F)

# OR and 95% CI
exp(cbind(OR = coef(glm_delsumcol_adjusted), confint(glm_delsumcol_adjusted)))

###############################################################################
### The survival analysis: the association of deletion sum mismatch to acute  
### rejection event

# Survival analysis of deletion sum mismatch
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ Del_sum_collision, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Kapplan-Meier plot with risk table for the manuscript
del_sum_col <- ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ Del_sum_collision, data = R_dos_pheno_dels_collision), 
  xlab = "Time to event (months)", 
  ylab = "Probability of graft survival without rejection",
  font.x = c(face = "bold"),
  font.y = c(face = "bold"),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Del sum 0", "Del sum 1", "Del sum 2", "Del sum 3", "Del sum 4", "Del sum 5", 
                  "Del sum 6", "Del sum 7", "Del sum 8", "Del sum 9"),
  legend = c(0.1,0.25))

del_sum_col$plot <- del_sum_col$plot+ 
  ggplot2::annotate("text", 
                    x = 100, y = 0.25, # x and y coordinates of the text
                    label = "HR 1.01, 95% CI 0.93-1.10 \n P-value 0.8",
                    size = 5)
del_sum_col

# HR and CI 95% for unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ Del_sum_collision, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# HR and CI 95% for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ Del_sum_collision + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + 
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Cox proportional hazards model for adjusted data
cox_delsum_collision <- coxph(Surv(Time_to_event_months, Rejection) ~ Del_sum_collision + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + 
                                Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision)
summary(cox_delsum_collision)
