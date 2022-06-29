# Performing permutation fdr test for deletion analysis
# installing packages
install.packages("tidyverse")
install.packages("data.table")
install.packages("magrittr")
install.packages("Amelia")
install.packages("broom")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("janitor")
install.packages("survival")
install.packages("survminer")
install.packages("gtsummary")

library(tidyverse)
library(data.table)
library(Amelia)
library(broom)
library(dplyr)
library(ggplot2)
library(janitor)
library(survival)
library(survminer)
library(gtsummary)
library(magrittr)
library(purrr)

set.seed(1234)

# Randomizing the rejection status with sample 
# and performing permutation test
permutation.test <- map(1:10000, function(y) {
  print(y)
  DATAFRAME <- R_covariates_deletions_collision
  DATAFRAME$Rejection <- sample(R_covariates_deletions_collision$Rejection)
  SURV_COEF <- map(colnames(DATAFRAME)[207:246], function(x) {
    DATA <- select(DATAFRAME, x, Time_to_event_months, Rejection, 
                   R_Age, D_Age, R_Gender, D_Gender, Cold_ischemia, 
                   Eplets_total_HLAI, Eplets_total_HLAII)
    COX <- coxph(Surv(Time_to_event_months, Rejection) ~., data = DATA)
    COEF <- summary(COX) %>% coef() %>% .[1,1]
    return(COEF)
  }) %>% unlist %>% max
  return(SURV_COEF)
}) %>% unlist

sample(R_covariates_deletions_collision$Rejection)
replicate(2, sample(R_covariates_deletions_collision$Rejection))

write.table(PERMUTATION_TEST, 
            file = "/results/PERMUTATION.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
