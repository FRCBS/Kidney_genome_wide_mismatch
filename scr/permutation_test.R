# Performing permutation fdr test for deletion analysis

cox <- coxph(Surv(Time_to_event_months, Rejection) ~ rs6943474_col + R_Age + D_Age + 
        R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_deletions_collision)

summary(cox) %>% coef() %>% .[1,1]

set.seed(1234)

function()
  
map(colnames(R_covariates_deletions_collision)[207:246], function(x) {
  DATA <- select(R_covariates_deletions_collision, x, Time_to_event_months, Rejection, 
                 R_Age, D_Age, R_Gender, D_Gender, Cold_ischemia, 
                 Eplets_total_HLAI, Eplets_total_HLAII)
  COX <- coxph(Surv(Time_to_event_months, Rejection) ~., data = DATA)
  COEF <- summary(COX) %>% coef() %>% .[1,1]
  return(COEF)
}) %>% unlist %>% max


# Randomize the rejection status, 
# Permutation test
permutation.test <- map(replicate(2, sample(R_covariates_deletions_collision[ ,3])), function(y) {
  DATAFRAME <- R_covariates_deletions_collision
  DATAFRAME[ ,3] <- y
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

# sample doesn't work with indexing!! -> works with dollar