####### n vary from 50 to 400, p vary from 10 to 40
####### 5% genomic markers have true large effect
set.seed(2025)
n = 50; p = 20; main_prop = 0.05
E = as.matrix(data.frame(continuous = rnorm(n), binary = rbinom(n, 2, 0.5)))
Z = matrix(rnorm(n * p), nrow = n, ncol = p)
beta_E = c(1, 1)
beta_Z = rep(0, p)
beta_Z[1:as.integer(round(p * main_prop))] = rnorm(p * main_prop)
Y = as.matrix(E) %*% beta_E + Z %*% beta_Z + rnorm(n)
library(iSKATtest)

########### OLS
iSKATtest::GESAT(Z, Y, E, ols = T, is_check_genotype = FALSE)

########### Ridge
iSKATtest::GESAT(Z, Y, E, is_check_genotype = FALSE)

########### Lasso OLS
iSKATtest::GESAT(Z, Y, E, lasso.ols = T, is_check_genotype = FALSE)
iSKATtest::GESAT(Z, Y, E, lasso.ols = T, lasso.criterion = "lambda.1se", is_check_genotype = FALSE)

########## Lasso ridge
iSKATtest::GESAT(Z, Y, E, lasso.select = T, is_check_genotype = FALSE)
iSKATtest::GESAT(Z, Y, E, lasso.select = T, lasso.criterion = "lambda.1se", is_check_genotype = FALSE)

library(dplyr)
results <- expand.grid(n = c(50, 100, 200, 300, 400), p = c(10, 20, 30, 40), sim = 1:5000, 
                       method = c("OLS", "Ridge", "Lasso OLS minMSE", "Lasso OLS 1se", "Lasso ridge minMSE", "Lasso ridge 1se")) %>% mutate(pval = NA_real_) 



n_values = c(50, 100, 200, 300, 400)
p_values = c(10, 20, 30, 40)
main_props = c(0.1, 0.2, 0.3)
n_sims = 5000
methods = c("OLS", "Ridge", "Lasso OLS minMSE", "Lasso OLS 1se", "Lasso ridge minMSE", "Lasso ridge 1se")

library(dplyr)
result <- tibble(n = integer(), p = integer(), sim = integer(), main_prop = numeric(), method = character(), pval = numeric())

for (n in n_values) {
  for (p in p_values) {
    for (main_prop in main_props) {
      for (sim in 1:n_sims) {
        set.seed(n + p + sim + main_prop * 10)
        E = as.matrix(data.frame(continuous = rnorm(n), binary = rbinom(n, 2, 0.5)))
        Z = matrix(rnorm(n * p), nrow = n, ncol = p)
        beta_E = c(5, 10)
        beta_Z = rep(0, p)
        #beta_Z[1:as.integer(round(p * main_prop))] = sqrt(2) * rnorm(p * main_prop)
        beta_Z[1:as.integer(round(p * main_prop))] = 40
        Y = as.matrix(E) %*% beta_E + Z %*% beta_Z + rnorm(n)
        if (n > p){
          OLS_p = iSKATtestSmallsample::GESAT(Z, Y, E, ols = T, is_check_genotype = FALSE)$pvalue
        } else{
          OLS_p = NA
        }
        result = bind_rows(
          result,
          tibble(
            n = n, p = p, sim = sim, main_prop = main_prop, method = methods,
            pval = c(OLS_p, # OLS
                     iSKATtestSmallsample::GESAT(Z, Y, E, lower = 1e-2, upper = 16, nintervals = 250, is_check_genotype = FALSE)$pvalue, # ridge
                     iSKATtestSmallsample::GESAT(Z, Y, E, lower = 1e-2, upper = 16, nintervals = 250, lasso.ols = T, is_check_genotype = FALSE)$pvalue, #lasso ridge
                     iSKATtestSmallsample::GESAT(Z, Y, E, lower = 1e-2, upper = 16, nintervals = 250, lasso.ols = T, lasso.criterion = "lambda.1se", is_check_genotype = FALSE)$pvalue, # lasso ridge select by 1se
                     iSKATtestSmallsample::GESAT(Z, Y, E, lower = 1e-2, upper = 16, nintervals = 250, lasso.select = T, is_check_genotype = FALSE)$pvalue, # lasso ridge select by minMSE
                     iSKATtestSmallsample::GESAT(Z, Y, E, lower = 1e-2, upper = 16, nintervals = 250, lasso.select = T, lasso.criterion = "lambda.1se", is_check_genotype = FALSE)$pvalue # lasso ridge select by 1se
            )
          )
        )
        #saveRDS(result, paste0("~/Documents/test_interaction/simulation_result/n", n, "_p", p, "_prop", main_prop, ".rds"))
      }
    }
  }
}
result %>% group_by(n, p, method) %>% summarize(reject_rate = mean(pval < 0.05), .groups = "drop")
