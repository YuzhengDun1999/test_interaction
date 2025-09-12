n = 100; p = 850; main_prop = 0.05
methods = c("OLS", "Ridge", "Lasso OLS minMSE")
result <- tibble(n = integer(), p = integer(), sim = integer(), main_prop = numeric(), method = character(), pval = numeric())

for (sim in 1:5000) {
  set.seed(n + p + sim + main_prop * 10)
  E = as.matrix(data.frame(continuous = rnorm(n), binary = rbinom(n, 2, 0.5)))
  Z = matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_E = c(5, 10)
  beta_Z = rep(0, p)
  beta_Z[1:as.integer(round(p * main_prop))] = 40
  V = NULL
  for(hhh in 1:ncol(E)){
    V = as.matrix(cbind(V, drop(E[, hhh]) * Z))
  }
  beta_V = rep(0, ncol(V))
  #beta_V[1] = 1
  Y = as.matrix(E) %*% beta_E + Z %*% beta_Z + V %*% beta_V+ rnorm(n)
  if (n > p){
    OLS_p = iSKATtest::GESAT(Z, Y, E, ols = T, is_check_genotype = FALSE)$pvalue
  } else{
    OLS_p = NA
  }
  result = bind_rows(
    result,
    tibble(
      n = n, p = p, sim = sim, main_prop = main_prop, method = methods,
      pval = c(OLS_p, # OLS
               iSKAT::GESAT(Z, Y, E, lower = 1e-1, upper = 16, nintervals = 250, is_check_genotype = FALSE)$pvalue, # ridge
               iSKATtest::GESAT(Z, Y, E, lasso.ols = T, lower = 1e-1, upper = 16, nintervals = 250, is_check_genotype = FALSE)$pvalue, #lasso OLS
               iSKATtest::GESAT(Z, Y, E, lasso.ols = T, lower = 1e-1, upper = 16, nintervals = 250, lasso.criterion = "lambda.1se", is_check_genotype = FALSE)$pvalue, # lasso OLS
               iSKATtest::GESAT(Z, Y, E, lasso.select = T, lower = 1e-1, upper = 16, nintervals = 250, is_check_genotype = FALSE)$pvalue, # lasso ridge select by minMSE
               iSKATtest::GESAT(Z, Y, E, lasso.select = T, lower = 1e-1, upper = 16, nintervals = 250, lasso.criterion = "lambda.1se", is_check_genotype = FALSE)$pvalue
      )
    )
  )
}

result %>% group_by(n, p, method) %>% summarize(reject_rate = mean(pval < 0.05), .groups = "drop")

