library(INLA)
library(dplyr)

inla_bym2_poisson <- function(y, X, adj) {
  
  n <- length(y)
  idarea = 1:n
  g_inla <- inla.read.graph(adj)
  data <- as.data.frame(cbind(y, X))
  colnames(data)[1] = "y"
  
  formula <- as.formula(paste0("y ~ 0 + ", 
                               paste(colnames(data)[-1], collapse = " + "), 
                               "+ f(idarea, model = 'bym2', graph = g_inla, hyper = prior)"))
  
  prior <- list(
    prec = list(
      prior = "logtnormal",
      param = c(0, 1)),
    phi = list(
      prior = "logitbeta",
      param = c(0.5, 0.5))
  )
  
  fit <- inla(formula, family = "poisson", data = data)
  
  beta <- fit$summary.fixed %>% pull(`0.5quant`)
  
  sig <- fit$summary.fixed %>%
    transmute(sig = !(`0.025quant` <= 0 & 0 <= `0.975quant`)) %>%
    pull(sig)
  
  alpha <- fit$summary.random$idarea[1:n,] %>% pull(`0.5quant`)
  
  cis <- fit$summary.fixed %>% 
    transmute(ci_lower = `0.025quant`, 
              ci_upper = `0.975quant`)
  
  output <- list(beta = beta,
                 sig = sig,
                 alpha = alpha,
                 cis = cis)
  
  return(output)
}


bym2_predict <- function(bym2_out, Xtest, Ln_full, train_ind) {
  
  nfull = nrow(Ln_full)
  test_ind = setdiff(1:nfull, train_ind)
  Ln_22 = Ln_full[test_ind, test_ind]
  Ln_21 = Ln_full[test_ind, train_ind]
  alphatrain = as.matrix(bym2_out$alpha, ncol = 1)
  alphatest = as.numeric(-ginv(Ln_22) %*% Ln_21 %*% alphatrain)
  testpred = exp(alphatest + as.numeric(Xtest %*% bym2_out$beta))
  
  return(testpred)
}
