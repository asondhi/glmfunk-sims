# glm-funk simulations
library(glmfunk)
library(glmnet)
library(MASS)
library(igraph)
library(Matrix)
library(scalreg)
library(lmtest)
library(sandwich)

library(doParallel)
registerDoParallel(cores=detectCores())
print(detectCores())

source("bym2_inla_fit.R")

library(tidyverse)

get_cov_prob <- function(cis, beta) {
  
  colnames(cis) = c("ci_lower", "ci_upper")
  
  cis %>%
    mutate(beta = beta,
           covered = 1*(beta <= ci_upper & beta >= ci_lower)) %>%
    pull(covered) %>%
    mean
  
}

# Parse arguments: p, s, rho, design, iter, method
# s = # nonzero features, should be even number
# rho = absolute value of effect size for betas
args = commandArgs(TRUE)
print(args)
for (i in 1:length(args)) {
  eval(parse(text = args[i]))
}

if (design == 1) {
  design = "inform" # fully informative nets
} else if (design == 2) {
  design = "uninform" # misspec sim
} 

if (method == 1) {
  method = "funkl2"
} else if (method == 2) {
  method = "funkl1"
} else if (method == 3) {
  method = "rnc-lasso"
} else if (method == 4) {
  method = "lasso"
} else if (method == 5) {
  method = "bym2"
} else if (method == 6) {
  method = "glm"
}

# load King County graph Laplacian
load("KC_graph.RData")
big_adj = as.matrix(bdiag(list(KC_adj, KC_adj, KC_adj, KC_adj, KC_adj, KC_adj,
                               KC_adj, KC_adj, KC_adj, KC_adj, KC_adj, KC_adj)))
g_big = graph_from_adjacency_matrix(big_adj)
Ln = diag(apply(big_adj, 1, sum)) - big_adj
adj_n = big_adj
colnames(adj_n) = NULL
rownames(adj_n) = NULL
n = nrow(adj_n) / 2

# Subject graph => alphas
Q <- 10 * Ln # scaled Laplacian
eigen_Q <- eigen(Q)
evals <- eigen_Q$values 
evals <- evals[evals > 1e-08]
z <- rnorm(length(evals), 0, 1 / evals)
E <- eigen_Q$vectors

do_sim_iter <- function(seed) {
  
  set.seed(seed)
  
  if (design == "inform") {
    alpha <- as.numeric(E[, 1:length(evals)] %*% z)
  } else {
    alpha = rnorm(2*n, mean = 0, sd = 0.24)
  }
  
  # Feature graph
  cs = s/2
  block = rbind(c(0, rep(1, cs - 1)), cbind(1, matrix(0, cs - 1, cs - 1)))
  blist = replicate(p/cs, block, FALSE)
  adj_p = as.matrix(bdiag(blist))
  if (design == "uninform") {
    for (ix in 1:(p/cs - 1)) {
      replace = adj_p[1:cs + (ix-1)*cs, (1 + (ix)*cs):p]
      adj_p[1:cs + (ix-1)*cs, (1 + (ix)*cs):p] = rbinom(nrow(replace)*ncol(replace), 1, 0.002)
    }
    adj_p[lower.tri(adj_p)] = t(adj_p)[lower.tri(adj_p)]
  }

  if (method == "funkl1" || method == "gracel1") {  
    g = graph_from_adjacency_matrix(adj_p, mode="undirected")
  	el = get.edgelist(g)
  	Lp = matrix(0, nrow = nrow(el), ncol = p)
  	for (i in 1:nrow(el)) {
  	  edge = el[i,]
  	  Lp[i, edge[1]] = 1
  	  Lp[i, edge[2]] = -1
  	}
  } else {
    Lp = laplacian(adj_p)
  }
  
  # Generate data
  beta = c(rep(rho, cs), rep(-rho, cs), rep(0, p - 2*cs))
  X1 = t(replicate(n, gen.one(cs)))
  X1 = scale(X1, scale = apply(X1, 2, function(x) max(x) - min(x)))
  X2 = t(replicate(n, gen.one(cs)))
  X2 = scale(X2, scale = apply(X2, 2, function(x) max(x) - min(x)))
  ind1 = sample(1:(2*n), size = n)
  ind2 = setdiff(1:(2*n), ind1)
  alpha1 = alpha[ind1]
  alpha2 = alpha[ind2]

  y1 = rpois(n, lambda = exp(alpha1 + as.numeric(X1%*%beta)))
  y2 = rpois(n, lambda = exp(alpha2 + as.numeric(X2%*%beta)))

  Ln1 = Ln[ind1, ind1]
  Ln2 = Ln[ind2, ind2]

  # Network information
  gn.info = sqrt(sum(((Ln1 %*% alpha1)^2)))
  gp.info = sum(abs(Lp %*% beta))

  if (p > n) {
    lambda_grid = sort(exp(seq(-5, 2, length = 50)), decreasing = TRUE)
    gamma_n_grid = sort(exp(seq(-5, 2, length = 50)), decreasing = TRUE)
    gamma_p_grid = sort(exp(seq(-5, 2, length = 50)), decreasing = TRUE)
  } else {
    lambda_grid = c(exp(-5))
    gamma_n_grid = c(exp(-5))
    gamma_p_grid = sort(exp(seq(-5, 2, length = 30)), decreasing = TRUE)
  }
  
  if (method == "funkl2") {
    param_list = list(gamma_n_grid, gamma_p_grid, lambda_grid)
    # CV
    cv.out = funkl2.cv.coorddesc(y1, X1, Ln1, Lp, param_list, model = "poisson", verb.cv = TRUE)
    cv.params = cv.out[[1]]
    cv.minerr = cv.out[[2]]
    cv.niter = cv.out[[3]]
    # Fit
    est = funkl2.fit(y1, X1, Ln1, Lp, cv.params, model = "poisson")
    # Inference
    Ts = debias.poisson.funk(y1, X1, c(est$alphahat, est$betahat), mu = NULL)
    t1er = mean(abs(Ts[(s+1):p]) > 1.96)
    pow = mean(abs(Ts[1:s]) > 1.96)
    
    cis = debias.poisson.funk.ci(y1, X1, c(est$alphahat, est$betahat), mu = NULL)
    cov_prob = get_cov_prob(cis, beta)
    
    # Test set error
    testpred = funkl2.predict(est, X2, Ln, ind1)
    test.err = poisson_loss(y2, testpred)
  }
  
  if (method == "funkl1") {
    param_list = list(gamma_n_grid, gamma_p_grid, lambda_grid)
    # CV
    cv.out = funkl1.cv.coorddesc(y1, X1, Ln1, Lp, param_list, model = "poisson", verb.cv = TRUE)
    cv.params = cv.out[[1]]
    cv.minerr = cv.out[[2]]
    cv.niter = cv.out[[3]]
    # Fit
    est = funkl1.fit(y1, X1, Ln1, Lp, cv.params, model = "poisson")
    # Inference
    Ts = debias.poisson.funk(y1, X1, c(est$alphahat, est$betahat), mu = NULL)
    t1er = mean(abs(Ts[(s+1):p]) > 1.96)
    pow = mean(abs(Ts[1:s]) > 1.96)

    cis = debias.poisson.funk.ci(y1, X1, c(est$alphahat, est$betahat), mu = NULL)
    cov_prob = get_cov_prob(cis, beta)
    
    # Test set error
    testpred = funkl1.predict(est, X2, Ln, ind1)
    test.err = poisson_loss(y2, testpred)
  }

  if (method == "rnc-lasso") {
    param_list = list(gamma_n_grid, lambda_grid)
    # CV
    cv.out = rnc.lasso.cv.coorddesc(y1, X1, Ln1, param_list, model = "poisson", verb.cv = TRUE)
    cv.params = cv.out[[1]]
    cv.minerr = cv.out[[2]]
    cv.niter = NA
    # Fit
    est = rnc.lasso.fit(y1, X1, Ln1, cv.params, model = "poisson")
    # Inference
    Ts = debias.poisson.funk(y1, X1, c(est$alphahat, est$betahat), mu = NULL)
    t1er = mean(abs(Ts[(s+1):p]) > 1.96)
    pow = mean(abs(Ts[1:s]) > 1.96)
    
    cis = debias.poisson.funk.ci(y1, X1, c(est$alphahat, est$betahat), mu = NULL)
    cov_prob = get_cov_prob(cis, beta)
    
    # Test set error
    testpred = rnc.lasso.predict(est, X2, Ln, ind1)
    test.err = poisson_loss(y2, testpred)
    cv.params = c(cv.params[1], NA, cv.params[2])
  }
  
  if (method == "lasso") {
    m = cv.glmnet(x = X1, y = y1, family = "poisson", lambda = lambda_grid)
    # Fit
    est = coef(m, s = "lambda.min")
    intcpt = as.numeric(est)[1]
    est = as.numeric(est)[-1]
    # Inference
    Ts = debias.poisson.lasso(y1, X1, est)
    t1er = mean(abs(Ts[(s+1):p]) > 1.96)
    pow = mean(abs(Ts[1:s]) > 1.96)
    cv.minerr = min(m$cvm)
    
    cis = debias.poisson.lasso.ci(y1, X1, est)
    cov_prob = get_cov_prob(cis, beta)
    
    # Test set error
    testpred = exp(intcpt + as.numeric(X2%*%est))
    test.err = poisson_loss(y2, testpred)
    cv.params = c(NA, NA, m$lambda.min)
    cv.niter = NA
  }
  
  if (method == "bym2") {
    # Fit
    adj_n1 = adj_n[ind1, ind1]
    m = inla_bym2_poisson(y = y1, X = X1, adj = as(adj_n1, "sparseMatrix"))
    # Inference
    t1er = mean(m$sig[(s+1):p])
    pow = mean(m$sig[1:s])
    
    cis = m$cis
    cov_prob = get_cov_prob(cis, beta)
    
    # Test set error
    testpred = bym2_predict(m, X2, Ln, ind1)
    test.err = poisson_loss(y2, testpred)
    cv.params = c(NA, NA, NA)
    cv.niter = NA
  }
  
  if (method == "glm") {
    glm_data <- data.frame(cbind(y1, X1))
    colnames(glm_data)[1] = "y"
    formula <- as.formula(paste0("y ~ ", 
                                 paste(colnames(glm_data)[-1], collapse = " + ")))
    # Fit
    m = glm(formula, data = glm_data, family = "poisson")
    # Inference
    Ts = broom::tidy(coeftest(m, vcov = sandwich))$statistic[-1]
    t1er = mean(abs(Ts[(s+1):p]) > 1.96)
    pow = mean(abs(Ts[1:s]) > 1.96)
    
    cis = broom::tidy(m, conf.int = TRUE, exponentiate = FALSE) %>%
      transmute(ci_lower = `conf.low`,
                ci_upper = `conf.high`)
    cov_prob = get_cov_prob(cis[-1,], beta)
    
    # Test set error
    testpred = exp(as.numeric(X2 %*% coef(m)[-1]) + coef(m)[1])
    test.err = poisson_loss(y2, testpred)
    cv.params = c(NA, NA, NA)
    cv.niter = NA
  }
  
  c(t1er, pow, cov_prob, cv.params[1], cv.params[2], cv.params[3], test.err, gn.info, gp.info, cv.niter)
}


results = foreach(i=1:iter, .combine='rbind') %dopar% {

  res <- tryCatch(do_sim_iter(seed = i), error = function(e) e)
  
  if (inherits(res, "error")) {
    NULL
  } else {
    res
  }
}


colnames(results) = c("t1er", "power", "cov_prob", "gamma_n", "gamma_p", "lambda", "test_err", "gn_info", "gp_info", "cv_niter")

results = as.data.frame(results)
filename = paste0(method, "_", design, "_p", p, "_s", s, "_rho", rho, "_iter", iter)
filename = paste0(filename, ".RData")
save(results, file = filename)

