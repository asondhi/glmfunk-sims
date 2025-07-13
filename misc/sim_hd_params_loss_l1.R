# glm-funk simulations
library(glmfunk)
library(glmnet)
library(MASS)
library(igraph)
library(Matrix)
library(scalreg)
library(lmtest)
library(sandwich)
library(tidyverse)

p = 300
s = 20
method = "funkl1"
rho = 0.4

# load King County graph Laplacian
load("KC_graph.RData")
Ln = KCL
adj_n = KC_adj
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

set.seed(1)
alpha <- as.numeric(E[, 1:length(evals)] %*% z)

# Feature graph
cs = s/2
block = rbind(c(0, rep(1, cs - 1)), cbind(1, matrix(0, cs - 1, cs - 1)))
blist = replicate(p/cs, block, FALSE)
adj_p = as.matrix(bdiag(blist))

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

gamma_n_grid = sort(exp(seq(-7, -3, length = 20)), decreasing = TRUE)
gamma_p_grid = sort(exp(seq(-5, 0, length = 20)), decreasing = TRUE)
lambda_grid = sort(exp(seq(-8, -3, length = 20)), decreasing = TRUE)

get_loss_at_params <- function(params) {
  
  est = funkl1.fit(y1, X1, Ln1, Lp, params, model = "poisson")
  testpred = funkl1.predict(est, X2, Ln, ind1)
  test.err = poisson_loss(y2, testpred)
  print(params)
  return(test.err)
}

default_gamma_n <- 0.005
default_gamma_p <- 0.1
default_lambda <- 0.01

# 1: vary gamma_n
# 2: vary gamma_p
# 3: vary lambda
grid <- rbind(cbind(gamma_n_grid, default_gamma_p, default_lambda, 1),
              cbind(default_gamma_n, gamma_p_grid, default_lambda, 2),
              cbind(default_gamma_n, default_gamma_p, lambda_grid, 3))

res <- apply(grid, 1, get_loss_at_params)

library(tidyverse)

plot_df <- data.frame(cbind(grid, res)) %>%
  as_tibble() %>%
  transmute(type = V4,
            gamma_n = gamma_n_grid,
            gamma_p = default_gamma_p,
            lambda = default_lambda,
            nll = res)

p1 <- plot_df %>%
  filter(type == 1) %>%
  ggplot(aes(x = gamma_n, y = nll)) +
  geom_line() +
  theme_bw()

p2 <- plot_df %>%
  filter(type == 2) %>%
  ggplot(aes(x = gamma_p, y = nll)) +
  geom_line() +
  theme_bw()

p3 <- plot_df %>%
  filter(type == 3) %>%
  ggplot(aes(x = lambda, y = nll)) +
  geom_line() +
  theme_bw()

ggpubr::ggarrange(p1, p2, p3, ncol = 3)

save.image(file = "nll_plot.RData")


