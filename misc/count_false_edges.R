# glm-funk simulations
library(glmfunk)
library(glmnet)
library(MASS)
library(igraph)
library(Matrix)
library(scalreg)
library(lmtest)
library(MASS)
library(sandwich)

p = 300
s = 20
design = "uninform"

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

# Feature graph
cs = s/2
block = rbind(c(0, rep(1, cs - 1)), cbind(1, matrix(0, cs - 1, cs - 1)))
blist = replicate(p/cs, block, FALSE)
adj_p = as.matrix(bdiag(blist))
adj_p_true = adj_p

counts = numeric(500)
for (i in 1:500) {
  
  set.seed(i)

  if (design == "uninform") {
    for (ix in 1:(p/cs - 1)) {
      replace = adj_p[1:cs + (ix-1)*cs, (1 + (ix)*cs):p]
      adj_p[1:cs + (ix-1)*cs, (1 + (ix)*cs):p] = rbinom(nrow(replace)*ncol(replace), 1, 0.002)
    }
    adj_p[lower.tri(adj_p)] = t(adj_p)[lower.tri(adj_p)]
  }

  g_true = graph_from_adjacency_matrix(adj_p_true, mode="undirected")
  g = graph_from_adjacency_matrix(adj_p, mode="undirected")

  counts[i] = length(E(g)) - length(E(g_true)) 
}

summary(counts) / 270
