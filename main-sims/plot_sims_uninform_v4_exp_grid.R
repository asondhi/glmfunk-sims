library(tidyverse)
library(ggpubr)
library(latex2exp)

collect = function(filename) {

  split = strsplit(filename, split = "_")[[1]]
  method = split[1]
  design = split[2]
  p = as.numeric(substring(split[3], 2))
  s = as.numeric(substring(split[4], 2))
  rho = as.numeric(strsplit(substring(split[5], 4), ".RData")[[1]])
  iter = as.numeric(strsplit(substring(split[6], 5), ".RData")[[1]])
  
  load(filename)
  results = data.frame(results)
  
  results <- results %>%
    mutate(cov_prob_signal = (p/s)*cov_prob - (1 - t1er)*(1 - s/p)*(p/s))
  
  if (!("power" %in% colnames(results))) { results$power = results$pow }

  out = data.frame(method = method, design = design, s = s, rho = rho, iter = iter,
             t1er = mean(results$t1er), 
             power = mean(results$power), 
             test_dev = mean(sqrt(results$test_err)),
             pow_se = sd(results$power)/sqrt(iter),
             t1er_se = sd(results$t1er)/sqrt(iter),
             test_se = sd(sqrt(results$test_err))/sqrt(iter),
             cov_prob = mean(results$cov_prob),
             cov_prob_se = sd(results$cov_prob)/sqrt(iter))
  
  return(out)           
}

res = NULL

setwd("simres_uninform_r2r")

for (file in list.files()) {
  res = rbind(res, collect(file))
}

setwd("..")

#####################

# Uninformative networks

cols <- c("funk l1" = "darkgreen", "funk l2" = "green3", "rnc-lasso" = "red4", "lasso" = "cornflowerblue")

p1 = res %>% 
  mutate_if(is.factor, as.character) %>%
  filter(design == "uninform") %>%
  filter(method %in% c("funkl1", "funkl2", "rnc-lasso", "lasso")) %>%
  select(-s, -design, -iter) %>%
  mutate(method = replace(method, method == "funkl2", "funk l2")) %>%
  mutate(method = replace(method, method == "funkl1", "funk l1")) %>%
  ggplot(aes(x = rho, y = power, color = method)) +
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin=power-pow_se, ymax=power+pow_se), width=0.01) +
  labs(y = "Power", x = TeX(r"( $\rho$ )")) +
  scale_color_manual(values = cols) +
  theme_bw()

p2 = res %>% 
  mutate_if(is.factor, as.character) %>%
  filter(design == "uninform") %>%
  filter(method %in% c("funkl1", "funkl2", "rnc-lasso", "lasso")) %>%
  select(-s, -design, -iter) %>%
  mutate(method = replace(method, method == "funkl2", "funk l2")) %>%
  mutate(method = replace(method, method == "funkl1", "funk l1")) %>%
  ggplot(aes(x = rho, y = t1er, color = method)) +
  geom_line() + 
  geom_point() +
  geom_errorbar(aes(ymin=t1er-t1er_se, ymax=t1er+t1er_se), width=0.01) +
  geom_hline(yintercept = 0.05) + 
  labs(y = "Type I error rate", x = TeX(r"( $\rho$ )")) +
  scale_color_manual(values = cols) +
  theme_bw()

p3 = res %>% 
  mutate_if(is.factor, as.character) %>%
  filter(design == "uninform") %>%
  filter(method %in% c("funkl1", "funkl2", "rnc-lasso", "lasso")) %>%
  select(-s, -design, -iter) %>%
  mutate(method = replace(method, method == "funkl2", "funk l2")) %>%
  mutate(method = replace(method, method == "funkl1", "funk l1")) %>%
  ggplot(aes(x = rho, y = test_dev, color = method)) +
  geom_line() + 
  geom_point() +
  geom_errorbar(aes(ymin=test_dev-test_se, ymax=test_dev+test_se), width=0.01) +
  labs(y = "Root mean squared error", x = TeX(r"( $\rho$ )")) +
  scale_color_manual(values = cols) + 
  theme_bw()

p4 = res %>% 
  mutate_if(is.factor, as.character) %>%
  filter(design == "uninform") %>%
  filter(method %in% c("funkl1", "funkl2", "rnc-lasso", "lasso")) %>%
  select(-s, -design, -iter) %>%
  mutate(method = replace(method, method == "funkl2", "funk l2")) %>%
  mutate(method = replace(method, method == "funkl1", "funk l1")) %>%
  ggplot(aes(x = rho, y = cov_prob, color = method)) +
  geom_line() + 
  geom_point() +
  geom_errorbar(aes(ymin=cov_prob-cov_prob_se, ymax=cov_prob+cov_prob_se), width=0.01) +
  labs(y = "Coverage probability", x = TeX(r"( $\rho$ )")) +
  scale_color_manual(values = cols) +
  theme_bw()


ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

ggsave("uninform_s20_hd_r2r.pdf", device=cairo_pdf, height=7, width=7)

