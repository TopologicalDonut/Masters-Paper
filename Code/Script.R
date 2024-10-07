# Romano-Wolf Correction Code ---------------------------------------------

# Computes adjusted p-values following the Romano-Wolf method.
# For a reference, see http://ftp.iza.org/dp12845.pdf page 8
#  t.orig: vector of t-statistics from original model 
#  t.boot: matrix of t-statistics from bootstrapped models
romano_wolf_correction <- function(t.orig, t.boot) {
  abs.t.orig <- abs(t.orig)
  abs.t.boot <- abs(t.boot)
  abs.t.sorted <- sort(abs.t.orig, decreasing = TRUE)
  
  max.order <- order(abs.t.orig, decreasing = TRUE)
  rev.order <- order(max.order)
  
  M <- nrow(t.boot)
  S <- ncol(t.boot)
  
  p.adj <- rep(0, S)
  p.adj[1] <- mean(apply(abs.t.boot, 1, max) > abs.t.sorted[1])
  for (s in seq(2, S)) {
    cur.index <- max.order[s:S]
    p.init <- mean(apply(abs.t.boot[, cur.index, drop=FALSE], 1, max) > abs.t.sorted[s])
    p.adj[s] <- max(p.init, p.adj[s-1])
  }
  p.adj[rev.order]
}

# Computes adjusted p-values for linear regression (lm) models.
summary_rw_lm <- function(model, indices=NULL, cov.type="HC2", num.boot=10000, seed=2020) {
  
  if (is.null(indices)) {
    indices <- 1:nrow(coef(summary(model)))
  }
  # Grab the original t values.
  summary <- coef(summary(model))[indices,,drop=FALSE]
  t.orig <- summary[, "t value"]
  
  # Null resampling.
  # This is a bit of trick to speed up bootstrapping linear models.
  # Here, we don't really need to re-fit linear regressions, which would be a bit slow.
  # We know that betahat ~ N(beta, Sigma), and we have an estimate Sigmahat.
  # So we can approximate "null t-values" by
  #  - Draw beta.boot ~ N(0, Sigma-hat)     note the 0 here, this is what makes it a *null* t-value.
  #  - Compute t.boot = beta.boot / sqrt(diag(Sigma.hat))
  Sigma.hat <- sandwich::vcovHC(model, type=cov.type)[indices, indices]
  se.orig <- sqrt(diag(Sigma.hat))
  num.coef <- length(se.orig)
  beta.boot <- MASS::mvrnorm(n=num.boot, mu=rep(0, num.coef), Sigma=Sigma.hat)
  t.boot <- sweep(beta.boot, 2, se.orig, "/")
  p.adj <- romano_wolf_correction(t.orig, t.boot)
  
  result <- cbind(summary[,c(1,2,4),drop=F], p.adj)
  colnames(result) <- c('Estimate', 'Std. Error', 'Orig. p-value', 'Adj. p-value')
  result
}

# Data Cleaning -----------------------------------------------------------

library(tidyverse)

data <- read.csv("data/oreopoulos-resume-study-replication-data-file.csv")

data_clean <- data %>%
  mutate(english_dummy = case_when(name_ethnicity == "British" ~ 1, name_ethnicity == "Canada" ~ 1, TRUE ~ 0),
         same_exp = ifelse(is.na(same_exp), 0, same_exp),
         reference = ifelse(is.na(reference), 0, reference),
         accreditation = ifelse(is.na(accreditation), 0, accreditation),
         legal = ifelse(is.na(legal), 0, legal),
         extracurricular_skills = ifelse(is.na(extracurricular_skills), 0, extracurricular_skills),
         type2 = case_when(type == 2 ~ 1, TRUE ~ 0),
         type3 = case_when(type == 3 ~ 1, TRUE ~ 0),
         type4 = case_when(type == 4 ~ 1, TRUE ~ 0)
  )

filtered_data <- data_clean %>%
  filter(type %in% c(0,1) & name_ethnicity != "Chn-Cdn")

# Causal Forest -----------------------------------------------------------

library(grf)
library(rpart)
library(lmtest)
library(sandwich)
library(ggplot2)

set.seed(1)

# accredition, reference, and legal need to be excluded from covariates since 
# none of the natives have these.
covariates <- c("ma", "female", "ba_quality", "exp_highquality", "language_skills",
                "extracurricular_skills", "same_exp")
treatment <- "english_dummy"
outcome <- "callback"

n <- nrow(filtered_data)

fmla <- formula(paste0("~ 0 + ", paste0(covariates, collapse="+")))
X <- model.matrix(fmla, filtered_data)
W <- filtered_data[,treatment]
Y <- filtered_data[,outcome]

# Number of rankings that the predictions will be ranking on
# (e.g., 2 for above/below median estimated CATE, 5 for estimated CATE quintiles, etc.)
num.rankings <- 3 

# Assign a fold number to each observation.x
# The argument 'clusters' in the causal forest will mimic K-fold cross-fitting.
num.folds <- 5
folds <- sample(rep(1:num.folds, length.out = n))

forest <- causal_forest(X, Y, W, W.hat = 0.5, clusters = folds, num.trees = 4000, tune.parameters = "all")

tau.hat <- predict(forest)$predictions

# Rank observations *within each fold* into quantiles according to their CATE predictions.
ranking <- rep(NA, n)
for (fold in seq(num.folds)) {
  tau.hat.quantiles <- quantile(tau.hat[folds == fold], probs = seq(0, 1, by=1/num.rankings))
  ranking[folds == fold] <- cut(tau.hat[folds == fold], tau.hat.quantiles, include.lowest=TRUE,labels=seq(num.rankings))
}

# Basic Predictions -------------------------------------------------------

ate_result <- average_treatment_effect(forest, target.sample = "all")
ate <- ate_result[1]
ate_se <- ate_result[2]
ci_lower <- ate - 1.96 * ate_se
ci_upper <- ate + 1.96 * ate_se

ggplot(data.frame(tau_hat = tau.hat), aes(x = tau_hat)) +
  geom_histogram(binwidth = (max(tau.hat) - min(tau.hat)) / 30, 
                 fill = "lightblue", color = "white") +
  geom_vline(xintercept = ate, color = "red", linewidth = 1) +
  geom_vline(xintercept = c(ci_lower, ci_upper), color = "red", 
             linetype = "dashed", linewidth = 0.5) +
  labs(x = "Treatment Effect",
       y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.1))

summary(tau.hat)
variable_importance(forest)

# Best Linear Projection and table ----------------------------------------

library(broom)
library(kableExtra)

blp_results <- best_linear_projection(forest, X)
blp_table <- tidy(blp_results)
blp_table

# Quantiles and Heatmap ---------------------------------------------------

e.hat <- forest$W.hat # P[W=1|X]
m.hat <- forest$Y.hat # E[Y|X]

# Estimating mu.hat(X, 1) and mu.hat(X, 0) for obs in held-out sample
# Note: to understand this, read equations 6-8 in this vignette:
# https://grf-labs.github.io/grf/articles/muhats.html
mu.hat.0 <- m.hat - e.hat * tau.hat        # E[Y|X,W=0] = E[Y|X] - e(X)*tau(X)
mu.hat.1 <- m.hat + (1 - e.hat) * tau.hat  # E[Y|X,W=1] = E[Y|X] + (1 - e(X))*tau(X)

# AIPW scores
aipw.scores <- tau.hat + W / e.hat * (Y -  mu.hat.1) - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0)
ols <- lm(aipw.scores ~ 0 + factor(ranking))
forest.ate <- data.frame("aipw", paste0("Q", seq(num.rankings)), coeftest(ols, vcov=vcovHC(ols, "HC2"))[,1:2])
colnames(forest.ate) <- c("method", "ranking", "estimate", "std.err")
rownames(forest.ate) <- NULL # just for display
forest.ate

# Plotting the point estimate of average treatment effect 
# and 95% confidence intervals around it.
ggplot(forest.ate) +
  aes(x = ranking, y = estimate) + 
  geom_point(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=estimate-2*std.err, ymax=estimate+2*std.err), width=.2, position=position_dodge(0.2)) +
  ylab("ATE estimate") + xlab("Ranking") +
  theme_minimal() +
  theme(legend.position="bottom", legend.title = element_blank())

# Using AIPW scores computed above
ols <- lm(aipw.scores ~ 1 + factor(ranking))
res <- summary_rw_lm(ols, indices=2:num.rankings)
rownames(res) <- paste("Rank", 2:num.rankings, "- Rank 1") # just for display
res

df <- mapply(function(covariate) {
  # Looping over covariate names
  # Compute average covariate value per ranking (with correct standard errors)
  fmla <- formula(paste0(covariate, "~ 0 + ranking"))
  ols <- lm(fmla, data=transform(filtered_data, ranking=factor(ranking)))
  ols.res <- coeftest(ols, vcov=vcovHC(ols, "HC2"))
  
  avg <- ols.res[,1]
  stderr <- ols.res[,2]
  
  data.frame(covariate, avg, stderr, ranking=paste0("Q", seq(num.rankings)),
             scaling=pnorm((avg - mean(avg))/sd(avg)),
             # We will order based on how much variation is 'explain' by the averages
             # relative to the total variation of the covariate in the data
             variation=sd(avg) / sd(filtered_data[,covariate]),
             labels=paste0(signif(avg, 3), "\n", "(", signif(stderr, 3), ")"))
}, covariates, SIMPLIFY = FALSE)
df <- do.call(rbind, df)

df$covariate <- factor(df$covariate, levels = unique(df$covariate[order(df$variation, decreasing = TRUE)]))

ggplot(df) +
  aes(ranking, covariate) +
  geom_tile(aes(fill = scaling)) + 
  geom_text(aes(label = labels)) +
  scale_fill_gradient(low = "#E1BE6A", high = "#40B0A6") +
  scale_y_discrete(limits = rev(levels(df$covariate)),
                   labels = c("female" = "Female",
                              "ba_quality" = "Top 200 University",
                              "extracurricular_skills" = "Extracurricular Skills Listed",
                              "ma" = "Canadian Master's Degree",
                              "exp_highquality" = "High Quality Experience",
                              "same_exp" = "Multinational Experience",
                              "language_skills" = "Fluent in Multiple Languages")) +
  theme_minimal() + 
  ylab("") + xlab("CATE estimate ranking") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 11))

# ATE Analysis ------------------------------------------------------------

subgroups <- c("extracurricular_skills", "female", "ba_quality")

compute_ate <- function(group, aipw_scores, data) {
  group_data <- split(aipw_scores, data[[group]])
  means <- sapply(group_data, mean)
  
  compute_se <- function(x) {
    n <- length(x)
    sd(x)/sqrt(n)
  }
  
  ses <- sapply(group_data, compute_se)
  
  data.frame(
    Group = names(means),
    ATE = means,
    SE = ses,
    CI_lower = means - 1.96 * ses,
    CI_upper = means + 1.96 * ses
  )
}

ate_results <- lapply(subgroups, function(group) {
  result <- compute_ate(group, aipw.scores, filtered_data)
  result$Subgroup <- group
  result
})

compute_diff <- function(ate_result) {
  diff <- ate_result$ATE[2] - ate_result$ATE[1]
  se_diff <- sqrt(sum(ate_result$SE^2))
  t_stat <- diff / se_diff
  c(diff = diff, se = se_diff, t_stat = t_stat)
}

diff_results <- lapply(ate_results, compute_diff)

t_stats <- sapply(diff_results, function(x) x["t_stat"])

# Apply Romano-Wolf correction
num_boot <- 20000
t_boot <- matrix(rnorm(length(t_stats) * num_boot), nrow = num_boot)
rw_pvalues <- romano_wolf_correction(t_stats, t_boot)

overview <- data.frame(
  Subgroup = subgroups,
  t_statistic = sapply(diff_results, function(x) x["t_stat"]),
  Original_p = sapply(diff_results, function(x) 2 * pt(abs(x["t_stat"]), df = nrow(filtered_data) - 2, lower.tail = FALSE)),
  RW_Adjusted_p = rw_pvalues
)

all_results <- do.call(rbind, ate_results)

subgroup_names <- c(
  "extracurricular_skills" = "Extracurricular Skills Listed",
  "female" = "Female",
  "ba_quality" = "Top 200 University"
)

all_results$Subgroup <- factor(all_results$Subgroup, 
                               levels = names(subgroup_names), 
                               labels = subgroup_names)

overview$Subgroup <- factor(overview$Subgroup, 
                            levels = names(subgroup_names), 
                            labels = subgroup_names)

p <- ggplot(all_results, aes(x = Group, y = ATE, color = Subgroup)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, 
                position = position_dodge(width = 0.5)) +
  facet_wrap(~ Subgroup, scales = "free_x") +
  geom_text(data = overview, 
            aes(x = 1.5, y = -Inf, 
                label = sprintf("p-value = %.3f", Original_p)),
            hjust = 0.5, vjust = -0.5, size = 3, color = "black") +
  theme_minimal() +
  labs(title = "ATEs and 95% Confidence Intervals by Subgroup",
       y = "Average Treatment Effect (ATE)",
       x = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

# Variable Importance -----------------------------------------------------
source("vimp_causal_forests.R")

vimp <- vimp_causal_forests(forest)
vimp

