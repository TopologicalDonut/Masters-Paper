library(tidyverse)
library(grf)
library(lmtest)
library(sandwich)
library(ggplot2)
library(broom)
library(MASS)
library(modelsummary)

source("code/functions/main_functions.R")
source("code/functions/vimp_causal_forests.R")

# ---- Data Cleaning ----

data <- read.csv("data/raw/oreopoulos-resume-study-replication-data-file.csv")

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

# ---- Main Analysis ----

set.seed(1)

covariates <- c("ma", "female", "ba_quality", "exp_highquality", "language_skills",
                "extracurricular_skills", "same_exp")
treatment <- "english_dummy"
outcome <- "callback"

forest_results <- fit_causal_forest(
  filtered_data, 
  covariates, 
  treatment, 
  outcome
)

ate_result <- average_treatment_effect(forest_results$forest, target.sample = "all")
ate <- ate_result[1]
ate_se <- ate_result[2]
ci_lower <- ate - 1.96 * ate_se
ci_upper <- ate + 1.96 * ate_se

# Analyze rankings
ranking_results <- analyze_rankings(
  forest_results$aipw_scores,
  forest_results$ranking,
  num_rankings = 3
)

# Best Linear Projection
blp_results <- best_linear_projection(forest_results$forest, forest_results$X)
modelsummary(blp_results,
             stars = TRUE)

# ---- Subgroup Analysis ----

subgroups <- c("extracurricular_skills", "female", "ba_quality")

ate_subgroup_results <- lapply(subgroups, function(group) {
  result <- compute_subgroup_ate(group, forest_results$aipw_scores, filtered_data)
  result$Subgroup <- group
  result
})

subgroup_diff_results <- lapply(ate_subgroup_results, compute_subgroup_ATE_diff)

# Apply Romano-Wolf correction
num_boot <- 20000
t_stats <- sapply(subgroup_diff_results, function(x) x["t_stat"])
t_boot <- matrix(rnorm(length(t_stats) * num_boot), nrow = num_boot)
rw_pvalues <- romano_wolf_correction(t_stats, t_boot)

overview <- data.frame(
  Subgroup = subgroups,
  t_statistic = sapply(subgroup_diff_results, function(x) x["t_stat"]),
  Original_p = sapply(subgroup_diff_results, function(x) 2 * pt(abs(x["t_stat"]), df = nrow(filtered_data) - 2, lower.tail = FALSE)),
  RW_Adjusted_p = rw_pvalues
)

all_results <- do.call(rbind, ate_subgroup_results)

# Define and apply nice labels
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

# ---- Visualization ----

# Treatment effects distribution
ggplot(data.frame(tau_hat = forest_results$tau_hat), aes(x = tau_hat)) +
  geom_histogram(binwidth = (max(forest_results$tau_hat) - min(forest_results$tau_hat)) / 30,
                 fill = "lightblue", color = "white") +
  geom_vline(xintercept = ate, color = "red", linewidth = 1) +
  geom_vline(xintercept = c(ci_lower, ci_upper), color = "red",
             linetype = "dashed", linewidth = 0.5) +
  labs(x = "Treatment Effect",
       y = "Count",
       title = "Distribution of Treatment Effects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.1))

# Quantiles plot
ggplot(ranking_results, aes(x = ranking, y = estimate)) +
  geom_point(position = position_dodge(0.2)) +
  geom_errorbar(
    aes(ymin = estimate - 2*std_err, ymax = estimate + 2*std_err),
    width = 0.2,
    position = position_dodge(0.2)
  ) +
  labs(
    title = "Treatment Effects by Ranking",
    y = "ATE estimate",
    x = "Ranking"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())

# Subgroup ATE plot
ggplot(all_results, aes(x = Group, y = ATE, color = Subgroup)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~ Subgroup, scales = "free_x") +
  geom_text(
    data = overview,
    aes(x = 1.5, y = -Inf, 
        label = sprintf("p-value = %.3f", Original_p)),
    hjust = 0.5, vjust = -0.5, size = 3, color = "black"
  ) +
  theme_minimal() +
  labs(
    title = "ATEs and 95% Confidence Intervals by Subgroup",
    y = "Average Treatment Effect (ATE)",
    x = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# ---- Variable Importance ----

# The output is a vector of numbers, they correspond to the covariates in order.
vimp <- vimp_causal_forests(forest_results$forest)
vimp

