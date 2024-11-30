#' Computes adjusted p-values following the Romano-Wolf method.
#' For reference, see http://ftp.iza.org/dp12845.pdf page 8
#' 
#' @param t_orig  Vector of t-statistics from original model 
#' @param t_boot  Matrix of t-statistics from bootstrapped models
#' 
#' @return        Vector of adjusted p-values
romano_wolf_correction <- function(t_orig, t_boot) {
  abs_t_orig <- abs(t_orig)
  abs_t_boot <- abs(t_boot)
  abs_t_sorted <- sort(abs_t_orig, decreasing = TRUE)
  
  max_order <- order(abs_t_orig, decreasing = TRUE)
  rev_order <- order(max_order)
  
  M <- nrow(t_boot)
  S <- ncol(t_boot)
  
  p_adj <- rep(0, S)
  p_adj[1] <- mean(apply(abs_t_boot, 1, max) > abs_t_sorted[1])
  
  for (s in seq(2, S)) {
    cur_index <- max_order[s:S]
    p_init <- mean(apply(abs_t_boot[, cur_index, drop=FALSE], 1, max) > abs_t_sorted[s])
    p_adj[s] <- max(p_init, p_adj[s-1])
  }
  
  p_adj[rev_order]
}

#' Computes adjusted p-values for linear regression (lm) models in particular.
#' 
#' @param model     The fitted lm model
#' @param indices   Indices of coefficients to adjust (default: all)
#' @param cov_type  Type of robust covariance estimation
#' @param num_boot  Number of bootstrap iterations
#' 
#' @return         Data frame with original and adjusted statistics
summary_rw_lm <- function(model, indices=NULL, cov_type="HC2", num_boot=10000) {
  if (is.null(indices)) {
    indices <- 1:nrow(coef(summary(model)))
  }
  
  # Get original t values
  summary <- coef(summary(model))[indices,,drop=FALSE]
  t_orig <- summary[, "t value"]
  
  # Null resampling
  # This is a bit of trick to speed up bootstrapping linear models.
  # Here, we don't really need to re-fit linear regressions, which would be a bit slow.
  # We know that betahat ~ N(beta, Sigma), and we have an estimate Sigmahat.
  # So we can approximate "null t-values" by
  #  - Draw beta.boot ~ N(0, Sigma-hat)     note the 0 here, this is what makes it a *null* t-value.
  #  - Compute t.boot = beta.boot / sqrt(diag(Sigma.hat))
  sigma_hat <- sandwich::vcovHC(model, type=cov_type)[indices, indices]
  se_orig <- sqrt(diag(sigma_hat))
  num_coef <- length(se_orig)
  
  beta_boot <- MASS::mvrnorm(n=num_boot, mu=rep(0, num_coef), Sigma=sigma_hat)
  t_boot <- sweep(beta_boot, 2, se_orig, "/")
  
  p_adj <- romano_wolf_correction(t_orig, t_boot)
  
  result <- cbind(summary[,c(1,2,4),drop=F], p_adj)
  colnames(result) <- c('Estimate', 'Std. Error', 'Orig. p-value', 'Adj. p-value')
  result
}

#' Fits a causal forest on a given set of data and specified covariates, treatment,
#' and outcome. Also ranks CATEs into quantiles and computes AIPW scores.
#'
#' @param data          Clean dataset for analysis.
#' @param covariates    Vector of covariate names.
#' @param treatment     Name of treatment variable.
#' @param outcome       Name of outcome variable.
#' @param num_rankings  Number of quantiles for ranking CATEs.
#' @param num_folds     Number of cross-validation folds.
#' @param num_trees     Number of trees in the forest.
#' @param w_hat         Optional propensity scores.
#' 
#' @return              List containing trained forest and other required outputs.
fit_causal_forest <- function(data, covariates, treatment, outcome, num_rankings = 3,
                              num_folds = 5, num_trees = 4000, W_hat = NULL) {
  
  # Prepare model matrices
  fmla <- formula(paste0("~ 0 + ", paste0(covariates, collapse="+")))
  X <- model.matrix(fmla, data)
  W <- data[,treatment]
  Y <- data[,outcome]
  n <- nrow(data)
  
  # Assign folds for cross-fitting
  folds <- sample(rep(1:num_folds, length.out = n))
  
  # Fit forest
  forest <- causal_forest(X, Y, W, 
                          W.hat = W_hat,
                          clusters = folds, 
                          num.trees = num_trees,
                          tune.parameters = "all")
  
  # Get predictions and compute AIPW scores
  tau_hat <- predict(forest)$predictions
  e_hat <- forest$W.hat
  m_hat <- forest$Y.hat
  
  mu_hat_0 <- m_hat - e_hat * tau_hat
  mu_hat_1 <- m_hat + (1 - e_hat) * tau_hat
  
  aipw_scores <- tau_hat + 
    W / e_hat * (Y - mu_hat_1) - 
    (1 - W) / (1 - e_hat) * (Y - mu_hat_0)
  
  # Rank observations within folds into quantiles according to their CATE predictions.
  ranking <- rep(NA, n)
  for (fold in seq(num_folds)) {
    tau_hat_quantiles <- quantile(tau_hat[folds == fold], 
                                  probs = seq(0, 1, by=1/num_rankings))
    ranking[folds == fold] <- cut(tau_hat[folds == fold], 
                                  tau_hat_quantiles, 
                                  include.lowest=TRUE,
                                  labels=seq(num_rankings))
  }
  
  list(
    forest = forest,
    tau_hat = tau_hat,
    ranking = ranking,
    X = X,
    W = W,
    Y = Y,
    aipw_scores = aipw_scores
  )
}

#' Analyze heterogeneous treatment effects by ranking
#'
#' @param aipw_scores   AIPW scores from forest analysis.
#' @param ranking       Ranking assignments from forest analysis.
#' @param num_rankings  Number of rankings used (same as input to `fit_causal_forest`).
#' 
#' @return              Data frame with results by ranking.
analyze_rankings <- function(aipw_scores, ranking, num_rankings) {
  ols <- lm(aipw_scores ~ 0 + factor(ranking))
  forest_ate <- data.frame(
    "method" = "aipw",
    "ranking" = paste0("Q", seq(num_rankings)), 
    coeftest(ols, vcov=vcovHC(ols, "HC2"))[,1:2]
  )
  colnames(forest_ate)[3:4] <- c("estimate", "std_err")
  forest_ate
}

#' Computes ATE for different subgroups.
#'
#' @param group        Subgroup variable name
#' @param aipw_scores  AIPW scores from forest analysis.
#' @param data         Dataset containing subgroup variables (usually same as forest dataset).
#' 
#' @return             Data frame with ATE results by subgroup.
compute_subgroup_ate <- function(group, aipw_scores, data) {
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

#' Computes difference in ATEs between groups
#' 
#' @param ate_subgroup_result  Data frame with ATE results for a subgroup
#' 
#' @return                     Named vector with difference, standard error, and t-statistic
compute_subgroup_ATE_diff <- function(ate_subgroup_result) {
  diff <- ate_subgroup_result$ATE[2] - ate_subgroup_result$ATE[1]
  se_diff <- sqrt(sum(ate_subgroup_result$SE^2))
  t_stat <- diff / se_diff
  c(diff = diff, se = se_diff, t_stat = t_stat)
}
