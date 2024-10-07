###############################################################################
## The original code is from https://gitlab.com/cbenard/vimp-causal-forests. ##
## A slight modification was made to track progress of the function.         ##
###############################################################################

### vimp function
#' Variable importance for causal forests
#'
#' @param c.forest A causal forest from the grf package.
#' @param variable.groups List of group of variables for group variable importance. Default is NULL, and each input variable forms a group.
#'
#' @return The list of importance values for all input variables or provided variable groups.
#' @export
#'
#' @examples
vimp_causal_forests <- function(c.forest, variable.groups = NULL){
  cat("Starting variable importance calculation...\n")
  
  # check input c.forest is a grf causal forest
  is.causal.forest <- all(class(c.forest) == c('causal_forest', 'grf'))
  if (!is.causal.forest){
    stop('c.forest must be a grf causal forest.')
  }
  
  # get initial data
  X <- c.forest$X.orig
  Y <- c.forest$Y.orig
  W <- c.forest$W.orig
  p <- ncol(X)
  n <- nrow(X)
  
  # get centered outcome and treatment assignment
  Y.centered <- Y - c.forest$Y.hat
  W.centered <- W - c.forest$W.hat
  
  # get oob predictions
  tau.hat <- c.forest$predictions
  
  # set variable groups
  if (is.null(variable.groups)){
    # when variable.groups is NULL, each input variable defines a group
    index.groups <- as.list(1:p)
  } else {
    # check provided variable groups are valid (non empty list, all variable names are contained in the data, and not all variables in one group)
    non.empty.list <- is.list(variable.groups) & length(variable.groups) > 0
    if (!non.empty.list){
      stop('variable.groups must be an non-empty list.')
    }
    names.valid <- all(unlist(variable.groups) %in% colnames(X))
    if (!names.valid){
      stop('Variable names provided in variable.groups not all found in input data.')
    }
    groups.valid <- max(sapply(variable.groups, function(group){length(unique(group))})) < p
    if (!groups.valid){
      stop('A group cannot contain all input variables.')
    }
    # get group of variable indices
    index.groups <- lapply(variable.groups, function(group){
      unique(sapply(group, function(j){which(colnames(X) == j)}))
    })
  }
  
  cat("Computing variable importance for each group...\n")
  total_groups <- length(index.groups)
  
  # compute vimp for all input variables using the settings of the initial causal forest
  In <- sapply(seq_along(index.groups), function(i) {
    index <- index.groups[[i]]
    cat(sprintf("Processing group %d of %d\n", i, total_groups))
    
    c.forest.drop.Xj <- causal_forest(X[, -index, drop=F], Y, W, Y.hat = c.forest$Y.hat, W.hat = c.forest$W.hat,
                                      num.trees = c.forest$`_num_trees`,
                                      sample.weights = c.forest$sample.weights,
                                      clusters = c.forest$clusters,
                                      equalize.cluster.weights = c.forest$equalize.cluster.weights,
                                      sample.fraction = c.forest$tunable.params$sample.fraction,
                                      mtry = min(c.forest$tunable.params$mtry, ncol(X) - length(index)),
                                      min.node.size = c.forest$tunable.params$min.node.size,
                                      honesty.fraction = c.forest$tunable.params$honesty.fraction,
                                      honesty.prune.leaves = c.forest$tunable.params$honesty.prune.leaves,
                                      alpha = c.forest$tunable.params$alpha,
                                      imbalance.penalty = c.forest$tunable.params$imbalance.penalty,
                                      ci.group.size = c.forest$ci.group.size)
    alpha <- get_forest_weights(c.forest.drop.Xj)
    vimp.Xj <- compute_vimp(alpha, Y.centered, W.centered, tau.hat)
    return(vimp.Xj)
  })
  
  cat("Variable importance computation complete.\n")
  cat("Computing retrain bias...\n")
  
  # compute retrain bias
  c.forest0 <- causal_forest(X, Y, W, Y.hat = c.forest$Y.hat, W.hat = c.forest$W.hat,
                             num.trees = c.forest$`_num_trees`,
                             sample.weights = c.forest$sample.weights,
                             clusters = c.forest$clusters,
                             equalize.cluster.weights = c.forest$equalize.cluster.weights,
                             sample.fraction = c.forest$tunable.params$sample.fraction,
                             mtry = c.forest$tunable.params$mtry,
                             min.node.size = c.forest$tunable.params$min.node.size,
                             honesty.fraction = c.forest$tunable.params$honesty.fraction,
                             honesty.prune.leaves = c.forest$tunable.params$honesty.prune.leaves,
                             alpha = c.forest$tunable.params$alpha,
                             imbalance.penalty = c.forest$tunable.params$imbalance.penalty,
                             ci.group.size = c.forest$ci.group.size)
  alpha <- get_forest_weights(c.forest0)
  In0 <- compute_vimp(alpha, Y.centered, W.centered, tau.hat)
  
  # compute debiased importance
  In <- In - In0
  
  cat("Variable importance calculation completed.\n")
  return(In)
}

# compute vimp from coefficients alpha of the retrained forest, centered outcome and treatment assignment, and original oob predictions
compute_vimp <- function(alpha, Y.centered, W.centered, tau.hat){
  W.bar <- alpha%*%W.centered
  W2.bar <- alpha%*%(W.centered^2)
  Y.bar <- alpha%*%Y.centered
  tau.bar <- alpha%*%tau.hat
  tau.hat.new <- (alpha%*%(W.centered*Y.centered - tau.hat*(W.centered^2)) - (W.bar*Y.bar - tau.bar*W2.bar))/(W2.bar - W.bar^2)
  vimp <- as.numeric(sum((tau.hat - tau.hat.new)^2)/(length(Y.centered)*var(tau.hat)))
  return(vimp)
}
