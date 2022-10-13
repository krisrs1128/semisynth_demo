
#' @param effect_estimates Length K vector of species effect estimates
#' @param true_indicators Length K binary vector of whether a species has a true
#'   effect
metrics <- function(effect_estimates, true_indicators) {
  K <- length(effect_estimates)
  precision <- vector(length = K)
  recall <- vector(length = K)
  ix <- order(effect_estimates)
  
  for (k in seq_len(K)) {
    precision[k] <- mean(true_indicators[1:ix[k]])
    recall[k] <- sum(true_indicators[1:ix[k]]) / sum(true_indicators)
  }
  
  list(precision = precision, recall = recall)
}