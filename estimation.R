
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

fsr_path <- function(intervals, truth, grid_ix = NULL) {
  f <- function(candidate_intervals, truth, ix) {
    n_false_sign <- sum(sign(candidate_intervals[, 1]) != sign(truth[ix]))
    n_false_sign / max(nrow(candidate_intervals), 1)
  }
  eval_path(intervals, truth, f, grid_ix)
}

power_path <- function(intervals, truth, grid_ix = NULL) {
  f <- function(candidate_intervals, truth, ix) {
    n_true_sign <- sum(sign(candidate_intervals[, 1]) == sign(truth[ix]))
    n_true_sign / max(sum(truth != 0), 1)
  }
  eval_path(intervals, truth, f, grid_ix)
}

eval_path <- function(intervals, truth, metric_fun, grid_ix = NULL) {
  if (is.null(grid_ix)) {
    grid_ix <- seq(0, max(abs(intervals)), length.out = n_grid)
  }
  
  result <- vector(length = length(grid_ix))
  for(i in seq_along(result)) {
    ix <- apply(abs(intervals), 1, min) > grid_ix[i] & as.logical(sign(intervals[, 1]) == sign(intervals[, 2]))
    result[i] <- metric_fun(intervals[ix, ], truth, ix)
  }
  result
}

postprocess_posterior <- function(posterior) {
  intervals <- t(apply(posterior, 2, function(x) quantile(x, c(0.05, 0.95)))) %>%
    as.data.frame() %>%
    add_rownames("parameter") %>%
    as_tibble() %>%
    mutate(
      feature = parameter %>%
        str_extract("[0-9]+,") %>%
        str_remove(",") %>%
        as.integer(),
      taxon = parameter %>%
        str_extract(",[0-9]+") %>%
        str_remove(",") %>%
        as.integer(),
    )
}