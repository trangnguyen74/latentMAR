

########################################
#### plugin_estimate

#' Compute the plug-in estimate of the CACE and NACE
#' @param pi.1,mu.11,mu.10,mu.01,mu.00 Vectors for estimated nuisance functions, computed at X values in the dataset
#' @export

plugin_estimate <- function(pi.1, mu.11, mu.10, mu.01, mu.00) {

    c(cace = weighted.mean(mu.11 - mu.01, nuis$pi.1),
      nace = weighted.mean(mu.10 - mu.00, 1-nuis$pi.1))
}


















