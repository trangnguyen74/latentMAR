

########################################
#### mixture_weights

#' Compute mixture weights
#'
#' Compute the mixture weights to be used in the outcome mixture equation based on the specific missingness assumption and estimated nuisance functions.
#' @param miss.assumption A character string indicating the specific missingness assumption to be used. Options are: "nSNR" (near stable noncomplier response), "nSCR" (near stable complier response), "rPI" (principal ignorability for response) and "rPO" (proportional response odds)
#' @param pi.1,varpi.11,varpi.10,lambda.0 Vectors for nuisance functions. Use the one(s) needed for the specific missingness assumption.
#' @export

mixture_weights <- function(miss.assumption, pi.1, varpi.11, varpi.10, lambda.0) {

    if(miss.assumption=="rPI")
        return(list(pi.01R = pi.1,
                    pi.00R = 1 - pi.1))


    epsilon <- 0.03
    trim_prob <- function(p, l, h) {
        (p<=l)*l + (p>=h)*h + (p>l & p<h)*p
    }

    lambda.0 <- trim_prob(lambda.0, epsilon, 1)


    if (miss.assumption=="nSNR") {
        varpi.10 <- trim_prob(varpi.10, epsilon, 1)

        varpi.01 <- (lambda.0 - (1-pi.1)*varpi.10) / pi.1
        varpi.01 <- trim_prob(varpi.01, epsilon, 1)

        pi.01R <- pi.1 * varpi.01 / lambda.0
        pi.00R <- 1 - pi.01R

    } else if(miss.assumption=="nSCR") {
        varpi.11 <- trim_prob(varpi.11, epsilon, 1)

        varpi.00 <- (lambda.0 - pi.1*varpi.11) / (1-pi.1)
        varpi.00 <- trim_prob(varpi.00, epsilon, 1)

        pi.00R <- (1-pi.1) * varpi.00 / lambda.0
        pi.01R <- 1 - pi.00R

    } else if(miss.assumption=="rPO") {
        varpi.11 <- trim_prob(varpi.11, epsilon, 1)
        varpi.10 <- trim_prob(varpi.10, epsilon, 1)

        varrho <- (varpi.11/(1-varpi.11))/(varpi.10/(1-varpi.10))

        gamma <- (pi.1+lambda.0)*(varrho-1)+1

        delta <- sqrt(gamma^2-4*pi.1*lambda.0*varrho*(varrho-1))

        varpi.01 <- (gamma-delta)/(2*(varrho-1)*pi.1)
        varpi.01 <- trim_prob(varpi.01, epsilon, 1)

        varpi.00 <- (lambda.0 - pi.1*varpi.01) / (1-pi.1)
        varpi.00 <- trim_prob(varpi.00, epsilon, 1)

        varpi.01 <- (lambda.0 - (1-pi.1)*varpi.00) / pi.1

        pi.01R <- pi.1 * varpi.01 / lambda.0
        pi.00R <- (1-pi.1) * varpi.00 / lambda.0

        pi.01R <- ifelse(varrho==1, pi.1,   pi.01R)
        pi.00R <- ifelse(varrho==1, 1-pi.1, pi.00R)
    }


    list(pi.01R = pi.01R,
         pi.00R = pi.00R)
}
