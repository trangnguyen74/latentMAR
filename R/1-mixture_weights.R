

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


########################################
#### contradiction_plot

#' Contradiction plot
#'
#' Make the diagnostic plot called "contradiction plot", which is a histogram of assumption-implied response probabilities under control, for compliers under the SNR/rER assumption and for noncompliers under the SCR assumption.
#' @param miss.assumption Options are: "SNR" or "rER" (for stable noncomplier response / response exclusion restriction -- two names of the same assumption), or "SCR" (stable complier response)
#' @param pi.1,varpi.11,varpi.10,lambda.0 Vectors for nuisance functions. If the missingness assumption is SNR/rER, provide \code{varpi.10}. If the missingness assumption is SCR, provide \code{varpi.11}.
#' @export

contradiction_plot <- function(miss.assumption, pi.1, varpi.11, varpi.10, lambda.0,
                               bin.width, x.step) {

    if (miss.assumption %in% c("SNR", "rER")) {
        probs <- (lambda.0 - (1-pi.1)*varpi.10) / pi.1
        wts   <- pi.1
    } else if(miss.assumption=="nSCR") {
        probs <- (lambda.0 - pi.1*varpi.11) / (1-pi.1)
        wts   <- 1 - pi.1
    }

    in.range <- ifelse(probs<=1 & probs>=0, "in", "out")

    probs <- ifelse(probs==0, .000001,
                    ifelse(probs==1, 1-.00001,
                           probs))

    wts <- wts / sum(wts)

    pdat <- data.frame(probs = probs,
                       wts = wts,
                       in.range = in.range)


    p <- ggplot2::ggplot(data = pdat,
                         mapping = ggplot2::aes(x = probs,
                                                weight = wts))

    break.points <- seq(floor(min(probs)/bin.width)*bin.width,
                        ceiling(max(probs)/bin.width)*bin.width,
                        bin.width)
    # x.ticks <- seq(ceiling(min(probs)/.1)*.1,
    #                floor(max(probs)/.1)*.1,
    #                .1)

    if (all(in.range=="in")) {
        p <- p + ggplot2::geom_histogram(breaks = break.points, fill = "black")

    } else if (all(in.range=="out")) {
        p <- p + ggplot2::geom_histogram(breaks = break.points, fill = "red")

    } else {
        p <- p +
            ggplot2::geom_histogram(breaks = break.points,
                                    mapping = ggplot2::aes(fill = in.range)) +
            scale_fill_manual(values = c("black", "red"), guide = "none")
    }

    p <- p +
        # scale_x_continuous(breaks = x.ticks) +
        labs(x = "response probability",
             y = "bin mass") +
        theme_bw()



    left.mass <- sum((probs<0)*wts)
    right.mass <- sum((probs>1)*wts)

    if (left.mass + right.mass > 0) {
        mass.line <- "Out of bounds mass (marked in red):"

        if (right.mass==0)       { mass.line <- paste(mass.line, signif(left.mass, 3), "on the left")
        } else if (left.mass==0) { mass.line <- paste(mass.line, signif(right.mass, 3), "on the right")
        } else                   { mass.line <- paste(mass.line, signif(left.mass, 3), "on the left;", signif(right.mass, 3), "on the right")
        }

        p <- p + labs(subtitle = mass.line)
    }

    if (!missing(x.step)) {
        x.ticks <- seq(floor(min(probs)/x.step)*x.step,
                       ceiling(max(probs)/x.step)*x.step,
                       x.step)
        p <- p + scale_x_continuous(breaks = x.ticks)
    }

    p


}
