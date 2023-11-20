

########################################
#### mus_under_control

#' @param principal.assumption A character string indicating the principal identification assumption to be used. Options are: "ER", "PI", "PIsens-GOR", "PIsens-MR", "PIsens-SMDe"
#' @param mix.wts A list of two vectors named pi.01R and pi.00R holding the mixture weights to be used in the outcome mixture equation. Output of function \code{mixture_weights()}.
#' @param mu.10,kappa.0R,varsigma.0R Vectors for nuisance functions. Use the one(s) needed for the principal identification assumption.
#' @param psi,rho,eta Sensitivity parameters required for the PIsens assumptions. Use psi for PIsens-GOR, rho for PIsens-MR and eta for PIsens-SMDe.
#' @param y.bounds A 2-length numeric vector holding the lower and upper bounds of the outcome. Needed if using PIsens-GOR and the outcome range is not [0,1].
#' @export


mus_under_control <- function(principal.assumption,
                              mix.wts, mu.10, kappa.0R, varsigma.0R,
                              psi, rho, eta,
                              y.bounds = NULL,
                              x.dat = NULL,
                              x.vars = NULL,
                              pi.1 = NULL) {


    if (principal.assumption=="PI")
        return(list(mu.01 = kappa.0R,
                    mu.00 = kappa.0R))


    if (principal.assumption=="ER") {
        mu.01 <- mu.10 + (kappa.0R - mu.10)/mix.wts$pi.01R

        # If want to deal with the issue that the ER assumption may estimate mu.01(X) that is out of range,
        # can use a bounded projection similar to that in Wang & Tchetgen Tchetgen.
        # (W & T bounded the effect in the [-1,1] range; this projection her bounds mu.01(X) within its natural bounds.)
        # (We don't do this in the paper -- keeping the focus on the missingness problem, not on primary PCE identification.)

        if (!is.null(y.bounds) & sum((mu.01>max(y.bounds)) + (mu.01<min(y.bounds))) > 0) {

            if (!is.null(x.dat)) {

                x.formula <- as.formula(paste("~", paste(x.vars, collapse = "+")))
                x.mat     <- model.matrix(x.formula, data = x.dat)

                project.obj <- .wang_bounded_projection(y.vec    = mu.01,
                                                        X.mat    = x.mat,
                                                        w.vec    = pi.1,
                                                        y.bounds = y.bounds)

                if (project.obj$convergence) {
                    cat("\n\n\nProjection model converged.\n\n\n")
                } else {
                    cat("\n\nProjection model did not converge.\n\n\n")
                }

                mu.01 <- project.obj$projection
            }
        }
        return(list(mu.00 = mu.10,
                    mu.01 = mu.01))
    }


    # This nearER is included here just to allow us to explore symmetric treatment of principal ID assumptions and missingness assumptions, out of curiosity.
    # We don't consider this as an assumption in the paper.
    if (principal.assumption=="nearER") {

        l <- min(y.bounds)
        h <- max(y.bounds)

        mu.01 <- mu.10 + (kappa.0R - mu.10)/mix.wts$pi.01R
        mu.01 <- (mu.01<=l)*l + (mu.01>=h)*h + (mu.01>l & mu.01<h)*mu.01
        mu.00 <- mu.01 + (kappa.0R - mu.01) / mix.wts$pi.00R
        return(list(mu.00 = mu.00,
                    mu.01 = mu.01))
    }


    if (principal.assumption=="PIsens-GOR") {
        psi.1 <- psi
        psi.0 <- 1/psi

        if (is.null(y.bounds)) {

            if (max(kappa.0R>1) | min(kappa.0R<0)) {
                stop(paste0(c("It looks like the outcome is not binary.",
                              "For PIsens-GOR, argument y.bounds needs to be specified")))
            } else {
                y.bounds <- c(0,1)
            }
        }

        kappa.diamond <- (kappa.0R-min(y.bounds)) / abs(diff(y.bounds))

        alpha.1 <- (mix.wts$pi.01R + kappa.diamond)*(psi.1-1) + 1
        alpha.0 <- (mix.wts$pi.00R + kappa.diamond)*(psi.0-1) + 1

        beta.1 <- sqrt(alpha.1^2 - 4*mix.wts$pi.01R*kappa.diamond*psi.1*(psi.1-1))
        beta.0 <- sqrt(alpha.0^2 - 4*mix.wts$pi.00R*kappa.diamond*psi.0*(psi.0-1))

        mu.01.diamond <- (alpha.1-beta.1) / (2*(psi.1-1)*mix.wts$pi.01R)
        mu.00.diamond <- (alpha.0-beta.0) / (2*(psi.0-1)*mix.wts$pi.00R)

        mu.01.diamond <- ifelse(mix.wts$pi.01R==0, 0, mu.01.diamond)
        mu.00.diamond <- ifelse(mix.wts$pi.00R==0, 0, mu.00.diamond)

        return(list(mu.01 = mu.01.diamond * abs(diff(y.bounds)) + min(y.bounds),
                    mu.00 = mu.00.diamond * abs(diff(y.bounds)) + min(y.bounds)))
    }


    if (principal.assumption=="PIsens-MR") {
        scale.0 <- 1 / ((rho-1)*mix.wts$pi.01R + 1)
        scale.1 <- scale.0*rho

        return(list(mu.01 = scale.1 * kappa.0R,
                    mu.00 = scale.0 * kappa.0R))
    }


    if (principal.assumption=="PIsens-SMDe") {
        denom <- sqrt(1 + eta^2 * mix.wts$pi.01R * mix.wts$pi.00R)

        return(list(mu.01 = kappa.0R + eta*(1-mix.wts$pi.01R)*varsigma.0R/denom,
                    mu.00 = kappa.0R - eta*(1-mix.wts$pi.00R)*varsigma.0R/denom))
    }
}
