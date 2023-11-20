


########################################
#### effects_all_assumptions

#' Compute effects under all choices of assumptions
#' This function computes effects under all choices of assumptions as in Tables 4 and 5 of the paper. The main use of this function is in bootstrapping.
#' @param nuis The list of estimated nuisance functions. Like output of \code{j2_nuisance_estimation()} or \code{ec_nuisance_estimation()}.
#' @param sens.type Type of sensitivity analysis (for PI violation) to be included.
#' @param sens.params Numeric vector including maximum and minimum values of the sensitivity parameter.
#' @param y.bounds Numeric vector holding the lower and upper bounds of the outcome. Needed only if \code{sens.type = "GOR"} and the outcome is not binary.
#' @export

effects_all_assumptions <- function(nuis,
                                    sens.type,
                                    sens.params,
                                    y.bounds = NULL,
                                    weights) {

    effects <- matrix(NA, ncol = 3, nrow = 13)
    colnames(effects) <- c("cace", "nace", "ate")
    rownames(effects) <- c("ER.nSNR", "ER.nSCR", "ER.rPI", "ER.rPO",
                           "PI",
                           "PIsens.nSNR.lo", "PIsens.nSNR.hi",
                           "PIsens.nSCR.lo", "PIsens.nSCR.hi",
                           "PIsens.rPI.lo", "PIsens.rPI.hi",
                           "PIsens.rPO.lo", "PIsens.rPO.hi")

    mix.wts.nSNR <- mixture_weights(miss.assumption = "nSNR",
                                    pi.1 = nuis$pi.1,
                                    varpi.10 = nuis$varpi.10,
                                    lambda.0 = nuis$lambda.0)

    mix.wts.nSCR <- mixture_weights(miss.assumption = "nSCR",
                                    pi.1 = nuis$pi.1,
                                    varpi.11 = nuis$varpi.11,
                                    lambda.0 = nuis$lambda.0)

    mix.wts.rPI <- mixture_weights(miss.assumption = "rPI",
                                   pi.1 = nuis$pi.1)

    mix.wts.rPO <- mixture_weights(miss.assumption = "rPO",
                                   pi.1 = nuis$pi.1,
                                   varpi.11 = nuis$varpi.11,
                                   varpi.10 = nuis$varpi.10,
                                   lambda.0 = nuis$lambda.0)

    plugin <- function(nuis, mus.0, wt) {
        c(cace = weighted.mean(nuis$mu.11 - mus.0$mu.01, wt * nuis$pi.1),
          nace = weighted.mean(nuis$mu.10 - mus.0$mu.00, wt * (1-nuis$pi.1)),
          ate  = weighted.mean(nuis$pi.1*(nuis$mu.11-mus.0$mu.01) + (1-nuis$pi.1)*(nuis$mu.10-mus.0$mu.00), wt))
    }


    effects["ER.nSNR",] <- plugin(nuis = nuis,
                                  mus.0 = mus_under_control(principal.assumption = "ER",
                                                            mix.wts = mix.wts.nSNR,
                                                            mu.10 = nuis$mu.10,
                                                            kappa.0R = nuis$kappa.0R),
                                  wt = weights)

    effects["ER.nSCR",] <- plugin(nuis = nuis,
                                  mus.0 = mus_under_control(principal.assumption = "ER",
                                                            mix.wts = mix.wts.nSCR,
                                                            mu.10 = nuis$mu.10,
                                                            kappa.0R = nuis$kappa.0R),
                                  wt = weights)

    effects["ER.rPI",] <- plugin(nuis = nuis,
                                 mus.0 = mus_under_control(principal.assumption = "ER",
                                                           mix.wts = mix.wts.rPI,
                                                           mu.10 = nuis$mu.10,
                                                           kappa.0R = nuis$kappa.0R),
                                 wt = weights)

    effects["ER.rPO",] <- plugin(nuis = nuis,
                                 mus.0 = mus_under_control(principal.assumption = "ER",
                                                           mix.wts = mix.wts.rPO,
                                                           mu.10 = nuis$mu.10,
                                                           kappa.0R = nuis$kappa.0R),
                                 wt = weights)

    # effects["nearER.nSNR",] <- plugin(nuis = nuis,
    #                               mus.0 = mus_under_control(principal.assumption = "nearER",
    #                                                         mix.wts = mix.wts.nSNR,
    #                                                         mu.10 = nuis$mu.10,
    #                                                         kappa.0R = nuis$kappa.0R,
    #                                                         y.bounds = y.bounds),
    #                               wt = weights)
    #
    # effects["nearER.nSCR",] <- plugin(nuis = nuis,
    #                               mus.0 = mus_under_control(principal.assumption = "nearER",
    #                                                         mix.wts = mix.wts.nSCR,
    #                                                         mu.10 = nuis$mu.10,
    #                                                         kappa.0R = nuis$kappa.0R,
    #                                                         y.bounds = y.bounds),
    #                               wt = weights)
    #
    # effects["nearER.rPI",] <- plugin(nuis = nuis,
    #                              mus.0 = mus_under_control(principal.assumption = "nearER",
    #                                                        mix.wts = mix.wts.rPI,
    #                                                        mu.10 = nuis$mu.10,
    #                                                        kappa.0R = nuis$kappa.0R,
    #                                                        y.bounds = y.bounds),
    #                              wt = weights)
    #
    # effects["nearER.rPO",] <- plugin(nuis = nuis,
    #                              mus.0 = mus_under_control(principal.assumption = "nearER",
    #                                                        mix.wts = mix.wts.rPO,
    #                                                        mu.10 = nuis$mu.10,
    #                                                        kappa.0R = nuis$kappa.0R,
    #                                                        y.bounds = y.bounds),
    #                              wt = weights)

    effects["PI",] <- plugin(nuis = nuis,
                             mus.0 = mus_under_control(principal.assumption = "PI",
                                                       kappa.0R = nuis$kappa.0R,
                                                       y.bounds = y.bounds),
                             wt = weights)



    sens.min <- min(sens.params)
    sens.max <- max(sens.params)


    if (sens.type=="GOR") {

        if (is.null(y.bounds)) {
            if (max(nuis$kappa.0R)>1 | min(nuis$kappa.0R)<0) {
                stop("y.bounds is required.")
            } else {
                y.bounds <- 0:1
            }
        }

        effects["PIsens.nSNR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       y.bounds = y.bounds,
                                                                       psi = sens.min),
                                             wt = weights)

        effects["PIsens.nSNR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       y.bounds = y.bounds,
                                                                       psi = sens.max),
                                             wt = weights)

        effects["PIsens.nSCR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       y.bounds = y.bounds,
                                                                       psi = sens.min),
                                             wt = weights)

        effects["PIsens.nSCR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       y.bounds = y.bounds,
                                                                       psi = sens.max),
                                             wt = weights)

        effects["PIsens.rPI.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      y.bounds = y.bounds,
                                                                      psi = sens.min),
                                            wt = weights)

        effects["PIsens.rPI.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      y.bounds = y.bounds,
                                                                      psi = sens.max),
                                            wt = weights)

        effects["PIsens.rPO.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      y.bounds = y.bounds,
                                                                      psi = sens.min),
                                            wt = weights)

        effects["PIsens.rPO.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      y.bounds = y.bounds,
                                                                      psi = sens.max),
                                            wt = weights)

    }



    if (sens.type=="SMDe") {

        effects["PIsens.nSNR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       varsigma.0R = nuis$varsigma.0R,
                                                                       eta = sens.min),
                                             wt = weights)

        effects["PIsens.nSNR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       varsigma.0R = nuis$varsigma.0R,
                                                                       eta = sens.max),
                                             wt = weights)

        effects["PIsens.nSCR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       varsigma.0R = nuis$varsigma.0R,
                                                                       eta = sens.min),
                                             wt = weights)

        effects["PIsens.nSCR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       varsigma.0R = nuis$varsigma.0R,
                                                                       eta = sens.max),
                                             wt = weights)

        effects["PIsens.rPI.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      varsigma.0R = nuis$varsigma.0R,
                                                                      eta = sens.min),
                                            wt = weights)

        effects["PIsens.rPI.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      varsigma.0R = nuis$varsigma.0R,
                                                                      eta = sens.max),
                                            wt = weights)

        effects["PIsens.rPO.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      varsigma.0R = nuis$varsigma.0R,
                                                                      eta = sens.min),
                                            wt = weights)

        effects["PIsens.rPO.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      varsigma.0R = nuis$varsigma.0R,
                                                                      eta = sens.max),
                                            wt = weights)

    }


    if (sens.type=="MR") {

        effects["PIsens.nSNR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       rho = sens.min),
                                             wt = weights)

        effects["PIsens.nSNR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       rho = sens.max),
                                             wt = weights)

        effects["PIsens.nSCR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       rho = sens.min),
                                             wt = weights)

        effects["PIsens.nSCR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       rho = sens.max),
                                             wt = weights)

        effects["PIsens.rPI.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      rho = sens.min),
                                            wt = weights)

        effects["PIsens.rPI.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      rho = sens.max),
                                            wt = weights)

        effects["PIsens.rPO.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      rho = sens.min),
                                            wt = weights)

        effects["PIsens.rPO.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      rho = sens.max),
                                            wt = weights)

    }

    effects
}







effects_all_assumptions_binary.y <- function(nuis,
                                             sens.type,
                                             sens.params,
                                             y.bounds = NULL,
                                             weights,
                                             data,
                                             x.vars) {

    effects <- matrix(NA, ncol = 2, nrow = 17)
    colnames(effects) <- c("cace", "nace")
    rownames(effects) <- c("ER.nSNR", "ER.nSCR", "ER.rPI", "ER.rPO",
                           "nearER.nSNR", "nearER.nSCR", "nearER.rPI", "nearER.rPO",
                           "PI",
                           "PIsens.nSNR.lo", "PIsens.nSNR.hi",
                           "PIsens.nSCR.lo", "PIsens.nSCR.hi",
                           "PIsens.rPI.lo", "PIsens.rPI.hi",
                           "PIsens.rPO.lo", "PIsens.rPO.hi")

    mix.wts.nSNR <- mixture_weights(miss.assumption = "nSNR",
                                    pi.1 = nuis$pi.1,
                                    varpi.10 = nuis$varpi.10,
                                    lambda.0 = nuis$lambda.0)

    mix.wts.nSCR <- mixture_weights(miss.assumption = "nSCR",
                                    pi.1 = nuis$pi.1,
                                    varpi.11 = nuis$varpi.11,
                                    lambda.0 = nuis$lambda.0)

    mix.wts.rPI <- mixture_weights(miss.assumption = "rPI",
                                   pi.1 = nuis$pi.1)

    mix.wts.rPO <- mixture_weights(miss.assumption = "rPO",
                                   pi.1 = nuis$pi.1,
                                   varpi.11 = nuis$varpi.11,
                                   varpi.10 = nuis$varpi.10,
                                   lambda.0 = nuis$lambda.0)

    plugin <- function(nuis, mus.0, wt) {
        c(cace = weighted.mean(nuis$mu.11-mus.0$mu.01, wt * nuis$pi.1),
          nace = weighted.mean(nuis$mu.10-mus.0$mu.00, wt * (1-nuis$pi.1)))
    }


    effects["ER.nSNR",] <- plugin(nuis = nuis,
                                  mus.0 = mus_under_control(principal.assumption = "ER",
                                                            mix.wts = mix.wts.nSNR,
                                                            mu.10 = nuis$mu.10,
                                                            kappa.0R = nuis$kappa.0R,
                                                            x.dat = data,
                                                            x.vars = x.vars,
                                                            pi.1 = nuis$pi.1 * weights),
                                  wt = weights)

    effects["ER.nSCR",] <- plugin(nuis = nuis,
                                  mus.0 = mus_under_control(principal.assumption = "ER",
                                                            mix.wts = mix.wts.nSCR,
                                                            mu.10 = nuis$mu.10,
                                                            kappa.0R = nuis$kappa.0R,
                                                            x.dat = data,
                                                            x.vars = x.vars,
                                                            pi.1 = nuis$pi.1 * weights),
                                  wt = weights)

    effects["ER.rPI",] <- plugin(nuis = nuis,
                                 mus.0 = mus_under_control(principal.assumption = "ER",
                                                           mix.wts = mix.wts.rPI,
                                                           mu.10 = nuis$mu.10,
                                                           kappa.0R = nuis$kappa.0R,
                                                           x.dat = data,
                                                           x.vars = x.vars,
                                                           pi.1 = nuis$pi.1 * weights),
                                 wt = weights)

    effects["ER.rPO",] <- plugin(nuis = nuis,
                                 mus.0 = mus_under_control(principal.assumption = "ER",
                                                           mix.wts = mix.wts.rPO,
                                                           mu.10 = nuis$mu.10,
                                                           kappa.0R = nuis$kappa.0R,
                                                           x.dat = data,
                                                           x.vars = x.vars,
                                                           pi.1 = nuis$pi.1 * weights),
                                 wt = weights)

    effects["nearER.nSNR",] <- plugin(nuis = nuis,
                                      mus.0 = mus_under_control(principal.assumption = "nearER",
                                                                mix.wts = mix.wts.nSNR,
                                                                mu.10 = nuis$mu.10,
                                                                kappa.0R = nuis$kappa.0R,
                                                                y.bounds = y.bounds),
                                      wt = weights)

    effects["nearER.nSCR",] <- plugin(nuis = nuis,
                                      mus.0 = mus_under_control(principal.assumption = "nearER",
                                                                mix.wts = mix.wts.nSCR,
                                                                mu.10 = nuis$mu.10,
                                                                kappa.0R = nuis$kappa.0R,
                                                                y.bounds = y.bounds),
                                      wt = weights)

    effects["nearER.rPI",] <- plugin(nuis = nuis,
                                     mus.0 = mus_under_control(principal.assumption = "nearER",
                                                               mix.wts = mix.wts.rPI,
                                                               mu.10 = nuis$mu.10,
                                                               kappa.0R = nuis$kappa.0R,
                                                               y.bounds = y.bounds),
                                     wt = weights)

    effects["nearER.rPO",] <- plugin(nuis = nuis,
                                     mus.0 = mus_under_control(principal.assumption = "nearER",
                                                               mix.wts = mix.wts.rPO,
                                                               mu.10 = nuis$mu.10,
                                                               kappa.0R = nuis$kappa.0R,
                                                               y.bounds = y.bounds),
                                     wt = weights)

    effects["PI",] <- plugin(nuis = nuis,
                             mus.0 = mus_under_control(principal.assumption = "PI",
                                                       kappa.0R = nuis$kappa.0R,
                                                       y.bounds = y.bounds),
                             wt = weights)



    sens.min <- min(sens.params)
    sens.max <- max(sens.params)


    if (sens.type=="GOR") {

        if (is.null(y.bounds)) {
            if (max(nuis$kappa.0R)>1 | min(nuis$kappa.0R)<0) {
                stop("y.bounds is required.")
            } else {
                y.bounds <- 0:1
            }
        }

        effects["PIsens.nSNR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       y.bounds = y.bounds,
                                                                       psi = sens.min),
                                             wt = weights)

        effects["PIsens.nSNR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       y.bounds = y.bounds,
                                                                       psi = sens.max),
                                             wt = weights)

        effects["PIsens.nSCR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       y.bounds = y.bounds,
                                                                       psi = sens.min),
                                             wt = weights)

        effects["PIsens.nSCR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       y.bounds = y.bounds,
                                                                       psi = sens.max),
                                             wt = weights)

        effects["PIsens.rPI.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      y.bounds = y.bounds,
                                                                      psi = sens.min),
                                            wt = weights)

        effects["PIsens.rPI.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      y.bounds = y.bounds,
                                                                      psi = sens.max),
                                            wt = weights)

        effects["PIsens.rPO.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      y.bounds = y.bounds,
                                                                      psi = sens.min),
                                            wt = weights)

        effects["PIsens.rPO.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-GOR",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      y.bounds = y.bounds,
                                                                      psi = sens.max),
                                            wt = weights)

    }



    if (sens.type=="SMDe") {

        effects["PIsens.nSNR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       varsigma.0R = nuis$varsigma.0R,
                                                                       eta = sens.min),
                                             wt = weights)

        effects["PIsens.nSNR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       varsigma.0R = nuis$varsigma.0R,
                                                                       eta = sens.max),
                                             wt = weights)

        effects["PIsens.nSCR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       varsigma.0R = nuis$varsigma.0R,
                                                                       eta = sens.min),
                                             wt = weights)

        effects["PIsens.nSCR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       varsigma.0R = nuis$varsigma.0R,
                                                                       eta = sens.max),
                                             wt = weights)

        effects["PIsens.rPI.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      varsigma.0R = nuis$varsigma.0R,
                                                                      eta = sens.min),
                                            wt = weights)

        effects["PIsens.rPI.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      varsigma.0R = nuis$varsigma.0R,
                                                                      eta = sens.max),
                                            wt = weights)

        effects["PIsens.rPO.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      varsigma.0R = nuis$varsigma.0R,
                                                                      eta = sens.min),
                                            wt = weights)

        effects["PIsens.rPO.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-SMDe",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      varsigma.0R = nuis$varsigma.0R,
                                                                      eta = sens.max),
                                            wt = weights)

    }


    if (sens.type=="MR") {

        effects["PIsens.nSNR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       rho = sens.min),
                                             wt = weights)

        effects["PIsens.nSNR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                       mix.wts = mix.wts.nSNR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       rho = sens.max),
                                             wt = weights)

        effects["PIsens.nSCR.lo",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       rho = sens.min),
                                             wt = weights)

        effects["PIsens.nSCR.hi",] <- plugin(nuis = nuis,
                                             mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                       mix.wts = mix.wts.nSCR,
                                                                       kappa.0R = nuis$kappa.0R,
                                                                       rho = sens.max),
                                             wt = weights)

        effects["PIsens.rPI.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      rho = sens.min),
                                            wt = weights)

        effects["PIsens.rPI.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                      mix.wts = mix.wts.rPI,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      rho = sens.max),
                                            wt = weights)

        effects["PIsens.rPO.lo",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      rho = sens.min),
                                            wt = weights)

        effects["PIsens.rPO.hi",] <- plugin(nuis = nuis,
                                            mus.0 = mus_under_control(principal.assumption = "PIsens-MR",
                                                                      mix.wts = mix.wts.rPO,
                                                                      kappa.0R = nuis$kappa.0R,
                                                                      rho = sens.max),
                                            wt = weights)

    }

    effects
}
