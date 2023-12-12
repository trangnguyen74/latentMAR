
nuis <- readRDS(here::here("analysis", "outputs", "EC", "nuisance.rds"))
y.bounds <- c(1,6)


# mixture weights under different specific missingness assumptions

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


# mu.10(X) and mu.11(X) under different principal identification assumptions
# combined with different specific missingness assumptions

mus.0c.ER.nSNR <- mus_under_control(principal.assumption = "ER",
                                    mix.wts = mix.wts.nSNR,
                                    mu.10 = nuis$mu.10,
                                    kappa.0R = nuis$kappa.0R)

mus.0c.ER.nSCR <- mus_under_control(principal.assumption = "ER",
                                    mix.wts = mix.wts.nSCR,
                                    mu.10 = nuis$mu.10,
                                    kappa.0R = nuis$kappa.0R)

mus.0c.ER.rPI <- mus_under_control(principal.assumption = "ER",
                                   mix.wts = mix.wts.rPI,
                                   mu.10 = nuis$mu.10,
                                   kappa.0R = nuis$kappa.0R)

mus.0c.ER.rPO <- mus_under_control(principal.assumption = "ER",
                                   mix.wts = mix.wts.rPO,
                                   mu.10 = nuis$mu.10,
                                   kappa.0R = nuis$kappa.0R)

mus.0c.PI <- mus_under_control(principal.assumption = "PI",
                               kappa.0R = nuis$kappa.0R)

mus.0c.PIsensSMDe.lo.nSNR <- mus_under_control(principal.assumption = "PIsens-SMDe",
                                               mix.wts = mix.wts.nSNR,
                                               kappa.0R = nuis$kappa.0R,
                                               varsigma.0R = nuis$varsigma.0R,
                                               eta = -.5)

mus.0c.PIsensSMDe.hi.nSNR <- mus_under_control(principal.assumption = "PIsens-SMDe",
                                               mix.wts = mix.wts.nSNR,
                                               kappa.0R = nuis$kappa.0R,
                                               varsigma.0R = nuis$varsigma.0R,
                                               eta = .5)

mus.0c.PIsensSMDe.lo.nSCR <- mus_under_control(principal.assumption = "PIsens-SMDe",
                                               mix.wts = mix.wts.nSCR,
                                               kappa.0R = nuis$kappa.0R,
                                               varsigma.0R = nuis$varsigma.0R,
                                               eta = -.5)

mus.0c.PIsensSMDe.hi.nSCR <- mus_under_control(principal.assumption = "PIsens-SMDe",
                                               mix.wts = mix.wts.nSCR,
                                               kappa.0R = nuis$kappa.0R,
                                               varsigma.0R = nuis$varsigma.0R,
                                               eta = .5)

mus.0c.PIsensSMDe.lo.rPI <- mus_under_control(principal.assumption = "PIsens-SMDe",
                                              mix.wts = mix.wts.rPI,
                                              kappa.0R = nuis$kappa.0R,
                                              varsigma.0R = nuis$varsigma.0R,
                                              eta = -.5)

mus.0c.PIsensSMDe.hi.rPI <- mus_under_control(principal.assumption = "PIsens-SMDe",
                                              mix.wts = mix.wts.rPI,
                                              kappa.0R = nuis$kappa.0R,
                                              varsigma.0R = nuis$varsigma.0R,
                                              eta = .5)

mus.0c.PIsensSMDe.lo.rPO <- mus_under_control(principal.assumption = "PIsens-SMDe",
                                              mix.wts = mix.wts.rPO,
                                              kappa.0R = nuis$kappa.0R,
                                              varsigma.0R = nuis$varsigma.0R,
                                              eta = -.5)

mus.0c.PIsensSMDe.hi.rPO <- mus_under_control(principal.assumption = "PIsens-SMDe",
                                              mix.wts = mix.wts.rPO,
                                              kappa.0R = nuis$kappa.0R,
                                              varsigma.0R = nuis$varsigma.0R,
                                              eta = .5)


# effect estimates from the plug-in estimator

eff.ER.nSNR <- plugin_estimate(pi.1 = nuis$pi.1,
                               mu.11 = nuis$mu.11,
                               mu.10 = nuis$mu.10,
                               mu.01 = mus.0c.ER.nSNR$mu.01,
                               mu.00 = mus.0c.ER.nSNR$mu.00) * abs(diff(y.bounds))

eff.ER.nSCR <- plugin_estimate(pi.1 = nuis$pi.1,
                               mu.11 = nuis$mu.11,
                               mu.10 = nuis$mu.10,
                               mu.01 = mus.0c.ER.nSCR$mu.01,
                               mu.00 = mus.0c.ER.nSCR$mu.00) * abs(diff(y.bounds))

eff.ER.rPI <- plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.ER.rPI$mu.01,
                              mu.00 = mus.0c.ER.rPI$mu.00) * abs(diff(y.bounds))
eff.ER.rPO <- plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.ER.rPO$mu.01,
                              mu.00 = mus.0c.ER.rPO$mu.00) * abs(diff(y.bounds))

eff.PI <- plugin_estimate(pi.1 = nuis$pi.1,
                          mu.11 = nuis$mu.11,
                          mu.10 = nuis$mu.10,
                          mu.01 = mus.0c.PI$mu.01,
                          mu.00 = mus.0c.PI$mu.00) * abs(diff(y.bounds))

eff.PIsensSMDe.nSNR <-
    list(lo = plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.PIsensSMDe.lo.nSNR$mu.01,
                              mu.00 = mus.0c.PIsensSMDe.lo.nSNR$mu.00) * abs(diff(y.bounds)),
         hi = plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.PIsensSMDe.hi.nSNR$mu.01,
                              mu.00 = mus.0c.PIsensSMDe.hi.nSNR$mu.00) * abs(diff(y.bounds)))

eff.PIsensSMDe.nSCR <-
    list(lo = plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.PIsensSMDe.lo.nSCR$mu.01,
                              mu.00 = mus.0c.PIsensSMDe.lo.nSCR$mu.00) * abs(diff(y.bounds)),
         hi = plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.PIsensSMDe.hi.nSCR$mu.01,
                              mu.00 = mus.0c.PIsensSMDe.hi.nSCR$mu.00) * abs(diff(y.bounds)))

eff.PIsensSMDe.rPI <-
    list(lo = plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.PIsensSMDe.lo.rPI$mu.01,
                              mu.00 = mus.0c.PIsensSMDe.lo.rPI$mu.00) * abs(diff(y.bounds)),
         hi = plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.PIsensSMDe.hi.rPI$mu.01,
                              mu.00 = mus.0c.PIsensSMDe.hi.rPI$mu.00) * abs(diff(y.bounds)))

eff.PIsensSMDe.rPO <-
    list(lo = plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.PIsensSMDe.lo.rPO$mu.01,
                              mu.00 = mus.0c.PIsensSMDe.lo.rPO$mu.00) * abs(diff(y.bounds)),
         hi = plugin_estimate(pi.1 = nuis$pi.1,
                              mu.11 = nuis$mu.11,
                              mu.10 = nuis$mu.10,
                              mu.01 = mus.0c.PIsensSMDe.hi.rPO$mu.01,
                              mu.00 = mus.0c.PIsensSMDe.hi.rPO$mu.00) * abs(diff(y.bounds)))


effect.estimates <- mget(c("eff.ER.nSNR",
                           "eff.ER.nSCR",
                           "eff.ER.rPI",
                           "eff.ER.rPO",
                           "eff.PI",
                           "eff.PIsensSMDe.nSNR",
                           "eff.PIsensSMDe.SCR",
                           "eff.PIsensSMDe.rPI",
                           "eff.PIsensSMDe.rPO"))


names(effect.estimates) <- substr(names(effect.estimates), 5, nchar(names(effect.estimates)))

saveRDS(effect.estimates,
        file = here::here("analysis", "outputs", "EC", "point_estimates.rds"))

