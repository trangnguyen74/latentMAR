
# First, we confirm that the functions that will be used in bootstrapping indeed
# yield the same point estimates with the original data

dat <- readRDS(here::here("analysis", "data", "jdat.rds"))

nuis <- j2_nuisance_estimation(dat)

effects <- effects_all_assumptions_binary.y(nuis        = nuis,
                                            sens.type   = "GOR",
                                            sens.params = c(1/2,2),
                                            weights     = 1,
                                            y.bounds    = 0:1,
                                            data = dat,
                                            x.vars = c("age", "sex", "race", "edu", "marital", "hh.kids", "hh.income",
                                                       "econ.hard", "occu", "wks.unemp", "part.motiv", "seek.motiv",
                                                       "seek.effi", "assertive", "depress"))
effects



# Now we bootstrap


samp.size <- nrow(dat)

boot.seed <- 12345
boot.num  <- 999



set.seed(boot.seed)
boot.wt <- samp.size * gtools::rdirichlet(n     = boot.num,
                                          alpha = rep(1, samp.size))


boot.ests <- sapply(1:boot.num, function(z) {

    nuis <- j2_nuisance_estimation(cbind(dat, wt = boot.wt[z,]))

    effects <- tryCatch(effects_all_assumptions_binary.y(nuis        = nuis,
                                                         sens.type   = "GOR",
                                                         sens.params = c(1/2,2),
                                                         y.bounds    = 0:1,
                                                         weights     = boot.wt[z,],
                                                         data = dat,
                                                         x.vars = c("age", "sex", "race", "edu", "marital", "hh.kids", "hh.income",
                                                                    "econ.hard", "occu", "wks.unemp", "part.motiv", "seek.motiv",
                                                                    "seek.effi", "assertive", "depress")),
                        error = function(e) return(NULL))

    effects

}, simplify = FALSE)

boot.ests <- boot.ests[!is.null(boot.ests)]
boot.ests <- simplify2array(boot.ests)

perc.ci <- apply(boot.ests, c(1,2), function(x) quantile(x, prob = c(0.025, 0.975), na.rm = TRUE))

effects <- abind::abind(effects, perc.ci[1,,], perc.ci[2,,], along = 3)
effects <- cbind(effects[,1,], effects[,2,])
colnames(effects) <- c("cace", "cace_2.5%", "cace_97.5%",
                       "nace", "nace_2.5%", "nace_97.5%")


saveRDS(effects,
        file = here::here("analysis", "outputs", "J2", "point_and_ci.rds"))
