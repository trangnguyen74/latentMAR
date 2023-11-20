


##############################
#### ec_nuisance_estimation

#' Experience Corps specific: Estimate nuisance functions
#'
#' @param dat Data frame to estimate nuisance function on. Could include a variable \code{wt} for weighting units, which is useful for bootstrapping. If not including \code{wt}, units are equally weighted.
#' @export

ec_nuisance_estimation <- function(dat) {

    if (!any(names(dat)=="wt")) dat$wt <- 1

    y.bounds <- c(1,6)
    y.scale <- abs(diff(y.bounds))
    y.min <- min(y.bounds)
    y.max <- max(y.bounds)

    dat$yin01       <- (dat$y - y.min) / y.scale
    dat$base.ylogit <- .generalized_logit_approximate(dat$base.y, bounds = y.bounds)

    X.vars <- c("cohort", "age", "sex", "race", "educ", "income",
                "major.morbidities", "depress", "base.ylogit")

    C.form <- paste("type ~",     paste(X.vars, collapse = " + "))
    R.form <- paste("response ~", paste(X.vars, collapse = " + "))
    Y.form <- paste("yin01 ~",    paste(X.vars, collapse = " + "))

    dat1 <- dat[dat$treat=="Intervention",]
    dat0 <- dat[dat$treat=="Control",]

    c.mod <- do.call("glm",
                     list(formula = C.form,
                          data    = dat1,
                          family  = "quasibinomial",
                          weights = dat1$wt))

    r1c.mod <- do.call("glm",
                       list(formula = R.form,
                            data    = dat1[dat1$type==1,],
                            family  = "quasibinomial",
                            weights = dat1[dat1$type==1,]$wt))
    r1n.mod <- do.call("glm",
                       list(formula = R.form,
                            data    = dat1[dat1$type==0,],
                            family  = "quasibinomial",
                            weights = dat1[dat1$type==0,]$wt))
    r0.mod  <- do.call("glm",
                       list(formula = R.form,
                            data    = dat0,
                            family  = "quasibinomial",
                            weights = dat0$wt))

    y1c.mod <- do.call("glm",
                       list(formula = Y.form,
                            data    = dat1[dat1$type==1,],
                            family  = "quasibinomial",
                            weights = dat1[dat1$type==1,]$wt))
    y1n.mod <- do.call("glm",
                       list(formula = Y.form,
                            data    = dat1[dat1$type==0,],
                            family  = "quasibinomial",
                            weights = dat1[dat1$type==0,]$wt))
    y0.mod  <- do.call("glm",
                       list(formula = Y.form,
                            data    = dat0,
                            family  = "quasibinomial",
                            weights = dat0$wt))

    pi.c <- predict(c.mod, newdata = dat, type ="response")

    varpi.1c <- predict(r1c.mod, newdata = dat, type = "response")
    varpi.1n <- predict(r1n.mod, newdata = dat, type = "response")
    lambda.0 <- predict(r0.mod,  newdata = dat, type = "response")

    mu.1c   <- predict(y1c.mod, newdata = dat, type = "response") * y.scale + y.min
    mu.1n   <- predict(y1n.mod, newdata = dat, type = "response") * y.scale + y.min
    kappa.0R <- predict(y0.mod,  newdata = dat, type = "response") * y.scale + y.min

    dispers.0R <- .extract_dispersion_wtd_glm(y0.mod)
    varsigma.0R <- sqrt(dispers.0R*(kappa.0R-y.min)*(y.max-kappa.0R))


    list(pi.1 = pi.c,
         varpi.11 = varpi.1c,
         varpi.10 = varpi.1n,
         lambda.0 = lambda.0,
         mu.11 = mu.1c,
         mu.10 = mu.1n,
         kappa.0R = kappa.0R,
         varsigma.0R = varsigma.0R)

}
