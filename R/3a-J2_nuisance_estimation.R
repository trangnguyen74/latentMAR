


##############################
#### j2_nuisance_estimation

#' Jobs II specific: Estimate nuisance functions
#'
#' @param dat Data frame to estimate nuisance function on. Could include a variable \code{wt} for weighting units, which is useful for bootstrapping. If not including \code{wt}, units are equally weighted.
#' @export

j2_nuisance_estimation <- function(dat) {

    if (!any(names(dat)=="wt")) dat$wt <- 1

    X.vars <- c("age", "sex", "race", "edu", "marital", "hh.kids", "hh.income",
                "econ.hard", "occu", "wks.unemp", "part.motiv", "seek.motiv",
                "seek.effi", "assertive", "depress")
    X.cont <- c("age", "econ.hard", "wks.unemp", "part.motiv", "seek.motiv",
                "seek.effi", "assertive", "depress")
    X.square <- paste0("I(", X.cont, "^2)")
    X.sqroot <- paste0("I(sqrt(", X.cont, "))")
    X.all <- c(X.vars, X.square, X.sqroot)

    C.form <- paste("complier ~", paste(X.all, collapse = " + "))
    R.form <- paste("response ~", paste(X.vars, collapse = " + "))
    Y.form <- paste("y ~",        paste(X.vars, collapse = " + "))

    dat1 <- dat[dat$treat==1,]
    dat0 <- dat[dat$treat==0,]

    c.mod <- do.call("glm",
                     list(formula = C.form,
                          data    = dat1,
                          family  = "quasibinomial",
                          weights = dat1$wt))

    r1c.mod <- do.call("glm",
                       list(formula = R.form,
                            data    = dat1[dat1$complier==1,],
                            family  = "quasibinomial",
                            weights = dat1[dat1$complier==1,]$wt))
    r1n.mod <- do.call("glm",
                       list(formula = R.form,
                            data    = dat1[dat1$complier==0,],
                            family  = "quasibinomial",
                            weights = dat1[dat1$complier==0,]$wt))
    r0.mod  <- do.call("glm",
                       list(formula = R.form,
                            data    = dat0,
                            family  = "quasibinomial",
                            weights = dat0$wt))

    y1c.mod <- do.call("glm",
                       list(formula = Y.form,
                            data    = dat1[dat1$complier==1,],
                            family  = "quasibinomial",
                            weights = dat1[dat1$complier==1,]$wt))
    y1n.mod <- do.call("glm",
                       list(formula = Y.form,
                            data    = dat1[dat1$complier==0,],
                            family  = "quasibinomial",
                            weights = dat1[dat1$complier==0,]$wt))
    y0.mod  <- do.call("glm",
                       list(formula = Y.form,
                            data    = dat0,
                            family  = "quasibinomial",
                            weights = dat0$wt))

    pi.c <- predict(c.mod, newdata = dat, type ="response")

    varpi.1c <- predict(r1c.mod, newdata = dat, type = "response")
    varpi.1n <- predict(r1n.mod, newdata = dat, type = "response")
    lambda.0 <- predict(r0.mod,  newdata = dat, type = "response")

    mu.1c   <- predict(y1c.mod, newdata = dat, type = "response")
    mu.1n   <- predict(y1n.mod, newdata = dat, type = "response")
    kappa.0R <- predict(y0.mod,  newdata = dat, type = "response")


    list(pi.1 = pi.c,
         varpi.11 = varpi.1c,
         varpi.10 = varpi.1n,
         lambda.0 = lambda.0,
         mu.11 = mu.1c,
         mu.10 = mu.1n,
         kappa.0R = kappa.0R)
}
