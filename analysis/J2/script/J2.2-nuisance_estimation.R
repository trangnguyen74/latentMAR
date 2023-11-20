
dat <- readRDS(here::here("analysis", "data", "jdat.rds"))

X.vars <- c("age", "sex", "race", "edu", "marital", "hh.kids", "hh.income",
            "econ.hard", "occu", "wks.unemp", "part.motiv", "seek.motiv",
            "seek.effi", "assertive", "depress")

C.form <- paste("complier ~", paste(X.vars, collapse = " + "))
R.form <- paste("response ~", paste(X.vars, collapse = " + "))
Y.form <- paste("y ~",        paste(X.vars, collapse = " + "))

dat1 <- dat[dat$treat==1,]
dat0 <- dat[dat$treat==0,]

c.mod <- glm(C.form, data = dat1, family = quasibinomial)

r1c.mod <- glm(R.form, data = dat1[dat1$complier==1,], family = quasibinomial)
r1n.mod <- glm(R.form, data = dat1[dat1$complier==0,], family = quasibinomial)
r0.mod  <- glm(R.form, data = dat0,                    family = quasibinomial)

y1c.mod <- glm(Y.form, data = dat1[dat1$complier==1,], family = quasibinomial)
y1n.mod <- glm(Y.form, data = dat1[dat1$complier==0,], family = quasibinomial)
y0.mod  <- glm(Y.form, data = dat0,                    family = quasibinomial)

pi.c <- predict(c.mod, newdata = dat, type ="response")

varpi.1c <- predict(r1c.mod, newdata = dat, type = "response")
varpi.1n <- predict(r1n.mod, newdata = dat, type = "response")
lambda.0 <- predict(r0.mod,  newdata = dat, type = "response")

mu.1c   <- predict(y1c.mod, newdata = dat, type = "response")
mu.1n   <- predict(y1n.mod, newdata = dat, type = "response")
kappa.0R <- predict(y0.mod,  newdata = dat, type = "response")


saveRDS(list(pi.1 = pi.c,
             varpi.11 = varpi.1c,
             varpi.10 = varpi.1n,
             lambda.0 = lambda.0,
             mu.11 = mu.1c,
             mu.10 = mu.1n,
             kappa.0R = kappa.0R),
        file = here::here("analysis", "outputs", "J2", "nuisance.rds"))
