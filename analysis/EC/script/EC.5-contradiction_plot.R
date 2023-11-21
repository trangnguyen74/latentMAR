
dat <- readRDS(here::here("analysis", "EC", "data", "edat.rds"))

y.bounds <- c(1,6)
dat$yin01       <- (dat$y - min(y.bounds)) / abs(diff(y.bounds))
dat$base.ylogit <- .generalized_logit_approximate(dat$base.y, bounds = y.bounds)

X.vars <- c("cohort", "age", "sex", "race", "educ", "income",
            "major.morbidities", "depress", "base.ylogit")

C.form <- paste("type ~",     paste(X.vars, collapse = " + "))
R.form <- paste("response ~", paste(X.vars, collapse = " + "))

dat1 <- dat[dat$treat=="Intervention",]
dat0 <- dat[dat$treat=="Control",]

c.mod <- glm(C.form, data = dat1, family = quasibinomial)

r1c.mod <- glm(R.form, data = dat1[dat1$type==1,], family = quasibinomial)
r1n.mod <- glm(R.form, data = dat1[dat1$type==0,], family = quasibinomial)
r0.mod  <- glm(R.form, data = dat0,                family = quasibinomial)


pi.c <- predict(c.mod, newdata = dat, type ="response")

varpi.1c <- predict(r1c.mod, newdata = dat, type = "response")
varpi.1n <- predict(r1n.mod, newdata = dat, type = "response")
lambda.0 <- predict(r0.mod,  newdata = dat, type = "response")


contradiction_plot(miss.assumption = "rER",
                   pi.1 = pi.c,
                   varpi.10 = varpi.1n,
                   lambda.0 = lambda.0,
                   bin.width = .05,
                   x.step = .2)

contradiction_plot(miss.assumption = "SCR",
                   pi.1 = pi.c,
                   varpi.11 = varpi.1c,
                   lambda.0 = lambda.0,
                   bin.width = .05,
                   x.step = .2)
