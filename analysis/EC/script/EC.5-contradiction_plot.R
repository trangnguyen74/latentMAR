
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


r1c.gbm = ps(as.formula(R.form), data = dat1[dat1$type==1,],
             n.trees           = 5000,
             interaction.depth = 1,
             shrinkage         = 0.01,
             stop.method       = "ks.max",
             n.minobsinnode    = 10,
             n.keep            = 1,
             n.grid            = 25,
             ks.exact          = NULL,
             verbose           = FALSE,
             keep.data         = TRUE)
r1n.gbm = ps(as.formula(R.form), data = dat1[dat1$type==0,],
             n.trees           = 20000,
             interaction.depth = 1,
             shrinkage         = 0.01,
             stop.method       = c("es.mean", "ks.max"),
             n.minobsinnode    = 10,
             n.keep            = 1,
             n.grid            = 25,
             ks.exact          = NULL,
             verbose           = FALSE)

varpi.1c.gbm <- predict(r1c.gbm$gbm.obj, newdata = dat, type = "response")
varpi.1n.gbm <- predict(r1n.gbm$gbm.obj, newdata = dat, type = "response")


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


ggplot(data = data.frame(probs = varpi.1n,
                         wts = 1-pi.c),
       aes(x = probs, weight = wts)) +
    geom_histogram(breaks = seq(0,1,.05))

ggplot(data = data.frame(probs = varpi.1c,
                         wts = pi.c),
       aes(x = probs, weight = wts)) +
    geom_histogram(breaks = seq(0,1,.05))

ggplot(data = data.frame(probs = lambda.0),
       aes(x = probs)) +
    geom_histogram(breaks = seq(0,1,.05))
