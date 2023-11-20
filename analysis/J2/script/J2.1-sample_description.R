
dat <- readRDS(here::here("analysis", "data", "jdat.rds"))

X.vars <- c("age", "sex", "race", "edu", "marital", "hh.kids", "hh.income",
            "econ.hard", "occu", "wks.unemp", "part.motiv", "seek.motiv",
            "seek.effi", "assertive", "depress")

dat$treat <- factor(dat$treat, levels = 0:1, labels = c("Control", "Treatment"))
dat$complier <- factor(dat$complier, levels = c(1,0), labels = c("Complier", "Noncomplier"))


tab.a <- tableone::CreateTableOne(vars = X.vars,
                                  strata = "treat",
                                  data = dat,
                                  addOverall = TRUE,
                                  test = FALSE)
tab.a

tab.b <- tableone::CreateTableOne(vars = X.vars,
                                  strata = "complier",
                                  data = dat[dat$treat=="Treatment",],
                                  addOverall = FALSE,
                                  test = FALSE)
tab.b

tab1.latex <- cbind(print(tab.a, printToggle = FALSE, noSpaces = TRUE),
                    print(tab.b, printToggle = FALSE, noSpaces = TRUE))
tab1.latex  <- xtable::xtable(tab1.latex)

saveRDS(tab1.latex,
        file = here::here("analysis", "outputs", "J2", "latex-table1.rds"))
