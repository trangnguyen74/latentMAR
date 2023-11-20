
dat <- readRDS(here::here("analysis", "data", "edat.rds"))

X.vars <- c("age", "sex", "race", "educ", "income",
            "major.morbidities", "depress", "base.y")


tab.a <- tableone::CreateTableOne(vars = X.vars,
                                  strata = "treat",
                                  data = dat,
                                  addOverall = TRUE,
                                  test = FALSE)
tab.a

tab.b <- tableone::CreateTableOne(vars = X.vars,
                                  strata = "type",
                                  data = dat[dat$treat=="Intervention",],
                                  addOverall = FALSE,
                                  test = FALSE)
tab.b

tab1.latex <- cbind(print(tab.a, printToggle = FALSE, noSpaces = TRUE),
                    print(tab.b, printToggle = FALSE, noSpaces = TRUE))
tab1.latex  <- xtable::xtable(tab1.latex)

saveRDS(tab1.latex,
        file = here::here("analysis", "outputs", "EC", "latex-table1.rds"))
