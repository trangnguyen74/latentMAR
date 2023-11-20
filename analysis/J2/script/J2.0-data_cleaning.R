
# The data here are from JOBS-II study,
# available from https://www.icpsr.umich.edu/web/ICPSR/studies/2739

# We use data file 02739-0011-Data.por within data item DS11, requested from
# https://www.icpsr.umich.edu/web/ICPSR/studies/2739/datadocumentation


library(tidyverse)


########################################
#### Get data

dat <- as.data.frame(
    memisc::as.data.set(
        memisc::spss.portable.file(
            here::here("analysis",
                       "data",
                       "J2",
                       "02739-0011-Data.por")))
) %>%
    filter(!is.na(V9106), RISK01=="high") %>%

    mutate(treat      = ifelse(V9106=="EXP", 1, 0),
           complier   = 1 * (V9107=="SHOW"),

           age        = V9002,
           sex        = V9001,
           race       = V1424,
           edu        = V1531,
           marital    = V1407,
           hh.kids    = V1408,
           hh.income  = V1532,
           econ.hard  = V1503,
           occu       = V1401A,
           wks.unemp  = V1405,
           part.motiv = V9101,
           seek.motiv = V1535,
           seek.effi  = V1510,
           assertive  = V1522,
           depress    = V1518,
           y          = V3461,
           earning    = V3404,
           id         = V1
    ) %>%
    select(id, treat, complier,
           age, sex, race, edu, marital, hh.kids,
           hh.income, econ.hard, occu, wks.unemp,
           part.motiv, seek.motiv, seek.effi, assertive,
           depress,
           y, earning)


########################################
#### Recode outcome

dat <- dat %>%
    mutate(y = ifelse(!is.na(earning) & earning>0, 1, y)) %>%
    select(-earning)


########################################
#### Recode covariates

dat <- dat %>%
    # recode education to 3 levels
    mutate(edu = ifelse(edu=="gradwk", 4, edu)) %>%
    # truncate 18 extreme values of wks.unemp
    mutate(wks.unemp = ifelse(wks.unemp>52, 52, wks.unemp)) %>%

    # truncate 16 large values of hh.kids
    mutate(hh.kids = ifelse(hh.kids>5, 5, hh.kids)) %>%

    # dichotomize race (collapsing small categories)
    mutate(race = ifelse(race=="WHITE", "white", "nonwhite"),
           race = factor(race)) %>%

    # collapse marital to three categories
    mutate(marital = ifelse(marital=="NEVMARR",
                            "nevmarr",
                            ifelse(marital=="MARRIED",
                                   "married",
                                   "divsepwid")),
           marital = factor(marital, levels = c("nevmarr",
                                                "married",
                                                "divsepwid"))) %>%

    mutate(sex = factor(sex, levels = c(0, 1), labels = c("male", "female")))

dat$edu <- factor(dat$edu, labels = c("lt-hs", "highsc", "somcol", "degree"))

levels(dat$occu) <- c("professional",
                      "managerial",
                      "clerical",
                      "sales",
                      "crafts/foremen",
                      "operative",
                      "labor/service")


########################################
#### Single imputation of covariates

mice::md.pattern(dat %>% select(-c(complier, y)),
                 rotate.names = TRUE)

dat$trio <- ifelse(is.na(dat$complier), 2, dat$complier)
dat$trio <- factor(dat$trio)

imp <- mice::mice(dat %>% select(-c(id, treat, complier)),
                  method = "cart",
                  seed = 33,
                  m = 1)

imp <- mice::complete(imp, 1)


dat <- cbind(dat %>% select(id, treat, complier, y, age, sex),
             imp %>% select(-c(age, sex, y, trio)))

dat$response <- 1 * (!is.na(dat$y))


saveRDS(dat, file = here::here("analysis", "data", "jdat.rds"))

