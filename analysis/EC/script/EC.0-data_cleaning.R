
ec.dir <- here::here("analysis", "EC")



####################
#### Get data



base_dat <- haven::read_dta(file.path(ec.dir, "data", "raw",
                                      "baseline_analytic_08_28_2013.dta"))
gen_dat  <- read.csv(file.path(ec.dir, "data", "raw",
                               "finalalldata3r13-red061614v3TG.csv"),
                     header = TRUE)
laq_dat  <- foreign::read.spss(file.path(ec.dir, "data", "raw",
                                         "Final dataset_CACE_8_21_13.sav"),
                               to.data.frame = TRUE,
                               use.value.labels = TRUE,
                               duplicated.value.labels="append")


gds <- base_dat[, paste0("GDS", 1:15)]
gds <- gds - 1
gds[, paste0("GDS", c(1, 5, 7, 11, 13))] <- 1 - gds[, paste0("GDS", c(1, 5, 7, 11, 13))]
base_dat$depress <- rowMeans(gds, na.rm = TRUE) * 15

base_dat$income <- ifelse(base_dat$DEM16 < 3.5, 1,
                          ifelse(base_dat$DEM16 < 5.5, 2,
                                 ifelse(base_dat$DEM16 < 8.5, 3,
                                        NA)))
base_dat$income <- factor(base_dat$income, labels = c("<15K", "15-34K", "35K+"))
base_dat <- base_dat[c("id", "income", "depress")]



gen_dat <- gen_dat[c("ID", "genach1", "genach24")]
colnames(gen_dat) <- c("id", "base.y", "y")
make_numeric <- function(x) {
    x <- ifelse(x==".", NA, x)
    as.numeric(as.character(x))
}
gen_dat$base.y <- make_numeric(gen_dat$base.y)
gen_dat$y      <- make_numeric(gen_dat$y)
rm(make_numeric)



laq_dat <- laq_dat[c("ID", "INTGRP", "cohort0", "age_base", "male_base",
                     "Education", "race_binary", "majmor", "sum24")]
colnames(laq_dat) <- c("id", "treat", "cohort", "age", "sex", "educ", "race",
                       "major.morbidities", "hours")
laq_dat$cohort    <- factor(laq_dat$cohort)
levels(laq_dat$educ) <- c("highschool or less", "some college")
levels(laq_dat$race) <- c("nonBlack", "Black")




dat = dplyr::full_join(laq_dat, base_dat,   by = "id")
dat = dplyr::full_join(dat, gen_dat, by = "id")
rm(base_dat, gen_dat, laq_dat)


####################
#### Subset to analysis sample

# remove 11 individuals who cross over from control to treatment arm
dat <- dat[!(dat$treat=="Control" & !is.na(dat$hours)), ]

# remove 68 individuals assigned to treatment who did not receive the
# intervention at the start (see CONSORT diagram in Gruenewald et al., 2016);
# in the dataset these are individuals without record of volunteering hours
# at 24 months
dat <- dat[dat$treat=="Control" | (dat$treat=="Intervention" & !is.na(dat$hours)),]


dat$type <- ifelse(dat$hours>730, 1, 0)
dat <- dat[,names(dat)!="hours"]

########################################
#### Single imputation of covariates

mice::md.pattern(dat %>% select(-c(type, y)),
                 rotate.names = TRUE)

dat$trio <- ifelse(is.na(dat$type), 2, dat$type)
dat$trio <- factor(dat$trio)

imp <- mice::mice(dat %>% select(-c(id, treat, type)),
                  method = "cart",
                  seed = 33,
                  m = 1)

imp <- mice::complete(imp, 1)


dat <- cbind(dat %>% select(id, treat, type, y, cohort, age, sex, educ, race, depress, base.y),
             imp %>% select(income, major.morbidities))

dat$response <- 1 * (!is.na(dat$y))

saveRDS(dat, file = file.path(ec.dir, "data", "edat.rds"))

