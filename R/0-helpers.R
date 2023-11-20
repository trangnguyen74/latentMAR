

########################################
#### .generlized_logit_approximate

#' Approximate generalized logit transformation of a variable
#' @param x Variable to be transformed.
#' @param bounds A 2-length numeric vector with the lower and upper bounds of the variable.
#' @keywords internal


.generalized_logit_approximate <- function(x, bounds) {

    x <- (x - min(bounds)) / abs(diff(bounds))
    x <- x * .99 + .005
    log(x/(1-x))
}





##############################
#### .extract_dispersion_wtd_glm

#' Extract the conditional distribution's dispersion param from a weighted GLM
#' (This function is borrowed from the package PIsens, almost verbatim.) Compute dispersion parameter from weighted models ignoring variance inflation due to weighting. Such variance inflation is desired when the dispersion parameter is used to compute SEs of regression coefficients. It is not desired when the dispersion parameter is used to estimate the conditional variance of the response variable, ie `var(Y|X)`.
#' @param mod The fitted linear model object from `lm` or `glm`
#' @keywords internal

.extract_dispersion_wtd_glm <- function(mod) {

    # MANUALLY compute squared Pearson residuals
    # Remember NOT TO USE mod$residuals or residuals(mod, "pearson")
    sq.pearson.resids <-
        (mod$y - mod$fitted.values)^2 / mod$family$variance(mod$fitted.values)

    weighted.mean(sq.pearson.resids, mod$prior.weights) *
        length(mod$fitted.values) / mod$df.residual
}




##############################
#### .wang_bounded_projection

#' Wang style


.wang_bounded_projection <- function(y.vec, X.mat, w.vec, y.bounds) {

    l <- min(y.bounds)
    h <- max(y.bounds)

    trans.to.01 <- function(x) (x - l) / (h - l)
    trans.back  <- function(x) x * (h-l) + l

    y.vec <- trans.to.01(y.vec)

    .expit <- function(x) 1/ (1 + exp(-x))



    Optimize <- function(objective, startpars, thres = 1e-6, max.step = 2000){

        Diff <- function(x,y) abs(x-y)/(x+thres)
        step <- 0
        value.old <- diff <- thres + 1

        while(diff > thres & value.old > thres & step < max.step){
            step <- step + 1
            opt  <- optim(startpars,objective,control=list(maxit=max.step))
            diff <- Diff(opt$value,value.old)
            value.old <- opt$value
            startpars <- opt$par
            if (step %% 10 == 0) {
                cat("This is the ", step, "th step. The optimum value is ",
                    opt$value," after ",opt$counts[1]," iterations \n",sep="")
            }
        }

        opt <- list(par         = opt$par,
                    convergence = (step < max.step),
                    value       = opt$value)
        return(opt)

    }

    objective <- function(beta) {
        obj <- t(X.mat) %*% ((y.vec - .expit(X.mat %*% beta)) * w.vec)
        return(sum(obj^2))
    }


    startpars <- rep(0,dim(X.mat)[2])
    opt       <- Optimize(objective, startpars)

    projs <- .expit(X.mat %*% opt$par)
    projs <- trans.back(projs)

    opt$projection <- projs

    return(opt)
}
