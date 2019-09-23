# boosted algorithms
SL.xgboost2 <- function(..., max_depth = 2){
	SL.xgboost(..., max_depth = max_depth)
}
SL.xgboost4 <- function(..., max_depth = 4){
	SL.xgboost(..., max_depth = max_depth)
}
SL.xgboost6 <- function(..., max_depth = 6){
	SL.xgboost(..., max_depth = max_depth)
}
SL.xgboost8 <- function(..., max_depth = 8){
	SL.xgboost(..., max_depth = max_depth)
}
descr_SL.xgboost <- "boosted regression trees with maximum depth of "
descr_SL.xgboost2 <- paste0(descr_SL.xgboost, 2)
descr_SL.xgboost4 <- paste0(descr_SL.xgboost, 4)
descr_SL.xgboost6 <- paste0(descr_SL.xgboost, 6)
descr_SL.xgboost8 <- paste0(descr_SL.xgboost, 8)

# random forests
SL.ranger.imp <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
    write.forest = TRUE, probability = family$family == "binomial", 
    min.node.size = ifelse(family$family == "gaussian", 5, 1), 
    replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), 
    num.threads = 1, verbose = TRUE, ...) {
    SuperLearner:::.SL.require("ranger")
    if (family$family == "binomial") {
        Y = as.factor(Y)
    }
    if (is.matrix(X)) {
        X = data.frame(X)
    }
    fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X), 
        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
        replace = replace, sample.fraction = sample.fraction, 
        case.weights = obsWeights, write.forest = write.forest, 
        probability = probability, num.threads = num.threads, 
        verbose = verbose, importance = "impurity")
    pred <- predict(fit, data = newX)$predictions
    if (family$family == "binomial") {
        pred = pred[, "1"]
    }
    fit <- list(object = fit, verbose = verbose)
    class(fit) <- c("SL.ranger")
    out <- list(pred = pred, fit = fit)
    return(out)
}
SL.ranger.reg <- function(..., X, mtry = floor(sqrt(ncol(X)))){
	SL.ranger.imp(..., mtry = mtry)
}

SL.ranger.small <- function(..., X, mtry = floor(sqrt(ncol(X)) * 1/2)){
	SL.ranger.imp(..., mtry = mtry)
}

SL.ranger.large <- function(..., X, mtry = floor(sqrt(ncol(X)) * 2)){
	SL.ranger.imp(..., mtry = mtry)
}
descr_SL.ranger <- "random forest with mtry equal to "
descr_SL.ranger.reg <- paste0(descr_SL.ranger, "square root of number of predictors")
descr_SL.ranger.small <- paste0(descr_SL.ranger, "one-half times square root of number of predictors")
descr_SL.ranger.large <- paste0(descr_SL.ranger, "two times square root of number of predictors")

# lasso
SL.glmnet.50 <- function(..., alpha = 0.5){
	SL.glmnet(..., alpha = alpha)
}
SL.glmnet.25 <- function(..., alpha = 0.25){
	SL.glmnet(..., alpha = alpha)
}
SL.glmnet.75 <- function(..., alpha = 0.75){
	SL.glmnet(..., alpha = alpha)
}

descr_SL.glmnet <- "GLMNET with lambda selected by CV and alpha equal to "
descr_SL.glmnet.50 <- paste0(descr_SL.glmnet, "0.5")
descr_SL.glmnet.25 <- paste0(descr_SL.glmnet, "0.25")
descr_SL.glmnet.75 <- paste0(descr_SL.glmnet, "0.75")
descr_SL.glmnet <- "GLMNET with lambda selected by CV and alpha equal to 0"

descr_SL.mean <- "intercept only regression"
descr_SL.glm <- "main terms generalized linear model"

default_library <- c("SL.mean","SL.xgboost2", "SL.xgboost4", "SL.xgboost6", "SL.xgboost8",
                     "SL.ranger.small", "SL.ranger.reg", "SL.ranger.large", 
                     "SL.glmnet", "SL.glmnet.25", "SL.glmnet.50", "SL.glmnet.75")

default_library_reduced <- c("SL.mean", "SL.glm")

#' Temporary fix for convex combination method mean squared error
#' Relative to existing implementation, we reduce the tolerance at which 
#' we declare predictions from a given algorithm the same as another
tmp_method.CC_LS <- function () 
{
    computeCoef = function(Z, Y, libraryNames, verbose, obsWeights, 
        errorsInLibrary = NULL, ...) {
        cvRisk <- apply(Z, 2, function(x) mean(obsWeights * (x - 
            Y)^2))
        names(cvRisk) <- libraryNames
        compute <- function(x, y, wt = rep(1, length(y))) {
            wX <- sqrt(wt) * x
            wY <- sqrt(wt) * y
            D <- crossprod(wX)
            d <- crossprod(wX, wY)
            A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
            bvec <- c(1, rep(0, ncol(wX)))
            fit <- tryCatch({quadprog::solve.QP(Dmat = D, dvec = d, Amat = A, 
                bvec = bvec, meq = 1)
          }, error = function(e){
            out <- list()
            class(out) <- "error"
            out
          })
            invisible(fit)
        }
        modZ <- Z
        naCols <- which(apply(Z, 2, function(z) {
            all(z == 0)
        }))
        anyNACols <- length(naCols) > 0
        if (anyNACols) {
            warning(paste0(paste0(libraryNames[naCols], collapse = ", "), 
                " have NAs.", "Removing from super learner."))
        }
        tol <- 4
        dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
        anyDupCols <- length(dupCols) > 0
        if (anyDupCols) {
            warning(paste0(paste0(libraryNames[dupCols], collapse = ", "), 
                " are duplicates of previous learners.", " Removing from super learner."))
        }
        if (anyDupCols | anyNACols) {
            rmCols <- unique(c(naCols, dupCols))
            modZ <- Z[, -rmCols, drop = FALSE]
        }
        fit <- compute(x = modZ, y = Y, wt = obsWeights)
        if(class(fit) != "error"){
          coef <- fit$solution
        }else{
          coef <- rep(0, ncol(Z))
          coef[which.min(cvRisk)] <- 1
        }
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] = 0
        }
        if(class(fit) != "error"){
          if (anyDupCols | anyNACols) {
              ind <- c(seq_along(coef), rmCols - 0.5)
              coef <- c(coef, rep(0, length(rmCols)))
              coef <- coef[order(ind)]
          }
          coef[coef < 1e-04] <- 0
          coef <- coef/sum(coef)
        }
        if (!sum(coef) > 0) 
            warning("All algorithms have zero weight", call. = FALSE)
        list(cvRisk = cvRisk, coef = coef, optimizer = fit)
    }
    computePred = function(predY, coef, ...) {
        predY %*% matrix(coef)
    }
    out <- list(require = "quadprog", computeCoef = computeCoef, 
        computePred = computePred)
    invisible(out)
}


#' Temporary fix for convex combination method negative log-likelihood loss
#' Relative to existing implementation, we reduce the tolerance at which 
#' we declare predictions from a given algorithm the same as another.
#' Note that because of the way \code{SuperLearner} is structure, one needs to 
#' install the optimization software separately.
tmp_method.CC_nloglik <- function () 
{
    computePred = function(predY, coef, control, ...) {
        if (sum(coef != 0) == 0) {
            stop("All metalearner coefficients are zero, cannot compute prediction.")
        }
        stats::plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*% 
            matrix(coef[coef != 0]))
    }
    computeCoef = function(Z, Y, libraryNames, obsWeights, control, 
        verbose, ...) {
        tol <- 4
        dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
        anyDupCols <- length(dupCols) > 0
        modZ <- Z
        if (anyDupCols) {
            warning(paste0(paste0(libraryNames[dupCols], collapse = ", "), 
                " are duplicates of previous learners.", " Removing from super learner."))
            modZ <- modZ[, -dupCols, drop = FALSE]
        }
        modlogitZ <- trimLogit(modZ, control$trimLogit)
        logitZ <- trimLogit(Z, control$trimLogit)
        cvRisk <- apply(logitZ, 2, function(x) -sum(2 * obsWeights * 
            ifelse(Y, stats::plogis(x, log.p = TRUE), stats::plogis(x, log.p = TRUE, 
                lower.tail = FALSE))))
        names(cvRisk) <- libraryNames
        obj_and_grad <- function(y, x, w = NULL) {
            y <- y
            x <- x
            function(beta) {
                xB <- x %*% cbind(beta)
                loglik <- y * stats::plogis(xB, log.p = TRUE) + (1 - 
                  y) * stats::plogis(xB, log.p = TRUE, lower.tail = FALSE)
                if (!is.null(w)) 
                  loglik <- loglik * w
                obj <- -2 * sum(loglik)
                p <- stats::plogis(xB)
                grad <- if (is.null(w)) 
                  2 * crossprod(x, cbind(p - y))
                else 2 * crossprod(x, w * cbind(p - y))
                list(objective = obj, gradient = grad)
            }
        }
        lower_bounds = rep(0, ncol(modZ))
        upper_bounds = rep(1, ncol(modZ))
        if (anyNA(cvRisk)) {
            upper_bounds[is.na(cvRisk)] = 0
        }
        r <- tryCatch({nloptr::nloptr(x0 = rep(1/ncol(modZ), ncol(modZ)), 
            eval_f = obj_and_grad(Y, modlogitZ), lb = lower_bounds, 
            ub = upper_bounds, eval_g_eq = function(beta) (sum(beta) - 
                1), eval_jac_g_eq = function(beta) rep(1, length(beta)), 
            opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08))
        }, error = function(e){
          out <- list()
          class(out) <- "error"
          out
        })
        if (r$status < 1 || r$status > 4) {
            warning(r$message)
        }
        if(class(r) != "error"){
          coef <- r$solution
        }else{
          coef <- rep(0, ncol(Z))
          coef[which.min(cvRisk)] <- 1
        }
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] <- 0
        }
        if (anyDupCols) {
            ind <- c(seq_along(coef), dupCols - 0.5)
            coef <- c(coef, rep(0, length(dupCols)))
            coef <- coef[order(ind)]
        }
        coef[coef < 1e-04] <- 0
        coef <- coef/sum(coef)
        out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
        return(out)
    }
    list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}