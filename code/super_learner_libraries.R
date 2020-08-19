# boosted algorithms

SL.xgboost.corrected <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), id, ntrees = 1000,
    max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(),
    nthread = 1, verbose = 0, save_period = NULL, ...)
{
    SuperLearner:::.SL.require("xgboost")
    if (packageVersion("xgboost") < 0.6)
        stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
    if (!is.matrix(X)) {
        X = model.matrix(~. - 1, X)
    }
    xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
    if (family$family == "gaussian") {
        model = xgboost::xgboost(data = xgmat, objective = "reg:squarederror",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, nthread = nthread,
            params = params, save_period = save_period)
    }
    if (family$family == "binomial") {
        model = xgboost::xgboost(data = xgmat, objective = "binary:logistic",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, nthread = nthread,
            params = params, save_period = save_period)
    }
    if (family$family == "multinomial") {
        model = xgboost::xgboost(data = xgmat, objective = "multi:softmax",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, num_class = length(unique(Y)),
            nthread = nthread, params = params, save_period = save_period)
    }
    if (!is.matrix(newX)) {
        newX = model.matrix(~. - 1, newX)
    }
    pred = predict(model, newdata = newX)
    fit = list(object = model)
    class(fit) = c("SL.xgboost")
    out = list(pred = pred, fit = fit)
    return(out)
}

SL.xgboost.2 <- function(..., max_depth = 2){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost.4 <- function(..., max_depth = 4){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost.6 <- function(..., max_depth = 6){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost.8 <- function(..., max_depth = 8){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
descr_SL.xgboost <- "boosted regression trees with maximum depth of "
descr_SL.xgboost.2 <- paste0(descr_SL.xgboost, 2)
descr_SL.xgboost.4 <- paste0(descr_SL.xgboost, 4)
descr_SL.xgboost.6 <- paste0(descr_SL.xgboost, 6)
descr_SL.xgboost.8 <- paste0(descr_SL.xgboost, 8)

# random forests
SL.ranger.imp <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), num.trees = 500, mtry = floor(sqrt(ncol(X))),
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
	SL.ranger.imp(..., X = X, mtry = mtry)
}

SL.ranger.small <- function(..., X, mtry = floor(sqrt(ncol(X)) * 1/2)){
	SL.ranger.imp(..., X  = X, mtry = mtry)
}

SL.ranger.large <- function(..., X, mtry = floor(sqrt(ncol(X)) * 2)){
	SL.ranger.imp(..., X = X, mtry = mtry)
}
descr_SL.ranger.imp <- "random forest with `mtry` equal to "
descr_SL.ranger.reg <- paste0(descr_SL.ranger.imp, "square root of number of predictors")
descr_SL.ranger.small <- paste0(descr_SL.ranger.imp, "one-half times square root of number of predictors")
descr_SL.ranger.large <- paste0(descr_SL.ranger.imp, "two times square root of number of predictors")

# function used to do smarter CV for glmnet
get_fold_id <- function(Y){
  fold_id <- rep(0, length(Y))
  wiY0 <- which(Y == 0)
  wiY1 <- which(Y == 1)
  #if <4 cases, no cv
  if(length(wiY1) == 4){
    #if exactly 4 cases, 4-fold cv
    #1 case per fold
    fold <- 1:4
    fold_id[sample(wiY1)] <- fold
    fold_id[sample(wiY0)] <- rep(fold, length = length(wiY0))
  }else{
    #if >=5 cases, 5 fold cv
    #cases split as evenly as possible
    fold <- 1:5
    fold_id[sample(wiY1)] <- rep(fold, length = length(wiY1))
    fold_id[sample(wiY0)] <- rep(fold, length = length(wiY0))
  }
  return(fold_id)
}


# function to have more robust behavior in SL.glmnet
SL.glmnet.0 <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), id, alpha = 1, nfolds = 5,
    nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
    SuperLearner:::.SL.require("glmnet")
    if (!is.matrix(X)) {
        X <- model.matrix(~-1 + ., X)
        newX <- model.matrix(~-1 + ., newX)
    }

    if(family$family == "binomial"){
        fold_id <- get_fold_id(Y)
        nfolds <- max(fold_id)
        if(nfolds != 0){
            fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                lambda = NULL, type.measure = loss, nfolds = nfolds,
                foldid = fold_id, family = family$family, alpha = alpha, nlambda = nlambda,
                ...)
            pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                "lambda.min", "lambda.1se"))
            fit <- list(object = fitCV, useMin = useMin)
            class(fit) <- "SL.glmnet"
        }else{
            # if fewer than 3 cases, just use mean
            meanY <- weighted.mean(Y, w = obsWeights)
            pred <- rep.int(meanY, times = nrow(newX))
            fit <- list(object = meanY)
            out <- list(pred = pred, fit = fit)
            class(fit) <- c("SL.mean")
        }
    }else{
        fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                lambda = NULL, type.measure = loss, nfolds = nfolds,
                family = family$family, alpha = alpha, nlambda = nlambda,
                ...)
            pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                "lambda.min", "lambda.1se"))
            fit <- list(object = fitCV, useMin = useMin)
            class(fit) <- "SL.glmnet"
    }
    out <- list(pred = pred, fit = fit)
    return(out)
}

# lasso
SL.glmnet.50 <- function(..., alpha = 0.5){
	SL.glmnet.0(..., alpha = alpha)
}
SL.glmnet.25 <- function(..., alpha = 0.25){
	SL.glmnet.0(..., alpha = alpha)
}
SL.glmnet.75 <- function(..., alpha = 0.75){
	SL.glmnet.0(..., alpha = alpha)
}

descr_SL.glmnet <- "elastic net with $\\lambda$ selected by 5-fold CV and $\\alpha$ equal to "
descr_SL.glmnet.50 <- paste0(descr_SL.glmnet, "0.5")
descr_SL.glmnet.25 <- paste0(descr_SL.glmnet, "0.25")
descr_SL.glmnet.75 <- paste0(descr_SL.glmnet, "0.75")
descr_SL.glmnet.0 <- "elastic net with $\\lambda$ selected by CV and $\\alpha$ equal to 0"

descr_SL.mean <- "intercept only regression"
descr_SL.glm <- "main terms generalized linear model"

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

# based on user-inputted options, make a vector describing the super learner library
make_sl_library_vector <- function(opts){
    default_library <- NULL
    # check if rf is requested
    if ("rf" %in% opts$learners) {
      if (opts$cvtune) {
        default_library <- c(default_library, "SL.ranger.small", "SL.ranger.reg", "SL.ranger.large")
      } else {
        if (length(opts$learners) == 1) {
            default_library <- "SL.ranger.reg"
        } else {
            default_library <- c(default_library, "SL.ranger.reg")
        }
      }
    }
    # check if xgboost is requested
    if("xgboost" %in% opts$learners){
      if (opts$cvtune) {
        default_library <- c(default_library, "SL.xgboost.2", "SL.xgboost.4", "SL.xgboost.6", "SL.xgboost.8")
      } else {
        if (length(opts$learners) == 1) {
            default_library <- "SL.xgboost.4"
        } else {
            default_library <- c(default_library, "SL.xgboost.4")
        }
      }
    }
    # check if elastic net is requested
    if("lasso" %in% opts$learners){
      if (opts$cvtune) {
        default_library <- c(default_library, "SL.glmnet.0", "SL.glmnet.25", "SL.glmnet.50", "SL.glmnet.75")
      } else {
        if (length(opts$learners) == 1) {
            default_library <- "SL.glmnet.0"
        } else {
            default_library <- c(default_library, "SL.glmnet.0")
        }
      }
    }
    # if fitting a super learner, throw in SL.mean
    if(length(opts$learners) > 1){
      default_library <- c(default_library, "SL.mean")
    }
    return(default_library)
}


#' function to run super learner and cv super learner on a single outcome
#' @param complete_dat the full dataset
#' @param outcome_name String name of outcome
#' @param pred_names Vector of string names of predictor variables
#' @param opts List of options outputted by get_global_options()
#' @param save_dir name of directory to save results to
#' @param fit_name name of fits (defaults to fit_<outcome_name>.rds)
#' @param cv_fit_name name of CV fits (defaults to cvfit_<outcome_name>.rds)
#' @param save_full_object Flag for whether or not to save the full fitted object, or just the fitted values
#' @param outer_folds a set of outer folds for VIM hypothesis testing
#' @param full_fit is it the full fit (TRUE) or a reduced fit (FALSE)?
#' @param SL.library the library to pass to SuperLearner
#' @param ... additional arguments to pass to individual algorithms or the SuperLearner
sl_one_outcome <- function(complete_dat, outcome_name,
                           pred_names,
                           opts,
                           save_dir = "/home/slfits/",
                           fit_name = paste0("fit_", outcome_name, ".rds"),
                           cv_fit_name = paste0("cvfit_", outcome_name, ".rds"),
                           save_full_object = TRUE,
                           outer_folds = rep(1, length(complete_dat[, outcome_name])),
                           full_fit = TRUE,
                           SL.library = "SL.mean",
                           ...){
  # three cases to worry about in terms of how to deal with missing data
  # 1. same_subset requested, but only studying ic80 or only studying ic50-derived outcomes,
  #    in which case, we will ignore your same_subset request because you're only studying
  #    endpoints based entirely on ic50 or entirely on ic80
  # 2. same_subset requested and studying mix of ic50(-derived) and ic80 endpoints, in which
  #    case we will subset down to sequences with both ic50 and ic80
  # 3. no same_subset requested, in which case we use all data available on each outcome
  if(!opts$same_subset | !(("ic80" %in% opts$outcomes | "iip" %in% opts$outcomes) & length(opts$outcomes) > 1)){
    complete_cases_idx <- complete.cases(complete_dat[,c(outcome_name,pred_names)])
  } else {
    complete_cases_idx <- complete.cases(complete_dat)
  }
  if (all(outer_folds == 1)) {
      outer_folds <- outer_folds[complete_cases_idx]
  }
  if (full_fit) {
      outer_bool <- outer_folds == 1
  } else {
      outer_bool <- outer_folds == 2
  }
  # subset data to only complete outcome and pred_names
  dat <- complete_dat[complete_cases_idx, ]
  newdat <- subset(dat, outer_bool)
  pred <- newdat[ , pred_names]
  # grab args
  L <- list(...)
  # make names for the full learner, cv learner, fitted values, cv fitted values, etc.
  learner_name <- gsub("fit_", "learner_", fit_name)
  fitted_name <- gsub("fit_", "fitted_", fit_name)
  cv_fitted_name <- gsub("cvfit_", "cvfitted_", cv_fit_name)
  cv_folds_name <- gsub("cvfitted_", "cvfolds_", gsub("cvfit_", "cvfolds_", cv_fit_name))
  cv_learner_name <- gsub("cvfit_", "cvlearner_", cv_fit_name)
  # unless we do not want to tune using CV nor evaluate performance using CV,
  # we will make a call to super learner
  if (!(opts$cvtune == FALSE & opts$cvperf == FALSE)) {
    fit <- SuperLearner(Y = newdat[ , outcome_name], X = pred, SL.library = SL.library, ...)
    # since cv super learner saves this as output and we need some parsimony later...
    fit$Y <- newdat[ , outcome_name]
    if (save_full_object) {
        saveRDS(fit, file = paste0(save_dir, fit_name))
        saveRDS(fit, file = paste0(save_dir, learner_name))
    }
    if (length(opts$learners) > 1) {
      # save super learner predictions
      saveRDS(fit$SL.predict, file = paste0(save_dir, fitted_name))
      # save super learner weights
      saveRDS(fit$coef, file = paste0(save_dir, "slweights_", fit_name))
    } else if (length(opts$learners) == 1) {
        # save learner predictions (these are cv-selected if cvtune = TRUE)
        saveRDS(fit$library.predict[, which.min(fit$cvRisk)], file = paste0(save_dir, fitted_name))
        # save learner
        if (save_full_object) {
            saveRDS(fit$fitLibrary[[which.min(fit$cvRisk)]]$object, file = paste0(save_dir, learner_name))
        }
    } else {
        # don't save anything
    }
  } else {
    # if we don't want to use CV at all, then use "default" learners
    if (length(opts$learners) > 1) {
        these_learners <- NULL
        if ("rf" %in% opts$learners) {
            these_learners <- c(these_learners, SL.library[grepl("ranger", SL.library)][1])
        }
        if ("lasso" %in% opts$learners) {
            these_learners <- c(these_learners, SL.library[grepl("glmnet", SL.library)][1])
        }
        if ("xgboost" %in% opts$learners) {
            these_learners <- c(these_learners, SL.library[grepl("xgboost", SL.library)][1])
        }
        these_learners <- c(these_learners, "SL.mean")
        fit <- SuperLearner(Y = newdat[ , outcome_name], X = pred, SL.library = these_learners, ...)
        fit$Y <- newdat[ , outcome_name]
        if (save_full_object) {
            saveRDS(fit, file = paste0(save_dir, fit_name))
            saveRDS(fit, file = paste0(save_dir, learner_name))
        }
        # save super learner predictions
        saveRDS(fit$SL.predict, file = paste0(save_dir, fitted_name))
        # save super learner weights
        saveRDS(fit$coef, file = paste0(save_dir, "slweights_", fit_name))
    } else {
        # in this case, directly call the wrapper function
        # note that it is the first instance of the first listed learner
        if (opts$learners[1] == "rf") {
            this_learner <- SL.library[grepl("ranger", SL.library)][1]
        } else if (opts$learners[1] == "lasso") {
            this_learner <- SL.library[grepl("glmnet", SL.library)][1]
        } else {
            this_learner <- SL.library[grepl("xgboost", SL.library)][1]
        }
        fit <- do.call(this_learner, args = c(L[!grepl("cvControl", names(L)) & !grepl("method", names(L))], list(Y = newdat[ , outcome_name], X = pred, newX = pred)))
        saveRDS(fit$pred, file = paste0(save_dir, fitted_name))
        # this will be an object with class native to what the individual learner is
        # i.e., if rf is desired, it'll be ranger object
        if (save_full_object) {
            saveRDS(fit$fit$object, file = paste0(save_dir, fit_name))
            saveRDS(fit$fit$object, file = paste0(save_dir, learner_name))
        }
    }
  }
  # now copy cvControl to innerCvControl to pass to CV.SL
  L$innerCvControl <- list(L$cvControl)
  new_arg_list <- c(list(Y = newdat[ , outcome_name], X = pred, SL.library = SL.library), L)
  if (length(opts$learners) == 1 & opts$cvtune & opts$cvperf) {
    # in this case, we want to report back only results for discrete super learner
    # since we're trying to pick out e.g., the best single random forest fit
    # note that if either opts$cvtune AND/OR opts$cvperf is FALSE then everything we need
    # will already be in the SuperLearner fit object.
    cv_fit <- do.call(CV.SuperLearner, new_arg_list)
    if (save_full_object) {
        saveRDS(cv_fit, file = paste0(save_dir, cv_fit_name))
        saveRDS(cv_fit, file = paste0(save_dir, cv_learner_name))
    }
    saveRDS(cv_fit$discreteSL.predict, file = paste0(save_dir, cv_fitted_name))
    saveRDS(cv_fit$folds, file = paste0(save_dir, cv_folds_name))
  } else if (length(opts$learners) > 1 & opts$cvperf) {
    # if multiple learners, then we are fitting a super learner so need CV superlearner
    # unless cvperf = FALSE
    cv_fit <- do.call(CV.SuperLearner, new_arg_list)
    if (save_full_object) {
        saveRDS(cv_fit, file = paste0(save_dir, cv_fit_name))
        saveRDS(cv_fit, file = paste0(save_dir, cv_learner_name))
    }
    saveRDS(cv_fit$SL.predict, file = paste0(save_dir, cv_fitted_name))
    saveRDS(cv_fit$folds, file = paste0(save_dir, cv_folds_name))
} else {
    # do nothing
}
  return(invisible(NULL))
}
