#' Calculate maximum score
#'
#' Optimizes the maximum score function defined by the data array. The objective
#' function is created using \code{makeScoreObjFun}.
#'
#' @section Optimization methods:
#' The optimization method is not bound to the problem. Any method that can
#' optimize a non-convex, non-smooth function is valid. However, we currently
#' only support the *Differential Evolution* method, as implemented in the
#' package **DEoptim**. In case the user wants to pass parameters to the solver,
#' we provide the parameter \code{optimParams}, which is forwarded to the
#' \code{control} parameter of the function \code{DEoptim::DEoptim}.
#' See \code{DEoptim::DEoptim.control} for more information.
#'
#' @section Sensitivity to parameter choices:
#' Results can be sensitive to parameter choices: changing the random seed or
#' optimization settings (e.g. bounds, \code{optimParams}) can lead to
#' substantially different estimates. We recommend setting a seed, using
#' \code{permuteInvariant = TRUE}, and \code{numRuns > 1} when stability or
#' reproducibility matters.
#'
#' @param dataArray The output of \code{CdataArray}.
#' @param bounds A list with elements \code{$lower}, \code{$upper} which are
#'   vectors defining lower and upper bounds for each variable in the objective
#'   function.
#' @param coefficient1 (optional) The first coefficient of the extended
#'   parameter vector. See \code{makeScoreObjFun} for more details.
#' @param method (optional) A string denoting the optimization method.
#'   Currently, only \code{"DEoptim"} (the default value) is supported.
#' @param optimParams (optional) A list of parameters to be passed to the
#'   optimization routine. Defaults to the empty list.
#'   See the section "Optimization methods" for more information.
#' @param getIneqSat (optional) A boolean indicating whether to include the
#'   \code{$ineqSat} member in the result. Defaults to \code{FALSE}.
#' @param permuteInvariant (optional) Whether to reorder the attributes rows of
#'   the data array before and after the optimization, such that the attribute
#'   whose row has the smallest standard deviation comes first. Note that
#'   setting this to `TRUE` affects only the internal procedures, as the initial
#'   order is restored in the results. Defaults to \code{TRUE}.
#' @param numRuns (optional) How many times to restart the optimization method. Useful with
#'   methods such as *Differential Evolution*, when restarting the method when the
#'   random seed has changed might yield a different result. Defaults to `1`.
#' @param progressUpdate (optional) How often to print progress. Defaults to `0` (never).
#'
#' @return A list with members:
#' \tabular{ll}{
#'   \code{$optVal}  \tab The optimal value of the objective function.\cr
#'   \code{$optArg}  \tab The argument vector which achieves that value.\cr
#'   \code{$ineqSat} \tab A vector of \code{0}s and \code{1}s. Each element
#'     corresponds to a single inequality, and is \code{1} if that inequality is
#'     satisfied for the optimal parameters, and \code{0} otherwise. Only
#'     present if \code{getIneqSat} is \code{TRUE}. \cr
#'   \code{$numSat} \tab The total number of satisfied inequalities. \cr
#'   \code{$bestRuns} \tab A list containing the results of optimization method
#'     that have yielded the best overall value. Each element is a list with
#'     members \code{$optVal}, \code{$optArg}, and \code{$ineqSat}. The first
#'     element of this list is returned as the default solution.
#' }
#' @export
optimizeScoreFunction <- function(
        dataArray, bounds,
        coefficient1 = NULL, method = NULL, optimParams = NULL,
        getIneqSat = FALSE, permuteInvariant = TRUE, numRuns = 1,
        progressUpdate = 0) {
    if (is.null(method)) {
        method <- "DEoptim"
    }
    stopifnot(numRuns >= 1)
    objDataArray <- dataArray
    if (permuteInvariant) {
        # TODO docs
        stdDevs <- apply(dataArray, 1, stats::sd)
        perm <- order(stdDevs[-1])
        # Now sort(stdDevs[-1]) == stdDevs[-1][perm].
        objDataArray[-1, ] <- objDataArray[-1, ][perm, ]
        # Now apply(objDataArray, 1, stats::sd)[-1] is sorted.
        invPerm <- rep(0, length(perm))
        for (i in 1:length(perm)) {
            invPerm[perm[i]] <- i
        }
    }
    switch (method,
        "DEoptim" = {
            objSign <- -1
            makeScoreObjFunArgs = list(dataArray = objDataArray, objSign = objSign)
            if (!is.null(coefficient1)) {
                makeScoreObjFunArgs$coefficient1 <- coefficient1
            }
            scoreObjFun <- do.call(makeScoreObjFun, makeScoreObjFunArgs)
            bestObjVal <- -Inf
            bestResult <- NULL
            runResults <- lapply(1:numRuns, function(runIdx) {
                runResult <- maximizeDEoptim(
                    scoreObjFun, bounds$lower, bounds$upper, optimParams)
                if (runResult$optVal > bestObjVal) {
                    bestResult <<- runResult
                    bestObjVal <<- runResult$optVal
                }
                if (progressUpdate > 0 && runIdx %% progressUpdate == 0) {
                    cat(sprintf("[optimizeScoreFunction] Iterations completed: %d\n", runIdx))
                }
                return(runResult)
            })
            bestIdxs <- lapply(
                runResults, function(runResult) { runResult$optVal }) == bestObjVal
            bestRuns <- runResults[bestIdxs]
            # Arbitrarily pick the first one as representative.
            result <- bestRuns[[1]]
        },
        stop(sprintf("optimizeScoreFunction: method %s is not implemented",
                      method))
    )
    if (permuteInvariant) {
        # Return parameters with the original order.
        result$optArg <- result$optArg[invPerm]
        makeScoreObjFunArgs$dataArray <- dataArray
        bestRuns <- lapply(bestRuns, function(runResult) {
            runResult$optArg <- runResult$optArg[invPerm]
            runResult
        })
    }
    if (getIneqSat) {
        scoreObjFunVec <- do.call(makeScoreObjFunVec, makeScoreObjFunArgs)
        result$ineqSat <- objSign * scoreObjFunVec(result$optArg)
        bestRuns <- lapply(bestRuns, function(runResult) {
            runResult$ineqSat <- objSign * scoreObjFunVec(runResult$optArg)
            runResult
        })
    }
    result$optArg <- unname(result$optArg)
    bestRuns <- lapply(bestRuns, function(runResult) {
        runResult$optArg <- unname(runResult$optArg)
        runResult
    })
    result$bestRuns <- bestRuns
    result$numSat <- result$optVal * dim(dataArray)[2]
    return(result)
}

#' Calculate bootstrap estimates
#'
#' Optimizes the objective function created be \code{makeBootstrapObjFun}. Its
#' argmax is used internally in the implementation of Cattaneo's bootstrap
#' method.
#'
#' @section Optimization methods:
#' The optimization method is not bound to the problem. Any method that can
#' optimize a non-convex, non-smooth function is valid. However, we currently
#' only support the *Differential Evolution* method, as implemented in the
#' package **DEoptim**. In case the user wants to pass parameters to the solver,
#' we provide the parameter \code{optimParams}, which is forwarded to the
#' \code{control} parameter of the function \code{DEoptim::DEoptim}.
#' See \code{DEoptim::DEoptim.control} for more information.
#'
#' @inheritParams optimizeScoreFunction
#' @inheritParams makeBootstrapObjFun
#'
#' @return A list with members:
#' \tabular{ll}{
#'   \code{$optVal}  \tab The optimal value of the objective function.\cr
#'   \code{$optArg}  \tab The argument vector which achieves that value.\cr
#' }
#' @export
optimizeBootstrapFunction <- function(
        fullDataArray, sampleDataArray, betaEst, H,
        bounds,
        coefficient1 = NULL, method = NULL, optimParams = NULL,
        useCorrectionFactor = TRUE) {
    if (is.null(method)) {
        method <- "DEoptim"
    }
    switch (method,
            "DEoptim" = {
                objSign <- -1
                makeBootstrapObjFunArgs = list(
                    fullDataArray = fullDataArray,
                    sampleDataArray = sampleDataArray,
                    betaEst = betaEst, H = H, objSign = objSign,
                    useCorrectionFactor = useCorrectionFactor)
                if (!is.null(coefficient1)) {
                    makeBootstrapObjFunArgs$coefficient1 <- coefficient1
                }
                # This function also returns some extra information relevant to
                # the calculation of the bootstrap function, useful for
                # experimenting.
                bootstrapExtraFunction <- do.call(
                    makeBootstrapExtraObjFun, makeBootstrapObjFunArgs)
                result <- maximizeDEoptim(
                    function(b) { bootstrapExtraFunction(b)$val },
                    bounds$lower, bounds$upper, optimParams)
                bootstrapEvalInfo <- bootstrapExtraFunction(result$optArg)
                result$bootstrapEvalInfo <- bootstrapEvalInfo
            },
            stop(sprintf("optimizeBootstrapFunction: method %s is not implemented",
                         method))
    )
    result$optArg <- unname(result$optArg)
    return(result)
}

#' Helper function for Differential Evolution
#'
#' An adapter function over \code{DEoptim::DEoptim}. When other methods are
#' implemented, they should follow this signature (for the three required
#' arguments).
#'
#' @param objective The objective function, produced by \code{makeScoreObjFun}.
#' @param lower,upper The lower and upper bounds in the box constraints.
#' @param control A list passed to \code{DEoptim::DEoptim}. See also
#'   \code{DEoptim::DEoptim.control}.
#' @return A list with elements:
#' \tabular{ll}{
#'   \code{$optVal} \tab The optimal value of the objective function. \cr
#'   \code{$optArg} \tab The argument for which the objective function attains
#'     that value, as a vector.
#' }
#' @keywords internal
maximizeDEoptim <- function(objective, lower, upper, control = list()) {
    outDEoptim <- DEoptim::DEoptim(objective, lower, upper, control = control)
    optArg <-  outDEoptim$optim$bestmem
    optVal <- -outDEoptim$optim$bestval
    return(list(optVal = optVal, optArg = optArg))
}

#' Get default values for the \code{optimParams} argument of \code{optimizeScoreFunction}
#'
#' @inheritParams optimizeScoreFunction
#' @return A list of parameters.
#' @export
getDefaultOptimParams <- function(method = NULL) {
    if (is.null(method)) {
        method = "DEoptim"
    }
    if (method == "DEoptim") {
        return(list(NP=50, F=0.6, CR=0.5, itermax=100, trace=FALSE, reltol=1e-3))
    }
    stop()
}

#' Create box constraints
#'
#' Creates a pair of vectors to be used as box constraints in
#' \code{optimizeScoreFunction}.
#'
#' @param numAttrs The total number of attributes in the data.
#' @param hi The upper bound for all attributes.
#' @param lo (optional) The lower bound for all attributes. If not given,
#'   defaults to \code{-hi}.
#' @return A list with elements \code{$upper}, \code{$lower}.
#' @export
#' @examples
#' makeBounds(5, 10)
#' makeBounds(5, 3, -4)
makeBounds <- function(numAttrs, hi, lo = NULL) {
    stopifnot(numAttrs >= 2, hi > 0, is.null(lo) || lo < 0)
    n <- numAttrs - 1
    if (is.null(lo)) {
        lo <- -hi
    }
    upper <- rep(hi, n)
    lower <- rep(lo, n)
    return(list(lower = lower, upper = upper))
}

#' Calculate per market statistics
#'
#' @param ineqSat The element \code{$ineqSat} from the return value of the
#'   function \code{optimizeScoreFunction}.
#' @param groupIDs A vector assigning each inequality to its market. See the
#'   function \code{makeGroupIDs}.
#' @return An array of dimension \code{(noM, 4)}. Its rows correspond to
#'   markets, and its columns are:
#'   \itemize{
#'     \item The market index.
#'     \item The total number of inequalities for that market.
#'     \item The number of inequalities satisfied by the calculated solution.
#'     \item The ratio of the two values above.
#'   }
#' @export
calcPerMarketStats <- function(ineqSat, groupIDs) {
    calcRow <- function (mIdx) {
        ineqIdxs <- which(groupIDs == mIdx)
        v <- ineqSat[ineqIdxs]
        numTotIneqs <- length(v)
        numSatIneqs <- sum(v)
        row <- c(mIdx, numTotIneqs, numSatIneqs, numSatIneqs/numTotIneqs)
        return(row)
    }
    mIdxs <- unique(groupIDs)
    numM <- length(mIdxs)
    perMarketStats <- t(sapply(mIdxs, calcRow))
    colnames(perMarketStats) <- c("Market ID", "Total inequalities",
                                  "Satisfied inequalities",
                                  "Satisfied/Total ratio")
    return(perMarketStats)
}
