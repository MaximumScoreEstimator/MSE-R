# Smoke tests for pointIdentifiedCR (and optionally cubeRootBootstrapCR)

test_that("pointIdentifiedCR runs and returns expected structure", {
    filename <- system.file("extdata", "precomp_testdata.dat", package = "maxscoreest")
    skip_if_not(nzchar(filename))
    matchedData <- importMatched(filename)
    ineqmembers <- Cineqmembers(matchedData$mate)
    dataArray <- CdataArray(matchedData$distanceMatrices, ineqmembers)
    bounds <- makeBounds(matchedData$noAttr, 100)
    optimParams <- getDefaultOptimParams()
    optimParams$itermax <- 20
    optimParams$NP <- 20
    set.seed(2)
    optResult <- optimizeScoreFunction(
        dataArray = dataArray,
        bounds = bounds,
        optimParams = optimParams,
        getIneqSat = FALSE,
        permuteInvariant = TRUE,
        numRuns = 1)
    groupIDs <- makeGroupIDs(ineqmembers)
    optimizeScoreArgs <- list(
        bounds = bounds,
        optimParams = optimParams,
        getIneqSat = FALSE,
        permuteInvariant = TRUE,
        numRuns = 1)
    set.seed(3)
    ssSize <- min(2, matchedData$noM)
    numSubsamples <- 2
    confidenceLevel <- 0.95
    cr <- pointIdentifiedCR(
        dataArray, groupIDs, optResult$optArg,
        ssSize, numSubsamples, confidenceLevel,
        optimizeScoreArgs,
        options = list(progressUpdate = 0))
    expect_true("crSymm" %in% names(cr))
    expect_true("crAsym" %in% names(cr))
    expect_true("estimates" %in% names(cr))
    expect_true("samples" %in% names(cr))
    expect_equal(nrow(cr$estimates), numSubsamples)
    expect_equal(ncol(cr$estimates), length(optResult$optArg))
    expect_equal(dim(cr$samples), c(numSubsamples, ssSize))
})
