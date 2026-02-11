# optimizeScoreFunction: optVal in [0,1], length(optArg) = noAttr-1, numSat consistent

test_that("optimizeScoreFunction returns valid structure and values", {
    filename <- system.file("extdata", "precomp_testdata.dat", package = "maxscoreest")
    skip_if_not(nzchar(filename))
    matchedData <- importMatched(filename)
    ineqmembers <- Cineqmembers(matchedData$mate)
    dataArray <- CdataArray(matchedData$distanceMatrices, ineqmembers)
    bounds <- makeBounds(matchedData$noAttr, 100)
    optimParams <- getDefaultOptimParams()
    # Use few iterations so test stays fast
    optimParams$itermax <- 20
    optimParams$NP <- 20
    set.seed(1)
    optResult <- optimizeScoreFunction(
        dataArray = dataArray,
        bounds = bounds,
        optimParams = optimParams,
        getIneqSat = TRUE,
        permuteInvariant = TRUE,
        numRuns = 1)
    expect_true(optResult$optVal >= 0 && optResult$optVal <= 1)
    expect_equal(length(optResult$optArg), matchedData$noAttr - 1)
    expect_equal(optResult$numSat, optResult$optVal * ncol(dataArray))
    expect_equal(length(optResult$ineqSat), ncol(dataArray))
    expect_true(all(optResult$ineqSat %in% c(0, 1)))
    expect_true("bestRuns" %in% names(optResult))
    expect_equal(length(optResult$bestRuns), 1)
})
