# Tests for modifydata.R — R port of MSE-Mathematica/modifydata.m.

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

# Compare two mate data.tables without depending on data.table internals
# (.internal.selfref, key ordering, list-column NULL vs integer(0)).
expectMateEqual <- function(actual, expected) {
    testthat::expect_equal(
        as.integer(actual$UpStream),
        as.integer(expected$UpStream))
    normalize <- function(col) lapply(col, function(v) {
        as.integer(if (is.null(v)) integer(0) else v)
    })
    testthat::expect_equal(
        normalize(actual$DownMates),
        normalize(expected$DownMates))
}

# Small synthetic state with hand-traceable values.
#
# Market 1 (noU = 3, noD = 4): 1-to-1 diagonal matches
#     u=1 ↔ d=1 ; u=2 ↔ d=2 ; u=3 ↔ d=3 ; d=4 unmatched
#         quotasU = (1, 1, 1), quotasD = (1, 1, 1, 0)
#
# Market 2 (noU = 3, noD = 4): many-to-many
#     u=1 ↔ {d=1, d=2} ; u=2 ↔ d=3 ; u=3 unmatched ; d=4 unmatched
#         quotasU = (2, 1, 0), quotasD = (1, 1, 1, 0)
makeSyntheticState <- function(withQuotas = TRUE, withPayoffs = FALSE,
                               beta = c(0.5)) {
    dist1 <- array(as.numeric(seq_len(2 * 4 * 3)), dim = c(2, 4, 3))
    dist2 <- array(100 + as.numeric(seq_len(2 * 4 * 3)), dim = c(2, 4, 3))
    mm1 <- matrix(0L, 4, 3)
    mm1[1, 1] <- 1L; mm1[2, 2] <- 1L; mm1[3, 3] <- 1L
    mm2 <- matrix(0L, 4, 3)
    mm2[1, 1] <- 1L; mm2[2, 1] <- 1L; mm2[3, 2] <- 1L
    state <- list(
        distanceMatrices = list(dist1, dist2),
        matchMatrices    = list(mm1, mm2),
        mate             = list(mateFromMatchMatrix(mm1),
                                mateFromMatchMatrix(mm2)),
        noM    = 2L,
        noU    = c(3L, 3L),
        noD    = c(4L, 4L),
        noAttr = 2L)
    if (withQuotas) state <- addQuotasFromMatches(state)
    if (withPayoffs) {
        state$payoffMatrices <- evaluatePayoffMatrices(
            state$distanceMatrices, beta)
    }
    return(state)
}

# ----------------------------------------------------------------------------
# storeState / restoreState
# ----------------------------------------------------------------------------

test_that("storeState preserves all members", {
    state <- makeSyntheticState()
    snap <- storeState(state)
    expect_equal(snap$noU, state$noU)
    expect_equal(snap$noD, state$noD)
    expect_equal(snap$distanceMatrices, state$distanceMatrices)
    expect_equal(snap$matchMatrices,    state$matchMatrices)
})

test_that("storeState deep-copies the mate data.tables", {
    state <- makeSyntheticState()
    snap <- storeState(state)
    # Modify the snapshot's mate in place; the original must be unaffected.
    data.table::set(snap$mate[[1]], 1L, "DownMates", list(list(99L)))
    expect_false(identical(
        snap$mate[[1]]$DownMates[[1]],
        state$mate[[1]]$DownMates[[1]]))
    expect_equal(state$mate[[1]]$DownMates[[1]], 1L)
})

test_that("restoreState round-trips an unmodified snapshot", {
    state <- makeSyntheticState()
    recovered <- restoreState(storeState(state))
    expect_equal(recovered$noU,           state$noU)
    expect_equal(recovered$matchMatrices, state$matchMatrices)
    expectMateEqual(recovered$mate[[1]], state$mate[[1]])
})

# ----------------------------------------------------------------------------
# mateFromMatchMatrix
# ----------------------------------------------------------------------------

test_that("mateFromMatchMatrix on a diagonal match matrix", {
    mm <- matrix(0L, 3, 3); diag(mm) <- 1L
    expected <- data.table::data.table(
        UpStream  = 1:3,
        DownMates = list(1L, 2L, 3L))
    data.table::setkey(expected, UpStream)
    expectMateEqual(mateFromMatchMatrix(mm), expected)
})

test_that("mateFromMatchMatrix preserves many-to-many relationships", {
    mm <- matrix(0L, 4, 3)
    mm[1, 1] <- 1L; mm[2, 1] <- 1L; mm[3, 2] <- 1L  # u=3 unmatched
    expected <- data.table::data.table(
        UpStream  = 1:3,
        DownMates = list(c(1L, 2L), 3L, integer(0)))
    data.table::setkey(expected, UpStream)
    expectMateEqual(mateFromMatchMatrix(mm), expected)
})

test_that("mateFromMatchMatrix on all-zeros returns empty DownMates", {
    mm <- matrix(0L, 3, 4)
    result <- mateFromMatchMatrix(mm)
    expect_equal(nrow(result), 4)
    expect_equal(unname(sapply(result$DownMates, length)), c(0, 0, 0, 0))
})

test_that("mateFromMatchMatrix on a zero-upstream matrix is empty", {
    mm <- matrix(integer(0), nrow = 3, ncol = 0)
    result <- mateFromMatchMatrix(mm)
    expect_equal(nrow(result), 0)
})

test_that("mateFromMatchMatrix agrees with importMatched on bundled data", {
    filename <- system.file("extdata", "precomp_testdata.dat",
                            package = "maxscoreest")
    state <- importMatched(filename)
    for (m in seq_len(state$noM)) {
        expectMateEqual(
            mateFromMatchMatrix(state$matchMatrices[[m]]),
            state$mate[[m]])
    }
})

test_that("mateFromMatchMatrices applies pointwise", {
    state <- makeSyntheticState()
    rebuilt <- mateFromMatchMatrices(state$matchMatrices)
    expect_equal(length(rebuilt), state$noM)
    for (m in seq_len(state$noM)) {
        expectMateEqual(rebuilt[[m]], state$mate[[m]])
    }
})

# ----------------------------------------------------------------------------
# quotaFromMatchMatrix / addQuotasFromMatches
# ----------------------------------------------------------------------------

test_that("quotaFromMatchMatrix returns row and column sums", {
    mm <- matrix(0L, 4, 3)
    mm[1, 1] <- 1L; mm[2, 1] <- 1L; mm[3, 2] <- 1L
    q <- quotaFromMatchMatrix(mm)
    expect_equal(q$quotaU, c(2L, 1L, 0L))
    expect_equal(q$quotaD, c(1L, 1L, 1L, 0L))
    expect_equal(sum(q$quotaU), sum(q$quotaD))
})

test_that("addQuotasFromMatches fills $quotasU and $quotasD for all markets", {
    state <- makeSyntheticState(withQuotas = FALSE)
    expect_null(state$quotasU)
    expect_null(state$quotasD)
    state <- addQuotasFromMatches(state)
    expect_equal(length(state$quotasU), state$noM)
    expect_equal(length(state$quotasD), state$noM)
    expect_equal(state$quotasU[[1]], c(1L, 1L, 1L))
    expect_equal(state$quotasD[[1]], c(1L, 1L, 1L, 0L))
    expect_equal(state$quotasU[[2]], c(2L, 1L, 0L))
    expect_equal(state$quotasD[[2]], c(1L, 1L, 1L, 0L))
})

# ----------------------------------------------------------------------------
# modifyMarket — unmatch
# ----------------------------------------------------------------------------

test_that("modifyMarket(unmatch=TRUE) clears exactly one (d, u) cell", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 1, u = 2L, d = 2L,
                      unmatch = TRUE, remove = FALSE)
    expect_equal(s$matchMatrices[[1]][2, 2], 0)
    # Original cell untouched
    expect_equal(state$matchMatrices[[1]][2, 2], 1)
    # Other cells in the same row/col untouched
    expect_equal(s$matchMatrices[[1]][1, 1], 1)
    expect_equal(s$matchMatrices[[1]][3, 3], 1)
})

test_that("modifyMarket(unmatch=TRUE) rebuilds the mate consistently", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 1, u = 2L, d = 2L,
                      unmatch = TRUE, remove = FALSE)
    expectMateEqual(s$mate[[1]], mateFromMatchMatrix(s$matchMatrices[[1]]))
})

test_that("modifyMarket(unmatch=TRUE) errors on length mismatch", {
    state <- makeSyntheticState()
    expect_error(
        modifyMarket(state, m = 1, u = c(1L, 2L), d = 1L,
                     unmatch = TRUE, remove = FALSE),
        "same length")
})

test_that("modifyMarket(unmatch=TRUE) warns on non-matched pair", {
    state <- makeSyntheticState()
    # In market 1, (u=1, d=2) is not a match.
    expect_warning(
        modifyMarket(state, m = 1, u = 1L, d = 2L,
                     unmatch = TRUE, remove = FALSE),
        "not matched")
})

test_that("modifyMarket(unmatch=TRUE) updates only the requested quota side", {
    state <- makeSyntheticState()
    # Market 2: unmatch (u=1, d=1) — u=1 still has d=2, d=1 has no mates left.
    s <- modifyMarket(state, m = 2, u = 1L, d = 1L,
                      unmatch = TRUE, remove = FALSE,
                      quotaUpdateUpstream = TRUE,
                      quotaUpdateDownstream = FALSE)
    expect_equal(s$quotasU[[2]], c(1L, 1L, 0L))
    expect_equal(s$quotasD[[2]], state$quotasD[[2]])   # not updated
})

test_that("modifyMarket(unmatch=TRUE) updates both quotas when requested", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 2, u = 1L, d = 1L,
                      unmatch = TRUE, remove = FALSE,
                      quotaUpdateUpstream = TRUE,
                      quotaUpdateDownstream = TRUE)
    expect_equal(s$quotasU[[2]], c(1L, 1L, 0L))
    expect_equal(s$quotasD[[2]], c(0L, 1L, 1L, 0L))
})

test_that("modifyMarket(unmatch=TRUE, quotaUpdate*) errors without quotas in state", {
    state <- makeSyntheticState(withQuotas = FALSE)
    expect_error(
        modifyMarket(state, m = 1, u = 1L, d = 1L,
                     unmatch = TRUE, remove = FALSE,
                     quotaUpdateUpstream = TRUE),
        "quotasU")
    expect_error(
        modifyMarket(state, m = 1, u = 1L, d = 1L,
                     unmatch = TRUE, remove = FALSE,
                     quotaUpdateDownstream = TRUE),
        "quotasD")
})

# ----------------------------------------------------------------------------
# modifyMarket — remove with quotaReset = TRUE (preserve shape)
# ----------------------------------------------------------------------------

test_that("modifyMarket(remove, quotaReset=TRUE) preserves shape", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 2, u = c(1L, 3L),
                      remove = TRUE, quotaReset = TRUE)
    expect_equal(s$noU[2], state$noU[2])
    expect_equal(dim(s$matchMatrices[[2]]),    dim(state$matchMatrices[[2]]))
    expect_equal(dim(s$distanceMatrices[[2]]), dim(state$distanceMatrices[[2]]))
})

test_that("modifyMarket(remove, quotaReset=TRUE) zeros listed columns", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 2, u = c(1L, 3L),
                      remove = TRUE, quotaReset = TRUE)
    expect_true(all(s$matchMatrices[[2]][, c(1, 3)] == 0))
    expect_equal(s$quotasU[[2]][c(1, 3)], c(0L, 0L))
    # The column of u=2 is unaffected
    expect_equal(s$matchMatrices[[2]][, 2], state$matchMatrices[[2]][, 2])
})

test_that("modifyMarket(remove, quotaReset=TRUE) zeros listed rows for downstreams", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 2, d = 1L,
                      remove = TRUE, quotaReset = TRUE)
    expect_true(all(s$matchMatrices[[2]][1, ] == 0))
    expect_equal(s$quotasD[[2]][1], 0L)
})

test_that("modifyMarket(remove, quotaReset=TRUE, quotaUpdate) decrements opposite", {
    state <- makeSyntheticState()
    # Market 2, remove u=1 (mated with d=1, d=2): decrement quotasD[1] and [2].
    s <- modifyMarket(state, m = 2, u = 1L,
                      remove = TRUE, quotaReset = TRUE, quotaUpdate = TRUE)
    expect_equal(s$quotasD[[2]],
                 state$quotasD[[2]] - c(1L, 1L, 0L, 0L))
})

# ----------------------------------------------------------------------------
# modifyMarket — remove with quotaReset = FALSE (slice)
# ----------------------------------------------------------------------------

test_that("modifyMarket(remove, slice) shrinks upstream arrays", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 2, u = c(1L, 3L),
                      remove = TRUE, quotaReset = FALSE)
    expect_equal(s$noU[2], state$noU[2] - 2)
    expect_equal(dim(s$matchMatrices[[2]])[2], 1)
    expect_equal(dim(s$distanceMatrices[[2]])[3], 1)
    expect_equal(length(s$quotasU[[2]]), 1)
    # The surviving column is the original u=2's column.
    expect_equal(s$matchMatrices[[2]][, 1], state$matchMatrices[[2]][, 2])
    expect_equal(s$quotasU[[2]], state$quotasU[[2]][2])
})

test_that("modifyMarket(remove, slice) shrinks downstream arrays", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 2, d = c(2L, 4L),
                      remove = TRUE, quotaReset = FALSE)
    expect_equal(s$noD[2], state$noD[2] - 2)
    expect_equal(dim(s$matchMatrices[[2]])[1], 2)
    expect_equal(dim(s$distanceMatrices[[2]])[2], 2)
    expect_equal(length(s$quotasD[[2]]), 2)
})

test_that("modifyMarket(remove, slice) operates on payoffMatrices when present", {
    state <- makeSyntheticState(withPayoffs = TRUE)
    origDim <- dim(state$payoffMatrices[[1]])
    s <- modifyMarket(state, m = 1, u = 1L,
                      remove = TRUE, quotaReset = FALSE)
    expect_equal(dim(s$payoffMatrices[[1]]),
                 c(origDim[1], origDim[2] - 1))
})

test_that("modifyMarket(remove, slice, quotaUpdate) decrements opposite quotas", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 2, u = 1L,
                      remove = TRUE, quotaReset = FALSE, quotaUpdate = TRUE)
    expect_equal(s$quotasD[[2]],
                 state$quotasD[[2]] - c(1L, 1L, 0L, 0L))
    # u=1 was sliced out → quotasU shrinks to old positions 2, 3.
    expect_equal(s$quotasU[[2]], state$quotasU[[2]][c(2, 3)])
})

test_that("modifyMarket(remove, slice) handles u and d together", {
    # Market 2 trace:
    #   step 1: quotaUpdate before u removal:
    #     quotasD -= rowSums(mm2[, u=1]) = (1,1,0,0)
    #              = (1,1,1,0) - (1,1,0,0) = (0,0,1,0)
    #   step 2: slice u=1 → matchMatrix cols become (old u=2, old u=3)
    #     quotasU becomes (1, 0)
    #   step 3: quotaUpdate before d removal:
    #     in the sliced matrix, d=3 was mated with (old u=2) → colSums = (1,0)
    #     quotasU -= (1, 0) → (0, 0)
    #   step 4: slice d=3 → quotasD becomes (0, 0, 0)
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 2, u = 1L, d = 3L,
                      remove = TRUE, quotaReset = FALSE, quotaUpdate = TRUE)
    expect_equal(s$noU[2], 2)
    expect_equal(s$noD[2], 3)
    expect_equal(s$quotasU[[2]], c(0L, 0L))
    expect_equal(s$quotasD[[2]], c(0L, 0L, 0L))
})

test_that("modifyMarket(remove, slice) handles removal of all upstreams", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 1, u = c(1L, 2L, 3L),
                      remove = TRUE, quotaReset = FALSE)
    expect_equal(s$noU[1], 0)
    expect_equal(dim(s$matchMatrices[[1]])[2], 0)
    expect_equal(dim(s$distanceMatrices[[1]])[3], 0)
    expect_equal(length(s$quotasU[[1]]), 0)
    expect_equal(nrow(s$mate[[1]]), 0)
})

# ----------------------------------------------------------------------------
# modifyMarket — rematch
# ----------------------------------------------------------------------------

test_that("modifyMarket(rematch=TRUE) re-solves and respects quotas", {
    state <- makeSyntheticState(withPayoffs = TRUE)
    s <- modifyMarket(state, m = 1, remove = FALSE, rematch = TRUE)
    # Other markets unchanged.
    expect_equal(s$matchMatrices[[2]], state$matchMatrices[[2]])
    # Market 1's match matrix respects the current quotas.
    expect_true(all(colSums(s$matchMatrices[[1]]) <= s$quotasU[[1]]))
    expect_true(all(rowSums(s$matchMatrices[[1]]) <= s$quotasD[[1]]))
})

test_that("modifyMarket(rematch=TRUE) errors when payoffMatrices is missing", {
    state <- makeSyntheticState(withPayoffs = FALSE)
    expect_error(
        modifyMarket(state, m = 1, remove = FALSE, rematch = TRUE),
        "payoffMatrices")
})

test_that("modifyMarket(rematch=TRUE) errors when quotas are missing", {
    state <- makeSyntheticState(withQuotas = FALSE, withPayoffs = TRUE)
    expect_error(
        modifyMarket(state, m = 1, remove = FALSE, rematch = TRUE),
        "quotas")
})

# ----------------------------------------------------------------------------
# modifyMarket — no-op, validation, immutability
# ----------------------------------------------------------------------------

test_that("modifyMarket with empty u and d and no flags is a no-op", {
    state <- makeSyntheticState()
    s <- modifyMarket(state, m = 1)
    expect_equal(s$noU, state$noU)
    expect_equal(s$noD, state$noD)
    expect_equal(s$matchMatrices,    state$matchMatrices)
    expect_equal(s$distanceMatrices, state$distanceMatrices)
})

test_that("modifyMarket errors on invalid market index", {
    state <- makeSyntheticState()
    expect_error(modifyMarket(state, m = 0))
    expect_error(modifyMarket(state, m = 99))
    expect_error(modifyMarket(state, m = c(1, 2)))
})

test_that("modifyMarket does not mutate the caller's state", {
    state <- makeSyntheticState()
    origMM      <- state$matchMatrices[[1]]
    origQuotasU <- state$quotasU[[1]]
    origNoU     <- state$noU
    origMate    <- data.table::copy(state$mate[[1]])
    modifyMarket(state, m = 1, u = 2L, d = 2L,
                 unmatch = TRUE, remove = FALSE,
                 quotaUpdateUpstream = TRUE)
    expect_equal(state$matchMatrices[[1]], origMM)
    expect_equal(state$quotasU[[1]],       origQuotasU)
    expect_equal(state$noU,                origNoU)
    expectMateEqual(state$mate[[1]],       origMate)
})

# ----------------------------------------------------------------------------
# removeUpstreams / removeDownstreams / removeStreams
# ----------------------------------------------------------------------------

test_that("removeUpstreams slices everything and updates noU", {
    state <- makeSyntheticState()
    s <- removeUpstreams(state, m = 1, u = c(1L, 2L))
    expect_equal(s$noU[1], 1)
    expect_equal(dim(s$matchMatrices[[1]])[2], 1)
    expect_equal(dim(s$distanceMatrices[[1]])[3], 1)
    expect_equal(length(s$quotasU[[1]]), 1)
    expect_equal(nrow(s$mate[[1]]), 1)
})

test_that("removeUpstreams updateMatchMatrix=FALSE leaves matchMatrix and mate", {
    state <- makeSyntheticState()
    s <- removeUpstreams(state, m = 1, u = 1L, updateMatchMatrix = FALSE)
    # noU and distanceMatrices were still sliced...
    expect_equal(s$noU[1], 2)
    expect_equal(dim(s$distanceMatrices[[1]])[3], 2)
    # ...but matchMatrix and mate are untouched.
    expect_equal(s$matchMatrices[[1]], state$matchMatrices[[1]])
    expectMateEqual(s$mate[[1]], state$mate[[1]])
})

test_that("removeUpstreams with empty u is a no-op", {
    state <- makeSyntheticState()
    s <- removeUpstreams(state, m = 1, u = integer(0))
    expect_equal(s$noU,           state$noU)
    expect_equal(s$matchMatrices, state$matchMatrices)
})

test_that("removeUpstreams handles payoffMatrices if present", {
    state <- makeSyntheticState(withPayoffs = TRUE)
    s <- removeUpstreams(state, m = 1, u = 1L)
    expect_equal(dim(s$payoffMatrices[[1]])[2], state$noU[1] - 1)
})

test_that("removeDownstreams slices everything and updates noD", {
    state <- makeSyntheticState()
    s <- removeDownstreams(state, m = 1, d = c(2L, 4L))
    expect_equal(s$noD[1], 2)
    expect_equal(dim(s$matchMatrices[[1]])[1], 2)
    expect_equal(dim(s$distanceMatrices[[1]])[2], 2)
    expect_equal(length(s$quotasD[[1]]), 2)
})

test_that("removeDownstreams updateMatchMatrix=FALSE leaves matchMatrix and mate", {
    state <- makeSyntheticState()
    s <- removeDownstreams(state, m = 1, d = 4L, updateMatchMatrix = FALSE)
    expect_equal(s$noD[1], 3)
    expect_equal(s$matchMatrices[[1]], state$matchMatrices[[1]])
})

test_that("removeStreams composes removeUpstreams and removeDownstreams", {
    state <- makeSyntheticState()
    expected <- removeDownstreams(
                    removeUpstreams(state, m = 2, u = 3L),
                    m = 2, d = 4L)
    actual <- removeStreams(state, m = 2, u = 3L, d = 4L)
    expect_equal(actual$noU,           expected$noU)
    expect_equal(actual$noD,           expected$noD)
    expect_equal(actual$matchMatrices, expected$matchMatrices)
})

test_that("removeStreams with both empty is a no-op", {
    state <- makeSyntheticState()
    s <- removeStreams(state, m = 1)
    expect_equal(s$noU,           state$noU)
    expect_equal(s$noD,           state$noD)
    expect_equal(s$matchMatrices, state$matchMatrices)
})

# ----------------------------------------------------------------------------
# Integration on bundled data
# ----------------------------------------------------------------------------

test_that("integration: modify+rematch on bundled precomp_testdata", {
    filename <- system.file("extdata", "precomp_testdata.dat",
                            package = "maxscoreest")
    state <- importMatched(filename)
    state <- addQuotasFromMatches(state)
    state$payoffMatrices <- evaluatePayoffMatrices(
        state$distanceMatrices, c(0.5, 1.0, 1.5, 2.0))

    origNoU1 <- state$noU[1]
    origNoM  <- state$noM

    s <- modifyMarket(state, m = 1, u = 1L,
                      remove = TRUE, quotaReset = FALSE,
                      quotaUpdate = TRUE, rematch = TRUE)

    # Market 1 lost an upstream; other markets are untouched.
    expect_equal(s$noU[1], origNoU1 - 1)
    expect_equal(s$noM,    origNoM)
    expect_equal(s$matchMatrices[[2]], state$matchMatrices[[2]])
    expect_equal(s$matchMatrices[[3]], state$matchMatrices[[3]])

    # Rematched market 1 respects the (decremented) quotas.
    expect_true(all(colSums(s$matchMatrices[[1]]) <= s$quotasU[[1]]))
    expect_true(all(rowSums(s$matchMatrices[[1]]) <= s$quotasD[[1]]))

    # The mate table is consistent with the rebuilt match matrix.
    expectMateEqual(s$mate[[1]], mateFromMatchMatrix(s$matchMatrices[[1]]))
})
