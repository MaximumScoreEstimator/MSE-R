# Routines for manipulating markets. R port of MSE-Mathematica/modifydata.m.
#
# The Mathematica module mutates a fixed set of top-level globals (matchMatrix,
# payoffMatrix, quota, mate, noU, noD, ...). The R package passes data
# explicitly, so we operate on a "market state" list — typically the return
# value of `importMatched`, optionally augmented with `$quotasU`, `$quotasD`,
# and `$payoffMatrices` for quota-aware and rematch operations.

#' Snapshot the current market state
#'
#' Returns a deep copy of the given state, suitable for restoring later.
#' Useful before calling \code{\link{modifyMarket}} (or another destructive
#' helper): assign the snapshot to a variable, perform modifications, then
#' call \code{\link{restoreState}} to recover the original. This mirrors the
#' \code{store}/\code{restore} pair in the Mathematica library.
#'
#' R uses value semantics for most types, so reassigning a state list already
#' protects the original. The exception is the \code{$mate} member: it is a
#' list of \code{data.table} objects, which carry reference semantics. This
#' function calls \code{data.table::copy} on each mate table to break the
#' aliasing.
#'
#' @section Market state structure:
#' The market state is a list with at least the following members (matching
#' the return value of \code{\link{importMatched}}):
#' \tabular{ll}{
#'   \code{$distanceMatrices} \tab A list of distance arrays, one per market;
#'     element \code{m} has dimension \code{(noAttr, noD[m], noU[m])}. \cr
#'   \code{$matchMatrices}    \tab A list of \code{0/1} match matrices, one per
#'     market; element \code{m} has dimension \code{(noD[m], noU[m])}. \cr
#'   \code{$mate}             \tab A list of mate \code{data.table} objects,
#'     one per market. See \code{\link{extractMate}}. \cr
#'   \code{$noM}, \code{$noU}, \code{$noD}, \code{$noAttr}
#'                            \tab Shape metadata.
#' }
#' Quota operations and rematching additionally use:
#' \tabular{ll}{
#'   \code{$quotasU}, \code{$quotasD} \tab Per-market quota vectors of length
#'     \code{noU[m]} and \code{noD[m]} respectively. These can be derived from
#'     the current match matrices with \code{\link{addQuotasFromMatches}}. \cr
#'   \code{$payoffMatrices}           \tab A list of evaluated payoff matrices
#'     (one per market, dimension \code{(noD[m], noU[m])}), required only when
#'     rematching. See \code{\link{evaluatePayoffMatrices}}.
#' }
#'
#' @param state The current market state.
#' @return A deep copy of \code{state}.
#' @seealso \code{\link{restoreState}}, \code{\link{modifyMarket}}
#' @export
storeState <- function(state) {
    snapshot <- state
    if (!is.null(state$mate)) {
        snapshot$mate <- lapply(state$mate, data.table::copy)
    }
    return(snapshot)
}

#' Restore a previously stored market state
#'
#' Returns a deep copy of the snapshot. Functionally identical to
#' \code{\link{storeState}}, exposed under a separate name to mirror the
#' Mathematica API.
#'
#' @param snapshot The state previously returned by \code{\link{storeState}}.
#' @return A deep copy of \code{snapshot}.
#' @seealso \code{\link{storeState}}, \code{\link{modifyMarket}}
#' @export
restoreState <- function(snapshot) {
    return(storeState(snapshot))
}

#' Compute a mate table from a single match matrix
#'
#' Inverse of \code{\link{extractMatchMatrices}}: rebuilds the mate
#' \code{data.table} (with columns \code{UpStream} and \code{DownMates}, see
#' \code{\link{extractMate}}) from a single market's match matrix. This is the
#' R analogue of the Mathematica function \code{Cmate}.
#'
#' The output uses the row indices of the input as upstream identifiers, so if
#' the match matrix has previously been sliced, the returned upstream indices
#' are renumbered consecutively from \code{1}.
#'
#' @param matchMatrix A 2-D array of dimension \code{(noD, noU)}, with values
#'   in \code{\{0, 1\}}.
#' @return A \code{data.table} keyed on \code{UpStream}, with one row per
#'   upstream. The \code{DownMates} column is a list column whose elements are
#'   integer vectors of downstream indices (\code{integer(0)} for unmatched
#'   upstreams).
#' @import data.table
#' @export
mateFromMatchMatrix <- function(matchMatrix) {
    UpStream <- NULL
    dims <- dim(matchMatrix)
    numU <- if (is.null(dims) || length(dims) < 2) 0L else dims[2]
    if (numU == 0) {
        return(data.table::data.table(
            UpStream = integer(0),
            DownMates = list()))
    }
    downMates <- lapply(seq_len(numU), function(uIdx) {
        as.integer(which(matchMatrix[, uIdx] == 1))
    })
    dt <- data.table::data.table(
        UpStream = seq_len(numU),
        DownMates = downMates)
    data.table::setkey(dt, UpStream)
    return(dt)
}

#' Compute mate tables for all markets
#'
#' List version of \code{\link{mateFromMatchMatrix}}.
#'
#' @param matchMatrices A list of match matrices, one per market.
#' @return A list of mate tables, one per market.
#' @export
mateFromMatchMatrices <- function(matchMatrices) {
    return(lapply(matchMatrices, mateFromMatchMatrix))
}

#' Derive quotas from a match matrix
#'
#' R analogue of the Mathematica function \code{Cquota} applied to a single
#' market. Computes per-agent match counts. In the modify routines these
#' counts double as the maximum allowed match counts during a subsequent
#' rematch (see \code{\link{generateAssignmentMatrix}}).
#'
#' @param matchMatrix A 2-D array of dimension \code{(noD, noU)}.
#' @return A list with two integer vectors:
#' \tabular{ll}{
#'   \code{$quotaU} \tab Column sums of \code{matchMatrix} (length
#'     \code{noU}). \cr
#'   \code{$quotaD} \tab Row sums of \code{matchMatrix} (length \code{noD}).
#' }
#' @export
quotaFromMatchMatrix <- function(matchMatrix) {
    return(list(
        quotaU = as.integer(colSums(matchMatrix)),
        quotaD = as.integer(rowSums(matchMatrix))))
}

#' Add per-market quotas to a state, derived from match matrices
#'
#' Convenience helper for using \code{\link{modifyMarket}} with quota
#' operations: derives \code{$quotasU} and \code{$quotasD} from the existing
#' \code{$matchMatrices} and adds them to the state.
#'
#' @inheritSection storeState Market state structure
#'
#' @param state The current market state.
#' @return The same state with members \code{$quotasU} and \code{$quotasD}
#'   added (each a list of per-market integer vectors).
#' @export
addQuotasFromMatches <- function(state) {
    quotas <- lapply(state$matchMatrices, quotaFromMatchMatrix)
    state$quotasU <- lapply(quotas, function(q) q$quotaU)
    state$quotasD <- lapply(quotas, function(q) q$quotaD)
    return(state)
}

#' Modify a market by unmatching agents, removing them, and/or rematching
#'
#' Port of the Mathematica function \code{modify}. Edits a single market in
#' one of two mutually exclusive ways (\code{unmatch} or \code{remove}), then
#' optionally re-solves the assignment problem for that market.
#'
#' @section Operations:
#' Exactly one of \code{unmatch} or \code{remove} controls the main operation
#' (\code{unmatch} takes precedence). If both are \code{FALSE}, only
#' \code{rematch} (if requested) takes effect.
#'
#' \describe{
#'   \item{\code{unmatch}}{
#'     For each \code{i}, sets \code{matchMatrices[[m]][d[i], u[i]]} to
#'     \code{0}. A warning is issued for any pair that was not actually
#'     matched. Vectors \code{u} and \code{d} must have the same length.
#'
#'     If \code{quotaUpdateUpstream} is \code{TRUE}, the upstream quotas for
#'     market \code{m} are recomputed from the new match matrix (i.e. column
#'     sums). \code{quotaUpdateDownstream} does the same for downstream
#'     quotas (row sums).
#'   }
#'   \item{\code{remove} with \code{quotaReset = FALSE} (the default)}{
#'     Slices the rows of \code{matchMatrices[[m]]} corresponding to
#'     upstreams \code{u} (and the columns for downstreams \code{d}); also
#'     slices the corresponding entries of \code{distanceMatrices[[m]]},
#'     \code{payoffMatrices[[m]]} (if present), and the quota vectors.
#'     Decrements \code{noU[m]} and \code{noD[m]} accordingly.
#'
#'     Indices for remaining agents are renumbered consecutively from
#'     \code{1}.
#'   }
#'   \item{\code{remove} with \code{quotaReset = TRUE}}{
#'     Preserves shape: zeros out the rows/columns and quotas for the listed
#'     agents instead of slicing them out. \code{noU[m]} and \code{noD[m]} do
#'     not change.
#'   }
#'   \item{\code{rematch}}{
#'     Re-solves the assignment problem for market \code{m} using the current
#'     \code{payoffMatrices[[m]]} and quotas, via
#'     \code{\link{generateAssignmentMatrix}}. The state must contain
#'     \code{$payoffMatrices}, \code{$quotasU}, and \code{$quotasD}.
#'   }
#' }
#'
#' When \code{remove = TRUE}, the flag \code{quotaUpdate} additionally
#' decrements the opposite stream's quotas to reflect each match that the
#' removed agents participated in. For example, removing an upstream that was
#' matched to two downstreams decrements those two downstream quotas by 1
#' each; if both upstreams \code{u_1} and \code{u_2} were matched to
#' downstream \code{d}, then \code{quotasD[[m]][d]} is decremented twice.
#'
#' The mate table for market \code{m} is rebuilt after each modification.
#'
#' @inheritSection storeState Market state structure
#'
#' @param state The current market state.
#' @param m The market index.
#' @param u,d Integer vectors of upstream and downstream indices to operate on.
#'   When \code{unmatch = TRUE}, the two must have the same length and
#'   correspond pairwise: \code{(u[i], d[i])} is the \code{i}-th match to
#'   unmatch. When \code{remove = TRUE}, they list the agents to remove (or
#'   reset). Both default to empty.
#' @param unmatch,remove,quotaReset,quotaUpdate,quotaUpdateUpstream,quotaUpdateDownstream,rematch
#'   See "Operations". All default to \code{FALSE}, except \code{remove},
#'   which defaults to \code{TRUE}.
#' @return The modified state. The original is not mutated.
#' @seealso \code{\link{removeUpstreams}}, \code{\link{removeDownstreams}},
#'   \code{\link{removeStreams}}, \code{\link{addQuotasFromMatches}},
#'   \code{\link{storeState}}
#' @export
modifyMarket <- function(state, m, u = integer(0), d = integer(0),
                         unmatch = FALSE, remove = TRUE,
                         quotaReset = FALSE, quotaUpdate = FALSE,
                         quotaUpdateUpstream = FALSE,
                         quotaUpdateDownstream = FALSE,
                         rematch = FALSE) {
    stopifnot(
        is.numeric(m), length(m) == 1, m >= 1, m <= state$noM,
        is.numeric(u) || length(u) == 0,
        is.numeric(d) || length(d) == 0)
    u <- as.integer(u)
    d <- as.integer(d)
    state <- storeState(state)

    if (unmatch) {
        if (length(u) != length(d)) {
            stop(sprintf(
                "modifyMarket: 'u' and 'd' must have the same length when unmatch=TRUE (got %d and %d)",
                length(u), length(d)))
        }
        matchMat <- state$matchMatrices[[m]]
        for (i in seq_along(u)) {
            if (matchMat[d[i], u[i]] != 1) {
                warning(sprintf(
                    "modifyMarket: in market %d, upstream %d is not matched with downstream %d",
                    m, u[i], d[i]))
            }
            matchMat[d[i], u[i]] <- 0L
        }
        state$matchMatrices[[m]] <- matchMat
        state$mate[[m]] <- mateFromMatchMatrix(matchMat)
        if (quotaUpdateUpstream) {
            if (is.null(state$quotasU)) stop(
                "modifyMarket: quotaUpdateUpstream=TRUE requires state$quotasU. ",
                "Call addQuotasFromMatches(state) first.")
            state$quotasU[[m]] <- as.integer(colSums(matchMat))
        }
        if (quotaUpdateDownstream) {
            if (is.null(state$quotasD)) stop(
                "modifyMarket: quotaUpdateDownstream=TRUE requires state$quotasD. ",
                "Call addQuotasFromMatches(state) first.")
            state$quotasD[[m]] <- as.integer(rowSums(matchMat))
        }
    } else if (remove) {
        # Process upstreams first; the resulting matrices feed into the
        # downstream step (matching the Mathematica order).
        if (length(u) > 0) {
            if (quotaUpdate) {
                if (is.null(state$quotasD)) stop(
                    "modifyMarket: quotaUpdate=TRUE requires state$quotasD. ",
                    "Call addQuotasFromMatches(state) first.")
                removedPerD <- rowSums(state$matchMatrices[[m]][, u, drop = FALSE])
                state$quotasD[[m]] <- state$quotasD[[m]] - as.integer(removedPerD)
            }
            if (quotaReset) {
                if (is.null(state$quotasU)) stop(
                    "modifyMarket: quotaReset=TRUE requires state$quotasU. ",
                    "Call addQuotasFromMatches(state) first.")
                state$quotasU[[m]][u] <- 0L
                state$matchMatrices[[m]][, u] <- 0L
                state$mate[[m]] <- mateFromMatchMatrix(state$matchMatrices[[m]])
            } else {
                state$distanceMatrices[[m]] <-
                    state$distanceMatrices[[m]][, , -u, drop = FALSE]
                if (!is.null(state$payoffMatrices)) {
                    state$payoffMatrices[[m]] <-
                        state$payoffMatrices[[m]][, -u, drop = FALSE]
                }
                state$matchMatrices[[m]] <-
                    state$matchMatrices[[m]][, -u, drop = FALSE]
                state$mate[[m]] <- mateFromMatchMatrix(state$matchMatrices[[m]])
                if (!is.null(state$quotasU)) {
                    state$quotasU[[m]] <- state$quotasU[[m]][-u]
                }
                state$noU[m] <- state$noU[m] - length(u)
            }
        }
        if (length(d) > 0) {
            if (quotaUpdate) {
                if (is.null(state$quotasU)) stop(
                    "modifyMarket: quotaUpdate=TRUE requires state$quotasU. ",
                    "Call addQuotasFromMatches(state) first.")
                removedPerU <- colSums(state$matchMatrices[[m]][d, , drop = FALSE])
                state$quotasU[[m]] <- state$quotasU[[m]] - as.integer(removedPerU)
            }
            if (quotaReset) {
                if (is.null(state$quotasD)) stop(
                    "modifyMarket: quotaReset=TRUE requires state$quotasD. ",
                    "Call addQuotasFromMatches(state) first.")
                state$quotasD[[m]][d] <- 0L
                state$matchMatrices[[m]][d, ] <- 0L
                state$mate[[m]] <- mateFromMatchMatrix(state$matchMatrices[[m]])
            } else {
                state$distanceMatrices[[m]] <-
                    state$distanceMatrices[[m]][, -d, , drop = FALSE]
                if (!is.null(state$payoffMatrices)) {
                    state$payoffMatrices[[m]] <-
                        state$payoffMatrices[[m]][-d, , drop = FALSE]
                }
                state$matchMatrices[[m]] <-
                    state$matchMatrices[[m]][-d, , drop = FALSE]
                state$mate[[m]] <- mateFromMatchMatrix(state$matchMatrices[[m]])
                if (!is.null(state$quotasD)) {
                    state$quotasD[[m]] <- state$quotasD[[m]][-d]
                }
                state$noD[m] <- state$noD[m] - length(d)
            }
        }
    }

    if (rematch) {
        if (is.null(state$payoffMatrices)) stop(
            "modifyMarket: rematch=TRUE requires state$payoffMatrices. ",
            "Call evaluatePayoffMatrices(state$distanceMatrices, beta) and ",
            "store the result as state$payoffMatrices first.")
        if (is.null(state$quotasU) || is.null(state$quotasD)) stop(
            "modifyMarket: rematch=TRUE requires state$quotasU and state$quotasD. ",
            "Call addQuotasFromMatches(state) first.")
        state$matchMatrices[[m]] <- generateAssignmentMatrix(
            state$payoffMatrices[[m]],
            state$quotasU[[m]],
            state$quotasD[[m]])
        state$mate[[m]] <- mateFromMatchMatrix(state$matchMatrices[[m]])
    }

    return(state)
}

#' Remove upstreams from a market
#'
#' Lightweight helper for the common case of slicing upstreams out of a single
#' market. R port of the Mathematica function \code{removeU}.
#'
#' Equivalent to \code{modifyMarket(state, m, u = u, remove = TRUE,
#' quotaReset = FALSE)}, but with an additional flag to leave
#' \code{matchMatrices[[m]]} untouched — useful in the precomputed-distances
#' workflow where the match matrix is not actively maintained.
#'
#' @inheritSection storeState Market state structure
#' @inheritParams modifyMarket
#' @param updateMatchMatrix If \code{TRUE} (default), also removes the
#'   corresponding columns from \code{state$matchMatrices[[m]]} and rebuilds
#'   \code{state$mate[[m]]}. If \code{FALSE}, the match matrix and mate are
#'   left intact (their shape will then no longer agree with the other
#'   arrays).
#' @return The modified state.
#' @seealso \code{\link{removeDownstreams}}, \code{\link{removeStreams}},
#'   \code{\link{modifyMarket}}
#' @export
removeUpstreams <- function(state, m, u, updateMatchMatrix = TRUE) {
    stopifnot(
        is.numeric(m), length(m) == 1, m >= 1, m <= state$noM,
        is.numeric(u) || length(u) == 0)
    u <- as.integer(u)
    if (length(u) == 0) return(state)
    state <- storeState(state)
    state$distanceMatrices[[m]] <- state$distanceMatrices[[m]][, , -u, drop = FALSE]
    if (!is.null(state$payoffMatrices)) {
        state$payoffMatrices[[m]] <- state$payoffMatrices[[m]][, -u, drop = FALSE]
    }
    if (!is.null(state$quotasU)) {
        state$quotasU[[m]] <- state$quotasU[[m]][-u]
    }
    if (updateMatchMatrix) {
        state$matchMatrices[[m]] <- state$matchMatrices[[m]][, -u, drop = FALSE]
        state$mate[[m]] <- mateFromMatchMatrix(state$matchMatrices[[m]])
    }
    state$noU[m] <- state$noU[m] - length(u)
    return(state)
}

#' Remove downstreams from a market
#'
#' R port of the Mathematica function \code{removeD}. See
#' \code{\link{removeUpstreams}} for the upstream variant.
#'
#' @inheritSection storeState Market state structure
#' @inheritParams modifyMarket
#' @param updateMatchMatrix See \code{\link{removeUpstreams}}.
#' @return The modified state.
#' @seealso \code{\link{removeUpstreams}}, \code{\link{removeStreams}},
#'   \code{\link{modifyMarket}}
#' @export
removeDownstreams <- function(state, m, d, updateMatchMatrix = TRUE) {
    stopifnot(
        is.numeric(m), length(m) == 1, m >= 1, m <= state$noM,
        is.numeric(d) || length(d) == 0)
    d <- as.integer(d)
    if (length(d) == 0) return(state)
    state <- storeState(state)
    state$distanceMatrices[[m]] <- state$distanceMatrices[[m]][, -d, , drop = FALSE]
    if (!is.null(state$payoffMatrices)) {
        state$payoffMatrices[[m]] <- state$payoffMatrices[[m]][-d, , drop = FALSE]
    }
    if (!is.null(state$quotasD)) {
        state$quotasD[[m]] <- state$quotasD[[m]][-d]
    }
    if (updateMatchMatrix) {
        state$matchMatrices[[m]] <- state$matchMatrices[[m]][-d, , drop = FALSE]
        state$mate[[m]] <- mateFromMatchMatrix(state$matchMatrices[[m]])
    }
    state$noD[m] <- state$noD[m] - length(d)
    return(state)
}

#' Remove upstreams and downstreams from a market
#'
#' Convenience helper combining \code{\link{removeUpstreams}} and
#' \code{\link{removeDownstreams}}. R port of the Mathematica function
#' \code{removeUD}.
#'
#' @inheritSection storeState Market state structure
#' @inheritParams modifyMarket
#' @inheritParams removeUpstreams
#' @return The modified state.
#' @seealso \code{\link{removeUpstreams}}, \code{\link{removeDownstreams}},
#'   \code{\link{modifyMarket}}
#' @export
removeStreams <- function(state, m, u = integer(0), d = integer(0),
                          updateMatchMatrix = TRUE) {
    state <- removeUpstreams(state, m, u, updateMatchMatrix = updateMatchMatrix)
    state <- removeDownstreams(state, m, d, updateMatchMatrix = updateMatchMatrix)
    return(state)
}
