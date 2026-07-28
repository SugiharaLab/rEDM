#------------------------------------------------------------------------
# Phase 0 fixture harness.
#
# Replicates the pyEDM validation-driver comparison:
#     round(reference, digits).equals( round(result, digits) )
# with three fidelity details that cause spurious failures if omitted:
#   (1) NaN-as-equal    : same-position NA compare equal (pandas .equals)
#   (2) row slicing     : compare only a caller-specified row range
#   (3) name reconcile  : e.g. pyEDM 'theta' vs rEDM 'Theta'
#
# Base R only (no RANN, no testthat needed by the core comparator), so the
# Phase 0 acceptance test runs before any kernel exists.
#
# Naming follows the rEDM convention: functions and constants begin with an
# uppercase letter; variables begin with a lowercase letter.
#------------------------------------------------------------------------

# Locate the fixtures directory whether run under testthat or via source().
FixtureDir <- function() {
  cand <- c(
    if (requireNamespace("testthat", quietly = TRUE))
      tryCatch(testthat::test_path("fixtures"), error = function(e) NULL),
    file.path("tests", "testthat", "fixtures"),
    file.path("fixtures"),
    Sys.getenv("EDM_FIXTURE_DIR", unset = NA)
  )
  cand <- cand[!is.na(cand) & nzchar(cand)]
  hit <- cand[dir.exists(cand)]
  if (length(hit) == 0L)
    stop("FixtureDir: fixtures directory not found; tried: ",
         paste(cand, collapse = ", "))
  hit[1]
}

# Load a fixture CSV by base name (with or without the _valid suffix).
# Empty fields and nan/NaN/NA are read as NA; strings are not factors.
LoadFixture <- function(name, dir = FixtureDir()) {
  cand <- file.path(dir, c(paste0(name, ".csv"),
                           paste0(name, "_valid.csv")))
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L)
    stop("LoadFixture: no CSV for '", name, "' in ", dir)
  read.csv(hit[1], header = TRUE, stringsAsFactors = FALSE,
           check.names = FALSE,
           na.strings = c("NA", "nan", "NaN", "NAN", ""))
}

#------------------------------------------------------------------------
# Core comparator.
#
#   result, reference : data.frames
#   digits            : rounding digits for numeric columns
#   rows              : optional integer row indices to compare (both sides)
#   reconcile         : named char vector c(referenceName = resultName);
#                       renames the given result column(s) to the reference
#                       label before comparison (e.g. c(theta = "Theta"))
#
# Returns a structured list so later phases get actionable failure output.
#------------------------------------------------------------------------
CompareEDM <- function(result, reference, digits = 6,
                       rows = NULL, reconcile = NULL) {
  out <- list(equal = FALSE, digits = digits,
              dimResult = dim(result), dimReference = dim(reference),
              nameDiff = NULL, ncolMismatch = FALSE,
              nrowMismatch = FALSE, diffs = NULL)

  # (3) name reconciliation on the result side
  if (!is.null(reconcile)) {
    for (refName in names(reconcile)) {
      resName <- reconcile[[refName]]
      j <- match(resName, colnames(result))
      if (!is.na(j)) colnames(result)[j] <- refName
    }
  }

  # column-name / arity check
  if (!identical(colnames(result), colnames(reference))) {
    out$nameDiff <- list(result = colnames(result),
                         reference = colnames(reference))
    out$ncolMismatch <- ncol(result) != ncol(reference)
    return(out)
  }

  # (2) optional row slice, applied to both sides
  if (!is.null(rows)) {
    result    <- result[rows, , drop = FALSE]
    reference <- reference[rows, , drop = FALSE]
    out$dimResult    <- dim(result)
    out$dimReference <- dim(reference)
  }

  if (nrow(result) != nrow(reference)) {
    out$nrowMismatch <- TRUE
    return(out)
  }

  diffs <- NULL
  for (col in colnames(reference)) {
    a <- result[[col]]
    b <- reference[[col]]

    if (is.numeric(a) && is.numeric(b)) {
      roundA <- round(a, digits)
      roundB <- round(b, digits)
      naA <- is.na(roundA); naB <- is.na(roundB)
      # (1) NaN-as-equal: mismatch only if NA-masks differ, or both finite & unequal
      bad <- (naA != naB) | (!naA & !naB & roundA != roundB)
    } else {
      charA <- as.character(a); charB <- as.character(b)
      naA <- is.na(charA); naB <- is.na(charB)
      bad <- (naA != naB) | (!naA & !naB & charA != charB)
    }

    if (any(bad)) {
      w <- which(bad)
      diffs <- rbind(diffs, data.frame(
        row = w, column = col,
        result    = as.character(a[w]),
        reference = as.character(b[w]),
        stringsAsFactors = FALSE))
    }
  }

  out$diffs <- diffs
  out$equal <- is.null(diffs)
  out
}

# Thin testthat wrapper (only used when testthat is loaded).
ExpectEDMEqual <- function(result, reference, digits = 6,
                           rows = NULL, reconcile = NULL) {
  cmp <- CompareEDM(result, reference, digits, rows, reconcile)
  if (!cmp$equal && !is.null(cmp$diffs)) {
    msg <- sprintf("%d cell mismatch(es) at %d digits; first: row %s col '%s' result=%s ref=%s",
                   nrow(cmp$diffs), digits, cmp$diffs$row[1], cmp$diffs$column[1],
                   cmp$diffs$result[1], cmp$diffs$reference[1])
  } else if (!cmp$equal) {
    msg <- "structural mismatch (column names / dimensions)"
  } else {
    msg <- "equal"
  }
  testthat::expect(cmp$equal, msg)
  invisible(cmp)
}

#------------------------------------------------------------------------
# Sample-data loader (works installed via system.file or from source tree).
#------------------------------------------------------------------------
SampleDataDir <- function() {
  cand <- c(
    if (requireNamespace("rEDM", quietly = TRUE))
      system.file("extdata", package = "rEDM"),
    file.path("inst", "extdata"),
    file.path("..", "..", "inst", "extdata"),
    Sys.getenv("EDM_DATA_DIR", unset = NA)
  )
  cand <- cand[!is.na(cand) & nzchar(cand)]
  hit <- cand[dir.exists(cand)]
  if (length(hit) == 0L) stop("SampleDataDir: extdata not found")
  hit[1]
}

LoadSampleData <- function(name) {
  read.csv(file.path(SampleDataDir(), name), check.names = FALSE)
}
