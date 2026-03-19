# tests/run_tests.R
# Run the full test suite.
#
# Usage:
#   Rscript tests/run_tests.R
#   or from an R session:
#   source("tests/run_tests.R")

if (!requireNamespace("testthat", quietly = TRUE))
  stop("Install testthat first: install.packages('testthat')")

library(testthat)
library(here)

results <- testthat::test_dir(
  here("tests", "testthat"),
  reporter = testthat::default_reporter()
)

if (any(as.data.frame(results)$failed > 0 | as.data.frame(results)$error > 0)) {
  message("\n[FAIL] Some tests failed. See output above.")
  quit(status = 1)
} else {
  message("\n[PASS] All tests passed.")
}
