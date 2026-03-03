library(here)
library(dplyr)

check_read <- function(path) {
  message("Checking: ", path)
  tryCatch({
    df <- read.csv(path)
    message("  Success. Rows: ", nrow(df), " Cols: ", ncol(df))
  }, error = function(e) {
    message("  FAILED: ", e$message)
  })
}

check_read(here("Study 3/Cleaning/output/fullTest_fixed.csv"))
check_read(here("Study 3/Cleaning/output/fullTrain_fixed.csv"))
check_read(here("Combined/input/adjacencyMatrix_p.csv"))
