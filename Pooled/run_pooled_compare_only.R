# Compare all 4 pooled models from saved LOO files (no fitting).
library(loo)
library(here)

all_model_names <- c("bias", "symmetric", "sym_lambda", "asym_lambda")

loo_list <- setNames(
  lapply(all_model_names, function(m) readRDS(here("Fits", paste0("loo_pooled_", m, ".rds")))),
  all_model_names
)

diag_list <- lapply(all_model_names, function(m) {
  f <- here("Results", paste0("diag_pooled_", m, ".csv"))
  if (file.exists(f)) read.csv(f) else NULL
})
diag_list <- Filter(Negate(is.null), diag_list)

comp <- loo_compare(loo_list)
print(comp)

write.csv(as.data.frame(comp), here("Results", "pooled_model_comparison_results.csv"), row.names = FALSE)
if (length(diag_list) > 0)
  write.csv(do.call(rbind, diag_list), here("Results", "pooled_model_diagnostics.csv"), row.names = FALSE)

message("Pooled model comparison complete.")
