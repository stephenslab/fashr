test_that("fdr_control computes cumulative FDR from sorted lfdr", {
  lfdr <- c(0.9, 0.01, 0.3, 0.02)
  names(lfdr) <- paste0("Dataset_", 1:4)
  obj <- structure(list(lfdr = lfdr), class = "fash")
  res <- suppressMessages(fdr_control(obj, alpha = 0.05, sort = TRUE))
  expect_equal(res$fdr_results$lfdr, sort(lfdr), ignore_attr = TRUE)
  expect_equal(res$fdr_results$FDR, cumsum(sort(lfdr)) / 1:4, ignore_attr = TRUE)
  # Only the first two pass FDR <= 0.05: (0.01, 0.015, ...)
  expect_setequal(res$significant_units, c("Dataset_2", "Dataset_4"))
  # sort = FALSE returns rows in the original dataset order
  res2 <- suppressMessages(fdr_control(obj, alpha = 0.05, sort = FALSE))
  expect_equal(res2$fdr_results$index, 1:4)
})

test_that("fdr_control requires lfdr in the fash object", {
  obj <- structure(list(), class = "fash")
  expect_error(fdr_control(obj), "lfdr")
})
