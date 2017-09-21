
library(genBaRcode)
context("raw data processing")

test_that("get.wobble.pos is number of N's within the regular expression", {
  expect_equal(.getWobblePos("xxxNNyyyNNNaaa"), c(4, 5, 9, 10, 11))
  expect_equal(.getWobblePos("xxxyyyNzNppp"), c(7, 9))
  expect_equal(.getWobblePos("xNx"), 2)
})
