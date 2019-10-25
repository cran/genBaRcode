library(genBaRcode)
context("helper functions")

test_that("helper functions - com_pair", {

  dat1 <- asBCdat(dat = data.frame(read_count = c(251, 20, 30, 12, 10, 4, 1),
                                   barcode = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                               "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", "AAAACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                               "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT", "AAAAACCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                               "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT")), label = "dat1")

  dat2 <- asBCdat(dat = data.frame(read_count = c(51, 40, 20, 16, 10, 1),
                                   barcode = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                               "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", "AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                               "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT", "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT"))
                  , label = "dat2")

  res <- list()
  res$shared <- data.frame(read_count = c(251, 30, 20, 10, 1), barcode = factor(x = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                 "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                                                                 "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                                                 "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT",
                                                                                 "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT"),
                                                                                levels = c("AAAAACCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                           "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                                                                           "AAAACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                           "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                                                           "ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                           "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                                                           "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT")),
                           read_count = c(51, 20, 40, 10,  1),
                           barcode = factor(x = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                       "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                       "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                       "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT",
                                       "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT"),
                                       levels = c("AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                  "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                                  "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                  "ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                  "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                  "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT")),
                           reads_diff = c(200,  10, -20,   0,   0))
  colnames(res$shared) <- c("read_count", "barcode", "read_count", "barcode", "reads_diff")
  attributes(res$shared)$row.names <- as.integer(c(1:3, 5, 7))

  res$unique_sample1 <- data.frame(read_count = c(12, 4), barcode = factor(x = c("AAAACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                 "AAAAACCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT"),
                                                                           levels = c("AAAAACCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                      "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                                                                      "AAAACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                      "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                                                      "ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                      "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                                                      "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT")))
  attributes(res$unique_sample1)$row.names <- as.integer(c(4, 6))

  res$unique_sample2 <- data.frame(read_count = 16, barcode = factor(x = "AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                     levels = c("AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                                                                "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                                                "ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                                                "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                                                "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT")))
  attributes(res$unique_sample2)$row.names <- as.integer(4)

  expect_error(com_pair(BC_dat1 = dat1), "# Two BCdat-objects required!")
  expect_error(com_pair(BC_dat1 = dat1, BC_dat2 = 1:5), "# Two BCdat-objects required!")

  test_dat <- com_pair(BC_dat1 = dat1, BC_dat2 = dat2)

  expect_equal(test_dat$shared, res$shared)
  expect_equal(test_dat$unique_sample1, res$unique_sample1)
  expect_equal(test_dat$unique_sample2, res$unique_sample2)

})

test_that("helper functions - rest", {

  expect_error(asBCdat(data.frame(), label = "empty", BC_backbone = "none", resDir = getwd()), "# Data object has to have two columns.")
  expect_error(asBCdat(data.frame(x = letters[1:5], y = letters[1:5]), label = "empty", BC_backbone = "none", resDir = getwd()), "# First column contains number of read counts therefore has to be of type numeric.")


  dat <- asBCdat(dat = data.frame(read_count = c(51, 40, 20, 16, 10, 1),
                                   barcode = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                               "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", "AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                               "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT", "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT"))
                  , label = "dat")

  expect_equal(.getMinDist(BC_dat = dat,
              ori_BCs = c("AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", "AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT"),
              m = "hamming"), c(4, 9, 0, 0, 6, 7))

  expect_equal(.getMinDist_one_sample(BC_dat = dat, m = "hamming"), c(3, 2, 3, 6, 10))

  #.revComp()

})
