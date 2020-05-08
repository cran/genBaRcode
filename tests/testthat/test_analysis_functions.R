
context("analysis functions")

test_that("analysis functions", {

  dat1 <- asBCdat(dat = data.frame(read_count = c(10, 7, 2, 1, 1,  12, 3, 1, 30, 20, 10, 2, 19, 14, 5, 1),
                                   barcode = c("AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                               "AAAAACCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                               "AAAAACCCGGGGATTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                               "AAAACCCCAGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                               "AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                               "AAAACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                               "AAAAACCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                               "ATTACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                               "ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                               "ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT",
                                               "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT",
                                               "ATTAGCCCTGGGATTTAGGACTTCCGGGTTTTAAAACCCCGGGGATTT",
                                               "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                               "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACGCCGGGGATTT",
                                               "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACGCCGGGGATTT",
                                               "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT")
                          ), label = "dat1"
                      )

  dat2 <- asBCdat(dat = data.frame(read_count = c(251, 20, 30, 12, 10, 4, 1),
                                    barcode = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", "AAAACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT", "AAAAACCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT")), label = "dat2")

  dat3 <- asBCdat(dat = data.frame(read_count = c(51, 40, 20, 16, 10, 1),
                                    barcode = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                                "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", "AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                                "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT", "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT"))
                  , label = "dat3")

  seqs <- c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT",
            "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACGCCGGGGATTT",
            "AAAACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT", "AAAAACCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "AAAAACCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "AAAAACCCGGGGATTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "ATTAGCCCTGGGATTTAGGACTTCCGGGTTTTAAAACCCCGGGGATTT", "AAAACCCCAGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
            "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT")

  c1 <- c(30, 20, 19, 14, 12, 10, 10,  7,  3,  2,  2,  1,  1,  1,  1)
  c2 <- c(251,   0,  20,   0,  12,  30,  10,   0,   4,   0,   0,   0,   0,   0,   1)
  c3 <- c(51,  0, 40,  0,  0, 20, 10,  0,  0,  0,  0,  0, 16,  0,  1)

  res <- matrix(c(c1, c2, c3), ncol = 3)
  colnames(res) <- c("dat1", "dat2", "dat3")
  rownames(res) <- seqs

  expect_error(generateTimeSeriesData(BC_dat_list = list(dat1, 1, "2", NA)), "# All list elements need to be of type BCdat.")

  expect_error(generateTimeSeriesData(BC_dat_list = list()), "# data object needs to be at least a list of length = 2")
  expect_error(generateTimeSeriesData(BC_dat_list = list(dat1), "# data object needs to be at least a list of length = 2"))

  expect_equal(generateTimeSeriesData(BC_dat_list = list(dat1, dat2, dat3)), res)

})
