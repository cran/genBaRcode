
context("errorCorrection")

test_that("errorCorrection", {

  dat <- asBCdat(
          data.frame(read_count = c(10, 7, 2, 1, 1,  12, 3, 1, 30, 20, 10, 2, 19, 14, 5, 1),
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
                                "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT"), stringsAsFactors = TRUE
                    )
        )

  res <- data.frame(read_count = c(61, 41, 20, 16),
                    barcode = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                              "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", "AAAACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT"),
                    stringsAsFactors = TRUE)
  test_res <- getReads(errorCorrection_single_variation(BC_dat = dat, maxDist = 2, save_it = FALSE, m = "hamming", EC_analysis = FALSE, nt = 1))
  expect_equal(res, test_res)

  res <- data.frame(read_count = c(51, 40, 20, 16, 10, 1),
                    barcode = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", "AAAACCCCGGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT", "AAAAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT"), stringsAsFactors = TRUE)
  test_res <- getReads(errorCorrection_single_connections(BC_dat = dat, maxDist = 2, save_it = FALSE, m = "hamming", EC_analysis = FALSE, nt = 1))
  expect_equal(res, test_res)

  res <- data.frame(read_count = c(97, 41),
                    barcode = c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT", "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT"), stringsAsFactors = TRUE)
  test_res <- getReads(errorCorrection_single_graphComp(BC_dat = dat, maxDist = 2, save_it = FALSE, m = "hamming", EC_analysis = FALSE, nt = 1))
  expect_equal(res, test_res)

  res <- matrix(c(c(51, 41, 20, 16, 10),
                  c("ATTAGCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                "ATTAGCCCTGGGATTTAGGACTTCCGGGTATTAAAACCCCGGGGATTT",
                                "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
                                "AAAACCCCTGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGATTT",
                                "ATTAGCCGCGGGATTTAAAACCCCCGGGTTTTAAAACCCCGGGGAATT")), ncol = 2)
  test_res <- getReads(errorCorrection_single_clustering_absolute(BC_dat = dat, maxDist = 2, save_it = FALSE, m = "hamming", EC_analysis = FALSE, nt = 1))
  test_res <- matrix(c(as.character(test_res$read_count), as.character(test_res$barcode)), ncol = 2)
  expect_equal(res, test_res)

  })
