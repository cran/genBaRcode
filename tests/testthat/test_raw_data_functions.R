
library(genBaRcode)
library(ShortRead)

context("raw data processing")

test_that("processingRawData - bc_backbone and file_name test", {

  eM <- "# Backbones are only valid if consisting of IUPAC-nucleotide-code symbols (backbone 1)"
  expect_error(processingRawData(file_name = "test_data.fastq",
                                source_dir = "./",
                                results_dir = "./",
                                bc_backbone = "X"), eM, fixed = TRUE)

  eM <- "# Backbones are only valid if consisting of IUPAC-nucleotide-code symbols (backbone 2)"
  expect_error(processingRawData(file_name = "test_data.fasta",
                                source_dir = "./",
                                results_dir = "./",
                                bc_backbone = c("ACTN", "ACTNL")), eM, fixed = TRUE)

  eM <- "# Unknown file format (test_data), so far only fastq and fasta files are supported."
  expect_error(processingRawData(file_name = "test_data",
                                source_dir = "./",
                                results_dir = "./",
                                bc_backbone = c("ACTN", "ACTNK")), eM, fixed = TRUE)

  eM <- "# Unknown file format (test_data.xy), so far only fastq and fasta files are supported."
  expect_error(processingRawData(file_name = "test_data.xy",
                                source_dir = "./",
                                results_dir = "./",
                                bc_backbone = c("ACTN", "ACTN")), eM, fixed = TRUE)
})

test_that("processingRawData - extractBarcodes", {

  bb <- "ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN"
  fastq <- c("ACGACTTTCGATTCTTTTCGTTTCTTTTGGATTCTATTACTTTCGATTCATTTCGATTCTTTTGGATTCTATTACTTTCGATTGCA", # 2x mismatch
             "CCCACTAACGAAACTTAACGAAACTTAAGGAAACTAAAACTAACGAAACTTAACGAAACTTAAGGAAACTAAAACTAACGAAAGGG", # AAA
             "AAAACTTTCGATTCTTTTCGATTCTTTTGGATTCTATTACTTTCGATTCTTTTCGATTCTTTTGGATTCTATTACTTTCGATTCCC", # TTT
             "AAAACTGGCGAGGCTTGGCGAGGCTTGGGGAGGCTAGGACTGGCGAGGCTTGGCGAGGCTTGGGGAGGCTAGGACTGGCGAGGTTT", # GGG
             "TTTACTCCCGACCCTTCCCGACCCTTCCGGACCCTACCACTCCCGACCCTTCCCGACCCTTCCGGACCCTACCACTCCCGACCAAA", # CCC
             "TCGACTTTCGATTGTTTTCGATTCTTTTGGATTCTATTACTTTCGATTCTTTTCGATTCTTTTGGATTCGATTACTTTCGATTTCG", # 2x mismatch
             "AAAACTGGCGAGGCTTGGCGAGGCTTGGGGAGGGTAGGACTCCCGACCCTTCCCGACCCTTCCGGACCCTACCACTCCCGACCTTT", # 1x mismatch
             "GTAAGTACCCCTGCTTTGCAGGCAAAGCAGTAGCTCCGAGGAGTAACACCACATCTGTCGACCCAATTTGTGACAACTCTAGGCAC") # random
  ids <- c("2MMa", "AA", "TT", "GG", "CC", "2MMb", "1MM", "random")

  dat <- methods::new(Class = "ShortRead",
                      sread = Biostrings::DNAStringSet(fastq),
                      id = Biostrings::BStringSet(ids))

  eM <- "# At the moment it is only possible to analyse the data either with no backbone at all or with one or multiple actual backbones."
  expect_error(extractBarcodes(dat = dat, label = "test", results_dir = "",
                  mismatch = 0, indels = FALSE,
                  bc_backbone = c("ACN", "none")), eM, fixed = TRUE)

  eM <- "# Barcodes can currently only be extracted from ShortRead data objects."
  expect_error(extractBarcodes(dat = data.frame(), label = "test", results_dir = "",
                  mismatch = 0, indels = FALSE,
                  bc_backbone = "ACN"), eM, fixed = TRUE)

  # maybe after "compilation"
  # wM <- "# mismatch needs to be an integer >= 0"
  # expect_warning(extractBarcodes(dat = dat, label = "test", results_dir = "",
  #                              mismatch = -1, indels = FALSE,
  #                              bc_backbone = "ACN"), wM, fixed = TRUE)

  res <- evaluate_promise(
            extractBarcodes(dat, label = "test", results_dir = "", mismatch = -1, indels = FALSE, bc_backbone = bb, wobble_extraction = FALSE)
  )

  expect_type(res$result, "list")
  expect_equal(dim(res$result[[1]])[1], 4)
  expect_equal(dim(res$result[[1]])[2], 2)
  expect_equal(sum((res$result[[1]])[, 2]), 4)

  wM <- "# mismatch needs to be an integer >= 0"
  expect_equal(res$warnings, wM)

  res <- extractBarcodes(dat, label = "test", results_dir = "", mismatch = 0, indels = FALSE, bc_backbone = bb, wobble_extraction = FALSE)

  expect_type(res, "list")
  expect_equal(dim(res[[1]])[1], 4)
  expect_equal(dim(res[[1]])[2], 2)
  expect_equal(sum((res[[1]])[, 2]), 4)

  res <- extractBarcodes(dat, label = "test", results_dir = "", mismatch = 1, indels = FALSE, bc_backbone = bb, wobble_extraction = FALSE)

  expect_type(res, "list")
  expect_equal(dim(res[[1]])[1], 5)
  expect_equal(dim(res[[1]])[2], 2)
  expect_equal(sum((res[[1]])[, 2]), 5)

  res <- extractBarcodes(dat, label = "test", results_dir = "", mismatch = 2, indels = FALSE, bc_backbone = bb, wobble_extraction = FALSE)

  expect_type(res, "list")
  expect_equal(dim(res[[1]])[1], 7)
  expect_equal(dim(res[[1]])[2], 2)
  expect_equal(sum((res[[1]])[, 2]), 7)

  res <- extractBarcodes(dat, label = "test", results_dir = "", mismatch = 3, indels = FALSE, bc_backbone = bb, wobble_extraction = FALSE)

  expect_type(res, "list")
  expect_equal(dim(res[[1]])[1], 7)
  expect_equal(dim(res[[1]])[2], 2)
  expect_equal(sum((res[[1]])[, 2]), 7)

})

test_that("processingRawData - fixBorderlineBCs and extractWobble", {

  bb <- "ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN"
  fastq <- c("CTAACGAAACTTAACGAAACTTAAGGAAACTAAAACTAACGAAACTTAACGAAACTTAAGGAAACTAAAACTAACGAAAGGGCGTAG",    # AAA - 1mm in the front
             "CCCACTAACGAAACTTAACGAAACTTAAGGAAACTAAAACTAACGAAACTTAACGAAACTTAAGGAAACTAAAACTAACGAAAGGGT",    # AAA
             "AGCGCTAGCCCACTAACGAAACTTAACGAAACTTAAGGAAACTAAAACTAACGAAACTTAACGAAACTTAAGGAAACTAAAACTAAC",    # AAA - 4mm in the back
             "ACTTTCGATTCTTTTCGATTCTTTTGGATTCTATTACTTTCGATTCTTTTCGATTCTTTTGGATTCTATTACTTTCGATTCCCAGGG",    # TTT
             "AAAACTGGCGAGGCTTGGCGAGGCTTGGGGAGGCTAGGACTGGCGAGGCTTGGCGAGGCTTGGGGAGGCTAGGACTGGCGAGGTTTA",    # GGG
             "TGCTCAAAACTGGCGAGGCTTGGCGAGGCTTGGGGAGGCTAGGACTGGCGAGGCTTGGCGAGGCTTGGGGAGGCTAGGACTGGCGAG",    # GGG - 1mm in the back
             "AGTCCCCAAAACTGGCGAGGCTTGGCGAGGCTTGGGGAGGCTAGGACTGGCGAGGCTTGGCGAGGCTTGGGGAGGCTAGGACTGGCG",    # GGG - 1 fixed pos mm in the back (3mm total)
             "TTTACTCCCGACCCTTCCCGACCCTTCCGGACCCTACCACTCCCGACCCTTCCCGACCCTTCCGGACCCTACCACTCCCGACCAAAG",    # CCC
             "TCCCGACCCTTCCCGACCCTTCCGGACCCTACCACTCCCGACCCTTCCCGACCCTTCCGGACCCTACCACTCCCGACCAAACACGGT",    # CCC - 2mm in the front
             "TTCGATTCTTTTCGATTCTTTTGGATTCTATTACTTTCGATTCTTTTCGATTCTTTTGGATTCTATTACTTTCGATTCCCGTGGACT",    # TTT - 3mm in the front
             "GTAAGTACCCCTGCTTTGCAGGCAAAGCAGTAGCTCCGAGGAGTAACACCACATCTGTCGACCCAATTTGTGACAACTCTAGGCACC")    # random

  ids <- c("AA-1mmF", "AA", "AA-3mmB", "TT", "GG", "GG-1mmB", "GG-3mmB", "CC", "CC-2mmF", "TT-3mmF", "random")

  dat <- methods::new(Class = "ShortRead",
                      sread = Biostrings::DNAStringSet(fastq),
                      id = Biostrings::BStringSet(ids))
  tmp <- list(list(), ShortRead::sread(dat))

  read_length <- nchar(fastq)[1]

  test_res <- c(6, 7, 9)
  test2_res <- list()
  test2_res[[1]] <- data.frame('dat...1.' = c("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGN", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
                                n = c(2, rep(1, 4)))
  test2_res[[2]] <- data.frame('dat...1.' = c("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGN", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
                                n = c(rep(2, 2), rep(1, 3)))
  test2_res[[3]] <- data.frame('dat...1.' = c("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                                              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGN", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGNN"),
                                n = c(rep(2, 3), rep(1,3)))

  for (mm in 1:3) {
    # test fixBorderlineBCs
    match_index <- Biostrings::vmatchPattern(pattern = bb,
                                             subject = Biostrings::DNAStringSet(fastq),
                                             max.mismatch = mm, min.mismatch = 0, with.indels = FALSE, fixed = FALSE, algorithm = "auto")
    match_matrix <- as.data.frame(match_index)
    BCs <- substr(as.character(fastq)[match_matrix$group], start = match_matrix$start, stop = match_matrix$end)
    res <- fixBorderlineBCs(BCs, mismatch = mm, match_matrix, read_length)

    expect_equal(length(res), test_res[mm])

    # test extractWobble
    res <- BC_matching(tmp, "test", results_dir = "", mismatch = mm, indels = FALSE, bc_backbone = bb, full_output = FALSE)
    res <- data.frame(extractWobble(i = 1, tmp = res, bc_backbone = bb, indels = FALSE, label = ""))

    expect_equal(res, test2_res[[mm]])
  }

})
