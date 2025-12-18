library(testthat)
local_edition(3)

test_that("BHiCect", {
  # basic example using chromosome 22 data from human GM12878 data (Rao et al. 2014)
  data(chr_dat_l)

  # manual entry for resolutions
  res_set <- c("1Mb", "500kb", "100kb")
  res_num <- c(1e6L, 5e5L, 1e5L)
  names(res_num) <- res_set

  BHiCect_results <- BHiCect(res_set, res_num, chr_dat_l,'smpl.cl', 4)

  # TODO implement some basic validation on the produced output
  cl_res_count <- table(purrr::map_chr(strsplit(names(BHiCect_results$cl_member),'_'),function(x) x[1]))
  expect_true(cl_res_count['1Mb'] < cl_res_count['100kb'])
})
