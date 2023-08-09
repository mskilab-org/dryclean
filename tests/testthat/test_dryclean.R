
context("unit testing dryclean operations")

library(dryclean)
library(GenomicRanges)

sample.1.path = system.file("extdata", "samp1.rds", package = 'dryclean')
sample.2.path = system.file("extdata", "samp2.rds", package = 'dryclean')
sample.3.path = system.file("extdata", "samp3.rds", package = 'dryclean')

decomp.1.path = system.file("extdata", "decomp1.rds", package = 'dryclean')
decomp.2.path = system.file("extdata", "decomp2.rds", package = 'dryclean')
decomp.3.path = system.file("extdata", "decomp3.rds", package = 'dryclean')

normal_table.path  = system.file("extdata", "normal_table.rds", package = 'dryclean')

normal_table = data.table(sample = c("samp1", "samp2", "samp3"), normal_cov = c(sample.1.path, sample.2.path, sample.3.path), decomposed_cov = c(decomp.1.path, decomp.2.path, decomp.3.path))


detergent.path = system.file("extdata", "detergent.rds", package = 'dryclean')



batch_outputs.path = system.file("extdata", "batch_outputs.rds", package = 'dryclean')
U.path = system.file("extdata", "U.rds", package = 'dryclean')
m.vec.path = system.file("extdata", "m.vec.rds", package = 'dryclean')

## Tests


test_that("prep_cov", {
  field = "reads.corrected"
  sample.1 = readRDS(sample.1.path)
  sample.1 = sample.1[, field] %>% gr2dt() %>% setnames(., field, "signal") %>% dt2gr()
  pcov = prep_cov(m.vec = sample.1)
  expect_identical(class(pcov)[1], "data.table")
  expect_equal(dim(pcov)[1], 50)
})


test_that("apg_project", {
  m.vec = readRDS(m.vec.path)
  U = readRDS(U.path)
  proj = apg_project(m.vec = m.vec, U = U, lambda1 = 0.0001, lambda2 = 0.0001)
  expect_equal(dim(proj[[1]])[1], dim(U)[2])
  expect_equal(dim(proj[[2]])[1], dim(U)[1])
})


test_that("wash_cycle", {
  batch_outputs = readRDS(batch_outputs.path)
  m.vec = readRDS(m.vec.path)
  samp.test = wash_cycle(m.vec = m.vec, L.burnin = batch_outputs$L, S.burnin = batch_outputs$S, r = batch_outputs$k, U.hat = batch_outputs$U.hat, V.hat = batch_outputs$V.hat, sigma.hat = batch_outputs$sigma.hat, verbose = FALSE)
  expect_equal(dim(samp.test[[1]])[1], dim(m.vec)[1])
  expect_equal(length(samp.test[[2]]), dim(m.vec)[1])
})


test_that("prepare_detergent", {
  dryclean_object = dryclean$new(
    pon_path = detergent.path, 
    cov_path = sample.1.path)
  expect_error(dryclean_object$prepare_detergent())
  
  dryclean_object = dryclean$new(
    normal_table_path = normal_table.path, 
    cov_path = sample.1.path)
  expect_error(dryclean_object$prepare_detergent(save.pon = TRUE))

  dryclean_object = dryclean$new(
    pon_path = detergent.path,
    normal_table_path = normal_table.path, 
    cov_path = sample.1.path)
  expect_message(dryclean_object$prepare_detergent(save.pon = TRUE))
  
  dryclean_object = dryclean$new(
    normal_table_path = normal_table.path, 
    cov_path = sample.1.path)
  dryclean_object$prepare_detergent(tolerance = 0.4, target_resolution = 3)
  all.smp = dryclean_object$get_pon()
  expect_true(identical(names(all.smp), c("L", "S","err", "k", "U.hat", "V.hat", "sigma.hat", "template", "inf_germ")))
  dryclean_object$prepare_detergent(use.all = FALSE, choose.by.clustering = TRUE, number.of.samples = 3, tolerance = 0.4, target_resolution = 3)
  hclust.smp = dryclean_object$get_pon()
  expect_true(identical(names(hclust.smp), c("L", "S","err", "k", "U.hat", "V.hat", "sigma.hat", "template", "inf_germ")))
})


test_that("start_wash_cycle", {
  dryclean_object = dryclean$new(
    pon_path = detergent.path, 
    cov_path = sample.1.path)
  
  dryclean_object$start_wash_cycle()

  strt = dryclean_object$get_drycleaned_cov()
  expect_true(identical(colnames(values(strt)), c("background.log", "foreground.log","input.read.counts", "median.chr", "foreground", "background", "log.reads")))
})

test_that("initialize", {
  expect_error(dryclean$new(
    cov_path = sample.1.path
  ))
  expect_error(dryclean$new(
    pon_path = "WRONG_PATH",
    cov_path = sample.1.path
  ))
})

test_that("get_normal_table_path", {
  dryclean_object = dryclean$new(
    normal_table_path = normal_table.path, 
    cov_path = sample.1.path)
  
  expect_true(identical(dryclean_object$get_normal_table_path(),normal_table.path))
})


test_that("get_normal_table", {
  dryclean_object = dryclean$new(
    normal_table_path = normal_table.path, 
    cov_path = sample.1.path)
  
  expect_true(identical(dryclean_object$get_normal_table(),readRDS(normal_table.path)))
})

test_that("get_cov_path", {
  dryclean_object = dryclean$new(
    pon_path = detergent.path,
    cov_path = sample.1.path)
  
  expect_true(identical(dryclean_object$get_cov_path(), sample.1.path))
})


test_that("get_pon_path", {
  dryclean_object = dryclean$new(
    pon_path = detergent.path,
    cov_path = sample.1.path)
  
  expect_true(identical(dryclean_object$get_pon_path(), detergent.path))
})

test_that("get_history", {
  dryclean_object = dryclean$new(
    normal_table_path = normal_table.path,
    cov_path = sample.1.path)
  
  dryclean_object$prepare_detergent()
  dryclean_object$start_wash_cycle()
  
  expect_true(any(grepl("Created dryclean object", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Loaded normal table from", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Started PON preparation", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Loaded coverage", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Created new PON from the normal table", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Loaded coverage from", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Finished drycleaning the coverage file", capture.output(dryclean_object$get_history()))))
})


test_that("save_drycleaned_cov", {
  dryclean_object = dryclean$new(
    pon_path = detergent.path,
    cov_path = sample.1.path)
  
  expect_error(dryclean_object$save_drycleaned_cov())
  
  dryclean_object$start_wash_cycle()
  
  expect_error(dryclean_object$save_drycleaned_cov())
})

test_that("save_cbs_drycleaned_cov", {
  dryclean_object = dryclean$new(
    pon_path = detergent.path,
    cov_path = sample.1.path)
  
  expect_error(dryclean_object$save_cbs_drycleaned_cov())
  
  dryclean_object$start_wash_cycle()
  
  expect_error(dryclean_object$save_cbs_drycleaned_cov())
  
  dryclean_object$cbs()
  
  expect_error(dryclean_object$save_cbs_drycleaned_cov())
  
})



test_that("cbs", {
  dryclean_object = dryclean$new(
    pon_path = detergent.path,
    cov_path = sample.1.path)
  
  expect_error(dryclean_object$cbs())
  
  dryclean_object$start_wash_cycle()
  
  dryclean_object$cbs()  
  
  cbs <- dryclean_object$get_cbs_drycleaned_cov()@seqinfo@seqlengths

  expect_true(identical(as.numeric(cbs[1]),101))
  
  expect_true(any(grepl("segments produced", capture.output(dryclean_object$cbs()))))
  expect_true(any(grepl("finished segmenting", capture.output(dryclean_object$cbs()))))
  
  expect_true(any(grepl("Applied CBS correction to the drycleaned coverage file", capture.output(dryclean_object$get_history()))))
  
  expect_true(any(grepl("The drycleaned coverage file corrected by CBS", capture.output(dryclean_object$check_status()))))
})

test_that("check_status", {
  dryclean_object = dryclean$new(
    pon_path = detergent.path,
    cov_path = sample.1.path)
  
  expect_true(any(grepl("The dryclean class object is created.", capture.output(dryclean_object$check_status()))))
})


test_that("get_cbs_drycleaned_cov", {
  dryclean_object = dryclean$new(
    pon_path = detergent.path,
    cov_path = sample.1.path)
  
  dryclean_object$start_wash_cycle()
  
  dryclean_object$cbs()  
  
  expect_true(identical(as.vector(dryclean_object$get_cbs_drycleaned_cov()@seqnames)[1],"22"))
})
