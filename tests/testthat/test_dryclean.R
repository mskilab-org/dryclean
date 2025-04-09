
context("unit testing dryclean operations")

library(dryclean)
library(GenomicRanges)
library(dplyr)

sample.1.path = system.file("extdata", "samp1.rds", package = 'dryclean')
sample.2.path = system.file("extdata", "samp2.rds", package = 'dryclean')
sample.3.path = system.file("extdata", "samp3.rds", package = 'dryclean')

decomp.1.path = system.file("extdata", "decomp1.rds", package = 'dryclean')
decomp.2.path = system.file("extdata", "decomp2.rds", package = 'dryclean')
decomp.3.path = system.file("extdata", "decomp3.rds", package = 'dryclean')

normal_vector = c(sample.1.path,sample.2.path,sample.3.path)

normal_table = data.table(sample = c("samp1", "samp2", "samp3"), normal_cov = c(sample.1.path, sample.2.path, sample.3.path), decomposed_cov = c(decomp.1.path, decomp.2.path, decomp.3.path))

detergent.path = system.file("extdata", "detergent.rds", package = 'dryclean')

detergent_test.path = system.file("extdata", "detergent_test.rds", package = 'dryclean')

batch_outputs.path = system.file("extdata", "batch_outputs.rds", package = 'dryclean')
U.path = system.file("extdata", "U.rds", package = 'dryclean')
m.vec.path = system.file("extdata", "m.vec.rds", package = 'dryclean')

## Tests

test_that("prep_cov", {
  field = "reads.corrected"
  sample.1 = readRDS(sample.1.path)
  sample.1 = sample.1[, field] %>% gr2dt() %>% setnames(., field, "signal") %>% dt2gr()
  pcov = prep_cov(m.vec = sample.1)
  expect_identical(class(pcov)[1], "GRanges")
  expect_equal(length(pcov)[1], 50)
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


test_that("initialize", {
  
  expect_error(
    pon_object = pon$new()
  )
  
  expect_error(
    pon_object = pon$new(pon_path = "WRONG_PATH")
  )
  
  expect_error(
    pon_object = pon$new(create_new_pon = TRUE, save_pon = TRUE)
  )
  
  expect_error(
    normal_samples = c("WRONG_FILE"),
    pon_object = pon$new(create_new_pon = TRUE, normal_vector = normal_samples)
  )
  
  expect_error(
    normal_samples = c(),
    pon_object = pon$new(create_new_pon = TRUE, normal_vector = normal_samples)
  )
  
  pon_object = pon$new(
    create_new_pon = TRUE, 
    normal_vector = normal_vector,
    field = "reads.corrected",
    tolerance = 0.4, target_resolution = 3
    )
  
  expect_equal(dim(pon_object$get_L())[1], 50)
  expect_equal(dim(pon_object$get_S())[1], 50)
  expect_equal(length(pon_object$get_err()),2)
  expect_equal((pon_object$get_k())[1],2)
  expect_equal(dim(pon_object$get_U_hat())[1], 50)
  expect_equal(dim(pon_object$get_V_hat())[1], 2)
  expect_equal(length(pon_object$get_sigma_hat()), 2)
  expect_equal(names(pon_object$get_seqlengths()), "22")
  expect_equal(length(pon_object$get_template()), 50)
  
  dryclean_object = dryclean$new(pon = pon_object)
  expect_true(any(grepl("Created dryclean object", capture.output(dryclean_object$get_history()))))
  
  pon_object = pon$new(
    create_new_pon = TRUE, 
    normal_vector = normal_vector,
    save_pon = TRUE,
    pon_path = detergent_test.path,
    field = "reads.corrected",
    tolerance = 0.4, target_resolution = 3
  )

  expect_equal(readRDS(detergent_test.path)$L[1], 0.06105891, tolerance = 0.00001)
  
})


test_that("clean", {
  
  pon_object = pon$new(pon_path = detergent.path)
  
  dryclean_object <- dryclean$new(pon = pon_object)
  
  a <- dryclean_object$clean(cov = sample.1.path, testing = TRUE)
  expect_true(identical(colnames(values(a)), c("background.log", "foreground.log", "signal", "input.read.counts", "center.all", "median.chr", "foreground", "background", "log.reads")))
  
  expect_equal(a$background.log[1], 0.0508, tolerance = 0.001)
})


test_that("get_history", {
  
  pon_object = pon$new(
    create_new_pon = TRUE, 
    normal_vector = normal_vector,
    field = "reads.corrected",
    tolerance = 0.4, target_resolution = 3
  )
  
  dryclean_object <- dryclean$new(pon = pon_object)
  
  a <- dryclean_object$clean(cov = sample.1.path, testing = TRUE)
  
  expect_true(any(grepl("Created dryclean object", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Loaded PON", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Loaded coverage from", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Started drycleaning the coverage file", capture.output(dryclean_object$get_history()))))
  expect_true(any(grepl("Finished drycleaning the coverage file", capture.output(dryclean_object$get_history()))))
  
  expect_true(any(grepl("Created pon object", capture.output(pon_object$get_history()))))
  expect_true(any(grepl("Loaded normal vector", capture.output(pon_object$get_history()))))
  expect_true(any(grepl("Started PON preparation", capture.output(pon_object$get_history()))))
  expect_true(any(grepl("Created new PON from the normal samples", capture.output(pon_object$get_history()))))
})


test_that("cbs", {
  
  pon_object = pon$new(pon_path = detergent.path)
  
  dryclean_object <- dryclean$new(pon = pon_object)
  
  a <- dryclean_object$clean(cov = sample.1.path, testing = TRUE, cbs = TRUE)
  
  cbs <- a@seqinfo@seqlengths

  expect_true(identical(as.numeric(cbs[1]),101))
  
  expect_true(any(grepl("Applied CBS correction to the drycleaned coverage file", capture.output(dryclean_object$get_history()))))
  
})


test_that("clustering", {
  expect_error(
    pon_object = pon$new(
      create_new_pon = TRUE, 
      normal_vector = normal_vector,
      field = "reads.corrected",
      tolerance = 0.4, target_resolution = 3,
      choose.by.clustering = TRUE,
      number.of.samples = 2
    )
  )
  
  pon_object = pon$new(
    create_new_pon = TRUE, 
    normal_vector = normal_vector,
    field = "reads.corrected",
    tolerance = 0.4, target_resolution = 3,
    use.all = FALSE,
    choose.by.clustering = TRUE,
    number.of.samples = 2
  )

  expect_equal(pon_object$get_L()[5], 0.004731082, tolerance = 0.00001) 
  
})

test_that("infer.germline", {

  pon_object = pon$new(
    create_new_pon = TRUE, 
    normal_vector = normal_vector,
    field = "reads.corrected",
    tolerance = 0.4, target_resolution = 3,
    infer.germline = TRUE
  )
  
  expect_equal(length(pon_object$get_inf_germ()), 50) 
})
