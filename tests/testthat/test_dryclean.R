
context("unit testing dryclean operations")

library(dryclean)
library(GenomicRanges)


sample.1.path = system.file("extdata", "samp1.rds", package = 'dryclean')
sample.2.path = system.file("extdata", "samp2.rds", package = 'dryclean')
sample.3.path = system.file("extdata", "samp3.rds", package = 'dryclean')


normal_table.path  = system.file("extdata", "normal_table.rds", package = 'dryclean')

batch_outputs.path = system.file("extdata", "batch_outputs.rds", package = 'dryclean')
U.path = system.file("extdata", "U.rds", package = 'dryclean')
m.vec.path = system.file("extdata", "m.vec.rds", package = 'dryclean')

## Tests

batch_outputs = list()
batch_outputs$L = matrix(runif(100*50, -1, 1), 100, 50)
batch_outputs$S = matrix(runif(100*50, -1, 1), 100, 50)
batch_outputs$k = 10
batch_outputs$U.hat = matrix(runif(100*50, -1, 1), 100, 10)
batch_outputs$V.hat = matrix(runif(100*50, -1, 1), 10, 50)
batch_outputs$sigma.hat = sort(runif(10, 0, 100), decreasing = T)

U = matrix(runif(100*10), 100, 10)
m.vec = matrix(runif(100, -1, 1))

samp = GRanges(1, IRanges(c(rep(1, 100)), c(rep(10, 100))), strand=c(rep("*", 100)), reads.corrected = c(runif(100, 0, 5)))

normal_table = system.file("data", "normal_table.rds", package = 'dryclean')

samples_path = system.file("data", package = 'dryclean')


test_that("prep_cov", {
    sample.1 = readRDS(sample.1.path)
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
    expect_error(prepare_detergent(normal_table.path))
    all.smp = prepare_detergent(normal.table.path = normal_table.path, path.to.save = "https://github.com/mskilab/dryclean/tree/master/inst/extdata/")
    expect_true(identical(names(all.smp), c("L", "S","err", "k", "U.hat", "V.hat", "sigma.hat")))
    hclust.smp = prepare_detergent(normal.table.path = normal_table.path, path.to.save = "https://github.com/mskilab/dryclean/tree/master/inst/extdata/", use.all = FALSE, choose.by.clustering = TRUE, number.of.samples = 3)
    expect_true(identical(names(hclust.smp), c("L", "S","err", "k", "U.hat", "V.hat", "sigma.hat")))
})


test_that("identify_germline", {
    expect_error(identify_germline(normal_table.path))
    all.smp = identify_germline(normal.table.path = normal_table.path, path.to.save = "https://github.com/mskilab/dryclean/tree/master/inst/extdata/")
})


test_that("start_wash_cycle", {
    sample.1 = readRDS(sample.1.path)
    expect_error(start_wash_cycle(sample.1))
    strt = start_wash_cycle(cov = sample.1, detergent.pon.path = "https://github.com/mskilab/dryclean/tree/master/inst/extdata/", whole_genome = TRUE, chr = NA, germline.filter = FALSE)
    expect_true(identical(colnames(values(strt)), c("L", "reads.corrected.log","reads.corrected.FC", "median.chr", "reads.corrected", "L1", "log.reads")))
    expect_error(start_wash_cycle(samp))
    strt = start_wash_cycle(cov = samp, burnin.samples.path = samples_path, whole_genome = FALSE, chr = "1")
})






