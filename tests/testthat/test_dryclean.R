
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

print(normal_table)

saveRDS(normal_table, normal_table.path)

detergent.path = system.file("extdata", "detergent.rds", package = 'dryclean')



batch_outputs.path = system.file("extdata", "batch_outputs.rds", package = 'dryclean')
U.path = system.file("extdata", "U.rds", package = 'dryclean')
m.vec.path = system.file("extdata", "m.vec.rds", package = 'dryclean')

## Tests


test_that("prep_cov", {
    sample.1 = readRDS(sample.1.path)
    pcov = prep_cov(m.vec = sample.1, build = "hg19")
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
    expect_error(prepare_detergent(normal_table.path, save.pon = TRUE))
    all.smp = prepare_detergent(normal.table.path = normal_table.path, tolerance = 0.4)
    expect_true(identical(names(all.smp), c("L", "S","err", "k", "U.hat", "V.hat", "sigma.hat")))
    hclust.smp = prepare_detergent(normal.table.path = normal_table.path, use.all = FALSE, choose.by.clustering = TRUE, number.of.samples = 3, tolerance = 0.4)
    expect_true(identical(names(hclust.smp), c("L", "S","err", "k", "U.hat", "V.hat", "sigma.hat")))
})


test_that("identify_germline", {
    expect_error(identify_germline(normal_table.path, save.grm = TRUE))
    all.smp = identify_germline(normal.table.path = normal_table.path)
    expect_true(colnames(values(all.smp)) == "germline.status")
})


test_that("start_wash_cycle", {
    sample.1 = readRDS(sample.1.path)
    expect_error(start_wash_cycle(sample.1))
    strt = start_wash_cycle(cov = sample.1, detergent.pon.path = detergent.path, whole_genome = TRUE, chr = NA, germline.filter = FALSE)
    ##print(colnames(values(strt)))
    expect_true(identical(colnames(values(strt)), c("background.log", "foreground.log","input.read.counts", "median.chr", "foreground", "background", "log.reads")))
    expect_error(start_wash_cycle(samp))
})






