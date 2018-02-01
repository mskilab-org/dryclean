library(dryclean)
library(GenomicRanges)


context("unit testing dryclean operations")

batch_outputs = list()
batch_outputs$L = matrix(runif(100*50, -1, 1), 100, 50)
batch_outputs$S = matrix(runif(100*50, -1, 1), 100, 50)
batch_outputs$k = 10
batch_outputs$U.hat = matrix(runif(100*50, -1, 1), 100, 10)
batch_outputs$V.hat = matrix(runif(100*50, -1, 1), 10, 50)
batch_outputs$sigma.hat = sort(runif(10, 0, 100), decreasing = T)

samp = GRanges(1, IRanges(c(rep(1, 100)), c(rep(10, 100))), strand=c(rep("*", 100)), reads.corrected = c(runif(100, 0, 100)))


U = matrix(runif(100*10), 100, 10)
m.vec = as.matrix(runif(100, -1, 1))



test_that("prep_cov", {
    pcov = prep_cov(samp)
    expect_identical(class(pcov)[1], "data.table")
    expect_equal(dim(pcov)[1], 100)
})

test_that("apg_project", {
    proj = apg_project(m.vec = m.vec, U = U, lambda1 = 0.0001, lambda2 = 0.0001)
    expect_equal(dim(proj[[1]])[1], dim(U)[2])
    expect_equal(dim(proj[[2]])[1], dim(U)[1])
})

test_that("wash_cycle", {
    samp.test = wash_cycle(m.vec = m.vec, L.burnin = batch_outputs$L, S.burnin = batch_outputs$S, r = batch_outputs$k, U.hat = batch_outputs$U.hat, V.hat = batch_outputs$V.hat, sigma.hat = batch_outputs$sigma.hat, verbose = FALSE)
    expect_equal(dim(samp.test[[1]])[1], dim(m.vec)[1])
    expect_equal(length(samp.test[[2]]), dim(m.vec)[1])
})

test_that("start_wash_cycle", {
    expect_error(start_wash_cycle(samp))
    strt = start_wash_cycle(cov = samp, burnin.samples.path = "~/git/dryclean/data", whole_genome = FALSE, chr = "1")
})



    
    
