library(testthat)
library(dryclean)

test_check("dryclean")
test_check(devtools::load_all())
