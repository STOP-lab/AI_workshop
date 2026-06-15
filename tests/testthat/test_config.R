test_that("experiment template can be read", {
  source("R/config.R")
  cfg <- read_experiment_config("config/experiment_template.yml")
  expect_equal(length(cfg$conditions), 2)
})
