library(testthat)
library(assertthat)
library(tidyverse)
source("R/gdid.R")

data <- data.frame(id = c(1, 1, 2, 2),
                   time = c(1, 2, 1, 2),
                   post = c(0, 1, 0, 1),
                   outcome = c(0, 1, 0, 1),
                   treatment = c(0, 1, 0, 1),
                   cluster = c(1, 1, 2, 2))

data2 <- data.frame(id = c(1, 1, 1, 2, 2, 2),
                   time = c(1, 2, 3, 1, 2, 3),
                   post = c(0, 1, 1, 0, 1, 1),
                   outcome = c(0, 1, 1, 0, 1, 1),
                   treatment = c(0, 0, 0, 1, 1, 1),
                   cluster = c(1, 1, 1, 2, 2, 2))

data3 <- data.frame(id = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
                    time = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
                    post = c(0, 0, 0, 1, 1, 0, 0, 0, 1, 1),
                    outcome = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                    treatment = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
                    cluster = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
                    Xt = seq(10, 19, 1),
                    X  = c(rep(5, 5), rep(10, 5)))

data4 <- data.frame(id = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
                    time = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
                    post = c(0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
                    outcome = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                    treatment = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0),
                    cluster = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3))

# Test that function returns a list

test_that("GDiD.DataProcess returns a list", {
  expect_is(GDiD.DataProcess(data,
                             "id",
                             "outcome",
                             "treatment",
                             time_invariant_vars = NULL,
                             number_lags = 1,
                             periods_post = 1), "list")})

# Test that function returns correct number of dataframes
test_that("GDiD.DataProcess returns correct number of dataframes",
          {expect_length(GDiD.DataProcess(data,
                                          "id",
                                          "outcome",
                                          "treatment",
                                          time_invariant_vars = NULL,
                                          number_lags = 1,
                                          periods_post = 1), 2)
            expect_length(GDiD.DataProcess(data,
                                           "id",
                                           "outcome",
                                           "treatment",
                                           time_invariant_vars = NULL,
                                           number_lags = 1,
                                           periods_post = 2), 2)})

# Test that function returns balanced panel
test_that("GDiD.DataProcess returns a balanced panel", {
  expect_error(GDiD.DataProcess(data %>% mutate(id = c(1, 2, 3, 3)),
                                "id",
                                "outcome",
                                "treatment",
                                cluster_var = "cluster",
                                time_invariant_vars = NULL,
                                number_lags = 1,
                                periods_post = 1),
               "Dataframe must be a balanced panel")})

# Test that function returns dataframes with correct columns
test_that("GDiD.DataProcess returns dataframes with correct columns", {
  expect_equal(colnames(GDiD.DataProcess(data, "id", "outcome", "treatment",
                                         time_invariant_vars = NULL,
                                         number_lags = 1,
                                         periods_post = 1)[[1]]),
               c("id", "cluster_id", "Y", "A"))
  expect_equal(colnames(GDiD.DataProcess(data, "id", "outcome", "treatment",
                                         time_invariant_vars = NULL,
                                         number_lags = 1,
                                         periods_post = 1)[[2]]),
               c("id", "cluster_id", "Y", "A"))})

# Test that function returns dataframes with correct number of rows
test_that("GDiD.DataProcess returns dataframes with correct number of rows",
          {expect_equal(nrow(GDiD.DataProcess(data, "id", "outcome", "treatment",
                                              time_invariant_vars = NULL,
                                              number_lags = 1,
                                              periods_post = 1)[[1]]), 2)
           expect_equal(nrow(GDiD.DataProcess(data, "id", "outcome", "treatment",
                                              time_invariant_vars = NULL,
                                              number_lags = 1,
                                              periods_post = 1)[[2]]), 2)})

# Test that function returns dataframes with correct data types
test_that("GDiD.DataProcess returns dataframes with correct data types",
          {expect_is(GDiD.DataProcess(data, "id", "outcome", "treatment",
                                      time_invariant_vars = NULL,
                                      number_lags = 1,
                                      periods_post = 1)[[1]]$id, "numeric")
           expect_is(GDiD.DataProcess(data, "id", "outcome", "treatment",
                                      time_invariant_vars = NULL,
                                      number_lags = 1,
                                      periods_post = 1)[[1]]$cluster_id, "numeric")
           expect_is(GDiD.DataProcess(data, "id", "outcome", "treatment",
                                      time_invariant_vars = NULL,
                                      number_lags = 1,
                                      periods_post = 1)[[1]]$Y, "numeric")
           expect_is(GDiD.DataProcess(data, "id", "outcome", "treatment",
                                      time_invariant_vars = NULL,
                                      number_lags = 1,
                                      periods_post = 1)[[1]]$A, "numeric")
           })


# Test that GDID.DataProcess returns correct lags
test_that("GDID.DataProcess returns correct lags for first post period",
          {pd3 = GDiD.DataProcess(data3,
                                  "id",
                                  "outcome",
                                  "treatment",
                                  time_invariant_vars = NULL,
                                  time_varying_vars = "outcome",
                                  number_lags = 2,
                                  periods_post = 1) %>%
            map(~select(.x, id, outcome_lag1, outcome_lag2)) %>%
            invoke(rbind, .) %>%
            arrange(id)
           expect_equal(pd3$outcome_lag1, c(2, 1, 7, 6))
           expect_equal(pd3$outcome_lag2, c(1, 0, 6, 5))})

test_that("GDID.DataProcess returns correct lags for second post period",
          {d3b = GDiD.DataProcess(data3,
                                  "id",
                                  "outcome",
                                  "treatment",
                                  time_invariant_vars = NULL,
                                  time_varying_vars = "outcome",
                                  number_lags = 1,
                                  periods_post = 2) %>%
            map(~select(.x, id, outcome_lag2)) %>%
            invoke(rbind, .) %>%
            arrange(id)
          expect_equal(pd3b$outcome_lag2, c(2, 0, 7, 5))})

test_that("GDID.DataProcess returns equal number of lags for each
          time-varying covariate",
          {pd = GDiD.DataProcess(data3,
                                 "id",
                                 "outcome",
                                 "treatment",
                                 time_invariant_vars = "X",
                                 time_varying_vars = c("outcome", "Xt"),
                                 number_lags = 1,
                                 periods_post = 1)

          expect_equal(c(2, 2), map_dbl(pd, ~length(grep("lag", names(.x)))))
          })

test_that("GDID.DataProcess returns same value of time-invariant covariate in
          each dataframe",
          {pd = GDiD.DataProcess(data3,
                                  "id",
                                  "outcome",
                                  "treatment",
                                  time_invariant_vars = "X",
                                  time_varying_vars = c("outcome", "Xt"),
                                  number_lags = 1,
                                  periods_post = 1)
            expect_equal(pd[[1]]$X, pd[[2]]$X)})

test_that("GDiD.StaggeredDataProcess returns the correct number of dataframes",
          {d4 = map(c(1, 2, 5), ~GDiD.StaggeredDataProcess(data4,
                                          id_var = "id",
                                          cluster_var = "cluster",
                                          outcome_var = "outcome",
                                          treatment_var = "treatment",
                                          time_varying_vars = "outcome",
                                          time_invariant_vars = NULL,
                                          number_bias_estimates = 1,
                                          number_lags = .x))
          expect_equal(length(d4[[1]]), 2)
          expect_equal(length(d4[[1]][[1]]), 2)
          expect_equal(length(d4[[1]][[2]]), 1)
          expect_equal(length(d4[[2]][[1]]), 1)
          expect_equal(length(d4[[2]][[2]]), 1)
          expect_type(d4[[3]][[1]], "character")
          expect_type(d4[[3]][[2]], "character")
          })


