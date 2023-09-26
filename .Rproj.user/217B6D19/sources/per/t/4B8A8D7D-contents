setwd("~/GDiD/gdid/R")
source("gdid.R")
library(dtplyr)
library(tidyverse)

# Functions to simulate data and check CI coverage -----------------------------------
GenerateFakeData <- function(number_units,
                             number_time,
                             time_treated,
                             groups = FALSE,
                             units_per_group = NULL) {

  if (groups == FALSE) {
    units_per_group = 1
    data = tibble(
      id = rep(1:number_units, each = number_time),
      cluster_id = rep(1:(number_units / units_per_group), each = number_time * units_per_group),
      time = rep(c(1:number_time), number_units),
    ) %>%
      group_by(id) %>%
      mutate(theta = rnorm(1, 0, 1)) %>%
      group_by(cluster_id) %>%
      mutate(U = rbinom(1, 1, 0.5),
             A = ifelse(U == 1, 1, 0)) %>%
      ungroup() %>%
      mutate(X = rnorm(number_units*number_time, 0, 1),
             mu0 = -U + X + theta, #+ rnorm(total_number_units*number_time, 0, 1),
             mu1 = 2 * (X + theta + U), #+ rnorm(total_number_units*number_time, 0, 1),
             post = ifelse(time >= time_treated, 1, 0),
             Y = case_when(
               post == 1 ~ A * mu1 + (1 - A) * mu0,
               post == 0 ~ mu0,
             ))
  }

  if (groups == TRUE) {

    units_per_group = floor(rexp(number_units, 1 / units_per_group)) + 1
    total_number_units = sum(units_per_group)

    data = tibble(
      id = rep(1:total_number_units, each = number_time),
      cluster_id = unlist(map2(1:number_units, units_per_group, ~rep(.x, .y * number_time))),
      time = rep(c(1:number_time), total_number_units)
    ) %>%
      group_by(id) %>%
      mutate(theta = rnorm(1, 0, 1)) %>%
      group_by(cluster_id) %>%
      mutate(U = rbinom(1, 1, 0.5),
             A = ifelse(U == 1, 1, 0)) %>%
      ungroup() %>%
      mutate(X = rnorm(total_number_units*number_time, 0, 1),
             mu0 = -U + X + theta, #+ rnorm(total_number_units*number_time, 0, 1),
             mu1 = 2 * (X + theta + U), #+ rnorm(total_number_units*number_time, 0, 1),
             post = ifelse(time >= time_treated, 1, 0),
             Y = case_when(
               post == 1 ~ A * mu1 + (1 - A) * mu0 + rnorm(total_number_units*number_time, 0, 1),
               post == 0 ~ mu0 + rnorm(total_number_units*number_time, 0, 1),
             ))
  }
  data
}

CoverageCheck <- function(finite_sample = TRUE,
                          weights = NULL,
                          N = 1000,
                          number_bias_estimates,
                          units_per_group,
                          seed = 10,
                          datasets = NULL,
                          aux_data = NULL,
                          sl_lib = c("SL.glm", "SL.ranger")) {
  if (finite_sample == TRUE) {
    A = datasets[[1]] %>%
      nest(data = c(-cluster_id, -A)) %>%
      .$A

    Anew = sample(A, size = length(A), replace = FALSE)

    datasets = map(datasets, ~nest(.x, data = c(-cluster_id, -A)) %>%
                     mutate(A = Anew) %>% unnest(cols = c(data))) %>%
      map2(aux_data, ~left_join(.x, .y, by = c("id", "cluster_id"))) %>%
      map2(c(1, rep(0, 2)), ~mutate(.x, post = .y)) %>%
      map(~mutate(.x,
                  Y = case_when(
                    post == 1 ~ A * Y1 + (1 - A) * Y0,
                    post == 0 ~ Y0)))

    estimands = unlist(map(datasets, ~with(.x, mean(Y1 - Y0))))

    estimands = c(estimands[1],
                  estimands[1],
                  0, 0, 0)

    datasets = map(datasets, ~select(.x, -Y1, -Y0, -post))

  }

  if (finite_sample == FALSE) {

    data = GenerateFakeData(number_units = N,
                            number_time = 10,
                            time_treated = 7,
                            units_per_group)

    datasets = GDiD.DataProcess(data, "id", "cluster_id", "X", "theta",
                                number_lags = 2,
                                periods_post = 1,
                                number_bias_estimates = number_bias_estimates)

    estimands = rep(0, 5)
  }

  res = GenerateEIFs(datasets, sl_lib = sl_lib)

  res.sp = OutputEstimates(res, variance_estimate = "SP",
                  cluster_id = "cluster_id",
                  bias_weights = rep(1 / number_bias_estimates, number_bias_estimates)) %>%
    mutate(cis = "SP")

  res.fs = OutputEstimates(res, variance_estimate = "FS",
                           cluster_id = "cluster_id",
                           bias_weights = rep(1 / number_bias_estimates, number_bias_estimates)) %>%
    mutate(cis = "FS")

  bind_rows(res.sp, res.fs) %>%
    mutate(estimand = rep(estimands, 2))
}

# Check a single dataset -------------------------------------------------------------
data = GenerateFakeData(number_units = 1000,
                        number_time = 6,
                        time_treated = 3,
                        groups = TRUE,
                        units_per_group = 5)

staggered_data = data %>%
  group_by(cluster_id) %>%
  mutate(cohort = rbinom(1, 1, 0.5)) %>%
  ungroup() %>%
  mutate(post = case_when(
    cohort == 1 & time >= 3 ~ 1,
    cohort == 1 & time < 3  ~ 0,
    cohort == 0 & time >= 4 ~ 1,
    cohort == 0 & time < 4  ~ 0)) %>%
  mutate(
    mu1 = case_when(
      cohort == 0 & time == 4 ~ 1 + mu1,
      cohort == 0 & time == 5 ~ 2 + mu1,
      cohort == 0 & time == 6 ~ 3 + mu1,
      cohort == 0 ~ mu1,
      TRUE ~ mu1
    )
  ) %>%
  mutate(
    Y = case_when(
      A == 1 & post == 1 ~ mu1,
      A == 0 | post == 0 ~ mu0)
  ) %>%
  mutate(Y = Y + rnorm(nrow(.), 0, 1))

test1 = StaggeredAdoptionEffects(staggered_data %>% select(-cohort),
                                 id_var = "id",
                                 outcome_var = "Y",
                                 treatment_var = "A",
                                 time_varying_vars = "Y",
                                 time_invariant_vars = "theta",
                                 number_lags = 1,
                                 sl_lib = c("SL.glm"),
                                 number_bias_estimates = 1,
                                 outer_folds = 1)

StaggeredSim <- function(groups, units_per_group, number_units = 1000) {

  data = GenerateFakeData(number_units = number_units,
                          number_time = 6,
                          time_treated = 5,
                          groups = groups,
                          units_per_group = units_per_group)

  staggered_data = data %>%
    group_by(cluster_id) %>%
    mutate(cohort = rbinom(1, 1, 0.5)) %>%
    ungroup() %>%
    mutate(post = case_when(
      cohort == 1 & time >= 3 ~ 1,
      cohort == 1 & time < 3  ~ 0,
      cohort == 0 & time >= 4 ~ 1,
      cohort == 0 & time < 4  ~ 0)) %>%
    mutate(
      mu1 = case_when(
        cohort == 0 & time == 4 ~ 1 + mu1,
        cohort == 0 & time == 5 ~ 2 + mu1,
        cohort == 0 & time == 6 ~ 3 + mu1,
        cohort == 0 ~ mu1,
        TRUE ~ mu1
      )
    ) %>%
    mutate(
      Y = case_when(
        A == 1 & post == 1 ~ mu1,
        A == 0 | post == 0 ~ mu0)
    ) %>%
    mutate(Y = Y + rnorm(nrow(.), 0, 1)) %>%
    select(-cohort)

  test1 = StaggeredAdoptionEffects(staggered_data,
                                   id_var = "id",
                                   outcome_var = "Y",
                                   treatment_var = "A",
                                   time_varying_vars = NULL,
                                   time_invariant_vars = "theta",
                                   number_lags = 1,
                                   number_bias_estimates = 2,
                                   sl_folds = 2,
                                   outer_folds = 1,
                                   sl_lib = "SL.glm")

  list(test1 = test1)
}

# -- run simulations
test.ss1 = map(1:250, ~StaggeredSim(groups = TRUE, number_units = 1000, units_per_group = 5))
test.ss2 = map(1:250, ~StaggeredSim(groups = FALSE, number_units = 1000))

saveRDS(list(test.ss1, test.ss2), "R/unit-tests/all-sims.rds")

# ---
sim_res <- readRDS("R/unit-tests/all-sims.rds")

AnalyzeSims <- function(sim_res_list) {
  map(sim_res_list, ~.x$test1) %>%
    invoke(rbind, .) %>%
    group_by(estimator, periods_post, cohort) %>%
    mutate(truth = case_when(
      grepl("Bias", estimator) ~ -1,
      grepl("^Ign", estimator) & cohort  == 1 ~ 2,
      grepl("^gIgn", estimator) & cohort == 1 ~ 3,
      grepl("^Ign", estimator) & cohort == 1 & periods_post == "average" ~ 2,
      grepl("gIgn", estimator) & cohort == 1 & periods_post == "average" ~ 3,
      grepl("^Ign", estimator)  & cohort == 2 & periods_post == 1 ~ 3,
      grepl("^Ign", estimator)  & cohort == 2 & periods_post == 2 ~ 4,
      grepl("^Ign", estimator)  & cohort == 2 & periods_post == 3 ~ 5,
      grepl("^gIgn", estimator) & cohort == 2 & periods_post == 1 ~ 4,
      grepl("^gIgn", estimator) & cohort == 2 & periods_post == 2 ~ 5,
      grepl("^gIgn", estimator) & cohort == 2 & periods_post == 3 ~ 6,
      grepl("^Ign", estimator)  & cohort == 2 & periods_post == "average" ~ 4,
      grepl("gIgn", estimator) & cohort == 2 & periods_post == "average" ~ 5,
      grepl("^gIgn", estimator) & cohort == "average" & periods_post == "average" ~ (4/7)*3 + (3/7)*5,
      grepl("^Ign", estimator) & cohort == "average" & periods_post == "average" ~ (4/7)*2 + (3/7)*4)) %>%
    reframe(covered = mean(lci < truth & uci > truth),
            estimate.avg = mean(estimate),
            cohort_weight = mean(cohort_weight),
            var.ratio = mean(varest) / var(estimate)) %>%
    arrange(estimator, periods_post, cohort)
}

print.data.frame(AnalyzeSims(sim_res[[2]])) %>%
  arrange(cohort, periods_post) %>%
  arrange(covered) %>%
  .$covered %>%
  summary()

