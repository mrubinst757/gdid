library(tidyverse)
library(assertthat)
library(tidyverse)
library(SuperLearner)

################################################################################
################## SIMULTANEOUS ADOPTION PERIOD FUNCTIONS ######################
################################################################################

#' Take a long-dataset and convert it to a list of wide datasets by
#' pre and post treatment periods and for each time period. Input dataset
#' must contain a time variable called "time" and a post-treatment indicator
#' called "post"
#'
#' @param data dataframe in long format. should contain a time identifier,
#' a unit identifier, a post-treatment indicator, treatment, covariates
#'
#' @param id_var a unique identifying variable for each unit
#'
#' @param outcome_var a string specifying the outcome variable name
#'
#' @param treatment_var a string specifying an indicator for whether the
#' unit was *ever* treated. dataframe must include some units that were
#' never treated
#'
#' @param time_varying_vars character vector of time-varying covariates (including
#' pre-treatment outcomes)
#'
#' @param time_invariant_vars character vector of time-invariant covariates
#'
#' @param number_lags number of lags to include for each time-varying covariate
#'
#' @param periods_post how many periods since treatment began are we estimating
#' the treatment effect for (starting with 1)
#'
#' @param number_bias_estimates how many estimates of the bias should be calculated
#' using pre-period data
#'
#' @return a list of pre-treatment and post-treatment dataframes

GDiD.DataProcess <- function(data,
                             id_var,
                             outcome_var,
                             treatment_var,
                             cluster_var = NULL,
                             time_varying_vars = NULL,
                             time_invariant_vars,
                             number_lags,
                             periods_post,
                             number_bias_estimates = 1) {

  id_counts = as.numeric(table(data[[id_var]]))

  missing_counts = apply(data, 2, function(x) sum(is.na(x)))

  agrp = sort(as.numeric(unique(tapply(data$A, data$id, mean))))

  assert_that("time" %in% names(data), msg = "Variable called 'time' must be
              included in data")

  assert_that("post" %in% names(data), msg = "Variable called 'post' must be
              included in data")

  assert_that(identical(c(0, 1), as.numeric(sort(unique(data[["post"]])))), msg = "Variable 'post' must
              contain values of 0 and 1")

  assert_that(identical(c(0, 1), as.numeric(sort(unique(data[[treatment_var]])))),
              msg = "Ever treated indicator must contain values of 0 and 1")

  assert_that(identical(c(0, 1), agrp),
              msg = "Treatment value must be the same for all unique IDs")

  assert_that(!is.na(id_var), msg = "Dataframe must contain a unique identifier")

  assert_that(length(unique(id_counts)) == 1, msg = "Dataframe must be a balanced panel")

  assert_that(all(missing_counts == 0), msg = "Dataframe must not contain missing data")

  assert_that(is.numeric(data[[id_var]]), msg = "Unit ID variable must be numeric")

  if (!is.null(cluster_var)) {
    assert_that(is.numeric(data[[cluster_var]]), msg = "Cluster ID variable must be numeric")
  }

  if (is.null(cluster_var)) {
    data[["cluster_id"]] = data[[id_var]]
    cluster_var = "cluster_id"
  }

  # Identify treatment implementation time
  T_0 = max(data[["time"]][data[["post"]] == 0])


  if (!is.null(time_varying_vars)) {
    # Create lagged variables
    for (i in 1:(number_lags + periods_post + number_bias_estimates)) {
      for (var in time_varying_vars) {
        data = data %>%
          group_by(!!sym(id_var), !!sym(cluster_var)) %>%
          mutate(!!sym(paste0(var, "_lag", i)) := lag(!!sym(var), i)) %>%
          ungroup()
      }
    }
  }

  # Filter data for time T_0 + m and T_0 - k to T_0
  data_t0_m <- data[data$time == T_0 + periods_post, ]


  # Extract relevant columns for the first dataframe
  lag_start = 1 + (periods_post - 1)
  lag_stop  = number_lags + (periods_post - 1)

  if (!is.null(time_varying_vars)) {
    df1_vars <- c(id_var, cluster_var, time_invariant_vars,
                  paste0(rep(time_varying_vars, each = number_lags), "_lag",
                         rep(lag_start:lag_stop, times = length(time_varying_vars))),
                  outcome_var, treatment_var)
  }

  if (is.null(time_varying_vars)) {
    df1_vars <- c(id_var, cluster_var, time_invariant_vars, outcome_var, treatment_var)
  }

  df1 <- data_t0_m[, df1_vars]

  # Extract relevant columns for the second dataframe

  df_pre_list <- list()

  for (i in 1:number_bias_estimates) {

    lag_start = 1 + (periods_post - 1)
    lag_stop  = number_lags + (periods_post - 1)

    data_t0 <- data[data$time == T_0 - (i - 1), ]

    if (!is.null(time_varying_vars)) {
      df_vars <- c(id_var, cluster_var, time_invariant_vars,
                   paste0(rep(time_varying_vars, each = number_lags),
                          "_lag", rep(lag_start:lag_stop, times = length(time_varying_vars))),
                   outcome_var, treatment_var)
    }
    if (is.null(time_varying_vars)) {
      df_vars <- c(id_var, cluster_var, time_invariant_vars, outcome_var, treatment_var)
    }

    df_pre_list[[i]] = data_t0[, df_vars]
  }

  append(list(df1), df_pre_list) %>%
    map(~rename(.x,
                id = !!sym(id_var),
                cluster_id = !!sym(cluster_var),
                Y = !!sym(outcome_var),
                A = !!sym(treatment_var)))
}

#' Create influence function values for the ATT using an input dataset with
#' estimated nuisance parameters
#'
#' @param data dataset containing estimates of the propensity score (pi),
#' the regression models (mu_1 and mu_0)
#'
#' @return dataframes with EIFs for the ATT and vectors to obtain
#' variance estimates for both finite-sample and super-population inferences
ConstructEIF.ATT <- function(data) {

  # EIF estimates for calculating the mean
  p = mean(data$A)

  data$eif.y1 = with(data, A * Y) / p

  data$eif.y0 = with(data, ((1 - A) * pi * Y + (A - pi) * mu_0) / (1 - pi)) / p

  data$eif.att = with(data, eif.y1 - eif.y0)

  data
}


#' Estimate nuisance functions for double robust estimator of ATT.
#'
#' @param outcome_tx_dat Dataset with outcome and treatment variables.
#'
#' @param X_dat Covariate data.
#'
#' @param outer_folds Number of folds for cross-validation and splitting
#'
#' @param sl_lib List of SuperLearner libraries for SuperLearner
#'
#' @param sl_folds Number of folds for the SuperLearner
#'
#' @return Dataset with estimated nuisance functions
EstimateNuisanceFunctions <- function(outcome_tx_dat,
                                      X_dat,
                                      outer_folds,
                                      sl_lib,
                                      sl_folds,
                                      folds_seed = NULL) {

  #######################
  ### Data Cleaning
  #######################

  ### Sort data
  X_dat = arrange(X_dat, cluster_id, id)
  outcome_tx_dat = arrange(outcome_tx_dat, cluster_id, id)
  clusters = unique(X_dat$cluster_id)
  number_per_cluster = as.numeric(table(X_dat$cluster_id))

  ### Create folds
  set.seed(folds_seed)
  folds_agg <- sample(outer_folds, length(clusters), replace = TRUE)
  folds <- unlist(map2(folds_agg, number_per_cluster, ~rep(.x, .y)))
  set.seed(NULL)

  ### Turn covariate data into model matrix
  X_dat <- X_dat %>%
    stats::model.matrix(~ . -1, data = .) %>%
    as.data.frame(.)

  ### Column naming hack to make models work with SuperLearner
  colnames(X_dat) <- paste0("V", seq(1, ncol(X_dat)))

  ### Initialize nuisance function columns

  dat <- dplyr::select(outcome_tx_dat, id)
  dat$mu_1 <- NA; dat$mu_0 <- NA; dat$pi <- NA;

  a_vals <- sort(unique(outcome_tx_dat$A))

  ###################################
  ### Nuisance Function Estimation
  ###################################

  for (FOLD in seq(1, outer_folds, 1)) {

    #cat("\n------------------------------------------------",
    #    "\n-- Estimating fold ", FOLD, " out of ", outer_folds,
    #    "\n------------------------------------------------")

    TEST_ROWS <- folds == FOLD
    TRAIN_ROWS <- !TEST_ROWS

    if (outer_folds == 1) {
      TRAIN_ROWS <- TEST_ROWS
    }

    ########################
    ### Pi: P(A = 1 | X) ###
    ########################

    #cat("\nEstimating Pi\n")

    outcome <- outcome_tx_dat$A[TRAIN_ROWS]

    Xdat.train = X_dat[TRAIN_ROWS, -c(1:2)]
    Xdat.test = X_dat[TEST_ROWS, -c(1:2)]

    if (is.numeric(Xdat.train)) {
      Xdat.train = tibble(X = Xdat.train)
    }

    if (is.numeric(Xdat.test)) {
      Xdat.test = tibble(X = Xdat.test)
    }

    mod_pi <-
      SuperLearner::SuperLearner(
        Y = outcome,   # Outcome
        X = Xdat.train, # X value training data (remove id)
        SL.library = sl_lib,                 # Algorithms to use
        family = binomial(),                   # Binomial outcome
        cvControl = SuperLearner.CV.control(V = sl_folds, stratifyCV = TRUE)
      )

    dat[["pi"]][TEST_ROWS] <- predict(mod_pi, newdata = Xdat.test)$pred[,1]

    ############################
    ### Mu: P(Y | X, A = a) ####
    ############################

    ### Indicator for each treatment value
    for (a in a_vals) {
      outcome <- outcome_tx_dat$Y[outcome_tx_dat$A == a & TRAIN_ROWS]

      outcome_type = if_else(length(unique(outcome)) == 2, "Binary", "Continuous")

      family_specification <- switch(outcome_type,
                                     "Binary" = binomial(),
                                     "Continuous" = gaussian())

      Xdat.train = X_dat[outcome_tx_dat$A == a & TRAIN_ROWS, -c(1:2)]
      Xdat.test  = X_dat[TEST_ROWS, -c(1:2)]

      if (is.numeric(Xdat.train)) {
        Xdat.train = tibble(X = Xdat.train)
      }

      if (is.numeric(Xdat.test)) {
        Xdat.test = tibble(X = Xdat.test)
      }

      mod_mu <-
        SuperLearner::SuperLearner(
          Y = outcome,
          X = Xdat.train,
          SL.library = sl_lib,
          family = family_specification,
          cvControl = SuperLearner.CV.control(V = sl_folds)
        )
      dat[[paste0("mu_", a)]][TEST_ROWS] <- predict(mod_mu, newdata = Xdat.test)$pred[,1]
    }
  }
  outcome_tx_dat <- inner_join(outcome_tx_dat, dat, by = "id")
  assert_that(nrow(outcome_tx_dat) == nrow(X_dat), msg = "We lost people!?")

  outcome_tx_dat
}

#' Estimate influence functions on list of datasets to estimate both the ignorability
#' estimators and to obtain estimates of the bias
#'
#' @param data_list list of datasets output from the GDiD.DataProcess function
#'
#' @param sl_lib SuperLearner libraries; default is SL.glm and SL.ranger
#'
#' @param outer_folds the number of outer folds for cross-fitting; default is 2
#'
#' @param sl_folds the number of cross-validation folds for the SuperLearner ensemble
#' weights; the default is 2
#'
#' @return list of dataframes with estimated influence functions for the average treatment
#' effect on the treated

GenerateEIFs <- function(data_list, sl_lib = c("SL.glm", "SL.ranger"), outer_folds = 2, sl_folds = 2,
                         folds_seed = NULL) {

  if (is.null(folds_seed)) {
    folds_seed = runif(1, 1, 2e6) # necessary to keep the folds constant for each time period
  }

  # separate dataframes into covariate and outcome information
  Ydat = map(data_list, ~.x %>% select(c("Y", "A", "id", "cluster_id")))

  Xdat = map(data_list, ~.x %>% select(-c("Y", "A")) %>%
               select(cluster_id, id, everything()))

  res  = map2(Ydat, Xdat,
              ~EstimateNuisanceFunctions(.x, .y, outer_folds, sl_lib, sl_folds, folds_seed) %>%
                ConstructEIF.ATT())
  res
}

#' Return a dataframe of treatment effect estimates using the output of ConstructEIF.ATT
#'
#' @param eif_data dataframe produced by ConstructEIF.ATT
#'
#' @param variance_estimate which type of variance estimator to use. SP indicates that the variance
#' estimate is with respect to a superpopulation parameter while FS is a conservative estimate for
#' the finite-sample variance with respect only to the randomization of treatment assignment. Note that
#' the SP variance estimate is valid for the finite sample estimand, but may be unnecessarily conservative.
#' The FS variance estimate is also in general conservative.
#'
#' @param cluster_id a vector of cluster identifiers; will affect variance calculation but not the
#' point estimate calculation
#'
#' @param bias_weights a vector of weights for how to average the pre-treatment bias estimates. if
#' NULL then the weights are set to be equal
#'
#' @param alpha Type I Error rate for associated hypothesis test
#'
#' @return a dataframe containing estimates of the debiased estimator, the
#' ignorability estimator, the bias, and the individual bias estimates along
#' with corresponding variance estimates and confidence intervals

OutputEstimates <- function(eif_data,
                            cluster_id = NULL,
                            bias_weights = NULL,
                            alpha = 0.05) {

  if (!is.null(bias_weights)) {
    assert_that(length(bias_weights) == (length(eif_data) - 1),
                msg = "Bias weights must have same number of
                elements as pre-treatment bias estimates")
  }

  # If no cluster ID given, assign each unit to it's own cluster
  if (is.null(cluster_id)) {
    cluster_id = 1:nrow(eif_data[[1]])
  }

  # If no bias weights are given, set to equal
  if (is.null(bias_weights)) {
    bias_weights = rep(1 / (length(eif_data) - 1), length(eif_data) - 1)
  }

  dat1 = eif_data[[1]]
  dat0 = eif_data[-1]

  # Point estimates -----------------------------------------------

  # Calculate Ignorability Estimate
  base_estimate = dat1 %>%
    ungroup() %>%
    summarize(est = mean(eif.att)) %>%
    .$est

  # Calculate Bias Estimate
  bias_estimates = map(dat0, ~ .x %>%
                         ungroup() %>%
                         summarize(est = mean(eif.att)) %>%
                         .$est)

  bias_estimate = map2(bias_estimates, bias_weights, ~.x * .y) %>%
    reduce(`+`)

  # Calculate debiased estimate
  final_estimate = base_estimate - bias_estimate

  # Variance estimates -------------------------------------------
  base_var_all  = EIFVarianceEstimate(eif_data[[1]],
                                      bias_estimate = FALSE)

  bias_var_all  = EIFVarianceEstimate(eif_data[-1],
                                      bias_weights,
                                      bias_estimate = TRUE)

  final_var_all = EIFVarianceEstimate(eif_data,
                                      bias_weights,
                                      bias_estimate = FALSE)

  base_var_eif  = base_var_all[[2]]
  final_var_eif = final_var_all[[2]]
  bias_var_eif  = bias_var_all[[2]]

  base_var  = base_var_all[[1]]
  bias_var  = bias_var_all[[1]]
  final_var = final_var_all[[1]]

  bias_vars = map_dbl(eif_data[-1], ~EIFVarianceEstimate(.x)[[1]])

  threshold = qnorm(1 - alpha / 2)

  res = tibble(
    estimator = c("gIgn", "Ign", "Bias", sprintf("Bias-t%s", 1:length(dat0))),
    estimate  = c(final_estimate, base_estimate, bias_estimate, unlist(bias_estimates)),
    weights   = c(1, 1/2, 1/2, 0.5 * bias_weights),
    varest = c(final_var, base_var, bias_var, bias_vars),
    lci = estimate - threshold * sqrt(varest),
    uci = estimate + threshold * sqrt(varest)
  )

  list(res = res, eif_data = list(base_var = base_var_eif,
                                  gign_var = final_var_eif,
                                  bias_var = bias_var_eif))
}

#' Calculate variance estimates for various functionals
#'
#' @param eif_data dataframe or list of dataframes containing output from
#' ConstructEIF.ATT (invoked in the GenerateEIFs)
#'
#' @param bias_weights vector of weights for each bias estimate
#'
#' @param variance_estimate "SP" for super-population; "FS" for finite-sample
#'
#' @param bias_estimate TRUE if the estimate is for an estimate of the bias
#' FALSE otherwise
#'
#' @return a one number estimate of the input variance
EIFVarianceEstimate <- function(eif_data,
                                bias_weights = NULL,
                                bias_estimate = FALSE) {

  if (!is.data.frame(eif_data)) {

    col_match = if_else(bias_estimate, "V[1-9][0-9]?", "V[2-9][0-9]?")

    cluster_id = eif_data[[1]]$cluster_id

    id = eif_data[[1]]$id

    variance_estimate = map(eif_data, ~select(.x, eif.att)) %>%
      invoke(cbind, .) %>%
      set_names(sprintf("V%s", 1:ncol(.))) %>%
      cbind(tibble(cluster_id = cluster_id, id = id)) %>%
      rowwise() %>%
      mutate(eif.b = sum(c_across(matches(col_match)) * bias_weights)) %>%
      ungroup()

    sample_size = length(id)

    if (bias_estimate == FALSE) {

      variance_estimate_data = variance_estimate %>%
        mutate(eif_final = V1 - eif.b) %>%
        select(eif_final, id, cluster_id)

      tau = mean(variance_estimate_data$eif_final)

      variance_estimate = variance_estimate_data %>%
        group_by(cluster_id) %>%
        summarize(eif = sum(eif_final - tau)) %>%
        summarize(variance_estimate = sum(eif^2) / sample_size^2) %>%
        .$variance_estimate
    }

    if (bias_estimate == TRUE) {

      tau = mean(variance_estimate$eif.b)

      variance_estimate_data = variance_estimate

      variance_estimate = variance_estimate_data %>%
        group_by(cluster_id) %>%
        summarize(eif = sum(eif.b - tau)) %>%
        summarize(variance_estimate = sum(eif^2) / sample_size^2) %>%
        .$variance_estimate
    }
  }

  if (is.data.frame(eif_data)) {

    sample_size = length(eif_data$id)

    variance_estimate_data = eif_data %>%
      select(eif.att, id, cluster_id)

    tau = mean(variance_estimate_data$eif.att)

    variance_estimate = variance_estimate_data %>%
      group_by(cluster_id) %>%
      summarize(eif = sum(eif.att - tau)) %>%
      summarize(variance_estimate = sum(eif^2) / sample_size^2) %>%
      .$variance_estimate
  }

  list(variance_estimate = variance_estimate,
       variance_estimate_data = variance_estimate_data)
}

################################################################################
####################### STAGGERED ADOPTION  FUNCTIONS ##########################
################################################################################

#' StaggeredAdoptionEffects: returns list of cohort by time, cohort, and overall
#' effect estimates by time. right now the overall estimates by time cannot account
#' for clustering; however, the cohort-averages and cohort/time averages will
#' account for clustering if present. to generate cohort-specific estimates, all
#' post-treatment time-periods are equally weighted (within cohort) by default
#'
#' @param data dataset containing variables including "time", "post", and a treatment
#' indicator, optional covariates, and an outcome
#' @param id_var unique id
#' @param outcome_var outcome variable
#' @param treatment_var treatment variable
#' @param time_varying_vars time-varying covariates
#' @param time_invariant_vars time-invariant covariates
#' @param number_lags number of pre-treatment outcomes / time-varying variables to control for
#' @param sl_lib SuperLearner libraries
#' @param number_bias_estimates number of bias estimates to calculate in pre-treatment period
#' @param alpha type I error rate
#' @param cohort_weights weights to put on each cohort. defaults to being proportional to the
#' number of post-treatment time periods given by each cohort
#' @param outer_folds number of splits for cross-fitting
#' @param sl_folds number of folds to learn SuperLearner weights
#'
#' @return dataframe of effect estimates and uncertainty estimates
StaggeredAdoptionEffects <- function(data,
                                     id_var,
                                     cluster_var = NULL,
                                     outcome_var,
                                     treatment_var,
                                     time_varying_vars,
                                     time_invariant_vars,
                                     number_lags,
                                     sl_lib = "SL.glm",
                                     number_bias_estimates = 1,
                                     alpha = 0.05,
                                     cohort_weights = NULL,
                                     outer_folds = 2,
                                     sl_folds = 2) {

  time_weights = NULL

  if (!is.null(cohort_weights)) {
    assert_that(all(cohort_weights) > 0)
    cohort_weights = cohort_weights / sum(cohort_weights)
  }

  procdat = GDiD.StaggeredDataProcess(data,
                                      id_var = id_var,
                                      cluster_var = cluster_var,
                                      outcome_var = outcome_var,
                                      treatment_var = treatment_var,
                                      time_varying_vars = time_varying_vars,
                                      time_invariant_vars = time_invariant_vars,
                                      number_lags = number_lags,
                                      number_bias_estimates = number_bias_estimates)

  res_list = list(); fin_list = list(); resdat = list(); res = list()

  cohort_ids = 1:length(procdat)
  dropped = map_lgl(procdat, is.character)
  cohort_ids = cohort_ids[dropped == FALSE]
  procdat = discard(procdat, is.character)

  assert_that(length(procdat) > 0, msg = "Not possible to compute estimate given
              input number of lags, bias estimates, and available data")

  if (is.null(cohort_weights)) {
    number_post_estimates = map_dbl(procdat, length)
    cohort_weights = number_post_estimates / sum(number_post_estimates)
  }

  for (i in 1:length(procdat)) {

    folds_seed = runif(1, 1, 2e6) # necessary to keep the folds constant for each time period

    res_list[[i]] = map(procdat[[i]],
                        ~GenerateEIFs(.x,
                                      sl_lib = sl_lib,
                                      outer_folds = outer_folds,
                                      sl_folds = sl_folds,
                                      folds_seed = folds_seed))

    fin_list = map(res_list[[i]],
                   ~OutputEstimates(.x,
                                    alpha = alpha,
                                    cluster_id = "cluster_id"))

    cohort_res = CohortAverageByTime(fin_list,
                                     alpha = alpha,
                                     time_weights = time_weights)

    resdat[[i]] = cohort_res$resdat

    res[[i]] = cohort_res$res %>%
      mutate(periods_post = "average") %>%
      bind_rows(
        map(fin_list, ~.x$res) %>%
          map2(1:length(fin_list),
               ~mutate(.x, periods_post = as.character(.y))) %>%
          invoke(rbind, .)
      ) %>%
      mutate(cohort = cohort_ids[i], cohort_weight = cohort_weights[i])
  }

  res = res %>%
    invoke(rbind, .) %>%
    select(-weights)

  if (length(procdat) > 1) {

    ign  = map(resdat, ~.x$res.ign)
    gign = map(resdat, ~.x$res.gign)
    bias = map(resdat, ~.x$res.bias)

    threshold = qnorm(1 - alpha / 2)

    res = res %>%
      bind_rows(
        res %>%
          filter(periods_post == "average") %>%
          group_by(estimator) %>%
          reframe(estimate = sum(estimate * cohort_weight)) %>%
          mutate(varest = c(CohortAverageVariance(bias,  cohort_weights),
                            CohortAverageVariance(ign, cohort_weights),
                            CohortAverageVariance(gign, cohort_weights)),
                 lci = estimate - threshold * sqrt(varest),
                 uci = estimate + threshold * sqrt(varest),
                 periods_post = "average")
      ) %>%
      mutate(cohort = replace_na(as.character(cohort), "average")) %>%
      arrange(estimator, cohort, periods_post)
  }
  return(res)
}

#' GDiD.StaggeredDataProcess
#'
#' Take a long-dataset and convert it to a list of wide datasets by
#' pre and post treatment periods and for each time period. Input dataset
#' must contain a time variable called "time" and a post-treatment indicator
#' called "post." Each element of the list corresponds to different cohorts
#' and within each cohort there is a list of datasets for each post-treatment
#' time-period
#'
#' @param data dataframe in long format. should contain a time identifier,
#' a unit identifier, a post-treatment indicator, treatment, covariates
#'
#' @param id_var a unique identifying variable for each unit
#'
#' @param cluster_var a unique identifying variable for each cluster. if NULL (
#' by default) is set equal to id_var
#'
#' @param outcome_var specifies the outcome of interest
#'
#' @param treatment_var specifies the treatment of interest
#'
#' @param time_varying_vars character vector of time-varying covariates (including
#' pre-treatment outcomes)
#'
#' @param time_invariant_vars character vector of time-invariant covariates
#'
#' @param number_lags number of lags to include for each time-varying covariate
#'
#' @param number_bias_estimates how many estimates of the bias should be calculated
#' using pre-period data
#'
#' @return a list of pre-treatment and post-treatment dataframes

GDiD.StaggeredDataProcess <- function(data,
                                      id_var,
                                      cluster_var = NULL,
                                      outcome_var,
                                      treatment_var,
                                      time_varying_vars,
                                      time_invariant_vars,
                                      number_lags,
                                      number_bias_estimates) {

  id_counts = as.numeric(table(data[[id_var]]))

  missing_counts = apply(data, 2, function(x) sum(is.na(x)))

  agrp = sort(as.numeric(unique(tapply(data$A, data$id, mean))))

  assert_that("time" %in% names(data), msg = "Variable called 'time' must be
              included in data")

  assert_that("post" %in% names(data), msg = "Variable called 'post' must be
              included in data")

  assert_that(!is.na(id_var), msg = "Dataframe must contain a unique identifier")

  assert_that(length(unique(id_counts)) == 1, msg = "Dataframe must be a balanced panel")

  assert_that(all(missing_counts == 0), msg = "Dataframe must not contain missing data")

  assert_that(is.numeric(data[[id_var]]), msg = "Unit ID variable must be numeric")

  assert_that(identical(c(0, 1), as.numeric(sort(unique(data[[treatment_var]])))),
              msg = "There must exist units that were never treated during the study period")

  if (!is.null(cluster_var)) {
    assert_that(is.numeric(data[[cluster_var]]), msg = "Cluster ID variable must be numeric")
  }

  if (is.null(cluster_var)) {
    data[["cluster_id"]] = data[[id_var]]
    cluster_var = "cluster_id"
  }

  print("Generating cohort indicators indexed from earliest to latest treated")

  data = GenerateCohortIndicator(data,
                                 treatment_var = treatment_var,
                                 cluster_var   = cluster_var)

  # extract cohort values for treated cohorts
  cohort_vals = sort(unique(data[["cohort"]]))[-1]

  # create seperate dataset for each cohort
  cohort_datalist = map(cohort_vals, ~filter(data, !!sym(treatment_var) == 0 | cohort == .x))

  data_list = rep(list(NULL), length(cohort_datalist))

  # apply GDiD.DataProcess to each cohort for each post-treatment time period
  for (i in 1:length(cohort_datalist)) {

    post_xwalk = cohort_datalist[[i]] %>%
      filter(cohort != 0) %>%
      distinct(post, time)

    data_input = cohort_datalist[[i]] %>%
      select(-post) %>%
      left_join(post_xwalk, by = "time")

    if (!is.null(time_varying_vars)) {

        stop = 1

        while (stop == 1) {

          T0 = max(data_input$time[data_input$post == 0])
          TF = max(data_input$time)
          K = TF - T0

          if (K > 0) {

            test = data_input %>%
              filter(time <= T0) %>%
              group_by(!!sym(id_var), !!sym(cluster_var)) %>%
              mutate(across(one_of(time_varying_vars), ~lag(., number_lags + (K - 1)))) %>%
              ungroup() %>%
              filter(time > T0 - number_bias_estimates) %>%
              select(one_of(time_varying_vars)) %>%
              filter_all(any_vars(is.na(.))) %>%
              nrow()

            if (test > 0)  {

              data_input = filter(data_input, time < TF)

              print(paste0("Effect for ", K, "th post-treatment time period for cohort ", i, " could not be calculated due to insufficient pre-period data. This cohort-x-time period was dropped from the analysis."))
            }
          }

          if (test == 0 | K == 0) {
            stop = 0
          }
        }

        if (K > 0) {

          post_values = data_input %>%
            filter(!!sym(treatment_var) != 0) %>%
            distinct(post, time) %>%
            filter(post == 1) %>%
            .$time

          for (j in 1:length(post_values)) {

            data_list[[i]][[j]] = GDiD.DataProcess(data = data_input,
                                                   id_var = id_var,
                                                   cluster_var = cluster_var,
                                                   outcome_var = outcome_var,
                                                   treatment_var = treatment_var,
                                                   time_varying_vars = time_varying_vars,
                                                   time_invariant_vars = time_invariant_vars,
                                                   number_lags = number_lags,
                                                   periods_post = j,
                                                   number_bias_estimates = number_bias_estimates)
          }
        }
        if (K == 0) {
          data_list[[i]] = paste0("The ", i, "-th cohort effect could not be calculated
                                  due to insufficient pre-period data given the input number
                                  of lags and time-varying variables.")
        }
      }

    if (is.null(time_varying_vars)) {

      post_values = data_input %>%
        filter(!!sym(treatment_var) != 0) %>%
        distinct(post, time) %>%
        filter(post == 1) %>%
        .$time

      for (j in 1:length(post_values)) {

        data_list[[i]][[j]] = GDiD.DataProcess(data = data_input,
                                               id_var = id_var,
                                               cluster_var = cluster_var,
                                               outcome_var = outcome_var,
                                               treatment_var = treatment_var,
                                               time_varying_vars = time_varying_vars,
                                               time_invariant_vars = time_invariant_vars,
                                               number_lags = number_lags,
                                               periods_post = j,
                                               number_bias_estimates = number_bias_estimates)
      }
    }
  }
  data_list
}

#' GenerateCohortIndicator: takes a dataframe with staggered adoption patterns
#' and creates indicators for each adoption-time cohort
#'
#' @param data a dataframe containing post-treatment indicator, time variable,
#' and cluster identifier
#'
#' @return dataframe with a cohort indicator; 0 for never-treated, the rest
#' are in order of treatment time
GenerateCohortIndicator <- function(data,
                                    cluster_var,
                                    treatment_var) {

  cohort_xwalk = data %>%
    group_by(time) %>%
    reframe(first_post = ifelse(sum(post) > 0, 1, 0)) %>%
    mutate(cohort = cumsum(ifelse(first_post >= 1, 1, 0)))

  res = data %>%
    group_by(!!sym(cluster_var)) %>%
    mutate(first_post = ifelse(cumsum(post) == 1, 1, 0)) %>%
    left_join(cohort_xwalk, by = c("time", "first_post")) %>%
    replace_na(list(cohort = 0)) %>%
    select(-first_post) %>%
    group_by(!!sym(cluster_var)) %>%
    mutate(cohort = ifelse(sum(!!sym(treatment_var)) > 0, max(cohort), 0)) %>%
    ungroup() %>%
    mutate(post = ifelse(!!sym(treatment_var) == 0, 0, post))

  res
}

#' CohortAverageByTime: returns the average effect estimate for a treatment
#' cohort across all post-treatment time periods
#'
#' @param input_list output from ConstructEIF.ATT
#' @param alpha Type I error rate
#' @param time_weights weights to average estimates across time. note: if they
#' do not sum to 1 they will be standardized to sum to 1
#'
#' @return list of effect estimates and confidence intervals averaged over time
CohortAverageByTime <- function(input_list,
                          alpha = 0.05,
                          time_weights = NULL) {

  if (is.null(time_weights)) {
    time_weights = rep(1 / length(input_list), length(input_list))
  }

  assert_that(all(time_weights > 0),
              msg = "Time weights must be greater than zero")

  time_weights = time_weights / sum(time_weights)

  eif_data  = map(input_list, ~.x$eif_data)
  ign_data  = map(eif_data, ~.x$base_var)
  gign_data = map(eif_data, ~.x$gign_var)
  bias_data = map(eif_data, ~.x$bias_var)

  res.ign  = CohortTimeVariance(ign_data,  time_weights, "Ign", alpha)
  res.gign = CohortTimeVariance(gign_data, time_weights, "gIgn", alpha)
  res.bias = CohortTimeVariance(bias_data, time_weights, "Bias", alpha)

  res = bind_rows(
    res.ign$res, res.gign$res, res.bias$res
  )

  resdat = list(res.ign  = res.ign$resdat,
                res.gign = res.gign$resdat,
                res.bias = res.bias$resdat)

  return(list(res = res, resdat = resdat))
}

#' CohortTimeVariance estimate the variance of the (weighted) average of all
#' post-treatment time periods for a single cohort
#'
#' @param eif_data output from ConstructEIF.ATT
#' @param time_weights weights for each post-treatment time-period in the average
#' @param estimator whether variance estimate is for the ignorability (Ign) or the
#'  Debiased-Ignorability estimator (gIgn)
#'
#' @return dataframe with the estimate and variance estimate associated
#' with the cohort average treatment effect estimate
CohortTimeVariance <- function(eif_data,
                               time_weights,
                               estimator,
                               alpha = 0.05) {

  eif_col = grep("eif", names(eif_data[[1]]), value = TRUE)

  final_data = map2(eif_data, time_weights,
                    ~mutate(.x, eif_final = !!sym(eif_col) * .y)) %>%
    reduce(left_join, by = c("id", "cluster_id")) %>%
    rowwise() %>%
    mutate(eif_final = sum(across(starts_with("eif_final")), na.rm = T)) %>%
    ungroup()

  estimate = mean(final_data$eif_final)

  sample_size = nrow(final_data)

  varest = final_data %>%
    group_by(cluster_id) %>%
    summarize(eif = sum(eif_final - estimate)) %>%
    summarize(variance_estimate = sum(eif^2) / sample_size^2) %>%
    .$variance_estimate

  threshold = qnorm(1 - alpha / 2)

  res = tibble(
    estimate  = estimate,
    estimator = estimator,
    varest = varest,
    lci = estimate - threshold * sqrt(varest),
    uci = estimate + threshold * sqrt(varest)
  )

  resdat = final_data %>%
    select(eif_final, id, cluster_id)

  list(res = res, resdat = resdat)
}

#' CohortAverageVariance: calculate the variance of the average of cohort-specific
#' estimates
#'
#' @param result_list list of results
#' @param cohort_weights weight to put on each cohort effect
#'
#' @return an estimate of the variance of the weighted average
CohortAverageVariance2 <- function(result_list,
                                  cohort_weights) {

  full_data = reduce(result_list, full_join, by = c("id", "cluster_id")) %>%
    mutate_at(vars(matches("eif")), ~replace_na(., 0))

  cohort_weights = (nrow(full_data) / map_dbl(result_list, nrow)) * cohort_weights

  sample_size = nrow(full_data)

  variance_estimate_data = full_data %>%
    rowwise() %>%
    mutate(eif = sum(c_across(matches("eif")) * cohort_weights))

  tau = mean(variance_estimate_data$eif)

  variance_estimate = variance_estimate_data %>%
    group_by(cluster_id) %>%
    summarize(eif = sum(eif - tau)) %>%
    summarize(varest = sum(eif^2) / sample_size^2) %>%
    .$varest

  variance_estimate
}

CohortAverageVariance <- function(result_list,
                                  cohort_weights) {

  cov_terms = expand.grid(
    indices1 = 1:length(result_list),
    indices2 = 1:length(result_list)
  ) %>%
    filter(indices1 > indices2)

  sample_sizes = map(result_list, ~nrow(.x))

  tau_list = map(result_list, ~mean(.x$eif_final))

  var_p1 = pmap(list(result_list, tau_list, sample_sizes),
                ~..1 %>%
                  group_by(cluster_id) %>%
                  summarize(eif = sum(eif_final - ..2)) %>%
                  summarize(var_est = sum(eif^2) / ..3^2)) %>%
    map2(cohort_weights, ~.x * .y^2) %>%
    reduce(`+`)

  var_p2 = map2(cov_terms$indices1, cov_terms$indices2,
                 ~CalculateCovariance(result_list[[.x]],
                                      result_list[[.y]],
                                      cohort_weights[c(.x, .y)])) %>%
    reduce(`+`)

  as.numeric(var_p1 + 2 * var_p2)
}

#' CalculateCovariance calculates the covariance between EIFs among shared units between two dataframes
#'
#' @param data1 dataframe with eif_final as a column
#' @param data2 dataframe with eif_final as a column
#' @param weights what weights are put on each column in the final average
#'
#' @return an estimate of two times the weighted covariance

CalculateCovariance <- function(data1, data2, weights) {
  all_data = full_join(data1, data2, by = c("id", "cluster_id")) %>%
    replace_na(list(eif_final.x = 0, eif_final.y = 0))

  sample_size = nrow(all_data)

  shared_data = all_data %>%
    group_by(cluster_id) %>%
    summarize(eif_final.x = sum(eif_final.x),
              eif_final.y = sum(eif_final.y)) %>%
    ungroup()

  number_groups = nrow(shared_data)

  cov.est = number_groups * (weights[1] * weights[2]) * cov(shared_data$eif_final.x, shared_data$eif_final.y) / sample_size^2

  cov.est
}

