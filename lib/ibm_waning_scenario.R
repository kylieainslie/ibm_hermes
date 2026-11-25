############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: generate example datasets for testing time-varying vaccine
#       effectiveness (VE) estimation methods. Unlike the lockdown/exit
#       scenario in ibm_interventions.R, this scenario:
#         - uses a slowed-down epidemic (lower transmission_prob) so the
#           outbreak unfolds over a long enough period for vaccination to
#           reach a meaningful share of the population before it burns out,
#         - rolls out vaccination gradually in many small increments rather
#           than a couple of large pulses, so vaccinated individuals are
#           staggered across many different calendar days, and
#         - uses a time-varying (ramp-then-wane) vaccine effectiveness, via
#           compute_time_varying_ve() in ibm_core.R, so there is a genuine
#           waning signal in the data for a VE(t) model to try to recover.
#
# Usage (from the root of the ibm_hermes repository):
#
#   source('lib/ibm_core.R')
#   source('lib/ibm_waning_scenario.R')
#   waning    <- simulate_waning_scenario(ve_floor = 0.1)
#   no_waning <- simulate_waning_scenario(ve_floor = 0.8)  # flat VE after ramp
#   export_waning_scenario(waning, "output/ibm_export_waning")
#
# Output tables follow the three-file schema (demographic_data,
# vaccination_data, incidence_data_reduced) used by the surviver R package
# and the group-D VE vignette in the Lorentz workshop repository.
#
# Requires ibm_core.R to have been sourced first.
############################################################################ #

#' Simulate the waning-VE test scenario
#'
#' Phase 1: initial spread without vaccination for `initial_days` days.
#' Phase 2: `n_chunks` calls to rerun_ibm(), each `chunk_days` long,
#' vaccinating `chunk_coverage` of the remaining eligible pool per chunk with
#' a distinct RNG seed (rng_seed + chunk index).
#'
#' Setting `ve_floor` equal to `vaccine_effectiveness` gives a scenario with
#' no waning (VE ramps up to its peak and then stays constant).
#'
#' @return A list with `demographic_data`, `vaccination_data`,
#'   `incidence_data_reduced` (data frames), `params` (the phase-1 parameter
#'   list, which records all scenario settings) and `study_end` (days).
simulate_waning_scenario <- function(ve_floor              = 0.1,
                                     vaccine_effectiveness = 0.8,
                                     ve_ramp_days          = 14,
                                     ve_half_life          = 45,
                                     pop_size              = 1e4,
                                     transmission_prob     = 0.02,
                                     num_infected_seeds    = 10,
                                     initial_days          = 20,
                                     n_chunks              = 16,
                                     chunk_days            = 10,
                                     chunk_coverage        = 0.03,
                                     output_dir            = tempfile("ibm_waning_"),
                                     verbose               = FALSE) {

  # ---- PHASE 1: INITIAL SPREAD (NO VACCINATION) ----
  params <- get_default_parameters()

  params$output_dir             <- output_dir
  params$pop_size               <- pop_size
  params$num_days               <- initial_days
  params$num_infected_seeds     <- num_infected_seeds
  params$bool_add_baseline      <- FALSE
  params$bool_show_demographics <- FALSE
  # NULL switches on the separate age-based natural mortality process
  # (ibm_pop_death.R); disease deaths are identified by time_of_death
  params$general_mortality_rate <- NULL

  # Slowed epidemic relative to the default (0.1). A sweep at pop_size = 1e4
  # found 0.035 peaks at day 41 and is extinct by day ~126 (before most of
  # the campaign), while 0.02 peaks mid-campaign (day ~88) and is still
  # transmitting at day 180, giving breakthrough events across the full
  # range of days since vaccination.
  params$transmission_prob      <- transmission_prob

  # vaccination is off during phase 1; VE settings apply once vaccinated
  params$vaccine_coverage       <- 0
  params$vaccine_effectiveness  <- vaccine_effectiveness
  params$ve_ramp_days           <- ve_ramp_days
  params$ve_half_life           <- ve_half_life
  params$ve_floor               <- ve_floor

  if (verbose) print_model_parameters(params)

  run_ibm(params, verbose = verbose)

  pop_data         <- readRDS(file.path(params$output_dir, "pop_data.rds"))
  health_time_data <- readRDS(file.path(params$output_dir, "health_time.rds"))

  # ---- PHASE 2: ROLLING VACCINATION CAMPAIGN ----
  for (i in seq_len(n_chunks)) {

    chunk_params <- params
    chunk_params$num_days           <- chunk_days
    chunk_params$num_infected_seeds <- 0
    chunk_params$vaccine_coverage   <- chunk_coverage
    # independent random stream per chunk
    chunk_params$rng_seed           <- params$rng_seed + i

    rerun_ibm(
      pop_data               = pop_data,
      prev_log_health_matrix = health_time_data$log_health,
      params                 = chunk_params,
      intervention           = i,
      verbose                = verbose
    )

    pop_data <- readRDS(file.path(
      chunk_params$output_dir, paste0("pop_data_intervention", i, ".rds")))
    health_time_data <- readRDS(file.path(
      chunk_params$output_dir, paste0("health_time_intervention", i, ".rds")))
  }

  study_end <- initial_days + n_chunks * chunk_days

  # ---- INCIDENCE DATA (WITH REPORTING DELAY) ----
  delay_params            <- get_default_delay_parameters()
  delay_params$output_dir <- params$output_dir
  # get_default_delay_parameters() has no rng_seed, so obtain_incidence_data()
  # would call set.seed(NULL) and draw different reporting delays/missingness
  # on every run; seed it explicitly so the scenario is reproducible
  delay_params$rng_seed   <- params$rng_seed
  incidence_data          <- obtain_incidence_data(pop_data, delay_params)$incidence_data

  # ---- THREE-FILE VE ANALYSIS SCHEMA ----
  demographic_data <- data.frame(
    ind_id       = rownames(pop_data),
    age          = pop_data$age,
    sex          = pop_data$sex,
    hh_id        = pop_data$hh_id,
    classroom_id = pop_data$classroom_id,
    workplace_id = pop_data$workplace_id,
    time         = 0,
    stringsAsFactors = FALSE
  )

  is_vaccinated <- !is.na(pop_data$time_of_vaccination)
  vaccination_data <- data.frame(
    ind_id              = rownames(pop_data)[is_vaccinated],
    time_of_vaccination = pop_data$time_of_vaccination[is_vaccinated],
    time                = study_end,
    stringsAsFactors = FALSE
  )

  incidence_data_reduced <- data.frame(
    ind_id            = rownames(incidence_data),
    time_of_reporting = incidence_data$time_of_reporting,
    time              = study_end,
    stringsAsFactors = FALSE
  )

  list(
    demographic_data       = demographic_data,
    vaccination_data       = vaccination_data,
    incidence_data_reduced = incidence_data_reduced,
    params                 = params,
    study_end              = study_end
  )
}

#' Write a simulated scenario to disk in the three-file schema
#'
#' Writes demographic_data.rds, vaccination_data.rds,
#' incidence_data_reduced.rds and scenario_params.rds (a list of `params`
#' and `study_end`) to `export_dir`.
export_waning_scenario <- function(scenario, export_dir) {
  if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
  saveRDS(scenario$demographic_data,
          file.path(export_dir, "demographic_data.rds"))
  saveRDS(scenario$vaccination_data,
          file.path(export_dir, "vaccination_data.rds"))
  saveRDS(scenario$incidence_data_reduced,
          file.path(export_dir, "incidence_data_reduced.rds"))
  saveRDS(list(params = scenario$params, study_end = scenario$study_end),
          file.path(export_dir, "scenario_params.rds"))
  invisible(export_dir)
}
