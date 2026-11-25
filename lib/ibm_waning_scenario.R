############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: generate an example dataset for testing time-varying vaccine
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
# Output is exported in the same three-file schema (demographic_data.rds,
# vaccination_data.rds, incidence_data_reduced.rds) used by the group-D
# VE-estimation vignette, so it can be dropped in as a drop-in replacement
# for the existing example data.
############################################################################ #

# ------------------------------------------------------------------------ -
# DEPENDENCIES ----
# ------------------------------------------------------------------------ -

rm(list = ls())

source('lib/ibm_core.R')
library(ggplot2)

# ------------------------------------------------------------------------ -
# PHASE 1: INITIAL SPREAD (NO VACCINATION) ----
# ------------------------------------------------------------------------ -

params <- get_default_parameters()

params$output_dir             <- "output/ibm_waning"
params$pop_size                <- 1e5
params$num_days                 <- 20
params$num_infected_seeds       <- 10
params$bool_add_baseline        <- FALSE
params$bool_show_demographics   <- FALSE
params$general_mortality_rate   <- NULL

# slow the epidemic down relative to the default/legacy scenario (0.07-0.1)
# so it stays sustained across the whole rolling vaccination campaign below,
# rather than peaking and going extinct before most people are vaccinated.
# A sweep at pop_size=1e4 found 0.035 peaks at day 41 and is extinct by day
# ~126 (i.e. before, not during, most of the campaign; attack rate 71%),
# while 0.02 peaks later (day ~88, mid-campaign) and sustains a 200-330
# case plateau from ~day 55 through day 180 without going extinct (attack
# rate ~32%) -- giving breakthrough events across the full range of
# days-since-vaccination instead of concentrating them in the first ~45 days.
params$transmission_prob        <- 0.02

# vaccination is off during the initial spread phase
params$vaccine_coverage         <- 0
params$vaccine_effectiveness    <- 0.8   # peak effectiveness once vaccinated
params$ve_ramp_days             <- 14
params$ve_half_life             <- 45
params$ve_floor                 <- 0.1

print_model_parameters(params)

ibm_results <- run_ibm(params, verbose = FALSE)

pop_data_file    <- file.path(ibm_results$params$output_dir, 'pop_data.rds')
health_time_file <- file.path(ibm_results$params$output_dir, 'health_time.rds')

pop_data         <- readRDS(pop_data_file)
health_time_data <- readRDS(health_time_file)

# ------------------------------------------------------------------------ -
# PHASE 2: ROLLING VACCINATION CAMPAIGN ----
# ------------------------------------------------------------------------ -

# roll out vaccination in many small increments (rather than 1-2 large
# pulses) so vaccinated individuals are staggered across many calendar
# days, giving the VE(t) model covariate variation to work with
n_chunks        <- 16
chunk_days      <- 10
chunk_coverage  <- 0.03   # fraction of remaining eligible pool per chunk
base_rng_seed   <- params$rng_seed

for (i in seq_len(n_chunks)) {

  chunk_params <- params
  chunk_params$num_days           <- chunk_days
  chunk_params$num_infected_seeds <- 0
  chunk_params$vaccine_coverage   <- chunk_coverage
  # vary the seed per chunk so each rolling-campaign call draws an
  # independent random stream rather than replaying the same one
  chunk_params$rng_seed           <- base_rng_seed + i

  rerun_results <- rerun_ibm(
    pop_data                = pop_data,
    prev_log_health_matrix  = health_time_data$log_health,
    params                  = chunk_params,
    intervention            = i,
    verbose                 = FALSE
  )

  pop_data_file    <- file.path(chunk_params$output_dir, paste0('pop_data_intervention', i, '.rds'))
  health_time_file <- file.path(chunk_params$output_dir, paste0('health_time_intervention', i, '.rds'))

  pop_data         <- readRDS(pop_data_file)
  health_time_data <- readRDS(health_time_file)
}

study_end <- params$num_days + n_chunks * chunk_days

# ------------------------------------------------------------------------ -
# CREATE FINAL POPULATION DATA ----
# ------------------------------------------------------------------------ -

final_pop_data <- pop_data

saveRDS(list(log_health = health_time_data$log_health, params = params),
        file = file.path(params$output_dir, "final_health_time.rds"))
saveRDS(final_pop_data,
        file = file.path(params$output_dir, "final_pop_data.rds"))

# ------------------------------------------------------------------------ -
# DATASET: INCIDENCE DATA (WITH REPORTING DELAY) ----
# ------------------------------------------------------------------------ -

delay_params <- get_default_delay_parameters()
delay_params$output_dir <- params$output_dir

incidence_results   <- obtain_incidence_data(final_pop_data, delay_params)
final_incidence_data <- incidence_results$incidence_data

saveRDS(final_incidence_data,
        file = file.path(delay_params$output_dir, "final_incidence_data.rds"))

# ------------------------------------------------------------------------ -
# EXPORT: DEMOGRAPHIC / VACCINATION / INCIDENCE FILES FOR VE ANALYSIS ----
# ------------------------------------------------------------------------ -
# schema matches analysis/group-D/data_wrangling.R in the workshop repo

export_dir <- file.path(params$output_dir, "for_ve_analysis")
if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

demographic_data <- data.frame(
  ind_id       = rownames(final_pop_data),
  age          = final_pop_data$age,
  sex          = final_pop_data$sex,
  hh_id        = final_pop_data$hh_id,
  classroom_id = final_pop_data$classroom_id,
  workplace_id = final_pop_data$workplace_id,
  time         = 0,
  stringsAsFactors = FALSE
)

is_vaccinated <- !is.na(final_pop_data$time_of_vaccination)
vaccination_data <- data.frame(
  ind_id              = rownames(final_pop_data)[is_vaccinated],
  time_of_vaccination = final_pop_data$time_of_vaccination[is_vaccinated],
  time                = study_end,
  stringsAsFactors = FALSE
)

incidence_data_reduced <- data.frame(
  ind_id            = rownames(final_incidence_data),
  time_of_reporting = final_incidence_data$time_of_reporting,
  time              = study_end,
  stringsAsFactors = FALSE
)

saveRDS(demographic_data,        file = file.path(export_dir, "demographic_data.rds"))
saveRDS(vaccination_data,        file = file.path(export_dir, "vaccination_data.rds"))
saveRDS(incidence_data_reduced,  file = file.path(export_dir, "incidence_data_reduced.rds"))

cat("\nExported VE analysis input files to:", export_dir, "\n")
cat("Study length (days):     ", study_end, "\n")
cat("Population size:         ", nrow(demographic_data), "\n")
cat("Vaccinated:              ", nrow(vaccination_data),
    sprintf("(%.1f%%)", 100 * nrow(vaccination_data) / nrow(demographic_data)), "\n")
cat("Reported cases:          ", nrow(incidence_data_reduced),
    sprintf("(%.1f%%)", 100 * nrow(incidence_data_reduced) / nrow(demographic_data)), "\n")

# ------------------------------------------------------------------------ -
# DIAGNOSTICS ----
# ------------------------------------------------------------------------ -

fig_dir <- file.path(params$output_dir, "final_figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# epidemic curve
epi_curve_df <- data.frame(
  day = seq_len(nrow(health_time_data$log_health)),
  I   = health_time_data$log_health$I
)

p_epi <- ggplot(epi_curve_df, aes(x = day, y = I)) +
  geom_line() +
  labs(x = "Day", y = "Proportion infectious", title = "Epidemic curve (slowed scenario)")

ggsave(file.path(fig_dir, "waning_scenario_epidemic_curve.png"), p_epi, width = 6, height = 4)

# cumulative vaccination coverage over calendar time
vacc_by_day <- data.frame(day = sort(unique(vaccination_data$time_of_vaccination)))
vacc_by_day$cum_vaccinated <- sapply(vacc_by_day$day, function(d) {
  sum(vaccination_data$time_of_vaccination <= d)
})
vacc_by_day$cum_coverage <- vacc_by_day$cum_vaccinated / nrow(demographic_data)

p_vacc <- ggplot(vacc_by_day, aes(x = day, y = cum_coverage)) +
  geom_step() +
  labs(x = "Day", y = "Cumulative coverage", title = "Rolling vaccination campaign")

ggsave(file.path(fig_dir, "waning_scenario_vaccination_coverage.png"), p_vacc, width = 6, height = 4)

# events among vaccinated, by days since vaccination (checks whether there
# is enough vaccinated-arm follow-up to identify waning out to ve_half_life)
vax_merged <- merge(incidence_data_reduced, vaccination_data, by = "ind_id")
vax_merged$days_since_vacc <- vax_merged$time_of_reporting - vax_merged$time_of_vaccination
vax_merged <- vax_merged[vax_merged$days_since_vacc >= 0, ]

cat("\nVaccinated-arm breakthrough events by days-since-vaccination bucket:\n")
print(table(cut(vax_merged$days_since_vacc,
                 breaks = c(0, 7, 14, 30, 45, 60, 90, Inf),
                 right = FALSE)))
