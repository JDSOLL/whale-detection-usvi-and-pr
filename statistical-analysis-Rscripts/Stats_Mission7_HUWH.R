##########################
# Author: Joshua Soll, Ocean Glider Lab, University of the Virgin Islands; joshuasoll440@gmail.com
# Filename: Stats_Mission7_HUWH.R
# Date last modified: Feb 9, 2026
# Description: Version six of a statistical analysis of the effect of cyclical time on manually verified humpback whale call detections for Mission 7, AKA UVI2_AroundSTX_20250306.
# This script includes code for:
#   A hurdle model analysis (zero: binomial; count: truncated NB) with sin/cos(hour) only, and day-block bootstrap CIs (95%) for the cyclical predictions
#   A Rayleigh test for clustering of detections across the 24-hr diel cycle.
##########################

##### Load packages #####
library(tidyverse)
library(lubridate) #parsing and working with date-times (hour(), with_tz()).
library(circular) #circular statistics (Rayleigh test & mean direction).
library(pscl) #hurdle() model (zero + count components).
library(suncalc) #sunrise/sunset times to classify bins as day vs night.
library(ggplot2)
library(cowplot) #combine multiple ggplots into a single figure


##### 1) Load Acoustics Data #####
## Reads the manually-annotated humpback detection spreadsheet for Mission .
# Each row should represent one analysis time bin (e.g., an hour) with start_date / start_time (bin start), and num_true_positives_pitchtracked (count).
# stringsAsFactors = FALSE prevents automatic conversion of text columns to factors.
df <- read.csv("annotated-detection-spreadsheets/Mission6_HUWH.csv", stringsAsFactors = FALSE)

##### 2) Convert datetime to local time (from UTC 0) + assigning day/night  #####
# Build a proper POSIXct datetime in UTC since the glider recorded data in UTC 0, then convert to local AST/UTC-4 (St. Thomas).
# Doing this because the diel/cyclical analyses must be in local time to align with sunlight/night.

df <- df %>%
  mutate(
    start_datetime_utc = as.POSIXct(paste(start_date, start_time),
                                    tz = "UTC", format = "%m/%d/%y %H:%M:%S"),
    start_datetime_ast = with_tz(start_datetime_utc, tzone = "America/St_Thomas")
  )

# Extract unique local dates so we only query sunrise/sunset once per date.
unique_dates <- unique(as.Date(df$start_datetime_ast))

# Compute nautical dawn/dusk for each day at the study location.
# Using nautical twilight is often a better proxy for 'usable light' than sunrise/sunset.
sun_times <- getSunlightTimes(
  date = unique_dates,
  lat = 18.15,
  lon = -64.8,
  keep = c("nauticalDawn", "nauticalDusk"),
  tz = "America/St_Thomas"
)

# Join the dawn/dusk table back onto the main data by local date, so each row/bin can be classified as day vs night.
df <- df %>%
  mutate(date_local = as.Date(start_datetime_ast)) %>%
  left_join(
    sun_times[, c("date", "nauticalDawn", "nauticalDusk")],
    by = c("date_local" = "date")
  ) %>%
  mutate(
    day_night = ifelse(start_datetime_ast >= nauticalDawn &
                         start_datetime_ast < nauticalDusk, "day", "night")
  )


##### 3) Add circular hour variables #####
# Convert local time-of-day to a circular representation on unit circle [0, 2π).
# sin_hour/cos_hour let the model represent a smooth 24h cycle without a hard break at midnight, where midnight is closest to 1:00 a.m., not furthest.
df <- df %>%
  mutate(
    hour_local = hour(start_datetime_ast),
    hour_rad   = 2 * pi * hour_local / 24,
    sin_hour   = sin(hour_rad),
    cos_hour   = cos(hour_rad)
  )

##### 4) Presence/absence column #####
# Presence is 1 if any manually verified/true detections occurred in the bin, else 0.
# This is used for the 'zero' component of the hurdle model and for plotting raw 0/1 jitter.
df <- df %>%
  mutate(presence = ifelse(num_true_positives_pitchtracked > 0, 1, 0))


##### 5) Keep full 24-hour days only #####
# The deployment and recovery day of all these glider missions have partial coverage (missing hours).
# For diel analyses, partial days can bias hourly patterns so we drop the first and last local day, keeping only complete interior days.
first_full_day <- as.Date(min(df$start_datetime_ast)) + 1
last_full_day  <- as.Date(max(df$start_datetime_ast)) - 1

df_full_days <- df %>%
  filter(date_local >= first_full_day & date_local <= last_full_day)

##### 6) Hurdle model: cyclical time (sin_hour + cos_hour) + day-block bootstrap for confidence intervals #####
# Hurdle model rationale:
#   Many bins have zero detections (zero-inflation), and positive bins have counts and are overdispersed.
#   A hurdle model fits:
#     1) a binomial model for Pr(Y > 0) (the 'zero' part)
#     2) a truncated count model for E[Y | Y > 0] (the 'count' part)
#     Then combines them to get E[Y] overall.
# Predictors: sin_hour + cos_hour capture the cyclical (24h) effect only. sin_hour and cos_hour can individually be significant.


# Fit model on full data
hurdle_model_cyclical <- hurdle(
  num_true_positives_pitchtracked ~ sin_hour + cos_hour,
  data = df_full_days,
  dist = "negbin",
  zero.dist = "binomial"
)

summary(hurdle_model_cyclical)

# Prediction grid (hourly resolution over 0–24 hrs)
# Create a fine-resolution time grid (0 to <24h) for smooth predicted curves.
# 'by = 0.1' produces ~240 points (every 6 minutes). Increase/decrease as desired.
pred_grid_cyc <- tibble(hour_local = seq(0, 23.9, by = 0.1)) %>%
  mutate(
    sin_hour = sin(2*pi*hour_local/24),
    cos_hour = cos(2*pi*hour_local/24)
  )

# Helper function to compute Pr(Y>0) robustly from hurdle via P(Y=0)
# 'type='response'' returns the overall expected value E[Y].
# 'type='count'' returns the positive mean E[Y | Y>0].
# For a hurdle model, E[Y] = Pr(Y>0) * E[Y|Y>0], so Pr(Y>0) = E[Y] / E[Y|Y>0].
# Clamp to [0,1] and drop non-finite values for numerical safety.
pred_p_nonzero <- function(mod, newdata) {
  mu_overall <- as.numeric(predict(mod, newdata = newdata, type = "response"))  # E[Y]
  mu_pos     <- as.numeric(predict(mod, newdata = newdata, type = "count"))     # E[Y|Y>0]
  p          <- mu_overall / mu_pos                                             # P(Y>0)
  p[!is.finite(p)] <- NA_real_
  p[p < 0] <- 0; p[p > 1] <- 1
  p
}


# Point estimates on original fit
pred_cyc <- pred_grid_cyc %>%
  mutate(
    p_nonzero  = pred_p_nonzero(hurdle_model_cyclical, pred_grid_cyc),  # Pr(Y>0)
    mu_pos     = predict(hurdle_model_cyclical, newdata = pred_grid_cyc, type = "count"),    # E[Y|Y>0]
    mu_overall = predict(hurdle_model_cyclical, newdata = pred_grid_cyc, type = "response")  # E[Y]
  )


##### Day-block bootstrap for ribbons #####
# The day-block bootstrap resamples whole day bins (not individual rows) to preserve within-day structure and quantifies uncertainty in the cyclical prediction curves.
# B <- 1000 replicates 95% pointwise ribbons via 2.5% and 97.5% quantiles.
set.seed(2025)
B <- 1000

days <- sort(unique(df_full_days$date_local))
n_days <- length(days)
n_x    <- nrow(pred_grid_cyc)

# The matrices store bootstrap predictions for each grid point (rows) across bootstrap replicates (cols).
# Later we take quantiles across columns to form confidence ribbons.
p_nonzero_mat  <- matrix(NA_real_, nrow = n_x, ncol = B)
mu_pos_mat     <- matrix(NA_real_, nrow = n_x, ncol = B)
mu_overall_mat <- matrix(NA_real_, nrow = n_x, ncol = B)

# Each iteration builds a bootstrap dataset by sampling the entire days with replacement.
# This keeps the autocorrelation/structure within a day intact.
# resamples days with replacement and rebuild dataset (with duplicates)
for (b in seq_len(B)) {
  samp_days <- sample(days, size = n_days, replace = TRUE)
  df_b <- bind_rows(lapply(samp_days, function(d) {
    df_full_days %>% filter(date_local == d)
  }))
  
  # Fit hurdle; skip replicate on failure
  # 'dist='negbin'' sets the count distribution for positive counts (negative binomial).
  # 'zero.dist='binomial'' models the zero hurdle (presence/absence).
  fit_b <- try(
    hurdle(num_true_positives_pitchtracked ~ sin_hour + cos_hour,
           data = df_b, dist = "negbin", zero.dist = "binomial"),
    silent = TRUE
  )
  if (inherits(fit_b, "try-error")) next
  
  # Predict on grid
  p_nonzero_mat[, b]  <- pred_p_nonzero(fit_b, pred_grid_cyc)
  mu_pos_mat[, b]     <- predict(fit_b, newdata = pred_grid_cyc, type = "count")
  mu_overall_mat[, b] <- predict(fit_b, newdata = pred_grid_cyc, type = "response")
}

# 'qfun' computes the 2.5% and 97.5% quantiles (95% interval).
# 'apply(..., 1, qfun)' computes these intervals at each hour_grid point.
qfun <- function(v) { quantile(v, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE) }

# Rowwise quantiles for ribbons
pn_lo_hi  <- t(apply(p_nonzero_mat,  1, qfun))
mp_lo_hi  <- t(apply(mu_pos_mat,     1, qfun))
mo_lo_hi  <- t(apply(mu_overall_mat, 1, qfun))

pred_cyc <- pred_cyc %>%
  mutate(
    p_nonzero_lo  = pn_lo_hi[,1],  p_nonzero_hi  = pn_lo_hi[,2],
    mu_pos_lo     = mp_lo_hi[,1],  mu_pos_hi     = mp_lo_hi[,2],
    mu_overall_lo = mo_lo_hi[,1],  mu_overall_hi = mo_lo_hi[,2]
  )

##### Observed summaries by hour (for dots) #####
# These summaries are optional 'dots' to show observed data by hour_bin.
# 'prop_present': observed Pr(Y>0) by hour (binned).
# 'mean_pos': mean count among positive bins only, per hour (binned).
# 'mean_all': mean count including zeros, per hour (binned).
obs_by_hour_all <- df_full_days %>%
  mutate(hour_bin = floor(hour_local)) %>%
  group_by(hour_bin) %>%
  summarise(
    prop_present = mean(presence == 1),
    mean_pos     = if (any(num_true_positives_pitchtracked > 0)) {
      mean(num_true_positives_pitchtracked[num_true_positives_pitchtracked > 0])
    } else NA_real_,
    mean_all     = mean(num_true_positives_pitchtracked),
    .groups = "drop"
  )

##### option to toggle between raw or binned dots #####
# 'USE_RAW <- TRUE' plots individual-hour bins (raw points/jitter) for transparency.
# Setting FALSE would plot only binned/hourly summaries, which would mean less clutter but less detail, not best for publication.
USE_RAW <- TRUE

if (USE_RAW) {
  
  # ZERO part: raw 0/1 presence by hour (light vertical jitter)
  plt_cyc_pnonzero <- ggplot() +
    geom_ribbon(data = pred_cyc,
                aes(hour_local, ymin = p_nonzero_lo, ymax = p_nonzero_hi),
                alpha = 0.25, fill = "grey70") +
    geom_jitter(
      data = df_full_days,
      aes(x = hour_local, y = presence),
      height = 0.05, width = 0, alpha = 0.15, size = 1
    ) +
    geom_line(data = pred_cyc, aes(hour_local, p_nonzero),
              size = 1.2, color = "tomato2") +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0,1), labels = c("No","Yes")) +
    scale_x_continuous(breaks = 0:23) +
    labs(x = "Hour of Day", y = "Any detection (raw 0/1)") +
    theme_classic()
  
  # COUNT part: raw positive counts by hour
  plt_cyc_meanpos <- ggplot() +
    geom_ribbon(data = pred_cyc,
                aes(hour_local, ymin = mu_pos_lo, ymax = mu_pos_hi),
                alpha = 0.25, fill = "grey70") +
    geom_point(
      data = dplyr::filter(df_full_days, presence == 1),
      aes(x = hour_local, y = num_true_positives_pitchtracked),
      alpha = 0.2, size = 1.2
    ) +
    geom_line(data = pred_cyc, aes(hour_local, mu_pos),
              size = 1.2, color = "tomato2") +
    scale_x_continuous(breaks = 0:23) +
    labs(x = "Hour of Day", y = "Calls (when present, raw)") +
    theme_classic()
  
  # Derived part: raw counts (including zeros)
  plt_cyc_overall <- ggplot() +
    geom_ribbon(data = pred_cyc,
                aes(hour_local, ymin = mu_overall_lo, ymax = mu_overall_hi),
                alpha = 0.25, fill = "grey70") +
    geom_point(
      data = df_full_days,
      aes(x = hour_local, y = num_true_positives_pitchtracked),
      alpha = 0.15, size = 1.1
    ) +
    geom_line(data = pred_cyc, aes(hour_local, mu_overall),
              size = 1.2, color = "tomato2") +
    scale_x_continuous(breaks = 0:23) +
    labs(x = "Hour of Day", y = "Calls per hour (raw)") +
    theme_classic()
  
} else {
  # binned dots here as alternative (if using FALSE, though not as good, don't recommend)
  plt_cyc_pnonzero <- plt_cyc_pnonzero
  plt_cyc_meanpos  <- plt_cyc_meanpos
  plt_cyc_overall  <- plt_cyc_overall
}

##### Prepare assembly of hurdle figure overall #####
# Combine the three hurdle outputs into one multi-panel figure:
#   A) Pr(Y>0) across hour
#   B) Mean count when present, E[Y|Y>0]
#   C) Overall mean count, E[Y]
# cowplot::plot_grid aligns and labels panels.
cyc_combo <- cowplot::plot_grid(
  plt_cyc_pnonzero, plt_cyc_meanpos, plt_cyc_overall,
  labels = c("A","B","C"), ncol = 3, align = "hv"
)

print(cyc_combo)

ggsave(
  "/insert/desired/directory.png",
  plot   = cyc_combo,
  width  = 18,
  height = 6,
  dpi    = 600
)

##### 7) Rayleigh test (only using days where sound was recorded all 24 hours, essentially cutting off deployment day and recovery day of glider) #####
# Tests whether angles (hours) are uniformly distributed around the clock.
# Here we test whether detections cluster at particular times of day.
# We build a circular vector of hours (in radians), then optionally weight by call counts.

rayleigh_only_when_present <- TRUE
# If TRUE, the test uses only bins with detections (presence==1).
# If FALSE, it uses all bins (including zeros) which changes the interpretation. Using TRUE here for this work.

df_rayleigh <- if (rayleigh_only_when_present) {
  df_full_days %>% filter(presence == 1)
} else {
  df_full_days
}

# Build a circular 'angle' object from hour_rad (radians), using a 24-hour clock template.
# This is required for circular::rayleigh.test().
hour_circ <- circular(df_rayleigh$hour_rad, type = "angles",
                      units = "radians", template = "clock24")

# Weighting: replicate each hour angle by the number of calls in that bin.
# This makes hours with higher call counts contribute more to the test and mean direction.
hour_circ_weighted <- rep(hour_circ, times = round(df_rayleigh$num_true_positives_pitchtracked))

if (length(hour_circ_weighted) > 0) {
  rayleigh_result <- rayleigh.test(hour_circ_weighted)
  print(rayleigh_result)
  mean_dir_deg <- as.numeric(mean(hour_circ_weighted)) * 180 / pi
  cat("Mean direction (degrees):", round(mean_dir_deg, 1), "\n")
} else {
  cat("Rayleigh test skipped: no non-zero bins to test.\n")
}

#Rayleigh visualization
# Visualization: polar bar chart of detections (or percent of detections) by hour bin.
# The red line indicates the mean direction (peak clustering time, and is adjusted manually after getting mean direction result for interpretability).
# Assume df_full_days exists and has hour_local and num_true_positives_pitchtracked
df_present <- df_full_days %>%
  filter(num_true_positives_pitchtracked > 0) %>%
  mutate(hour_bin = floor(hour_local)) %>%
  count(hour_bin, name = "n_calls") %>%
  mutate(
    hour_bin = factor(hour_bin, levels = 0:23, labels = sprintf("%02d:00", 0:23)),
    percent = n_calls / sum(n_calls) * 100
  )

library(ggplot2)
library(dplyr)

# Using existing df_present with hour_bin as factor from 0:23
# Convert mean direction (in degrees) to x-position on discrete scale
mean_deg <- 175.5
mean_x <- mean_deg / 14.35  # 15 degrees per hour, but adjusting for visuals because the rayleigh plot in this code was shifted to put 00:00 at the top of the plot. 14.35 deg to appear as 11:42 a.m.
# 'mean_x' is a manual adjustment to map degrees to the discrete hour factor scale.

# Maximum radius for line (max percent)
max_radius <- max(df_present$percent) * 1.1  # a bit beyond max bar

# 'coord_polar' transforms the bar chart into a circular clock-like plot.
# 'start = -pi/24' rotates so that 00:00 sits at the top (adjust if you want a different orientation).
rayleigh_plot <- ggplot(df_present, aes(x = hour_bin, y = percent)) +
  geom_col(fill = "skyblue", color = "black", size = 0.1, width = 1) +
  coord_polar(start = -pi/24, direction = 1) +
  scale_y_sqrt() +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 9),
    legend.position = "bottom",
    plot.background = element_rect(fill = alpha("cornsilk", 0.5), color = NA),
    plot.margin = unit(c(10, 10, 10, 10), "pt"),
    plot.title = element_text(size = 14, hjust = 0.5, vjust = 5)
  ) +
  
  # Add red mean direction line with mapped color for legend
  geom_segment(aes(x = mean_x, xend = mean_x, y = 0, yend = max_radius, color = "Mean Direction"),
               size = 1, lineend = "round") +
  
  # Manually define legend color
  scale_color_manual(name = NULL, values = c("Mean Direction" = "tomato2"))

print(rayleigh_plot)

#Save plot
ggsave("insert/desired/directory.pdf", plot = rayleigh_plot, width = 10, height = 8, dpi = 300)

##### 8) Summary stats #####
# Quick mission-level summary numbers for reporting.
# precision_score compares true detections to algorithm detections (as provided in df$sum).
# Ensure df$sum exists and represents total algorithm detections per bin.
summary_stats <- list(
  total_days               = length(unique(df$date_local)),
  days_with_detections     = length(unique(df$date_local[df$num_true_positives_pitchtracked > 0])),
  total_hours              = nrow(df),
  hours_with_detections    = sum(df$num_true_positives_pitchtracked > 0),
  total_true_detections    = sum(df$num_true_positives_pitchtracked, na.rm = TRUE),
  total_algorithm_detections = sum(df$sum, na.rm = TRUE)
)
summary_stats$precision_score <- with(summary_stats,
                                      total_true_detections / total_algorithm_detections)
summary_stats

######### End of script #########