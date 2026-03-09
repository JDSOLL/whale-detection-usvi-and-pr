##########################
# Author: Joshua Soll, Ocean Glider Lab, University of the Virgin Islands; joshuasoll440@gmail.com
# Filename: Effect_of_Depth_on_Detections.R
# Date last modified: Mar 09, 2026
# Description: Effects of depth on whale detections using a hurdle model only.
#   Computes hourly mean depth
#   Fits hurdle model: occurrence (Bernoulli/logit) + intensity given presence (truncated NB)
#   Plots P(Y>0), E[Y|Y>0], and overall E[Y] vs depth
#   OVERLAY = RAW HOURLY OBSERVATIONS (not bin means)
##########################

######### Libraries #########
# Purpose: Load only the packages needed for this analysis, while keeping the console clean.
# Why suppressPackageStartupMessages(): reduces clutter so warnings/errors are easier to spot.
# General package roles:
#   dplyr/tidyr/tibble: data wrangling + tidy data frames.
#   lubridate: timezone-safe datetime parsing + hourly binning.
#   ggplot2/cowplot: plotting + multi-panel figure layouts.
#   readr: fast/consistent CSV reading (even though read.csv is used below).
#   ncdf4: read glider netCDF time series (depth/time/lat/lon).
#   pscl: hurdle() models (zero part + positive-count part).
#   abind: stack bootstrap replicate matrices into a 3D array.
######### Libraries #########
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(readr)
  library(ncdf4)
  library(pscl)
  library(cowplot)
  library(tibble)
  library(abind)   
})

######### User Toggles (edit these for mission/species) #########
# Purpose: Central place to point the script at a specific mission + species and set output labels.
# Why this matters: most downstream joins and plots assume these paths/times are correct.
# Note: This code can be used for all four datasets. Keep only ONE whale_df read.csv() line active at a time (comment out the others).

# Whale detections (hourly)
#whale_df <- read.csv("annotated-detection-spreadsheets/Mission7_HUWH.csv", stringsAsFactors = FALSE)
whale_df <- read.csv("annotated-detection-spreadsheets/Mission7_MIWH.csv", stringsAsFactors = FALSE)
#whale_df <- read.csv("annotated-detection-spreadsheets/Mission6_HUWH.csv", stringsAsFactors = FALSE)
#whale_df <- read.csv("annotated-detection-spreadsheets/Mission6_MIWH.csv", stringsAsFactors = FALSE)

# Glider netCDF path
nc_path <- "mission-nc-files/Mission7_TimeSeries.nc"
#nc_path <- "mission-nc-files/Mission6_TimeSeries.nc"

# Labels for filenames
mission_label <- "Mission 7"
species_label <- "Minke Whale"

# Timezone (USVI)
tz_loc <- "America/St_Thomas"

# Define start/end in local time (adjust per mission); toggle comment out/in the associated mission start/end time as needed
# Mission 6
#start_time <- as.POSIXct("2025-02-11 00:00:00", tz = tz_loc)
#end_time   <- as.POSIXct("2025-03-04 23:59:59", tz = tz_loc)

# Mission 7
start_time <- as.POSIXct("2025-03-06 00:00:00", tz = tz_loc)
end_time   <- as.POSIXct("2025-03-25 23:59:59", tz = tz_loc)


# Depth post-processing toggles
scale_depth_by_10 <- TRUE # Mission 7 depth data needs to be multiplied by 10 because of format the data was recorded in. Set as FALSE for Mission 6
cap_depth_max     <- 1100 # filter out extreme values (also catches Mission 6 outliers)

# Output directories
out_dir_base <- "/Insert/Directory/Here"
out_dir_hz   <- file.path(out_dir_base, "hurdle")
dir.create(out_dir_hz, recursive = TRUE, showWarnings = FALSE)

######### 1) Whale detections (hourly) #########
# Purpose: Convert the detection CSV's date+time columns into a single POSIXct timestamp,
#   localize to USVI time, and create an hourly bin key used for joining with glider depth.
# Key idea: hour_bin is the "join key" that aligns detections with the *hourly mean* depth.
whale_df <- whale_df %>%
  mutate(
    start_datetime_utc = mdy_hms(paste(start_date, start_time), tz = "UTC"),
    start_datetime_loc = with_tz(start_datetime_utc, tzone = tz_loc),
    hour_bin = floor_date(start_datetime_loc, unit = "hour")
  )

######### 2) Glider time series #########
# Purpose: Read glider depth/time/lat/lon from the netCDF time series and subset to the mission time window (in local time).
# Why: whale detections are binned hourly; we need depth aligned to those same hours.
nc <- nc_open(nc_path)  # Open the netCDF file connection (required before ncvar_get()).
on.exit(try(nc_close(nc), silent = TRUE), add = TRUE)  # Safety: ensures the file handle closes even if the script errors.

glider_time <- ncvar_get(nc, "time") # netCDF time vector (seconds since origin).
depth <- ncvar_get(nc, "depth") # Depth measurements (units depend on mission export; may need scaling).
latitude <- ncvar_get(nc, "latitude")  # Latitude (not used in model; retained for possible QC/plots).
longitude <- ncvar_get(nc, "longitude") # Longitude (not used in model; retained for possible QC/plots).
nc_close(nc) # Closes the .nc file

glider_df <- data.frame(
  time_utc  = as.POSIXct(glider_time, origin = "1970-01-01", tz = "UTC"),
  depth     = as.numeric(depth),
  latitude  = as.numeric(latitude),
  longitude = as.numeric(longitude)
) %>%
  mutate(time_loc = with_tz(time_utc, tzone = tz_loc)) %>%
  filter(time_loc >= start_time & time_loc <= end_time) %>%
  mutate(
    depth = if (scale_depth_by_10) pmax(depth * 10, 0) else pmax(depth, 0),  # Optional unit fix (some exports are in 0.1 m).
    depth = if (!is.na(cap_depth_max)) pmin(depth, cap_depth_max) else depth  # Optional cap to drop extreme/outlier depths.
  )

# Hourly mean depth
# Purpose: Collapse high-frequency glider depth samples into a single mean depth per hour.
# Why mean(): gives a representative depth for each detection bin; reduces noise/oversampling.
glider_hourly <- glider_df %>%
  mutate(hour_bin = floor_date(time_loc, unit = "hour")) %>%
  group_by(hour_bin) %>%
  summarise(
    depth_avg = pmax(mean(depth, na.rm = TRUE), 0),
    .groups = "drop"
  )

######### 3) Merge whale + glider datasets; build cyclical time #########
# Purpose: Join hourly detections with hourly mean depth and add diel-cycle covariates.
# Diel cycle encoding: sin/cos transforms represent a 24-hr cycle without a discontinuity (plotting on a unit circle) at midnight (i.e., hour 23 is closest to hour 0 in circular time, not the furthest).
ts_df <- whale_df %>%
  left_join(glider_hourly, by = "hour_bin") %>%
  filter(!is.na(depth_avg)) %>%
  mutate(
    depth_avg = pmax(depth_avg, 0),
    hour_local = hour(hour_bin) + minute(hour_bin) / 60,
    sin_hour = sin(2 * pi * hour_local / 24),
    cos_hour = cos(2 * pi * hour_local / 24)
  )

######### 4) Hurdle model (mean depth only) #########
# Purpose: Fit a hurdle model for hourly call detections vs depth while controlling for diel cycle.
# Model structure (pscl::hurdle):
#   (1) Zero/occurrence part: Pr(Y > 0) ~ logit(depth + sin + cos)
#   (2) Positive-count part: E[Y | Y > 0] ~ truncated NegBin(depth + sin + cos)
# Why hurdle model: call counts often have many zeros (zero inflation) + overdispersed positive counts.
hurdle_mean <- hurdle(
  num_true_positives_pitchtracked ~ depth_avg + sin_hour + cos_hour,
  data = ts_df, dist = "negbin"
)
cat("\n--- Hurdle (mean-only) ---\n"); print(summary(hurdle_mean))

######### 5) Predictions over depth grid #########
# Purpose: Generate smooth predicted curves across the observed depth range.
# Approach:
#   1) Create a depth grid (depth_seq).
#   2) Hold time-of-day effects fixed at their sample means (sin_bar/cos_bar) so the curves represent an "average" diel condition rather than a specific hour.
#   3) Extract:
#      - p_nonzero: Pr(Y > 0)
#      - mu_pos:    E[Y | Y > 0]
#      - mu_overall:E[Y]
# Note on "type="prob"": returns Pr(Y = k) for k in 0:ymax, so we back out Pr(Y > 0) = 1 - Pr(Y=0).
depth_seq <- seq(min(ts_df$depth_avg, na.rm = TRUE),
                 max(ts_df$depth_avg, na.rm = TRUE),
                 length.out = 200)

# Hold time effects at sample means
sin_bar <- mean(ts_df$sin_hour, na.rm = TRUE)
cos_bar <- mean(ts_df$cos_hour, na.rm = TRUE)

newdata_depth <- tibble(
  depth_avg = depth_seq,
  sin_hour  = sin_bar,
  cos_hour  = cos_bar
)

# Hurdle predictions
ymax <- max(ts_df$num_true_positives_pitchtracked, na.rm = TRUE)
prob_mat <- predict(hurdle_mean, newdata = newdata_depth, type = "prob", at = 0:ymax)  # Pr(Y=k) for k=0..ymax.
p0 <- prob_mat[, "0"]
p_nonzero <- 1 - p0
mu_pos <- predict(hurdle_mean, newdata = newdata_depth, type = "count")     # E[Y|Y>0]
mu_overall <- predict(hurdle_mean, newdata = newdata_depth, type = "response")  # E[Y]

pred_df <- newdata_depth %>%
  transmute(depth_avg,
            p_nonzero = p_nonzero,
            mu_pos = mu_pos,
            mu_overall = mu_overall)

######### 6) Day-block bootstrap ribbons #########
# Purpose: Add uncertainty bands (95% bootstrap CIs) around the depth-effect predictions.
# Why "day-block" bootstrap: preserves within-day autocorrelation/structure by resampling whole days (blocks) rather than individual hours (which would overstate effective sample size).
# Output: pred_df gets lower/upper ribbons for each of the three quantities:
#   Pr(Y>0), E[Y|Y>0], and E[Y].
set.seed(42)
ts_df <- ts_df %>% mutate(day_block = as.Date(hour_bin))

boot_once_block <- function() {
  # Purpose: One bootstrap replicate:
  #   1) Resamples days with replacement.
  #   2) Refits the hurdle model on the resampled dataset.
  #   3) Predicts on the SAME depth grid (newdata_depth) so replicates are comparable.
  # Returns: a 200 x 3 matrix (p_nonzero, mu_pos, mu_overall) for this replicate, or NULL if fit fails.
  samp_days <- sample(unique(ts_df$day_block),
                      replace = TRUE,
                      size = length(unique(ts_df$day_block)))
  dat_b <- dplyr::bind_rows(lapply(samp_days, function(d) ts_df[ts_df$day_block == d, ]))
  fit_b <- try(pscl::hurdle(
    num_true_positives_pitchtracked ~ depth_avg + sin_hour + cos_hour,
    data = dat_b, dist = "negbin"
  ), silent = TRUE)
  if (inherits(fit_b, "try-error")) return(NULL)
  ymax_b <- max(dat_b$num_true_positives_pitchtracked, na.rm = TRUE)
  if (!is.finite(ymax_b) || ymax_b < 0) ymax_b <- 0
  prob_b   <- predict(fit_b, newdata = newdata_depth, type = "prob", at = 0:ymax_b)
  if (is.null(colnames(prob_b)) && ncol(prob_b) >= 1) {
    p0_b <- prob_b[, 1]
  } else {
    p0_b <- prob_b[, "0"]
  }
  cbind(
    p_nonzero_b = 1 - p0_b,
    mu_pos_b    = predict(fit_b, newdata = newdata_depth, type = "count"),
    mu_over_b   = predict(fit_b, newdata = newdata_depth, type = "response")
  )
}

B <- 1000  # for good plot quality, do 1000+
boot_list <- replicate(B, boot_once_block(), simplify = FALSE)
boot_list <- boot_list[!vapply(boot_list, is.null, logical(1))]

if (length(boot_list) < 10)
  warning("Few successful bootstrap fits; CI may be unstable (consider reducing B or revisiting blocks).")

boot_arr <- abind::abind(boot_list, along = 3)

q_lo <- 0.025; q_hi <- 0.975
p_nonzero_lo  <- apply(boot_arr[, 1, , drop = FALSE], 1, quantile, q_lo, na.rm = TRUE)
p_nonzero_hi  <- apply(boot_arr[, 1, , drop = FALSE], 1, quantile, q_hi, na.rm = TRUE)
mu_pos_lo     <- apply(boot_arr[, 2, , drop = FALSE], 1, quantile, q_lo, na.rm = TRUE)
mu_pos_hi     <- apply(boot_arr[, 2, , drop = FALSE], 1, quantile, q_hi, na.rm = TRUE)
mu_overall_lo <- apply(boot_arr[, 3, , drop = FALSE], 1, quantile, q_lo, na.rm = TRUE)
mu_overall_hi <- apply(boot_arr[, 3, , drop = FALSE], 1, quantile, q_hi, na.rm = TRUE)

pred_df <- pred_df %>%
  mutate(
    p_nonzero_lo  = pmax(0, pmin(1, p_nonzero_lo)),
    p_nonzero_hi  = pmax(0, pmin(1, p_nonzero_hi)),
    mu_pos_lo     = pmax(0, mu_pos_lo),
    mu_pos_hi     = pmax(0, mu_pos_hi),
    mu_overall_lo = pmax(0, mu_overall_lo),
    mu_overall_hi = pmax(0, mu_overall_hi)
  )

######### 7) p-values (kept for reporting, not plotted) #########
# Purpose: Pull depth term p-values from the model summary for quick reporting.
# Note: hurdle models have separate coefficients/p-values for the zero and count components.
sum_hz <- summary(hurdle_mean)
pval_hz_zero  <- tryCatch(coef(sum_hz)$zero["depth_avg", "Pr(>|z|)"], error = function(e) NA_real_)
pval_hz_count <- tryCatch(coef(sum_hz)$count["depth_avg", "Pr(>|z|)"], error = function(e) NA_real_)
cat("Zero-part p(depth): ", pval_hz_zero, "\n")
cat("Count-part p(depth):", pval_hz_count, "\n")

######### 8) Figures of raw hourly overlays #########
# Purpose: Plot modeled depth relationships while showing the raw hourly data underneath.
# Panels:
#   A) Occurrence: Pr(Calls > 0) vs depth  (raw shown as 0/1 jittered points)
#   B) Intensity: Calls when present (>0) vs depth
#   C) Overall: Calls per hour (including zeros) vs depth

# Zero plot: raw 0/1 hourly detections (jittered vertically) + curve + ribbon
plt_hz_p <- ggplot() +
  geom_ribbon(data = pred_df,
              aes(depth_avg, ymin = p_nonzero_lo, ymax = p_nonzero_hi),
              alpha = 0.25) +
  geom_line(data = pred_df, aes(depth_avg, p_nonzero), linewidth = 1.1) +
  geom_jitter(data = ts_df,
              aes(depth_avg, as.numeric(num_true_positives_pitchtracked > 0)),
              height = 0.05, width = 0, alpha = 0.18, size = 1.1) +
  scale_y_continuous(breaks = c(0,1), labels = c("No","Yes"), limits = c(-0.05, 1.05)) +
  labs(x = "Hourly Average Depth (m)", y = "Pr(Calls > 0)") +
  theme_classic(base_size = 13) + theme(plot.title = element_text(hjust = 0.5))

# Count plot: raw hourly counts WHEN PRESENT (>0) + curve + ribbon
plt_hc_mu <- ggplot() +
  geom_ribbon(data = pred_df,
              aes(depth_avg, ymin = mu_pos_lo, ymax = mu_pos_hi),
              alpha = 0.25) +
  geom_line(data = pred_df, aes(depth_avg, mu_pos), linewidth = 1.1) +
  geom_point(data = dplyr::filter(ts_df, num_true_positives_pitchtracked > 0),
             aes(depth_avg, num_true_positives_pitchtracked),
             alpha = 0.20, size = 1.2) +
  labs(x = "Hourly Average Depth (m)", y = "Calls when present (raw)") +
  theme_classic(base_size = 13)

# Overall plot: raw hourly counts (including zeros) + curve + ribbon
plt_hd_overall <- ggplot() +
  geom_ribbon(data = pred_df,
              aes(depth_avg, ymin = mu_overall_lo, ymax = mu_overall_hi),
              alpha = 0.25) +
  geom_line(data = pred_df, aes(depth_avg, mu_overall), linewidth = 1.1) +
  geom_point(data = ts_df,
             aes(depth_avg, num_true_positives_pitchtracked),
             alpha = 0.15, size = 1.1) +
  labs(x = "Hourly Average Depth (m)", y = "Calls per hour (raw)") +
  theme_classic(base_size = 13)

# Compose the three plots together
hurdle_combo <- plot_grid(plt_hz_p, plt_hc_mu, plt_hd_overall,
                          labels = c("A", "B", "C"), ncol = 3, align = "hv") +
  theme(plot.background = element_rect(fill = "white", color = NA))

print(hurdle_combo)

# Save plot using desired directory, created above (hz_path)
# Keep dpi high (600) for higher quality images
hz_path <- file.path(out_dir_hz,
                     paste0("Hurdle_Raw", gsub(" ", "", species_label), "_", gsub(" ", "", mission_label), "_COMBINED.png"))
ggsave(hz_path, plot = hurdle_combo, width = 18, height = 6, dpi = 600)

######### 9) Effect-size helpers (reporting) #########
# Purpose: Convert model coefficients into interpretable effect sizes.
#   Zero part (logit): exponentiate coefficient to get an odds ratio.
#   Count part (log): exponentiate coefficient to get a multiplicative change; convert to percent.
# Scale: per +100 m depth makes the magnitude easier to interpret than per 1 m.
# Zero part: odds ratio per +100 m depth
b_zero <- suppressWarnings(coef(summary(hurdle_mean))$zero["depth_avg", "Estimate"])
or_zero_100m <- exp(100 * b_zero)
cat("Hurdle zero-part odds ratio per +100 m depth:", sprintf("%.2f", or_zero_100m), "\n")

# Count part: percent change in E[Calls | >0] per +100 m
b_cnt <- suppressWarnings(coef(summary(hurdle_mean))$count["depth_avg", "Estimate"])
pc_cnt_100m <- (exp(100 * b_cnt) - 1) * 100
cat("Hurdle count-part % change in E[Calls|>0] per +100 m:", sprintf("%.1f%%", pc_cnt_100m), "\n")

######### End of script #########
