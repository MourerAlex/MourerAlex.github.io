# ============================================================
# part5_run.R
# Executable companion to Part 5: Implementation — From Formula to Code
#
# Separable effects for the DES prostate cancer trial
# G-formula and IPW estimators (Stensrud et al. 2020)
# ============================================================


# ── Section 0: Libraries + Data Loading ──────────────────────

library(Hmisc)      # getHdata() for prostate dataset
library(survival)   # survSplit for person-time expansion
library(dplyr)      # readability (case_when, mutate, group_by)

a_y <- 1   # A_Y intervention: DES testosterone suppression
a_d <- 0   # A_D intervention: placebo cardiovascular rate


# ── Load and preprocess ──

getHdata(prostate)

# Strip Hmisc labels (they conflict with dplyr mutate)
prostate <- as.data.frame(prostate)
prostate$rx     <- as.character(prostate$rx)
prostate$status <- as.character(prostate$status)
prostate$pf     <- as.character(prostate$pf)

subject_df <- prostate %>%
  filter(rx %in% c("placebo", "5.0 mg estrogen")) %>%
  mutate(
    A = as.integer(rx == "5.0 mg estrogen"),              # DES = 1, placebo = 0
    event_type = case_when(                                # 0=censored, 1=Y(prostate), 2=D(other)
      status == "dead - prostatic ca" ~ 1L,
      status == "alive"               ~ 0L,
      TRUE                            ~ 2L),
    tstart     = 0,                                        # month 0 = enrollment for everyone
    event_time = dtime + 1,                                # shift +1: "died month 0" → month 1
    normal_act = as.integer(pf == "normal activity"),
    age_cat    = cut2(age, c(0, 60, 70, 80, 100)),
    cv_hist    = hx,
    hemo_bin   = as.integer(hg < 12)
  ) %>%
  select(id = patno, A, event_time, event_type, tstart,
         normal_act, age_cat, cv_hist, hemo_bin)

cat("═══ Section 0: Data Summary ═══\n")
cat("n =", nrow(subject_df), "\n")
cat("Treatment: DES =", sum(subject_df$A == 1),
    ", Placebo =", sum(subject_df$A == 0), "\n")
cat("Events: Y =", sum(subject_df$event_type == 1),
    ", D =", sum(subject_df$event_type == 2),
    ", Censored =", sum(subject_df$event_type == 0), "\n")
cat("Event time range:", range(subject_df$event_time), "\n\n")


# ── Section 1: Person-Time Expansion ─────────────────────────

# event_indicator = 1L for ALL subjects, including censored.
# survSplit uses this to flag each subject's terminal row (= 1 on last row,
# 0 on all earlier rows). If censored subjects had 0, their last row
# wouldn't be flagged and we couldn't derive c_event = 1 from it.
subject_df <- subject_df %>%
  mutate(event_indicator = 1L)

cut_times <- 1:75   # event_time range is 1–76 after +1 shift

# Expand: one row per subject per interval [tstart, tstop)
pt_df <- survSplit(
  data  = subject_df,
  cut   = cut_times,
  start = "tstart",
  end   = "event_time",
  event = "event_indicator"
)

pt_df <- pt_df %>%
  rename(tstop = event_time) %>%   # after survSplit these are interval boundaries
  mutate(k = tstart)               # k = interval index (month 0, 1, 2, ...)

# Derive event indicators from event_indicator + event_type.
# event_indicator == 1 only on the terminal row; event_type says which event.
pt_df <- pt_df %>%
  mutate(
    y_event = as.integer(event_indicator == 1 & event_type == 1),  # prostate death
    d_event = as.integer(event_indicator == 1 & event_type == 2),  # other death
    c_event = as.integer(event_indicator == 1 & event_type == 0)   # censored/alive
  )

# Temporal ordering (C_k, D_k, Y_k) — set unobserved outcomes to NA.
# glm drops NA rows automatically, defining the correct risk sets.
pt_df <- pt_df %>%
  mutate(
    d_event = ifelse(c_event == 1, NA, d_event),     # censored → D unobserved
    y_event = ifelse(c_event == 1, NA, y_event),     # censored → Y unobserved
    y_event = ifelse(d_event == 1, NA, y_event),     # D occurred → Y unobserved (sequential in mutate)
    A_y = A,                                          # Y-hazard model sees A_y
    A_d = A                                           # D-hazard model sees A_d
  )

baseline_df <- pt_df %>% filter(k == 0)   # one row per subject
n <- n_distinct(pt_df$id)

cat("═══ Section 1: Person-Time Data ═══\n")
cat("Total rows:", nrow(pt_df), ", n =", n, ", k range:", range(pt_df$k), "\n")
cat("y_event:"); print(table(pt_df$y_event, useNA = "always"))
cat("d_event:"); print(table(pt_df$d_event, useNA = "always"))
cat("c_event:"); print(table(pt_df$c_event, useNA = "always"))
cat("\n")


# ── Section 2: G-formula ─────────────────────────────────────
#
# Target: P(Y^{a_Y, a_D, c̄=0}_{K+1} = 1)
#
# Identification (Stensrud 2020, eq. 8/9): the counterfactual cumulative
# incidence equals a function of observable cause-specific hazards, with
# Y-hazard evaluated at A=a_Y and D-hazard at A=a_D (the cross-arm trick).
#
# Model specification: additive treatment (A_y + k + ...) — constant effect
# over time. The workshop uses interactions (rx * (dtime + ...)) allowing
# time-varying effects. Both valid; results differ by modelling choice.

# ── 2.1 Fit hazard models ──

fit_y_haz <- glm(
  y_event ~ A_y + k + I(k^2) + I(k^3) + normal_act + age_cat + cv_hist + hemo_bin,
  data = pt_df, family = binomial()
)

fit_d_haz <- glm(
  d_event ~ A_d + k + I(k^2) + I(k^3) + normal_act + age_cat + cv_hist + hemo_bin,
  data = pt_df, family = binomial()
)

cat("═══ Section 2: G-formula ═══\n")
cat("Y-hazard coefficients:\n"); print(round(coef(fit_y_haz), 4))
cat("\nD-hazard coefficients:\n"); print(round(coef(fit_d_haz), 4))

# ── 2.2 Create cloned datasets ──
# For each (a_Y, a_D), expand baseline over all time points.
# Every subject "lives" through all times — the c̄=0 world.
# The standardization Σ_l f(l)P(L=l) becomes (1/n) Σ_i f(l_i):
# each subject carries their own covariates, the average IS the standardization.

make_clone <- function(baseline_df, a_y, a_d, times) {
  n_subj <- nrow(baseline_df)
  K      <- length(times)
  clone       <- baseline_df[rep(seq_len(n_subj), each = K), ]
  clone$k     <- rep(times, n_subj)
  clone$A_y   <- a_y
  clone$A_d   <- a_d
  return(clone)
}

time_points <- c(0, cut_times)   # all k values: 0, 1, ..., 75

clone_11 <- make_clone(baseline_df, a_y = 1, a_d = 1, time_points)  # DES
clone_00 <- make_clone(baseline_df, a_y = 0, a_d = 0, time_points)  # Placebo
clone_10 <- make_clone(baseline_df, a_y = 1, a_d = 0, time_points)  # Modified

# ── 2.3 Predict hazards (cross-arm trick) ──
# fit_y_haz sees A_y, fit_d_haz sees A_d.
# In clone_10: Y-model gets A_y=1 (DES), D-model gets A_d=0 (placebo).

predict_hazards <- function(clone, fit_y, fit_d) {
  clone$haz_y <- predict(fit_y, newdata = clone, type = "response")
  clone$haz_d <- predict(fit_d, newdata = clone, type = "response")
  # Joint event-free probability within ONE interval (not cumulative).
  # Assumes Y and D conditionally independent within each interval.
  clone$surv  <- (1 - clone$haz_y) * (1 - clone$haz_d)
  return(clone)
}

clone_11 <- predict_hazards(clone_11, fit_y_haz, fit_d_haz)
clone_00 <- predict_hazards(clone_00, fit_y_haz, fit_d_haz)
clone_10 <- predict_hazards(clone_10, fit_y_haz, fit_d_haz)

# ── 2.4 Compute cumulative incidence ──
# Parametric g-formula (Stensrud 2020 eq. 8/9).
# Must compute per-subject first, then average — cumprod doesn't commute
# with mean (unlike the nonparametric Aalen-Johansen which works at the
# population level and could average first).
compute_cum_inc <- function(clone_df) {
  clone_df %>%
    group_by(id) %>%
    arrange(k) %>%
    mutate(
      cum_surv     = cumprod(surv),               # S_i(s): cumulative event-free survival
      lag_cum_surv = lag(cum_surv, default = 1),   # S_i(s-1); default=1 avoids workshop time-0 bug
      inc          = haz_y * (1 - haz_d) * lag_cum_surv   # subdensity increment
    ) %>%
    ungroup() %>%
    # Average over subjects at each k, then cumulate over time.
    # Equivalent to cumsum-per-subject-then-mean (Σ_s Σ_i = Σ_i Σ_s),
    # but collapses n×K → K rows first — more memory-efficient.
    group_by(k) %>%
    summarise(mean_inc = mean(inc), .groups = "drop") %>%
    arrange(k) %>%
    mutate(cum_inc = cumsum(mean_inc))
}

cum_inc_11 <- compute_cum_inc(clone_11)   # ν(1,1): DES
cum_inc_00 <- compute_cum_inc(clone_00)   # ν(0,0): Placebo
cum_inc_10 <- compute_cum_inc(clone_10)   # ν(1,0): Modified

K_end <- max(cum_inc_11$k)
cat("\n── G-formula cumulative incidence at month", K_end, "──\n")
cat("  ν(1,1) DES:     ", round(tail(cum_inc_11$cum_inc, 1), 4), "\n")
cat("  ν(0,0) Placebo: ", round(tail(cum_inc_00$cum_inc, 1), 4), "\n")
cat("  ν(1,0) Modified:", round(tail(cum_inc_10$cum_inc, 1), 4), "\n\n")


# ── Section 3: IPW Estimator ─────────────────────────────────
#
# IPW Representation 1 (Part 4a Section 5.6):
# Stand in the A=a_Y arm. Y hazard and Y event-free are already correct.
# W_D swaps the D distribution from a_Y to a_D. W_C corrects for censoring.
# No Y hazard model needed — only D hazard (for W_D) and censoring (for W_C).

# ── 3.1 Fit censoring model ──
# Uses observed A (not overridden a_Y): W_C corrects for censoring as it
# actually operated. In the a_Y arm, A = a_y already.

fit_cens <- glm(
  c_event ~ A + k + I(k^2) + I(k^3) + normal_act + age_cat + cv_hist + hemo_bin,
  data = pt_df, family = binomial()
)

cat("═══ Section 3: IPW ═══\n")
cat("Censoring model coefficients:\n"); print(round(coef(fit_cens), 4))

# ── 3.2 Predict cause-specific D-free probability at both arms ──
# predict(type="response") returns the hazard P(D=1|...).
# d_free_k = 1 - hazard = P(D=0|...) at this ONE interval (not cumulative).
pt_df <- pt_df %>%
  mutate(
    d_free_k_ad  = 1 - predict(fit_d_haz,                       # rate we WANT (numerator of W_D)
                                newdata = pt_df %>% mutate(A_d = a_d),
                                type = "response"),
    d_free_k_ay  = 1 - predict(fit_d_haz,                       # rate to CANCEL (denominator of W_D)
                                newdata = pt_df %>% mutate(A_d = a_y),
                                type = "response"),
    uncens_k     = 1 - predict(fit_cens, newdata = pt_df,        # P(not censored) this interval
                                type = "response")
  )

# ── 3.3 Compute weights ──
# W_D(s) = cs_surv_d_ad / cs_surv_d_ay  (Stensrud 2020, eq. 11)
#   cs_surv = cause-specific survival = cumprod of interval D-free probabilities.
#   Conditional on being event-free (Y-free AND D-free) at each prior interval
#   (the conditioning on Y-free is implicit — the hazard model is fitted on
#   the risk set that excludes Y events).
#   W_D > 1 when D is more likely under a_Y than a_D: surviving D in the a_Y
#   arm is "too rare", so we upweight those who made it.
#
# W_C(s) = 1 / cs_surv_cens  (standard IPCW: upweight uncensored subjects)
pt_df <- pt_df %>%
  group_by(id) %>%               # cumprod MUST be per-subject
  arrange(k) %>%
  mutate(
    cs_surv_d_ad = cumprod(d_free_k_ad),    # cumulative D-free under a_D
    cs_surv_d_ay = cumprod(d_free_k_ay),    # cumulative D-free under a_Y
    cs_surv_cens = cumprod(uncens_k),        # cumulative uncensored
    w_d     = cs_surv_d_ad / cs_surv_d_ay,   # swap D distribution
    w_cens  = 1 / cs_surv_cens,              # correct for censoring
    w_total = w_d * w_cens                   # combined weight
  ) %>%
  ungroup()

cat("\nWeight summary (A = a_y arm only):\n")
ay_weights <- pt_df %>% filter(A == a_y)
cat("  W_D:     "); print(summary(ay_weights$w_d))
cat("  W_C:     "); print(summary(ay_weights$w_cens))
cat("  W_total: "); print(summary(ay_weights$w_total))

# ── 3.4 Weighted estimator ──
# At each time s in the a_Y arm:
#   1. Risk set = uncensored, D-free, Y-observable subjects
#   2. Among those, find Y events
#   3. Sum their weights / n at BASELINE (not current risk set — weights
#      already correct for attrition, this is Horvitz-Thompson)
#   4. Cumulate over time → cumulative incidence
#
# Key difference from g-formula: no Y model, no clones, no per-subject
# prediction. Just observed Y events in one arm, reweighted. Trade-off:
# doesn't need Y model correct, but needs D and censoring models correct.
estimate_sep_eff <- function(pt_df, a_y_val, max_k) {
  ay_df <- pt_df %>% filter(A == a_y_val)
  n_ay  <- sum(ay_df$k == 0)           # baseline count, NOT current risk set
  times <- sort(unique(ay_df$k))
  times <- times[times > 0 & times <= max_k]

  event_vec <- numeric(length(times))
  for (i in seq_along(times)) {
    s <- times[i]
    riskset_s  <- ay_df %>%             # alive, uncensored, D-free at time s
      filter(k == s, !is.na(y_event), !is.na(d_event),
             c_event == 0, d_event == 0)
    eventset_s <- riskset_s %>%         # Y events at time s
      filter(y_event == 1)
    event_vec[i] <- sum(eventset_s$w_total) / n_ay   # weighted incidence increment
  }
  data.frame(k = times, cum_inc = cumsum(event_vec))
}

max_k <- max(pt_df$k)

# Observed arms (a_Y = a_D): W_D = 1, only censoring weights needed
cum_inc_11_ipw <- estimate_sep_eff(
  pt_df %>% mutate(w_total = w_cens), a_y_val = 1, max_k = max_k
)
cum_inc_00_ipw <- estimate_sep_eff(
  pt_df %>% mutate(w_total = w_cens), a_y_val = 0, max_k = max_k
)

# Cross-arm (a_Y=1, a_D=0): full W_D × W_C
cum_inc_10_ipw <- estimate_sep_eff(pt_df, a_y_val = 1, max_k = max_k)

cat("\n── IPW cumulative incidence at month", max_k, "──\n")
cat("  ν(1,1) DES:     ", round(tail(cum_inc_11_ipw$cum_inc, 1), 4), "\n")
cat("  ν(0,0) Placebo: ", round(tail(cum_inc_00_ipw$cum_inc, 1), 4), "\n")
cat("  ν(1,0) Modified:", round(tail(cum_inc_10_ipw$cum_inc, 1), 4), "\n\n")


# ── Section 4: Causal Contrasts ──────────────────────────────

gf_11 <- tail(cum_inc_11$cum_inc, 1)
gf_00 <- tail(cum_inc_00$cum_inc, 1)
gf_10 <- tail(cum_inc_10$cum_inc, 1)

ipw_11 <- tail(cum_inc_11_ipw$cum_inc, 1)
ipw_00 <- tail(cum_inc_00_ipw$cum_inc, 1)
ipw_10 <- tail(cum_inc_10_ipw$cum_inc, 1)

cat("═══ Section 4: Causal Contrasts ═══\n\n")

cat("── Cumulative incidence at month", K_end, "──\n")
print(data.frame(
  arm      = c("(1,1) DES", "(0,0) Placebo", "(1,0) Modified"),
  gformula = round(c(gf_11, gf_00, gf_10), 4),
  ipw      = round(c(ipw_11, ipw_00, ipw_10), 4)
), row.names = FALSE)

cat("\n── Risk differences ──\n")
print(data.frame(
  estimand = c("Total effect", "Sep. indirect (A_D)", "Sep. direct (A_Y)"),
  gformula = round(c(gf_11 - gf_00, gf_11 - gf_10,
                      (gf_11 - gf_00) - (gf_11 - gf_10)), 4),
  ipw      = round(c(ipw_11 - ipw_00, ipw_11 - ipw_10,
                      (ipw_11 - ipw_00) - (ipw_11 - ipw_10)), 4)
), row.names = FALSE)

cat("\nG-formula and IPW should roughly agree if models are correctly specified.\n")
cat("Discrepancies suggest model misspecification.\n")
cat("\nDone.\n")
