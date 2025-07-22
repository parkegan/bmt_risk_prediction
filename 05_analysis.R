
# Hospitalization-level analysis for BMT cohort poster/manuscript

# ------------------------------------------------------------------------------

# libraries
library(data.table)
library(tidyverse)
library(tidytable)
library(collapse)
library(janitor)
library(here)
library(pROC)
library(dplyr)
library(tidyr)
library(caret)


# Read in data
#-------------------------------------------------------------------------------

df <- arrow::read_parquet(here("clean", "analysis_250721.parquet"))

# Pull in HCT dates
cohort <- fread(here("clean", "cohort_bmt_250721.csv")) |>
  fselect(mrn, date_of_hct) |>
  funique()

df <-
  left_join(df, cohort, by = "mrn") |>
  ftransform(hct_01 = case_when(t_adt <= date_of_hct & date_of_hct <= t_dc   ~ 1,
                                TRUE                                         ~ 0))

# Create time series df --------------------------------------------------------

ts_df <-
  df |>
  mutate(time_block = floor(as.double(difftime(time, t_ward, units = "hours")/4))) |>
  group_by(mrn, enc, time_block) |>
  arrange(mrn, enc, time_block, time) |>
  filter(row_number() == 1) |>
  ungroup() |>
  fmutate(outcome_24h = case_when(as.double(difftime(t_outcome, time, units = "hours")) <= 24 
                                  & outcome_01 == 1                                                ~ 1,
                                  TRUE                                                             ~ 0)) |>
  funique()

# Examine outcomes -------------------------------------------------------------


cohort_outcomes <-
  df |>
  select(mrn, enc, age_yrs, female_01, auto_01, outcome_01, dead_01, hospice_01, 
         wicutx_01, hct_01, starts_with("min"), starts_with("max")) |>
  funique()

cohort_outcomes |>
  select(
    auto_01, 
    age_yrs,
    female_01,
    outcome_01,
    wicutx_01,
    dead_01,
    hospice_01,
    hct_01
  ) |>
  gtsummary::tbl_summary(by      = outcome_01,
                         missing = "no",
                         type    = list()) |>
  gtsummary::add_p()


icu_ever <-
  cohort_outcomes |>
  filter(wicutx_01 == 1)


# For people who were seen in the ICU, their time of event is the first
# Ward-ICU admission. Censor any SOI metrics after that point.

df <- 
  df |>
  filter(time < t_outcome & 
           time < t_dc & 
           t_ward <= time) 

only_outcomes <- 
  df |>
  filter(outcome_01 == 1)


# See how many people met each threshold ---------------------------------------

inc_df <-
  df |>
  mutate(edi_threshm   = as.factor(case_when(37.4 <= max_edi    ~ 1,
                                             TRUE                ~ 0)),
         edi_threshh   = as.factor(case_when(68.8 <= max_edi    ~ 1,
                                             TRUE               ~ 0)),
         sirs_thresh  = as.factor(case_when(2 <= max_sirs      ~ 1,
                                            TRUE                ~ 0)),
         news_thresh  = as.factor(case_when(7 <= max_news      ~ 1,
                                            TRUE                ~ 0)),
         mews_thresh  = as.factor(case_when(5 <= max_mews      ~ 1,
                                            TRUE                ~ 0)),
         qsofa_thresh = as.factor(case_when(2 <= max_qsofa     ~ 1,
                                            TRUE                ~ 0)))

# How many people never met any threshold?
never_alarm <-
  inc_df |>
  filter(outcome_01 == 1 &
           max_edi < 37.4 &
           max_qsofa < 2 &
           max_mews < 5 &
           max_news < 7) |>
  fselect(mrn, enc, t_adt, t_outcome, wicutx_01, outcome_01, max_edi, 
          max_sirs, max_qsofa, max_mews, max_news) |>
  funique()


# Cumulative incidence of positive scores
inc_df |>
  collapse::fselect(
    enc,
    edi_threshm,
    edi_threshh,
    sirs_thresh,
    news_thresh,
    mews_thresh,
    qsofa_thresh
  ) |>
  collapse::funique() |>
  collapse::fselect(-enc) |>
  gtsummary::tbl_summary() 


# Lead times for positive scores

df_lead <-
  df |>
  filter(outcome_01 == 1) |>
  ftransform(sirs_lead = case_when(sirs >= 2 ~ as.double(difftime(t_outcome, time, units = "days")),
                                   TRUE      ~ NA_real_),
             qsofa_lead = case_when(qsofa >= 2 ~ as.double(difftime(t_outcome, time, units = "days")),
                                    TRUE      ~ NA_real_),
             mews_lead = case_when(mews >= 5 ~ as.double(difftime(t_outcome, time, units = "days")),
                                   TRUE      ~ NA_real_),
             news_lead = case_when(news >= 7 ~ as.double(difftime(t_outcome, time, units = "days")),
                                   TRUE      ~ NA_real_),
             edi_lead = case_when(edi >= 37.4 ~ as.double(difftime(t_outcome, time, units = "days")),
                                  TRUE      ~ NA_real_)) |>
  group_by(mrn, enc) |>
  mutate(sirs_lead = fmax(sirs_lead),
         qsofa_lead = fmax(qsofa_lead),
         mews_lead = fmax(mews_lead),
         news_lead = fmax(news_lead),
         edi_lead = fmax(edi_lead)) |>
  ungroup() |>
  fselect(mrn, enc, sirs_lead, qsofa_lead, mews_lead, news_lead, edi_lead) |>
  funique()

df_lead |>
  collapse::fselect(-mrn, -enc) |>
  gtsummary::tbl_summary() 

# make tables and figures
# --------------------------------------------------------------------------

# table
table_1 <-
  df |>
  collapse::fselect(
    enc,
    age_yrs,
    female_01,
    race_cat,
    auto_01,
    outcome_01,
    wicutx_01
  ) |>
  collapse::funique() |>
  collapse::fselect(-enc) |>
  gtsummary::tbl_summary(by = outcome_01) |>
  gtsummary::add_p()

table_1

# maximum score per hospitalization --------------------------------------------
# grouped by outcome status

c_edi  = "#57B147"
c_news = "#5E97C9"
c_mews = "#FFC939"
c_qsof = "#5D3FD3"
c_sirs = "#C34D36"
c_gray = "#D3D3D3"
scols  = c(c_sirs, c_qsof, c_mews, c_news, c_edi)
sorder = c("SIRS", "qSOFA", "MEWS", "NEWS", "EDI")
olabs  = c("No Deterioration", "Clinical Deterioration")

# Dataframe with maximum scores

df2 <-
  df |>
  collapse::fsubset(!is.na(auto_01)) |>
  collapse::fgroup_by(enc) |>
  collapse::fsummarize(
    auto       = fmax(auto_01),
    qsofa      = fmax(qsofa),
    sirs       = fmax(sirs),
    mews       = fmax(mews),
    news       = fmax(news),
    edi        = fmax(edi),
    outcome    = fmax(outcome_01)
  ) |>
  ungroup() |>
  ftransform(outcome_f = factor(outcome, labels = olabs))

## plot individual histograms
hs <- 
  ggplot(df2, aes(x = sirs, fill = outcome_f)) +
  geom_histogram(binwidth = 1L, color = "black", position = "identity") +
  scale_y_continuous(limits = c(0, 801), breaks = seq(0, 800, 200)) +
  scale_fill_manual(values = c(c_gray, c_sirs)) +
  theme_bw(base_size = 16L) +  
  labs(
    x    = "Highest SIRS",
    y    = "N",
    fill = "Outcome"
  ) +
  theme(legend.position = "none")

hq <- 
  ggplot(df2, aes(x = qsofa, fill = outcome_f)) +
  geom_histogram(binwidth = 1L, color = "black", position = "identity") +
  scale_y_continuous(limits = c(0, 820), breaks = seq(0, 800, 200)) +
  scale_fill_manual(values = c(c_gray, c_qsof)) +
  theme_bw(base_size = 16L) +  
  labs(
    x    = "Highest qSOFA",
    y    = "N",
    fill = "Outcome"
  ) +
  theme(legend.position = "none")

hm <- 
  ggplot(df2, aes(x = mews, fill = outcome_f)) +
  geom_histogram(binwidth = 1L, color = "black", position = "identity") +
  scale_x_continuous(breaks = seq(0L, 10L, 2L)) +
  scale_y_continuous(limits = c(0, 801), breaks = seq(0, 800, 200)) +
  scale_fill_manual(values = c(c_gray, c_mews)) +
  theme_bw(base_size = 16L) +  
  labs(
    x    = "Highest MEWS",
    y    = "N",
    fill = "Outcome"
  ) +
  theme(legend.position = "none")

hn <- 
  ggplot(df2, aes(x = news, fill = outcome_f)) +
  geom_histogram(binwidth = 1L, color = "black", position = "identity") +
  scale_x_continuous(breaks = seq(0L, 1L, 4L)) +
  scale_y_continuous(limits = c(0, 801), breaks = seq(0, 800, 200)) +
  scale_fill_manual(values = c(c_gray, c_news)) +
  theme_bw(base_size = 16L) +  
  labs(
    x    = "Highest NEWS",
    y    = "N",
    fill = "Outcome"
  ) +
  theme(legend.position = "none")

he <- 
  ggplot(df2, aes(x = edi, fill = outcome_f)) +
  geom_histogram(binwidth = 2.5, color = "black", position = "identity") +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  scale_y_continuous(limits = c(0, 125)) +
  scale_fill_manual(values = c(c_gray, c_edi)) +
  theme_bw(base_size = 16L) +  
  labs(
    x    = "Highest EDI",
    y    = "N",
    fill = "Outcome"
  ) +
  theme(legend.position = "bottom")

f1top = hs + hq + hm + hn + plot_layout(nrow = 1, axes = "collect")

ggsave(
  here("figs", paste0("score_hist_top.pdf")), 
  width = 12, height = 6, units = "in"
)

f1bottom = he

ggsave(
  here("figs", paste0("score_hist_bottom.pdf")), 
  width = 12, height = 6, units = "in"
)

f1 = f1top / f1bottom

ggsave(
  here("figs", paste0("f01_hist.pdf")), 
  width = 12, height = 12, units = "in"
)

rm(hs, hq, hm, hn, he)

# time to positivity by each score ---------------------------------------------

# requires 4h discrete time blocks for ease of curve plotting
df4 <-
  df |>
  collapse::fmutate(
    hours   = as.numeric(difftime(time, t_adt, units = "hours")),
    time_4h = 4*floor(hours/4)
  ) |>
  tidytable::select(
    enc,
    auto_01,
    t_ward,
    t_outcome,
    time_4h,
    qsofa,
    sirs,
    mews,
    news,
    edi,
    outcome_01
  ) |>
  # need worst values per time block
  collapse::fgroup_by(enc, time_4h) |>
  collapse::fmutate(across(qsofa:edi, fmax)) |>
  collapse::fungroup() |>
  collapse::funique()

# check to make sure each block has just one row per
#df4 |> get_dupes(enc, time_4h)

#  balance  time series up to 1 week (x axis will go 7 days)

df4 <-
  df4 |>
  group_by(enc, time_4h) |>
  tidytable::mutate(
    block_4_id    = tidytable::cur_group_id()
  )

# make threshold dummies
df4 <-
  df4 |>
  collapse::ftransform(
    thresh_sirs  = fifelse(sirs  >= 2,      as.integer(1), NA_integer_),
    thresh_qsofa = fifelse(qsofa >= 2,      as.integer(1), NA_integer_),
    thresh_mews5 = fifelse(mews  >= 5,      as.integer(1), NA_integer_),
    thresh_news8 = fifelse(news  >= 7,      as.integer(1), NA_integer_),
    # doi:10.1001/jamanetworkopen.2023.24176
    thresh_edim  = fifelse(edi   >= 37.4,   as.integer(1), NA_integer_), 
    thresh_edih  = fifelse(edi   >= 68.8,   as.integer(1), NA_integer_)) |>
  tidytable::arrange(block_4_id) |>
  tidytable::fill(everything(), .direction = "down", .by = enc) |>
  tidytable::mutate(
    tidytable::across(
      .cols = starts_with("thresh"),
      .fns  = ~fifelse(is.na(.x), 0, .x)
    ))

# summarize into proportion-positive by each 4h time block
df5 <-
  df4 |>
  tidytable::group_by(time_4h) |>
  tidytable::summarize(
    total_met_sirs_threshold  = sum(thresh_sirs),
    total_met_qsofa_threshold = sum(thresh_qsofa),
    total_met_mews_threshold  = sum(thresh_mews5),
    total_met_news8_threshold = sum(thresh_news8),
    total_met_edim_threshold  = sum(thresh_edim),
    total_patients_in_block   = tidytable::n()
  ) |>
  collapse::ftransform(
    prop_sirs  = total_met_sirs_threshold   / total_patients_in_block,
    prop_qsofa = total_met_qsofa_threshold  / total_patients_in_block,
    prop_mews  = total_met_mews_threshold   / total_patients_in_block,
    prop_news8 = total_met_news8_threshold  / total_patients_in_block,
    prop_edim =  total_met_edim_threshold   / total_patients_in_block,
    days       = time_4h/24
  ) |>
  tidytable::select(-starts_with("total"))

f2 <-
  df5 |>
  tidytable::pivot_longer(
    cols         = starts_with("prop"),
    names_to     = "score",
    values_to    = "prop",
    names_prefix = "prop_"
  ) |>
  collapse::fmutate(
    score = factor(
      score,
      levels = c("sirs", "qsofa", "mews", "news8", "edim"),
      labels = c("SIRS (>=2)", "qSOFA (>= 2)", "MEWS (>= 5)", "NEWS (>= 7)", "EDI (>= 37.4)")
    )
  ) |>
  ggplot(aes(x = days, y = prop, group = score)) +
  geom_smooth(aes(color = score), se = F, linewidth = 1.7, alpha = 0.9) +
  scale_x_continuous(
    limits = c(0, 7),
    breaks = c(0:7),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, .6),
    breaks = seq(0, .6, 0.1),
    expand = c(0.01, 0.01)
  ) +
  scale_color_manual(values = score_colors) +
  labs(
    x     = "Days in Hospital",
    y     = "Cumulative Incidence",
    color = "Score",
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position       = c(0.2, 0.8),
    legend.box.background = element_rect(),
    legend.box.margin     = margin(6, 6, 6, 6)
  )

f2

ggsave(here("figs", "figure_2.jpg"), height = 8, width = 6, units = "in")

# --------------------------------------------------------------------------

# combine plot panels

library(patchwork)

f1 + f2

ggsave(here("figs", "figure_combined.jpg"), height = 10, width = 6, units = "in")


# Visual abstract figures ------------------------------------------------------

q = ci.auc(roc(df2$outcome, df2$qsofa), method = "bootstrap")
s = ci.auc(roc(df2$outcome, df2$sirs),  method = "bootstrap")
m = ci.auc(roc(df2$outcome, df2$mews),  method = "bootstrap")
n = ci.auc(roc(df2$outcome, df2$news),  method = "bootstrap")
e = ci.auc(roc(df2$outcome, df2$edi),   method = "bootstrap")

f2a = 
  tidytable(
    score = sorder,
    auroc = c(s[[2]], q[[2]], m[[2]], n[[2]], e[[2]]),
    ci_lo = c(s[[1]], q[[1]], m[[1]], n[[1]], e[[1]]),
    ci_hi = c(s[[3]], q[[3]], m[[3]], n[[3]], e[[3]])
  ) |> 
  ftransform(score = factor(score, levels = sorder)) |>
  ggplot(aes(x = auroc, y = score, group = score)) +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = 1.2) +
  geom_point(aes(color = score), size = 4) +
  scale_x_continuous(
    limits = c(0.5, 0.86),
    breaks = seq(0.5, 0.85, 0.05)
  ) +
  scale_color_manual(values = scols) +
  labs(
    x     = "AUROC for Clinical Deterioration",
    y     = "",
    title = "A. Hospitalization-Level Discrimination"
  ) +
  guides(color = "none") +
  theme(aspect.ratio = 0.4)

f2a

ts <-
  probably::threshold_perf(
    .data       = df2,
    truth       = outcome_f,
    estimate    = sirs,
    thresholds  = seq(0, 4, 1),
    na_rm       = TRUE,
    event_level = "second"
  ) |>
  ftransform(score = "sirs")

tq <-
  probably::threshold_perf(
    .data       = df2,
    truth       = outcome_f,
    estimate    = qsofa,
    thresholds  = seq(0, 3, 1),
    na_rm       = TRUE,
    event_level = "second"
  ) |>
  ftransform(score = "qsofa")

tm <-
  probably::threshold_perf(
    .data       = df2,
    truth       = outcome_f,
    estimate    = mews,
    thresholds  = seq(0, 10, 1),
    na_rm       = TRUE,
    event_level = "second"
  ) |>
  ftransform(score = "mews")

tn <-
  probably::threshold_perf(
    .data       = df2,
    truth       = outcome_f,
    estimate    = news,
    thresholds  = seq(0, 17, 1),
    na_rm       = TRUE,
    event_level = "second"
  ) |>
  ftransform(score = "news")

te <-
  probably::threshold_perf(
    .data       = df2,
    truth       = outcome_f,
    estimate    = edi,
    thresholds  = seq(14, 85, 0.1),
    na_rm       = TRUE,
    event_level = "second"
  ) |>
  ftransform(score = "edi")

thresh <-
  ts |>
  bind_rows(tq) |>
  bind_rows(tm) |>
  bind_rows(tn) |>
  bind_rows(te) |>
  fmutate(
    .estimate = signif(.estimate, 4),
    .estimate = case_when(
      .estimate > 0.999 ~ 0.999,
      .estimate < 0.001 ~ 0.001,
      TRUE              ~ .estimate
    )
  ) |>
  pivot_wider(
    names_from  = .metric,
    values_from = .estimate
  ) |>
  select(score, .threshold, sens = sensitivity) 

f2b =  
  df2 |>
  select(-starts_with("outcome")) |>
  pivot_longer(-enc, names_to = "score", values_to = "maxval") |>
  join(thresh, how = "left", multiple = T) |>
  ftransform(
    pos_01 = if_else(maxval >= .threshold, 1L, 0L),
    one    = 1L
  ) |>
  fgroup_by(score, sens) |>
  fsummarize(
    pos = fsum(pos_01),
    n   = fsum(one)
  ) |>
  ftransform(
    pct   = pos/n,
    score = factor(
      score, 
      levels = c("sirs", "qsofa", "mews", "news", "edi"), 
      labels = c("SIRS", "qSOFA", "MEWS", "NEWS", "EDI")
    )
  ) |>
  ggplot(aes(x = sens, y = pct, color = score)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "gray") +
  geom_line(linewidth = 1.7, alpha = 0.9) +
  labs(
    x     = "Sensitivity",
    y     = "Percent of Patients Crossing Threshold",
    color = "Score",
    title = "B. Hospitalization-Level Efficiency"
  ) +
  scale_color_manual(values = scols) +
  guides(color = "none") +
  theme(aspect.ratio = 1)

# 80% sensitivity thresholds
dthresh <-
  thresh |>
  fsubset(sens >= 0.8) |>
  fgroup_by(score) |>
  fsummarize(threshold = fmax(.threshold))

rm(ts, tq, tm, tn, te, thresh)

figure2 = f2a / f2b

figure2 + plot_layout(heights = c(1, 2)) & theme_bw(base_size = 14)

ggsave(
  here("figs", paste0("f02_square_", today, ".pdf")), 
  width = 8, height = 8, units = "in"
)


# Model performance comparisons ------------------------------------------------

# Build hospitalization-level and time series models for each predictor

analyze_performance <- function(data, 
                                predictor_name, 
                                outcome_name, 
                                conf_level = 0.95, 
                                boot_n = 2000) {
  
  # Check if required packages are installed
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' needed for this function to work. Please install it.")
  }
  
  prevalence <- sum(data[[outcome_name]]) / nrow(data)
  
  formula <- as.formula(paste(outcome_name, "~", predictor_name))
  
  model_fit <- glm(formula, data = data, family = binomial)
  
  model_summary <- summary(model_fit)
  
  # Generate predictions
  pred_values <- predict(model_fit, newdata = data, type = "response")
  
  # ROC and AUC
  roc_obj <- pROC::roc(data[[outcome_name]], pred_values)
  auc_value <- pROC::auc(roc_obj)
  
  # AUC CIs
  auc_ci <- pROC::ci.auc(roc_obj, conf.level = conf_level,
                         method = c("delong", "bootstrap"),
                         boot.n = boot_n, boot.stratified = TRUE,
                         reuse.auc = TRUE,
                         progress = getOption("pROCProgress")$name,
                         parallel = FALSE)
  
  list(
    prevalence = prevalence,
    model_summary = model_summary,
    roc = roc_obj,
    auc = auc_value,
    auc_ci = auc_ci,
    predictor = predictor_name
  )
}

predictors <- c("edi", "sirs", "qsofa", "mews", "news")

hospital_level <- lapply(predictors, 
                         function(p) analyze_performance(df2, p, "outcome"))

time_series <- lapply(predictors, 
                      function(p) analyze_performance(ts_df, p, "outcome_24h")) 

names(hospital_level) <- predictors
names(time_series) <- predictors

# Compare performance of predictors between hospitalization-level and 
# time series models 

compare_outcomes <- 
  function(hospital_results, timeseries_results) {
    # Check if pROC is available
    if (!requireNamespace("pROC", quietly = TRUE)) {
      stop("Package 'pROC' needed for ROC comparisons. Please install it.")
    }
    
    comparisons <- list()
    
    for (pred in names(hospital_results)) {
      roc_hospital <- hospital_results[[pred]]$roc
      roc_timeseries <- timeseries_results[[pred]]$roc
      
      roc_comparison <- tryCatch(
        {
          pROC::roc.test(
            roc1 = roc_hospital,
            roc2 = roc_timeseries,
            method = "delong",  
            paired = FALSE      
          )
        },
        error = function(e) {
          warning(paste("ROC comparison failed for", pred, ":", e$message))
          return(NULL)
        }
      )
      
      # Store all comparison metrics
      comparisons[[pred]] <- list(
        hospital_auc = as.numeric(hospital_results[[pred]]$auc),
        timeseries_auc = as.numeric(timeseries_results[[pred]]$auc),
        auc_diff = as.numeric(timeseries_results[[pred]]$auc - hospital_results[[pred]]$auc),
        roc_test = roc_comparison,
        hospital_prevalence = hospital_results[[pred]]$prevalence,
        timeseries_prevalence = timeseries_results[[pred]]$prevalence,
        hospital_auc_ci = hospital_results[[pred]]$auc_ci,
        timeseries_auc_ci = timeseries_results[[pred]]$auc_ci
      )
    }
    
    return(comparisons)
  }

outcome_comparisons <- compare_outcomes(hospital_level, time_series)

# Compare performance between auto/allo groups ---------------------------------
#  (hospitalization level)

compare_groups <- function(data, 
                           predictor_names, 
                           outcome_name = "outcome", 
                           auto_var = "auto", 
                           conf_level = 0.95, 
                           boot_n = 2000) {
  
  # Check if required packages are installed
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' needed for this function to work. Please install it.")
  }
  
  results <- list()
  
  data_auto0 <- data[data[[auto_var]] == 0, ]
  data_auto1 <- data[data[[auto_var]] == 1, ]
  
  analyze_predictor <- function(d, pred, out) {
    prevalence <- sum(d[[out]]) / nrow(d)
    
    formula <- as.formula(paste(out, "~", pred))
    
    model_fit <- glm(formula, data = d, family = binomial)
    model_summary <- summary(model_fit)
    
    pred_values <- predict(model_fit, newdata = d, type = "response")
    
    # ROC and AUC
    roc_obj <- pROC::roc(d[[out]], pred_values)
    auc_value <- pROC::auc(roc_obj)
    
    # AUC CI
    auc_ci <- pROC::ci.auc(roc_obj, conf.level = conf_level,
                           method = c("delong", "bootstrap"),
                           boot.n = boot_n, boot.stratified = TRUE,
                           reuse.auc = TRUE,
                           progress = getOption("pROCProgress")$name,
                           parallel = FALSE)
    
    list(
      prevalence = prevalence,
      model_summary = model_summary,
      roc = roc_obj,
      auc = auc_value,
      auc_ci = auc_ci
    )
  }
  
  for (predictor in predictor_names) {
    auto0_results <- analyze_predictor(data_auto0, predictor, outcome_name)
    
    auto1_results <- analyze_predictor(data_auto1, predictor, outcome_name)
    
    roc_test <- pROC::roc.test(
      auto0_results$roc,
      auto1_results$roc,
      method = "delong",
      paired = FALSE
    )
    
    results[[predictor]] <- list(
      auto0 = auto0_results,
      auto1 = auto1_results,
      comparison = list(
        roc_test = roc_test,
        auc_diff = as.numeric(auto1_results$auc - auto0_results$auc),
        prevalence_diff = auto1_results$prevalence - auto0_results$prevalence
      )
    )
  }
  
  return(results)
}

comparison_results <- 
  compare_groups(
    data = df2,
    predictor_names = predictors,
    outcome_name = "outcome",
    auto_var = "auto"
  )

summary_table <- 
  do.call(rbind, lapply(names(comparison_results), function(pred) {
    data.frame(
      Predictor = pred,
      Auto0_AUC = as.numeric(comparison_results[[pred]]$auto0$auc),
      Auto1_AUC = as.numeric(comparison_results[[pred]]$auto1$auc),
      AUC_Difference = comparison_results[[pred]]$comparison$auc_diff,
      P_Value = comparison_results[[pred]]$comparison$roc_test$p.value,
      Auto0_Prevalence = comparison_results[[pred]]$auto0$prevalence,
      Auto1_Prevalence = comparison_results[[pred]]$auto1$prevalence
    )
  }))

print(summary_table)


# End -------------------------- -----------------------------------------------