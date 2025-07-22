
# Data prep for poster/manuscript

# setup ------------------------------------------------------------------------

## libraries
library(tidytable)
library(collapse)
library(arrow)
library(here)
library(pROC)

## helpers
today = format(Sys.Date(), "%y%m%d")

# data -------------------------------------------------------------------------

df  = fread(here("clean/cohort_bmt.csv"))
adt = read_parquet(here("clean/adt.parquet")) 
edi = read_parquet(here("clean", "epic_scores.parquet"))
vs  = read_parquet(here("../000_data/ohsu/clean/vs_dynam_ohsu_240126.parquet"))
ms  = read_parquet(here("../000_data/ohsu/clean/vs_mentalstatus_ohsu_240218.parquet"))
ox  = read_parquet(here("../000_data/ohsu/clean/fs_resp_ohsu_240216.parquet"))
lab = read_parquet(here("../000_data/ohsu/intermediate/lab_ohsu.parquet"))

# outcomes ---------------------------------------------------------------------

## remove short (i.e., fake) ICU stays
## pgl reviewed on 24-05-25 and stays < 8h are for procedures consistently

short <-
  read_parquet(here("clean", "short_icu_stays_240425.parquet")) |>
  fsubset(time_diff_hours < 8) |> 
  ftransform(
    day     = lubridate::date(time_start),
    drop    = 1L,
    loc_cat = "ICU"
  ) |>
  select(mrn, enc, day, drop, loc_cat)

adt <- 
  adt |>
  ftransform(day  = lubridate::date(time)) |>
  join(short, how = "left") |>
  fsubset(is.na(drop)) |>
  select(-drop, -day) |>
  fsubset(loc_cat %in% c("ED", "ICU", "Wards")) |>
  join(df |> fselect(enc), how = "inner", multiple = T) |>
  roworder(time) |>
  mutate(
    w_icu_tx = if_else(loc_cat == "ICU" & lag(loc_cat == "Wards"), 1L, 0L),
    .by      = enc
  ) 

t_ward = 
  fsubset(adt, loc_cat == "Wards") |> 
  fgroup_by(enc) |> 
  fsummarize(
    t_adt  = fmin(admit_dttm_adt),
    t_dc   = fmax(disch_dttm_adt),
    t_ward = fmin(time)
  )

t_icutx = 
  fsubset(adt, w_icu_tx == 1) |>
  fgroup_by(enc) |>
  fsummarize(t_icutx = fmin(time))

df <-
  df |>
  ftransform(
    auto_01    = if_else(donor_type   == "Auto", 1L, 0L),
    outcome_01 = if_else(dischg_dispo %in% c("dead", "hospice"), 1L, 0L),
    dead_01    = if_else(dischg_dispo == "dead",    1L, 0L),
    hospice_01 = if_else(dischg_dispo == "hospice", 1L, 0L)
  ) |>
  select(mrn, enc, age_yrs, ends_with("cat"), ends_with("1"), ends_with("tm")) |>
  join(t_ward,  how = "left") |>
  join(t_icutx, how = "left") |>
  ftransform(
    outcome_01 = if_else(is.na(t_icutx), outcome_01, 1L),
    wicutx_01  = if_else(is.na(t_icutx), 0L, 1L),
    t_outcome  = pmin(t_dc, t_icutx, na.rm = T)
  ) |>
  relocate(starts_with("t_"), .before = outcome_01) |>
  select(-t_icutx, -ends_with("tm")) |>
  fsubset(!is.na(t_ward))

rm(adt, t_ward, t_icutx)

# score parameters -------------------------------------------------------------

vs <-
  vs |>
  join(
    df |> fselect(mrn, enc, t_adt, t_dc), 
    how      = "inner",
    multiple = T
  ) |>
  fsubset(time >= t_adt & time <= t_dc) |>
  ftransform(
    spo2 = if_else(spo2 < 60, NA_integer_, spo2),
    sbp  = if_else(is.na(sbp_ni), sbp_art, sbp_ni),
    tmp  = if_else(is.na(temp_f), 32 + (9*temp_c/5), temp_f)
  ) |>
  select(enc, time, hr, rr, sbp, spo2, tmp)

ms <-
  ms |>
  join(
    df |> fselect(mrn, enc, t_adt, t_dc), 
    how      = "inner",
    multiple = T
  ) |>
  fsubset(time >= t_adt & time <= t_dc) |>
  ftransform(
    gcs = case_when(
      !is.na(gcs) ~ as.integer(gcs),
      avpu_cat == "alert" ~ 15L,
      avpu_cat == "voice" ~ 13L,
      avpu_cat == "pain"  ~ 08L,
      !is.na(avpu_cat)    ~ 03L,
      TRUE                ~ NA_integer_
    )
  ) |>
  fselect(enc, time, gcs)

ox <-
  ox |>
  fselect(-enc) |>
  join(
    df |> fselect(mrn, enc, t_adt, t_dc), 
    how      = "inner",
    multiple = T
  ) |>
  fsubset(time >= t_adt & time <= t_dc) |>
  fsubset(device_name != "room air") |>
  fsubset(o2_flow > 0) |>
  ftransform(o2_01 = 1L) |>
  fselect(enc, time, o2_01)

lab <-
  lab |>
  fselect(-enc, -order_dttm) |>
  join(
    df |> fselect(mrn, enc, t_adt, t_dc), 
    how      = "inner",
    multiple = T
  ) |>
  ftransform(t_adt = t_adt - lubridate::dhours(12)) |>
  fsubset(result_dttm >= t_adt & result_dttm <= t_dc) |>
  ftransform(
    test = case_when(
      common_nm == "WHITE CELL COUNT" ~ "wbc",
      common_nm == "PH, VENOUS"       ~ "ph_ven",
      common_nm == "PH, ARTERIAL"     ~ "ph_art",
      TRUE                            ~ NA_character_
    )
  ) |>
  fsubset(!is.na(test)) |>
  select(enc, time = result_dttm, test, result) |>
  funique() |>
  pivot_wider(
    names_from  = test,
    values_from = result
  ) |>
  fmutate(
    wbc = stringr::str_remove(wbc, "<"),
    wbc = as.numeric(wbc),
    ph_ven = as.numeric(ph_ven),
    ph_art = as.numeric(ph_art)
  )

soi <-
  vs |> 
  join(ms,  how = "full", multiple = T) |>
  join(lab, how = "full", multiple = T) |>
  select(-starts_with("ph")) |>
  roworder(time) |>
  fill(everything(), .direction = "down", .by = enc) |>
  join(ox,  how = "full", multiple = T) |>
  ftransform(time = lubridate::ceiling_date(time, "20 minutes")) |>
  fgroup_by(enc, time) |>
  fsummarize(
    hr    = fmax(hr),
    rr    = fmax(rr),
    sbp   = fmin(sbp),
    spo2  = fmin(spo2),
    tmp   = fmax(tmp),
    gcs   = fmin(gcs),
    wbc   = fmin(wbc),
    o2_01 = fmax(o2_01)
  )

rm(vs, ms, ox)

# SIRS -------------------------------------------------------------------------

soi <-
  soi |>
  ftransform(
    sirs_wbc = case_when(
      wbc > 12             ~ 1L,
      wbc < 04             ~ 1L,
      wbc >= 4 & wbc <= 12 ~ 0L,
      TRUE                 ~ NA_integer_
    ),
    sirs_resp = case_when(
      rr > 20             ~ 1L,
      rr <= 20            ~ 0L,
      TRUE                ~ NA_integer_
    ),
    sirs_pulse = case_when(
      hr > 90             ~ 1L,
      hr <= 90            ~ 0L,
      TRUE                ~ NA_integer_
    ),
    sirs_temp = case_when(
      tmp > 100.4         ~ 1L,
      tmp < 96.8          ~ 1L,
      tmp >= 96.8 & 
        tmp <= 100.4      ~ 0L,
      TRUE                ~ NA_integer_
    )) |>
  group_by(enc) |>
  fill(starts_with("sirs"), .direction = "down") |>
  ungroup() |>
  mutate(
    sirs = 
      rowSums(across(c(sirs_wbc, sirs_resp, sirs_pulse, sirs_temp)), na.rm = T)) |>
  fselect(-sirs_wbc, -sirs_resp, -sirs_pulse, -sirs_temp)

# qSOFA ------------------------------------------------------------------------

soi <-
  soi |>
  ftransform(
    qsofa_bp = case_when(
      sbp <= 100 ~ 1L,
      sbp > 100  ~ 0L,
      TRUE       ~ NA_integer_
    ),
    qsofa_resp = case_when(
      rr >= 22  ~ 1L,
      rr <  22  ~ 0L,
      TRUE      ~ NA_integer_
    ),
    qsofa_gcs = case_when(
      gcs < 15   ~ 1L,
      gcs >= 15  ~ 0L,
      TRUE       ~ NA_integer_
    )) |>
  group_by(enc) |>
  fill(starts_with("qsofa"), .direction = "down") |>
  ungroup() |>
  mutate(
    qsofa = rowSums(across(c(qsofa_bp, qsofa_resp, qsofa_gcs)), na.rm = TRUE)) |>
  fselect(-qsofa_bp, -qsofa_resp, -qsofa_gcs)

# MEWS -------------------------------------------------------------------------

soi <-
  soi |>
  mutate(
    mews_temp = case_when(
      tmp > 101.2    ~ 2L,
      tmp > 100.4    ~ 1L,
      tmp > 96.8     ~ 0L,
      tmp > 95       ~ 1L,
      tmp <= 95      ~ 2L,
      TRUE           ~ NA_integer_
    ),
    mews_resp = case_when(
      rr >= 30       ~ 3L,
      rr >= 21       ~ 2L,
      rr >= 15       ~ 1L,
      rr >= 9        ~ 0L,
      rr <= 8        ~ 2L,
      TRUE           ~ NA_integer_
    ),
    mews_pulse = case_when(
      hr >= 130      ~ 3L,
      hr >= 111      ~ 2L,
      hr >= 101      ~ 1L,
      hr >= 51       ~ 0L,
      hr >= 41       ~ 1L,
      hr <= 40       ~ 2L,
      TRUE           ~ NA_integer_
    ),
    mews_sbp = case_when(
      sbp >= 200    ~ 2L,
      sbp >= 101    ~ 0L,
      sbp >= 81     ~ 1L,
      sbp >= 71     ~ 2L,
      sbp <= 70     ~ 3L,
      TRUE            ~ NA_integer_
    ),
    # AVPU not present in our data - use GCS as per PMID 14687096
    mews_gcs = case_when(
      gcs == 15    ~ 0L,
      gcs >= 13    ~ 1L,
      gcs >= 8     ~ 2L,
      gcs >= 0     ~ 3L,
      TRUE         ~ NA_integer_
    )) |>
  group_by(enc) |>
  fill(starts_with("mews"), .direction = "down") |>
  ungroup() |>
  mutate(
    mews = rowSums(across(c(mews_temp, mews_resp, mews_pulse,
                            mews_sbp, mews_gcs)), na.rm = TRUE)) |>
  fselect(-mews_temp, -mews_resp, -mews_pulse, -mews_sbp, -mews_gcs)

# NEWS -------------------------------------------------------------------------

soi <-
  soi |>
  ftransform(
    news_resp = case_when(
      rr >= 25         ~ 3L,
      rr >= 21         ~ 2L,
      rr >= 12         ~ 0L,
      rr >= 9          ~ 1L,
      rr <= 8          ~ 3L,
      TRUE             ~ NA_integer_
    ),
    news_spo2 = case_when(
      spo2 >= 96       ~ 0L,
      spo2 >= 94       ~ 1L,
      spo2 >= 92       ~ 2L,
      spo2 <= 91       ~ 3L,
      TRUE             ~ NA_integer_
    ),
    news_o2 = case_when(
      o2_01 == 1       ~ 2L,
      o2_01 == 0       ~ 0L,
      TRUE             ~ NA_integer_
    ),
    news_temp = case_when(
      tmp >= 102.3    ~ 2L,
      tmp > 100.4     ~ 1L,
      tmp > 96.8      ~ 0L,
      tmp > 95        ~ 1L,
      tmp <= 95       ~ 3L,
      TRUE            ~ NA_integer_
    ),
    news_sbp = case_when(
      sbp >= 220      ~ 3L,
      sbp >= 111      ~ 0L,
      sbp >= 101      ~ 1L,
      sbp >= 91       ~ 2L,
      sbp <= 90       ~ 3L,
      TRUE            ~ NA_integer_
    ),
    news_pulse = case_when(
      hr >= 131       ~ 3L,
      hr >= 111       ~ 2L,
      hr >= 91        ~ 1L,
      hr >= 51        ~ 0L,
      hr >= 41        ~ 1L,
      hr <= 40        ~ 3L,
      TRUE            ~ NA_integer_
    ),
    news_gcs = case_when(
      gcs < 15        ~ 3L,
      gcs == 15       ~ 0L,
      TRUE            ~ NA_integer_
    )) |>
  group_by(enc) |>
  fill(starts_with("news"), .direction = "down") |>
  ungroup() |>
  mutate(
    news = rowSums(across(c(news_resp, news_spo2, news_temp, news_gcs,
                            news_sbp, news_pulse)), na.rm = TRUE)) |>
  mutate(
    news = case_when(
      news < 7 & news_resp  == 3 ~ 7L,
      news < 7 & news_spo2  == 3 ~ 7L,
      news < 7 & news_temp  == 3 ~ 7L,
      news < 7 & news_sbp   == 3 ~ 7L,
      news < 7 & news_pulse == 3 ~ 7L,
      news < 7 & news_gcs   == 3 ~ 7L,
      TRUE                       ~ news
    )) |>
  fselect(-news_resp, -news_spo2, -news_temp, -news_gcs, -news_sbp, -news_pulse)

# save scores ------------------------------------------------------------------

edi <-
  edi |>
  fselect(-enc) |>
  join(
    df |> fselect(mrn, enc, t_adt, t_dc), 
    how      = "inner",
    multiple = T
  ) |>
  fsubset(time >= t_adt & time <= t_dc) |>
  ftransform(time = lubridate::round_date(time, "20 minutes")) |>
  fsubset(!is.na(edi)) |>
  fgroup_by(enc, time) |>
  fsummarize(edi = fmax(edi))

soi <-
  soi |>
  join(edi, how = "full", multiple = T) |>
  fgroup_by(enc, time) |>
  fmax() |>
  roworder(time) |>
  fill(everything(), .by = enc, .direction = "down") 

df <-
  df |>
  join(soi, how = "left", multiple = T) |>
  fsubset(time >= t_ward & time < t_outcome) |>
  ftransform(
    sirs  = if_else(is.na(sirs),  0L, sirs),
    qsofa = if_else(is.na(qsofa), 0L, qsofa),
    mews  = if_else(is.na(mews),  0L, mews),
    news  = if_else(is.na(news),  0L, news),
    o2_01 = if_else(is.na(o2_01), 0L, o2_01)
  ) |>
  fgroup_by(enc) |>
  fmutate(
    max_hr    = fmax(hr),
    max_rr    = fmax(rr),
    min_sbp   = fmin(sbp),
    min_spo2  = fmin(spo2),
    max_temp  = fmax(tmp),
    min_gcs   = fmin(gcs),
    min_wbc   = fmin(wbc),
    any_o2    = fmax(o2_01),
    max_sirs  = fmax(sirs),
    max_qsofa = fmax(qsofa),
    max_mews  = fmax(mews),
    max_news  = fmax(news),
    max_edi   = fmax(edi)
  ) |>
  fungroup() |> 
  fsubset(!is.na(max_edi)) |>
  fselect(-hr, -rr, -sbp, -spo2, -tmp, -gcs, -wbc, -o2_01)

rm(soi, edi)


# Save df ----------------------------------------------------------------------

df |> write_parquet(here("clean", paste0("analysis_", today, ".parquet")))

# End --------------------------------------------------------------------------
