
# This file identifies all patients who underwent BMT on or after 1/1/2019
# and their subsequent hospitalizations

# libraries
library(tidytable)
library(collapse)
library(here)

today = format(Sys.Date(), "%y%m%d")

# start with all OHSU patients -------------------------------------------------

# load all ADT data
adt <-
  arrow::read_parquet(here("../000_data/ohsu/clean", "adt_ohsu.parquet")) |>
  ## keep only the hospitalizations (NOT ED-only visits)
  ftransform(
    admitted_01 = case_when(
      loc_cat == "icu"   ~ 1,
      loc_cat == "wards" ~ 1,
      TRUE               ~ 0
    )) |>
  fgroup_by(enc) |>
  fmutate(admitted_01 = fmax(admitted_01)) |>
  fungroup() |>
  fsubset(admitted_01 == 1) |>
  select(enc, admit_dttm_adt, disch_dttm_adt, time = time_start, loc_cat) |>
  funique()

# flag important times
adt <-
  adt |>
  roworder(time) |>
  fgroup_by(enc, loc_cat) |>
  fmutate(mintime = fmin(time)) |>
  fungroup() |>
  ftransform(
    mintime_wards = if_else(loc_cat == "wards", mintime, NA),
    mintime_icu   = if_else(loc_cat == "icu",   mintime, NA)
  ) |>
  fgroup_by(enc) |>
  fmutate(
    loc_first     = ffirst(loc_cat),
    mintime_wards = fmin(mintime_wards),
    mintime_icu   = fmin(mintime_icu)
  ) |>
  fungroup() |>
  fselect(-time, -loc_cat, -mintime) |>
  funique()

# check for problems with these times
# (i.e., no patient should be missing both mintime_wards and mintime_icu)
adt |> fsubset(is.na(mintime_wards) & is.na(mintime_icu)) # no problems!

# subset the BMT patients ------------------------------------------------------
df <-
  fread(here("../009_bmt_riskprediction/pre/bmt_chart_review.csv")) |>
  fselect(mrn, date_of_hct, donor_type) |>
  inner_join(adt) |>
  rename(disch_dttm = disch_dttm_adt) |>
  ftransform(admit_dttm = pmin(admit_dttm_enc, admit_dttm_adt)) |>
  fselect(-last_name, -first_name, -admit_dttm_enc, -admit_dttm_adt, -pat_id)

# exclude admission from prior to transplant
df |>
  fsubset(date_of_hct > disch_dttm) |>
  fselect(enc, admit_dttm, disch_dttm, date_of_hct, dischg_dispo) |>
  ftransform(
    dc_to_hct    = as.numeric(difftime(date_of_hct, disch_dttm), "days"),
    dischg_dispo = qF(dischg_dispo)
  ) |>
  summary()

# save df ----------------------------------------------------------------------

df |>
  fsubset(date_of_hct <= disch_dttm) |>
  relocate(admit_dttm, disch_dttm, loc_first, .before = mintime_wards) |>
  fwrite(here("clean", paste0("cohort_bmt_", today, ".csv")))

rm(adt, df)

# End --------------------------------------------------------------------------