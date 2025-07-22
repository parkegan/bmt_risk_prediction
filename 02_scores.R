
# This file extracts EDI scores for all patients in the study cohort

# libraries
library(tidytable)
library(collapse)
library(here)


# load  scores -----------------------------------------------------------------

bmt_mrns <- fread(here("clean", "cohort_bmt.csv")) |> fselect(mrn) |> funique()


df <-
  arrow::read_parquet(
    "../000_data/ohsu/intermediate/scores_epic_inpatient.parquet",
    as_data_frame = F
  ) |>
  dplyr::filter(mrn %in% bmt_mrns$mrn & disp_name == "Deterioration Index") |>
  dplyr::collect() |>
  # Clean midnight time stamps
  ftransform(
    recorded_time_date = case_when(
      stringr::str_ends(recorded_time_date, "M") ~ recorded_time_date,
      TRUE                                       ~ paste0(recorded_time_date, " 00:00:01 AM")
    )) |>
  ftransform(
    time       = lubridate::mdy_hms(recorded_time_date),
    meas_value = as.numeric(meas_value)
  ) |>
  rename(edi = meas_value) |>
  fselect(-recorded_time_date)


# save files -------------------------------------------------------------------

fread(here("clean", "cohort_bmt.csv")) |>
  left_join(df) |>
  fsubset(time >= admit_dttm & time <= disch_dttm) |>
  fselect(mrn, enc, time, edi) |>
  arrow::write_parquet(here("clean", "epic_scores.parquet"))
