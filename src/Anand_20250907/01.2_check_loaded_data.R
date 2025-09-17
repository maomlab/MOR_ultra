

data_path <- "data/Anand_20250907"

intermediate_path <- "intermediate/Anand_20250907"
if (!dir.exists(intermediate_path)) {
  cat("Creating output_path '", intermediate_path, "'\n", sep = "")
  dir.create(intermediate_path)
}

date_code <- "20250915"

load(paste0(intermediate_path, "/MOR_binding.Rdata"))

############################################
# check file names are parsed consistently #
############################################
check_data <- MOR_binding |>
  dplyr::filter(
    pzfx_fname == paste0(
      data_path, "/", researcher, "/", date_code, "/", assay_id, ".pzfx")) |>
  dplyr::distinct(pzfx_fname, .keep_all = TRUE)
assertthat::assert_that(
  nrow(check_data) == 0,
  msg = paste0(nrow(check_data), " pzfx_fname values are not parsed correctly."))

###########################
# Manually Check Targets: #
###########################
MOR_binding |>
  dplyr::distinct(
    pzfx_fname,
    table_name,
    protocol_id,
    researcher,
    date_code,
    assay_id,
    run_id,
    target) |>
  readr::write_tsv(
    file = paste0(intermediate_path, "/MOR_binding_check_parse_targets_", date_code, ".tsv"))

#######################################
# check that substance_ids are not NA #
#######################################
check_data <- MOR_binding |>
  dplyr::distinct(pzfx_fname, table_name, .keep_all = TRUE) |>
  dplyr::filter(is.na(substance_id) | is.na(replica))
assertthat::assert_that(
  nrow(check_data) == 0,
  msg = paste0(nrow(check_data), " substance_id or replica values are NA"))

# check that log10_dose values are reasonable
check_data <- MOR_binding |>
  dplyr::distinct(pzfx_fname, table_name, log10_dose, .keep_all = TRUE) |>
  dplyr::filter((log10_dose < -20) | (log10_dose > 0))
assertthat::assert_that(
  nrow(check_data) == 0,
  msg = paste0(nrow(check_data), " log10_dose < 20 or log10_dose > 0"))

###############################################
# check that log10_dose values are reasonable #
###############################################
check_data <- MOR_binding |>
  dplyr::distinct(pzfx_fname, table_name, log10_dose, .keep_all = TRUE) |>
  dplyr::filter((log10_dose < -12)) |>
  dplyr::filter(treatment_type == "treatment") |>
  dplyr::count(protocol_id, researcher, date_code, assay_id, run_id)
assertthat::assert_that(
  nrow(check_data) == 0,
  msg = paste0(nrow(check_data), " log10_dose < 20 or log10_dose > 0"))

##############################################
# Check values are within a reasonable range #
##############################################
# out of reasonable range can happen e.g. if the difference in activity between
# the positive and negative controls are small
check_data <- MOR_binding |>
  dplyr::filter(
    treatment_type == "treatment",
    value_normalized > 3)
assertthat::assert_that(
  nrow(check_data) == 0,
  msg = paste0(nrow(check_data), " value_normalized < 300%"))

check_data <- MOR_binding |>
  dplyr::filter(
    treatment_type == "treatment",
    value_normalized < -1)
assertthat::assert_that(
  nrow(check_data) == 0,
  msg = paste0(nrow(check_data), " value_normalized > -100%"))


