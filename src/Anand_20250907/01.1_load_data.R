
# this script loads the binding and GTPyS MOR data produced by
# trainees working with Jess Anand

# The data are organized by as follows
# data/Anand_20250907/
#  <protocol_id>/
#    <researcher>/
#      <date_code>/
#        <assay_id>.pzfx
#          <run_id>
#            <substance_id>
#              <replica>
#
# Where
#   <protocol_id>: "MOR Binding Ki Data", "MOR GTPyS Data"
#   <researcher>:  "Allyson", "Arturo", etc.
#   <date_code>: "3H 2019/MOR/4.1.19", "Competition Binding 2023/2.17.23", etc.
#      Note the <date_code> can be multiple nested folders
#      all assays in in the same date_code folder are assumed to be done with the same batch?
#   <assay_id>: "DRC25.21.20.KOR. 3H. 2.7.23" "DRC9-DRC10-DRC3-DRC8 hDOR competition binding 02.22.22", etc.
#      Note the <assay_id> is the filename before the .pzfx
#   <run_id>: "01_Competition Binding 17421", "03_competition bindin—DRC4"
#      Note the run_id is <table_index>_<table_name> in the pzfx, where table_index is 2-digit zero padded
#      The table_index values are not necessarily sequential as some tables may be excluded
#   <substance_id>: "17421", "phenylephrine", etc.
#.     Note: substance_ids are normalized by the following rules
#        * no spaces (use dashes if needed)
#        * no nicknames (e.g. U69593 instead of just U69)
#        * Common acronyms are used (e.g. DAMGO for [D-Ala2, N-MePhe4, Gly-ol]-enkephalin))
#        * generic drug names are lower case (e.g. carfentanil)
#        * organization drug ids use the organization convention e.g. (U69593)
#   <replica>: "1", "2"
#      Note: the ycolumn header in the pzfx data table is parsed as <substance_id>_<replica>

# for assay_id: "MOR Binding Ki Data", the measured data are
#   <log10_dose>: -15, -14, etc.
#   <value>:
#      readout, maybe fluorescence intensity? (e.g. 154)
#      data points[1,2] are assumed to be negative controls (i.e. bind at 0%)
#      data points[3,4] are assumed to be positive controls (i.e. bind at 100%)
#
#   computed values:
#      <value_normalized>: values normalized to the positive and negative controls
#          bottom_mean = (value[1] + value[2]) / 2
#          top_mean = (value[3] + value[4]) / 2
#          value_normalized = (value - bottom_mean) / (top_mean - bottom_mean)

library(pzfx)

# in the pzfx files there are some empty columns that mess up parsing
# this script strips them from the xml before parsing with the pzfx package
source(here::here("src/Anand_20250907/pzfx_remove_empty_ycolumns.R"))



data_path <- "data/Anand_20250907"

intermediate_path <- "intermediate/Anand_20250907"
if (!dir.exists(intermediate_path)) {
    cat("Creating output_path '", intermediate_path, "'\n", sep = "")
    dir.create(intermediate_path)
}

output_path <- "product/Anand_20250907"
if (!dir.exists(output_path)) {
    cat("Creating output_path '", output_path, "'\n", sep = "")
    dir.create(output_path)
}
dataset_name <- "Anand Data 20250907"
date_code <- "20250907"




pzf_files <- data.frame(
  pzf_fname = data_path |>
    list.files(pattern = "*pzf$", full.names = TRUE, recursive = TRUE))

# convert pzf to pzfx
# this uses an applescript to open the pzf file in Prism 10 and save it again as a .pzfx file
# Only tested on MacOS with Prism10
# It takes a few seconds per file so ~10-20 minutes overall
pzf_files |>
  dplyr::rowwise() |>
  dplyr::do({
    data <- .
    pzf_fname <- data$pzf_fname[1]
    pzfx_fname <- pzf_fname |> paste0("x")

    cat("Converting '", pzf_fname, "' => '", pzfx_fname, "' ...\n", sep = "")
    if(file.exists(pzfx_fname)) {
      cat(pzfx_fname, " already exists, skipping!\n")
    } else {
      cmd <- paste0(
      "osascript ",
        here::here("src/Anand_20250907/prism_pzf_to_pzfx.scpt"), " ",
        shQuote(here::here(pzf_fname)), " ",
        shQuote(here::here(pzfx_fname)), " ")
      cat(cmd, "\n")
      system(cmd)
      cat("\n\n")
    }
    data.frame()
  })


# check if all pzf files now have matching pzfx files
pzf_files |>
  dplyr::mutate(
    pzfx_fname = pzf_fname |> paste0("x")) |>
  dplyr::do({
    pzfx_fname <- .$pzfx_fname[1]
    if(!file.exists(pzfx_fname)) {
      cat("WARNING: ", pzfx_fname, " does not exist.\n", sep = "")
    }
    data.frame()
  })

# PARSE MOR binding Data
# Read each .pzfx file and concatenate into a single data table

protocol_id_label = "MOR Binding Ki Data"
MOR_binding <- data.frame(
  pzfx_fname = paste0(data_path, "/", protocol_id_label) |>
      list.files(pattern = "*pzfx$", full.names = TRUE, recursive = TRUE)) |>
  dplyr::rowwise() |>
  dplyr::do({
    pzfx_fname <- .$pzfx_fname[1]
      cat("Reading pzfx file: '", pzfx_fname, "'\n", sep = "")

      # remove empty y-columns and save to a temporary .pzfx file
      pzfx_fname_clean <- tempfile("prism_temp_", fileext = ".pzfx")
      cat("Cleaning empty columns and saving to '", pzfx_fname_clean, "'\n", sep = "")
      pzfx_fname |> pzfx_remove_empty_ycolumns(pzfx_fname_clean)

      data <- pzfx::pzfx_tables(pzfx_fname_clean) |>
      purrr::imap_dfr(function(table_name, table_index) {
          # skip reading tables where the table name is "std "
          run_id <- paste0(
            table_index |> stringr::str_pad(width = 2, side = "left", pad = "0"),
            "_",
            table_name)

          if (table_name %in% c("std ")) {
              cat("  Skipping table '", run_id, "'\n", sep = "")

          } else {
              cat("  Reading table: '", run_id, "'\n", sep = "")
              tryCatch({

                # there seem to be a few columns with no data with columns headers that
                # start with "α" or "_" => remove them
                data <- pzfx::read_pzfx(
                    path = pzfx_fname_clean,
                    table = table_name) |>
                  dplyr::select(
                    -tidyselect::starts_with("α"),
                    -tidyselect::starts_with("_"))

                # There appear to be two different formats
                # parse tables with 3 columns: [log10_dose, substance_id, value]
                if (ncol(data) == 3) {
                  data <- data |>
                    dplyr::rename(log10_dose = 1) |>
                    dplyr::mutate(
                      ROWTITLE = NA,
                      .before = 1)

                  data <- data |>
                    dplyr::mutate(
                      pzfx_fname = pzfx_fname,
                      table_name = table_name,
                      .before = 1)

                  data <-  data |>
                    tidyr::pivot_longer(
                      cols = -c("pzfx_fname", "table_name", "ROWTITLE", "log10_dose"),
                      names_to = "substance_id",
                      values_to = "value") |>
                    dplyr::mutate(run_id = run_id)

                # parse tables with 4 columns: [log10_dose, substance_id, value]
                } else if (ncol(data) == 4) {
                  data <- data |>
                    dplyr::rename(log10_dose = 2)

                  data <- data |>
                    dplyr::mutate(
                      pzfx_fname = pzfx_fname,
                      table_name = table_name,
                      .before = 1)

                  data <-  data |>
                    tidyr::pivot_longer(
                      cols = -c("pzfx_fname", "table_name", "ROWTITLE", "log10_dose"),
                      names_to = "substance_id",
                      values_to = "value") |>
                    dplyr::mutate(
                      run_id = run_id)

                # raise error when there are fewer than 3 or more than 4 columns
                } else {
                  # break if we hit this edge case to inspect
                  cat("error reading table '", table_name, "' in '", pzfx_fname, "' ...\n", sep = "")
                  data <- data.frame(
                    pzfx_fname = pzfx_fname,
                    table_name = table_name,
                    run_id = run_id,
                    ROWTITLE = NA,
                    log10_dose = NA,
                    substance_id = NA,
                    value = NA)
                }
              }, error = function(msg) {
                cat("error reading table '", table_name, "' in '", pzfx_fname, "' ...\n", sep = "")
                data <- data.frame(
                  pzfx_fname = pzfx_fname,
                  table_name = table_name,
                  run_id = run_id,
                  ROWTITLE = NA,
                  log10_dose = NA,
                  substance_id = NA,
                  value = NA)
              })

          }
      })

      # some data files don't actually have any data
      if( data |> dplyr::filter(!is.na(value)) |> nrow() == 0) {
        data <- data.frame()
      }
      data
    })

MOR_binding <- MOR_binding |>
  dplyr::mutate(
    protocol_id = pzfx_fname |> stringr::str_extract(paste0(data_path, "/(MOR Binding Ki Data)/([^/]+)/(.+)/([^/]+pzfx)$"), group = 1),
    researcher = pzfx_fname |> stringr::str_extract(paste0(data_path, "/(MOR Binding Ki Data)/([^/]+)/(.+)/([^/]+pzfx)$"), group = 2),
    date_code = pzfx_fname |> stringr::str_extract(paste0(data_path, "/(MOR Binding Ki Data)/([^/]+)/(.+)/([^/]+pzfx)$"), group = 3),
    assay_id = pzfx_fname |> stringr::str_extract(paste0(data_path, "/(MOR Binding Ki Data)/([^/]+)/(.+)/([^/]+)[.]pzfx$"), group = 4))


MOR_binding <- MOR_binding |>
  dplyr::mutate(
    target = dplyr::case_when(
      assay_id == ("lobelineMOR-lobelineDOR-lobelineKOR 3H competition binding 07.20.22") &
        run_id == "02_competition binding lobeline" ~ "OPRM_HUMAN",
      assay_id == ("lobelineMOR-lobelineDOR-lobelineKOR 3H competition binding 07.20.22") &
        run_id == "03_competition binding lobeline" ~ "OPRD_HUMAN",
      assay_id == ("lobelineMOR-lobelineDOR-lobelineKOR 3H competition binding 07.20.22") &
        run_id == "04_competition binding lobeline" ~ "OPRK_HUMAN",

      assay_id == ("lobelinemor-lobelinekor-phenylephrinemor-phenylephrinekor 3H competition binding 7.22.22") &
        run_id == "02_competition binding lobeline" ~ "OPRM_HUMAN",
      assay_id == ("lobelinemor-lobelinekor-phenylephrinemor-phenylephrinekor 3H competition binding 7.22.22") &
        run_id == "03_competition binding lobeline" ~ "OPRK_HUMAN",
      assay_id == ("lobelinemor-lobelinekor-phenylephrinemor-phenylephrinekor 3H competition binding 7.22.22") &
        run_id == "04_competition binding phenylephrine" ~ "OPRM_HUMAN",
      assay_id == ("lobelinemor-lobelinekor-phenylephrinemor-phenylephrinekor 3H competition binding 7.22.22") &
        run_id == "05_competition binding phenylephrine" ~ "OPRK_HUMAN",

      assay_id == ("phenylephrine MOR.DOR.KOR competition binding 01.31.23") &
        run_id == "02_competition binding lobeline" ~ "OPRM_HUMAN",
      assay_id == ("phenylephrine MOR.DOR.KOR competition binding 01.31.23") &
        run_id == "03_competition binding lobeline" ~ "OPRD_HUMAN",
      assay_id == ("phenylephrine MOR.DOR.KOR competition binding 01.31.23") &
        run_id == "04_competition binding lobeline" ~ "OPRK_HUMAN",

      assay_id == ("20859-22807-14719-BINDING CHO HUMAN DOR(NLX at MOR) 01.16.20") &
        run_id == "02_competition binding template 20859" ~ "OPRD_HUMAN",
      assay_id == ("20859-22807-14719-BINDING CHO HUMAN DOR(NLX at MOR) 01.16.20") &
        run_id == "03_competition binding template 22807" ~ "OPRD_HUMAN",
      assay_id == ("20859-22807-14719-BINDING CHO HUMAN DOR(NLX at MOR) 01.16.20") &
        run_id == "04_competition binding template 14719" ~ "OPRD_HUMAN",
      assay_id == ("20859-22807-14719-BINDING CHO HUMAN DOR(NLX at MOR) 01.16.20") &
        run_id == "05_competition binding template NLX" ~ "OPRM_HUMAN",

      assay_id |> stringr::str_detect(
        "(hMOR| MOR | mor |[(]MOR[)]|MOR[.]|[.]MOR| mU |Human Mor|HUMAN MOR|Human MOR|human-mor)") ~ "OPRM_HUMAN",
      assay_id |> stringr::str_detect(
        "(hKOR|hkOR| KOR |[.]KOR|KOR[.]| KOR$| kor |[(]KOR[)]|Human Kor|Human KOR)") ~ "OPRK_HUMAN",
      assay_id |> stringr::str_detect(
        "(hDOR| DOR | dor |[.]DOR|DOR[.]|[(]DOR[)]|Human Dor|Human DOR|Human dor|Human DELTA)") ~ "OPRD_HUMAN",
      assay_id |> stringr::str_detect("rMOR") ~ "OPRM_RAT",
      assay_id |> stringr::str_detect("rKOR") ~ "OPRK_RAT",
      assay_id |> stringr::str_detect("rDOR") ~ "OPRD_RAT",


      run_id |> stringr::str_detect("MOR$") ~ "OPRM_HUMAN",
      run_id |> stringr::str_detect("DOR$") ~ "OPRD_HUMAN",
      run_id |> stringr::str_detect("KOR$") ~ "OPRK_HUMAN",

      TRUE ~ "OPRM_HUMAN"),
    .before = 8)


MOR_binding <- MOR_binding |>
  dplyr::mutate(ROWTITLE = ROWTITLE |> stringr::str_trim()) |>
  dplyr::mutate(
    ROWTITLE = dplyr::case_when(
      ROWTITLE == "nlx" ~ "NLX",
      ROWTITLE == "TRIS" ~ "Tris",
      ROWTITLE == "tris" ~ "Tris",
      TRUE ~ ROWTITLE))

MOR_binding <- MOR_binding |>
    tidyr::separate_wider_delim(
        cols = substance_id,
        names = c("substance_id", "replica"), delim = "_")

MOR_binding <- MOR_binding |>
  dplyr::mutate(
    substance_id_orig = substance_id)

# For many datasets, Rasha only put the compound ID in the table name
MOR_binding <- MOR_binding |>
  dplyr::mutate(
    substance_id = ifelse(
      (researcher == "Rasha") & is.na(substance_id),
      run_id |> stringr::str_extract("(U69|DRC[0-9]+)"),
      substance_id))

MOR_binding <- MOR_binding |>
  dplyr::mutate(
    substance_id = substance_id_orig |> stringr::str_replace_all(" ", "")) |>
  dplyr::mutate(
    substance_id = dplyr::case_when(
      substance_id == "Alfentanil" ~ "alfentanil",
      substance_id == "Bremazocine" ~ "bremazocine",
      substance_id == "Buprenorphine" ~ "buprenorphine",
      substance_id == "CARFENTANIL" ~ "carfentanil",
      substance_id == "Carfentanil" ~ "carfentanil",
      substance_id == "Cirazoline" ~ "cirazoline",
      substance_id == "Clonidine" ~ "clonidine",
      substance_id == "DRC-37" ~ "DRC37",
      substance_id == "DRC-38" ~ "DRC38",
      substance_id == "DRC-39" ~ "DRC39",
      substance_id == "DRC-40" ~ "DRC40",
      substance_id == "DRC2-EA" ~ "DRC2EA",
      substance_id == "DRC2-LA" ~ "DRC2LA",
      substance_id == "Damgo" ~ "DAMGO",
      substance_id == "Fentalog1" ~ "fentalog1",
      substance_id == "Fentalog2" ~ "fentalog2",
      substance_id == "Fentanyl" ~ "fentanyl",
      substance_id == "Morhpine" ~ "morphine",
      substance_id == "Morphine" ~ "morphine",
      substance_id == "NLXM" ~ "NLX-M",
      substance_id == "NLXMEI" ~ "NLX-MeI",
      substance_id == "NLXMEL" ~ "NLX-MEL",
      substance_id == "NLXMeI" ~ "NLX-MeI",
      substance_id == "Phenylephrine" ~ "phenylephrine",
      substance_id == "Prazosin" ~ "prazosin",
      substance_id == "RMethadone" ~ "R-methadone",
      substance_id == "R-Methadone" ~ "R-methadone",
      substance_id == "Sufentail" ~ "sufentanil",
      substance_id == "Sufentanil" ~ "sufentanil",
      substance_id == "U69" ~ "U69593",
      substance_id == "damgo" ~ "DAMGO",
      substance_id == "damgo1" ~ "DAMGO1",
      substance_id == "damgo2" ~ "DAMGO2",
      substance_id == "drc2" ~ "DRC2",
      substance_id == "drc4" ~ "DRC4",
      substance_id == "mmp220" ~ "MMP220",
      substance_id == "mmp2200" ~ "MMP2200",
      substance_id == "nlxmel" ~ "NXL-MEL",
      TRUE ~ substance_id))

MOR_binding |>
  dplyr::group_by(pzfx_fname, table_name, date_code, assay_id, run_id) |>
  dplyr::mutate(index = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::filter(index <= 8) |>
  dplyr::count(index, ROWTITLE)

MOR_binding <- MOR_binding |>
  dplyr::group_by(protocol_id, researcher, date_code, assay_id, run_id) |>
  dplyr::mutate(
    treatment_type = dplyr::case_when(
      dplyr::row_number() %in% 1:4 ~ "negative_control",
      dplyr::row_number() %in% 5:8 ~ "positive_control",
      TRUE ~ "treatment")) |>
    dplyr::mutate(
        bottom_mean = value[1:4] |> mean(na.rm = TRUE),
        top_mean = value[5:8] |> mean(na.rm = TRUE),
        value_normalized = (value - bottom_mean) / (top_mean - bottom_mean)) |>
    dplyr::ungroup()

MOR_binding <- MOR_binding |>
  dplyr::select(
    pzfx_fname,
    table_name,
    protocol_id,
    researcher,
    date_code,
    assay_id,
    run_id,
    target,
    treatment_type,
    bottom_mean,
    top_mean,
    substance_id_orig,
    substance_id,
    replica,
    log10_dose,
    value,
    value_normalized)

save(MOR_binding, file = paste0(intermediate_path, "/MOR_binding.Rdata"))
