
library(rcrossref)


intermediate_path <- "intermediate/01.1_gather_KeBOB"
if (!dir.exists(intermediate_path)) {
  cat("Creating ", intermediate_path, "\n", sep = "")
  dir.create(intermediate_path)
}

# References
KeBOB_references <- googlesheets4::read_sheet(
  "https://docs.google.com/spreadsheets/d/18Bptm8MlUOZMUB9GPsEgJzbcGCvtEg9d7PWpomdJTs8/edit?gid=2071146896#gid=2071146896",
  sheet = "References") |>
  dplyr::mutate(
    DOI = Reference |> stringr::str_extract("[(][^,]+, [^,]+, (.+)[)]$", group = 1),
    .before = "Title")


KeBOB_journal_info <- KeBOB_references |>
    dplyr::filter(!is.na(DOI)) |>
    dplyr::rowwise() |>
    dplyr::do({
        DOI <- .$DOI[1]
        cat("Getting Journal info for '", DOI, "' ...\n", sep = "")
        Sys.sleep(0.1) # wait to prevent hammering the crossref server
        tryCatch({
            cr <- rcrossref::cr_works(dois = DOI)
            tibble::tibble(
                DOI = DOI,
                journal_full = cr$data$container.title,
                journal_abbrev = cr$data$short.container.title)
        }, error = function(err) {
            cat("ERROR: Failed to get DOI for '", DOI, "' \n", sep = "")
            data.frame()
        })
    })

KeBOB_references <- KeBOB_references |>
    dplyr::left_join(
        KeBOB_journal_info,
        by = "DOI")

KeBOB_references <- KeBOB_references |>
    dplyr::select(
        index,
        Reference,
        Year,
        DOI,
        journal_full,
        journal_abbrev,
        Title,
        Corresponding,
        `Article Type`,
        `Curation of Compounds`,
        `Curation of Activities`,
        Notes,
        MOR,
        DOR,
        KOR,
        Other,
        Allostery,
        Metal,
        Mutation,
        Priming,
        Simulation,
        `In Vitro`,
        FRET,
        `In Vivo`,
        `Endogenous Ligands`,
        Nitazenes,
        Fentanyls)



KeBOB_references |>
  readr::write_tsv(
      paste0(intermediate_path, "/KeBOB_references.tsv"))
save(KeBOB_references, file = paste0(intermediate_path, "/KeBOB_references.Rdata"))

# Receptors
KeBOB_receptors <- googlesheets4::read_sheet(
  "https://docs.google.com/spreadsheets/d/18Bptm8MlUOZMUB9GPsEgJzbcGCvtEg9d7PWpomdJTs8/edit?gid=2071146896#gid=2071146896",
  sheet = "Receptors")
KeBOB_receptors |>
  readr::write_tsv(
    paste0(intermediate_path, "/KeBOB_receptors.tsv"))
save(KeBOB_receptors, file = paste0(intermediate_path, "/KeBOB_receptors.Rdata"))


# Substances
KeBOB_substances <- googlesheets4::read_sheet(
  "https://docs.google.com/spreadsheets/d/18Bptm8MlUOZMUB9GPsEgJzbcGCvtEg9d7PWpomdJTs8/edit?gid=2071146896#gid=2071146896",
  sheet = "Substances")
KeBOB_substances |>
  readr::write_tsv(
    paste0(intermediate_path, "/KeBOB_substances.tsv"))
save(KeBOB_substances, file = paste0(intermediate_path, "/KeBOB_substances.Rdata"))
