


intermediate_path <- "intermediate/01.3_gather_chembl"
if (!dir.exists(intermediate_path)) {
  cat("Creating ", intermediate_path, "\n", sep = "")
  dir.create(intermediate_path)
}

chembl_con <- DBI::dbConnect(
  RSQLite::SQLite(),
  dbname = "data/ChEMBL35/chembl_35/chembl_35_sqlite/chembl_35.db")


chembl_KeBOB_receptors <- readr::read_tsv(
  "intermediate/01.1_gather_KeBOB/KeBOB_receptors.tsv",
  show_col_types = FALSE) |>
  dplyr::rename(chembl_id = "ChEMBL ID") |>  
  dplyr::distinct(chembl_id, .keep_all = TRUE)
          
dplyr::copy_to(
  chembl_con,
  chembl_KeBOB_receptors,
  overwrite = TRUE)

chembl_con |> DBI::dbListTables()

target_dictionary_tbl <- chembl_con |> dplyr::tbl("target_dictionary") |>
    dplyr::rename(target_chembl_id = chembl_id)    

KeBOB_receptors_tbl <- chembl_con |> dplyr::tbl("chembl_KeBOB_receptors") |>
    dplyr::rename(target_chembl_id = chembl_id)

assays_tbl <- chembl_con |> dplyr::tbl("assays") |>
    dplyr::rename(assay_chembl_id = chembl_id)    

docs_tbl <- chembl_con |> dplyr::tbl("docs") |>
    dplyr::rename(doc_chembl_id = chembl_id)

activities_tbl <- chembl_con |> dplyr::tbl("activities")

compound_id_tbl <- chembl_con |>
    dplyr::tbl("chembl_id_lookup") |>
    dplyr::filter(entity_type == "COMPOUND") |>
    dplyr::select(
        molregno = entity_id,
        chemical_chembl_id = chembl_id)

compound_structures_tbl <- chembl_con |> dplyr::tbl("compound_structures") |>
    dplyr::select(
        molregno,
        standard_inchi,
        standard_inchi_key,
        canonical_smiles)

molecule_synonyms <- chembl_con |>
    dplyr::tbl("molecule_synonyms") |>
    dplyr::collect(n = Inf)

chembl_KeBOB_assays <- target_dictionary_tbl |>
  dplyr::semi_join(KeBOB_receptors_tbl, by = "target_chembl_id") |>
  dplyr::rename(target_chembl_id = chembl_id) |>
  dplyr::left_join(assays_tbl, by = "tid") |>
  dplyr::collect(n = Inf)

chembl_KeBOB_documents <- target_dictionary_tbl |>
  dplyr::semi_join(KeBOB_receptors_tbl, by = "target_chembl_id") |>
  dplyr::left_join(assays_tbl, by = "tid") |>
  dplyr::distinct(doc_id) |>
  dplyr::select(doc_id) |>
  dplyr::left_join(docs_tbl, by = "doc_id") |>
  dplyr::collect(n = Inf)

chembl_KeBOB_documents |>
  readr::write_tsv(
    paste0(intermediate_path, "/chembl_KeBOB_documents.tsv"))
save(chembl_KeBOB_documents, file = paste0(intermediate_path, "/chembl_KeBOB_documents.Rdata"))



chembl_KeBOB_activities <- target_dictionary_tbl |>
  dplyr::semi_join(KeBOB_receptors_tbl, by = "target_chembl_id") |>
  dplyr::left_join(assays_tbl, by = "tid") |>
  dplyr::left_join(activities_tbl, by = c("assay_id", "doc_id", "src_id")) |>
  dplyr::left_join(docs_tbl, by = c("doc_id", "src_id")) |>
  dplyr::left_join(compound_id_tbl, by = "molregno") |>
  dplyr::left_join(compound_structures_tbl, by = "molregno") |>
  dplyr::collect(n = Inf)


KeBOB_molecule_synonyms <- molecule_synonyms |>
    dplyr::filter(syn_type %in% c("INN", "FDA",  "RESEARCH_CODE", "OTHER")) |>
    dplyr::semi_join(chembl_KeBOB_activities, by = "molregno") |>
    dplyr::group_by(molregno, syn_type) |>
    dplyr::summarize(synonyms = paste(synonyms, collapse = "|"), .groups = "drop") |>
    tidyr::pivot_wider(id_cols = "molregno", names_from = "syn_type", values_from = "synonyms") |>
    dplyr::mutate(
        compound_name = dplyr::case_when(
            !is.na(INN) ~ INN,
            !is.na(FDA) ~ FDA,
            !is.na(OTHER) ~ RESEARCH_CODE,
            TRUE ~ OTHER)) |>
    dplyr::select(
        molregno, compound_name, INN, FDA, RESEARCH_CODE, OTHER)

chembl_KeBOB_activities <- chembl_KeBOB_activities |>
    dplyr::left_join(
        KeBOB_molecule_synonyms,
        by = "molregno")

chembl_KeBOB_activities |>
  readr::write_tsv(
    paste0(intermediate_path, "/chembl_KeBOB_activities.tsv"))


chembl_KeBOB_substances <- chembl_KeBOB_activities |>
    dplyr::distinct(molregno, .keep_all = TRUE) |>
    dplyr::transmute(
        molregno,
        substance_id =  chemical_chembl_id, 
        chemical_chembl_id,
        standard_inchi,
        standard_inchi_key,
        canonical_smiles,
        compound_name,
        INN,
        FDA,
        RESEARCH_CODE,
        OTHER)

chembl_KeBOB_substances |>
  readr::write_tsv(
    paste0(intermediate_path, "/chembl_KeBOB_substances.tsv"))
