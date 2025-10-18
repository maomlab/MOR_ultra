
proteins <- data.frame(
    uniprot_entry = c(
        "OPRM_HUMAN",
        "OPRD1_HUMAN",
        "OPRK1_HUMAN",
        "OPRL1_HUMAN",
        "AA2AR_HUMAN")) |>
    dplyr::rowwise() |>
    dplyr::do({
        uniprot_entry <- .$uniprot_entry
        path <- paste0("data/uniprot/", uniprot_entry, ".fasta")
        fasta <- path |> seqinr::read.fasta(seqtype = "AA")
        fasta <- fasta[[1]] |>  paste0(collapse = "")
        data.frame(
            uniprot_entry = uniprot_entry,
            fasta = fasta)
    })

activity_data <- readxl::read_excel(
  "data/MOR_ligands_20250606.xlsx",
  skip = 2) |>
  dplyr::filter(!(Ligand |> stringr::str_detect("^[*]")))

prediction_path <-
    "intermediate/boltz2_prediction_dimer_20250712"

inputs_path <- paste0(prediction_path, "/inputs")
if (!dir.exists(inputs_path)) {
    cat("Creating prediction path: ", inputs_path, "\n", sep = "")
    dir.create(inputs_path, recursive = TRUE)
}

runs <- dplyr::cross_join(
    proteins |>
    dplyr::rename(
        uniprot_entry_1 = uniprot_entry,
        fasta_1 = fasta),
    proteins |>
    dplyr::rename(
        uniprot_entry_2 = uniprot_entry,
        fasta_2 = fasta)) |>
    dplyr::cross_join(activity_data)

cat("Defined ", nrow(runs), " runs\n", sep = "")


runs |>
    dplyr::rowwise() |>
    dplyr::do({
        run_data <- .
        input_params_path <- paste0(
            inputs_path, "/",
            run_data$uniprot_entry_1, ";",
            run_data$uniprot_entry_2, ";",            
            run_data$Ligand, ".yaml")
        cat("Writing inputs file to ", input_params_path, "\n", sep = "")
        input_params <- list(
            "sequences" = list(
                list("protein" = list(
                    "id" = "A",
                    "sequence" = run_data$fasta_1 |>
                        paste0(collapse = ""))),
                list("protein" = list(
                    "id" = "B",
                    "sequence" = run_data$fasta_2 |>
                        paste0(collapse = ""))),
                list("ligand" = list(
                    "id" = "C",
                    "smiles" = run_data$SMILES))),
            "properties" = list(
                list("affinity" = list(
                    "binder" = "C"))))
                    
        input_params |>
            yaml::write_yaml(file = input_params_path)


        data.frame()
    })


cmd <- paste0(
    "boltz predict \'", inputs_path, "\' ",
    "--out_dir \'", prediction_path, "\' ",
    "--num_workers 1 ",
    "--use_msa_server ",
#    "--recycling_steps 20 ",
#    "--diffusion_samples 200 ",
#    "--step_scale 1.3 ",
    "--max_parallel_samples 1 ",
#    "--override ",
    "--seed 123 ",
    "--no_trifast", sep = "")
cat(cmd, "\n", sep = "")
system(cmd)


