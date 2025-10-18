

mor_fasta <- seqinr::read.fasta(
    file = "data/OPRM_HUMAN.fasta",
    seqtype = "AA")

activity_data <- readxl::read_excel(
  "data/MOR_ligands_20250606.xlsx",
  skip = 2) |>
  dplyr::filter(!(Ligand |> stringr::str_detect("^[*]")))

prediction_path <-
    "intermediate/boltz2_prediction_extra1_20250606/OPRM_HUMAN"

inputs_path <- paste0(prediction_path, "/inputs")
if (!dir.exists(inputs_path)) {
    cat("Creating prediction path: ", inputs_path, "\n", sep = "")
    dir.create(inputs_path, recursive = TRUE)
}

activity_data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand_data <- .
        SMILES <- ligand_data$SMILES
        input_params_path <- paste0(
            inputs_path, "/", ligand_data$Ligand, ".yaml")
        cat("Writing inputs file to ", input_params_path, "\n", sep = "")
        input_params <- list(
            "sequences" = list(
                list("protein" = list(
                    "id" = "A",
                    "sequence" = mor_fasta[[1]] |>
                        paste0(collapse = ""))),
                list("ligand" = list(
                    "id" = "B",
                    "smiles" = SMILES))),
            "properties" = list(
                list("affinity" = list(
                    "binder" = "B"))))
                    
        input_params |>
            yaml::write_yaml(file = input_params_path)


        data.frame()
    })


cmd <- paste0(
    "boltz predict \'", inputs_path, "\' ",
    "--out_dir \'", prediction_path, "\' ",
    "--num_workers 1 ",
    "--use_msa_server ",
    "--recycling_steps 20 ",
    "--diffusion_samples 200 ",
    "--step_scale 1.3 ",
    "--max_parallel_samples 1 ",
#    "--override ",
    "--seed 123 ",
    "--no_trifast", sep = "")
cat(cmd, "\n", sep = "")
system(cmd)


