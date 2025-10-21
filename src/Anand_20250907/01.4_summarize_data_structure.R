data_path <- "data/Anand_20250907"

intermediate_path <- "intermediate/Anand_20250907"
if (!dir.exists(intermediate_path)) {
  cat("Creating intermediate_path '", intermediate_path, "'\n", sep = "")
  dir.create(intermediate_path)
}

output_path <- "product/Anand_20250907"
if (!dir.exists(output_path)) {
  cat("Creating output_path '", output_path, "'\n", sep = "")
  dir.create(output_path)
}
date_code <- "20251021"

load(paste0(intermediate_path, "/MOR_binding.Rdata"))

MOR_treatments <- MOR_binding |>
    dplyr::filter(treatment_type == "treatment")

cat(
    "Data Report\n",
    "\n",
    "N Protocols: ", MOR_treatments |> dplyr::distinct(protocol_id) |> nrow(), "\n",
    "    e.g. '", MOR_treatments$protocol_id[1], "'\n",
    "N Targets: ", MOR_treatments |> dplyr::distinct(target) |> nrow(), "\n",
    "    e.g. '", MOR_treatments$target[1], "'\n",    
    "N Substances: ", MOR_treatments |> dplyr::distinct(substance_id) |> nrow(), "\n",
    "    e.g. '", MOR_treatments$substance_id[1], "'\n",
    "N Researchers: ", MOR_treatments |> dplyr::distinct(researcher) |> nrow(), "\n",
    "    e.g. '", MOR_treatments$researcher[1], "'\n",
    "N Date Codes: ", MOR_treatments |> dplyr::distinct(researcher, date_code) |> nrow(), "\n",
    "   e.g. '", MOR_treatments$date_code[1], "'\n",
    "   per Researcher range: [",
    MOR_treatments |>
        dplyr::distinct(researcher, date_code) |>
        dplyr::count(researcher) |>    
        dplyr::pull("n") |>
        min(), ", ",
    MOR_treatments |>
        dplyr::distinct(researcher, date_code) |>
        dplyr::count(researcher) |>    
        dplyr::pull("n") |>
        max(), "]\n",
    "N Assays: ", MOR_treatments |> dplyr::distinct(researcher, date_code, assay_id) |> nrow(), "\n",
    "   e.g. '", MOR_treatments$assay_id[1], "'\n",    
    "   per Date Code range: [", 
    MOR_treatments |>
        dplyr::distinct(researcher, date_code, assay_id) |>
        dplyr::count(researcher, date_code) |>    
        dplyr::pull("n") |>
        min(), ", ",
    MOR_treatments |>
        dplyr::distinct(researcher, date_code, assay_id) |>
        dplyr::count(researcher, date_code) |>    
        dplyr::pull("n") |>
        max(), "]\n",
    "N Runs: ", MOR_treatments |> dplyr::distinct(researcher, date_code, assay_id, run_id) |> nrow(), "\n",
    "   e.g. '", MOR_treatments$run_id[1], "'\n",    
    "   per Assay range: [", 
    MOR_treatments |>
        dplyr::distinct(researcher, date_code, assay_id, run_id) |>
        dplyr::count(researcher, date_code, assay_id) |>    
        dplyr::pull("n") |>
        min(), ", ",
    MOR_treatments |>
        dplyr::distinct(researcher, date_code, assay_id, run_id) |>
        dplyr::count(researcher, date_code, assay_id) |>    
        dplyr::pull("n") |>
        max(), "]\n",
    "N Replicas: ", MOR_treatments |> dplyr::distinct(replica) |> nrow(), "\n",
    "   e.g. '", MOR_treatments$replica[1], "'\n",    
    "N Log10[Dose]: ", MOR_treatments |> dplyr::distinct(log10_dose) |> nrow(), "\n",
    "   e.g. '", MOR_treatments$log10_dose[1], "'\n",
    sep = "",
    file = paste0(output_path, "/data_structure_summary_", date_code, ".txt"))
