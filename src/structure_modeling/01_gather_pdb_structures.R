

pdb_summary <- readr::read_csv(
  "data/pdb_structures/OPRM1/OPRM1_pdb_summary_20250531.csv",
  show_col_types = FALSE)

pdb_summary |>
  dplyr::rowwise() |>
  dplyr::do({
    pdb_entry <- .$`Entry ID`
    cmd <- paste0("cd data/pdb_structures/OPRM1 && curl -O https://files.rcsb.org/download/", pdb_entry, ".cif")
    cat(cmd, "\n", sep = "")
    system(cmd)
    data.frame()
  })