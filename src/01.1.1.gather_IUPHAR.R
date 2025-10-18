

# Version 2025.2 of the IUPHAR/BPS
# Guide to Pharmacology database was released on 18th June 2025
# 2025.2

staging_directory <- "data/IUPHAR/dump"

if (!dir.exists(staging_directory)) {
  cat("Creating staging directory: ", staging_directory, "\n", sep = "")
  dir.create(staging_directory)
}

#########################
HGNC_mapping_fname <- paste0(staging_directory, "/GtP_to_HGNC_mapping.csv")
httr::GET(
  url = "http://www.guidetopharmacology.org/DATA/GtP_to_HGNC_mapping.csv",
  httr::write_disk(HGNC_mapping_fname, overwrite=TRUE))

uniprot_mapping_fname <- paste0(staging_directory, "/GtP_to_UniProt_mapping.csv")
httr::GET(
  url = "http://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv",
  httr::write_disk(uniprot_mapping_fname, overwrite=TRUE))


targets_fname <- paste0(staging_directory, "/targets_and_families.csv")
httr::GET(
  url = "http://www.guidetopharmacology.org/DATA/targets_and_families.csv",
  httr::write_disk(targets_fname, overwrite=TRUE))

ligands_fname <- paste0(staging_directory, "/ligands.csv")
httr::GET(
  url = "http://www.guidetopharmacology.org/DATA/ligands.csv",
  httr::write_disk(ligands_fname, overwrite=TRUE))

ligand_id_mapping_fname <- paste0(staging_directory, "/ligand_id_mapping.csv")
httr::GET(
  url = "https://www.guidetopharmacology.org/DATA/ligand_id_mapping.csv",
  httr::write_disk(ligand_id_mapping_fname, overwrite=TRUE))


peptides_fname <- paste0(staging_directory, "/peptides.csv")
httr::GET(
  url = "http://www.guidetopharmacology.org/DATA/peptides.csv",
  httr::write_disk(peptides_fname, overwrite=TRUE))

interactions_fname <- paste0(staging_directory, "/interactions.csv")
httr::GET(
  url = "http://www.guidetopharmacology.org/DATA/interactions.csv",
  httr::write_disk(interactions_fname, overwrite=TRUE))



gtp_rdf_fname <- paste0(staging_directory, "/gtp-rdf.n3")
httr::GET(
  url = "https://www.guidetopharmacology.org/DATA/rdf/2025.2/gtp-rdf.n3",
  httr::write_disk(gtp_rdf_fname, overwrite = TRUE))


gtpo <- rdflib::rdf_parse("data/IUPHAR/dump/gtp-rdf.n3", format = "turtle")

opioid_targets <- gtpo |> rdflib::rdf_query("
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX gtpo: <https://rdf.guidetopharmacology.org/ns/gtpo#>

  SELECT ?targetLabel
  WHERE {
  	?target gtpo:hasTargetFamily ?targetFamily .
  	?targetFamily gtpo:targetFamilyName \"Opioid receptors\" .

  	?target rdfs:label ?targetLabel
  }")

interactions <- gtpo |> rdflib::rdf_query("
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX gtpo: <https://rdf.guidetopharmacology.org/ns/gtpo#>

  SELECT ?targetLabel ?ligandLabel ?reference ?action ?medianAffinity ?affinityUnits
  WHERE {
  	?interaction gtpo:hasTarget ?target .
  	?interaction gtpo:hasLigand ?ligand .
  	?interaction gtpo:hasAffinity ?affinity .

  	?target gtpo:hasTargetFamily ?targetFamily .
  	?targetFamily gtpo:targetFamilyName \"Opioid receptors\" .

  	?target rdfs:label ?targetLabel .
  	?ligand rdfs:label ?ligandLabel .

  	?interaction gtpo:hasReference ?reference .
  	?interaction gtpo:hasAction ?action .

  	?affinity gtpo:hasMedianValue ?medianAffinity .
  	?affinity gtpo:hasUnits ?affinityUnits
  }")

