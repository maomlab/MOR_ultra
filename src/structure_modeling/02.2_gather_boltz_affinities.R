




mor_fasta <- seqinr::read.fasta(
    file = "data/OPRM_HUMAN.fasta",
    seqtype = "AA")

activity_data <- readxl::read_excel(
  "data/MOR_ligands_20250606.xlsx",
  skip = 2) |>
  dplyr::filter(!(Ligand |> stringr::str_detect("^[*]")))



##############################
# Chai-1 Affinity Predictions #
##############################
prediction_path <-
    "data/chai1_predictions_Wtemplate_20250602"

chai1_prediction_data <- activity_data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand_data <- .

        z <- data.frame(prediction_rank = 0:4) |>
            dplyr::rowwise() |>
            dplyr::do({
                prediction_rank <- .$prediction_rank                            
                affinities_path <- paste0(
                    prediction_path, "/", ligand_data$Ligand, "/",
                    "scores.rank_", prediction_rank, ".json")
                affinities <- rjson::fromJSON(file = affinities_path)
                ligand_data |>
                    data.frame() |>
                    dplyr::mutate(
                        prediction_rank = prediction_rank,
                        aggregate_score = affinities$aggregate_score,
                        ptm = affinities$ptm,
                        iptm = affinities$iptm,
                        receptor_ptm = affinities$per_chain_ptm[[1]][1],
                        ligand_ptm = affinities$per_chain_ptm[[1]][2])
            })
    }) |>
    dplyr::ungroup()

chai1_prediction_data |>
    readr::write_tsv("product/chai1_predictions_Wtemplate_20250602.tsv")

chai1_prediction_data_top <- chai1_prediction_data |>
    dplyr::filter(prediction_rank == 0)
cor(
    x = chai1_prediction_data_top$MOR.Emax,
    y = chai1_prediction_data_top$aggregate_score,
    method = "spearman")




##############################
# Boltz Affinity Predictions #
##############################
prediction_path <-
    "intermediate/boltz2_prediction_20250606/OPRM_HUMAN"

boltz_prediction_data <- activity_data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand_data <- .
        results_path <- paste0(
            prediction_path, "/", ligand_data$Ligand, "/boltz_results_", ligand_data$Ligand, "/",
            "predictions/", ligand_data$Ligand)
        confidences_path <- paste0(results_path, "/confidence_", ligand_data$Ligand, "_model_0.json")
        confidences <- rjson::fromJSON(file = confidences_path)

        affinities_path <- paste0(results_path, "/affinity_", ligand_data$Ligand, ".json")
        affinities <- rjson::fromJSON(file = affinities_path)

        ligand_data |>
            data.frame() |>
            dplyr::mutate(
                prediction_rank = 0,
                confidence_score = confidences$confidence_score,
                ptm = confidences$ptm,
                iptm = confidences$iptm,
                ligand_iptm = confidences$ligand_iptm,
                protein_iptm = confidences$protein_iptm,
                complex_plddt = confidences$complex_plddt,
                complex_iplddt = confidences$complex_iplddt,
                complex_pde = confidences$complex_pde,
                complex_ipde = confidences$complex_ipde,
                receptor_ptm = confidences$chains_ptm$`0`,
                ligand_ptm = confidences$chains_ptm$`1`,
                affinity_pred_value = affinities$affinity_pred_value,
                affinity_probability_binary = confidences$affinity_probability_binary,
                affinity_pred_value1 = confidences$affinity_pred_value1,
                affinity_probability_binary1 = confidences$affinity_probability_binary1,
                affinity_pred_value2 = confidences$affinity_pred_value2,
                affinity_probability_binary2 = confidences$affinity_probability_binary2)
    }) |>
    dplyr::ungroup()

boltz_prediction_data |>
    readr::write_tsv("product/boltz_affinity_predictions_20250606.tsv")



##############################
# Boltz Affinity Predictions #
# Extra1                     #
##############################
prediction_path <-
    "intermediate/boltz2_prediction_extra1_20250606/OPRM_HUMAN"

boltz_prediction_data <- activity_data |>
    dplyr::rowwise() |>
    dplyr::do({
        ligand_data <- .
        results_path <- paste0(
            prediction_path, "/boltz_results_inputs/",
            "predictions/", ligand_data$Ligand)
        affinities_path <- paste0(results_path, "/affinity_", ligand_data$Ligand, ".json")
        affinities <- rjson::fromJSON(file = affinities_path)

        z <- data.frame(prediction_rank = 0:199) |>
            dplyr::rowwise() |>
            dplyr::do({
                prediction_rank <- .$prediction_rank                            

                confidences_path <- paste0(results_path, "/confidence_", ligand_data$Ligand, "_model_", prediction_rank, ".json")
                confidences <- rjson::fromJSON(file = confidences_path)

                ligand_data |>
                    data.frame() |>
                    dplyr::mutate(
                        prediction_rank = prediction_rank,
                        confidence_score = confidences$confidence_score,
                        ptm = confidences$ptm,
                        iptm = confidences$iptm,
                        ligand_iptm = confidences$ligand_iptm,
                        protein_iptm = confidences$protein_iptm,
                        complex_plddt = confidences$complex_plddt,
                        complex_iplddt = confidences$complex_iplddt,
                        complex_pde = confidences$complex_pde,
                        complex_ipde = confidences$complex_ipde,
                        receptor_ptm = confidences$chains_ptm$`0`,
                        ligand_ptm = confidences$chains_ptm$`1`,
                        affinity_pred_value = affinities$affinity_pred_value,
                        affinity_probability_binary = confidences$affinity_probability_binary,
                        affinity_pred_value1 = confidences$affinity_pred_value1,
                        affinity_probability_binary1 = confidences$affinity_probability_binary1,
                        affinity_pred_value2 = confidences$affinity_pred_value2,
                        affinity_probability_binary2 = confidences$affinity_probability_binary2)
                })
    }) |>
    dplyr::ungroup()

boltz_prediction_data |>
    readr::write_tsv("product/boltz_affinity_predictions_20250606.tsv")

cor(
    x = prediction_data$MOR.Emax,
    y = prediction_data$confidence_score,
    method = "spearman")



z <- chai1_prediction_data_top |>
    dplyr::select(
        Ligand, MOR.Emax, chai1_score = aggregate_score) |>
    dplyr::left_join(
        boltz_prediction_data |>
        dplyr::select(
            Ligand,
            boltz_score = confidence_score),
        by = "Ligand")
cor(
    x = z$chai1_score,
    y = z$MOR.Emax,
    method = "spearman")



#############################33
# Structure features

structure_features <- readr::read_tsv(
    "product/structure_features_20250608.tsv",
    show_col_types = FALSE)


boltz_prediction_data <- boltz_prediction_data |>
    dplyr::left_join(
        structure_features |>
        dplyr::filter(source == "boltz2_prediction_extra1_20250606/OPRM_HUMAN") |>
        dplyr::transmute(
            Ligand = structure_id,
            prediction_rank = structure_index |>
                stringr::str_replace(".cif", "") |>
                as.numeric(),
            openness),
        by = c("Ligand", "prediction_rank"))


plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        data = boltz_prediction_data,
        mapping = ggplot2::aes(
            x = confidence_score,
            y = openness,
            color = MOR.Emax)) +
    ggplot2::facet_wrap(
        facets = dplyr::vars(Ligand))

ggplot2::ggsave(
    filename = "product/boltz_extra1_openness_by_confidence_20250608.pdf",
    plot = plot,
    width = 10,
    height = 10,
    useDingbats = FALSE)



plot_data <- boltz_prediction_data |>
    dplyr::mutate(
        open_state = ifelse(openness > 12.8, "open", "closed")) |>
    dplyr::group_by(Ligand, open_state) |>
    dplyr::summarize(
        mean_confidence = mean(confidence_score),
        MOR.Emax = MOR.Emax[1],
        .groups = "drop") |>
    tidyr::pivot_wider(
        id_cols = c(Ligand, MOR.Emax),
        names_from = open_state,
        values_from = mean_confidence) |>
    dplyr::mutate(
        open_close_diff = open - closed)

plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        data = plot_data,
        mapping = ggplot2::aes(
            x = open_close_diff,
            y = MOR.Emax))
ggplot2::ggsave(
    filename = "product/boltz_extra1_Emax_by_open_close_confidence_diff_20250608.pdf",
    plot = plot,
    width = 5,
    height = 5,
    useDingbats = FALSE)




plot_data <- boltz_prediction_data |>
    dplyr::group_by(Ligand) |>
    dplyr::arrange(dplyr::desc(confidence_score)) |>
    dplyr::filter(dplyr::row_number() < 100) |>
    dplyr::ungroup() |>
    dplyr::mutate(
        open_state = ifelse(openness > 12.8, "open", "closed")) |>
    dplyr::group_by(Ligand, open_state) |>
    dplyr::summarize(
        n_state = dplyr::n(),
        MOR.Emax = MOR.Emax[1],
        .groups = "drop") |>
    tidyr::pivot_wider(
        id_cols = c(Ligand, MOR.Emax),
        names_from = open_state,
        values_from = n_state) |>
    dplyr::mutate(
        open_close_diff = open - closed)

plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        data = plot_data,
        mapping = ggplot2::aes(
            x = open_close_diff,
            y = MOR.Emax))
ggplot2::ggsave(
    filename = "product/boltz_extra1_Emax_by_open_close_count_diff_20250608.pdf",
    plot = plot,
    width = 5,
    height = 5,
    useDingbats = FALSE)



plot_data <- plot_data |>
    dplyr::mutate(
        class = dplyr::case_when(
            Ligand
