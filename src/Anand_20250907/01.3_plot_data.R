
library(GeomIndicator)

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

load(paste0(intermediate_path, "/MOR_binding.Rdata"))


#########
protocol_id_label <- "MOR Binding Ki Data"
researcher_label <- "Allyson"
date_code_label <- "3H 2022/August/8.1.22"
assay_id_label <- "9003717, 9003719, 9003720, 9003721 hDOR Binding 8.1.22"
taget_label <- "OPRD_HUMAN"

plot_data <- MOR_binding |>
  dplyr::filter(treatment_type == "treatment") |>
  dplyr::filter(
    protocol_id == protocol_id_label,
    researcher == researcher_label,
    date_code == date_code_label,
    assay_id == assay_id_label,
    target == target_label)

plot_data <- plot_data |>
  dplyr::left_join(
    plot_data |>
      dplyr::group_by(run_id) |>
      dplyr::distinct(substance_id) |>
      dplyr::mutate(substance_index = dplyr::row_number() |> as.factor()) |>
      dplyr::ungroup(),
    by = c("run_id", "substance_id"))


common_layers <- list(
  ggplot2::theme_bw(),
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.8,
    color = "grey90"),
  ggplot2::geom_hline(
    yintercept = 1,
    linewidth = 0.8,
    color = "grey90"),
  ggplot2::geom_line(
    mapping = ggplot2::aes(
      x = log10_dose,
      y = value_normalized,
      color = substance_index,
      group = paste0(substance_index, replica, sep = "_"))),
  GeomIndicator::geom_indicator(
    mapping = ggplot2::aes(
      indicator = substance_id,
      color = substance_index,
      group = substance_index),
    size = 2,
    ypos = "top",
    xpos = "right"),
  ggplot2::scale_x_continuous(expression("Log"[10]*"[Dose]")),
  ggplot2::facet_wrap(
    facets = dplyr::vars(run_id)),
  ggplot2::theme(legend.position="none"))


n_facets <- plot_data |> dplyr::distinct(date_code) |> nrow()
n_cols <- n_facets |> sqrt() |> ceiling()
n_rows <- ceiling(n_facets / n_cols)

plot <- ggplot2::ggplot(data = plot_data) +
  common_layers +
  ggplot2::ggtitle(
    label = paste0(dataset_name),
    subtitle = paste0(
      protocol_id_label, "\n",
      researcher_label, "\n",
      date_code_label, "\n",
      assay_id_label)) +
  ggplot2::scale_y_continuous(
    expression("% Displacement of Reference Ligand"),
    breaks = c(0, 0.5, 1, 1.5),
    labels = c("0%", "50%", "100%", "150%"))

ggplot2::ggsave(
  filename = paste0(
    output_path,
    "/data_summary_",
    protocol_id_label, "_",
    researcher_label, "_",
    date_code_label |> stringr::str_replace_all("/", "|"), "_",
    assay_id_label, "_",
    date_code, ".pdf"),
  width = n_cols * 3,
  height = n_rows * 3,
  useDingbats = FALSE)
plot_data |>
  readr::write_tsv(
    file = paste0(
      output_path,
      "/data_summary_",
      protocol_id_label, "_",
      researcher_label, "_",
      date_code_label |> stringr::str_replace_all("/", "|"), "_",
      assay_id_label, "_",
      date_code, ".tsv"))



###########

common_layers <- list(
  ggplot2::theme_bw(),
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.8,
    color = "grey90"),
  ggplot2::geom_hline(
    yintercept = 1,
    linewidth = 0.8,
    color = "grey90"),
  ggplot2::geom_line(
    mapping = ggplot2::aes(
      x = log10_dose,
      y = value_normalized,
      color = substance_index,
      group = paste0(substance_index, replica, sep = "_"))),
  GeomIndicator::geom_indicator(
    mapping = ggplot2::aes(
      indicator = substance_id,
      color = substance_index,
      group = substance_index),
    size = 2,
    ypos = "top",
    xpos = "right"),
  ggplot2::scale_x_continuous(expression("Log"[10]*"[Dose]")),
  ggplot2::facet_wrap(
    facets = dplyr::vars(date_code, assay_id),
    scales = "free_y"),
  ggplot2::theme(legend.position="none"))


MOR_binding |>
  dplyr::group_by(protocol_id, researcher, target) |>
  dplyr::do({
    data <- .
    protocol_id_label <- data$protocol_id[1]
    researcher_label <- data$researcher[1]
    target_label <- data$target[1]

  cat("Plotting data for ", protocol_id_label,
      " and ", researcher_label, " at ", target_label, "\n", sep = "")
  plot_data <- MOR_binding |>
    dplyr::filter(treatment_type == "treatment") |>
    dplyr::filter(
      protocol_id == protocol_id_label,
      researcher == researcher_label,
      target == target_label)

  plot_data <- plot_data |>
    dplyr::left_join(
      plot_data |>
        dplyr::group_by(date_code, assay_id) |>
        dplyr::distinct(substance_id) |>
        dplyr::mutate(substance_index = dplyr::row_number() |> as.factor()) |>
        dplyr::ungroup(),
      by = c("date_code", "assay_id", "substance_id"))

  n_facets <- plot_data |> dplyr::distinct(date_code) |> nrow()
  n_cols <- n_facets |> sqrt() |> ceiling()
  n_rows <- ceiling(n_facets / n_cols)

  plot <- ggplot2::ggplot(data = plot_data) +
    common_layers +
    ggplot2::ggtitle(
      label = paste0(dataset_name),
      subtitle = paste0(protocol_id_label, " ", researcher_label, " at ", target_label)) +
    ggplot2::scale_y_continuous(
      expression("% Displacement of Reference Ligand"),
      breaks = c(0, 0.5, 1, 1.5),
      labels = c("0%", "50%", "100%", "150%"))

  ggplot2::ggsave(
    filename = paste0(
      output_path,
      "/data_summary_",
      protocol_id_label, "_",
      researcher_label, "_",
      target_label, "_",
      date_code, ".pdf"),
    width = n_cols * 2,
    height = n_rows * 2,
    useDingbats = FALSE)
  plot_data |>
    readr::write_tsv(
      file = paste0(
        output_path,
        "/data_summary_",
        protocol_id_label, "_",
        researcher_label, "_",
        target_label, "_",
        date_code, ".tsv"))

  data.frame()
})


##################

common_layers <- list(
  ggplot2::theme_bw(),
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.8,
    color = "grey90"),
  ggplot2::geom_hline(
    yintercept = 1,
    linewidth = 0.8,
    color = "grey90"),
  ggplot2::geom_line(
    mapping = ggplot2::aes(
      x = log10_dose,
      y = value_normalized,
      color = researcher_index,
      group = paste0(researcher_index, date_code, assay_id, run_id, replica, sep = "_")),
    alpha = 0.8),
  GeomIndicator::geom_indicator(
    mapping = ggplot2::aes(
      indicator = researcher,
      color = researcher_index,
      group = researcher_index),
    size = 2,
    ypos = "top",
    xpos = "right"),
  ggplot2::scale_x_continuous(expression("Log"[10]*"[Dose]")),
  ggplot2::facet_wrap(
    facets = dplyr::vars(substance_id)),
  ggplot2::theme(legend.position="none"))




plot_data <- MOR_binding |>
  dplyr::filter(treatment_type == "treatment")

plot_data <- plot_data |>
  dplyr::left_join(
    plot_data |>
      dplyr::group_by(substance_id) |>
      dplyr::distinct(researcher) |>
      dplyr::mutate(researcher_index = dplyr::row_number() |> as.factor()) |>
      dplyr::ungroup(),
    by = c("substance_id", "researcher"))

n_facets <- plot_data |> dplyr::distinct(substance_id) |> nrow()
n_cols <- n_facets |> sqrt() |> ceiling()
n_rows <- ceiling(n_facets / n_cols)

plot <- ggplot2::ggplot(data = plot_data) +
  common_layers +
  ggplot2::ggtitle(
    label = paste0(dataset_name)) +
  ggplot2::scale_y_continuous(
    expression("% Displacement of Reference Ligand"),
    breaks = c(0, 0.5, 1, 1.5),
    labels = c("0%", "50%", "100%", "150%"),
    limits = c(-.3, 2))

ggplot2::ggsave(
  filename = paste0(
    output_path,
    "/data_summary_", protocol_id_label, "_by_substance_",
    date_code, ".pdf"),
  width = n_cols * 2,
  height = n_rows * 2,
  useDingbats = FALSE)
plot_data |>
  readr::write_tsv(
    file = paste0(
      output_path,
      "/data_summary_", protocol_id_label, "_by_substance_",
      date_code, ".tsv"))

