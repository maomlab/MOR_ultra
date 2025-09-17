


library(pzfx)
library(GeomIndicator)

data_path <- "data/AnandLab/Nitazene_20250818"

intermediate_path <- "intermediate/AnandLab_Nitazene_20250818"
if (!dir.exists(intermediate_path)) {
    cat("Creating output_path '", intermediate_path, "'\n", sep = "")
    dir.create(intermediate_path)
}


output_path <- "product/AnandLab_Nitazene_20250818"
if (!dir.exists(output_path)) {
    cat("Creating output_path '", output_path, "'\n", sep = "")
    dir.create(output_path)
}

dataset_name <- "Nitazene Analogs 20250818"

date_code <- "20250822"


data <- data_path |>
    list.files(pattern = "*pzfx", full.names = TRUE, recursive = TRUE) |>
    purrr::map_dfr(function(pzfx_fname){
        cat("Reading pzfx file: '", pzfx_fname, "'\n", sep = "")
        pzfx::pzfx_tables(pzfx_fname) |>
        purrr::map_dfr(function(table_name) {
            if (table_name %in% c("std ")) {
                cat("  Skipping table '", table_name, "'\n", sep = "")
                
            } else {
                cat("  Reading table: '", table_name, "'\n", sep = "")
                data <- pzfx::read_pzfx(
                    path = pzfx_fname,
                    table = table_name) |>
                    dplyr::mutate(
                        pzfx_fname = pzfx_fname,
                        table_name = table_name,
                        .before = 1)
                if (!"ROWTITLE" %in% names(data)) {
                    data <- data |> dplyr::mutate(ROWTITLE = NA)
                }
                data <-  data |>
                    tidyr::pivot_longer(
                        cols = -c("pzfx_fname", "table_name", "ROWTITLE", "log10_dose"),
                        names_to = "substance_id",
                        values_to = "value")
            }
        })
    })

data <- data |>
    tidyr::separate_wider_delim(
        cols = substance_id,
        names = c("substance_id", "replica"), delim = "_")

data <- data |> 
    dplyr::mutate(
        assay = pzfx_fname |> stringr::str_extract(".*[/]([^/]+)[/][^/]+", group = 1),
        date_code = pzfx_fname |>
            stringr::str_extract("(([0-9]+[.][0-9]+[.][0-9]+|[0-9]+[-][0-9]+[-][0-9]+))") |>
            stringr::str_replace_all("[-]", "."))

data <- data |>
    dplyr::group_by(pzfx_fname, table_name, date_code, assay) |>
    dplyr::mutate(
        bottom_mean = value[1:4] |> mean(na.rm = TRUE),
        top_mean = value[5:8] |> mean(na.rm = TRUE),
        value_normalized = (value - bottom_mean) / (top_mean - bottom_mean)) |>
    dplyr::ungroup()
save(data, file = paste0(intermediate_path, "/data.Rdata"))


plot_data <- data |>
    dplyr::group_by(date_code, assay, substance_id) |>
    dplyr::filter(dplyr::row_number() > 8) |>
    dplyr::ungroup() |>
    dplyr::left_join(
        data |>
            dplyr::group_by(date_code, assay) |>
            dplyr::distinct(substance_id) |>
            dplyr::mutate(substance_index = dplyr::row_number() |> as.factor()) |>
            dplyr::ungroup(),
        by = c("date_code", "assay", "substance_id"))

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
        facets = dplyr::vars(date_code)),
    ggplot2::theme(legend.position="none"))


assay_type <- "Binding"
plot_data_assay <- plot_data |> dplyr::filter(assay == assay_type)
n_facets <- plot_data_assay |>
    dplyr::distinct(date_code) |>
    nrow()
n_cols <- n_facets |> sqrt() |> ceiling()
n_rows <- ceiling(n_facets / n_cols)

plot <- ggplot2::ggplot(data = plot_data_assay) +
    common_layers +
    ggplot2::ggtitle(
        label = paste0(dataset_name),
        subtitle = assay_type) +
    ggplot2::scale_y_continuous(
        expression("% Displacement of Reference Ligand"),
        breaks = c(0, 0.5, 1, 1.5),
        labels = c("0%", "50%", "100%", "150%"))

ggplot2::ggsave(
    filename = paste0(output_path, "/data_summary_", assay_type, "_", date_code, ".pdf"),
    width = n_cols * 2,
    height = n_rows * 2,
    useDingbats = FALSE)
plot_data |>
    readr::write_tsv(
        paste0(output_path, "/data_summary_", assay_type, "_", date_code, ".tsv"))



assay_type <- "GTPyS"
plot_data_assay <- plot_data |> dplyr::filter(assay == assay_type)
n_facets <- plot_data_assay |>
    dplyr::distinct(date_code) |>
    nrow()
n_cols <- n_facets |> sqrt() |> ceiling()
n_rows <- ceiling(n_facets / n_cols)

plot <- ggplot2::ggplot(data = plot_data_assay) +
    common_layers +
    ggplot2::ggtitle(
        label = paste0(dataset_name),
        subtitle = assay_type) +
    ggplot2::scale_y_continuous(
        expression("Scaled ["^"35"*"S]GTP"*gamma*"S"),
        breaks = c(0, 0.5, 1, 1.5),
        labels = c("0%", "50%", "100%", "150%"))

ggplot2::ggsave(
    filename = paste0(output_path, "/data_summary_", assay_type, "_", date_code, ".pdf"),
    width = n_cols * 2,
    height = n_rows * 2,
    useDingbats = FALSE)
plot_data |>
    readr::write_tsv(
        paste0(output_path, "/data_summary_", assay_type, "_", date_code, ".tsv"))

