
activity_data <- readxl::read_excel(
  "data/MOR_ligands_20250531.xlsx",
  skip = 2)

rmsds_long <- readr::read_tsv(
  "intermediate/rmsd_2.tsv",
  show_col_types = FALSE)

rmsds_wide <- rmsds_long |>
  dplyr::mutate(ref_pdb_id = ref_path |> stringr::str_extract("([^/]+).cif", group = 1)) |>
  dplyr::select(-ref_path) |>
  tidyr::pivot_wider(
    id_cols = target_path,
    names_from = ref_pdb_id,
    values_from = rmsd)

structure_info <- rmsds_wide |>
  dplyr::select(target_path) |>
  tidyr::extract(
    target_path,
    into = c("source", "type", "params", "date_code", "ligand", "index"),
    regex = "^([^_]+)_([^_]+)_([^_]+)_(\\d{8})/(.+)/pred\\.rank_(\\d+)\\.cif$",
    remove = FALSE
  ) |>
  dplyr::mutate(index = as.integer(index))

umap_embedding <- uwot::umap(
  X = rmsds_wide |> dplyr::select(-target_path),
  n_neighbors = 20,
  min_dist = .5,
  verbose = TRUE)
    
features <- structure_info |>
  dplyr::mutate(
    UMAP1 = umap_embedding[,1],
    UMAP2 = umap_embedding[,2]) |>
  dplyr::bind_cols(rmsds_wide |> dplyr::select(-target_path)) |>
  dplyr::left_join(
    activity_data,
    by = c("ligand" = "Ligand"))

save(features, file = "intermediate/rmsd_features.Rdata")
    
plot <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::coord_fixed() +
  ggplot2::geom_point(
    data = features,
    mapping = ggplot2::aes(
      x = UMAP1,
      y = UMAP2,
      color = `MOR Emax`),
    size = 1.2,
    alpha = 0.9,
    shape = 16) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  filename = "product/02_rmsd_Emax_embedding.pdf",
  width = 5,
  height = 5)

############################


fit <- lm(
  formula = `MOR Emax` ~ `8E0G` + `9MQH` + `9BJK` +   `5C1M` + `4DKL` + `9MQJ` + `7UL4` + `8QOT` + `9MQI` + `8EF5` + `8EF6` + `8EFB` + `7T2H` + `8EFL` + `8EFO` + `7SBF` + `8EFQ` + `8F7R` + `7T2G` + `8K9L` + `6DDE` + `7SCG` + `8F7Q` + `6DDF` + `8K9K` + `9BQJ` + `7U2K` + `7U2L`,
  data = features |> dplyr::filter(params == "Wtemplate"))

##################3

plot <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::coord_fixed() +
  ggplot2::geom_point(
    data = features,
    mapping = ggplot2::aes(
      x = `8QOT`,
      y = `7SBF`,
      color = `MOR Emax`),
    size = 1.2,
    alpha = 0.9,
    shape = 16) +
  ggplot2::facet_wrap(facets = dplyr::vars(params)) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  filename = "product/02_rmsd_8QOT_7SBF_20250603.pdf",
  width = 8,
  height = 4)


########################3

plot <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    data = features,
    mapping = ggplot2::aes(
      x = `7SBF`,
      y = `MOR Emax`),
    size = 1.2,
    alpha = 0.9,
    shape = 16) +
  ggplot2::geom_smooth(
    data = features,
    mapping = ggplot2::aes(
      x = `7SBF`,
      y = `MOR Emax`),
    method = "lm",
    formula = y ~ x,
    se = FALSE) +
  ggplot2::facet_wrap(facets = dplyr::vars(params)) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  filename = "product/02_rmsd_Emax_7SBF_20250603.pdf",
  width = 8,
  height = 4)


####################

plot <- ggplot2::ggplot(
  data = features |>
    dplyr::filter(index = 0) |>
    dplyr::filter(params == "Wtemplate") |>
    dplyr::mutate(
      pred_Emax = predict(fit, data = features))) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = pred_Emax,
      y = `MOR Emax`,
      color = index),
    size = 1.2,
    alpha = 0.9,
    shape = 16) +
  ggplot2::geom_smooth(
    mapping = ggplot2::aes(
      x = pred_Emax,
      y = `MOR Emax`),
    method = "lm",
    formula = y ~ x,
    se = FALSE) +
  ggplot2::facet_wrap(facets = dplyr::vars(params)) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_x_continuous("Predicted Emax")

ggplot2::ggsave(
  filename = "product/02_rmsd_Emax_pred_20250603.pdf",
  width = 8,
  height = 4)
