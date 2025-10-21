
intermediate_path <- "intermediate/Anand_20250907"
if (!dir.exists(intermediate_path)) {
  cat("Creating output_path '", intermediate_path, "'\n", sep = "")
  dir.create(intermediate_path)
}


output_path <- "product/Anand_20250907/model_fits"
if (!dir.exists(output_path)) {
  cat("Creating output_path '", output_path, "'\n", sep = "")
  dir.create(output_path)
}

dataset_name <- "Anand Data 20250907"

date_code <- "20250907"

load(paste0(intermediate_path, "/MOR_binding.Rdata"))

###########################################
# 4 substances tested by Allyson 8/1/2022 #
###########################################

protocol_id_label <- "MOR Binding Ki Data"
researcher_label <- "Allyson"
date_code_label <- "3H 2022/August/8.1.22"
assay_id_label <- "9003717, 9003719, 9003720, 9003721 hDOR Binding 8.1.22"

model_data <- MOR_binding |>
  dplyr::filter(treatment_type == "treatment") |>
  dplyr::filter(
    protocol_id == protocol_id_label,
    researcher == researcher_label,
    date_code == date_code_label,
    assay_id == assay_id_label) |>
    dplyr::select(
        substance_id,
        log_dose = log10_dose,
        response = value_normalized)


model <- BayesPharma::sigmoid_model(
    data = model_data,
    formula = BayesPharma::sigmoid_antagonist_formula(
        response_units = "% Displacement of Reference Ligand",
        predictor = 0 + substance_id),
    prior = BayesPharma::sigmoid_antagonist_prior(
       ic50 = brms::prior(prior = normal(-8.5, 3.5), nlpar = "ic50"),
       hill = brms::prior(prior = normal(-1, 1), nlpar = "hill", ub = -0.01),
       top = brms::prior(prior = normal(1, 0.5), nlpar = "top"),
       bottom = brms::prior(prior = normal(0, 0.5), nlpar = "bottom")),
    init = BayesPharma::sigmoid_antagonist_init(
       ic50 = function() stats::runif(n = 1, min = -8.5 - 3, max = -8.5 + 3),
       hill = function() stats::runif(n = 1, min = -1.2, max = -0.8),
       top = function() stats::runif(n = 1, min = 0.8, max = 1.2),
       bottom = function() stats::runif(n = 1, min = -0.2, max = 0.2)),
    cores = 4,
    backend = "cmdstanr")


plot <- BayesPharma::plot_posterior_draws(model)

n_facets <- model_data |>
    dplyr::distinct(substance_id) |>
    nrow()
n_cols <- n_facets |> sqrt() |> ceiling()
n_rows <- ceiling(n_facets / n_cols)

ggplot2::ggsave(
    filename = paste0(
        output_path, "/", protocol_id_label, "_",
        researcher_label, "_",
        date_code |> stringr::str_replace_all("/", "|"), "_",
        assay_id_label, "_",
        date_code, ".pdf"),
    plot = plot,
    width = n_cols * 2,
    height = n_rows * 2,
    useDingbats = FALSE)


###########################################
# 4 substances tested by Allyson 8/1/2022 #
###########################################

protocol_id_label <- "MOR Binding Ki Data"
researcher_label <- "Allyson"
date_code_label <- "3H 2022/August/8.1.22"
assay_id_label <- "9003717, 9003719, 9003720, 9003721 hDOR Binding 8.1.22"

model_data <- MOR_binding |>
  dplyr::filter(treatment_type == "treatment")

model_data <- model_data |>
  dplyr::semi_join(
    model_data |>
      dplyr::distinct(substance_id, researcher) |>
      dplyr::count(substance_id) |>
      dplyr::filter(n >= 4),
    by = c("substance_id")) |>
  dplyr::select(
    researcher,
    substance_id,
    log_dose = log10_dose,
    response = value_normalized)

models <- model_data |>
  dplyr::filter(substance_id == "20018") |>
  dplyr::group_by(substance_id) |>
  dplyr::do({
    data <- .
    substance_id_label <- data$substance_id[1]
    cat("Fitting model for substance ", substance_id_label, "\n", sep = "")

    model <- BayesPharma::sigmoid_model(
      data = model_data,
      formula = BayesPharma::sigmoid_antagonist_formula(
        response_units = "% Displacement of Reference Ligand",
        predictor = 0 + (1|researcher)),
      prior = BayesPharma::sigmoid_antagonist_prior(
        ic50 = brms::prior(prior = normal(-8.5, 3.5), nlpar = "ic50"),
        hill = brms::prior(prior = normal(-1, 1), nlpar = "hill", ub = -0.01),
        top = brms::prior(prior = normal(1, 0.5), nlpar = "top"),
        bottom = brms::prior(prior = normal(0, 0.5), nlpar = "bottom")),
      init = BayesPharma::sigmoid_antagonist_init(
        ic50 = function() stats::runif(n = 1, min = -8.5 - 3, max = -8.5 + 3),
        hill = function() stats::runif(n = 1, min = -1.2, max = -0.8),
        top = function() stats::runif(n = 1, min = 0.8, max = 1.2),
        bottom = function() stats::runif(n = 1, min = -0.2, max = 0.2)),
      cores = 4,
      backend = "cmdstanr")

    plot <- BayesPharma::plot_posterior_draws(model)

    n_facets <- model_data |>
      dplyr::distinct(substance_id) |>
      nrow()
    n_cols <- n_facets |> sqrt() |> ceiling()
    n_rows <- ceiling(n_facets / n_cols)

    ggplot2::ggsave(
      filename = paste0(
        output_path, "/", protocol_id_label, "_",
        substance_id_label, "_",
        date_code, ".pdf"),
      plot = plot,
      width = n_cols * 2,
      height = n_rows * 2,
      useDingbats = FALSE)
  tibble::tibble(
    substance_id = substance_id_label,
    model = model)
})


##################
# All Substances #
##################


protocol_id_label <- "MOR Binding Ki Data"

model_data <- MOR_binding |>
  dplyr::filter(treatment_type == "treatment") |>
  dplyr::filter(
    protocol_id == protocol_id_label) |>
  dplyr::select(
    substance_id,
    log_dose = log10_dose,
    response = value_normalized)


model <- BayesPharma::sigmoid_model(
  data = model_data,
  formula = BayesPharma::sigmoid_antagonist_formula(
    response_units = "% Displacement of Reference Ligand",
    predictor = 0 + substance_id),
  prior = BayesPharma::sigmoid_antagonist_prior(
    ic50 = brms::prior(prior = normal(-8.5, 3.5), nlpar = "ic50"),
    hill = brms::prior(prior = normal(-1, 1), nlpar = "hill", ub = -0.01),
    top = brms::prior(prior = normal(1, 0.5), nlpar = "top"),
    bottom = brms::prior(prior = normal(0, 0.5), nlpar = "bottom")),
  init = BayesPharma::sigmoid_antagonist_init(
    ic50 = function() stats::runif(n = 1, min = -8.5 - 3, max = -8.5 + 3),
    hill = function() stats::runif(n = 1, min = -1.2, max = -0.8),
    top = function() stats::runif(n = 1, min = 0.8, max = 1.2),
    bottom = function() stats::runif(n = 1, min = -0.2, max = 0.2)),
  cores = 4,
  backend = "cmdstanr")


plot <- BayesPharma::plot_posterior_draws(model)

n_facets <- model_data |>
  dplyr::distinct(substance_id) |>
  nrow()
n_cols <- n_facets |> sqrt() |> ceiling()
n_rows <- ceiling(n_facets / n_cols)

ggplot2::ggsave(
  filename = paste0(
    output_path, "/", protocol_id_label, "_",
    researcher_label, "_",
    date_code |> stringr::str_replace_all("/", "|"), "_",
    assay_id_label, "_",
    date_code, ".pdf"),
  plot = plot,
  width = n_cols * 2,
  height = n_rows * 2,
  useDingbats = FALSE)

