

library(BayesPharma)

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

date_code <- "20250822"

load(paste0(intermediate_path, "/data.Rdata"))

#########################
#GPTyS model individual #
#########################
assay_type = "GTPyS"
model_data <- data |>
    dplyr::group_by(date_code, assay, substance_id) |>
    dplyr::filter(dplyr::row_number() > 8) |>
    dplyr::ungroup() |>    
    dplyr::filter(assay == assay_type) |>
    dplyr::select(
        substance_id,
        log_dose = log10_dose,
        response = value_normalized)


model_1 <- BayesPharma::sigmoid_model(
    data = model_data,
    formula = BayesPharma::sigmoid_agonist_formula(
        response_units = "% Reference GTPyS Stimulation",
        predictor = 0 + substance_id),
    prior = BayesPharma::sigmoid_agonist_prior(
       ec50 = brms::prior(prior = normal(-8.5, 2.5), nlpar = "ec50"),
       hill = brms::prior(prior = normal(1, 1), nlpar = "hill", lb = -0.01),
       top = brms::prior(prior = normal(1, 0.5), nlpar = "top"),
       bottom = brms::prior(prior = normal(0, 0.5), nlpar = "bottom")),
    init = BayesPharma::sigmoid_agonist_init(
       ec50 = function() stats::runif(n = 1, min = -8.5 - 2, max = -8.5 + 2),
       hill = function() stats::runif(n = 1, min = 0.8, max = 1.2),
       top = function() stats::runif(n = 1, min = 0.8, max = 1.2),
       bottom = function() stats::runif(n = 1, min = -0.2, max = 0.2)),
    cores = 4,
    backend = "cmdstanr")


plot <- BayesPharma::plot_posterior_draws(model_1)

n_facets <- model_data |>
    dplyr::distinct(substance_id) |>
    nrow()
n_cols <- n_facets |> sqrt() |> ceiling()
n_rows <- ceiling(n_facets / n_cols)

ggplot2::ggsave(
    filename = paste0(output_path, "/GPTyS_model_1_posterior_draws_", date_code, ".pdf"),
    plot = plot,
    width = n_cols * 2,
    height = n_rows * 2,
    useDingbats = FALSE)



draws <- model_1 |>
    posterior::as_draws_df() |>
    dplyr::select(
        .chain, .iteration,
        b_ec50_substance_id30278,
        b_top_substance_id30278,
        b_ec50_substance_id30279,
        b_top_substance_id30279) |>
    dplyr::filter(.iteration > 2000) |>
    tidyr::pivot_longer(
        cols = tidyselect::starts_with("b_"),
        names_to = "parameter",
        values_to = "value") |>
    dplyr::mutate(
        parameter_type = parameter |>
            stringr::str_extract("b_([^_]+)", group = 1),
        substance_id = parameter |>
            stringr::str_extract("substance_id([^_]+)$", group = 1)) |>
    dplyr::select(-parameter) |>
    tidyr::pivot_wider(
        id_cols = c(".chain", ".iteration", "substance_id"),
        names_from = "parameter_type",
        values_from = "value")


plot <- ggplot2::ggplot(data = draws) +
    ggplot2::ggtitle("Bayesian Posterior Estimates") +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = ec50,
            y = top),
        alpha = 0.2,
        size = 0.4,
        shape = 16) +
    ggplot2::facet_wrap(
        facets = dplyr::vars(substance_id))
    
ggplot2::ggsave(
    filename = paste0(
        output_path,
        "/GTPyS_model_1_example_posterior_distributions_",
        date_code, ".pdf"),
    plot = plot,
    width = 5,
    height = 3,
    useDingbats = FALSE)


##############################
# GTPyS Unnormalized model_2 #
##############################
assay_type = "GTPyS"
model_data <- data |>
    dplyr::filter(assay == assay_type) |>    
    dplyr::group_by(date_code, assay, substance_id) |>
    dplyr::mutate(
        is_NC = dplyr::row_number() %in% c(5:8),
        is_PC = dplyr::row_number() %in% c(5:8),
        substance_id = dplyr::case_when(
            dplyr::row_number() %in% c(1:8) ~ NA,
            TRUE ~ substance_id)) |>
    dplyr::ungroup() |>
    dplyr::select(
        date_code,
        substance_id,
        log_dose = log10_dose,
        response = value)

model_formula <- brms::brmsformula(
    response ~ baseline + radioactivity * sigmoid(ec50, hill, top, bottom, log_dose))
    baseline ~ date_code,
    radioactivity ~ date_code,
    brms::bf(ec50 + hill + top + bottom ~ isPC * 1 + isNC * 0 + isTreatment * substance_id,
    nl = TRUE)
model_formula$bayes_pharma_info <- list(
    formula_type = "sigmoid_agonist",
    treatment_variable = "log_dose",
    treatment_units = "Log[Molar]",
    response_variable = "response",
    response_units = "GTPyS Stimulation")
class(model_formula) <- c("bpformula", class(model_formula))

model_prior = BayesPharma::sigmoid_agonist_prior(
    baseline = brms::prior(prior = normal(50, 100), nlpar = "baseline"),
    radioactivity = brms::prior(prior = normal(300, 100), nlpar = "radioactivity"),
    ec50 = brms::prior(prior = normal(-8.5, 2.5), nlpar = "ec50"),
    hill = brms::prior(prior = normal(1, 1), nlpar = "hill", lb = -0.01),
    top = brms::prior(prior = normal(1, 0.5), nlpar = "top"),
    bottom = brms::prior(prior = normal(0, 0.5), nlpar = "bottom"),
    

    
model_init = BayesPharma::sigmoid_agonist_init(
    baseline = function() stats::runif(n = 1, min = 50 - 20, max = 50 + 20),
    radioactivity = function() stats::runif(n = 1, min = 300 - 100, max = 300 + 100),
    ec50 = function() stats::runif(n = 1, min = -8.5 - 2, max = -8.5 + 2),
    hill = function() stats::runif(n = 1, min = 0.8, max = 1.2),
    top = function() stats::runif(n = 1, min = 0.8, max = 1.2),
    bottom = function() stats::runif(n = 1, min = -0.2, max = 0.2)),




model_2 <- BayesPharma::sigmoid_model(
    data = model_data,
    formula = model_formula,
    cores = 4,
    backend = "cmdstanr")

#############33


assay_type = "GTPyS"
model_data <- data |>
    dplyr::group_by(date_code, assay, substance_id) |>
    dplyr::mutate(
        isNC = dplyr::row_number() %in% c(1:4),
        isPC = dplyr::row_number() %in% c(5:8),
        isTreatment = dplyr::row_number() > 0) |>
    dplyr::ungroup() |>    
    dplyr::filter(assay == assay_type) |>
    dplyr::select(
        date_code,
        isPC, isNC, isTreatment,
        substance_id,
        log_dose = log10_dose,
        response = value_normalized)


model_formula <- brms::brmsformula(
    brms::bf(response | subset(isPC) ~ baseline + radioactivity * 0) +
    brms::bf(response | subset(isNC) ~ baseline + radioactivity * 1) +
    brms::bf(response | subset(isTreatment) ~
                 baseline + radioactivity * sigmoid(ec50, hill, top, bottom, log_dose)) +
    set_rescor(FALSE),
    baseline ~ date_code,
    radioactivity ~ date_code,
    ec50 + hill + top + bottom ~ substance_id,
    nl = TRUE)
model_formula$bayes_pharma_info <- list(
    formula_type = "sigmoid_agonist",
    treatment_variable = "log_dose",
    treatment_units = "Log[Molar]",
    response_variable = "response",
    response_units = "GTPyS Stimulation")
class(model_formula) <- c("bpformula", class(model_formula))


model_2 <- BayesPharma::sigmoid_model(
    data = model_data,
    formula = model_formula,
    prior = BayesPharma::sigmoid_agonist_prior(
        baseline = brms::prior(prior = normal(50, 100), nlpar = "baseline"),
        radioactivity = brms::prior(prior = normal(300, 100), nlpar = "radioactivity"),
        ec50 = brms::prior(prior = normal(-8.5, 2.5), nlpar = "ec50"),
        hill = brms::prior(prior = normal(1, 1), nlpar = "hill", lb = -0.01),
        top = brms::prior(prior = normal(1, 0.5), nlpar = "top"),
        bottom = brms::prior(prior = normal(0, 0.5), nlpar = "bottom")),
    init = BayesPharma::sigmoid_agonist_init(
        baseline = function() stats::runif(n = 1, min = 50 - 20, max = 50 + 20),
        radioactivity = function() stats::runif(n = 1, min = 300 - 100, max = 300 + 100),
        ec50 = function() stats::runif(n = 1, min = -8.5 - 2, max = -8.5 + 2),
        hill = function() stats::runif(n = 1, min = 0.8, max = 1.2),
        top = function() stats::runif(n = 1, min = 0.8, max = 1.2),
        bottom = function() stats::runif(n = 1, min = -0.2, max = 0.2)),
    cores = 4,
    backend = "cmdstanr")


