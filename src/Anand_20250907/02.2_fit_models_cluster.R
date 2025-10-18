

---
title: "02.1 fit models cluster"
output: html_document
---
```{r load-packages}
library(future)
library(future.batchtools)
library(BayesPharma)
library(cmdstanr)
```

```{r prepare-future-plan}
# submitting jobs will use the SLURM queue
future.batchtools::batchtools_slurm(
    template = here::here("src/hpc/batchtools.greatlakes.tmpl"),
    resources = list(
        ncpus = 4,
        memory = 5000,
        walltime = 1000,
        account = "maom99",
        partition = "standard")) |>
    future::plan()
```

```{r run-params}
intermediate_path <- "intermediate/Anand_20250907"
if (!dir.exists(intermediate_path)) {
  cat("Creating output_path '", intermediate_path, "'\n", sep = "")
  dir.create(intermediate_path, recursive = TRUE)
}


output_path <- "product/Anand_20250907"
if (!dir.exists(output_path)) {
  cat("Creating output_path '", output_path, "'\n", sep = "")
  dir.create(output_path, recursive = TRUE)
}

date_code <- "20250922"
```

```{r load-data}
load(paste0(intermediate_path, "/MOR_binding.Rdata"))
```

```{r fit-model-function}
fit_model <- function(data, output_path) {
  model_data <- data |>
    dplyr::filter(treatment_type == "treatment") |>
    dplyr::select(
      substance_id,
      log_dose = log10_dose,
      response = value_normalized)

  cat(
    "Fitting model for:\n",
    "  protocol_id: ", data$protocol_id[1], "\n",
    "  researcher:  ", data$researcher[1], "\n",
    "  date_code:   ", data$date_code[1], "\n",
    "  assay_id:    ", data$assay_id[1], "\n",
    "  target:   ", data$target[1], "\n",
    "  run_id: ", data$run_id[1], "\n")


    model <- BayesPharma::sigmoid_model(
        data = model_data,
        formula = BayesPharma::sigmoid_antagonist_formula(
            response_units = "% Displacement of Reference Ligand",
            predictor = 1),
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
        backend = "cmdstanr",
        expose_functions = FALSE)

     result <- tibble::tibble(
       data = list(data),
       model = list(model))

     save(result, file = output_path)
}



```{r fit-per-run-models}

n_jobs <- MOR_binding |>
  dplyr::distinct(protocol_id, researcher, date_code, assay_id) |>
  nrow()

n_models <- MOR_binding |>
  dplyr::distinct(protocol_id, researcher, date_code, assay_id, target, run_id) |>
  nrow()

cat("Submitting ", n_jobs, " jobs to fit ", n_models, " models\n", sep = "")

job_index <- 0
models <- MOR_binding |>
  dplyr::group_by(protocol_id, researcher, date_code, assay_id) |> dplyr::do({
    job_index <<- job_index + 1
    assay_data <- .

    n_models_per_assay <- assay_data |> dplyr::distinct(target, run_id) |> nrow()

    cat("Submitting job ", job_index, " to fit ", n_models_per_assay, " models\n", sep = "")
    job <- future::future({
      model_index <- 0
      assay_data |> dplyr::group_by(target, run_id) |> dplyr::do({
        run_data <- .
        model_index <<- model_index + 1
        fit_model(
          data = run_data,
          output_path = paste0(
            intermediate_path, "/model_per_run_", date_code, "_",
            "j", job_index, "_m", model_index, ".Rdata"))
        data.frame()
      })
    })
    data.frame()
  })
```


```{r gather-fit-models}

models <- data.frame(
    path = list.files(
        path = intermediate_path,
        pattern = "model_per_run_.+.Rdata",
        full.names = TRUE)) |>
        dplyr::rowwise() |>
        dplyr::do({
            path <- .$path
            cat("loading ", path, "\n", sep = "")
            load(path)
            result
        })
```

```{r summarize-model-fits}


model_summaries <- models |>
    dplyr::rowwise() |>
    dplyr::do({
        row <- .
        data <- row$data
        model <- row$model

        model |>
             posterior::summarize_draws() |>
             dplyr::mutate(
                 protocol_id = data$protocol_id[1],
                 researcher = data$researcher[1],
                 date_code = data$date_code[1],
                 assay_id = data$assay_id[1],
                 target = data$target[1],
                 run_id = data$run_id[1],
                 .before = 1)
    }) |>
    dplyr::ungroup()


model_summaries <- model_summaries |>
  dplyr::left_join(
      MOR_binding |>
        dplyr::distinct(
          protocol_id, researcher, date_code, target, run_id, substance_id)
      by = c("protocol_id", "researcher", "date_code", "target", "run_id"))

model_summaries |>
    readr::write_tsv(
        file = paste0(output_path, "/model_by_run_fit_summaries_", date_code, ".tsv"))

```


```{r plot-estimates}

plot_data <- model_summaries |>
    dplyr::filter(variable == "b_ic50_Intercept")

plot_data <- plot_data |>
    dplyr::semi_join(
        plot_data |>
        dplyr::count(substance_id) |>
        dplyr::filter(n >= 20),
        by = "substance_id")

plot <- ggplot2::ggplot(data = plot_data) +
    ggplot2::ggtitle(
        "Estiamted Binding Affinities") +
    ggplot2::theme_bw() +
    ggplot2::geom_abline(
        intercept = 0,
        slope = 1,
        linewidth = 1.2) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = mean,
            y = q95 - q5,
            color = researcher,
            shape = target),
        alpha = .8) +
    ggplot2::facet_wrap(
        facets = dplyr::vars(substance_id)) +
    ggplot2::scale_x_continuous(
        "Estiamted IC50 of binding",
        breaks = c(-9, -7, -5),
        labels =c("1nM", "100nM", "10uM")) +
    ggplot2::scale_y_continuous(
        "Log10 of the [5%, 95%] Credible Interval Width")

ggplot2::ggsave(
    filename = paste0(output_path, "/ic50_confidence_by_substance_id_", date_code, ".pdf"),
    plot = plot,
    width = 10,
    height = 8,
    useDingbats = FALSE)

plot_data |>
    readr::write_tsv(
        file = paste0(output_path, "/ic50_confidence_by_substance_id_", date_code, ".tsv"))


```
