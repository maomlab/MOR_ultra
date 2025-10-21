

---
title: "02.2.1 fit one model per run on cluster"
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
