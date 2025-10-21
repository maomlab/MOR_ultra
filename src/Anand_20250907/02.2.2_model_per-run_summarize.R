
---
title: "02.2 fit one model per run on cluster"
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
