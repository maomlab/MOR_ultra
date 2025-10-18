library(tidymodels)

load("data/rmsd_features.Rdata")

# Filter your data first
data <- features |>
  dplyr::filter(params == "Wtemplate")

# Define your model formula
model_formula <- as.formula(
  `MOR Emax` ~ `8E0G` + `9MQH` + `9BJK` + `5C1M` + `4DKL` + `9MQJ` + `7UL4` +
    `8QOT` + `9MQI` + `8EF5` + `8EF6` + `8EFB` + `7T2H` + `8EFL` + `8EFO` +
    `7SBF` + `8EFQ` + `8F7R` + `7T2G` + `8K9L` + `6DDE` + `7SCG` + `8F7Q` +
    `6DDF` + `8K9K` + `9BQJ` + `7U2K` + `7U2L`
)

# Set up leave-ligand-out CV
cv_splits <- data |>
   rsample::group_vfold_cv(
     group = ligand,
     v = 5) #dplyr::n_distinct(data$ligand))

# Define model spec
lm_spec <- parsnip::linear_reg() |>
  parsnip::set_engine("lm")

# Fit model with resampling
fit_resamples <- tune::fit_resamples(
  lm_spec,
  model_formula,
  resamples = cv_splits,
  control = tune::control_resamples(save_pred = TRUE))

# Collect predictions on the held-out ligands
predictions <- collect_predictions(fit_resamples) |>
  dplyr::arrange(.row)

# View predictions
print(predictions)


###############

plot <- ggplot2::ggplot(
  data = predictions) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = .pred,
      y = `MOR Emax`),
    size = 1.2,
    alpha = 0.9,
    shape = 16) +
  ggplot2::geom_smooth(
    mapping = ggplot2::aes(
      x = .pred,
      y = `MOR Emax`),
    method = "lm",
    formula = y ~ x,
    se = FALSE) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_x_continuous("Predicted Emax")

ggplot2::ggsave(
  filename = "product/02_rmsd_Emax_cv_pred_20250603.pdf",
  width = 4,
  height = 4)
