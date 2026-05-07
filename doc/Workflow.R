## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
devtools::load_all()
script_path <- system.file("vignette_helpers", "db_functions.R", package = "curveRmetrics")
if (file.exists(script_path)) {
  source(script_path)
} else {
  warning("Script file not found at expected path: ", script_path)
}
conn <- get_db_connection()

freq_curves_df  <- get_freq_curves(conn)
bayes_curves_df <- get_bayes_curves(conn)
bayes_param_ci_df <- get_bayes_param_ci(conn)
freq_param_ci_df <- get_freq_parameters(conn)

## -----------------------------------------------------------------------------
freq_curves_df  <- compute_inflection_point(freq_curves_df)
bayes_curves_df <- compute_inflection_point(bayes_curves_df)

head(freq_curves_df)

## -----------------------------------------------------------------------------
freq_sens <-compute_assay_sensitivity(freq_curves_df)
bayes_sens <- compute_assay_sensitivity(bayes_curves_df)

head(freq_sens)

## -----------------------------------------------------------------------------
lods <-generate_lods(freq_param_ci_df)
bayes_lod <- generate_lods(bayes_param_ci_df)
head(bayes_lod)

## -----------------------------------------------------------------------------
freq_rdl_tbl <- compute_mdc_rdl(freq_param_ci_df, lods = lods)

bayes_rdl_tbl <- compute_mdc_rdl(bayes_param_ci_df, lods = bayes_lod)

head(freq_rdl_tbl)

## -----------------------------------------------------------------------------
d2_freq <- compute_second_deriv_df(curves_df = freq_curves_df)

d2_bayes <- compute_second_deriv_df(curves_df = bayes_curves_df)

## -----------------------------------------------------------------------------
freq_loqs <- compute_loqs(freq_curves_df, second_deriv_df = d2_freq)

bayes_loqs <- compute_loqs(bayes_curves_df, second_deriv_df = d2_bayes)

head(bayes_loqs)

## -----------------------------------------------------------------------------
qc_freq <-attach_quality_metrics(freq_curves_df, lods, rdls = freq_rdl_tbl, 
                                 sensitivity = freq_sens,
                                 loqs = freq_loqs)

qc_bayes <-attach_quality_metrics(bayes_curves_df, lods = bayes_lod, 
                                  rdls = bayes_rdl_tbl, 
                                  sensitivity = bayes_sens, 
                                  loqs = bayes_loqs)

head(qc_freq)

