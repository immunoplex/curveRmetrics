library(DBI)
library(dplyr)

get_db_connection <- function() {
  dbConnect(RPostgres::Postgres(),
            dbname = Sys.getenv("db"),
            host = Sys.getenv("db_host"),
            port = Sys.getenv("db_port"),
            user = Sys.getenv("db_userid_x"),
            password = Sys.getenv("db_pwd_x"),
            sslmode = 'require',
            options = "-c search_path=madi_results"
  )
}

conn <- get_db_connection()


get_freq_curves <- function(con = get_db_connection()) {
  sql <- "
    SELECT
      freq_curves_id,
      curve_id,
      'frequentist' AS method,
      iter,
      status,
      model_name,
      a, b, c, d, g,
      dfresidual,
      nobs,
      rsquare_fit,
      aic,
      bic,
      loglik,
      mse,
      cv,
      bkg_method,
      is_log_response,
      is_log_x,
      formula
    FROM madi_results.freq_curves;
  "

  # Pull the whole result set into a tibble (data.frame)
  DBI::dbGetQuery(con, sql) %>% as_tibble()
}

# ------------------------------------------------------------------
# 2️⃣  bayes_curves table (only the columns you asked for)
# ------------------------------------------------------------------
get_bayes_curves <- function(con = get_db_connection()) {
  sql <- "
    SELECT
      curve_id,
      'bayesian'    AS method,
      curve_family AS model_name,
      a,
      b,
      c,
      d,
      g,
      apply_prozone
    FROM madi_results.bayes_curves;
  "

  DBI::dbGetQuery(con, sql) %>% as_tibble()
}


get_bayes_param_ci  <- function(con = get_db_connection()) {
  sql <- "WITH base AS (
  SELECT *
  FROM madi_results.bayes_ensemble
  WHERE is_global_best = true
)
SELECT
  curve_id,
  'bayesian' AS method,
  family AS model_name,
  TRUE AS converged,
  v.parameter,
  v.estimate,
  v.conf_lower,
  v.conf_upper,
  TRUE AS is_best_model
FROM base
CROSS JOIN LATERAL (
  VALUES
    ('a', a, a_lower, a_upper),
    ('b', b, b_lower, b_upper),
    ('c', c, c_lower, c_upper),
    ('d', d, d_lower, d_upper),
    ('g', g, g_lower, g_upper)
) AS v(parameter, estimate, conf_lower, conf_upper);"
DBI::dbGetQuery(con, sql) %>% as_tibble()
}

get_freq_parameters <- function(con = get_db_connection()) {
  sql <- "
SELECT  curve_id, 'frequentist' AS method, model_name, converged, parameter, estimate, conf_lower, conf_upper, is_best_model
	FROM madi_results.freq_candidate_parameters
	WHERE is_best_model = true;
  "
  DBI::dbGetQuery(con, sql) %>% as_tibble()

}

get_freq_second_derivative <- function(con = get_db_connection()) {
  sql <- "SELECT freq_second_derivative_id, curve_id, model_name, concentration, d2x_y
	FROM madi_results.freq_second_derivative;"

  DBI::dbGetQuery(con, sql) %>% as_tibble()
}
