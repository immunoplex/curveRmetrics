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


get_freq_curves <- function(con = get_db_connection(), study_accession) {
  sql <- glue::glue_sql("
    SELECT
      freq_curves.curve_id,
      'frequentist' AS method,
      model_name,
      a, b, c, d, g
    FROM madi_results.freq_curves
    INNER JOIN madi_results.curve_lookup
    ON curve_lookup.curve_id = freq_curves.curve_id
    WHERE curve_lookup.study_accession = {study_accession};
  ",  .con = con)
  
  # Pull the whole result set into a tibble (data.frame)
  DBI::dbGetQuery(con, sql) %>% as_tibble()
}

# ------------------------------------------------------------------
# 2.  bayes_curves table
# ------------------------------------------------------------------
get_bayes_curves <- function(con = get_db_connection(), study_accession) {
  sql <- glue::glue_sql("
    SELECT
      bayes_curves.curve_id,
      'bayesian'    AS method,
      curve_family AS model_name,
      a,
      b,
      c,
      d,
      g
    FROM madi_results.bayes_curves
    INNER JOIN madi_results.curve_lookup
    ON curve_lookup.curve_id = bayes_curves.curve_id
    WHERE curve_lookup.study_accession = {study_accession};
  ",  .con = con)
  
  DBI::dbGetQuery(con, sql) %>% as_tibble()
}


get_freq_parameters <- function(con = get_db_connection()) {
  sql <- "
SELECT  curve_id, 'frequentist' AS method, model_name, parameter, estimate, conf_lower, conf_upper
	FROM madi_results.freq_candidate_parameters
	WHERE is_best_model = true;
  "
  DBI::dbGetQuery(con, sql) %>% as_tibble()
  
}

get_bayes_param_ci  <- function(con = get_db_connection(), study_accession) {
  sql <- glue::glue_sql("SELECT curve_id,
    family AS model_name,
    'bayesian' as method,
    parameter,
    estimate,
    conf_lower,
    conf_upper

	FROM  (
SELECT DISTINCT ON (plateid, parameter)
    plateid,
    plate_elpd,
    curve_id,
    family,
    parameter,
    estimate,
    conf_lower,
    conf_upper
FROM madi_results.bayes_ensemble
CROSS JOIN LATERAL (
    VALUES
        ('a', a, a_lower, a_upper),
        ('b', b, b_lower, b_upper),
        ('c', c, c_lower, c_upper),
        ('d', d, d_lower, d_upper),
        ('g', g, g_lower, g_upper)
) AS params(parameter, estimate, conf_lower, conf_upper)
WHERE study_accession = {study_accession}
ORDER BY plateid, parameter, plate_elpd DESC

	) AS a
;", ,  .con = con)
  #   sql <- glue::glue_sql("WITH base AS (
  #   SELECT *
  #   FROM madi_results.bayes_ensemble
  #   WHERE is_global_best = true
  # )
  # SELECT
  #   base.curve_id,
  #   'bayesian' AS method,
  #   family AS model_name,
  #   TRUE AS converged,
  #   v.parameter,
  #   v.estimate,
  #   v.conf_lower,
  #   v.conf_upper,
  #   TRUE AS is_best_model
  # FROM base
  # INNER JOIN madi_results.curve_lookup
  #     ON curve_lookup.curve_id = base.curve_id
  # CROSS JOIN LATERAL (
  #   VALUES
  #     ('a', a, a_lower, a_upper),
  #     ('b', b, b_lower, b_upper),
  #     ('c', c, c_lower, c_upper),
  #     ('d', d, d_lower, d_upper),
  #     ('g', g, g_lower, g_upper)
  # ) AS v(parameter, estimate, conf_lower, conf_upper)
  # WHERE curve_lookup.study_accession = {study_accession};",  .con = con)
  # 
  
  DBI::dbGetQuery(con, sql) %>% as_tibble()
}


get_freq_parameters <- function(con = get_db_connection(), study_accession) {
  sql <- glue::glue_sql("
SELECT  base.curve_id, 'frequentist' AS method, model_name, parameter, estimate, conf_lower, conf_upper
	FROM madi_results.freq_candidate_parameters as base
	INNER JOIN madi_results.curve_lookup
    ON curve_lookup.curve_id = base.curve_id
	WHERE is_best_model = true
	AND curve_lookup.study_accession = {study_accession};
  ",  .con = con)
  DBI::dbGetQuery(con, sql) %>% as_tibble()
  
}


get_standards <- function(con = get_db_connection(), study_accession) {
  #   sql <- glue::glue_sql("
  #   SELECT curve_lookup.curve_id, well, sampleid,
  # (1.0 / NULLIF(xmap_standard.dilution, 0)) * af.standard_curve_concentration AS concentration,
  # xmap_standard.antigen, antibody_mfi AS assay_response, xmap_standard.feature
  # 	FROM madi_results.xmap_standard
  # 	INNER JOIN madi_results.curve_lookup
  # 	ON curve_lookup.project_id = xmap_standard.project_id
  # 	AND curve_lookup.study_accession = xmap_standard.study_accession
  # 	AND curve_lookup.experiment_accession = xmap_standard.experiment_accession
  # 	AND curve_lookup.plate = xmap_standard.plate
  # 	AND curve_lookup.nominal_sample_dilution = xmap_standard.nominal_sample_dilution
  # 	AND curve_lookup.source = xmap_standard.source
  # 	AND curve_lookup.wavelength = xmap_standard.wavelength
  # 	AND curve_lookup.antigen = xmap_standard.antigen
  # 	AND curve_lookup.feature = xmap_standard.feature
  #     LEFT JOIN madi_results.xmap_antigen_family AS af
  #     ON xmap_standard.antigen               = af.antigen
  #     AND xmap_standard.experiment_accession = af.experiment_accession
  #     AND xmap_standard.project_id           = af.project_id
  # 	WHERE curve_lookup.study_accession = {study_accession}
  # 	;
  #   ",  .con = con)
  sql <- glue::glue_sql("
  SELECT  curve_lookup.curve_id, well, sampleid,
(1.0 / NULLIF(xmap_standard.dilution, 0)) * af.standard_curve_concentration AS concentration,
xmap_standard.antigen, antibody_mfi AS assay_response, xmap_standard.feature
	FROM madi_results.xmap_standard
	INNER JOIN madi_results.curve_lookup
	ON curve_lookup.project_id = xmap_standard.project_id
	AND curve_lookup.study_accession = xmap_standard.study_accession
	AND curve_lookup.experiment_accession = xmap_standard.experiment_accession
	AND curve_lookup.plate = xmap_standard.plate
	AND curve_lookup.nominal_sample_dilution = xmap_standard.nominal_sample_dilution
	AND curve_lookup.source = xmap_standard.source
	AND curve_lookup.wavelength = xmap_standard.wavelength
	AND curve_lookup.antigen = xmap_standard.antigen
	AND curve_lookup.feature = xmap_standard.feature
   LEFT JOIN madi_results.xmap_antigen_family AS af
    ON xmap_standard.antigen               = af.antigen
	AND xmap_standard.feature              = af.feature
	AND xmap_standard.study_accession      = af.study_accession
    AND xmap_standard.experiment_accession = af.experiment_accession
    AND xmap_standard.project_id           = af.project_id
	WHERE curve_lookup.study_accession = {study_accession}
	;
  ",  .con = con)
  
  
  
  DBI::dbGetQuery(con, sql) %>% as_tibble()
}

