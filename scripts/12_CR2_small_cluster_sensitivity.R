# ============================================================
# CR2 SMALL-CLUSTER SENSITIVITY — FIXED VERSION
#
# Endpoints:
#   1. MAC M2
#   2. FIB M1
#   3. FIB M3
#   4. MAC-FIB spatial neighborhood
#
# Method:
#   patient-balanced weighted lm
#   CR2 cluster-robust covariance
#   HTZ small-sample joint Wald test
#
# Cluster = Patient
# ============================================================


# ============================================================
# 0. Install / load package
# ============================================================

if (!requireNamespace("clubSandwich", quietly = TRUE)) {

  install.packages(
    "clubSandwich",
    repos = "https://cloud.r-project.org"
  )
}


library(clubSandwich)
library(splines)


# ============================================================
# 1. Paths
# ============================================================

# Repository-friendly working directory.
# Set the environment variable CGN_PROJECT_DIR to override.
BASE <- Sys.getenv(
  "CGN_PROJECT_DIR",
  unset = file.path(getwd(), "workspace")
)

MODULE_SCORE_PATH <- file.path(
  BASE,
  "figure5_frozen_module_scores.csv"
)

SPATIAL_PATH <- file.path(
  BASE,
  "figure4_driver_neighbor_data_smoothed.csv"
)

OUT_DIR <- file.path(
  BASE,
  "reproduction",
  "results"
)


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


if (!file.exists(MODULE_SCORE_PATH)) {

  stop(
    paste(
      "Missing:",
      MODULE_SCORE_PATH
    )
  )
}


if (!file.exists(SPATIAL_PATH)) {

  stop(
    paste(
      "Missing:",
      SPATIAL_PATH
    )
  )
}


cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\nCR2 SMALL-CLUSTER SENSITIVITY\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")


# ============================================================
# 2. Load data
# ============================================================

module_scores <- read.csv(
  MODULE_SCORE_PATH,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


spatial <- read.csv(
  SPATIAL_PATH,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = 1
)


cat(
  "\nModule-score table:",
  nrow(module_scores),
  "x",
  ncol(module_scores),
  "\n"
)


cat(
  "Spatial table:",
  nrow(spatial),
  "x",
  ncol(spatial),
  "\n"
)


# ============================================================
# 3. Required-column check
# ============================================================

required_module_cols <- c(
  "module_label",
  "module_score",
  "Disease",
  "Patient",
  "PC1"
)


required_spatial_cols <- c(
  "MAC_to_FIB",
  "Disease",
  "Patient",
  "PC1"
)


missing_module <- setdiff(
  required_module_cols,
  colnames(module_scores)
)


missing_spatial <- setdiff(
  required_spatial_cols,
  colnames(spatial)
)


if (length(missing_module) > 0) {

  stop(
    paste(
      "Module table missing:",
      paste(
        missing_module,
        collapse = ", "
      )
    )
  )
}


if (length(missing_spatial) > 0) {

  stop(
    paste(
      "Spatial table missing:",
      paste(
        missing_spatial,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 4. Prepare endpoint
# ============================================================

prepare_endpoint <- function(
  dat,
  outcome
) {

  cols <- c(
    outcome,
    "Disease",
    "Patient",
    "PC1"
  )


  d <- dat[
    complete.cases(
      dat[, cols, drop = FALSE]
    ),
    cols,
    drop = FALSE
  ]


  # LN vs anti-GBM
  d <- d[
    d$Disease %in% c(
      "SLE",
      "GBM"
    ),
    ,
    drop = FALSE
  ]


  d$Disease <- factor(
    d$Disease,
    levels = c(
      "SLE",
      "GBM"
    )
  )


  d$Patient <- as.character(
    d$Patient
  )


  d$PC1 <- as.numeric(
    d$PC1
  )


  d[[outcome]] <- as.numeric(
    d[[outcome]]
  )


  # ----------------------------------------------------------
  # Patient-balanced weights
  # ----------------------------------------------------------

  patient_counts <- table(
    d$Patient
  )


  d$patient_weight <- 1 /
    as.numeric(
      patient_counts[
        d$Patient
      ]
    )


  return(d)
}


# ============================================================
# 5. CR2 + HTZ function
# ============================================================

fit_CR2 <- function(
  d,
  outcome,
  endpoint_name
) {

  formula_text <- paste0(
    outcome,
    " ~ ",
    "bs(PC1, df = 3, degree = 3, intercept = FALSE)",
    " * Disease"
  )


  f <- as.formula(
    formula_text
  )


  # ----------------------------------------------------------
  # Weighted lm
  # ----------------------------------------------------------

  fit <- lm(
    f,
    data = d,
    weights = patient_weight
  )


  coef_names <- names(
    coef(fit)
  )


  cat("\n")
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")
  cat(endpoint_name)
  cat("\n")
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")


  cat(
    "N ROI:",
    nrow(d),
    "\n"
  )


  cat(
    "Patients:",
    length(
      unique(
        d$Patient
      )
    ),
    "\n"
  )


  cat(
    "LN patients:",
    length(
      unique(
        d$Patient[
          d$Disease == "SLE"
        ]
      )
    ),
    "\n"
  )


  cat(
    "anti-GBM patients:",
    length(
      unique(
        d$Patient[
          d$Disease == "GBM"
        ]
      )
    ),
    "\n"
  )


  # ----------------------------------------------------------
  # Identify 3 Disease × spline interaction coefficients
  # ----------------------------------------------------------

  interaction_idx <- which(
    grepl(
      "bs\\(PC1",
      coef_names
    )
    &
    grepl(
      "DiseaseGBM",
      coef_names
    )
  )


  interaction_names <- coef_names[
    interaction_idx
  ]


  cat(
    "\nInteraction coefficients:\n"
  )


  print(
    interaction_names
  )


  if (length(interaction_idx) != 3) {

    stop(
      paste0(
        endpoint_name,
        ": expected 3 interaction coefficients, found ",
        length(interaction_idx)
      )
    )
  }


  # ----------------------------------------------------------
  # CR2 variance
  # ----------------------------------------------------------

  V_CR2 <- vcovCR(
    fit,
    cluster = d$Patient,
    type = "CR2",
    inverse_var = FALSE
  )


  # ----------------------------------------------------------
  # Joint H0:
  # all 3 interaction coefficients = 0
  # ----------------------------------------------------------

  constraints <- constrain_zero(
    interaction_idx,
    coefs = coef(fit)
  )


  WT <- Wald_test(
    fit,
    constraints = constraints,
    vcov = V_CR2,
    test = c(
      "HTZ",
      "chi-sq"
    ),
    tidy = TRUE
  )


  cat(
    "\nCR2 Wald tests:\n"
  )


  print(
    WT
  )


  HTZ_row <- WT[
    WT$test == "HTZ",
    ,
    drop = FALSE
  ]


  CHI_row <- WT[
    WT$test == "chi-sq",
    ,
    drop = FALSE
  ]


  if (nrow(HTZ_row) != 1) {

    stop(
      paste(
        endpoint_name,
        ": HTZ result not found."
      )
    )
  }


  result <- data.frame(

    endpoint =
      endpoint_name,

    n_ROI =
      nrow(d),

    n_patients =
      length(
        unique(
          d$Patient
        )
      ),

    n_LN_patients =
      length(
        unique(
          d$Patient[
            d$Disease == "SLE"
          ]
        )
      ),

    n_GBM_patients =
      length(
        unique(
          d$Patient[
            d$Disease == "GBM"
          ]
        )
      ),

    HTZ_F =
      as.numeric(
        HTZ_row$Fstat
      ),

    HTZ_df_num =
      as.numeric(
        HTZ_row$df_num
      ),

    HTZ_df_denom =
      as.numeric(
        HTZ_row$df_denom
      ),

    CR2_HTZ_P =
      as.numeric(
        HTZ_row$p_val
      ),

    CR2_chisq_P =
      if (nrow(CHI_row) == 1) {
        as.numeric(
          CHI_row$p_val
        )
      } else {
        NA_real_
      },

    stringsAsFactors = FALSE
  )


  return(result)
}


# ============================================================
# 6. Build endpoint datasets
# ============================================================

selected_modules <- c(
  "MAC M2",
  "FIB M1",
  "FIB M3"
)


endpoint_data <- list()


for (module_name in selected_modules) {

  temp <- module_scores[
    module_scores$module_label == module_name,
    ,
    drop = FALSE
  ]


  endpoint_data[[module_name]] <- prepare_endpoint(
    temp,
    "module_score"
  )
}


# ------------------------------------------------------------
# MAC-FIB spatial
# ------------------------------------------------------------

spatial_temp <- spatial[
  spatial$Disease %in% c(
    "SLE",
    "GBM"
  ),
  ,
  drop = FALSE
]


endpoint_data[["MAC-FIB"]] <- prepare_endpoint(
  spatial_temp,
  "MAC_to_FIB"
)


# ============================================================
# 7. Endpoint audit
# ============================================================

audit_rows <- list()


for (endpoint_name in names(endpoint_data)) {

  d <- endpoint_data[[endpoint_name]]


  audit_rows[[endpoint_name]] <- data.frame(

    endpoint =
      endpoint_name,

    n_ROI =
      nrow(d),

    n_patients =
      length(
        unique(
          d$Patient
        )
      ),

    n_LN_patients =
      length(
        unique(
          d$Patient[
            d$Disease == "SLE"
          ]
        )
      ),

    n_GBM_patients =
      length(
        unique(
          d$Patient[
            d$Disease == "GBM"
          ]
        )
      ),

    stringsAsFactors = FALSE
  )
}


endpoint_audit <- do.call(
  rbind,
  audit_rows
)


rownames(
  endpoint_audit
) <- NULL


cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\nENDPOINT AUDIT\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")


print(
  endpoint_audit,
  row.names = FALSE
)


# ============================================================
# 8. Run four CR2 models
# ============================================================

results <- list()


results[["MAC M2"]] <- fit_CR2(
  endpoint_data[["MAC M2"]],
  outcome = "module_score",
  endpoint_name = "MAC M2"
)


results[["FIB M1"]] <- fit_CR2(
  endpoint_data[["FIB M1"]],
  outcome = "module_score",
  endpoint_name = "FIB M1"
)


results[["FIB M3"]] <- fit_CR2(
  endpoint_data[["FIB M3"]],
  outcome = "module_score",
  endpoint_name = "FIB M3"
)


results[["MAC-FIB"]] <- fit_CR2(
  endpoint_data[["MAC-FIB"]],
  outcome = "MAC_to_FIB",
  endpoint_name = "MAC-FIB"
)


CR2_results <- do.call(
  rbind,
  results
)


rownames(
  CR2_results
) <- NULL


# ============================================================
# 9. Add original clustered P values
# ============================================================

original_P <- c(

  "MAC M2" =
    7.462782e-05,

  "FIB M1" =
    1.267494e-06,

  "FIB M3" =
    3.961392e-08,

  "MAC-FIB" =
    1.794172e-02
)


CR2_results$original_clustered_P <-
  original_P[
    CR2_results$endpoint
  ]


# ============================================================
# 10. BH correction across four sensitivity endpoints
# ============================================================

CR2_results$CR2_HTZ_FDR_across_4 <-
  p.adjust(
    CR2_results$CR2_HTZ_P,
    method = "BH"
  )


CR2_results$CR2_HTZ_nominal_lt_005 <-
  CR2_results$CR2_HTZ_P < 0.05


CR2_results$CR2_HTZ_FDR_lt_005 <-
  CR2_results$CR2_HTZ_FDR_across_4 < 0.05


# ============================================================
# 11. Interpretation
# ============================================================

CR2_results$interpretation <- ifelse(

  CR2_results$CR2_HTZ_FDR_across_4 < 0.05,

  "Robust after CR2-HTZ and BH correction",

  ifelse(

    CR2_results$CR2_HTZ_P < 0.05,

    paste0(
      "Nominally significant with CR2-HTZ; ",
      "not FDR<0.05 across four sensitivity endpoints"
    ),

    "Attenuated under CR2-HTZ small-sample inference"
  )
)


# ============================================================
# 12. Full results
# ============================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\nFINAL CR2 SMALL-CLUSTER RESULTS\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")


print(
  CR2_results,
  row.names = FALSE,
  digits = 6
)


# ============================================================
# 13. Compact reviewer table
# ============================================================

compact <- CR2_results[
  ,
  c(
    "endpoint",
    "n_patients",
    "n_GBM_patients",
    "original_clustered_P",
    "CR2_chisq_P",
    "HTZ_df_denom",
    "CR2_HTZ_P",
    "CR2_HTZ_FDR_across_4",
    "interpretation"
  )
]


cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\nREVIEWER-ORIENTED SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")


print(
  compact,
  row.names = FALSE,
  digits = 6
)


# ============================================================
# 14. Save
# ============================================================

RESULT_PATH <- file.path(
  OUT_DIR,
  "CR2_small_cluster_results.csv"
)


AUDIT_PATH <- file.path(
  OUT_DIR,
  "CR2_small_cluster_endpoint_audit.csv"
)


write.csv(
  CR2_results,
  RESULT_PATH,
  row.names = FALSE
)


write.csv(
  endpoint_audit,
  AUDIT_PATH,
  row.names = FALSE
)


# ============================================================
# 15. Final counts
# ============================================================

n_nominal <- sum(
  CR2_results$CR2_HTZ_nominal_lt_005
)


n_fdr <- sum(
  CR2_results$CR2_HTZ_FDR_lt_005
)


cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\nCR2 ANALYSIS COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")


cat(
  "\nCR2-HTZ nominal P < 0.05:",
  n_nominal,
  "/ 4\n"
)


cat(
  "CR2-HTZ FDR < 0.05 across four:",
  n_fdr,
  "/ 4\n"
)


cat(
  "\nSaved:\n"
)


cat(
  RESULT_PATH,
  "\n"
)


cat(
  AUDIT_PATH,
  "\n"
)


cat(
  paste(rep("=", 80), collapse = ""),
  "\n"
)
