library(data.table)

## ============================================================================
## 0. Setup + variable groups
## ============================================================================
cohort    <- "Cohort_name"   # one of: "GenR", "MOBA", "PREDO", "Raine"
data_path <- "path/to/dataset/final_dataset_cohort_name.csv"   # <-- point to the PREDO file when cohort = "PREDO"

results_path <- "path/to/results/results"

## one subfolder per cohort, e.g. .../results/descriptives_GenR
out_dir <- file.path(results_path, paste0("descriptives_", cohort))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(name) file.path(out_dir, paste0(name, "_", cohort, ".", 
                                                     ifelse(name == "distributions", "pdf", "csv")))

dataset <- as.data.frame(fread(data_path))

prs_vars <- c("prs_PPD_c","prs_PPD_m","prs_PPD_f",
              "prs_PPD_unique_c","prs_PPD_unique_m","prs_PPD_unique_f",
              "prs_MDD_c","prs_MDD_m","prs_MDD_f",
              "prs_MDD_female_c","prs_MDD_female_m","prs_MDD_female_f")
pc_vars  <- paste0("PC", 1:10)

## Continuous vs categorical (full set, as used for GenR and MoBa)
continuous_vars <- c("age_pregnancy_m", "dep_pre_m","dep_6m_m","dep_36m_m",
                     "BMI_m", "birthweight_c", "age_gestational_c",
                     "intern_3y_c","intern_9y_c","intern_14y_c",
                     "age_3y_c","age_9y_c","age_14y_c", prs_vars, pc_vars)
categorical_vars <- c("gender_c", "smoke_m", "history_dep_m")

outcome_vars  <- c("intern_3y_c","intern_9y_c","intern_14y_c")
mediator_vars <- c("history_dep_m","dep_pre_m","dep_6m_m","dep_36m_m",
                   "BMI_m","birthweight_c","smoke_m")
sex_var       <- "gender_c"          # child sex; confirm coding with table(dataset[[sex_var]])

## ----------------------------------------------------------------------------
## Cohort-specific conditions
##   GenR, MOBA : run everything as is
##   PREDO      : no history_dep_m, no internalizing at 14y (so no age_14y_c)
##   Raine      : no dep_pre_m and no dep_36m_m
## ----------------------------------------------------------------------------
cohort_key <- toupper(cohort)

if (cohort_key %in% c("GENR", "MOBA")) {
  drop_vars <- character(0)
} else if (cohort_key == "PREDO") {
  drop_vars <- c("history_dep_m", "intern_14y_c", "age_14y_c")
} else if (cohort_key == "RAINE") {
  drop_vars <- c("dep_pre_m", "dep_36m_m")
} else {
  stop("Unknown cohort '", cohort, "'. Use GenR, MOBA, PREDO or Raine.")
}

continuous_vars  <- setdiff(continuous_vars,  drop_vars)
categorical_vars <- setdiff(categorical_vars, drop_vars)
mediator_vars    <- setdiff(mediator_vars,    drop_vars)
outcome_vars     <- setdiff(outcome_vars,     drop_vars)

if (length(drop_vars) > 0)
  message("Cohort ", cohort, ": excluding ", paste(drop_vars, collapse = ", "))

## keep only variables that exist in this cohort (safety net)
continuous_vars  <- intersect(continuous_vars,  names(dataset))
categorical_vars <- intersect(categorical_vars, names(dataset))
mediator_vars    <- intersect(mediator_vars,    names(dataset))
outcome_vars     <- intersect(outcome_vars,     names(dataset))

## ============================================================================
## 1. Helper summaries
## ============================================================================
.skew <- function(x){ x <- x[!is.na(x)]; n <- length(x)
if (n < 3) return(NA_real_); m <- mean(x); s <- sd(x)
if (s == 0) return(NA_real_); (sum((x - m)^3)/n) / s^3 }
.kurt <- function(x){ x <- x[!is.na(x)]; n <- length(x)
if (n < 4) return(NA_real_); m <- mean(x); s <- sd(x)
if (s == 0) return(NA_real_); (sum((x - m)^4)/n) / s^4 - 3 }

describe_continuous <- function(dt, vars, cohort){
  do.call(rbind, lapply(vars, function(v){
    x <- suppressWarnings(as.numeric(dt[[v]]))
    data.frame(cohort = cohort, variable = v, type = "continuous",
               n = sum(!is.na(x)), pct_missing = round(mean(is.na(x)) * 100, 1),
               mean = round(mean(x, na.rm = TRUE), 3), sd = round(sd(x, na.rm = TRUE), 3),
               median = round(median(x, na.rm = TRUE), 3),
               q1 = round(quantile(x, .25, na.rm = TRUE), 3),
               q3 = round(quantile(x, .75, na.rm = TRUE), 3),
               min = round(min(x, na.rm = TRUE), 3), max = round(max(x, na.rm = TRUE), 3),
               skewness = round(.skew(x), 2), kurtosis = round(.kurt(x), 2),
               nonnormal_flag = abs(.skew(x)) > 2 | abs(.kurt(x)) > 7,
               row.names = NULL)
  }))
}

describe_categorical <- function(dt, vars, cohort){
  do.call(rbind, lapply(vars, function(v){
    x <- dt[[v]]; tab <- table(x, useNA = "no"); tot <- sum(tab)
    data.frame(cohort = cohort, variable = v, level = names(tab),
               n = as.integer(tab), pct = round(as.integer(tab)/tot * 100, 1),
               n_missing = sum(is.na(x)), pct_missing = round(mean(is.na(x)) * 100, 1),
               row.names = NULL)
  }))
}

## ============================================================================
## 2. Table 1 — overall + by child sex (with a group-comparison p)
## ============================================================================
fmt_meansd  <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
fmt_np      <- function(x, lvl) { n <- sum(x == lvl, na.rm=TRUE)
sprintf("%d (%.1f%%)", n, 100*n/sum(!is.na(x))) }

make_table1 <- function(dt, cont_vars, cat_vars, group_var, cohort){
  g   <- dt[[group_var]]; grps <- sort(unique(g[!is.na(g)])); rows <- list()
  
  for (v in cont_vars){
    x <- suppressWarnings(as.numeric(dt[[v]]))
    p <- tryCatch(t.test(x ~ g)$p.value, error = function(e) NA_real_)
    r <- data.frame(variable = v, level = "mean (SD)", overall = fmt_meansd(x))
    for (gr in grps) r[[paste0("g_", gr)]] <- fmt_meansd(x[g == gr])
    r$p <- round(p, 4); rows[[length(rows)+1]] <- r
  }
  for (v in cat_vars){
    x <- dt[[v]]
    p <- tryCatch(suppressWarnings(chisq.test(table(x, g))$p.value), error = function(e) NA_real_)
    for (lvl in sort(unique(x[!is.na(x)]))){
      r <- data.frame(variable = v, level = as.character(lvl), overall = fmt_np(x, lvl))
      for (gr in grps) r[[paste0("g_", gr)]] <- fmt_np(x[g == gr], lvl)
      r$p <- round(p, 4); rows[[length(rows)+1]] <- r
    }
  }
  out <- do.call(rbind, rows); out$cohort <- cohort
  out[, c("cohort", setdiff(names(out), "cohort"))]
}

## ============================================================================
## 3. Missing-data summary (justifies FIML)
## ============================================================================
missing_summary <- function(dt, vars, cohort){
  vars <- intersect(vars, names(dt))
  data.frame(cohort = cohort, variable = vars,
             n_total   = nrow(dt),
             n_present = sapply(vars, function(v) sum(!is.na(dt[[v]]))),
             n_missing = sapply(vars, function(v) sum(is.na(dt[[v]]))),
             pct_missing = round(sapply(vars, function(v) mean(is.na(dt[[v]]))) * 100, 1),
             row.names = NULL)
}

## ============================================================================
## 4. Sample-size flow (only the outcome waves this cohort has)
## ============================================================================
prs_present <- intersect(prs_vars, names(dataset))
n_genotyped <- if (length(prs_present) > 0)
  sum(rowSums(!is.na(dataset[prs_present])) > 0) else NA_integer_

outcome_labels <- c(intern_3y_c  = "Internalizing 3y",
                    intern_9y_c  = "Internalizing 9y",
                    intern_14y_c = "Internalizing 14y")

flow <- data.frame(cohort = cohort,
                   step = c("Total records", "With any PRS (genotyped)",
                            outcome_labels[outcome_vars]),
                   n = c(nrow(dataset), n_genotyped,
                         sapply(outcome_vars, function(v) sum(!is.na(dataset[[v]])))),
                   row.names = NULL)

## ============================================================================
## 5. Correlation matrices
## ============================================================================
## safe correlation: needs at least 2 variables, otherwise returns NULL
safe_cor <- function(dt, vars, label){
  vars <- intersect(vars, names(dt))
  if (length(vars) < 2){
    warning(label, ": fewer than 2 variables found in the dataset (",
            length(vars), "), correlation matrix skipped.")
    return(NULL)
  }
  num <- as.matrix(as.data.frame(lapply(dt[vars], function(x) suppressWarnings(as.numeric(x)))))
  colnames(num) <- vars
  round(cor(num, use = "pairwise.complete.obs"), 3)
}

missing_prs <- setdiff(prs_vars, names(dataset))
if (length(missing_prs) > 0)
  message("PRS variables not found in ", cohort, ": ", paste(missing_prs, collapse = ", "))

cor_mat <- safe_cor(dataset, c(mediator_vars, outcome_vars, prs_vars), "Main correlation matrix")
prs_cor <- safe_cor(dataset, prs_vars, "PRS correlation matrix")

## ============================================================================
## 6. Write descriptive outputs
## ============================================================================
write.csv(describe_continuous(dataset, continuous_vars, cohort),
          out_file("desc_continuous"), row.names = FALSE)
write.csv(describe_categorical(dataset, categorical_vars, cohort),
          out_file("desc_categorical"), row.names = FALSE)
write.csv(make_table1(dataset, continuous_vars, categorical_vars, sex_var, cohort),
          out_file("table1_bysex"), row.names = FALSE)
write.csv(missing_summary(dataset, c(continuous_vars, categorical_vars), cohort),
          out_file("missingness"), row.names = FALSE)
write.csv(flow, out_file("sample_flow"), row.names = FALSE)
if (!is.null(cor_mat))
  write.csv(data.frame(var = rownames(cor_mat), cor_mat),
            out_file("correlation_matrix"), row.names = FALSE)
if (!is.null(prs_cor))
  write.csv(data.frame(var = rownames(prs_cor), prs_cor),
            out_file("prs_correlations"), row.names = FALSE)

## Distribution check: histograms of every continuous variable to one PDF
pdf(out_file("distributions"), width = 10, height = 7)
op <- par(mfrow = c(3, 3))
for (v in continuous_vars){
  x <- suppressWarnings(as.numeric(dataset[[v]]))
  hist(x, main = v, xlab = v, breaks = 30, col = "grey80", border = "white")
}
par(op); dev.off()

message("Descriptives complete for cohort: ", cohort, "\nResults saved in: ", out_dir)
