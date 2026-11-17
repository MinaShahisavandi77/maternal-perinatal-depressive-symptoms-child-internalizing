################################################################################
# Script for analysing data in MOBA cohort
# written by M.Shahisavandi 
# analist : Adrian Dahl Askelund
## Depression-PRS mediation models 
## Depression PRS (PPD / MDD / MDD-female/PPD_unique) -> perinatal mediator(s)
##   -> child internalizing
## STEP 1 (this script, per cohort): fit all models, extract the UNSTANDARDIZED
##   SEM parameter set + fit measures + vcov, and write cohort-tagged files.
## STEP 2 stratified analysis: fit all models, extract the UNSTANDARDIZED SEM
##   parameter set + fit measures + vcov, and write cohort-tagged files.
## step 3 adjust for history of depression and investigate the intrauterine effects of consequential mediation  
##
################################################################################
# we start with loading the cran libraries 

cran_pkgs <- c("remotes", "lavaan", "data.table")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
# Please have lavaan vesrion of (>= 0.6-12)
#======loading library===================
library(remotes)
remotes::install_github("MinaShahisavandi77/genmediation")
#=======================================
library(lavaan)
library(genmediation)
library(data.table)
## ============================================================================
## STEP 1 — FIT + EXTRACT (per cohort)
## ============================================================================
data_MOBA_std<-fread("/test/final_dataset_MOBA.csv")
results<- "/test/result"
## ---- 0. Cohort + shared settings --------------------------------------------
cohort    <- "MOBA"   
dataset <- standardise_vars(data_MOBA_std, vars = c("prs_PPD_unique_c","prs_PPD_unique_m","prs_PPD_unique_f","prs_PPD_c","prs_PPD_m","prs_PPD_f","prs_MDD_c","prs_MDD_m","prs_MDD_f","prs_MDD_female_c","prs_MDD_female_m","prs_MDD_female_f","dep_pre_m","dep_6m_m","dep_36m_m","intern_3y_c","intern_9y_c","intern_14y_c"))

### residual the prs and depression for mother 

# Define PC sets
pcs_c <- paste0("PC", 1:10, "_c")
pcs_m <- paste0("PC", 1:10, "_m")
pcs_f <- paste0("PC", 1:10, "_f")

# Helper function: regress var on PCs, return residuals (NA-safe, preserves row order)
residualize <- function(data, var, pcs) {
  form <- as.formula(paste(var, "~", paste(pcs, collapse = " + ")))
  model <- lm(form, data = data, na.action = na.exclude)
  resid(model)
}

# --- Depression variables (mother), residualized on mother's PCs ---
dep_vars <- c("dep_pre_m", "dep_6m_m", "dep_36m_m")

for (v in dep_vars) {
  dataset[[paste0(v, "_res")]] <- residualize(dataset, v, pcs_m)
}

# --- PRS variables, residualized on their own group's PCs (by _c/_m/_f suffix) ---
prs_vars <- c("prs_MDD_c", "prs_PPD_c", "prs_MDD_m", "prs_PPD_m", "prs_MDD_f", "prs_PPD_f",
              "prs_PPD_unique_f", "prs_PPD_unique_m", "prs_PPD_unique_c",
              "prs_MDD_female_c", "prs_MDD_female_m", "prs_MDD_female_f")

for (v in prs_vars) {
  suffix <- sub(".*_", "", v)          # last part after final underscore: c, m, or f
  pcs <- switch(suffix,
                c = pcs_c,
                m = pcs_m,
                f = pcs_f,
                stop(paste("Unrecognized suffix for variable:", v)))
  dataset[[paste0(v, "_res")]] <- residualize(dataset, v, pcs)
}

## Depression PRS traits, each with child / mother / father columns
prs_sets <- list(
  PPD       = list(c = "prs_PPD_c_res",  p1 = "prs_PPD_m_res",  p2 = "prs_PPD_f_res"),
  MDD       = list(c = "prs_MDD_c_res",  p1 = "prs_MDD_m_res",  p2 = "prs_MDD_f_res"),
  MDDfemale = list(c = "prs_MDD_female_c_res", p1 = "prs_MDD_female_m_res", p2 = "prs_MDD_female_f_res"),
  PDDunique=list(c="prs_PPD_unique_c_res",p1="prs_PPD_unique_m_res",p2="prs_PPD_unique_f_res")
)

## Mediator blocks (PRS comes from prs_sets, not from here)
exposures <- list(
  depression  = c(w1 = "dep_pre_m_res", w2 = "dep_6m_m_res", w3 = "dep_36m_m_res"),
  bmi         = c(w1 = "BMI_m"),
  birthweight = c(w1 = "birthweight_c"),
  smoke=c(w1="smoke_m"),
  gestational=c(w1="age_gestational_c")
)
## Outcome variable + its age covariate, per timepoint
outcomes <- list(
  y3y  = list(var = "intern_3y_c",age = "age_3y_c"),
  y9y  = list(var = "intern_9y_c",age = "age_9y_c"),
  y14y = list(var = "intern_14y_c",age = "age_14y_c")
)

cov_mediator_base <- c("age_pregnancy_m", "gender_c")   
cov_outcome_base  <- c("gender_c", "age_pregnancy_m") 

setwd(results)

## ---- 1. Extraction helpers --------------------------------------------------

## run_mediation may return the lavaan object directly or a wrapper list.
## Adjust this if your wrapper stores it under a different slot.
get_lav <- function(x) {
  if (inherits(x, "lavaan"))      x
  else if (!is.null(x$fit))       x$fit
  else if (!is.null(x$lavaan))    x$lavaan
  else stop("Can't locate a lavaan object inside the run_mediation result.")
}

## UNSTANDARDIZED parameter table only: est, SE (delta-method), z, p, CIs,
## R-square rows, N, convergence/SE flags, and a cross-cohort key.
extract_sem_full <- function(fit, tag, cohort, prs, exposure, outcome) {
  lav <- get_lav(fit)
  
  converged <- isTRUE(lavInspect(lav, "converged"))
  se_opt    <- lavInspect(lav, "options")$se
  if (!converged)               warning(tag, ": model did NOT converge — SEs will be NA")
  if (identical(se_opt, "none")) warning(tag, ": se='none' — no SEs were computed")
  
  ## Unstandardized parameter table ONLY (no standardized solution).
  out <- parameterEstimates(lav, ci = TRUE, standardized = FALSE, rsquare = TRUE)
  
  out$param_id  <- ifelse(nzchar(out$label), out$label, paste(out$lhs, out$op, out$rhs))
  out$N         <- lavInspect(lav, "ntotal")
  out$ngroups   <- lavInspect(lav, "ngroups")
  out$converged <- converged
  out$se_ok     <- !is.na(out$se)
  out$meta_use  <- out$op %in% c("~", ":=")   # rows you actually pool in meta-analysis
  
  cbind(cohort = cohort, model = tag, prs = prs,
        exposure = exposure, outcome = outcome, out, row.names = NULL)
}

## All global fit indices (long format so models with differing measures rbind).
extract_fit_measures <- function(fit, tag, cohort, prs, exposure, outcome) {
  lav <- get_lav(fit)
  fm  <- fitMeasures(lav)
  data.frame(cohort = cohort, model = tag, prs = prs,
             exposure = exposure, outcome = outcome,
             measure = names(fm), value = as.numeric(fm),
             row.names = NULL, stringsAsFactors = FALSE)
}

## ---- 2. Fit every PRS x mediator x outcome combination ----------------------
params      <- list()
fit_indices <- list()

for (prs_name in names(prs_sets)) {
  ps <- prs_sets[[prs_name]]
  for (exp_name in names(exposures)) {
    meds <- exposures[[exp_name]]
    for (out_name in names(outcomes)) {
      out <- outcomes[[out_name]]
      tag <- paste(prs_name, exp_name, out_name, sep = "__")
      message("Fitting: ", cohort, " | ", tag)
      
      model_syntax <- build_mediation_model(
        outcome             = out$var,
        mediator_vars       = meds,
        prs_child           = ps$c,
        prs_parent1         = ps$p1,
        prs_parent2         = ps$p2,
        covariates_mediator = cov_mediator_base,
        covariates_outcome  = c(cov_outcome_base, out$age)
      )
      
      fit <- tryCatch(
        run_mediation(data = dataset, model = model_syntax,missing = "fiml.x",
                      cluster = "family_id",
                      seed    = 123),
        error = function(e) { warning("FAILED ", tag, ": ", conditionMessage(e)); NULL }
      )
      if (is.null(fit)) next
      
      params[[tag]]      <- extract_sem_full(fit, tag, cohort, prs_name, exp_name, out_name)
      fit_indices[[tag]] <- extract_fit_measures(fit, tag, cohort, prs_name, exp_name, out_name)
      
      ## parameter covariance matrix (for multivariate meta later)
      saveRDS(vcov(get_lav(fit)), paste0("vcov_", cohort, "_", tag, ".rds"))
    }
  }
}

## ---- 3. Write cohort-tagged outputs -----------------------------------------
params_all <- do.call(rbind, params)
fit_all    <- do.call(rbind, fit_indices)
# if you get this error is harmless " mean cluster size is 1.xxx.This means that many clusters only contain a single observation".
write.csv(params_all, paste0("sem_params_",      cohort, ".csv"), row.names = FALSE)
write.csv(fit_all,    paste0("sem_fitmeasures_", cohort, ".csv"), row.names = FALSE)

## ============================================================================
## STEP 2 — SENSITIVITY ANALYSIS: sex-stratified (boys = 1, girls = 2)
## ============================================================================
cov_mediator_sex <- setdiff(cov_mediator_base, "gender_c")
cov_outcome_sex  <- setdiff(cov_outcome_base,  "gender_c")

sex_groups <- c(boys = 1, girls = 2)

for (grp in names(sex_groups)) {
  
  keep    <- !is.na(dataset$gender_c) & dataset$gender_c == sex_groups[[grp]]
  dat_sex <- dataset[keep, , drop = FALSE]
  message("\n=== Sensitivity: ", grp, " (n = ", nrow(dat_sex), ") ===")
  
  params_s      <- list()
  fit_indices_s <- list()
  
  for (prs_name in names(prs_sets)) {
    ps <- prs_sets[[prs_name]]
    for (exp_name in names(exposures)) {
      meds <- exposures[[exp_name]]
      for (out_name in names(outcomes)) {
        out <- outcomes[[out_name]]
        tag <- paste(prs_name, exp_name, out_name, sep = "__")
        message("Fitting [", grp, "]: ", cohort, " | ", tag)
        
        model_syntax <- build_mediation_model(
          outcome             = out$var,
          mediator_vars       = meds,
          prs_child           = ps$c,
          prs_parent1         = ps$p1,
          prs_parent2         = ps$p2,
          covariates_mediator = cov_mediator_sex,
          covariates_outcome  = c(cov_outcome_sex, out$age)
        )
        
        fit <- tryCatch(
          run_mediation(data = dat_sex, model = model_syntax,missing = "fiml.x",
                        cluster = "family_id",
                        seed    = 123),
          error = function(e) { warning("FAILED ", grp, " ", tag, ": ", conditionMessage(e)); NULL }
        )
        if (is.null(fit)) next
        
        pr <- extract_sem_full(fit, tag, cohort, prs_name, exp_name, out_name)
        pr$group <- grp
        fm <- extract_fit_measures(fit, tag, cohort, prs_name, exp_name, out_name)
        fm$group <- grp
        params_s[[tag]]      <- pr
        fit_indices_s[[tag]] <- fm
        
        saveRDS(vcov(get_lav(fit)),
                paste0("vcov_", cohort, "_", grp, "_", tag, ".rds"))
      }
    }
  }
  
  write.csv(do.call(rbind, params_s),
            paste0("sem_params_",      cohort, "_", grp, ".csv"), row.names = FALSE)
  write.csv(do.call(rbind, fit_indices_s),
            paste0("sem_fitmeasures_", cohort, "_", grp, ".csv"), row.names = FALSE)
}

## ============================================================================
## STEP 3 — EXTENDED CHAIN: upstream depression mediator(s)
## If `history_dep` exists, it enters as the first mediator (before prenatal
## depression). If it does NOT, fall back to prenatal depression as the upstream
## mediator (bmi / birthweight blocks only).
## ============================================================================
if ("history_dep_m" %in% names(dataset)) {
  dataset <- standardise_vars(dataset, vars = "history_dep_m")
  
  exposures_hist <- list(
    depression  = c(w1 = "history_dep_m", w2 = "dep_pre_m_res", w3 = "dep_6m_m_res", w4 = "dep_36m_m_res")
  )
  analysis_tag <- "history_dep_m"
  message("STEP 3: history_dep_m available — extended chain WITH history_dep.")
  
} 

params_h      <- list()
fit_indices_h <- list()

for (prs_name in names(prs_sets)) {
  ps <- prs_sets[[prs_name]]
  for (exp_name in names(exposures_hist)) {
    meds <- exposures_hist[[exp_name]]
    for (out_name in names(outcomes)) {
      out <- outcomes[[out_name]]
      tag <- paste(prs_name, exp_name, out_name, sep = "__")
      message("Fitting [histdep]: ", cohort, " | ", tag)
      
      model_syntax <- build_mediation_model(
        outcome             = out$var,
        mediator_vars       = meds,
        prs_child           = ps$c,
        prs_parent1         = ps$p1,
        prs_parent2         = ps$p2,
        covariates_mediator = cov_mediator_base,
        covariates_outcome  = c(cov_outcome_base, out$age)
      )
      
      fit <- tryCatch(
        run_mediation(data = dataset, model = model_syntax,missing = "fiml.x",
                      cluster = "family_id",
                      seed    = 123),
        error = function(e) { warning("FAILED histdep ", tag, ": ", conditionMessage(e)); NULL }
      )
      if (is.null(fit)) next
      
      pr <- extract_sem_full(fit, tag, cohort, prs_name, exp_name, out_name)
      pr$analysis <- analysis_tag
      fm <- extract_fit_measures(fit, tag, cohort, prs_name, exp_name, out_name)
      fm$analysis <- analysis_tag
      params_h[[tag]]      <- pr
      fit_indices_h[[tag]] <- fm
      
      saveRDS(vcov(get_lav(fit)),
              paste0("vcov_", cohort, "_histdep_", tag, ".rds"))
    }
  }
}

write.csv(do.call(rbind, params_h),
          paste0("sem_params_",      cohort, "_histdep.csv"), row.names = FALSE)
write.csv(do.call(rbind, fit_indices_h),
          paste0("sem_fitmeasures_", cohort, "_histdep.csv"), row.names = FALSE)



