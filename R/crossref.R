# Table that cross references litmz to provide putative IDs
# Author: Benjamin R. Gordon
# Data: 2019-03-11



table_crossref <- function(ppm = 50,
                           ref.type = c("generic", "endnote"),
                           save.impvars = FALSE,
                           save.matches = FALSE) {
  
  ref.type <- match.arg(ref.type)
  
  # retrieve important variables -----------------------------------------------
  mzpls_mod <- readRDS("./inst/extdata/mzpls_model.rds")
  mzrf_mod <- readRDS("./inst/extdata/mzrf_model.rds")
  
  mzpls_impvars <- caret::varImp(mzpls_mod, scale = TRUE)
  mzpls_impvars <- tibble::as_tibble(mzpls_impvars$importance)
  mzpls_impvars <-
    mzpls_impvars %>%
    mutate(mz = gsub("mz_", "", rownames(mzpls_impvars))) %>%
    rowwise() %>%
    transmute(
      mz = as.numeric(mz),# warnings arise from features with two decimal points
      importance = max(control, eT, eCO2, eCO2eT)) %>%
    arrange(desc(importance)) %>%
    slice(1:20)
  
  mzrf_impvars <- caret::varImp(mzrf_mod, scale = TRUE)
  mzrf_impvars <- tibble::as_tibble(mzrf_impvars$importance)
  mzrf_impvars <-
    mzrf_impvars %>%
    mutate(mz = gsub("mz_", "", rownames(mzrf_impvars))) %>%
    rowwise() %>%
    transmute(
      mz = as.numeric(mz),# warnings arise from features with two decimal points
      importance = max(control, eT, eCO2, eCO2eT)) %>%
    arrange(desc(importance)) %>%
    slice(1:20)
  
  # Add variables for ppm error ranges -----------------------------------------
  mzpls_impvars <-
    mzpls_impvars %>%
    rowwise() %>%
    mutate(
      model = "pls",
      mz_neutral = mz - 1.007276,
      mz_low = mz_neutral - (mz * ppm/10^6),
      mz_high = mz_neutral + (mz * ppm/10^6)) %>%
    select(model, everything()) %>%
    ungroup()
  
  mzrf_impvars <-
    mzrf_impvars %>%
    rowwise() %>%
    mutate(
      model = "rf",
      mz_neutral = mz - 1.007276,
      mz_low = mz_neutral - (mz * ppm/10^6),
      mz_high = mz_neutral + (mz * ppm/10^6)) %>%
    select(model, everything()) %>%
    ungroup()
  
  pls_rf_impvars <-
    bind_rows(mzpls_impvars, mzrf_impvars)
  
  # Cross reference with litmz -------------------------------------------------
  pls_rf_matches <-
    pls_rf_impvars %>%
    mutate(dummy = TRUE) %>%
    left_join(litmz %>% mutate(dummy = TRUE))  %>%
    filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
    select(-dummy,
           -mz_neutral,
           -mz_low,
           -mz_high)
  
  # References -----------------------------------------------------------------
  if (ref.type == "generic") {
    pls_rf_matches <-
      pls_rf_matches %>%
      mutate(ref = stringr::str_replace(endnote_ref, ",", ";")) %>%
      select(-endnote_ref)
  }
  
  else if (ref.type == "endnote") {
    pls_rf_matches <-
      pls_rf_matches %>%
      mutate(endnote_ref = stringr::str_replace(endnote_ref, ",", ";")) %>%
      select(-ref)
  }
  
  # Output and saves -----------------------------------------------------------
  if (save.impvars) {
    saveRDS(pls_rf_impvars, "./inst/extdata/pls_rf_impvars.rds")
  }
  
  if (save.matches) {
    readr::write_csv(pls_rf_matches, "./tables/important_variable_matches.txt")
  }
  
  impvars_matches <- list(imp_vars = pls_rf_impvars,
                          lit_matches = pls_rf_matches)
  
  message("duplicate mz values have .1 appended to the value and may produce
          some warnings")
  
  impvars_matches
  
}
