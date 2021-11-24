# This script provides plotting, risk modeling and prediction tools for the app

# tools ----

  require(plyr)
  require(tidyverse)
  require(cowplot)
  require(gridExtra)
  require(stringi)

  ## auxiliary scripts

  scrpt_paths <- c('./tools/modeling_tools.R', 
                   './tools/km_toolbox.R')

  src_result <- try(scrpt_paths %>% 
                      walk(source), 
                    silent = T) ## if launched via app

  if(any(class(src_result) == 'try-error')) {
    
    scrpt_paths <- stri_replace(scrpt_paths, 
                                fixed = './', 
                                replacement = './PAH Scorer/') ## if run standalone
    
    src_result <- try(scrpt_paths %>% 
                        walk(source), 
                      silent = T)
    
  }
  
# reading the tables with values of both prediction scores and outcomes in the establishment cohort -----
  
  score_tbl <- try(read_tsv('./data/score_tbl.tsv'), 
                   silent = T) ## if run in the app
  
  if(any(class(score_tbl) == 'try-error')) {
    
    score_tbl <- try(read_tsv('./PAH Scorer/data/score_tbl.tsv'), 
                     silent = T) ## if run standalone
    
  }

# functions ----
  
  ## general functions
  
  count_feature <- function(inp_tbl, var_to_count, remove_na = T) {
    
    ## calculates the percentage and number of participants with/without the given feature
    
    if(remove_na) {
      
      count_tbl <- inp_tbl %>% 
        filter(!is.na(.data[[var_to_count]]), 
               .data[[var_to_count]] != 'no_answer')
      
    } else {
      
      count_tbl <- inp_tbl
      
    }
    
    feature_counts <- count_tbl %>% 
      count(.data[[var_to_count]]) %>% 
      mutate(percent = n/sum(n) * 100, 
             total_n = sum(n))
    
    return(feature_counts)
    
  }
  
  ## score calculation functions

  calculate_risk <- function(inp_data, inp_model, type = 'risk', ...) {
    
    ## calculates HR or other response for the given table with cox model variables 
    ## and the input cox or glm model object
    
    mod_prediction <- predict(inp_model, 
                              type = type, 
                              newdata = inp_data, 
                              ...)
    
    return(mod_prediction)
    
  }
  
  risk_stats_tbl <- function(score_val, risk_class, predict_res) {
    
    ## creates a table with 5-year risk modeling summary
    
    return(tibble(score = score_val, 
                  risk_class = risk_class, 
                  mort_risk = predict_res$fit, 
                  lower_ci = predict_res$fit + predict_res$se.fit* qnorm(0.025), 
                  upper_ci = predict_res$fit + predict_res$se.fit* qnorm(0.975)))
    
  }
  
  calculate_sign2525 <- function(Gender, age_fc, cardiac_index, mPAP, score_only = F, ...) {
    
    ## calculates the model_10967 score (HR from the respective Cox model)
    ## and the predicted 5-year mortality with 95% CI obtained from normal distribution
    
    age_strat <- ifelse(age_fc > 60, 2.1009330, 0)
    
    Gender_strat <- ifelse(Gender == 'female', -0.7506657, 0)

    mPAP_strat <- ifelse(mPAP >= 49, 0.9143866, 0)

    ci_strat <- cut(cardiac_index, 
                    c(-Inf, 2, 2.5, Inf), 
                    c(0, -0.8219363, -1.2468615)) %>% 
      as.character %>% 
      as.numeric
    
    inp_tbl <- tibble(sign2525 = Gender_strat + age_strat + ci_strat + mPAP_strat)
    
    if(score_only) {
      
      return(inp_tbl$sign2525)
      
    }
    
    risk <- calculate_risk(inp_data = inp_tbl, 
                           inp_model = glm_models$sign2525, 
                           type = 'response', 
                           se.fit = T)
    
    return(risk_stats_tbl(score_val = inp_tbl$sign2525, 
                          risk_class = NA, 
                          predict_res = risk))
    
  }
  
  calculate_enh_french <- function(FRENCH3p, age_fc, RAA, ...) {
    
    ## calculates the enchanced FPHR3p score as described by Sonnweber et al.
    ## and predicts 5-year mortality
    
    age_strat <- cut(age_fc, 
                     c(-1, 40, 65, 120), 
                     c(0, 1, 2)) %>% 
      as.character %>% 
      as.numeric
    
    RAA_strat <- cut(RAA, 
                     c(-1, 18, 26,  100), 
                     c(0, 1, 2)) %>% 
      as.character %>% 
      as.numeric
    
    inp_tbl <- tibble(french_enh = FRENCH3p + age_strat + RAA_strat)
    
    risk <- calculate_risk(inp_data = inp_tbl, 
                           inp_model = glm_models$french_enh, 
                           type = 'response', 
                           se.fit = T)

    return(risk_stats_tbl(score_val = inp_tbl$french_enh, 
                          risk_class = NA, 
                          predict_res = risk))
    
  }
  
  calculate_french3p <- function(WHOFc, SMWD, NTproBNP, ...) {
    
    ## calculates FRENCH3p score (number of low risk criteria missed) 
    ## and 5-year mortality stats
    
    WHOFc_strat <- ifelse(WHOFc %in% c(1, 2), 0, 1)
    
    SMWD_strat <- ifelse(SMWD <= 440, 1, 0)
    
    NTproBNP_strat <- ifelse(NTproBNP < 300, 0, 1)
    
    inp_tbl <- tibble(FRENCH3p = WHOFc_strat + SMWD_strat + NTproBNP_strat)

    risk <- calculate_risk(inp_data = inp_tbl, 
                           inp_model = glm_models$FRENCH3p, 
                           type = 'response', 
                           se.fit = T)
    
    risk_class <- paste(inp_tbl$FRENCH3p, 'low risk criteria missed')
    
    return(risk_stats_tbl(score_val = inp_tbl$FRENCH3p, 
                          risk_class = risk_class, 
                          predict_res = risk))
    
  }
  
  calculate_french4p <- function(WHOFc, SMWD, mRAP, cardiac_index, ...) {
    
    ## calculates FRENCH4p score (number of low risk criteria missed) 
    ## and 5-year mortality stats
    
    WHOFc_strat <- ifelse(WHOFc %in% c(1, 2), 0, 1)
    
    SMWD_strat <- ifelse(SMWD > 440, 0, 1)
    
    mRAP_strat <- ifelse(mRAP < 8, 0, 1)
    
    CI_strat <- ifelse(cardiac_index > 2.5, 0, 1)
    
    inp_tbl <- tibble(FRENCH4p = WHOFc_strat + SMWD_strat + mRAP_strat + CI_strat)
    
    risk <- calculate_risk(inp_data = inp_tbl, 
                           inp_model = glm_models$FRENCH4p, 
                           type = 'response', 
                           se.fit = T)
    
    risk_class <- paste(inp_tbl$FRENCH4p, 'low risk criteria missed')
    
    
    return(risk_stats_tbl(score_val =inp_tbl$FRENCH4p, 
                          risk_class = risk_class, 
                          predict_res = risk))
    
  }
  
  calculate_compera <- function(WHOFc, SMWD, NTproBNP, mRAP, cardiac_index, SvO2, ...) {
    
    ## calculates FRENCH4p score (number of low risk criteria missed) 
    ## and 5-year mortality stats
    
    WHOFc_strat <- cut(WHOFc, 
                       c(-1, 2.5, 3.5, 5), 
                       c(1, 2, 3)) %>% 
      as.character %>% 
      as.numeric
    
    SMWD_strat <- cut(SMWD, 
                      c(-1, 165, 440, 10000), 
                      c(3, 2, 1), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    NTproBNP_strat <- cut(NTproBNP, 
                          c(-1, 300, 1400, 300000), 
                          c(1, 2, 3), 
                          right = F) %>% 
      as.character %>% 
      as.numeric
    
    mRAP_strat <- cut(mRAP, 
                      c(-1, 8, 14, 100), 
                      c(1, 2, 3), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    mRAP_strat <- cut(mRAP, 
                      c(-1, 8, 14, 100), 
                      c(1, 2, 3), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    CI_strat <- cut(cardiac_index, 
                    c(-1, 2, 2.5, 100), 
                    c(3, 2, 1), 
                    right = F) %>% 
      as.character %>% 
      as.numeric
    
    SvO2_strat <- cut(SvO2, 
                      c(-1, 60, 65, 110), 
                      c(3, 2, 1), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    inp_tbl <- tibble(Compera = mean(c(WHOFc_strat, 
                                       SMWD_strat, 
                                       NTproBNP_strat, 
                                       mRAP_strat, 
                                       CI_strat, 
                                       SvO2_strat)) %>% 
                        round)
    
    risk <- calculate_risk(inp_data = inp_tbl, 
                           inp_model = glm_models$Compera, 
                           type = 'response', 
                           se.fit = T)
    
    return(risk_stats_tbl(score_val = inp_tbl$Compera, 
                          risk_class = cut(inp_tbl$Compera, 
                                           c(0, 1.5, 2.5, 4), 
                                           c('Low', 'Intermediate', 'High')), 
                          predict_res = risk))
    
  }
  
  calculate_spahr <- function(WHOFc, SMWD, NTproBNP, mRAP, cardiac_index, SvO2, RAA, pericardial_effusion, ...) {
    
    ## calculates FRENCH4p score (number of low risk criteria missed) 
    ## and 5-year mortality stats
    
    WHOFc_strat <- cut(WHOFc, 
                       c(-1, 2.5, 3.5, 5), 
                       c(1, 2, 3)) %>% 
      as.character %>% 
      as.numeric
    
    SMWD_strat <- cut(SMWD, 
                      c(-1, 165, 440, 10000), 
                      c(3, 2, 1), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    NTproBNP_strat <- cut(NTproBNP, 
                          c(-1, 300, 1400, 300000), 
                          c(1, 2, 3), 
                          right = F) %>% 
      as.character %>% 
      as.numeric
    
    mRAP_strat <- cut(mRAP, 
                      c(-1, 8, 14, 100), 
                      c(1, 2, 3), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    mRAP_strat <- cut(mRAP, 
                      c(-1, 8, 14, 100), 
                      c(1, 2, 3), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    CI_strat <- cut(cardiac_index, 
                    c(-1, 2, 2.5, 100), 
                    c(3, 2, 1), 
                    right = F) %>% 
      as.character %>% 
      as.numeric
    
    SvO2_strat <- cut(SvO2, 
                      c(-1, 60, 65, 110), 
                      c(3, 2, 1), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    RAA_strat <- cut(RAA, 
                     c(-1, 18, 26, 200), 
                     c(1, 2, 3), 
                     right = F) %>% 
      as.character %>% 
      as.numeric
    
    pericardial_eff_rec <- ifelse(pericardial_effusion == 'no', 0, 2)
    
    inp_tbl <- tibble(SPAHR = mean(c(WHOFc_strat, 
                                     SMWD_strat, 
                                     NTproBNP_strat, 
                                     mRAP_strat, 
                                     CI_strat, 
                                     SvO2_strat, 
                                     RAA_strat, 
                                     pericardial_eff_rec)) %>% 
                        round)
    
    risk <- calculate_risk(inp_data = inp_tbl, 
                           inp_model = glm_models$SPAHR, 
                           type = 'response', 
                           se.fit = T)
    
    return(risk_stats_tbl(score_val =inp_tbl$SPAHR, 
                          cut(inp_tbl$SPAHR, 
                              c(0, 1.5, 2.5, 4), 
                              c('Low', 'Intermediate', 'High')), 
                          predict_res = risk))
    
  }
  
  calculate_mrasp <- function(WHOFc, SMWD, NTproBNP, RAA, ...) {
    
    ## calculates risk classes of mRASP and the 5-year mortality risk
    
    WHOFc_strat <- cut(WHOFc, 
                       c(-1, 2.5, 3.5, 5), 
                       c(0, 1, 2)) %>% 
      as.character %>% 
      as.numeric
    
    SMWD_strat <- cut(SMWD, 
                      c(-1, 165, 440, 10000), 
                      c(2, 1, 0), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    NTproBNP_strat <- cut(NTproBNP, 
                          c(-1, 300, 1400, 300000), 
                          c(0, 1, 2), 
                          right = F) %>% 
      as.character %>% 
      as.numeric
    
    RAA_strat <- cut(RAA, 
                     c(-1, 18, 26, 200), 
                     c(0, 1, 2), 
                     right = F) %>% 
      as.character %>% 
      as.numeric
    
    inp_tbl <- tibble(mRASP_score = WHOFc_strat + 
                        SMWD_strat + 
                        NTproBNP_strat + 
                        RAA_strat) %>% 
      mutate(mRASP = cut(mRASP_score, 
                         c(-1, 2.5, 5.5, 10), 
                         c(0, 1, 2)) %>% 
               as.character %>% 
               as.numeric)
    
    risk <- calculate_risk(inp_data = inp_tbl, 
                           inp_model = glm_models$mRASP, 
                           type = 'response', 
                           se.fit = T)
    
    return(risk_stats_tbl(score_val = inp_tbl$mRASP_score, 
                          risk_class = cut(inp_tbl$mRASP_score, 
                                           c(-1, 2.5, 5.5, 10), 
                                           c('Low', 'Intermediate', 'High')), 
                          predict_res = risk))
    
    
  }
  
  calculate_reveal_lite <- function(WHOFc, SMWD, NTproBNP, eGFR, SBP, HR, ...) {
    
    ## calculated REVEAL Lite
    
    WHOFc_strat <- cut(WHOFc, 
                       c(-1, 1.5, 2.5, 3.5, 5), 
                       c(-1, 0, 1, 2)) %>% 
      as.character %>% 
      as.numeric
    
    SMWD_strat <- cut(SMWD, 
                      c(-1, 165, 320, 440, 10000), 
                      c(1, 0, -1, -2), 
                      right = F) %>% 
      as.character %>% 
      as.numeric
    
    NTproBNP_strat <- cut(NTproBNP, 
                          c(-1, 300, 1100, 300000), 
                          c(-2, 0, 2), 
                          right = F) %>% 
      as.character %>% 
      as.numeric
    
    eGFR_strat <- ifelse(eGFR < 60, 1, 0)
    
    SBP_strat <- ifelse(SBP < 110, 1, 0) 
    
    HR_strat <- ifelse(HR > 96, 1, 0)
    
    
    inp_tbl <- tibble(reveal_lite_score = WHOFc_strat + 
                        SMWD_strat + 
                        NTproBNP_strat + 
                        eGFR_strat + 
                        SBP_strat + 
                        HR_strat + 6) %>% 
      mutate(reveal_lite = cut(reveal_lite_score, 
                                      c(-1, 5.5, 7.5, 16), 
                                      c(1, 2, 3)) %>% 
               as.character %>% 
               as.numeric)
    
    risk <- calculate_risk(inp_data = inp_tbl, 
                           inp_model = glm_models$reveal_lite, 
                           type = 'response', 
                           se.fit = T)
    
    return(risk_stats_tbl(score_val = inp_tbl$reveal_lite_score, 
                          risk_class = cut(inp_tbl$reveal_lite, 
                                           c(-1, 1.5, 2.5, 10), 
                                           c('Low', 'Intermediate', 'High')), 
                          predict_res = risk))
    
    
  }
  
  make_input_tbl <- function(PatientID = '', Gender, mPAP, age_fc, RAA, WHOFc, SMWD, 
                             NTproBNP, mRAP, cardiac_index, SvO2, pericardial_effusion, eGFR, SBP, HR) {
    
    ## makes a table with input parameters
    
    inp_param_tbl <- tibble(Paramater = c('Patient ID', 
                                          'Sex', 
                                          'Age', 
                                          'Six Minute Walking Distance, m', 
                                          'WHO Functional Class', 
                                          'Systolic Blood Pressure, mmHg', 
                                          'Heart Rate, bpm', 
                                          'Paricardial Effusion', 
                                          'Right Atrial Area, sq. cm', 
                                          'NT-proBNP, pg/mL', 
                                          'eGFR, mL/min/1.73 sq m', 
                                          'Mixed Venous Oxygen Saturation, %', 
                                          'Cardiac Index', 
                                          'Mean Right Atrial Pressure, mmHg', 
                                          'Mean Pulmonary Artherial Pressure, mmHg'), 
                            Value = list(PatientID, 
                                         Gender, 
                                         age_fc, 
                                         SMWD,
                                         car::recode(as.character(WHOFc), "'1' = 'I'; '2' = 'II'; '3' = 'III'; '4' = 'IV'"), 
                                         SBP, 
                                         HR, 
                                         pericardial_effusion, 
                                         RAA, 
                                         NTproBNP, 
                                         eGFR, 
                                         SvO2, 
                                         cardiac_index, 
                                         mRAP, 
                                         mPAP) %>% 
                              map(function(x) if(class(x) == 'numeric') signif(x, 2) else x) %>% 
                              map_chr(as.character))
    
    
    return(inp_param_tbl)
    
  }
  
  make_score_tbl <- function(Gender, mPAP, age_fc, RAA, WHOFc, SMWD, 
                             NTproBNP, mRAP, cardiac_index, SvO2, pericardial_effusion, eGFR, SBP, HR) {
    
    ## makes a table with score values, risk categories and predicted 5-year mortality
    ## for the two experimental scores (signature 2525 and Sonnweber Enhanced French)
    ## French3p, 4p, SAPHR, COMPERA and mRASP
    
    FRENCH3p <- calculate_french3p(WHOFc = WHOFc, 
                                   SMWD = SMWD, 
                                   NTproBNP = NTproBNP)

    other_results <- list(sign2525 = calculate_sign2525, 
                          enh_french = calculate_enh_french, 
                          FRENCH4p = calculate_french4p, 
                          SPAHR = calculate_spahr, 
                          Compera = calculate_compera, 
                          mRASP = calculate_mrasp, 
                          reveal_lite = calculate_reveal_lite) %>% 
      map(function(x) x(Gender = Gender, 
                        mPAP = mPAP, 
                        age_fc = age_fc, 
                        RAA = RAA, 
                        WHOFc = WHOFc, 
                        SMWD = SMWD, 
                        NTproBNP = NTproBNP, 
                        mRAP = mRAP, 
                        cardiac_index = cardiac_index, 
                        SvO2 = SvO2, 
                        pericardial_effusion = pericardial_effusion, 
                        FRENCH3p = FRENCH3p$score, 
                        eGFR = eGFR, 
                        SBP = SBP, 
                        HR = HR))
    
    res_tbl <- other_results %>% 
      c(list(FRENCH3p = FRENCH3p))
    
    res_tbl <- res_tbl[names(proj_globs$score_labs)] %>% 
      map2_dfr(., 
               names(.), 
               function(x, y) mutate(x, risk_scale = y)) %>% 
      mutate(score_label = proj_globs$score_labs[risk_scale], 
             five_mort_lab = paste(signif(mort_risk, 2) * 100, 
                                   ' [95%CI: ', 
                                   signif(lower_ci, 2) * 100, 
                                   ', ', 
                                   signif(upper_ci, 2) * 100, 
                                   ']', sep = ''), 
             reference = proj_globs$score_refs[risk_scale]) %>% 
      select(risk_scale, 
             score_label, 
             reference, 
             score, 
             risk_class, 
             mort_risk, 
             lower_ci, 
             upper_ci, 
             five_mort_lab)
    
    return(res_tbl)
    
  }
  
  ## report function
  
  render_report <- function(inp_param_tbl, score_summ_tbl, path_to_save) {
    
    ## makes a simple ggplot/cowplot report to be saved by the user
    ## score value and predicted 5-year mortality risk
    
    ## input parameter table 
    
    hj <- matrix(c(0, 1), 
                 ncol = 2, 
                 nrow = nrow(inp_param_tbl), 
                 byrow = TRUE) ## justification of the table head
    
    x <- matrix(c(0.05, 0.95), 
                ncol = 2, 
                nrow = nrow(inp_param_tbl), 
                byrow = TRUE) ## justification of the table head
    
    just_mtx <- cbind(rep(0, nrow(inp_param_tbl)), 
                      rep(1, nrow(inp_param_tbl))) ## justification of the table body
    
    x_mtx <- cbind(rep(0.05, nrow(inp_param_tbl)), 
                   rep(0.95, nrow(inp_param_tbl))) ## justification of the table body
    
    theme_1 <- ttheme_minimal(core = list(fg_params = list(hjust = as.vector(just_mtx), 
                                                           x = as.vector(x_mtx), 
                                                           fontsize = 9, 
                                                           fontface = 'plain')),
                              colhead = list(fg_params = list(fontsize = 9, 
                                                              fontface = 'bold', 
                                                              hjust = as.vector(hj), 
                                                              x = as.vector(x)), 
                                             bg_params = list(fill = 'gray96')))
    
    inp_grob <- tableGrob(inp_param_tbl, 
                          theme = theme_1, 
                          rows = NULL)
    
    ## summary table with score values
    
    nice_tbl <- score_summ_tbl %>% 
      select(score_label, 
             reference, 
             score, 
             risk_class, 
             five_mort_lab) %>% 
      mutate(score = signif(score, 2)) %>% 
      set_names(c('Risk assessment tool', 
                  'Reference', 
                  'Score', 
                  'Risk class', 
                  'Predicted 5-year mortality risk, %'))
    
    hj <- matrix(c(0, 1), 
                 ncol = 2, 
                 nrow = nrow(nice_tbl), 
                 byrow = TRUE) ## justification of the table head
    
    x <- matrix(c(0.05, 0.95), 
                ncol = 2, 
                nrow = nrow(nice_tbl), 
                byrow = TRUE) ## justification of the table head
    
    just_mtx <- cbind(rep(0, nrow(nice_tbl)), 
                      rep(1, nrow(nice_tbl)), 
                      rep(1, nrow(nice_tbl)), 
                      rep(1, nrow(nice_tbl)), 
                      rep(1, nrow(nice_tbl))) ## justification of the table body
    
    x_mtx <- cbind(rep(0.05, nrow(nice_tbl)), 
                   rep(0.95, nrow(nice_tbl)), 
                   rep(0.95, nrow(nice_tbl)), 
                   rep(0.95, nrow(nice_tbl)), 
                   rep(0.95, nrow(nice_tbl))) ## justification of the table body
    
    
    theme_1 <- ttheme_minimal(core = list(fg_params = list(hjust = as.vector(just_mtx), 
                                                           x = as.vector(x_mtx), 
                                                           fontsize = 9, 
                                                           fontface = 'plain')),
                              colhead = list(fg_params = list(fontsize = 9, 
                                                              fontface = 'bold', 
                                                              hjust = as.vector(hj), 
                                                              x = as.vector(x)), 
                                             bg_params = list(fill = 'gray96')))
    
    
    table_grob <- tableGrob(nice_tbl, 
                            theme = theme_1, 
                            rows = NULL)
    
    ## putting everything into a report
    
    report <- list()
    
    banner <- plot_grid(ggdraw() + 
                          draw_image('./www/mui_logo.png'), 
                        nrow = 1)
    
    report$title_panel <- plot_grid(ggdraw() + 
                                      draw_text('PAH Risk Assessment Tools', 
                                                size = 14, 
                                                fontface = 'bold', 
                                                color = '#800000', 
                                                x = 0.05, 
                                                hjust = 0), 
                                    ggdraw() + 
                                      draw_text('Score Calulation and Risk Prediction Summary', 
                                                size = 12, 
                                                x = 0.05, 
                                                hjust = 0),  
                                    ncol = 1, 
                                    rel_heights = c(1, 1)) %>% 
      plot_grid(.,
                banner, 
                rel_widths = c(1, 0.7))
    
    report$inp_paramaters <- plot_grid(ggdraw() + 
                                         draw_text('Input parameters', 
                                                   size = 11, 
                                                   fontface = 'bold', 
                                                   x = 0.05, 
                                                   hjust = 0), 
                                       ggdraw() + 
                                         draw_grob(inp_grob, 
                                                   x = 0.01, 
                                                   hjust = 0), 
                                       ncol = 1, 
                                       rel_heights = c(0.1, 1))
    
    report$summ_table <- plot_grid(ggdraw() + 
                                     draw_text('Score values and predicted 5-year mortality risk', 
                                               size = 11, 
                                               fontface = 'bold', 
                                               x = 0.05, 
                                               hjust = 0), 
                                   ggdraw() + 
                                     draw_grob(table_grob, 
                                               x = 0.01, 
                                               hjust = 0), 
                                   ggdraw() + 
                                     draw_text(proj_globs$method_text, 
                                               size = 9, 
                                               x = 0.05, 
                                               hjust = 0), 
                                   ncol = 1, 
                                   rel_heights = c(0.1, 1, 0.2))
    
    report$references <- plot_grid(ggdraw() + 
                                     draw_text('References', 
                                               size = 11, 
                                               fontface = 'bold', 
                                               x = 0.05, 
                                               hjust = 0), 
                                   ggdraw() + 
                                     draw_text(paste(proj_globs$score_citations_lite, 
                                                     collapse = '\n'), 
                                               size = 9, 
                                               x = 0.05, 
                                               hjust = 0), 
                                   ncol = 1, 
                                   rel_heights = c(0.05, 1))
    
    
    report_file <- plot_grid(report$title_panel + 
                               theme(plot.margin = margin(t = 15, 
                                                          l = 10, 
                                                          r = 10, 
                                                          b = 3, unit = 'mm')), 
                             report$inp_paramaters + 
                               theme(plot.margin = margin(t = 3, 
                                                          l = 10, 
                                                          r = 10, 
                                                          b = 3, unit = 'mm')), 
                             report$summ_table + 
                               theme(plot.margin = margin(t = 3, 
                                                          l = 10, 
                                                          r = 10, 
                                                          b = 3, unit = 'mm')), 
                             report$references + 
                               theme(plot.margin = margin(t = 3, 
                                                          l = 10, 
                                                          r = 10, 
                                                          b = 15, unit = 'mm')), 
                             ncol = 1, 
                             rel_heights = c(1, 4.1, 2.8, 1.7))

    ggsave(plot = report_file, 
           filename = path_to_save, 
           width = 210, 
           height = 290, 
           units = 'mm')
    
  }
  
  ## visualization of the modeling results
  
  plot_mod_results <- function(score_val, risk_scale, fontsize = 4) {
    
    ## makes a simple plot of the modeled risk value in comparison to the prevalence in the study population
    
    ## prevalence of the response in the study cohort
    
    comparator_tbl <- score_tbl %>% 
      filter(!is.na(.data[[risk_scale]])) %>% 
      count_feature('death_acute') %>% 
      filter(death_acute == 1)
    
    comparator_prevalence <- comparator_tbl$percent
    
    comparator_n <- comparator_tbl$total_n
    
    ## risk values predicted from the model
    
    glm_model <- glm_models[[risk_scale]]
    
    inp_data <- tibble(score = score_val) %>% 
      set_names(risk_scale)
    
    pred_res <- calculate_risk(inp_data = inp_data, 
                               inp_model = glm_model, 
                               type = 'response', 
                               se.fit = T)
    
    pred_tbl <- risk_stats_tbl(score_val = score_val, 
                               risk_class = NA, 
                               predict_res = pred_res) %>% 
      mutate(plot_lab = paste(signif(mort_risk, 2) * 100, 
                              ' [95%CI: ', 
                              signif(lower_ci, 2) * 100, 
                              ', ', 
                              signif(upper_ci, 2) * 100, 
                              ']', sep = ''))
    
    ## plot
    
    risk_plot <- pred_tbl %>% 
      ggplot(aes(x = mort_risk * 100, 
                 y = 1, 
                 fill = mort_risk * 100))
    
    ## adding the low (0 to 15%), intermediate (15 - 48%) and high risk (over 48%) ranges to the plot
    ## cutoffs from Kylhammar et al.
    
    risk_plot <- risk_plot + 
      annotate('rect', 
               xmin = -Inf, 
               xmax = 15, 
               ymin = -Inf, 
               ymax = Inf, 
               alpha = 0.15, 
               fill = 'steelblue4') + 
      annotate('rect', 
               xmin = 15, 
               xmax = 48, 
               ymin = -Inf, 
               ymax = Inf, 
               alpha = 0.15, 
               fill = 'white') + 
      annotate('rect', 
               xmin = 48, 
               xmax = Inf, 
               ymin = -Inf, 
               ymax = Inf, 
               alpha = 0.15, 
               fill = 'firebrick4') + 
      annotate('text', 
               label = 'Low risk', 
               size = fontsize, 
               x = 15/2, 
               y = 1.7, 
               color = 'steelblue4', 
               fontface = 'bold') + 
      annotate('text', 
               label = 'Intermediate risk', 
               size = fontsize, 
               x = (48 - 15)/2 + 15, 
               y = 1.7, 
               color = 'gray60', 
               fontface = 'bold') + 
      annotate('text', 
               label = 'High risk', 
               size = fontsize, 
               x = (100 - 48)/2 + 48, 
               y = 1.7, 
               color = 'firebrick4', 
               fontface = 'bold') + 
      geom_vline(xintercept = comparator_prevalence, 
                 linetype = 'dashed') + 
      scale_x_continuous(limits = c(0, 100), 
                         name = 'Predicted 5-year mortality risk, %')
    
    ## all other layers
    
    risk_plot <- risk_plot + 
      geom_point(data = tibble(p = comparator_prevalence), 
                 aes(x = p, 
                     y = 1), 
                 shape = 22, 
                 size = 4, 
                 fill = 'white') + 
      geom_errorbarh(aes(xmin = lower_ci * 100, 
                         xmax = upper_ci * 100), 
                     height = 0.25) + 
      geom_point(shape = 21, 
                 size = 4) + 
      geom_text(aes(label = plot_lab), 
                size = fontsize, 
                hjust = 0.5, 
                vjust = 0, 
                nudge_y = 0.25) + 
      proj_globs$common_theme + 
      theme(axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.line.y = element_blank()) + 
      scale_fill_gradient2(low = 'steelblue4', 
                           mid = 'white', 
                           high = 'firebrick4', 
                           midpoint = comparator_prevalence, 
                           limits = c(0, 100), 
                           name = 'Predicted\nrisk') + 
      scale_y_continuous(limits = c(0.5, 1.8)) + 
      guides(fill = F)
    
    ## adding the total n and mortality in the plot
    
    risk_plot <- risk_plot + 
      labs(tag = paste('\nTotal: n = ', 
                       comparator_n, 
                       ', 5-year mortality: ', 
                       signif(comparator_prevalence, 3),
                       '%', 
                       sep = ''))
    
    return(risk_plot)
    
  }
  
  ## prevalence plotting
  
  plot_score_event <- function(risk_scale, score_val, strata_breaks, strata_labels = NULL, 
                               fill_color = 'coral3', line_color = 'coral4', label_bars = T, 
                               x_lab = proj_globs$score_labs[risk_scale], y_lab = '5-year mortality within score strata', 
                               plot_title = NULL, plot_subtitle = NULL, 
                               plot_tag = NULL, fontsize = 4) {
    
    ## calculates a summary table with the percent positive response in the each score strata
    ## defined by strata_breaks.
    ## Generates and returns the respective bar plot with the strata containing the given score value
    ## (score_val) highlighted
    
    plotting_tbl <- score_tbl %>% 
      filter(!is.na(.data[[risk_scale]])) %>% 
      mutate(score_strata = cut(.data[[risk_scale]], 
                                breaks = strata_breaks, 
                                labels = strata_labels, 
                                include.lowest = T)) %>% 
      dlply(.(score_strata), count_feature, 'death_acute') %>% 
      map(function(x) if(any(x$death_acute == 1)) x else rbind(x, tibble(death_acute = 1, 
                                                                         n = 0, 
                                                                         percent = 0, 
                                                                         total_n = x$total_n))) %>% 
      map2_dfr(., names(.), function(x, y) mutate(x, score_strata = y)) %>% 
      filter(death_acute == 1)
    
    highlight_strata <- tibble(score = score_val) %>% 
      mutate(score_strata = cut(score,
                                breaks = strata_breaks, 
                                labels = strata_labels, 
                                include.lowest = T)) %>% 
      .$score_strata
    
    plotting_tbl <- plotting_tbl %>% 
      mutate(fill_var = ifelse(score_strata == highlight_strata, 
                               'highlight', 
                               'normal'))
    
    if(is.null(plot_tag)) {
      
      plot_tag <- paste('\nTotal: n =', 
                        sum(plotting_tbl$total_n), 
                        ', whole-cohort 5-year mortality: ', 
                        signif(sum(plotting_tbl$n)/sum(plotting_tbl$total_n) * 100, 
                               3), 
                        '%', sep = '')
      
    }
    
    if(label_bars) {
      
      ## should percent response and n number be displayed in the plot?
      
      plotting_tbl <- plotting_tbl %>% 
        mutate(plot_lab = paste(signif(percent, 2), 
                                '%\nn =',
                                total_n))
      
      score_plot <- plotting_tbl %>% 
        ggplot(aes(x = score_strata, 
                   y = percent)) + 
        geom_text(aes(label = plot_lab), 
                  size = fontsize, 
                  hjust = 0.5, 
                  vjust = 0, 
                  nudge_y = 1)
      
    } else {
      
      score_plot <- plotting_tbl %>% 
        ggplot(aes(x = score_strata, 
                   y = percent))
      
    }
    
    score_plot <- score_plot + 
      geom_bar(aes(fill = fill_var), 
               stat = 'identity', 
               color = line_color) + 
      proj_globs$common_theme + 
      labs(x = x_lab, 
           y = y_lab, 
           title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag) + 
      scale_fill_manual(values = c(highlight = 'firebrick4', 
                                   normal = fill_color)) + 
      guides(fill = F)
    
    if(!is.null(strata_labels)) {
      
      score_plot <- score_plot + 
        scale_x_discrete(limits = strata_labels)
      
    }
    
    return(score_plot)
    
  }

# globals ----
  
  proj_globs <- list()
  
  ## graphics
  
  proj_globs$common_text <- element_text(size = 14, 
                                         face = 'plain', 
                                         color = 'black')
  
  proj_globs$common_margin <- margin(t = 5, 
                                     l = 5, 
                                     r = 5, 
                                     unit = 'mm')
  
  proj_globs$common_theme <- theme_classic() + theme(axis.text = proj_globs$common_text, 
                                                     axis.title = proj_globs$common_text, 
                                                     plot.title = element_text(size = 16, 
                                                                               face = 'bold'), 
                                                     plot.subtitle = proj_globs$common_text, 
                                                     plot.tag = element_text(size = 14, 
                                                                             face = 'plain', 
                                                                             color = 'black', 
                                                                             hjust = 0, 
                                                                             vjust = 1), 
                                                     plot.tag.position = 'bottom', 
                                                     legend.text = proj_globs$common_text, 
                                                     legend.title = proj_globs$common_text, 
                                                     strip.text = proj_globs$common_text,
                                                     plot.margin = proj_globs$common_margin)
  
  ## score names, references and plotting colors
  
  proj_globs$score_labs <- c(sign2525 = 'Signature #2525', 
                             enh_french = 'Age/RAA-enhanced FPHR3p', 
                             FRENCH3p = 'FPHR3p', 
                             FRENCH4p = 'FPHR4p', 
                             SPAHR = 'SPAHR', 
                             Compera = 'Compera', 
                             mRASP = 'mRASP', 
                             reveal_lite = 'REVEAL 2.0 Lite')
  
  proj_globs$score_refs <- c(model_10967 = '', 
                             enh_french = '1', 
                             FRENCH3p = '2', 
                             FRENCH4p = '2', 
                             SPAHR = '3', 
                             Compera = '4', 
                             mRASP = '5', 
                             reveal_lite = '6')
  
  proj_globs$score_citations <- c('1' = '1. Sonnweber, T., et al. (2021). Risk assessment in precapillary pulmonary hypertension: a comparative analysis. Respir. Res. 22. doi:10.1186/s12931-021-01624-z.', 
                                  '2' = '2. Boucly, A., et al. (2017). Risk assessment, prognosis and guideline implementation in pulmonary arterial hypertension. Eur. Respir. J. 50, 1700889. doi:10.1183/13993003.00889-2017.', 
                                  '3' = '3. Kylhammar, D., et al. (2018). A comprehensive risk stratification at early follow-up determines prognosis in pulmonary arterial hypertension. Eur. Heart J. 39, 4175-4181. doi:10.1093/eurheartj/ehx257.', 
                                  '4' = '4. Hoeper, M., et al. (2017). Mortality in pulmonary arterial hypertension: Prediction by the 2015 European pulmonary hypertension guidelines risk stratification model. Eur. Respir. J. 50, 1700740. doi:10.1183/13993003.00740-2017.', 
                                  '5' = '5. Xiong, W., et al. (2018). A modified risk score in one-year survival rate assessment of group 1 pulmonary arterial hypertension. BMC Pulm. Med. 18. doi:10.1186/s12890-018-0712-7.', 
                                  '6' = '6. Benza, R. L., et.al. (2021). Development and Validation of an Abridged Version of the REVEAL 2.0 Risk Score Calculator, REVEAL Lite 2, for Use in Patients With Pulmonary Arterial Hypertension. Chest 159, 337-346. doi:10.1016/j.chest.2020.08.2069.')
  
  proj_globs$score_citations_lite <- c('1' = '1. Sonnweber, T., et al. (2021). doi:10.1186/s12931-021-01624-z.', 
                                       '2' = '2. Boucly, A., et al. (2017). doi:10.1183/13993003.00889-2017.', 
                                       '3' = '3. Kylhammar, D., et al. (2018). doi:10.1093/eurheartj/ehx257.', 
                                       '4' = '4. Hoeper, M., et al. (2017). doi:10.1183/13993003.00740-2017.', 
                                       '5' = '5. Xiong, W., et al. (2018). doi:10.1186/s12890-018-0712-7.', 
                                       '6' = '6. Benza, R. L., et.al. (2021). doi:10.1016/j.chest.2020.08.2069.')
  
  proj_globs$method_text <- paste('The 5-year mortality risk was estimated by logistic regression in a cohort of', 
                                  nrow(score_tbl), 
                                  'PAH patients\nrecruited by Innsbruck, Linz and Vienna PAH reference centers as described by Sonnweber et al. (pending).')
  
  ## study prevalence of 5-year deaths
  
  proj_globs$study_prevalence <- score_tbl %>% 
    count_feature('death_acute') %>% 
    filter(death_acute == 1)
  
  
# creating a series of logistic models for 5-year mortality, one for each risk scale ------
  
  
  score_tbl <- score_tbl %>% 
    mutate(reveal_lite = Reveal_lite2_3_cat)
  
  score_tbl$sign2525 <- calculate_sign2525(Gender = score_tbl$Gender, 
                                           age_fc = score_tbl$age_fc, 
                                           cardiac_index = score_tbl$cardiac_index, 
                                           mPAP = score_tbl$mPAP, 
                                           score_only = T)
  
  glm_models <- list(sign2525 = death_acute ~ sign2525, 
                     french_enh = death_acute ~ french_enh, 
                     FRENCH3p = death_acute ~ FRENCH3p, 
                     FRENCH4p = death_acute ~ FRENCH4p, 
                     Compera = death_acute ~ Compera, 
                     SPAHR = death_acute ~ SPAHR, 
                     mRASP = death_acute ~ mRASP, 
                     reveal_lite = death_acute ~ reveal_lite) %>% 
    map(glm, 
        family = 'binomial', 
        data = score_tbl)
  
# END ---