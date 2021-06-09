# This script provides the interface and server functions for an web app
# for calculation of the following risk scores: French3p, French4p, mRASP, Compera, 
# SPAHR, Sonnweber enhanced FPHR3p and #10967 signature
# In the scores and risk categories for the established risk scales are returned
# For the experimental tools (Sonnweber Score and #10967), only the points 
# and study data are displayed.
# Study data to present: modeled risk of overall mortality by Cox (HR) and 5-year mortality by GLM (OR)
# 5-year mortality in the score strata.

# data and tools ----

    library(plyr)
    library(tidyverse)
    library(stringi)
    library(shinyWidgets)

    source('app_tools.R')


# User interface -----

  ui <- fluidPage(
     
     ## some styling
     
     tags$head(
        tags$style(
           "
            .title {
                background: url('banner_pah.png');
                background-repeat: no-repeat; 
                background-size: cover;
                font-size: 55px;
                color: #800000;
                font-family: Bahnschrift, Verdana, Helvetica;
                text-shadow: 1px 1px #e2e2e2; 
                padding-left: 3%;
                padding-top: 1%; 
                padding-bottom: 0.05%
            }
            
            h2 {
               font-size: 30;
               font-weight: bold; 
               font-family: Bahnschrift, Verdana, Helvetica
               }
            
             h3 {
               font-size: 26;
               font-weight: bold; 
               font-family: Bahnschrift, Verdana, Helvetica
               }
            
            h4 {
               font-size: 22; 
               font-family: Bahnschrift, Verdana, Helvetica
            }
            
            .shiny-text-output {
              font_size: 18, 
              font-family: Bahnschrift, Verdana, Helvetica
            }
            
            .shiny-html-output {
              font_size: 18, 
              font-family: Bahnschrift, Verdana, Helvetica
            }
            
            "
        )
     ),
      
     ## Title panel with the logos and names
     
     titlePanel(title =  HTML("<div class = 'title'>
                                 <strong>PAH</strong>
                                 Risk Assessment Tools
                                 <img src = '' width = 50%>
                                 <img src = 'mui_logo.png' width = 7%><br/>
                                 <div/>
                                 <hr style = 'height: 5px'>"), 
                windowTitle = 'PAH Risk Assessment Tools'), 
     
    ## Side panel with user's entries.
    ## The parameters are grouped by demographics (age, sex), performance (SMWD and WHO Fc), 
    ## Laboratory values (NTproBNP, Hb), echo (pericardial effusion, RAA), 
    ## and cardiopulmonary readouts
    
    sidebarLayout(
      
      sidebarPanel(h4('Patient data'), 
                   textInput(inputId = 'pat_ID', 
                             label = 'Patient ID', 
                             value = 'Not specified', 
                             placeholder = 'To be included in the calculation report'), 
                   hr(), 
                   h4('Demographics'), 
                   numericInput(inputId = 'age', 
                                label = 'Age, years', 
                                value = 40, 
                                min = 18, max = 120), 
                   radioButtons(inputId = 'sex', 
                                label = 'Biological sex', 
                                choices = c('Male' = 'male', 
                                            'Female' = 'female')), 
                   hr(), 
                   h4('Performance'), 
                   numericInput(inputId = 'smwd', 
                                label = 'Six Minute Walkig Distance, m', 
                                value = 440, 
                                min = 0, max = 10000), 
                   radioButtons(inputId = 'who_fc', 
                                label = 'WHO Functional Class', 
                                choices = c('I' = 1, 
                                            'II' = 2, 
                                            'III' = 3, 
                                            'IV' = 4)), 
                   hr(), 
                   h4('Ultrasound'), 
                   radioButtons(inputId = 'perricard', 
                                label = 'Paricardial Effusion', 
                                choices = c('Absent or minimal' = 'no', 
                                            'Present' = 'yes')), 
                   numericInput(inputId = 'raa', 
                                label = 'Right Atrial Area, sq. cm', 
                                value = 16, 
                                min = 0, max = 40), 
                   hr(), 
                   h4('Blood Laboratory Parameters'), 
                   numericInput(inputId = 'nt_pro_bnp', 
                                label = 'NT-proBNP, pg/mL', 
                                value = 50, 
                                min = 0, max = 30000), 
                   numericInput(inputId = 'Hb', 
                                label = 'Blood Hemoglobin, g/L', 
                                value = 140, 
                                min = 40, max = 200), 
                   numericInput(inputId = 'SO2_RL', 
                                label = 'Oxygen Saturation, %', 
                                value = 100, 
                                min = 0, max = 100), 
                   hr(), 
                   h4('Cardiopulmonary variables'), 
                   numericInput(inputId = 'SvO2', 
                                label = 'Mixed Venous Oxygen Saturation, %', 
                                value = 100, 
                                min = 0, max = 100),
                   numericInput(inputId = 'card_index', 
                                label = 'Cardiac Index', 
                                value = 2.6, 
                                min = 0, 
                                max = 5), 
                   numericInput(inputId = 'mRAP', 
                                label = 'Mean Right Atrial Pressure, mmHg', 
                                value = 7, 
                                min = 0, max = 40), 
                   numericInput(inputId = 'mPAP', 
                                label = 'Mean Pulmonary Artherial Pressure, mmHg', 
                                value = 65, 
                                min = 0, max = 150), 
                   width = 3),
      
      ## Main panel to hold the dynamic output
      
      mainPanel(
         tabsetPanel(tabPanel('General information', 
                              h3('Welcome to PAH Risk Assessment Tools!'), 
                              hr(), 
                              p('The score calculator and risk modeling app was developed based on the pending 
                                publication by Sonnweber et al. and Sonnweber T et at. (doi:10.1186/s12931−021−01624−z) 
                                as well as seminal papers describing the development of the respecitive scoring tools 
                                (Boucly, A., et al. (2017). doi:10.1183/13993003.00889−2017, 
                                Kylhammar, D., et al. (2018). doi:10.1093/eurheartj/ehx257, 
                                Hoeper, M., et al. (2017). doi:10.1183/13993003.00740−2017 
                                and Xiong, W., et al. (2018). doi:10.1186/s12890−018−0712−7).'), 
                              br(), 
                              p(proj_globs$method_text), 
                              br(), 
                              p('The app developers and publication authors carry no responsibility 
                                for correctness of the risk predictions and the app source code. 
                                This tool may not be used for diagnostic purposes. For scientific use only.'), 
                              hr(), 
                              em('By using the application you accept', 
                                 a('the terms of use and licensing', 
                                   href = 'readme.txt')), 
                              br(),
                              HTML("<div style =  'text-align: right'>
                                  <img src = '' width = 80%>
                                  <p>Powered by </p>
                                  <a href = 'http://www.daas.tirol'>
                                  <img src = 'logo_large.png' width = 60 alt = 'daas.tirol'>
                                  </a>
                                  <img src = '' width = 30>
                                   <img src = 'shiny_logo.png' width = 60></div>")), 
            tabPanel('Score Calulation Summary', 
                     h3('Input Parameters'), 
                     hr(), 
                     tableOutput('inp_par_tbl'), 
                     br(), 
                     h3('Score Calculation and Risk Modeling Summary'), 
                     hr(), 
                     tableOutput('summ_table'), 
                     br(),
                     htmlOutput('methods'),
                     br(), 
                     h3('References'), 
                     hr(), 
                     htmlOutput('score_refs'), 
                     hr(), 
                     downloadButton('downloadReport', 
                                    label = 'Download report')), 
            tabPanel('Signature #10967', 
                     h3('Signature #10967 and study cohort mortality risk'), 
                     br(), 
                     h4('Modeling of the 5-year mortality risk'), 
                     hr(), 
                     plotOutput('model_10967_model', 
                                width = '80%', 
                                height = '200px'), 
                     br(), 
                     p('5-year mortality risk was modeled by logistic regression. 
                     Risk estimate with 95% confidence interval is shown, dashed line represents 
                     the 5-year mortality in the modeling data set. The limits of the 5-year risk 
                       mortality classes were defined as described by Kylhammar, D., et al. (2018). doi:10.1093/eurheartj/ehx257'), 
                     br(), 
                     h4('Real-life 5-year mortality'), 
                     hr(), 
                     plotOutput('model_10967_prev', 
                                width = '65%', 
                                height = '500px'), 
                     br(), 
                     p('5-year mortality in the study cohort stratified by quantiles of the Signature #10967 score.')), 
            tabPanel('Enhanced FPHR3p', 
                     h3('Age/RAA-enhanced FPHR3p and study cohort mortality risk'), 
                     br(), 
                     h4('Modeling of the 5-year mortality risk'), 
                     hr(), 
                     plotOutput('french_enh_model', 
                                width = '65%', 
                                height = '200px'), 
                     br(), 
                     p('5-year mortality risk was modeled by logistic regression. 
                     Risk estimate with 95% confidence interval is shown, dashed line represents 
                     the 5-year mortality in the modeling data set. The limits of the 5-year risk 
                       mortality classes were defined as described by Kylhammar, D., et al. (2018). doi:10.1093/eurheartj/ehx257'), 
                     br(), 
                     h4('Real-life 5-year mortality'), 
                     hr(), 
                     plotOutput('french_enh_prev', 
                                width = '65%', 
                                height = '500px'), 
                     br(), 
                     p('5-year mortality in the study cohort stratified by quantiles of the Age/RAA-enhanced FPHR3p score.')), 
            tabPanel('FPHR3p', 
                     h3('FPHR3p and study cohort mortality risk'), 
                     br(), 
                     h4('Modeling of the 5-year mortality risk'), 
                     hr(), 
                     plotOutput('french3p_model', 
                                width = '65%', 
                                height = '200px'), 
                     br(), 
                     p('5-year mortality risk was modeled by logistic regression. 
                     Risk estimate with 95% confidence interval is shown, dashed line represents 
                     the 5-year mortality in the modeling data set. The limits of the 5-year risk 
                       mortality classes were defined as described by Kylhammar, D., et al. (2018). doi:10.1093/eurheartj/ehx257'), 
                     br(), 
                     h4('Real-life 5-year mortality'), 
                     hr(), 
                     plotOutput('french3p_prev', 
                                width = '65%', 
                                height = '500px'), 
                     br(), 
                     p('5-year mortality in the study cohort stratified by values of the FPHR3p score (number of missed low-risk criteria)')), 
            tabPanel('FPHR4p', 
                     h3('FPHR4p and study cohort mortality risk'), 
                     br(), 
                     h4('Modeling of the 5-year mortality risk'), 
                     hr(), 
                     plotOutput('french4p_model', 
                                width = '65%', 
                                height = '200px'), 
                     br(), 
                     p('5-year mortality risk was modeled by logistic regression. 
                     Risk estimate with 95% confidence interval is shown, dashed line represents 
                     the 5-year mortality in the modeling data set. The limits of the 5-year risk 
                       mortality classes were defined as described by Kylhammar, D., et al. (2018). doi:10.1093/eurheartj/ehx257'), 
                     br(), 
                     h4('Real-life 5-year mortality'), 
                     hr(), 
                     plotOutput('french4p_prev', 
                                width = '65%', 
                                height = '500px'), 
                     br(), 
                     p('5-year mortality in the study cohort stratified by values of the FPHR4p score (number of missed low-risk criteria)')), 
            tabPanel('SPAHR', 
                     h3('SPAHR and study cohort mortality risk'), 
                     br(), 
                     h4('Modeling of the 5-year mortality risk'), 
                     hr(), 
                     plotOutput('spahr_model', 
                                width = '65%', 
                                height = '200px'), 
                     br(), 
                     p('5-year mortality risk was modeled by logistic regression. 
                     Risk estimate with 95% confidence interval is shown, dashed line represents 
                     the 5-year mortality in the modeling data set. The limits of the 5-year risk 
                       mortality classes were defined as described by Kylhammar, D., et al. (2018). doi:10.1093/eurheartj/ehx257'), 
                     br(), 
                     h4('Real-life 5-year mortality'), 
                     hr(), 
                     plotOutput('spahr_prev', 
                                width = '65%', 
                                height = '500px'), 
                     br(), 
                     p('5-year mortality in the study cohort stratified by SPAHR risk classes')), 
            tabPanel('COMPERA', 
                     h3('Compera and study cohort mortality risk'), 
                     br(), 
                     h4('Modeling of the 5-year mortality risk'), 
                     hr(), 
                     plotOutput('compera_model', 
                                width = '65%', 
                                height = '200px'), 
                     br(), 
                     p('5-year mortality risk was modeled by logistic regression. 
                     Risk estimate with 95% confidence interval is shown, dashed line represents 
                     the 5-year mortality in the modeling data set. The limits of the 5-year risk 
                       mortality classes were defined as described by Kylhammar, D., et al. (2018). doi:10.1093/eurheartj/ehx257'), 
                     br(), 
                     h4('Real-life 5-year mortality'), 
                     hr(), 
                     plotOutput('compera_prev', 
                                width = '65%', 
                                height = '500px'), 
                     br(), 
                     p('5-year mortality in the study cohort stratified by Compera risk classes')),
            tabPanel('mRASP', 
                     h3('mRASP and study cohort mortality risk'), 
                     br(), 
                     h4('Modeling of the 5-year mortality risk'), 
                     hr(), 
                     plotOutput('mrasp_model', 
                                width = '65%', 
                                height = '200px'), 
                     br(), 
                     p('5-year mortality risk was modeled by logistic regression. 
                     Risk estimate with 95% confidence interval is shown, dashed line represents 
                     the 5-year mortality in the modeling data set. The limits of the 5-year risk 
                       mortality classes were defined as described by Kylhammar, D., et al. (2018). doi:10.1093/eurheartj/ehx257'), 
                     br(), 
                     h4('Real-life 5-year mortality'), 
                     hr(), 
                     plotOutput('mrasp_prev', 
                                width = '65%', 
                                height = '500px'), 
                     br(), 
                     p('5-year mortality in the study cohort stratified by mRASP risk classes'))
         ), width = 9
      )
      
    )
  )

# Define server logic ----

  server <- function(input, output) {
    
    ## Table with input paramaters
    
    inp_param_tbl <- reactive({
      
      make_input_tbl(PatientID = input$pat_ID, 
                     Gender = input$sex, 
                     Hb = input$Hb, 
                     mPAP = input$mPAP, 
                     SO2_RL = input$SO2_RL, 
                     age_fc = input$age, 
                     RAA = input$raa, 
                     WHOFc = as.numeric(input$who_fc), 
                     SMWD = input$smwd, 
                     NTproBNP = input$nt_pro_bnp, 
                     mRAP = input$mRAP, 
                     cardiac_index = input$card_index, 
                     SvO2 = input$SvO2, 
                     pericardial_effusion = input$perricard)
      
      
    })
    
    output$inp_par_tbl <- renderTable({
      
      inp_param_tbl()
      
    }, hover = T, bordered = T)
    
    ## Score calculation summary tab
    
    score_calc_tbl <- reactive({
      
      ## generates a raw table with score values, risk classes and modeling of 5-year mortality 
      ## in the pooled cohort
      
      make_score_tbl(Gender = input$sex, 
                     Hb = input$Hb, 
                     mPAP = input$mPAP, 
                     SO2_RL = input$SO2_RL, 
                     age_fc = input$age, 
                     RAA = input$raa, 
                     WHOFc = as.numeric(input$who_fc), 
                     SMWD = input$smwd, 
                     NTproBNP = input$nt_pro_bnp, 
                     mRAP = input$mRAP, 
                     cardiac_index = input$card_index, 
                     SvO2 = input$SvO2, 
                     pericardial_effusion = input$perricard)
      
    })
    
    output$summ_table <- renderTable({
      
      ## a simlified summary output

      score_calc_tbl() %>% 
        select(score_label, 
               reference, 
               score, 
               risk_class, 
               five_mort_lab) %>% 
        set_names(c('Risk assessment tool', 
                    'Reference', 
                    'Score', 
                    'Risk class', 
                    'Predicted 5-year mortality risk, %'))
      
    }, hover = T, bordered = T)
    
    ## Text for the methods and references
    
    output$methods <- renderUI({
      
      HTML(proj_globs$method_text %>% 
             stri_replace(fixed = '\n', 
                          replacement = '<br>'))
      
    })
    
    output$score_refs <- renderUI({
      
      HTML(proj_globs$score_citations %>% 
             map(stri_replace, 
                 regex = '^\\d\\.', 
                 replacement = '') %>% 
             map2(.,
                  names(.), 
                  function(x, y) paste('<strong>', y, '.</strong> ', x, sep = '')) %>% 
             paste(collapse = '<br>'))
      
    })
    
    ## model_10967 graphics
    
    model_10967_plots <- reactive({
      
      score_val <- calculate_model_10967(Gender = input$sex, 
                                         Hb = input$Hb, 
                                         mPAP = input$mPAP, 
                                         SO2_RL = input$SO2_RL)
      
      forest_plot <- plot_mod_results(score_val = score_val$score, 
                                      risk_scale = 'model_10967', 
                                      fontsize = 4) + 
        theme(plot.title = element_blank(), 
              plot.subtitle = element_blank())
      
      prev_plot <- plot_score_event(risk_scale = 'model_10967', 
                                    score_val = score_val$score, 
                                    strata_breaks = quantile(score_tbl$model_10967, 
                                                             c(0, 0.25, 0.5, 0.75, 1), 
                                                             na.rm = T)) + 
        scale_x_discrete(limits = levels(cut(score_tbl$model_10967, 
                                             quantile(score_tbl$model_10967, 
                                                      c(0, 0.25, 0.5, 0.75, 1), 
                                                      na.rm = T), 
                                             include.lowest = T))) + 
        expand_limits(y = 65)
      
      return(list(forest = forest_plot, 
                  prevalence = prev_plot))
      
    })
    
    output$model_10967_model <- renderPlot({
      
      model_10967_plots()$forest
      
    })
    
    output$model_10967_prev <- renderPlot({
      
      model_10967_plots()$prevalence
      
    })
    
    ## French enhanced graphics
    
    french_enh_plots <- reactive({
      
      french_val <- calculate_french3p(WHOFc = as.numeric(input$who_fc), 
                                       SMWD = input$smwd, 
                                       NTproBNP = input$nt_pro_bnp)
      
      score_val <- calculate_enh_french(FRENCH3p = french_val$score, 
                                        age_fc = input$age, 
                                        RAA = input$raa)
      
      forest_plot <- plot_mod_results(score_val = score_val$score, 
                                      risk_scale = 'french_enh', 
                                      fontsize = 4) + 
        theme(plot.title = element_blank(), 
              plot.subtitle = element_blank())
      
      prev_plot <- plot_score_event(risk_scale = 'french_enh', 
                                    score_val = score_val$score, 
                                    strata_breaks = c(0, 3, 4, 5, 7)) + 
        expand_limits(y = 65) + 
        scale_x_discrete(limits = levels(cut(score_tbl$french_enh, 
                                             c(0, 3, 4, 5, 7), 
                                             include.lowest = T)))
      
      return(list(forest = forest_plot, 
                  prevalence = prev_plot))
      
    })
    
    output$french_enh_model <- renderPlot({
      
      french_enh_plots()$forest
      
    })
    
    output$french_enh_prev <- renderPlot({
      
      french_enh_plots()$prevalence
      
    })
   
    ## French3p graphics
    
    french3p_plots <- reactive({
      
      score_val <- calculate_french3p(WHOFc = as.numeric(input$who_fc), 
                                      SMWD = input$smwd, 
                                      NTproBNP = input$nt_pro_bnp)
  
      forest_plot <- plot_mod_results(score_val = score_val$score, 
                                      risk_scale = 'FRENCH3p', 
                                      fontsize = 4) + 
        theme(plot.title = element_blank(), 
              plot.subtitle = element_blank())
      
      prev_plot <- plot_score_event(risk_scale = 'FRENCH3p', 
                                    score_val = score_val$score, 
                                    strata_breaks = c(-1, 0.5, 1.5, 2.5, 3.5), 
                                    strata_labels = c('0', '1', '2', '3'), 
                                    x_lab = 'FPHR4p, # missed low-risk criteria') + 
        expand_limits(y = 65)
      
      return(list(forest = forest_plot, 
                  prevalence = prev_plot))
      
    })
    
    output$french3p_model <- renderPlot({
      
      french3p_plots()$forest
      
    })
    
    output$french3p_prev <- renderPlot({
      
      french3p_plots()$prevalence
      
    })
    
    ## French4p graphics
    
    french4p_plots <- reactive({
      
      score_val <- calculate_french4p(WHOFc = as.numeric(input$who_fc), 
                                      SMWD = input$smwd, 
                                      mRAP = input$mRAP, 
                                      cardiac_index = input$card_index)
      
      forest_plot <- plot_mod_results(score_val = score_val$score, 
                                      risk_scale = 'FRENCH4p', 
                                      fontsize = 4) + 
        theme(plot.title = element_blank(), 
              plot.subtitle = element_blank())
      
      prev_plot <- plot_score_event(risk_scale = 'FRENCH4p', 
                                    score_val = score_val$score, 
                                    strata_breaks = c(-1, 0.5, 1.5, 2.5, 3.5, 4.5), 
                                    strata_labels = c('0', '1', '2', '3', '4'), 
                                    x_lab = 'FPHR4p, # missed low-risk criteria') + 
        expand_limits(y = 65)
      
      return(list(forest = forest_plot, 
                  prevalence = prev_plot))
      
    })
    
    output$french4p_model <- renderPlot({
      
      french4p_plots()$forest
      
    })
    
    output$french4p_prev <- renderPlot({
      
      french4p_plots()$prevalence
      
    })
    
    ## SPAHR graphics
    
    spahr_plots <- reactive({
      
      score_val <- calculate_spahr(WHOFc = as.numeric(input$who_fc), 
                                   SMWD = input$smwd, 
                                   NTproBNP = input$nt_pro_bnp, 
                                   mRAP = input$mRAP,
                                   cardiac_index = input$card_index, 
                                   SvO2 = input$SvO2, 
                                   RAA = input$raa, 
                                   pericardial_effusion = input$perricard)
      
      forest_plot <- plot_mod_results(score_val = score_val$score, 
                                      risk_scale = 'SPAHR', 
                                      fontsize = 4) + 
        theme(plot.title = element_blank(), 
              plot.subtitle = element_blank())
      
      prev_plot <- plot_score_event(risk_scale = 'SPAHR', 
                                    score_val = score_val$score, 
                                    strata_breaks = c(0, 1.5, 2.5, 3.5), 
                                    strata_labels = c('low', 'intermediate', 'high'), 
                                    x_lab = 'SPAHR, risk class') + 
        expand_limits(y = 65)
      
      return(list(forest = forest_plot, 
                  prevalence = prev_plot))
      
    })
    
    output$spahr_model <- renderPlot({
      
      spahr_plots()$forest
      
    })
    
    output$spahr_prev <- renderPlot({
      
      spahr_plots()$prevalence
      
    })
    
    ## Compera graphics
    
    compera_plots <- reactive({
      
      score_val <- calculate_compera(WHOFc = as.numeric(input$who_fc), 
                                     SMWD = input$smwd, 
                                     NTproBNP = input$nt_pro_bnp, 
                                     mRAP = input$mRAP,
                                     cardiac_index = input$card_index, 
                                     SvO2 = input$SvO2)
      
      forest_plot <- plot_mod_results(score_val = score_val$score, 
                                      risk_scale = 'Compera', 
                                      fontsize = 4) + 
        theme(plot.title = element_blank(), 
              plot.subtitle = element_blank())
      
      prev_plot <- plot_score_event(risk_scale = 'Compera', 
                                    score_val = score_val$score, 
                                    strata_breaks = c(0, 1.5, 2.5, 3.5), 
                                    strata_labels = c('low', 'intermediate', 'high'), 
                                    x_lab = 'Compera, risk class') + 
        expand_limits(y = 65)
      
      return(list(forest = forest_plot, 
                  prevalence = prev_plot))
      
    })
    
    output$compera_model <- renderPlot({
      
      compera_plots()$forest
      
    })
    
    output$compera_prev <- renderPlot({
      
      compera_plots()$prevalence
      
    })
    
    ## Compera graphics
    
    mrasp_plots <- reactive({
      
      score_val <- calculate_mrasp(WHOFc = as.numeric(input$who_fc), 
                                   SMWD = input$smwd, 
                                   NTproBNP = input$nt_pro_bnp, 
                                   RAA = input$raa)
      
      mRASP_class = cut(score_val$score, 
                        c(-1, 2.5, 5.5, 10), 
                        c(0, 1, 2)) %>% 
        as.character %>% 
        as.numeric
      
      forest_plot <- plot_mod_results(score_val = mRASP_class, 
                                      risk_scale = 'mRASP', 
                                      fontsize = 4) + 
        theme(plot.title = element_blank(), 
              plot.subtitle = element_blank())
      
      prev_plot <- plot_score_event(risk_scale = 'mRASP', 
                                    score_val = mRASP_class, 
                                    strata_breaks = c(-1, 0.5, 1.5, 3.5), 
                                    strata_labels = c('low', 'intermediate', 'high'), 
                                    x_lab = 'mRASP, risk class') + 
        expand_limits(y = 65)
      
      return(list(forest = forest_plot, 
                  prevalence = prev_plot))
      
    })
    
    output$mrasp_model <- renderPlot({
      
      mrasp_plots()$forest
      
    })
    
    output$mrasp_prev <- renderPlot({
      
      mrasp_plots()$prevalence
      
    })
   
   ## download report
   
   output$downloadReport <- downloadHandler(
      
      ## defining the filename
      
      filename = function() {
         
         return(paste('PAH_Tool_report_', Sys.Date(), '.pdf', sep=''))
         
      },
      
      ## calling the saving function
      
      content = function(con) {
         
         render_report(inp_param_tbl = inp_param_tbl(), 
                       score_summ_tbl = score_calc_tbl(), 
                       path_to_save = con)
         
      }
   )
  
  }

# Run the app ----

  shinyApp(ui = ui, server = server)
  
  
  