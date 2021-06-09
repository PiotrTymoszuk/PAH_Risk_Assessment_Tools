# This script contains common tools used throughout the MS monocyte manuscript like plotting functions and modeling toolbox

# libraries ----

  library(survival)
  library(survminer)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)

    
# survival toolbox ----
    
    create_surv <- function(inp_table, time_variable, event_variable) {
      
      ## creates a survival object
      
      return(Surv(inp_table[[time_variable]], inp_table[[event_variable]]))
      
    }
    
    model_coxph <- function(inp_table, surv_object, dep_variable, return_simple = T) {
      
      ## models time to the event for as a function of a dependent variable
      ## with a one paremeter Cox model
      
      int_surv_object <- surv_object
      
      mod_formula <- paste('int_surv_object', 
                           dep_variable, 
                           sep = '~') %>% 
        as.formula
      
      #surv_table <- inp_table %>% 
        #mutate(dep_variable = .[[dep_variable]])
      
      #cox_object <- coxph(surv_object ~ dep_variable, data = surv_table)
      
      cox_object <- coxph(formula = mod_formula, 
                          data = inp_table)
      
      if(return_simple) {
        
        return(cox_object)
        
      } else {
        
        return(list(surv_table = surv_table, 
                    surv_object = surv_object, 
                    cox_object = cox_object, 
                    summary = summary(cox_object))) 
        
      }
      
    }
    
    model_km_cutpoint <- function(inp_table, surv_object, dep_variable, cutpoint = NULL, ...) {
      
      ## performs KM modeling for the given table, survial object, variable of interest and a distance cutpoint
      
      modeling_table <- inp_table %>% 
        mutate(dep_variable = .[[dep_variable]])
      
      min_var <- min(modeling_table$dep_variable, na.rm = T) - 1
      
      max_var <- max(modeling_table$dep_variable, na.rm = T)
      
      if(!is.null(cutpoint)) {
        
        modeling_table <- modeling_table %>% 
          mutate(var_strata = cut(dep_variable, c(min_var, cutpoint, max_var)))
        
      } else {
        
        modeling_table <- modeling_table %>% 
          mutate(var_strata = dep_variable)
        
      }
      
      km_fit <- survfit(surv_object ~ var_strata, data = modeling_table)
      
      test <- survdiff(surv_object ~ var_strata, data = modeling_table, ...)
      
      stats <- list(chi = test$chisq, df = length(test$n) - 1, n = test$n)
      
      p_value <- 1 - pchisq(test$chisq, length(test$n) - 1)
      
      return(list(test_table = modeling_table, 
                  km_fit = km_fit, 
                  test = test, 
                  cutpoint = cutpoint, 
                  stats = stats, 
                  p_value = p_value))
      
    }
    
    optimize_km <- function(inp_table, surv_object, dep_variable, step = 0.1, min_n = 0, ...) {
      
      ## finds am optimal cutpoint for the KM modeling of the survival in a function of the staratified distance
      ## the optional min_n argument defines the minimal group size for a dep_variable strata
      
      end_cutoff <- max(inp_table[[dep_variable]], na.rm = T)
      
      current_cutoff <- min(inp_table[[dep_variable]], na.rm = T)
      
      result_frame <- data.frame()
      
      while(current_cutoff < end_cutoff) {
        
        current_km_model <- try(model_km_cutpoint(inp_table, surv_object, 
                                                  dep_variable = dep_variable, 
                                                  cutpoint = current_cutoff, ...), silent = T)
        
        if(class(current_km_model) == 'try-error') {
          
          current_cutoff <- current_cutoff + step
          
        } else if(min(current_km_model$km_fit$n) < min_n) {
          
          current_cutoff <- current_cutoff + step
          
        } else {
          
          p_value <- current_km_model$p_value
          
          result_frame <- rbind(result_frame, c(current_cutoff, current_km_model$p_value, current_km_model$km_fit$n))
          
          current_cutoff <- current_cutoff + step
          
        }
        
      }
      
      result_frame <- result_frame %>% 
        set_names(c('cutoff', 'p_value', 'n1', 'n2')) %>% 
        as_tibble
      
      return(list(result_table = result_frame, 
                  optimal_cutoff = result_frame$cutoff[result_frame$p_value == min(result_frame$p_value)]))
      
    }
    
    km_optimization_plot <- function(optimize_km_results) {
      
      ## takes results of the optimize_km results and returns a basic diagnostic plot
      ## to assess the quality of KM optimization: group sizes (n) and p value as a function
      ## of cutoff
      
      opt_plot <- tcga_km_optimization_plot <- optimize_km_results$result_table %>% 
        mutate(neg_log10_p_value = -log10(p_value)) %>% 
        gather(key = 'parameter', value = 'par_value', p_value, n1, n2, neg_log10_p_value) %>% 
        ggplot(aes(x = cutoff, y = par_value)) + 
        geom_line() + 
        facet_grid(parameter ~ ., 
                   scales = 'free', 
                   labeller = labeller(.cols = as_labeller(c(n1 = 'n1',
                                                             n2 = 'n2', 
                                                             p_value = 'p value', 
                                                             neg_log10_p_value = '-log10 p value')))) + 
        theme(axis.title.y = element_blank())
      
      return(opt_plot)
      
    }
    
    km_summary <- function(km_model_object) {
      
      ## returns the summary of the km_model_object
      
      output <- km_model_object$stats[c('chi', 'df')] %>% 
        reduce(cbind)
      
      output <- km_model_object$stats$n %>% 
        reduce(cbind) %>% 
        cbind(output, .)
      
      output <- cbind(output, km_model_object$p_value)
      
      if(!is.null(km_model_object$cutpoint)) {
        
        output <- cbind(output, km_model_object$cutpoint) %>% 
          as_tibble
        
        names(output) <- c('ChiSq', 'Df', paste('N', 1:length(km_model_object$stats$n), sep = ''), 'P_value', 'Cutoff')
        
      } else {
        
        output <- output %>% 
          as_tibble
        
        names(output) <- c('ChiSq', 'Df', paste('N', 1:length(km_model_object$stats$n), sep = ''), 'P_value')
        
      }
      
      
      
      return(output)
      
    }
    
    plot_km <- function(km_model_object, surv_object, ...) {
      
      # plots the km model of interest with the ggsurvplot function
      # ... are additional arguments for the ggsurvplot function
      
      fit <- surv_fit(surv_object ~ var_strata, data = km_model_object$test_table)
      
      output <- ggsurvplot(fit, ggtheme = common_theme, ...)
      
      return(output)
      
    }
    
    save_km_plot <- function(km_plot, file_name, width, height) {
      
      ## saves a km plot on the disc together with the risk table
      
      cairo_pdf(file_name, width = cm_to_inch(width), height = cm_to_inch(height))
      
      print(km_plot, newpage = FALSE)
      
      dev.off()
      
    }

# varia ----
    
    cm_to_inch <- function(inp_cm){
      
      return(0.393701 * inp_cm)
      
    }
    
# END ---
    

    
    