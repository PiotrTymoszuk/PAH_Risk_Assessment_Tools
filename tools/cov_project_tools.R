# This script provides tools for the COV-Questionaire project


# libraries and data ----

   library(tidyr)
   library(tibble)

# globals ----

  proj_globs <- list()

  proj_globs$present_scale <- c(yes = 'firebrick4', 
                                no = 'steelblue4')
  
  proj_globs$reg_scale <- c(positive = 'firebrick3', 
                            negative = 'steelblue3', 
                            ns = 'gray60')
  
  proj_globs$common_text <- element_text(size = 8, 
                                         face = 'plain', 
                                         color = 'black')
  
  proj_globs$common_margin <- margin(t = 5, 
                                     l = 5, 
                                     r = 5, 
                                     unit = 'mm')
  
  proj_globs$common_theme <- theme_classic() + theme(axis.text = proj_globs$common_text, 
                                                     axis.title = proj_globs$common_text, 
                                                     plot.title = element_text(size = 10, 
                                                                               face = 'bold'), 
                                                     plot.subtitle = proj_globs$common_text, 
                                                     plot.tag = element_text(size = 8, 
                                                                             face = 'plain', 
                                                                             color = 'black', 
                                                                             hjust = 0, 
                                                                             vjust = 1), 
                                                     plot.tag.position = 'bottom', 
                                                     legend.text = proj_globs$common_text, 
                                                     legend.title = proj_globs$common_text, 
                                                     strip.text = proj_globs$common_text,
                                                     plot.margin = proj_globs$common_margin)

# variable summing and recoding functions ----
  
   yn_variable <- function(inp_vec) {
     
     ## recodes 0 to no, 1 to yes
     
     out_vec <- car::recode(inp_vec, 
                            "'0' = 'no'; 
                            '1' = 'yes'")
     
     return(out_vec)
     
  }

   sum_variables <- function(inp_tbl, vars_to_sum) {
      
      ## calculates a sum of the given variables
      ## coded as yes/no
      
      var_list <- vars_to_sum %>% 
        map(function(x) car::recode(inp_tbl[[x]], "'no' = 0; 'yes' = 1")) %>% 
        #map(function(x) ifelse(inp_tbl[[x]] == 'no' | is.na(inp_tbl[[x]]), 0, 1)) %>% 
        reduce(function(x, y) x + y)
      
      return(var_list)
      
   }
 
   binarize_variables <- function(inp_tbl, vars_to_recode) {
     
     ## recodes yes/no to 1/0
     
     tbl_list <- vars_to_recode %>% 
       map(function(x) tibble(ID = inp_tbl[['ID']], 
                              sympt_recoded = as.numeric(car::recode(inp_tbl[[x]], 
                                                                     "'no' = 0; 
                                                                      'yes' = 1")))) %>% 
       set_names(vars_to_recode)
     
     rec_tbl <- tbl_list %>% 
       map2(., names(.), 
            function(x, y) set_names(x, c('ID', y))) %>% 
       reduce(left_join, 
              by = 'ID')
     
     return(rec_tbl)
     
   }
   
   translate_var <- function(var_vector, dictionary = globals$risk_labels) {
      
      ## translates the variable svector to the labels
      
      out_vector <- dictionary[var_vector]
      
      return(out_vector)
      
   }
   
# feature percent calculation, plotting and modeling -----
   
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
   
   count_feature_lst <- function(inp_tbl, var_to_count_vec, remove_na = T, positive_only = T) {
      
      ## wrapper for variable vectors
      
      count_tbl <- var_to_count_vec %>% 
         map(count_feature, 
             inp_tbl = inp_tbl, 
             remove_na = remove_na) %>% 
         map(set_names, 
             c('feature_strata', 
               'n', 
               'percent', 
               'total_n')) %>% 
         set_names(var_to_count_vec) %>% 
         map2_dfr(., names(.), 
                  function(x, y) mutate(x, feature = y))
      
      if(positive_only) {
         
         count_tbl <- count_tbl %>% 
            filter(feature_strata == 'yes')
         
      }
      
      return(count_tbl)
      
   }
   
   pie_feature <- function(perc_tbl, plot_title = names(perc_tbl)[1], plot_subtitle = NULL, 
                           order_tbl = NULL, pie = T, repel = T, fill_colors = NULL) {
     
     ## plots the distribution of the given feature: takes an output of the count_feature function
     ## a tibble with order for plotting may be provided (names: plot_order, feature_strata)
     
     plotting_tbl <- perc_tbl %>% 
       set_names(c('feature_strata', 
                   'n', 
                   'percent', 
                   'total_n'))
     
     plot_tag <- paste('n =', plotting_tbl[['total_n']][1])
     
     if(!is.null(order_tbl)) {
       
       plotting_tbl <- plotting_tbl %>% 
         left_join(., 
                   order_tbl, 
                   by = 'feature_strata') %>% 
         arrange(desc(plot_order))
       
     } else {
       
       plotting_tbl <- plotting_tbl %>% 
         arrange(desc(feature_strata))
       
     }
     
     plotting_tbl <- plotting_tbl %>% 
       mutate(plot_lab = paste(signif(percent, 2), '%', sep = ''), 
              plot_y = cumsum(percent) - 0.5 * percent)
     
     if(!is.null(order_tbl)) {
       
       perc_plot <- plotting_tbl %>% 
         ggplot(aes(x = '', 
                    y = percent, plot_order, 
                    fill = reorder(feature_strata, plot_order)))
       
     } else {
       
       perc_plot <- plotting_tbl %>% 
         ggplot(aes(x = '', 
                    y = percent, 
                    fill = feature_strata))
       
     }
     
     perc_plot <- perc_plot + 
       geom_bar(stat = 'identity', 
                position = 'stack', 
                color = 'black') + 
       proj_globs$common_theme + 
       theme(axis.title.x = element_blank(), 
             axis.text.x = element_blank(), 
             axis.ticks.x = element_blank()) + 
       labs(title = plot_title, 
            subtitle = plot_subtitle, 
            tag = plot_tag, 
            y = '% of complete answers')
     
     if(repel) {
       
       perc_plot <- perc_plot + 
         geom_label_repel(aes(label = plot_lab, 
                              y = plot_y), 
                          size = 2.75, 
                          box.padding = 0.1, 
                          label.padding = 0.1, 
                          show.legend = F)
       
     } else {
       
       perc_plot <- perc_plot + 
         geom_label(aes(label = plot_lab, 
                       y = plot_y), 
                   size = 2.75, 
                   label.padding = 0.1, 
                   show.legend = F)
       
     }
     
     if(pie) {
       
       perc_plot <- perc_plot + 
         coord_polar(theta = 'y') + 
         theme(axis.line = element_blank(), 
               axis.title = element_blank(), 
               axis.text = element_blank(), 
               axis.ticks = element_blank())
       
     }
     
     if(!is.null(fill_colors)) {
       
       perc_plot <- perc_plot + 
         scale_fill_manual(values = fill_colors, 
                           name = '')
       
     }
     
     return(perc_plot)
     
     
   }
   


   
   plot_score_event <- function(score_tbl, strata_breaks, strata_labels = NULL, 
                                fill_color = 'coral3', label_bars = T, 
                                x_lab = 'Score', y_lab = '% cases within score strata', 
                                plot_title = NULL, plot_subtitle = NULL, 
                                plot_tag = NULL) {
      
      ## returns a summary table with the percent positive response in the each score strata
      ## and the respective bar plot
      
      plotting_tbl <- score_tbl %>% 
         mutate(score_strata = cut(score, 
                                   breaks = strata_breaks, 
                                   labels = strata_labels)) %>% 
         dlply(.(score_strata), count_feature, 'response') %>% 
         map2_dfr(., names(.), function(x, y) mutate(x, score_strata = y)) %>% 
         filter(response == 1)
      
      if(label_bars) {
         
         plotting_tbl <- plotting_tbl %>% 
            mutate(plot_lab = paste(signif(percent, 2), 
                                    '%\nn =',
                                    total_n))
         
         score_plot <- plotting_tbl %>% 
            ggplot(aes(x = score_strata, 
                       y = percent)) + 
            geom_text(aes(label = plot_lab), 
                      size = 2.75, 
                      hjust = 0.5, 
                      vjust = 0, 
                      nudge_y = 1)
         
      } else {
         
         score_plot <- plotting_tbl %>% 
            ggplot(aes(x = score_strata, 
                       y = percent))
         
      }
      
      score_plot <- score_plot + 
         geom_bar(stat = 'identity', 
                  color = 'black', 
                  fill = fill_color) + 
         proj_globs$common_theme + 
         labs(x = x_lab, 
              y = y_lab, 
              title = plot_title, 
              subtitle = plot_subtitle, 
              tag = plot_tag)
      
      if(!is.null(strata_labels)) {
         
         score_plot <- score_plot + 
            scale_x_discrete(limits = strata_labels)
         
      }
      
      return(score_plot)
      
   }
   
# END ----