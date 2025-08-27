################################################################################
# Function Rare Indval. 
################################################################################

# This function rarefies the abundance table using different seed values, and
# runs the IndVal analysis for each run. 
# The function stores the selected ASVs or genes for each iteration to determine
# those cases that are more frquently selected. 

################################################################################
# 1. Load libraries
################################################################################

library(tidyverse)
library(vegan)
library(indicspecies)

################################################################################
# 2. Define variable for dev
################################################################################

# Dev only
#########
# abundance_table <- ASV_ABUND_TBL_filt_wide
# sample_size <- rowSums(abundance_table) %>% min()
# semester <- SAMPLES2COMMUNITY %>%
#   dplyr::select(Date_formatted, Community) %>%
#   column_to_rownames("Date_formatted")
# 
# semester_ord <- semester[rownames(abundance_table),, drop = F]
# all(rownames(semester_ord) == rownames(abundance_table))
# vector_group <- semester_ord$Community
# nperm <- 999
# p_value <- 1e-2
# iterations <- 100

# iterations = 100
# abundance_table = abundance_table
# sample_size = sample_size
# vector_group = vector_group
# nperm = 9999
# p_value = 1e-2

#########

################################################################################
# 3. Define function
################################################################################

rare_indval <- function(iterations = 100, abundance_table, sample_size,
                        vector_group, nperm = 9999, p_value = 1e-3) {

  
  output_list <- foreach(i = 1:iterations) %dopar%   {
    
    # Rarefy the abundance table
    set.seed(i)
    rarefied_data <- rrarefy(abundance_table, sample_size)
    
    # Run IndVal analysis
    set.seed(i)
    indval_output <- multipatt(rarefied_data, 
                        vector_group, 
                        control = how(nperm = nperm), 
                        func = "IndVal.g", 
                        max.order = 1)
    
    # Extract indval results
    indval_output_df <- indval_output$sign
    
    indval_output_df_sig_s1 <- indval_output_df %>%
                               filter(p.value <=  p_value & s.S1 == 1) %>%
                               arrange(desc(stat)) %>%
                               mutate(semester = "S1", iteration = i) %>%
                               rownames_to_column("asv_id")
    
    indval_output_df_sig_s2 <- indval_output_df %>%
                               filter(p.value <= p_value & s.S2 == 1) %>%
                               arrange(desc(stat)) %>%
                               mutate(semester = "S2", iteration  = i) %>%
                               rownames_to_column("asv_id")
    
    # Save the results to a list or data frame
    output_list_i <- list()
    output_list_i[["S1"]] <- indval_output_df_sig_s1
    output_list_i[["S2"]] <- indval_output_df_sig_s2
    
    return(output_list_i)
    
  }
  
  # Extract output_list to a data frame
  output_df <- bind_rows(lapply(output_list, function(x) {
    bind_rows(x[["S1"]], x[["S2"]])
  }))
  
  # Count the number of times each ASV or gene was selected
  selected_counts_df <- output_df %>%
    group_by(asv_id, semester) %>%
    summarise(consistency = 100*n()/iterations,
              stat_mean = mean(stat),
              stat_sd = sd(stat)) %>%
    arrange(desc(consistency))
  
  return(selected_counts_df)
  
}

