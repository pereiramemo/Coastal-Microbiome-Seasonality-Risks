# abund_wide = ABUND_s2 
# abund_wide_ref = ABUND_s1
# seqs2range = asvs2range
# semester2range = s22range
# temp_increase = temp_increase

################################################################################
### 1. Function to update abundance matrix based on selected and remove seqs
################################################################################

update_abundance_matrix <- function(abund_wide, abund_wide_prev, seq_id_selected_i, seq_id_removed_i) {

  # identify seqs to select and remove
  abund_wide_selected_i <- abund_wide[,seq_id_selected_i,  drop = F]
  abund_wide_removed_i <- abund_wide[,seq_id_removed_i, drop = F]
  # set all removes seqs to 0 abundance
  abund_wide_removed_i[,] <- 0
  # combine selected and removed abundance matrices
  abund_wide_i <- cbind(abund_wide_selected_i, abund_wide_removed_i)
  # order columns
  abund_wide_i <- abund_wide_i[, colnames(abund_wide_prev)]
  
  return(abund_wide_i)

}

################################################################################
### 2. Function to compute diss by sample
################################################################################

compute_diss_by_sample <- function(abund_wide, abund_wide_ref) {
  
  output <- list()
  for (i in 1:nrow(abund_wide)) {
    
    for (j in 1:nrow(abund_wide_ref)) {
      
      x_ij <- rbind(abund_wide[i,], abund_wide_ref[j,])
      diss_ij <- vegdist(x_ij, method = "bray")
      output[[paste0(i, "_", j)]] <- as.numeric(diss_ij)
    }
  }
  return(unlist(output))
}

################################################################################
### 3. Function to compute similarity with prev iteration
################################################################################

compute_diss_with_prev <- function(abund_wide, abund_wide_prev) {

  output <- list()
  for (j in 1:nrow(abund_wide)) {
    
    x_ij <- rbind(abund_wide_prev[j,], abund_wide[j,])
    dis_j <- vegdist(x_ij, method = "bray") %>% as.numeric()
    sim_j <- 1 - dis_j
    output[as.character(j)] <- sim_j
    
  }
  return(unlist(output))
}
  
################################################################################
### 4. Function to simulate increasing temperature
################################################################################

increasing_temp_simulation <- function(temp_increase, abund_wide, abund_wide_ref, 
                                       seqs2range, semester2range) {

  n <- 0
  output_list <- list()
  for (i in temp_increase) {

    n <- n +1
    
    # Border case: initial community
    if (n == 1) {
      abund_wide_prev <- abund_wide
      diversity_init <- diversity(abund_wide_prev)
      rich_obs_init <- estimateR(abund_wide_prev)["S.obs",]
      rich_chao1_init <- estimateR(abund_wide_prev)["S.chao1",]
      seq_id_ref <-  colnames(abund_wide_ref[,colSums(abund_wide_ref)>0])
      
    } else {
      abund_wide_prev <- abund_wide_i
    }
    
    # Update minimum temperature value
    semester_min_temp_i <- semester2range$min_temp + i
    
    # Select ASVs capable of surviving
    seq_id_selected_i <- seqs2range %>%
                         filter(max_temp >= semester_min_temp_i) %>%
                         pull(seq_id)
    
    # Select ASVs that cannot withstand this temperature 
    seq_id_removed_i <- seqs2range %>%
                        filter(max_temp < semester_min_temp_i) %>%
                        pull(seq_id)
    
    # Update abundance matrix 
    abund_wide_i <- update_abundance_matrix(abund_wide = abund_wide, 
                                            abund_wide_prev = abund_wide_prev,
                                            seq_id_selected_i = seq_id_selected_i, 
                                            seq_id_removed_i = seq_id_removed_i)
    
    # Sanity check: ASVs are the same between current and previous abund tables
    if (all(colnames(abund_wide_i) != colnames(abund_wide_prev))) {
      # exit function
      return("ASVs are not the same")
    }
    
    # Compute diversity and richness
    diversity_i <- diversity(abund_wide_i)
    rich_obs_i <- estimateR(abund_wide_i)["S.obs",]
    rich_chao1_i <- estimateR(abund_wide_i)["S.chao1",]
    
    # Compute proportions of shared Seq IDs
    seq_id_selected_i_gt_zero_index <- colnames(abund_wide[,colSums(abund_wide) >0]) %in% seq_id_selected_i
    seq_id_selected_i_gt_zero <-  seq_id_selected_i[seq_id_selected_i_gt_zero_index]
    
    seq_id_removed_i_gt_zero_index <- colnames(abund_wide[,colSums(abund_wide) >0]) %in% seq_id_removed_i
    seq_id_removed_i_gt_zero <-  seq_id_removed_i[seq_id_removed_i_gt_zero_index]
    
    selected_shared_seq_id <- 100*sum(seq_id_ref %in% seq_id_selected_i_gt_zero)/length(seq_id_ref)
    removed_shared_seq_id <- 100*sum(seq_id_ref %in% seq_id_removed_i_gt_zero)/length(seq_id_ref)
    # selected_shared_seq_id <- 100*sum(seq_id_selected_i_gt_zero %in% seq_id_ref)/length(seq_id_selected_i_gt_zero)
    # removed_shared_seq_id <- 100*sum(seq_id_removed_i_gt_zero %in% seq_id_ref)/length(seq_id_removed_i_gt_zero)
    
    
    # Compute dissimilarity with i-1
    diss_with_prev <- compute_diss_with_prev(abund_wide = (abund_wide_i > 0) +0, 
                                             abund_wide_prev = (abund_wide_prev > 0) + 0) 
    
    # Compute dissimilarity with reference community
    diss_by_sample <- compute_diss_by_sample(abund_wide = (abund_wide_i > 0) +0, 
                                             abund_wide_ref = (abund_wide_ref + 0) + 0)  
    
    # Save outputs
    output_list[[as.character(i)]][["temp"]] <- i
    
    output_list[[as.character(i)]][["div"]] <- list(diversity = 100*diversity_i/diversity_init,
                                                    rich_obs = 100*rich_obs_i/rich_obs_init,
                                                    rich_chao1 = 100*rich_chao1_i/rich_chao1_init)
    
    output_list[[as.character(i)]][["shared_seq_id"]] <- list(selected = selected_shared_seq_id,
                                                              removed = removed_shared_seq_id)
    
    output_list[[as.character(i)]][["seq_id_selected"]] <- seq_id_selected_i
    output_list[[as.character(i)]][["seq_id_removed"]] <- seq_id_removed_i
    output_list[[as.character(i)]][["diss_with_prev"]] <- diss_with_prev 
    output_list[[as.character(i)]][["diss_with_ref"]] <- diss_by_sample
    
  }

  return(output_list)
  
}
