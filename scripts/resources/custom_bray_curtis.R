###############################################################################
### Custom Bray Curtis
###############################################################################

# It computes the Bray Curtis dissimilarity between samples, including
# the B (losses) and C (gains) components.
# Additionally, it keeps track of the contribution of all OTUs 
# to each component
# The output consists of three matrices and two data frames:
# 1. The standard Bray Curtis dissimilarity matrix
# 2. The Bray Curtis B (losses) component dissimilarity matrix
# 3. The Bray Curtis C (gains) component dissimilarity matrix
# 4. The species contributions to the B (losses) component
# 5. The species contributions to the c (gains) component

###############################################################################
#### custom_bray_curtis_subfun sub function: rows iteration
###############################################################################

custom_bray_curtis_subfun <- function(i, mtx, comp_log, unscale) {
  
  dimnames_list <- list(row = rownames(mtx)[i], col = rownames(mtx))
  D_diff_vector <- matrix(nrow = 1, ncol = dim(mtx)[1], dimnames = dimnames_list, NA) 
  B_loss_vector <- matrix(nrow = 1, ncol = dim(mtx)[1], dimnames = dimnames_list, NA) 
  C_gain_vector <- matrix(nrow = 1, ncol = dim(mtx)[1], dimnames = dimnames_list, NA) 
  output_list <- list()
  
  if (comp_log == T) {
    output_list_comp <- list()
    output_list_comp[["B_loss_sp"]] <- list()
    output_list_comp[["C_gain_sp"]] <- list()
  }
  
  # compute only n(n-1)/2 values  
  if (i < dim(mtx)[1]) {
    
    for (j in (i+1):(dim(mtx)[1])) {
      
      # reduce columns to iterate
      x <- mtx[i,] > 0
      y <- mtx[j,] > 0
      z <- (y & x)
      
      # compute A component
      A <- 0
      for (k in colnames(mtx[,z, drop = F])) {
        Aj <- min(mtx[i,k], mtx[j,k])
        A <- A + Aj
      }
      
      # record B comp
      if (comp_log == T) {
        B_loss_comp <- setNames(rep(NA, dim(mtx)[2]), colnames(mtx))
      }
      
      # compute B component
      B <- 0
      for (k in colnames(mtx[,x, drop = F])) {
        Bj_test <- mtx[i,k] - mtx[j,k]
        if (Bj_test > 0) {
          
          Bj <- Bj_test
          B <- B + Bj
          
          if (comp_log == T) {
            B_loss_comp[k] <- Bj
          }
          
        }
      }
      
      # record C comp
      if (comp_log == T) {
        C_gain_comp <- setNames(rep(NA, dim(mtx)[2]), colnames(mtx)) 
      }
      
      # compute C component
      C <- 0
      for (k in colnames(mtx[,y, drop = F])) {
        Cj_test <- mtx[j,k] - mtx[i,k]
        if (Cj_test > 0) {
          
          Cj <- Cj_test
          C <- C + Cj
          
          if (comp_log == T) {
            C_gain_comp[k] <- Cj
          }
          
        }
      }
      
      # load values in vector
      if (unscale == F) { 
        D_diff_ij <- (B + C)/(2*A + B + C)
        B_loss_ij <- (B)/(2*A + B + C)
        C_gain_ij <- (C)/(2*A + B + C)
        D_diff_vector[,j] <- D_diff_ij
        B_loss_vector[,j] <- B_loss_ij
        C_gain_vector[,j] <- C_gain_ij
      }
      
      if (unscale == T) { 
        D_diff_ij <- (B + C)
        B_loss_ij <- (B)
        C_gain_ij <- (C)
        D_diff_vector[,j] <- D_diff_ij
        B_loss_vector[,j] <- B_loss_ij
        C_gain_vector[,j] <- C_gain_ij
      }
        
      # load values in list
      if (comp_log == T) {
        comparison <- paste(rownames(mtx)[i], rownames(mtx)[j], sep = " vs ")
        output_list_comp[["B_loss_sp"]][[comparison]] <- B_loss_comp/(2*A + B + C)
        output_list_comp[["C_gain_sp"]][[comparison]] <- C_gain_comp/(2*A + B + C)
      }    
      
    } # end for loop in i 
    
    # add final value to match squared matrix
  } else { 
    
    j <- i
    
    if (comp_log == T) {
      B_loss_comp <- setNames(rep(NA, dim(mtx)[2]), colnames(mtx))  
      C_gain_comp <- setNames(rep(NA, dim(mtx)[2]), colnames(mtx))
      comparison <- paste(rownames(mtx)[i], rownames(mtx)[j], sep = " vs ")
      output_list_comp[["B_loss_sp"]][[comparison]] <- B_loss_comp
      output_list_comp[["C_gain_sp"]][[comparison]] <- C_gain_comp
    }  
    
  }# close if statement
  
  output_list[["mtx"]][[1]] <- D_diff_vector
  output_list[["mtx"]][[2]] <- B_loss_vector
  output_list[["mtx"]][[3]] <- C_gain_vector
  
  if (comp_log == T) {
    output_list[["comp"]][[1]] <- do.call(cbind.data.frame, output_list_comp[["B_loss_sp"]])
    output_list[["comp"]][[2]] <- do.call(cbind.data.frame, output_list_comp[["C_gain_sp"]])
  }
  
  return(output_list)
  
} # function

###############################################################################
#### output combine function: matrix rows rbind and data frame cbind
###############################################################################

combine_fun <- function(X1, X2) {
  
  output_list <- list()
  
  D_diff_mtx <- rbind(X1[["mtx"]][[1]], X2[["mtx"]][[1]])
  B_loss_mtx <- rbind(X1[["mtx"]][[2]], X2[["mtx"]][[2]])
  C_gain_mtx <- rbind(X1[["mtx"]][[3]], X2[["mtx"]][[3]])
  
  output_list[["mtx"]][["D_diff"]] <- D_diff_mtx
  output_list[["mtx"]][["B_loss"]] <- B_loss_mtx
  output_list[["mtx"]][["C_gain"]] <- C_gain_mtx
  
  if (("comp" %in% names(X1)) ==  T) {
    B_loss_comp_df <- cbind(X1[["comp"]][[1]], X2[["comp"]][[1]])
    C_loss_comp_df <- cbind(X1[["comp"]][[2]], X2[["comp"]][[2]])
    output_list[["comp"]][["B_loss_sp"]] <- B_loss_comp_df
    output_list[["comp"]][["C_gain_sp"]] <- C_loss_comp_df
  }
    
  return(output_list)
  
}

###############################################################################
#### custom_bray_curtis function: columns iterations
###############################################################################

custom_bray_curtis <- function(mtx, comp_log = F, unscale = F, ncores = 16) {
  
  library(doParallel)  
  
  doParallel::registerDoParallel(ncores)
  
  n <- (dim(mtx)[1])
  X <- foreach(i = 1:n, 
               .combine = "combine_fun", 
               .export="custom_bray_curtis_subfun") %dopar% custom_bray_curtis_subfun(i, mtx, comp_log, unscale) 
  
  return(X)
}  

###############################################################################
### tests function
###############################################################################
# 
# mtx <- abund_asvs

# mtx <- abund_opus_ord[1:2,]

# ABUND_total_dist_list <- custom_bray_curtis(mtx, comp_log = T, unscale = F)
# 
# ABUND_total_D_diff_mtx <- ABUND_total_dist_list[["mtx"]][["D_diff"]]
# ABUND_total_B_loss_mtx <- ABUND_total_dist_list[["mtx"]][["B_loss"]]
# ABUND_total_C_gain_mtx <- ABUND_total_dist_list[["mtx"]][["C_gain"]]
# 
# X_ref_vegan <- as.matrix(vegdist(mtx, method = "bray"))
# plot(X_ref_vegan[upper.tri(X_ref_vegan, diag = F)],
#      ABUND_total_D_diff_mtx[upper.tri(ABUND_total_D_diff_mtx, diag = F)])
# abline(a = 0, b = 1)
# mean(X_ref_vegan[upper.tri(X_ref_vegan, diag = F)]- ABUND_total_D_diff_mtx[upper.tri(ABUND_total_D_diff_mtx, diag = F)])
# 
# ABUND_total_D_diff_mtx_test <- ABUND_total_B_loss_mtx + ABUND_total_C_gain_mtx
# plot(ABUND_total_D_diff_mtx_test[upper.tri(ABUND_total_D_diff_mtx_test, diag = F)],
#      X_ref_vegan[upper.tri(X_ref_vegan, diag = F)])
# mean(ABUND_total_D_diff_mtx_test[upper.tri(ABUND_total_D_diff_mtx_test, diag = F)] -
#        X_ref_vegan[upper.tri(X_ref_vegan, diag = F)])
# 
# ###############################################################################
# ### test subfunction
# ###############################################################################
# 
# i <- 10
# X_test <- custom_bray_curtis_subfun(i, mtx, comp_log = F, unscale = F)
# ABUND_total_D_diff_mtx_test_row <- X_test[["mtx"]][[1]]
# sum(ABUND_total_D_diff_mtx_test_row[1,] == ABUND_total_D_diff_mtx[i,], na.rm = T)
# 
# ABUND_total_B_loss_mtx_test_row <- X_test[["mtx"]][[2]]
# sum(ABUND_total_B_loss_mtx_test_row[1,] == ABUND_total_B_loss_mtx[i,], na.rm = T)
# 
# ABUND_total_C_gain_mtx_test_row <- X_test[["mtx"]][[3]]
# sum(ABUND_total_C_gain_mtx_test_row[1,] == ABUND_total_C_gain_mtx[i,], na.rm = T)
# 
# dim(X_test[["mtx"]][[1]])
# sum(!is.na(X_test[["mtx"]][[1]]))
# 
# ###############################################################################
# ### test comp log
# ###############################################################################
#  
# ABUND_total_B_loss_sp_df <- ABUND_total_dist_list[["comp"]][["B_loss_sp"]]
# ABUND_total_C_gain_sp_df <- ABUND_total_dist_list[["comp"]][["C_gain_sp"]]
# 
# sum(ABUND_total_B_loss_sp_df[,"2018-03-14 vs 2018-03-20"], na.rm = T)
# sum(ABUND_total_C_gain_sp_df[,"2018-03-14 vs 2018-03-20"], na.rm = T)
# ABUND_total_B_loss_mtx["2018-03-14", "2018-03-20"]
# ABUND_total_C_gain_mtx["2018-03-14", "2018-03-20"]
# 
# sum(ABUND_total_B_loss_sp_df[,"2018-03-14 vs 2018-03-20"], na.rm = T) +
# sum(ABUND_total_C_gain_sp_df[,"2018-03-14 vs 2018-03-20"], na.rm = T)
# ABUND_total_D_diff_mtx["2018-03-14", "2018-03-20"]
# 
# ###############################################################################
# ### Sanity check
# ###############################################################################
# 
# date_1 <- rownames(mtx)[10]
# date_2 <- rownames(mtx)[24]
# 
# tbi_test <- adespatial::TBI(mat1 = mtx[date_1,,drop = F],
#                             mat2 = mtx[date_2,,drop = F],
#                             save.BC = T, 
#                             method = "%difference",
#                             test.BC = F,
#                             test.t.perm = F,
#                             nperm = 0)
# 
# tbi_test$BCD.mat
# ABUND_total_D_diff_mtx[date_1, date_2]
# ABUND_total_C_gain_mtx[date_1, date_2]
# ABUND_total_B_loss_mtx[date_1, date_2]
# X_ref_vegan[date_1, date_2]
# 
# ABUND_total_D_diff_mtx[date_1,date_2] - tbi_test$BCD.mat$`D=(B+C)/(2A+B+C)`
# ABUND_total_B_loss_mtx[date_1,date_2] - tbi_test$BCD.mat$`B/(2A+B+C)`
# ABUND_total_C_gain_mtx[date_1,date_2] - tbi_test$BCD.mat$`C/(2A+B+C)`
# 
# 
# 
# ABUND_total_norm_long_test <- mtx %>%
#                               as.data.frame() %>%
#                               rownames_to_column("Date_formatted") %>%
#                               gather(key = "asv_id", value = "abund", 2:(dim(mtx)[2] +1)) %>%
#                               filter(Date_formatted %in% c(date_1, date_2)) %>%
#                               mutate(time = if_else(Date_formatted == date_1, 1, 2),
#                                      abund = as.numeric(abund))
# 
# turnover_total_ref <- codyn::turnover(df = ABUND_total_norm_long_test,
#                                        time.var = "time",
#                                        species.var = "asv_id",
#                                        abundance.var = "abund",
#                                        metric = "total")
# 
# turnover_appear_ref <- codyn::turnover(df = ABUND_total_norm_long_test,
#                                         time.var = "time",
#                                         species.var = "asv_id",
#                                         abundance.var = "abund",
#                                         metric = "appearance")
# 
# turnover_disappear_ref <- codyn::turnover(df = ABUND_total_norm_long_test,
#                                            time.var = "time",
#                                            species.var = "asv_id",
#                                            abundance.var = "abund",
#                                            metric = "disappearance")
# 
# sp_rich_1 <- sum(mtx[date_1,,drop = F] > 0)
# sp_rich_2 <- sum(mtx[date_2,,drop = F] > 0)
# sp_rich_total <- sum( (mtx[date_1,,drop = F] + mtx[date_2,,drop = F]) > 0)
# sp_rich_appear <- sum(mtx[date_1,,drop = F] == 0 &  mtx[date_2,,drop = F] > 0)
# sp_rich_disappear <- sum(mtx[date_1,,drop = F]> 0 &  mtx[date_2,,drop = F] == 0)
# 
# (sp_rich_appear + sp_rich_disappear)/sp_rich_total
# (sp_rich_appear)/sp_rich_total
# (sp_rich_disappear)/sp_rich_total
# 
# turnover_total_ref
# turnover_appear_ref
# turnover_disappear_ref
# turnover_appear_ref[1,1] + turnover_disappear_ref[1,1] == turnover_total_ref[1,1]


