###############################################################################
### 1. Set env
###############################################################################

library(vegan)
library(tidyverse)

###############################################################################
### 2. Get Season function
###############################################################################

getSeason <- function(input.date){
  numeric.date <- 100*month(input.date)+day(input.date)
  ## input Seasons upper limits in the form MMDD in the "break =" option:
  cuts <- base::cut(numeric.date, breaks = c(0,320,0621,0922,1221,1231)) 
  # rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
  levels(cuts) <- c("Summer","Fall","Winter","Spring","Summer")
  return(cuts)
}

###############################################################################
### 3. Custom PCoA function
###############################################################################

custom_pcoa <- function(ABUND, dataset, normalize = TRUE) {
  
  output_list <- list()
  
  if (normalize == T) {
    ABUND <- decostand(ABUND, method = "hellinger")
  }
  
  ABUND_dist <- vegdist(ABUND, method = "bray")
  ABUND_pcoa <- cmdscale(ABUND_dist, k = 4, eig = T)
  
  X <- data.frame(ABUND_pcoa$points) %>%
       rownames_to_column("sample_name")
  
  colnames(X) <- c("sample_name", "PCo1", "PCo2", "PCo3", "PCo4")
  X$dataset <- dataset
  
  output_list[["X"]] <- X
  output_list[["eig"]] <- ABUND_pcoa$eig
  
  return(output_list)
  
}

###############################################################################
### 4. Custom PCoA plot function
###############################################################################

custom_pcoa_plot <- function(X_output_list, METADATA = METADATA, text_size = 13, 
                             pcoa_plot_margin = c(1,1,1,3),
                             boxplot_margin = c(0.1,14.75,0.2,0.55)) {
  
  output_list <- list()
  
  X_output_list_extend <- X_output_list[["X"]] %>% 
                          left_join(x = ., y = METADATA, by = "sample_name") 
  
  output_list[["X_extend"]] <- X_output_list_extend
  
  # set vars
  
  x_title <- paste("PC1 ", (X_output_list[["eig"]]/sum(X_output_list[["eig"]]))[1] %>% round(4) * 100, "% of variance")
  y_title <- paste("PC2 ", (X_output_list[["eig"]]/sum(X_output_list[["eig"]]))[2] %>% round(4) * 100, "% of variance")
  x_max <- max(X_output_list_extend$PCo1) + max(X_output_list_extend$PCo1)/10
  x_min <- min(X_output_list_extend$PCo1) - abs(min(X_output_list_extend$PCo1)/10)
  y_max <- max(X_output_list_extend$PCo2) + max(X_output_list_extend$PCo2)/10
  y_min <- min(X_output_list_extend$PCo2) - abs(min(X_output_list_extend$PCo2)/10)
  season_colors <- c("#c93f1b","#98482b", "#154360", "#3c7810")
  
  # PCoA plot
  
  pcoa_plot_pco1_vs_pco2 <- ggplot(X_output_list_extend, aes(x = PCo1, y = PCo2, 
                                                             # shape = interaction(Filter_names, Origin, sep = " "),
                                                             color = Season), alpha = 0.8) +
                            geom_point(size = 5) +
                            scale_color_manual(values = season_colors) +
                            # scale_shape_manual(name = "Size fraction", values = c(17,19,15)) +
                            geom_text(aes(label = Date_formatted), vjust = 2, size = 3, check_overlap = TRUE) +
                            xlab(x_title) +
                            ylab(y_title) +                      
                            theme_bw() +
                            theme(
                                 axis.title = element_text(size = text_size +4),
                                 axis.text = element_text(size = text_size),
                                 legend.text = element_text(size= text_size),
                                 legend.title = element_text(size = text_size),
                                 plot.margin = unit(pcoa_plot_margin, "lines")
                           ) +
                           xlim(c(x_min, x_max)) +
                           ylim(c(y_min, y_max))
  
  ### boxplot
  
  x_title_boxplot <- x_title
  y_title_boxplot <- "Season"
  X_output_list_extend$Season_rev <- factor(X_output_list_extend$Season, 
                                            levels = c("Spring", "Winter", "Fall", "Summer" ))
  
  boxplot_pco1_vs_season <- ggplot(X_output_list_extend, aes(x = PCo1, y = Season_rev, fill = Season_rev)) +
                            geom_boxplot() +
                            scale_fill_manual(values = rev(season_colors), name = "Season") +
                            xlab(x_title_boxplot) +
                            ylab(y_title_boxplot) +                      
                            theme_bw() +
                            theme(
                                  axis.title.x = element_text(size = text_size +4, margin = unit(c(4,0,0,0), "mm")),
                                  axis.title.y = element_text(size = text_size +4, margin = unit(c(0,2,0,0), "mm")),
                                  axis.text = element_text(size = text_size),
                                  legend.text = element_text(size= text_size),
                                  legend.title = element_text(size = text_size),
                            plot.margin = unit(boxplot_margin, "lines")
                          ) +
                          xlim(c(x_min, x_max))
  
  # join plots
  
  pcoa_plot_pco1_vs_pco2_and_boxplot <- ggarrange(pcoa_plot_pco1_vs_pco2 + rremove("x.title") +rremove("x.text"), 
                                                  boxplot_pco1_vs_season + rremove("legend"), 
                                                  ncol = 1, nrow = 2, 
                                                  heights = c(0.7, 0.3)) 
  
  output_list[["plot"]] <- pcoa_plot_pco1_vs_pco2_and_boxplot
  
  return(output_list)
    
}

###############################################################################
### Compute angles function
###############################################################################

angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

