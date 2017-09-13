

# Custom factor levels and palettes for plotting
pal_methods <- c("Summary" = "red", "Binary" = "#2f60af")

pal_indices <- c(ARI = "red",
                 NMI = "darkorange",
                 Completeness = "deepskyblue",
                 Homogeneity = "darkblue",
                 V_measure = "lightslateblue",
                 Purity = "darkgreen",
                 Entropy = "purple4",
                 Consensus = "black",
                 Consensus_shuffle = "black",
                 Random = "darkred")

lvls_idx <- c("Consensus", "ARI", "NMI", "V_measure", "Completeness",
              "Homogeneity", "Purity", "Entropy", "Random")

pal_strategy <- c("Method: Summary, Input: H3K4me3" = "red",
                  "Method: Summary, Input: ChromHMM" = "red",
                  "Method: Summary, Input: DNase I" = "red",
                  "Method: Binary, Input: H3K4me3" = "#2f60af",
                  "Method: Binary, Input: ChromHMM" = "#2f60af",
                  "Method: Binary, Input: DNase I" = "#2f60af",
                  "Method: Summary (shuffled), Input: H3K4me3" = "red",
                  "Method: Binary (shuffled), Input: H3K4me3" = "#2f60af",
                  "Method: Summary (shuffled), Input: ChromHMM" = "red",
                  "Method: Binary (shuffled), Input: ChromHMM" = "#2f60af",
                  "Method: Summary (shuffled), Input: DNase I" = "red",
                  "Method: Binary (shuffled), Input: DNase I" = "#2f60af",
                  "Method: Random predictor" = "gray70")

linetype_strategy <- c("Method: Summary, Input: H3K4me3" = 1,
                       "Method: Summary, Input: ChromHMM" = 2,
                       "Method: Summary, Input: DNase I" = 3,
                       "Method: Binary, Input: H3K4me3" = 1,
                       "Method: Binary, Input: DNase I" = 3,
                       "Method: Binary, Input: ChromHMM" = 2,
                       "Method: Summary (shuffled), Input: H3K4me3" = 1,
                       "Method: Binary (shuffled), Input: H3K4me3" = 1,
                       "Method: Summary (shuffled), Input: ChromHMM" = 2,
                       "Method: Binary (shuffled), Input: ChromHMM" = 2,
                       "Method: Summary (shuffled), Input: DNase I" = 3,
                       "Method: Binary (shuffled), Input: DNase I" = 3,
                       "Method: Random predictor" = 1)

linetype_strategy2 <- c("Method: Summary, Input: H3K4me3" = 1,
                        "Method: Summary, Input: ChromHMM" = 1,
                        "Method: Summary, Input: DNase I" = 1,
                        "Method: Binary, Input: H3K4me3" = 1,
                        "Method: Binary, Input: DNase I" = 1,
                        "Method: Binary, Input: ChromHMM" = 1,
                        "Method: Summary (shuffled), Input: H3K4me3" = 2,
                        "Method: Binary (shuffled), Input: H3K4me3" = 2,
                        "Method: Summary (shuffled), Input: ChromHMM" = 2,
                        "Method: Binary (shuffled), Input: ChromHMM" = 2,
                        "Method: Summary (shuffled), Input: DNase I" = 2,
                        "Method: Binary (shuffled), Input: DNase I" = 2,
                        "Method: Random predictor" = 1)

lvls_idx <- c("Consensus", "ARI", "NMI", "V_measure", "Completeness",
              "Homogeneity", "Purity", "Entropy", "Random")

lvls_idx_shuf <- c("Consensus", "ARI", "NMI", "V_measure", "Completeness",
                   "Homogeneity", "Purity", "Entropy", "Random",
                   "Consensus_shuffle",
                   "ARI_shuffle", "NMI_shuffle", "V_measure_shuffle",
                   "Completeness_shuffle", "Homogeneity_shuffle",
                   "Purity_shuffle", "Entropy_shuffle")

getPalette <- brewer.pal(6, name = "Greys") %>% colorRampPalette
pal_iter <- getPalette(100) %>% rev() %>% setNames(seq(1, 100))

# Plot the results of the grid search for tuning parameters
plotGrid <- function(grid) {
  
  rdbu <- brewer.pal(n = 8, name = "RdBu")
  pal_grid <- colorRampPalette(rdbu)(n=100) %>% rev()
  
  par(mfrow = c(1, 4), mar = c(1.5, 1.5, 1.5, 1.5), oma = c(5, 1, 1, 1))
  # palette <- colorRampPalette(c("blue", "white", "darkorchid4"))(n = 100)
  palette <- pal_grid
  titles <- c("F1 Score", "Precision", "Recall", "False Positive Rate")
  
  for (i in 1:4) {
    
    mat <- data.matrix(grid_scores[[i]])
    
    # Plot the false-colour image of the data frame
    image(mat, col = palette, zlim = c(0, 1), axes = FALSE)
    
    # Axis ticks
    axis(2, at = seq(0, 1, by = 0.165), labels = seq(0.2, 0.8, 0.1),
         tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.5, 0))
    axis(1, at = seq(0, 1, by = 0.165), labels = seq(200, 500, 50),
         tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.5, 0), line = 0.2)
    
    # Titles
    mtext(titles[i], outer = FALSE, side = 3, cex = 0.8, line = 0.4)
    
    # Display the value of the cells
    xlabs <- seq(0, 1, by = 0.165)
    ylabs <- seq(0, 1, by = 0.165)
    text(x = rep(xlabs, ncol(mat)), y = rep(ylabs, each = nrow(mat)),
         labels = round(c(mat), digits = 2), cex = 0.9, col = "white")
    
  }
  
  mtext("p", side = 1, outer = TRUE, col = "grey20", line = 0, cex = 0.8)
  mtext("gap", side = 2, outer = TRUE, col = "grey20", line = 0, cex = 0.8)
  
}

# Plot the ROC curves for the subsampling experiment
plotSubsamp <- function(df, col, method, input) {
  
  auc <- mean_auc %>% 
    filter(Method == method) %>% 
    filter(Input == input)
  
  df %>% ggplot(aes(x = FPR, y = TPR, colour = factor(iteration))) +
    geom_line() + 
    scale_colour_manual(values = pal_iter) +
    # geom_smooth(colour = col, method = "loess", span = 0.8, se = FALSE) +
    stat_summary(fun.y = mean, colour = col, geom = "line") +
    facet_wrap(~ dataset, ncol = 6) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    xlab("False Positive Rate") + ylab("True Positive Rate") +
    ggmin::theme_min() +
    theme(legend.position = "none",
          axis.text = element_text(size = rel(0.6)),
          strip.text.x = element_text(size = 8)) +
    geom_text(data = auc,
              aes(x = 0.6, y = 0.25, label = mean_AUC, group = NULL),
              size = 3,
              inherit.aes = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
  
}
