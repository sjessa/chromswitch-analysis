
# Make a new 'Label' column in a dataframe, combining the Method & Input
makeLabel <- function(df) {
  
  df %>% rowwise() %>% 
    mutate(Label = case_when(
      Method == "Summary" && Input == "H3K4me3" && index == "Consensus" ~
        "Method: Summary, Input: H3K4me3",
      Method == "Summary" && Input == "ChromHMM" && index == "Consensus" ~
        "Method: Summary, Input: ChromHMM",
      Method == "Summary" && Input == "DNase I" && index == "Consensus" ~
        "Method: Summary, Input: DNase I",
      Method == "Binary" && Input == "H3K4me3" && index == "Consensus" ~
        "Method: Binary, Input: H3K4me3",
      Method == "Binary" && Input == "ChromHMM" && index == "Consensus"~
        "Method: Binary, Input: ChromHMM",
      Method == "Binary" && Input == "DNase I" && index == "Consensus"~
        "Method: Binary, Input: DNase I",
      Method == "Summary" && index == "Consensus_shuffle" && Input == "H3K4me3" ~
        "Method: Summary (shuffled), Input: H3K4me3",
      Method == "Binary" && index == "Consensus_shuffle" && Input == "H3K4me3" ~
        "Method: Binary (shuffled), Input: H3K4me3",
      Method == "Summary" && index == "Consensus_shuffle" && Input == "ChromHMM" ~
        "Method: Summary (shuffled), Input: ChromHMM",
      Method == "Binary" && index == "Consensus_shuffle" && Input == "ChromHMM" ~
        "Method: Binary (shuffled), Input: ChromHMM",
      Method == "Summary" && index == "Consensus_shuffle" && Input == "DNase I" ~
        "Method: Summary (shuffled), Input: DNase I",
      Method == "Binary" && index == "Consensus_shuffle" && Input == "DNase I" ~
        "Method: Binary (shuffled), Input: DNase I",
      index == "Random" ~ "Method: Random predictor"
    )) %>%
    ungroup()
  
}

makeLabelSimple <- function(df) {
  
  df %>% rowwise() %>% 
    mutate(Label = case_when(
      Method == "Summary" && Input == "H3K4me3"  ~ "Method: Summary, Input: H3K4me3",
      Method == "Summary" && Input == "ChromHMM" ~ "Method: Summary, Input: ChromHMM",
      Method == "Summary" && Input == "DNase I"  ~ "Method: Summary, Input: DNase I",
      Method == "Binary"  && Input == "H3K4me3"  ~ "Method: Binary, Input: H3K4me3",
      Method == "Binary"  && Input == "ChromHMM" ~ "Method: Binary, Input: ChromHMM",
      Method == "Binary"  && Input == "DNase I"  ~ "Method: Binary, Input: DNase I"
    )) %>%
    ungroup()
  
}

# From the output of the grid search, calculate the types of decisions
# and calculate performance scores
makeGrid <- function(grid, input) {
  
  grid_results <- grid$result[map_lgl(grid$error, is_null)] %>% 
    map_df("stats") %>% 
    bind_rows() %>% 
    mutate(precision = tp/(tp + fp),
           recall = tp/(tp + fn),
           fpr = fp/(tn + fp),
           f1 = 2/(1/recall + 1/precision)) %>%
    select(-c(tp, tn, fp, fn))
  
  makeGrid <- function(col) {
    
    df <- grid_results %>%
      as.data.frame %>%
      select_("gap", "p", col) %>%
      spread_("p", col)
    
    rownames(df) <- df$gap
    df <- df %>% select(-gap)
    return(df)
  }
  
  eval_scores <- c("F1 Score", "Precision", "Recall", "False Positive Rate")
  
  grid_scores <- c("f1", "precision", "recall", "fpr") %>%
    lapply(makeGrid)
  
  save(grid_scores,
       file = paste0("figure_data/grid_scores.", input, ".RData"))
  
  grid_scores_flat <- lapply(seq_along(grid_scores),
                             function(i) mutate(grid_scores[[i]],
                                                score = eval_scores[i])) %>% 
    bind_rows() %T>%
    write_tsv(paste0("figure_data/grid_scores.", input, ".tsv"))
}