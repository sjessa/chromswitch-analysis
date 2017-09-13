
# Calculate ROC given external cluster validity stats for the 
# test set and the true classes
calcROC2 <- function(stats, class) {
  
  stats %>% 
    dplyr::select(-gene_name) %>% 
    map(~pROC::roc(response = class, predictor = as.numeric(.))) %>%
    map(`[`, c("sensitivities", "specificities")) %>%
    map(as.data.frame) %>%
    map2_df(names(.), function(df, index) mutate(df, index = index)) %>%
    mutate(TPR = sensitivities, FPR = 1 - specificities)
  
}


# Calculate ROC given external cluster validity stats for the 
# test set and the true classes
calcAUC2 <- function(stats, class) {
  
  stats %>% 
    dplyr::select(-gene_name) %>% 
    map(~pROC::roc(response = class, predictor = as.numeric(.))) %>%
    map_df("auc")
  
}


evalSubsamplePerformance <- function(subsample_out, n_other, n_brain) {
  
  roc <- subsample_out %>%
    lapply(calcROC2, benchmark$switch)
  
  roc_tidy <- seq_along(roc) %>%
    lapply(function(i) mutate(roc[[i]],
                              iteration = i,
                              n = n_other + n_brain,
                              n_other = n_other,
                              n_brain = n_brain)) %>%
    bind_rows()
  
  auc <- subsample_out %>%
    lapply(calcAUC2, benchmark$switch)
  
  auc_tidy <- seq_along(auc) %>%
    lapply(function(i) mutate(auc[[i]],
                              iteration = i,
                              n = n_other + n_brain,
                              n_other = n_other,
                              n_brain = n_brain)) %>%
    bind_rows()
  
  return(list(roc = roc_tidy, auc = auc_tidy))
  
}

# Repeat this subsampling for each region, using the feature
# matrices as input, 100 times, and return ROC and AUC info
subsampleFromFeatures <- function(saved_features, n_other, n_brain, idx_other) {
  
  runExp <- function(subsample_idx) {
    
    seq_along(saved_features) %>% 
      map_df(function(i) cluster(ft_mat = saved_features[[i]][subsample_idx, ],
                                 metadata = meta[subsample_idx, ],
                                 region = benchmark[i],
                                 optimal_clusters = TRUE,
                                 heatmap = FALSE)) %>% 
      dplyr::select(Purity, Entropy, ARI, NMI, Homogeneity, Completeness,
                    V_measure, Consensus, gene_name)
    
  }
  
  subsampleOnce <- function(n_other, n_brain) {
    
    subsample_idx <- c(sample(idx_other, n_other, replace = FALSE),
                       sample(idx_brain, n_brain, replace = FALSE))
    
    runExp(subsample_idx)
    
    
  }
  
  runForGivenIdx <- function(j) {
    
    which_other <- comb_other[all_combos[j, ][[1]]] %>% unlist()
    which_brain <- comb_brain[all_combos[j, ][[2]]] %>% unlist()
    subsample_idx <- c(which_other, which_brain)
    
    runExp(subsample_idx)
    
  }
  
  # The number of ways to subsample to these parameters is the number of
  # ways to subsample from one condition multiplied by the # of ways to subsample
  # from the other
  n_combos <- choose(length(idx_other), n_other) * choose(length(idx_brain), n_brain)
  
  if (n_combos < 100) {
    # Then, rather than sampling randomly, repeat the experiment for all possibilities
    
    comb_other <- combn(idx_other, n_other) %>%
      t() %>%
      as.data.frame() %>%
      split(seq(nrow(.))) %>%
      setNames(NULL) %>%
      lapply(setNames, NULL)
    
    
    comb_brain <- combn(idx_brain, n_brain) %>%
      t() %>%
      as.data.frame() %>%
      split(seq(nrow(.))) %>%
      setNames(NULL) %>%
      lapply(setNames, NULL)
    
    comb_other_idx <- seq(1, length(comb_other))
    comb_brain_idx <- seq(1, length(comb_brain))
    all_combos <- expand.grid(comb_other_idx, comb_brain_idx)
    
    subsample_out <- lapply(seq(1:nrow(all_combos)), runForGivenIdx)
    
  } else subsample_out <- replicate(100, subsampleOnce(n_other, n_brain), simplify = FALSE)
  
  evalSubsamplePerformance(subsample_out, n_other, n_brain)
}


callBinaryFromPk <- function(lpk, idx, meta, region, p) {
  
  lpk@samples <- chromswitch:::lpkSamples(lpk)[idx]
  lpk@peaks   <- chromswitch:::lpkPeaks(lpk)[idx]
  
  lpk %>%
    binarizePeaks(p) %>% 
    cluster(metadata = meta,
            region = region,
            heatmap = FALSE,
            optimal_clusters = TRUE,
            n_features = TRUE)
  
}

subsampleFromPeaks <- function(saved_peaks, n_other, n_brain, idx_other) {
  
  runExp <- function(subsample_idx) {
    
    subsamp_out <- seq_along(saved_peaks) %>% 
      map_df(function(i) try(callBinaryFromPk(lpk = saved_peaks[[i]],
                                                     idx = subsample_idx,
                                                     meta = meta[subsample_idx, ],
                                                     region = benchmark[i],
                                                     p = 0.4))) %>% 
      dplyr::select(Purity, Entropy, ARI, NMI, Homogeneity, Completeness,
                    V_measure, Consensus, gene_name, k, n_features)
    
  }
  
  subsampleOnce <- function(n_other, n_brain) {
    
    subsample_idx <- c(sample(idx_other, n_other, replace = FALSE),
                       sample(idx_brain, n_brain, replace = FALSE))
    
    runExp(subsample_idx)
  }
  
  
  runForGivenIdx <- function(j) {
    
    which_other <- comb_other[all_combos[j, ][[1]]] %>% unlist()
    which_brain <- comb_brain[all_combos[j, ][[2]]] %>% unlist()
    subsample_idx <- c(which_other, which_brain)
    
    runExp(subsample_idx)
    
  }
  
  # The number of ways to subsample to these parameters is the number of
  # ways to subsample from one condition multiplied by the # of ways to subsample
  # from the other
  n_combos <- choose(length(idx_other), n_other) * choose(length(idx_brain), n_brain)
  
  if (n_combos < 100) {
    # Then, rather than sampling randomly, repeat the experiment for all possibilities
    
    comb_other <- combn(idx_other, n_other) %>%
      t() %>%
      as.data.frame() %>%
      split(seq(nrow(.))) %>%
      setNames(NULL) %>%
      lapply(setNames, NULL)
    
    
    comb_brain <- combn(idx_brain, n_brain) %>%
      t() %>%
      as.data.frame() %>%
      split(seq(nrow(.))) %>%
      setNames(NULL) %>%
      lapply(setNames, NULL)
    
    comb_other_idx <- seq(1, length(comb_other))
    comb_brain_idx <- seq(1, length(comb_brain))
    all_combos <- expand.grid(comb_other_idx, comb_brain_idx)
    
    subsample_out <- lapply(seq(1:nrow(all_combos)), runForGivenIdx)
    
  } else subsample_out <- replicate(100, subsampleOnce(n_other, n_brain), simplify = FALSE)
  
  scores <- subsample_out %>% 
    lapply(dplyr::select, -k, -n_features) %>% 
    evalSubsamplePerformance(n_other, n_brain)
  
  return(list(out = subsample_out, scores = scores))
}