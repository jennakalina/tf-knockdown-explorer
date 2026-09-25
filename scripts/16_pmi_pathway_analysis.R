library(dplyr)
library(tidyr)
library(stringr)

# Read in
conc <- read.csv('results/concordance/conc_results.csv')
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_path_metonly.csv') # Switch between subcl and path
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_path_metonly.csv') # Switch between subcl and path
# Only focusing on significant results
conc <- conc %>% filter(padj < 0.05)
# Add weight, associated protein paths, associated protein paths
conc <- conc %>%
  mutate(weight = -log10(padj)) %>%
  left_join(metab_path_all %>% rename(pathway_met = pathways)) %>%
  left_join(prot_path_all %>% rename(pathway_prot = pathways))
# Save with new info
write.csv(conc, 'results/concordance/pmi_pathway_analysis/conc_res_sig_with_paths.csv', row.names = FALSE)

# Read in for background
conc_full <- read.csv('results/concordance/conc_results_full.csv')

# Function to run at any level, direction, p value threshold
get_indiv_tables <- function(level, dir, p_thresh) {
  # Read in data
  if (level == 'subclass') {
    metab_path_all <- read.csv('raw_data/kegg_mapping/metab_subcl_metonly.csv') 
    prot_path_all <- read.csv('raw_data/kegg_mapping/prot_subcl_metonly.csv') 
    conc <- read.csv('results/concordance/pmi_pathway_analysis/conc_res_sig_with_subclass.csv') 
  } else if (level == 'pathway') {
    metab_path_all <- read.csv('raw_data/kegg_mapping/metab_path_metonly.csv') 
    prot_path_all <- read.csv('raw_data/kegg_mapping/prot_path_metonly.csv') 
    conc <- read.csv('results/concordance/pmi_pathway_analysis/conc_res_sig_with_paths.csv') 
  }
  
  # perform separately for concordant and discordant
  conc <- conc %>% filter(direction == dir) %>% mutate(label = paste0(protein, '-', metabolite)) %>% filter(padj < p_thresh)
  
  # Initialize matrix of pathway x PMI
  pathways <- conc %>% separate_rows(pathway_met, sep = '; ') %>% filter(pathway_met != 'Unannotated') %>% pull(pathway_met) %>% unique()
  
  mat <- matrix(0, nrow = length(unique(conc$label)), ncol = length(pathways))
  colnames(mat) <- pathways
  rownames(mat) <- unique(conc$label)
  
  # Count matrix; 2 if both analytes map to it, 1 if one does, 0 if none
  for (pair in rownames(mat)) {
    prot <- str_split_fixed(pair, '-', n = 2)[1]
    met <- str_split_fixed(pair, '-', n = 2)[2]
    
    if (level == 'subclass') {
      for (path in colnames(mat)) {
        # If path is found in pathway dataframe, add 1 to count matrix
        if (grepl(path, metab_path_all$subclasses[metab_path_all$metabolite == met])) { 
          mat[pair, path] <- mat[pair, path] + 1
        }
        # Same for protein
        if (grepl(path, prot_path_all$subclasses[prot_path_all$protein == prot])) { 
          mat[pair, path] <- mat[pair, path] + 1
        }
      }
    } else if (level == 'pathway') {
      for (path in colnames(mat)) {
        # If path is found in pathway dataframe, add 1 to count matrix
        if (grepl(path, metab_path_all$pathways[metab_path_all$metabolite == met])) { 
          mat[pair, path] <- mat[pair, path] + 1
        }
        # Same for protein
        if (grepl(path, prot_path_all$pathways[prot_path_all$protein == prot])) { 
          mat[pair, path] <- mat[pair, path] + 1
        }
      }
    }
  }
  
  # Create weights matrix (1D array because one value for each pair)
  weights_mat <- conc %>% select(label, weight) %>% tibble::column_to_rownames('label')
  stopifnot(identical(rownames(mat), rownames(weights_mat))) # Ensure order
  weights <- weights_mat$weight
  weights_mat <- matrix(weights, nrow = length(weights),
                        ncol = length(pathways),
                        dimnames = list(names(weights), pathways))
  
  # Score matrix - counts x weight
  score_mat <- sweep(mat, MARGIN = 1, STATS = weights, FUN = '*')
  
  write.csv(mat, paste0('results/concordance/pmi_pathway_analysis/mats/', dir, '_', level, '_p', p_thresh, '_counts_mat.csv'))
  write.csv(weights_mat, paste0('results/concordance/pmi_pathway_analysis/mats/', dir, '_', level, '_p', p_thresh, '_weights_mat.csv'))
  write.csv(score_mat, paste0('results/concordance/pmi_pathway_analysis/mats/', dir, '_', level, '_p', p_thresh, '_score_mat.csv'))
  
  # Calculate PMI score (actual observed)
  pmi_scores_obs <- colSums(score_mat)
  
  # Permutation test
  all_pairs <- conc_full %>%
    mutate(label = paste0(protein, "-", metabolite)) %>%
    distinct(label, protein, metabolite)
  
  full_mat <- matrix(0, nrow = nrow(all_pairs), ncol = length(pathways), dimnames = list(all_pairs$label, pathways))
  # Calculate Crand 
  for (pair in rownames(full_mat)) {
    prot <- str_split_fixed(pair, '-', n = 2)[1]
    met <- str_split_fixed(pair, '-', n = 2)[2]
    
    if (level == 'subclass') {
      for (path in colnames(full_mat)) {
        # If path is found in pathway dataframe, add 1 to count matrix
        if (grepl(path, metab_path_all$subclasses[metab_path_all$metabolite == met])) { 
          full_mat[pair, path] <- full_mat[pair, path] + 1
        }
        # Same for protein
        if (grepl(path, prot_path_all$subclasses[prot_path_all$protein == prot])) { 
          full_mat[pair, path] <- full_mat[pair, path] + 1 
        }
        }
      } else if (level == 'pathway') {
        for (path in colnames(full_mat)) {
          # If path is found in pathway dataframe, add 1 to count matrix
          if (grepl(path, metab_path_all$pathways[metab_path_all$metabolite == met])) { 
            full_mat[pair, path] <- full_mat[pair, path] + 1
          }
          # Same for protein
          if (grepl(path, prot_path_all$pathways[prot_path_all$protein == prot])) { 
            full_mat[pair, path] <- full_mat[pair, path] + 1 
            }
        }
        }
  }
  
  # Perform permutation
  n_perm <- 10000
  n_obs <- length(weights)
  n_path <- length(pathways)
  
  # Set up matrix for PMI scores
  pmi_null_dist <- matrix(NA_real_, nrow = n_perm, ncol = n_path)
  colnames(pmi_null_dist) <- colnames(full_mat)
  
  for (i in seq_len(n_perm)) {
    # Randomly select n pairs from full list of 40k
    random_idx <- sample(seq_len(nrow(full_mat)),
                         size = n_obs,
                         replace = FALSE)
    
    # Get pathway counts for selected pairs
    crand <- full_mat[random_idx, , drop = FALSE]
    
    # Randomly assign weights
    wrand <- sample(weights, size = n_obs, replace = FALSE)
    
    # Calculate scores
    srand <- sweep(crand, MARGIN = 1, STATS = wrand, FUN = '*')
    
    # Calculate pPMI per pathway
    pmi_null_dist[i, ] <- colSums(srand)
  }
  
  # Calculate z score
  null_mean <- colMeans(pmi_null_dist)
  null_sd <- apply(pmi_null_dist, 2, sd)
  
  z_scores <- (pmi_scores_obs - null_mean) / null_sd
  
  # Calculate empirical permutation p-value
  empirical_p_val <- sapply(seq_along(pmi_scores_obs), function(j) {
    (sum(pmi_null_dist[, j] >= pmi_scores_obs[j]) + 1) / (n_perm + 1)
  })
  
  # Adjust with BH
  empirical_p_val_adj <- p.adjust(empirical_p_val, method = 'BH')
  
  # # Check if distribution is normal
  # par(mfrow = c(3, 4))
  # for (p in colnames(pmi_null_dist)) {
  #   hist(pmi_null_dist[, p], breaks = 50, main = p, xlab = "Permutation PMI score")
  #   abline(v = pmi_scores_obs[p], lwd = 2)
  # }
  
  # Calculate corresponding one-sided p-value
  # one_sided_p <- 1 - pnorm(z_scores)
  
  path_info_df <- data.frame(pathway = pathways, 
                             pPMI = pmi_scores_obs,
                             z.score = z_scores,
                             p.val = empirical_p_val,
                             p.adj = empirical_p_val_adj)
  
  write.csv(path_info_df, paste0('results/concordance/pmi_pathway_analysis/summary_tables/', dir, '_', level, '_p', p_thresh, '_table.csv'), row.names = FALSE)
}

# Run for all combinations
levels <- c('subclass', 'pathway')
dirs <- c('Concordant', 'Discordant')
p_thresholds <- c(0.05, 0.01, 1e-3, 1e-4)

grid <- expand.grid(level = levels,
                    dir = dirs,
                    p_thresh = p_thresholds)

set.seed(123)

# Run over grid
for (i in seq_len(nrow(grid))) {
  level_i <- grid$level[i]
  dir_i <- grid$dir[i]
  p_i <- grid$p_thresh[i]
  
  message("\n\nStarting combination - level: ", level_i, ", direction: ", dir_i, ", p-value: ", p_i)
  
  get_indiv_tables(level = level_i, dir = dir_i, p_thresh = p_i)
  
  message("Finished!")
}

### CAN START FROM HERE after running above for all combinations (takes a long time)
levels <- c('subclass', 'pathway')
p_thresholds <- c(0.05, 0.01, 1e-3, 1e-4)

# Now do stats for each level/p-value (combines conc + disc)
grid <- expand.grid(level = levels,
                    p_thresh = p_thresholds)

for (i in seq_len(nrow(grid))) {
  level_i <- grid$level[i]
  p_i <- grid$p_thresh[i]
  
  conc_res <- read.csv(paste0('results/concordance/pmi_pathway_analysis/summary_tables/Concordant_', level_i, '_p', p_i, '_table.csv'))
  disc_res <- read.csv(paste0('results/concordance/pmi_pathway_analysis/summary_tables/Discordant_', level_i, '_p', p_i, '_table.csv'))
  
  # Rename columns to preserve directionality
  conc_res <- conc_res %>% tibble::column_to_rownames('pathway')
  colnames(conc_res) <- paste0('concordant_', colnames(conc_res))
  disc_res <- disc_res %>% tibble::column_to_rownames('pathway')
  colnames(disc_res) <- paste0('discordant_', colnames(disc_res))
  
  # Ensure order then bind
  disc_res <- disc_res[rownames(conc_res),]
  res <- cbind(conc_res, disc_res)
  
  # Calculate additional stats
  res <- res %>%
    mutate(node.strength = concordant_pPMI + discordant_pPMI,
           node.balance = (concordant_pPMI - discordant_pPMI) / (concordant_pPMI + discordant_pPMI)) %>%
    tibble::rownames_to_column('pathway')
  res$Id <- res$Label <- res$pathway
  res[is.na(res)] <- 0
  
  # Filter out those with no significant results
  res <- res %>% filter(node.strength > 0)
  
  # Write out
  write.csv(res, paste0('results/concordance/pmi_pathway_analysis/tables_for_network/', level_i, '_p', p_i, '_node_table.csv'), row.names = FALSE)
}

# Update PMI table
conc_res <- read.csv('results/concordance/conc_results.csv')

for (i in seq_len(nrow(grid))) {
  level_i <- as.character(grid$level[i])
  p_i <- grid$p_thresh[i]
  
  conc_res_filt <- conc_res %>% filter(padj < p_i)
  
  # Read in mapping files depending on level
  if (level_i == 'subclass') {
    # Subclass 
    metab_path_all <- read.csv('raw_data/kegg_mapping/metab_subcl_metonly.csv') %>% rename(pathways = subclasses)
    prot_path_all <- read.csv('raw_data/kegg_mapping/prot_subcl_metonly.csv') %>% rename(pathways = subclasses)
  } else if (level_i == 'pathway') {
    # Pathway
    metab_path_all <- read.csv('raw_data/kegg_mapping/metab_path_metonly.csv') 
    prot_path_all <- read.csv('raw_data/kegg_mapping/prot_path_metonly.csv') 
  }
  
  # add weight and pathways
  conc_res_filt <- conc_res_filt %>% 
    mutate(weight = -log10(padj)) %>%
    left_join(metab_path_all %>% rename(pathway_met = pathways)) %>%
    left_join(prot_path_all %>% rename(pathway_prot = pathways))
  
  # Extend so pathways are unique
  conc_res_long <- conc_res_filt %>%
    separate_rows(pathway_met, sep = '; ') %>%
    separate_rows(pathway_prot, sep = '; ') %>%
    filter(pathway_met != pathway_prot, # Excluding those that are only annotated in the same pathway
           pathway_met != 'Unannotated',
           pathway_prot != 'Unannotated') 
  
  # Since we want order independent pathways, ignore met vs prot and just make it alphabetical
  pathway_pairs <- conc_res_long %>%
    rowwise() %>%
    mutate(pathway_1 = sort(c(pathway_met, pathway_prot))[1],
           pathway_2 = sort(c(pathway_met, pathway_prot))[2]) %>%
    ungroup() %>%
    group_by(pathway_1, pathway_2) %>% # For each unique pair
    summarise(N_conc = sum(direction == 'Concordant'),
              N_disc = sum(direction == 'Discordant'),
              Score_conc = sum(weight[direction == 'Concordant']),
              Score_disc = sum(weight[direction == 'Discordant']),
              .groups = 'drop') %>%
    mutate(edge.strength = Score_conc + Score_disc,
           edge.balance = (Score_conc - Score_disc) / (Score_conc + Score_disc),
           Weight = edge.strength) %>%
    rename(Source = pathway_1, Target = pathway_2)
  pathway_pairs[is.na(pathway_pairs)] <- 0
  
  # Make sure set of pathways match
  nodes_table <- read.csv(paste0('results/concordance/pmi_pathway_analysis/tables_for_network/', level_i, '_p', p_i, '_node_table.csv'))
  unique_paths <- nodes_table %>% pull(pathway) %>% unique()
  pathway_pairs <- pathway_pairs %>% filter(Source %in% unique_paths, Target %in% unique_paths)
  
  # Write out
  write.csv(pathway_pairs, paste0('results/concordance/pmi_pathway_analysis/tables_for_network/', level_i, '_p', p_i, '_edge_table.csv'), row.names = FALSE)
}

