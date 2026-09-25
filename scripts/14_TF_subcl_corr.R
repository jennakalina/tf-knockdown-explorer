library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(patchwork)

# Read in data
dams <- read.csv('results/diff_analysis/DAMs_all_discovery.csv')
daps <- read.csv('results/diff_analysis/DAPs_all_discovery.csv')
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv', check.names = FALSE)
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv('raw_data/kegg_mapping/kegg_pathways_manual.csv')
clust_mem_met <- read.csv('results/corr_analysis/clust_mem/cluster_mem_tf_metab_k7.csv')
clust_mem_prot <- read.csv('results/corr_analysis/clust_mem/cluster_mem_tf_prot_k6.csv')

# Take top RNAi per df (don't run to do RNAi level)
dams_top_rnai <- dams %>%
  filter(Reg != 'Not significant') %>%
  filter(!Target %in% c('Attp40', 'KH017-Attp2')) %>%
  count(RNAi, Target) %>%
  group_by(Target) %>%
  slice_max(n, n=1, with_ties = FALSE) %>%
  pull(RNAi)
dams <- dams %>% filter(RNAi %in% dams_top_rnai)

daps_top_rnai <- daps %>%
  filter(Reg != 'Not significant') %>%
  rename(Target = target) %>%
  count(RNAi, Target) %>%
  group_by(Target) %>%
  slice_max(n, n=1, with_ties = FALSE) %>%
  pull(RNAi)
daps <- daps %>% filter(RNAi %in% daps_top_rnai) %>% rename(Target = target)

# Unnest pathways
valid_sub <- path_classes %>% 
  filter(class == 'Metabolism',
         subclass != 'Global and overview maps') %>% 
  pull(subclass) %>%
  unique()

metab_path_all_long <- metab_path_all %>%
  mutate(subclasses = strsplit(subclasses, ";\\s*")) %>% 
  unnest(subclasses) %>%
  rename(subclass = subclasses) %>%
  select(metabolite, subclass) %>%
  filter(subclass %in% valid_sub)
prot_path_all_long <- prot_path_all %>%
  mutate(subclasses = strsplit(subclasses, ";\\s*")) %>% 
  unnest(subclasses) %>%
  rename(subclass = subclasses) %>%
  select(protein, subclass) %>%
  filter(subclass %in% valid_sub)

dams <- dams %>% 
  select(metabolite, Target, Reg) %>% # Switch for RNAi level
  left_join(metab_path_all_long, by = 'metabolite')
da_scores_met <- dams %>%
  group_by(Target, subclass) %>% # Switch for RNAi level
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
full_to_nick <- data.frame(subclass = unique(da_scores_met$subclass),
                           sub_short = c('Amino Acid', 'Carbohydrate', 'Energy', 'Glycan', 'Lipid', 
                                         'Cofactors/Vitamins', 'Nucleotide', 'Xenobiotics'))
da_scores_met <- da_scores_met %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

# Add pathway data to DAPs
daps <- daps %>% 
  mutate(Reg = ifelse(is.na(Reg), 'Not significant', Reg)) %>%
  select(protein, Target, Reg) %>% # Switch for RNAi
  left_join(prot_path_all_long, by = 'protein')
da_scores_prot <- daps %>%
  group_by(Target, subclass) %>% # Switch for RNAi
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
da_scores_prot <- da_scores_prot %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

# Calculate DA scores
da_scores_met <- da_scores_met %>%
  mutate(da_score_met = (n_up - n_down) / n_total) %>%
  select(Target, subclass, da_score_met)

da_scores_prot <- da_scores_prot %>%
  mutate(da_score_prot = (n_up - n_down) / n_total) %>%
  select(Target, subclass, da_score_prot)

tfs_keep <- intersect(unique(dams$Target), unique(daps$Target)) # Switch for RNAi
met_clusters <- unique(clust_mem_met$cluster)

# All in metab
for (i in 1:length(met_clusters)) {
  clust <- met_clusters[[i]]
  
  clust_tfs <- clust_mem_met %>% filter(cluster == clust) %>% pull(rnai) ## Include for per cluster
  tfs_keep_sub <- tfs_keep[tfs_keep %in% clust_tfs]
  
  plot_df <- da_scores_met %>% 
    left_join(da_scores_prot, by = c('Target', 'subclass')) %>% # Switcgh for RNAi
    filter(Target %in% tfs_keep_sub)
  
  subclasses <- unique(plot_df$subclass)
  plots <- list()
  
  for (i in 1:length(subclasses)) {
    sub <- subclasses[i]
    
    sub_df <- plot_df %>% filter(subclass == sub)
    
    p <- ggplot(sub_df, aes(x = da_score_met, y = da_score_prot)) +
      geom_vline(xintercept = 0, color = 'gray70') +
      geom_hline(yintercept = 0, color = 'gray70') +
      geom_point() +
      xlim(c(-0.2, 0.55)) + ylim(c(-0.25, 0.4)) +
      geom_smooth(method = 'lm', formula = 'y ~ x', fullrange = TRUE) +
      labs(x = 'Metabolite DA Score', y = 'Protein DA Score', title = paste0(sub, ' Correlation')) +
      theme_bw() 
    
    m <- lm(da_score_prot ~ da_score_met, sub_df)
    r2 <- format(summary(m)$r.squared, digits = 3)
    r2_lab <- paste0('R^2 == ', r2)
    
    p1 <- p + geom_text(x = 0.35, y = 0.3, label = r2_lab, parse = TRUE)
    
    plots[[i]] <- p1
  }
  
  
  p <- wrap_plots(plots, ncol = 4)
  ggsave(paste0('results/tf_systems/correlations/DAscore_corr_met_clust_', clust, '.png'), p, width = 12, height = 6)
}

# In prot
prot_clusters <- unique(clust_mem_prot$cluster)

for (i in 1:length(prot_clusters)) {
  clust <- prot_clusters[[i]]
  
  clust_tfs <- clust_mem_prot %>% filter(cluster == clust) %>% pull(rnai) ## Include for per cluster
  tfs_keep_sub <- tfs_keep[tfs_keep %in% clust_tfs]
  
  plot_df <- da_scores_met %>% 
    left_join(da_scores_prot, by = c('Target', 'subclass')) %>% # Switcgh for RNAi
    filter(Target %in% tfs_keep_sub)
  
  subclasses <- unique(plot_df$subclass)
  plots <- list()
  
  for (i in 1:length(subclasses)) {
    sub <- subclasses[i]
    
    sub_df <- plot_df %>% filter(subclass == sub)
    
    p <- ggplot(sub_df, aes(x = da_score_met, y = da_score_prot)) +
      geom_vline(xintercept = 0, color = 'gray70') +
      geom_hline(yintercept = 0, color = 'gray70') +
      geom_point() +
      xlim(c(-0.2, 0.55)) + ylim(c(-0.25, 0.4)) +
      geom_smooth(method = 'lm', formula = 'y ~ x', fullrange = TRUE) +
      labs(x = 'Metabolite DA Score', y = 'Protein DA Score', title = paste0(sub, ' Correlation')) +
      theme_bw() 
    
    m <- lm(da_score_prot ~ da_score_met, sub_df)
    r2 <- format(summary(m)$r.squared, digits = 3)
    r2_lab <- paste0('R^2 == ', r2)
    
    p1 <- p + geom_text(x = 0.35, y = 0.3, label = r2_lab, parse = TRUE)
    
    plots[[i]] <- p1
  }
  
  
  p <- wrap_plots(plots, ncol = 4)
  ggsave(paste0('results/tf_systems/correlations/DAscore_corr_prot_clust_', clust, '.png'), p, width = 12, height = 6)
}




## Repeat for DF Score
da_scores_met <- dams %>%
  group_by(Target, subclass) %>% # Switch for RNAi level
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
da_scores_met <- da_scores_met %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

# Add pathway data to DAPs
daps <- daps %>% 
  mutate(Reg = ifelse(is.na(Reg), 'Not significant', Reg)) %>%
  select(protein, Target, Reg) %>% # Switch for RNAi
  left_join(prot_path_all_long, by = 'protein')
da_scores_prot <- daps %>%
  group_by(Target, subclass) %>% # Switch for RNAi
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
da_scores_prot <- da_scores_prot %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

# Calculate DF scores
df_scores_met <- da_scores_met %>%
  mutate(df_score_met = (n_up + n_down) / n_total) %>%
  select(Target, subclass, df_score_met)

df_scores_prot <- da_scores_prot %>%
  mutate(df_score_prot = (n_up + n_down) / n_total) %>%
  select(Target, subclass, df_score_prot)

tfs_keep <- intersect(unique(dams$Target), unique(daps$Target)) # Switch for RNAi

plot_df <- df_scores_met %>% 
  left_join(df_scores_prot, by = c('Target', 'subclass')) %>% # Switcgh for RNAi
  filter(Target %in% tfs_keep)

subclasses <- unique(plot_df$subclass)
plots <- list()

for (i in 1:length(subclasses)) {
  sub <- subclasses[i]
  
  sub_df <- plot_df %>% filter(subclass == sub)
  
  p <- ggplot(sub_df, aes(x = df_score_met, y = df_score_prot)) +
    geom_vline(xintercept = 0, color = 'gray70') +
    geom_hline(yintercept = 0, color = 'gray70') +
    geom_point() +
    xlim(c(0, 0.55)) + ylim(c(0, 0.8)) +
    geom_smooth(method = 'lm', formula = 'y ~ x', fullrange = TRUE) +
    labs(x = 'Metabolite DF Score', y = 'Protein DF Score', title = paste0(sub, ' Correlation')) +
    theme_bw() 
  
  m <- lm(df_score_prot ~ df_score_met, sub_df)
  r2 <- format(summary(m)$r.squared, digits = 3)
  r2_lab <- paste0('R^2 == ', r2)
  
  p1 <- p + geom_text(x = 0.35, y = 0.7, label = r2_lab, parse = TRUE)
  
  plots[[i]] <- p1
}


p <- wrap_plots(plots, ncol = 4)
p
ggsave('results/tf_systems/correlations/DFscore_scatterplots.png', width = 12, height = 6)




