library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)

# Read in data
dams <- read.csv('results/diff_analysis/DAMs_all_discovery.csv')
daps <- read.csv('results/diff_analysis/DAPs_all_discovery.csv')
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv', check.names = FALSE)
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv('raw_data/kegg_mapping/kegg_pathways_manual.csv')
clust_mem <- read.csv('results/rCCA/clustering/cluster_membership.csv')

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
daps <- daps %>% filter(RNAi %in% dams_top_rnai) %>% rename(Target = target)

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
clusters <- unique(clust_mem$cluster)

# All clusters individually (will save PNGs if run)
for (i in 1:length(clusters)) {
  clust <- clusters[[i]]
  
  clust_tfs <- clust_mem %>% filter(cluster == clust) %>% pull(TF) ## Include for per cluster
  tfs_keep_sub <- tfs_keep[tfs_keep %in% clust_tfs]
  
  plot_df <- da_scores_met %>% 
    left_join(da_scores_prot, by = c('Target', 'subclass')) %>% # Switch for RNAi
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
  ggsave(paste0('results/rCCA/clustering/correlations/DAscore_corr_clust_', clust, '.png'), p, width = 12, height = 6)
}

### All in one
subclasses <- unique(da_scores_met$subclass)
plot_list <- list()
r2_df <- data.frame(subclass = NA,
                    cluster = NA,
                    n = NA,
                    r2 = NA)

for (sub in subclasses) {
  for (clust in clusters) {
    clust_tfs <- clust_mem %>%
      filter(cluster == clust) %>%
      pull(TF)
    
    tfs_keep_sub <- tfs_keep[tfs_keep %in% clust_tfs]
    
    sub_df <- da_scores_met %>%
      left_join(da_scores_prot,
                by = c("Target", "subclass")) %>%
      filter(Target %in% tfs_keep_sub,
             subclass == sub)
    
    if (nrow(sub_df) > 1) {
      
      m <- lm(da_score_prot ~ da_score_met, sub_df)
      r2 <- format(summary(m)$r.squared, digits = 3)
      
      # Save r2 to df
      r2_df <- bind_rows(r2_df,
                         tibble(subclass = sub,
                                cluster = clust,
                                n = nrow(sub_df),
                                r2 = r2))
      
      p <- ggplot(sub_df,
                  aes(x = da_score_met,
                      y = da_score_prot)) +
        geom_vline(xintercept = 0, color = "gray70") +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_point() +
        geom_smooth(method = "lm",
                    formula = y ~ x,
                    fullrange = TRUE) +
        coord_cartesian(xlim = c(-0.2, 0.55), ylim = c(-0.25, 0.4)) +
        annotate("text", x = 0.35, y = 0.3,
                 label = paste0("R^2 == ", r2),
                 parse = TRUE) +
        theme_bw() +
        labs(x = 'Metabolite DA Score', y = 'Protein DA Score')
      
      # Add column title to top row
      if (sub == subclasses[1])
        p <- p + ggtitle(paste0('Cluster ', clust))
      
      # Add row label to first column
      if (clust == clusters[1])
        p <- p + ylab(sub)
      else
        p <- p
      
    } else {
      r2_df <- bind_rows(r2_df,
                         tibble(subclass = sub,
                                cluster = clust,
                                n = nrow(sub_df),
                                r2 = NA_real_))
      
      p <- ggplot() + theme_void()
    }
    
    plot_list[[length(plot_list) + 1]] <- p
  }
}

final_plot <- wrap_plots(plot_list, ncol = length(clusters))

ggsave("results/rCCA/clustering/correlations/DAscore_corr_all.png", final_plot, width = 18, height = 20)

r2_mat <- r2_df %>% 
  filter(!is.na(subclass)) %>%
  mutate(cluster = paste0('Cluster ', as.character(cluster)),
         r2 = as.numeric(r2)) %>%
  select(subclass, cluster, r2) %>%
  pivot_wider(names_from = cluster, values_from = r2) %>%
  column_to_rownames('subclass') %>%
  as.matrix()

pdf('results/rCCA/clustering/correlations/r2_heatmap.pdf', width = 6, height = 6)
par(mar = c(10, 10, 10, 10))
ComplexHeatmap::pheatmap(r2_mat,
                         cluster_rows = F,
                         cluster_cols = F,
                         main = expression(R^2 ~ of ~ Cluster ~ per ~ Pathway ~ Subclass),
                         name = 'R Squared',
                         angle_col = '45')
dev.off()


