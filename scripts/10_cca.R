### Regularized canonical correlation analysis
library(mixOmics)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggtext)
library(ggrepel)
library(ggforce)
library(patchwork)

# Read in data
metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE)
metab_metadata <- read.csv('processed_data/filtered_metab_metadata.csv')
prot_data <- read.csv('processed_data/filtered_prot_data.csv', check.names = FALSE)
prot_metadata <- read.csv('processed_data/filtered_prot_metadata.csv')
prot_paths <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
metab_paths <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')

# Convert to matrix
prot_data <- prot_data %>%
  left_join(prot_metadata %>% dplyr::select(Sample, Sample.info1), by = 'Sample') %>%
  dplyr::select(-Sample) %>%
  relocate(Sample.info1) %>%
  rename(Sample = Sample.info1)

samps_keep <- intersect(prot_data$Sample, metab_data$Sample) # 159 overlapping samples

prot_data <- prot_data %>% 
  filter(Sample %in% samps_keep) %>% 
  dplyr::select(-RNAi) %>% arrange(Sample) %>% 
  column_to_rownames('Sample') %>% 
  as.matrix()
metab_data <- metab_data %>% 
  filter(Sample %in% samps_keep) %>% 
  dplyr::select(-RNAi) %>% arrange(Sample) %>% 
  column_to_rownames('Sample') %>% 
  as.matrix()


######### CCA ######### 
### Sample level
res.cca <- rcc(metab_data, prot_data, method = 'ridge', 
               lambda1 = 0.5, lambda2 = 0.05)

### RNAi Level (average RNAi replicates)
prot_data_rnai <- prot_data %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metab_metadata %>% select(Sample, RNAi), by = 'Sample') %>%
  select(-Sample) %>%
  group_by(RNAi) %>%
  summarise(across(everything(), ~ mean(.x))) %>%
  column_to_rownames('RNAi') %>%
  as.matrix()
metab_data_rnai <- metab_data %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metab_metadata %>% select(Sample, RNAi), by = 'Sample') %>%
  select(-Sample) %>%
  group_by(RNAi) %>%
  summarise(across(everything(), ~ mean(.x))) %>%
  column_to_rownames('RNAi') %>%
  as.matrix()

res.cca.rnai <- rcc(metab_data_rnai, prot_data_rnai, method = 'ridge', 
                    lambda1 = 0.5, lambda2 = 0.05)

### TARGET LEVEL
# Take most responsive RNAi per TF
all_samps <- data.frame(Sample = rownames(metab_data))
all_tfs <- all_samps %>% left_join(metab_metadata %>% dplyr::select(Sample, RNAi, Target), by = 'Sample') 

regsum_met <- read.csv('results/diff_analysis/reg_summary_metab_rnai.csv')

rnai_remove <- all_tfs %>%
  distinct(RNAi, Target) %>%
  group_by(Target) %>%
  filter(n() > 1) %>%
  left_join(regsum_met %>% dplyr::select(RNAi, sig), by = "RNAi") %>%
  filter(sig < max(sig)) %>%
  pull(RNAi)

prot_data_tf <- prot_data %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metab_metadata %>% dplyr::select(Sample, RNAi, Target), by = 'Sample') %>%
  filter(!RNAi %in% rnai_remove) %>%
  dplyr::select(-c(Sample, RNAi)) %>%
  group_by(Target) %>%
  summarise(across(everything(), ~ mean(.x))) %>%
  column_to_rownames('Target') %>%
  as.matrix()
metab_data_tf <- metab_data %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metab_metadata %>% dplyr::select(Sample, RNAi, Target), by = 'Sample') %>%
  filter(!RNAi %in% rnai_remove) %>%
  dplyr::select(-c(Sample, RNAi)) %>%
  group_by(Target) %>%
  summarise(across(everything(), ~ mean(.x))) %>%
  column_to_rownames('Target') %>%
  as.matrix()

res.cca.tf <- rcc(metab_data_tf, prot_data_tf, method = 'ridge', 
                  lambda1 = 0.5, lambda2 = 0.05)


######### PLOTTING ######### 
# Sample level
metab_variates <- res.cca$variates$X %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample') %>%
  left_join(metab_metadata %>% select(Sample, Target, RNAi), by = 'Sample')

# Sample level, TF_sample labeling
ggplot(metab_variates, aes(x = V1, y = V2, label = paste0(Target, '_', Sample))) +
  geom_text_repel(aes(colour = Target), size = 3, show.legend = FALSE,
                  box.padding = 0.01, max.overlaps = 50) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(x = 'Dimension 1', y = 'Dimension 2') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(. ~ "Sample Projection Plot - Sample Level") +
  theme(strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 14, face = 'bold'))
ggsave('results/rCCA/sample_projection_plots/spp_samp_level_samp_labels.png', width = 12, height = 12)

# Sample level, TF only labeling
ggplot(metab_variates, aes(x = V1, y = V2, label = Target)) +
  geom_text_repel(aes(colour = Target), size = 3, show.legend = FALSE,
                  box.padding = 0.01, max.overlaps = 50) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(x = 'Dimension 1', y = 'Dimension 2') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(. ~ "Sample Projection Plot - Sample Level, TF Labels") +
  theme(strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 14, face = 'bold'))
ggsave('results/rCCA/sample_projection_plots/spp_samp_level_tf_labels.png', width = 10, height = 10)

# RNAi level
metab_variates <- res.cca.rnai$variates$X %>% 
  as.data.frame() %>% 
  rownames_to_column('RNAi') %>%
  left_join(metab_metadata %>% select(Target, RNAi) %>% distinct(), by = 'RNAi')

# RNAi level, RNAi labeling
ggplot(metab_variates, aes(x = V1, y = V2, label = RNAi)) +
  geom_text_repel(aes(colour = Target), size = 3, show.legend = FALSE,
                  box.padding = 0.01, max.overlaps = 50) +  
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(x = 'Dimension 1', y = 'Dimension 2') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(. ~ "Sample Projection Plot - RNAi Level") +
  theme(strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 14, face = 'bold'))
ggsave('results/rCCA/sample_projection_plots/spp_rnai_level_rnai_labels.png', width = 10, height = 10)
ggsave('results/for_figures/rCCA/spp_rnai_level_rnai_labels.png', width = 8, height = 8)


# RNAi level, TF labeling
ggplot(metab_variates, aes(x = V1, y = V2, label = Target)) +
  geom_text_repel(aes(colour = Target), size = 3, show.legend = FALSE,
                  box.padding = 0.01, max.overlaps = 50) +  
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(x = 'Dimension 1', y = 'Dimension 2') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(. ~ "Sample Projection Plot - RNAi Level, TF Labels") +
  theme(strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 14, face = 'bold'))
ggsave('results/rCCA/sample_projection_plots/spp_rnai_level_tf_labels.png', width = 8, height = 8)

# TF level
metab_variates <- res.cca.tf$variates$X %>% 
  as.data.frame() %>% 
  rownames_to_column('Target')

ggplot(metab_variates, aes(x = V1, y = V2, label = Target)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_label_repel(size = 3, show.legend = FALSE, label.size = 0, colour = 'midnightblue',
                  box.padding = 0.01, max.overlaps = 50, fill = alpha("white", 0.65)) +  
  labs(x = 'Dimension 1', y = 'Dimension 2') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(. ~ "Sample Projection Plot - TF Level") +
  theme(strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 14, face = 'bold'))
ggsave('results/rCCA/sample_projection_plots/spp_tf_level.png', width = 8, height = 7)
write.csv(metab_variates, 'results/rCCA/tf_level_spp_coords.csv', row.names = FALSE)

### Correlation circle plot
# Function to plot CCPs at different levels (uncomment ggsaves to save)
plot_ccp <- function(level) {
  if (level == 'Sample') {
    res.cca.plot <- res.cca
    metab.data <- metab_data
    prot.data <- prot_data
  } else if (level == 'RNAi') {
    res.cca.plot <- res.cca.rnai
    metab.data <- metab_data_rnai
    prot.data <- prot_data_rnai
  } else if (level == 'TF') {
    res.cca.plot <- res.cca.tf
    metab.data <- metab_data_tf
    prot.data <- prot_data_tf
    }
  
  metab_corr <- cor(metab.data, res.cca.plot$variates$X) %>% 
    as.data.frame() %>%
    rownames_to_column('compound') %>%
    mutate(type = 'metabolite')
  prot_corr <- cor(prot.data, res.cca.plot$variates$Y) %>% 
    as.data.frame() %>%
    rownames_to_column('compound') %>%
    mutate(type = 'protein')
  ccp_df <- rbind(metab_corr, prot_corr)
  
  # Full labels
  p1 <- ggplot(ccp_df, aes(x = V1, y = V2, label = compound)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.35) + 
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.35) +
    geom_circle(aes(x0 = 0, y0 = 0, r = 0.5), colour = 'gray50') +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), colour = 'gray50') +
    geom_text_repel(aes(colour = type), size = 1, show.legend = FALSE,
                    box.padding = 0.001, max.overlaps = 500) +  
    labs(x = 'Component 1', y = 'Component 2') +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    facet_grid(. ~ paste0("Correlation Circle Plot - Full Grid, ", level, " Level")) +
    theme(strip.background = element_rect(fill = "grey90", color = "black"),
          strip.text = element_text(size = 12, face = 'bold'))
  #ggsave(paste0('results/rCCA/correlation_circle_plots/', tolower(level), '_level/ccp_full_labels.png'), p1, width=16,height=16)
  
  # Full points
  p2 <- ggplot(ccp_df, aes(x = V1, y = V2, label = compound)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.35) + 
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.35) +
    geom_circle(aes(x0 = 0, y0 = 0, r = 0.5), colour = 'gray50') +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), colour = 'gray50') +
    geom_point(aes(colour = type)) +  
    labs(x = 'Component 1', y = 'Component 2') +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    facet_grid(. ~ paste0("Correlation Circle Plot - Full Grid, ", level, " Level")) +
    theme(strip.background = element_rect(fill = "grey90", color = "black"),
          strip.text = element_text(size = 12, face = 'bold'))
  #ggsave(paste0('results/rCCA/correlation_circle_plots/', tolower(level), '_level/ccp_full_points.png'), p2, width=8,height=8)
  
  # Lipid metabolism only
  prots_keep <- prot_paths %>% filter(grepl('Lipid metabolism', subclasses)) %>% pull(protein)
  metabs_keep <- metab_paths %>% filter(grepl('Lipid metabolism', subclasses)) %>% pull(metabolite)
  cmpds_keep <- c(prots_keep, metabs_keep)
  
  ccp_df_sub <- ccp_df %>% filter(compound %in% cmpds_keep)
  
  p3 <- ggplot(ccp_df_sub, aes(x = V1, y = V2, label = compound)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.35) + 
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.35) +
    geom_circle(aes(x0 = 0, y0 = 0, r = 0.5), colour = 'gray50') +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), colour = 'gray50') +
    geom_text_repel(aes(colour = type), size = 3, show.legend = FALSE,
                    box.padding = 0.001, max.overlaps = 200) +  
    labs(x = 'Component 1', y = 'Component 2') +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    facet_grid(. ~ paste0("Correlation Circle Plot (", level, " Level, Lipid Metabolism only)")) +
    theme(strip.background = element_rect(fill = "grey90", color = "black"),
          strip.text = element_text(size = 12, face = 'bold'))
    #ggsave(paste0('results/rCCA/correlation_circle_plots/', tolower(level), '_level/ccp_lipid_met.png'), p3, width=10,height=10)
    
    # Amino acid metabolism only
    prots_keep <- prot_paths %>% filter(grepl('Amino acid metabolism', subclasses)) %>% pull(protein)
    metabs_keep <- metab_paths %>% filter(grepl('Amino acid metabolism', subclasses)) %>% pull(metabolite)
    cmpds_keep <- c(prots_keep, metabs_keep)
    
    ccp_df_sub <- ccp_df %>% filter(compound %in% cmpds_keep)
    
    p4 <- ggplot(ccp_df_sub, aes(x = V1, y = V2, label = compound)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.35) + 
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.35) +
      geom_circle(aes(x0 = 0, y0 = 0, r = 0.5), colour = 'gray50') +
      geom_circle(aes(x0 = 0, y0 = 0, r = 1), colour = 'gray50') +
      geom_text_repel(aes(colour = type), size = 3, show.legend = FALSE,
                      box.padding = 0.001, max.overlaps = 200) +  
      labs(x = 'Component 1', y = 'Component 2') +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      facet_grid(. ~ paste0("Correlation Circle Plot (", level, " Level, Amino Acid Metabolism only)")) +
      theme(strip.background = element_rect(fill = "grey90", color = "black"),
            strip.text = element_text(size = 12, face = 'bold'))
    #ggsave(paste0('results/rCCA/correlation_circle_plots/', tolower(level), '_level/ccp_aa_met.png'), p4, width=10,height=10)
    
    plots <- c(p1, p2, p3, p4)
    return(plots)
}

# Viewing plots (run with ggsaves commented out, return_plots uncommented)
plots_samp <- plot_ccp('Sample')
plots_rnai <- plot_ccp('RNAi')
plots_tf <- plot_ccp('TF')

# Saving plots (run with ggsaves uncommented, return_plots commented out)
plot_ccp('Sample')
plot_ccp('RNAi')
plot_ccp('TF')

### CIM plots (built in)
# RNAi level
png('results/rCCA/cim/rnai_level_cim.png', width = 12, height = 10)
cim(res.cca.rnai, comp = 1:2, xlab = 'Proteins', ylab = 'Metabolites')
dev.off()


### Now get clustering
df <- read.csv('results/rCCA/tf_level_spp_coords.csv')

# Look at silhouette width
coord_mat <- df %>% select(-Target) %>% as.matrix()
sil_width <- lapply(2:10, function(k) {
  km <- kmeans(coord_mat, centers = k, nstart = 50)
  sil <- silhouette(km$cluster, dist(coord_mat))
  data.frame(k = k, avg_sil_width = mean(sil[, 'sil_width']))
  }) %>%
  bind_rows()

ggplot(sil_width, aes(x = k, y = avg_sil_width)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 2:10) +
  theme_bw() +
  labs(title = 'Silhouette Width for k = 2 to k = 10',
       x = 'Number of clusters (k)',
       y = 'Silhouette Width')

# sil width shows k = 6 is best
k <- 6

# Cluster with k = 6
df <- df %>% tibble::column_to_rownames('Target')

km_res <- kmeans(df, centers = k)
df$cluster <- km_res$cluster

df <- df %>% tibble::rownames_to_column('TF')
ggplot(df, aes(x = V1, y = V2)) +
  geom_point(aes(colour = as.character(cluster))) +
  scale_color_brewer(palette = "Dark2") +
  geom_label_repel(aes(label = TF, colour = as.character(cluster)), show.legend = FALSE) +
  theme_bw() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(title = 'Sample Projection Plot by Cluster',
       x = 'Dimension 1',
       y = 'Dimension 2',
       color = 'Cluster')
ggsave('results/rCCA/clustering/SPP_clustered.png', width = 10, height = 8)

df_write <- df %>% select(TF, cluster) %>% arrange(cluster)
write.csv(df_write, 'results/rCCA/clustering/cluster_membership.csv', row.names = FALSE)

