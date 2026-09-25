### Correlation analysis on filtered metabolomics data 
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(KEGGREST)
library(dendextend)

# Read in data
metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE)
metadata <- read.csv('processed_data/filtered_metab_metadata.csv')
dam_df <- read.csv('results/diff_analysis/DAMs_all.csv')

# Metabolite-metabolite correlation
metab_data <- column_to_rownames(metab_data, var = 'Sample')
metab_mat <- as.matrix(metab_data %>% select(-RNAi)) 
metab_cor <- cor(metab_mat, method = 'pearson', use = 'pairwise.complete.obs')

## Plot
p <- pheatmap::pheatmap(
  metab_cor,
  #annotation_col = metab_annotation,
  #annotation_row = metab_annotation,
  #annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize_row = 8,
  fontsize_col = 8,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Metabolite-to-Metabolite Correlation Heatmap",
  annotation_legend=FALSE,
)

# Clustering the metab-metab correlation heatmap
row_tree <- p$tree_row

# Get best k with silhouette width
library(cluster)

# Convert correlation to distance
dist_mat <- as.dist(1 - metab_cor)
sil_widths <- sapply(2:15, function(k) {
  cl <- cutree(row_tree, k = k)
  mean(silhouette(cl, dist_mat)[, 3])
})
plot(2:15, sil_widths, type = "b",
     main = 'Clusters vs Silhouette for Metabolites',
     xlab = "Number of clusters",
     ylab = "Average Silhouette Width") # 2 is optimal, but could also use 4 visually

k <- 4
clusters <- cutree(row_tree, k = k)
# Convert to dataframe
cluster_df <- data.frame(
  metabolite = names(clusters),
  cluster = clusters
)

# Add compound and pathway info
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')
cluster_annotation <- cluster_df %>%
  left_join(metab_path_all,
            by = "metabolite") %>%
  arrange(cluster)
# write.csv(cluster_annotation, 'results/corr_analysis/clust_mem/cluster_mem_metab_k4.csv', row.names = FALSE)

clust_mem <- cluster_annotation

# Get cluster annotation
metab_annotation <- clust_mem[, 1:2] %>% distinct() %>% column_to_rownames('metabolite') %>% mutate(cluster = as.factor(cluster))

# Define annotation colors for ALL RNAi (keep full colors)
annotation_colors <- list(cluster = RColorBrewer::brewer.pal(n = length(unique(metab_annotation$cluster)), name = "Set2"))

# Name the colors
names(annotation_colors$cluster) <- unique(metab_annotation$cluster)

p_ann <- pheatmap::pheatmap(
  metab_cor,
  annotation_col = metab_annotation,
  annotation_row = metab_annotation,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize_row = 8,
  fontsize_col = 8,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Metabolite-to-Metabolite Correlation Heatmap",
  annotation_legend=TRUE,
)
ggsave('results/corr_analysis/metab_metab_cor.png', plot = p_ann, width=10, height=8)

## TF-TF Correlation Heatmap
# Average metabolite data by RNAi
metab_data$Sample <- rownames(metab_data)
averaged_data <- metab_data %>%
  group_by(RNAi) %>%
  summarise(across(-Sample, mean, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("RNAi") %>%
  as.matrix()

# Calculate RNAi-to-RNAi correlation
rnai_cor <- cor(t(averaged_data), use = "complete.obs")

# Create heatmap with custom labels but full annotations
p2 <- pheatmap::pheatmap(
  rnai_cor,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  labels_row = rownames(rnai_cor),          # Custom row labels (only controls)
  labels_col = colnames(rnai_cor),          # Custom column labels (only controls)
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "RNAi-to-RNAi Correlation Heatmap"
)

# Clustering
library(dendextend)
# Extract the row clustering and convert to dendrogram
row_clust <- p2$tree_row
row_dend <- as.dendrogram(row_clust)

k <- 6  # or any desired number of clusters
row_clusters <- cutree(row_clust, k = k)

cluster_table <- data.frame(
  RNAi_label = names(row_clusters),
  Cluster = row_clusters
)

# Build dendrogram
par(mar = c(8, 4, 4, 2))  # c(bottom, left, top, right) - increased bottom from ~5 to 8
dend_rnai <- plot(row_dend, main = "RNAi Clustering Dendrogram")

# Save
png(filename = "results/corr_analysis/rnai_rnai_clust_dendro.png", width = 20,
    height = 10, units = "in", res = 800)
par(mar = c(8, 4, 4, 2))  # bottom, left, top, right
plot(row_dend, main = "RNAi Clustering Dendrogram")
dev.off()

# Clustering the rnai-rnai correlation heatmap
row_tree2 <- p2$tree_row

# Get best k with silhouette width
library(cluster)

# Convert correlation to distance
dist_mat2 <- as.dist(1 - rnai_cor)
sil_widths2 <- sapply(2:20, function(k) {
  cl <- cutree(row_tree2, k = k)
  mean(silhouette(cl, dist_mat2)[, 3])
})
plot(2:20, sil_widths2, type = "b",
     main = 'Clusters vs Silhouette for RNAi',
     xlab = "Number of clusters",
     ylab = "Average Silhouette Width")

k <- 6
clusters2 <- cutree(row_tree2, k = k)

# Convert to dataframe
cluster_df2 <- data.frame(
  rnai = names(clusters2),
  cluster = clusters2
)
cluster_df2 <- cluster_df2 %>% arrange(cluster)
# write.csv(cluster_df2, 'results/corr_analysis/clust_mem/cluster_mem_rnai_metab_k6.csv', row.names = FALSE)

# Get cluster annotation
metab_annotation <- cluster_df2[, 1:2] %>% distinct() %>% select(-rnai) %>% mutate(cluster = as.factor(cluster))

# Define annotation colors for ALL RNAi (keep full colors)
annotation_colors <- list(cluster = RColorBrewer::brewer.pal(n = length(unique(metab_annotation$cluster)), name = "Set2"))

# Name the colors
names(annotation_colors$cluster) <- unique(metab_annotation$cluster)

## Plot
p2_save <- pheatmap::pheatmap(
  rnai_cor,
  annotation_col = metab_annotation,
  annotation_row = metab_annotation,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  labels_row = rownames(rnai_cor),          # Custom row labels (only controls)
  labels_col = colnames(rnai_cor),          # Custom column labels (only controls)
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "RNAi-to-RNAi Correlation Heatmap"
)

ggsave('results/corr_analysis/rnai_rnai_cor_metab.png', plot = p2_save, width=10, height=8)


## TF-metabolite Correlation Analysis 
# Get -log10 of adjusted p values, sign correlates with up/downregulation
dam_df <- dam_df %>%
  mutate(signed_log10padj = sign(log2FC) * -log10(p.adj))
padj_mat <- dam_df %>%
  select(RNAi, metabolite, signed_log10padj) %>%
  pivot_wider(names_from = metabolite,
              values_from = signed_log10padj) %>%
  column_to_rownames('RNAi') %>%
  as.matrix()

# Plot
library(ComplexHeatmap)
Heatmap(
    t(padj_mat),
    name = "Adjusted p-value",
    col = circlize::colorRamp2(c(-1, 0, 1), c("steelblue3", "white", "salmon2")),
    column_title = "68 Responsive Perturbations",
    row_title = "245 Metabolites",
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    show_row_names = FALSE,
    show_column_names = FALSE
  )
# Save
png(filename = "results/corr_analysis/tf_metab_padj_cor.png", width = 10,
    height = 8, units = "in", res = 800)
par(mar = c(8, 4, 4, 2))  # bottom, left, top, right
Heatmap(
  t(padj_mat),
  name = "Adjusted p-value",
  col = circlize::colorRamp2(c(-1, 0, 1), c("steelblue3", "white", "salmon2")),
  column_title = "101 Responsive Perturbations",
  row_title = "245 Metabolites",
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  show_row_names = FALSE,
  show_column_names = FALSE
)
dev.off()


##### SAME LAST TWO, BUT TF LEVEL
rm(list = ls())
# Read in data
metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE)
metadata <- read.csv('processed_data/filtered_metab_metadata.csv')
dam_df <- read.csv('results/diff_analysis/DAMs_all.csv')
reg_summary <- read.csv('results/diff_analysis/reg_summary_metab_rnai.csv')

# Get best RNAi per TF (highest number of DAMs)
best_rnai_per_tf <- reg_summary %>%
  group_by(Target) %>%
  slice_max(order_by = sig, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Target, RNAi = RNAi)
# Add controls to best rnai to keep
control <- data.frame(Target='control', RNAi='w1118_F')
best_rnai_per_tf <- rbind(best_rnai_per_tf, control)

# Filter data to only keep the best rnais
rnai_to_keep <- best_rnai_per_tf$RNAi
metab_data <- metab_data %>%
  filter(RNAi %in% rnai_to_keep) %>%
  left_join(best_rnai_per_tf, by='RNAi') %>%
  select(-RNAi)

# Average metabolite data by TF
averaged_data <- metab_data %>%
  group_by(Target) %>%
  summarise(across(-Sample, mean, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("Target") %>%
  as.matrix()

# Calculate RNAi-to-RNAi correlation
tf_cor <- cor(t(averaged_data), use = "complete.obs")

# Create heatmap with custom labels but full annotations
p2 <- pheatmap::pheatmap(
  tf_cor,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  labels_row = rownames(tf_cor),          # Custom row labels (only controls)
  labels_col = colnames(tf_cor),          # Custom column labels (only controls)
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "TF-to-TF Correlation Heatmap"
)

# Clustering
# Extract the row clustering and convert to dendrogram
row_clust <- p2$tree_row
row_dend <- as.dendrogram(row_clust)

# Clustering the rnai-rnai correlation heatmap
row_tree2 <- p2$tree_row

# Get best k with silhouette width
library(cluster)

# Convert correlation to distance
dist_mat2 <- as.dist(1 - tf_cor)
sil_widths2 <- sapply(2:20, function(k) {
  cl <- cutree(row_tree2, k = k)
  mean(silhouette(cl, dist_mat2)[, 3])
})
plot(2:20, sil_widths2, type = "b",
     main = 'Clusters vs Silhouette for RNAi',
     xlab = "Number of clusters",
     ylab = "Average Silhouette Width")

k <- 4
clusters2 <- cutree(row_tree2, k = k)

# Convert to dataframe
cluster_df2 <- data.frame(
  rnai = names(clusters2),
  cluster = clusters2
)
cluster_df2 <- cluster_df2 %>% arrange(cluster)
write.csv(cluster_df2, 'results/corr_analysis/clust_mem/cluster_mem_tf_metab_k4.csv', row.names = FALSE)

k <- 4  # or any desired number of clusters
row_clusters <- cutree(row_clust, k = k)

cluster_table <- data.frame(
  RNAi_label = names(row_clusters),
  Cluster = row_clusters
)

# Build dendrogram
par(mar = c(8, 4, 4, 2))  # c(bottom, left, top, right) - increased bottom from ~5 to 8
dend_rnai <- plot(row_dend, main = "TF Clustering Dendrogram")

# Save
png(filename = "results/corr_analysis/tf_to_tf/tf_tf_clust_dendro_metab.png", width = 20,
    height = 10, units = "in", res = 800)
par(mar = c(8, 4, 4, 2))  # bottom, left, top, right
plot(row_dend, main = "TF Clustering Dendrogram")
dev.off()

# Get cluster annotation
metab_annotation <- cluster_df2[, 1:2] %>% distinct() %>% select(-rnai) %>% mutate(cluster = as.factor(cluster))

# Define annotation colors for ALL RNAi (keep full colors)
annotation_colors <- list(cluster = RColorBrewer::brewer.pal(n = length(unique(metab_annotation$cluster)), name = "Set3"))

# Name the colors
names(annotation_colors$cluster) <- unique(metab_annotation$cluster)

## Plot
p2_save <- pheatmap::pheatmap(
  tf_cor,
  annotation_col = metab_annotation,
  annotation_row = metab_annotation,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  labels_row = rownames(tf_cor),          # Custom row labels (only controls)
  labels_col = colnames(tf_cor),          # Custom column labels (only controls)
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "TF-to-TF Correlation Heatmap (Metabolite Data)"
)

ggsave('results/corr_analysis/tf_to_tf/tf_tf_cor_metab_k4.png', plot = p2_save, width=10, height=8)
ggsave('results/for_figures/corr_analysis/tf_tf_cor_metab_k4.png', plot = p2_save, width=10, height=8)

# Get -log10 of adjusted p values, sign correlates with up/downregulation
dam_df <- dam_df %>%
  mutate(signed_log10padj = sign(log2FC) * -log10(p.adj)) %>%
  filter(RNAi %in% rnai_to_keep) %>%
  select(-RNAi) %>%
  relocate(Target)
padj_mat <- dam_df %>%
  select(Target, metabolite, signed_log10padj) %>%
  pivot_wider(names_from = metabolite,
              values_from = signed_log10padj) %>%
  column_to_rownames('Target') %>%
  as.matrix() %>%
  t()

# Label only 3 of interest
row_labels <- rownames(padj_mat)
row_labels[!row_labels %in% c('N6-N6-N6-Trimethyl-L-lysine', 'Citrate', 'L-serine')] <- ''

# Plot
p3_save <- pheatmap::pheatmap(
  padj_mat,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-4, 4, length.out = 101),
  labels_col = colnames(padj_mat),
  #labels_row = row_labels,
  #show_rownames = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "TF-to-Metabolite Correlation Heatmap"
)

ggsave('results/corr_analysis/tf_metab_cor.png', plot = p3_save, width=8, height=8)


## Correlations of just TML
rm(list=ls())
# Read in data
metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE)
metadata <- read.csv('processed_data/filtered_metab_metadata.csv')
dam_df <- read.csv('results/diff_analysis/DAMs_all.csv')

# Metabolite-metabolite correlation
metab_data <- column_to_rownames(metab_data, var = 'Sample')
metab_mat <- as.matrix(metab_data %>% select(-RNAi)) 
metab_cor <- cor(metab_mat, method = 'pearson', use = 'pairwise.complete.obs')

tml_cor <- metab_cor %>% 
  as.data.frame() %>% 
  select('N6-N6-N6-Trimethyl-L-lysine') %>%
  tibble::rownames_to_column('Metabolite') %>%
  rename(correlation = 'N6-N6-N6-Trimethyl-L-lysine') %>%
  filter(Metabolite != 'N6-N6-N6-Trimethyl-L-lysine') %>%
  mutate(abs_corr = abs(correlation),
         sign = ifelse(sign(correlation) == 1, 'Positive', 'Negative'),
         sign = factor(sign, levels = c('Positive', 'Negative')))

top_15_mets <- tml_cor %>% arrange(desc(correlation)) %>% slice_head(n = 15) %>% pull(Metabolite)
bot_15_mets <- tml_cor %>% arrange(correlation) %>% slice_head(n = 15) %>% pull(Metabolite)
label_mets <- c(bot_15_mets, top_15_mets)

tml_cor <- tml_cor %>%
  mutate(label = ifelse(Metabolite %in% label_mets, Metabolite, ''))

ggplot(tml_cor, aes(x = reorder(Metabolite, desc(correlation)), y = correlation)) +
  geom_point(aes(fill = sign), shape = 21) +
  scale_fill_manual(values = c('tomato3', 'cyan3')) +
  theme_bw() +
  geom_text_repel(aes(label = label), max.overlaps = 50, size = 3, 
                  box.padding = 0.5, segment.color = alpha('gray20', 0.5)) +
  geom_hline(yintercept = 0) +
  labs(title = 'Pearson Correlation of Metabolites with TML',
       x = 'Metabolite', y = 'Pearson Coefficient', fill = 'Sign of Pearson Coefficient') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave('results/corr_analysis/tml_correlations.png', width = 10, height = 6)

keep_mets <- c(label_mets, 'N6-N6-N6-Trimethyl-L-lysine')
tml_mets_df <- metab_cor %>% 
  as.data.frame() %>% 
  select(any_of(keep_mets)) %>%
  rownames_to_column('Metabolite') %>%
  filter(Metabolite %in% keep_mets) %>%
  column_to_rownames('Metabolite')
tml_mets_df <- tml_mets_df[colnames(tml_mets_df),]

top_corr_df <- data.frame(metabolite = character(),
                          top5 = character(),
                          direction = character())

for (i in 1:length(top_15_mets)) {
  metab <- top_15_mets[[i]]
  
  met_list <- tml_mets_df %>% 
    select(metab) %>% 
    rename(corr = metab) %>%
    arrange(desc(corr)) %>% 
    rownames_to_column('metabolite') %>%
    filter(metabolite != metab) %>%
    slice_head(n = 5) %>%
    pull(metabolite)
  
  top_corr_df <- top_corr_df %>% add_row(metabolite = metab,
                                         top5 = paste0(met_list, collapse = '; '),
                                         direction = 'Positive')
}

for (i in 1:length(bot_15_mets)) {
  metab <- bot_15_mets[[i]]
  
  met_list <- tml_mets_df %>% 
    select(metab) %>% 
    rename(corr = metab) %>%
    arrange(corr) %>% 
    rownames_to_column('metabolite') %>%
    filter(metabolite != metab) %>%
    slice_head(n = 5) %>%
    pull(metabolite)
  
  top_corr_df <- top_corr_df %>% add_row(metabolite = metab,
                                         top5 = paste0(met_list, collapse = '; '),
                                         direction = 'Negative')
}

write.csv(top_corr_df, 'results/corr_analysis/tml_correlations_top5.csv', row.names = FALSE)

