### Differential analysis on filtered proteomics data 
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(cluster)

# Read in data
prot_data <- read.csv('processed_data/filtered_prot_data.csv', check.names = FALSE)
metadata <- read.csv('processed_data/filtered_prot_metadata.csv')
dap_df <- read.csv('results/diff_analysis/DAPs_all.csv')
### NOTE: If running the first time, skip all clustering until after you calculate it later this script

## TF-TF Correlation Heatmap
# Average protein data by RNAi
averaged_data <- prot_data %>%
  group_by(RNAi) %>%
  summarise(across(-Sample, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("RNAi") %>%
  as.matrix()

# Calculate RNAi-to-RNAi correlation
rnai_cor <- cor(t(averaged_data), use = "complete.obs")

# Create heatmap with custom labels but full annotations
p2 <- pheatmap::pheatmap(
  rnai_cor,
  #annotation_row = rnai_annotation,        # Full annotation (colors for all)
  #annotation_col = rnai_annotation,        # Full annotation (colors for all)
  #annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  #labels_row = rownames(rnai_cor),
  #labels_col = colnames(rnai_cor),
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize_row = 7,
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

k <- 4  # or any desired number of clusters
row_clusters <- cutree(row_clust, k = k)

cluster_table <- data.frame(
  RNAi_label = names(row_clusters),
  Cluster = row_clusters
)

# Plot dendrogram
par(mar = c(8, 4, 4, 2))  # c(bottom, left, top, right) - increased bottom from ~5 to 8
dend_rnai <- plot(row_dend, main = "RNAi Clustering Dendrogram")

# Save
png(filename = "results/corr_analysis/rnai_rnai_clust_dendro_prot.png", width = 20,
    height = 10, units = "in", res = 800)
par(mar = c(8, 4, 4, 2))  # bottom, left, top, right
plot(row_dend, main = "RNAi Clustering Dendrogram")
dev.off()

# Clustering the rnai-rnai correlation heatmap
row_tree2 <- p2$tree_row

# Get best k with silhouette width
# Convert correlation to distance
dist_mat2 <- as.dist(1 - rnai_cor)
sil_widths2 <- sapply(2:15, function(k) {
  cl <- cutree(row_tree2, k = k)
  mean(silhouette(cl, dist_mat2)[, 3])
})
plot(2:15, sil_widths2, type = "b",
     main = 'Clusters vs Silhouette for RNAi',
     xlab = "Number of clusters",
     ylab = "Average Silhouette Width")

k <- 4
clusters2 <- cutree(row_tree2, k = k)

# Convert to dataframe
cluster_df2 <- data.frame(
  protein = names(clusters2),
  cluster = clusters2
)
cluster_df2 <- cluster_df2 %>% arrange(cluster)
write.csv(cluster_df2, 'results/corr_analysis/clust_mem/cluster_mem_rnai_prot_k4.csv', row.names = FALSE)

# Get cluster membership for annotation
clust_mem_rnairnai <- read.csv("results/corr_analysis/clust_mem/cluster_mem_rnai_prot_k4.csv") %>% column_to_rownames('protein')
clust_mem_rnairnai$cluster <- factor(clust_mem_rnairnai$cluster,
                                     levels = sort(unique(clust_mem_rnairnai$cluster)))

# Ensure row order matches correlation matrix
rnai_annotation <- clust_mem_rnairnai[rownames(rnai_cor), , drop = FALSE]

# Get cluster levels
cluster_levels <- levels(rnai_annotation$cluster)

# Generate colors
annotation_colors <- list(cluster = setNames(RColorBrewer::brewer.pal(n = length(unique(rnai_annotation$cluster)),name = 'Set2'), 
                                             cluster_levels))

p2_ann <- pheatmap::pheatmap(
  rnai_cor,
  annotation_row = rnai_annotation,        # Full annotation (colors for all)
  annotation_col = rnai_annotation,        # Full annotation (colors for all)
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  #labels_row = rownames(rnai_cor),
  #labels_col = colnames(rnai_cor),
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize_row = 7,
  fontsize_col = 8,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "RNAi-to-RNAi Correlation Heatmap"
)

ggsave('results/corr_analysis/rnai_rnai_cor_prot.png', plot = p2_ann, width=10, height=8)



## TF-protein Correlation Analysis 
# Get -log10 of adjusted p values, sign correlates with up/downregulation
dap_df <- dap_df %>%
  mutate(signed_log10padj = sign(log2FC) * -log10(p.adj))
padj_mat <- dap_df %>%
  select(RNAi, protein, signed_log10padj) %>%
  pivot_wider(names_from = protein,
              values_from = signed_log10padj) %>%
  column_to_rownames('RNAi') %>%
  as.matrix()

# Handle NAs for plotting
m <- t(padj_mat)
m[is.na(m)] <- 0
dist_rows <- dist(m)
dist_cols <- dist(t(m))

# Plot
library(ComplexHeatmap)
Heatmap(
  t(padj_mat),
  name = "Adjusted p-value (NAs in gray)",
  col = circlize::colorRamp2(c(-1, 0, 1), c("steelblue3", "white", "salmon2")),
  na_col = "gray90",
  column_title = "143 Responsive Perturbations",
  row_title = "755 Proteins",
  cluster_rows = hclust(dist_rows),
  cluster_columns = hclust(dist_cols),
  show_row_names = FALSE,
  show_column_names = FALSE
)
# Save
png(filename = "results/corr_analysis/tf_prot_padj_cor.png", width = 10,
   height = 8, units = "in", res = 800)
par(mar = c(8, 4, 4, 2))  # bottom, left, top, right
Heatmap(
  t(padj_mat),
  name = "Adjusted -log10 p-value (NAs in gray)",
  col = circlize::colorRamp2(c(-1, 0, 1), c("steelblue3", "white", "salmon2")),
  na_col = "gray90",
  column_title = "143 Responsive Perturbations",
  row_title = "755 Proteins",
  cluster_rows = hclust(dist_rows),
  cluster_columns = hclust(dist_cols),
  show_row_names = FALSE,
  show_column_names = FALSE
)
dev.off()


# Protein-protein correlation
prot_data <- column_to_rownames(prot_data, var = 'Sample')
prot_mat <- as.matrix(prot_data %>% select(-RNAi)) # Leave out RNAi column 
prot_cor <- cor(prot_mat, method = 'pearson', use = 'pairwise.complete.obs')

# Plot heatmap
library(pheatmap)

p <- pheatmap::pheatmap(
  prot_cor,
  #annotation_row = prot_annotation,        # Full annotation (colors for all)
  #annotation_col = prot_annotation,        # Full annotation (colors for all)
  #annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Protein-to-Protein Correlation Heatmap",
  annotation_legend = TRUE
)

# Clustering the prot-prot correlation heatmap
row_tree <- p$tree_row

# Get best k with silhouette width
# Convert correlation to distance
dist_mat <- as.dist(1 - prot_cor)
sil_widths <- sapply(2:15, function(k) {
  cl <- cutree(row_tree, k = k)
  mean(silhouette(cl, dist_mat)[, 3])
})
plot(2:15, sil_widths, type = "b",
     main = 'Clusters vs Silhouette for Protein',
     xlab = "Number of clusters",
     ylab = "Average Silhouette Width")

k <- 4
clusters <- cutree(row_tree, k = k)

# Convert to dataframe
cluster_df <- data.frame(
  protein = names(clusters),
  cluster = clusters
)

# Add compound and pathway info
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
cluster_annotation <- cluster_df %>%
  left_join(prot_path_all,
            by = "protein") %>%
  arrange(cluster)
write.csv(cluster_annotation, 'results/corr_analysis/clust_mem/cluster_mem_prot_k4.csv', row.names = FALSE)

# Ensure row order matches correlation matrix
prot_annotation <- cluster_annotation[, 1:2] %>% distinct() %>% column_to_rownames('protein') %>% mutate(cluster = as.factor(cluster))

# Define annotation colors for ALL RNAi (keep full colors)
annotation_colors <- list(cluster = RColorBrewer::brewer.pal(n = length(unique(prot_annotation$cluster)), name = "Set2"))

# Name the colors
names(annotation_colors$cluster) <- unique(prot_annotation$cluster)

p_ann <- pheatmap::pheatmap(
  prot_cor,
  annotation_row = prot_annotation,        # Full annotation (colors for all)
  annotation_col = prot_annotation,        # Full annotation (colors for all)
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Protein-to-Protein Correlation Heatmap",
  annotation_legend = TRUE
)

ggsave('results/for_figures/corr_analysis/prot_prot_cor.png', p_ann, width=10, height=8)
ggsave('results/corr_analysis/prot_prot_cor.png', p_ann, width=10, height=8)



##### SAME FIRST TWO, BUT TF LEVEL
rm(list = ls())
# Read in data
prot_data <- read.csv('processed_data/filtered_prot_data.csv', check.names = FALSE)
metadata <- read.csv('processed_data/filtered_prot_metadata.csv')
dap_df <- read.csv('results/diff_analysis/DAPs_all.csv')
reg_summary <- read.csv('results/diff_analysis/reg_summary_prot_rnai.csv')

# Get best RNAi per TF (highest number of DAMs)
best_rnai_per_tf <- reg_summary %>%
  group_by(target) %>%
  slice_max(order_by = sig, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(target, RNAi = RNAi)
# Add controls to best rnai to keep
control <- data.frame(target='control', RNAi='w1118_F')
best_rnai_per_tf <- rbind(best_rnai_per_tf, control)

# Filter data to only keep the best rnais
rnai_to_keep <- best_rnai_per_tf$RNAi
prot_data <- prot_data %>%
  filter(RNAi %in% rnai_to_keep) %>%
  left_join(best_rnai_per_tf, by='RNAi') %>%
  select(-RNAi)

# Average protein data by TF
averaged_data <- prot_data %>%
  group_by(target) %>%
  summarise(across(-Sample, mean, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("target") %>%
  as.matrix()

## UNCOMMENT TO FILTER DOWN TO METABOLISM CLASS
# prot_path <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
# proteins_metab_class <- prot_path %>% filter(grepl('Metabolism', classes)) %>% pull(protein) %>% unique()
# averaged_data <- averaged_data %>% as.data.frame() %>% select(any_of(proteins_metab_class)) %>% as.matrix()

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

k <- 6 # 5 for metabolism only
clusters2 <- cutree(row_tree2, k = k)

# Convert to dataframe
cluster_df2 <- data.frame(
  rnai = names(clusters2),
  cluster = clusters2
)
cluster_df2 <- cluster_df2 %>% arrange(cluster)
write.csv(cluster_df2, 'results/corr_analysis/clust_mem/cluster_mem_tf_prot_k6.csv', row.names = FALSE)

row_clusters <- cutree(row_clust, k = k)

cluster_table <- data.frame(
  RNAi_label = names(row_clusters),
  Cluster = row_clusters
)

# Build dendrogram
par(mar = c(8, 4, 4, 2))  # c(bottom, left, top, right) - increased bottom from ~5 to 8
dend_rnai <- plot(row_dend, main = "TF Clustering Dendrogram")

# Save
png(filename = "results/corr_analysis/tf_to_tf/tf_tf_clust_dendro_prot.png", width = 20,
    height = 10, units = "in", res = 800)
par(mar = c(8, 4, 4, 2))  # bottom, left, top, right
plot(row_dend, main = "TF Clustering Dendrogram")
dev.off()

# Get cluster annotation
prot_annotation <- cluster_df2[, 1:2] %>% distinct() %>% select(-rnai) %>% mutate(cluster = as.factor(cluster))

# Define annotation colors for ALL TF (keep full colors)
annotation_colors <- list(cluster = RColorBrewer::brewer.pal(n = length(unique(prot_annotation$cluster)), name = "Set3"))

# Name the colors
names(annotation_colors$cluster) <- unique(prot_annotation$cluster)

## Plot
p2_save <- pheatmap::pheatmap(
  tf_cor,
  annotation_col = prot_annotation,
  annotation_row = prot_annotation,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  labels_row = rownames(tf_cor),          # Custom row labels (only controls)
  labels_col = colnames(tf_cor),          # Custom column labels (only controls)
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 6,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "TF-to-TF Correlation Heatmap (Protein Data)"
)

ggsave('results/corr_analysis/tf_to_tf/tf_tf_cor_prot.png', plot = p2_save, width=10, height=8)
ggsave('results/for_figures/corr_analysis/tf_tf_cor_prot.png', plot = p2_save, width=10, height=8)


# TF-protein Correlation Analysis 
# Get -log10 of adjusted p values, sign correlates with up/downregulation
dap_df <- dap_df %>%
  mutate(target = ifelse(RNAi == 'attP40_F', 'attP40', 
                         ifelse(RNAi == 'attP2_F', 'attP2', target))) %>%
  mutate(signed_log10padj = sign(log2FC) * -log10(p.adj)) %>%
  filter(RNAi %in% rnai_to_keep) %>%
  select(-RNAi)
padj_mat <- dap_df %>%
  select(target, protein, signed_log10padj) %>%
  pivot_wider(names_from = protein,
              values_from = signed_log10padj) %>%
  column_to_rownames('target') %>%
  as.matrix()

# Handle NAs for plotting
m <- t(padj_mat)
m[is.na(m)] <- 0
dist_rows <- dist(m)
dist_cols <- dist(t(m))

# Plot
library(ComplexHeatmap)
Heatmap(
  t(padj_mat),
  name = "Adjusted p-value (NAs in gray)",
  col = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
  na_col = "gray90",
  column_title = "TF-to-Protein Correlation Heatmap",
  cluster_rows = hclust(dist_rows),
  cluster_columns = hclust(dist_cols),
  show_row_names = FALSE,
  show_column_names = TRUE
)
# Save
png(filename = "results/for_figures/corr_analysis/tf_prot_padj_cor.png", width = 10,
    height = 8, units = "in", res = 400)
par(mar = c(8, 4, 4, 2))  # bottom, left, top, right
Heatmap(
  t(padj_mat),
  name = " ",
  col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
  na_col = "gray90",
  column_title = "TF-to-Protein Correlation Heatmap",
  cluster_rows = hclust(dist_rows),
  cluster_columns = hclust(dist_cols),
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 6),
)
dev.off()


