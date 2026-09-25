### Differential analysis on preprocessed proteomics data 
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggrepel)
library(purrr)
library(limma)
library(glue)
library(pheatmap)
library(RColorBrewer)

# Read in data
prot_data <- read.csv("processed_data/preprocessed_prot_data.csv", check.names = FALSE)
rownames(prot_data) <- prot_data$Sample
prot_data$plate <- NULL
metadata <- read.csv('processed_data/preprocessed_prot_metadata.csv')
rownames(metadata) <- metadata$Sample


### Differential Analysis ###
prot_cols <- setdiff(colnames(prot_data), c("Sample", "RNAi")) # Get protein column names

# Average replicates for each RNAi and control
rnai_means <- prot_data %>%
  group_by(RNAi) %>%
  summarise(across(all_of(prot_cols), \(x) mean(x, na.rm = TRUE)))
control <- rnai_means %>%
  filter(RNAi %in% "w1118_F") %>%
  select(-RNAi) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

# Calculate log2FC with the averaged data per RNAi
log2fc_df <- rnai_means %>%
  filter(!RNAi %in% "w1118_F") %>%
  mutate(across(
    all_of(prot_cols),
    ~ .x - control[[cur_column()]]
  ))

# Function to run Welch's t-test for one RNAi
ttest_rnai <- function(data, rnai) {
  # Get RNAi and control data
  rnai_data <- data %>% filter(RNAi == rnai)
  control_data <- data %>% filter(RNAi %in% "w1118_F")
  
  # Map data frame across rows
  map_dfr(prot_cols, function(prot) {
    
    x <- rnai_data[[prot]]
    y <- control_data[[prot]]
    
    vx <- var(x, na.rm = TRUE)
    vy <- var(y, na.rm = TRUE)
    
    # Check for invalid variance
    if (is.na(vx) || is.na(vy) || vx == 0 || vy == 0) {
      return(tibble(RNAi = rnai,
                    protein = prot,
                    p.value = NA_real_))
    }
    
    t_res <- t.test(x, y, var.equal=F)
    
    tibble(RNAi = rnai,
           protein = prot,
           p.value = t_res$p.value)
  })
}

# Apply function to each RNAi
pval_df <- map_dfr(setdiff(unique(prot_data$RNAi), "w1118_F"),
                   ~ ttest_rnai(prot_data, .x))

# Adjust p values with BH
pval_df <- pval_df %>%
  group_by(RNAi) %>% mutate(p.adj = p.adjust(p.value, method='BH')) %>%
  ungroup()

# Reshape df
log2fc_long <- log2fc_df %>%
  pivot_longer(cols = all_of(prot_cols),
               names_to = 'protein',
               values_to = 'log2FC')

# Get target per RNAi
metadata_unique <- metadata %>%
  distinct(RNAi, target)

# Add TF data
results_df <- log2fc_long %>%
  left_join(pval_df, by = c("RNAi", "protein")) %>%
  left_join(metadata_unique, by = "RNAi") %>%
  arrange(p.adj)
results_df$Reg <- ifelse(results_df$log2FC > 0.6 & results_df$p.adj < 0.05, 'Upregulated', 
                         ifelse(results_df$log2FC < -0.6 & results_df$p.adj < 0.05, 'Downregulated', 'Not significant'))
# Discovery version with p < 0.10
results_df_disc <- log2fc_long %>%
  left_join(pval_df, by = c("RNAi", "protein")) %>%
  left_join(metadata_unique, by = "RNAi") %>%
  arrange(p.adj)
results_df_disc$Reg <- ifelse(results_df_disc$log2FC > 0.6 & results_df_disc$p.adj < 0.1, 'Upregulated', 
                         ifelse(results_df_disc$log2FC < -0.6 & results_df_disc$p.adj < 0.1, 'Downregulated', 'Not significant'))

results_df_sig <- results_df %>% filter(!Reg %in% c('Not significant', NA))

# Get number of up/down/ns per protein, RNAi, and TF
regulation_summary_tf <- results_df %>%
  dplyr::count(RNAi, Reg) %>%
  tidyr::pivot_wider(
    names_from = Reg,
    values_from = n,
    values_fill = 0) 
regulation_summary_tf$sig <- regulation_summary_tf$Downregulated + regulation_summary_tf$Upregulated
regulation_summary_tf <- regulation_summary_tf %>%
  left_join(metadata_unique, by = "RNAi")

# Histogram of log2FC
log_hist <- ggplot(regulation_summary_tf, aes(x=sig)) +
  geom_histogram(bins=30, fill = 'gray', color='black') +
  geom_vline(xintercept = 15, col = "firebrick", linetype = 'dashed') +
  xlab('Number of DAPs') +
  ylab('Frequency') +
  ggtitle('Count of DAPs per RNAi') +
  theme_minimal()
log_hist
ggsave('results/diff_analysis/plots_prot/log2fc_histogram_prot.png', height = 6, width = 8)
ggsave('results/for_figures/diff_analysis/log2fc_histogram_prot.png', height = 4, width = 6)

# Filter regulation summary and dataframes down to responsive TFs
# Set threshold to RNAi with > 15 significant proteins; leaves 133/144 of RNAi lines + control
threshold <- 15
rnai_to_keep <- regulation_summary_tf$RNAi[regulation_summary_tf$sig > threshold]
# Add controls to RNAi to keep
controls_keep <- c('w1118_F', 'attp40_F', 'attP2_F')
rnai_to_keep <- c(rnai_to_keep, controls_keep)

# Filter down data and metadata
prot_data_resp <- prot_data %>% filter(prot_data$RNAi %in% rnai_to_keep)
write.csv(prot_data_resp, 'processed_data/filtered_prot_data.csv', row.names = FALSE)
metadata_resp <- metadata %>% filter(metadata$RNAi %in% rnai_to_keep)
write.csv(metadata_resp, 'processed_data/filtered_prot_metadata.csv', row.names = FALSE)

# Filter down DAPs
results_df <- results_df %>% filter(RNAi %in% rnai_to_keep)
results_df_disc <- results_df_disc %>% filter(RNAi %in% rnai_to_keep)
results_df_sig <- results_df_sig %>% filter(RNAi %in% rnai_to_keep)
write.csv(results_df, 'results/diff_analysis/DAPs_all.csv', row.names = FALSE)
write.csv(results_df_disc, 'results/diff_analysis/DAPs_all_discovery.csv', row.names = FALSE)
write.csv(results_df_sig, 'results/diff_analysis/DAPs_sig.csv', row.names = FALSE)

# Filtered down reg summaries
regulation_summary <- results_df %>%
  dplyr::count(protein, Reg) %>%
  tidyr::pivot_wider(names_from = Reg,
                     values_from = n,
                     values_fill = 0) %>%
  mutate(sig = Downregulated + Upregulated)

regulation_summary_disc <- results_df_disc %>%
  group_by(protein) %>%
  summarise(
    Downregulated = sum(Reg == "Downregulated" & !is.na(p.adj)),
    Upregulated = sum(Reg == "Upregulated" & !is.na(p.adj)),
    Not.significant = sum(Reg == "Not significant" & !is.na(p.adj)),
    NAs = sum(is.na(p.adj)),
    sig = Downregulated + Upregulated,
    sig_RNAi = paste(RNAi[Reg %in% c("Upregulated","Downregulated")], collapse=";"),
    .groups = 'drop')

write.csv(regulation_summary, 'results/diff_analysis/reg_summary_prot.csv', row.names = FALSE)
write.csv(regulation_summary_disc, 'results/diff_analysis/reg_summary_prot_disc.csv', row.names = FALSE)

regulation_summary_tf <- regulation_summary_tf %>%
  filter(regulation_summary_tf$RNAi %in% rnai_to_keep)
write.csv(regulation_summary_tf, 'results/diff_analysis/reg_summary_prot_rnai.csv', row.names = FALSE)


# Add fractions of upregulated and downregulated genes
regulation_summary_tf$frac_up <- regulation_summary_tf$Upregulated / 
  (regulation_summary_tf$sig + regulation_summary_tf$`Not significant`)
regulation_summary_tf$frac_down <- regulation_summary_tf$Downregulated / 
  (regulation_summary_tf$sig + regulation_summary_tf$`Not significant`)
# Plot
p1 <- ggplot(regulation_summary_tf, aes(x=frac_up)) +
  geom_histogram(bins=30, fill = 'gray', color='black') +
  xlab('Fraction of Upregulated DAPs') +
  ylab('Number of Perturbations') +
  ggtitle('Fractions of Upregulated DAPs') +
  theme_minimal()
p2 <- ggplot(regulation_summary_tf, aes(x=frac_down)) +
  geom_histogram(bins=15, fill = 'gray', color='black') +
  xlab('Fraction of Downregulated DAPs') +
  ylab('Number of Perturbations') +
  coord_cartesian(xlim = c(0, 0.15)) +
  ggtitle('Fractions of Downregulated DAPs') +
  theme_minimal()
fractions <- p1 + p2
fractions
ggsave('results/diff_analysis/plots_prot/reg_fractions_prot.png', height = 5, width = 8)

regulation_summary_prot <- results_df_sig %>%
  count(protein, Reg) %>%
  tidyr::pivot_wider(
    names_from = Reg,
    values_from = n,
    values_fill = 0) 
regulation_summary_prot$sig <- regulation_summary_prot$Downregulated + regulation_summary_prot$Upregulated
# Plot
p3 <- ggplot(regulation_summary_prot, aes(x=sig)) +
  geom_histogram(bins=30, fill = 'gray', color='black') +
  xlab('Number of TFs Affecting a Protein') +
  ylab('Frequency') +
  ggtitle('Histogram of Number of TFs Affecting Proteins') +
  theme_minimal()
ggsave('results/diff_analysis/plots_prot/tfs_affecting_prot.png', p3, height = 6, width = 8)
ggsave('results/for_figures/diff_analysis/tfs_affecting_prot.png', p3, height = 4, width = 6)

### Volcano plot
# Get top 20 DAM labels
results_df <- results_df %>% arrange(p.adj)
# Plot
volcano_plot <- ggplot(data = results_df, aes(x=log2FC, y=-log10(p.adj), col=Reg)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(xlim = c(-8, 10)) +
  labs(color = 'Regulation',
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  ggtitle('Differentially Abundant Proteins')
volcano_plot
ggsave('results/for_figures/diff_analysis/volcanoplot.png', width = 6, height = 4)


### RNAi-RNAi correlation within targets V2
expr_mat <- prot_data_resp %>% select(-c(RNAi, Sample)) %>% as.matrix()
expr_mat <- expr_mat[metadata_resp$Sample, ]

# Get control mean
metadata_resp$target[is.na(metadata_resp$target)] <- metadata_resp$Genotype[is.na(metadata_resp$target)]
ctrl_samps <- metadata_resp$Sample[metadata_resp$target == 'control']
ctrl_means <- colMeans(expr_mat[ctrl_samps, ])

# Center matrix to controls
expr_mat_centered <- sweep(expr_mat, 2, ctrl_means, '-')

# Add RNAi and target information back
expr_df <- as.data.frame(expr_mat_centered) %>% 
  tibble::rownames_to_column('Sample') %>%
  left_join(metadata_resp %>% select(Sample, RNAi, target), by = 'Sample')

# Group by target and calculate correlation
tf_cor_list <- expr_df %>%
  filter(target != "w1118") %>%
  group_split(target)

tf_cor_results <- map_df(tf_cor_list, function(df_target) {
  reps <- df_target$Sample
  target <- unique(df_target$target)
  
  mat_subset <- expr_mat_centered[reps, , drop = FALSE]
  cor_mat <- cor(t(mat_subset)) 
  
  cor_vals <- cor_mat[upper.tri(cor_mat)]
  
  data.frame(Target = target,
             Correlation = mean(cor_vals))
})

# Full correlation matrix
full_cor_mat <- cor(t(expr_mat_centered))

# Get all pair combinations
rnai_info <- expr_df[, c("RNAi", "target")]

bg_values <- c()

# Get the background mean and standard deviation from all other pairs
for (i in 1:(nrow(full_cor_mat)-1)) {
  for (j in (i+1):nrow(full_cor_mat)) {
    target_i <- rnai_info$target[i]
    target_j <- rnai_info$target[j]
    
    # Calculate background corr for every pair where the target is not the same
    if (target_i != target_j) {
      bg_values <- c(bg_values, full_cor_mat[i, j])
    }
  }
}

bg_mean <- mean(bg_values) # 0.22
bg_sd <- sd(bg_values) # 0.13

# Calculating z-score and p-value
tf_cor_results <- tf_cor_results %>%
  mutate(z.score = (Correlation - bg_mean) / bg_sd,
         p.value = 1 - pnorm(z.score),
         p.adj = p.adjust(p.value, method = 'BH')) %>% 
  arrange(p.adj)
write.csv(tf_cor_results, 'results/diff_analysis/null_dist/tf_null_dist_prot.csv', row.names = FALSE)
