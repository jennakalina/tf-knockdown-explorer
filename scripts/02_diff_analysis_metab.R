### Differential analysis on preprocessed metabolomics data 
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggrepel)
library(purrr)
library(limma)
library(glue)

# Read in data
metab_data <- read.csv("processed_data/preprocessed_metab_data.csv", check.names = FALSE)
metadata <- read.csv('processed_data/preprocessed_metab_metadata.csv')

### Differential Analysis ###
metab_cols <- setdiff(colnames(metab_data), c("Sample", "RNAi")) # Get metabolite column names

# Average replicates for each RNAi and control
rnai_means <- metab_data %>%
  group_by(RNAi) %>%
  summarise(across(all_of(metab_cols), mean, na.rm=T))
control <- rnai_means %>%
  filter(RNAi %in% "w1118_F") %>%
  select(-RNAi) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Calculate log2FC with the averaged data per RNAi
log2fc_df <- rnai_means %>%
  filter(!RNAi %in% "w1118_F") %>%
  mutate(across(all_of(metab_cols), ~ .x - control[[cur_column()]]))

# Function to run Welch's t-test for one RNAi
ttest_rnai <- function(data, rnai) {
  # Get RNAi and control data
  rnai_data <- data %>% filter(RNAi == rnai)
  control_data <- data %>% 
    filter(RNAi %in% "w1118_F")
  
  # Map data frame across rows
  map_dfr(metab_cols, function(metab) {
    t_res <- t.test(rnai_data[[metab]],
                    control_data[[metab]],
                    var.equal=F)
    tibble(RNAi = rnai,
           metabolite = metab,
           p.value = t_res$p.value)
  })
}

# Apply function to each RNAi
pval_df <- map_dfr(setdiff(unique(metab_data$RNAi), "w1118_F"),
                   ~ ttest_rnai(metab_data, .x))

# Adjust p values with BH
pval_df <- pval_df %>%
  group_by(RNAi) %>% mutate(p.adj = p.adjust(p.value, method='BH')) %>%
  ungroup()

# Reshape df
log2fc_long <- log2fc_df %>%
  pivot_longer(cols = all_of(metab_cols),
               names_to = 'metabolite',
               values_to = 'log2FC')

# Get target per RNAi
metadata_unique <- metadata %>%
  distinct(RNAi, Target)

# Add TF data
results_df <- log2fc_long %>%
  left_join(pval_df, by = c("RNAi", "metabolite")) %>%
  left_join(metadata_unique, by = "RNAi")
results_df$Reg <- ifelse(results_df$log2FC > 0.6 & results_df$p.adj < 0.05, 'Upregulated', 
                         ifelse(results_df$log2FC < -0.6 & results_df$p.adj < 0.05, 'Downregulated', 'Not significant'))

# Make discovery version
results_df_disc <- log2fc_long %>%
  left_join(pval_df, by = c("RNAi", "metabolite")) %>%
  left_join(metadata_unique, by = "RNAi")
results_df_disc$Reg <- ifelse(results_df_disc$log2FC > 0.6 & results_df_disc$p.adj < 0.1, 'Upregulated', 
                         ifelse(results_df_disc$log2FC < -0.6 & results_df_disc$p.adj < 0.1, 'Downregulated', 'Not significant'))

# Get number of up/down/ns per RNAi
regulation_summary <- results_df %>%
  count(RNAi, Reg) %>%
  tidyr::pivot_wider(
    names_from = Reg,
    values_from = n,
    values_fill = 0) 
regulation_summary$sig <- regulation_summary$Downregulated + regulation_summary$Upregulated

# Histogram of log2FC
log_hist <- ggplot(regulation_summary, aes(x=sig)) +
  geom_histogram(bins=30, fill = 'gray', color='black') +
  geom_vline(xintercept = 2, col = "firebrick", linetype = 'dashed') +
  xlab('Number of DAMs') +
  ylab('Frequency') +
  ggtitle('Count of DAMs per RNAi') +
  theme_minimal()
log_hist
ggsave('results/diff_analysis/log2fc_histogram_metab.png', width = 8, height = 6)
ggsave('results/for_figures/diff_analysis/log2fc_histogram_metab.png', width = 6, height = 4)

# Filter regulation summary and dataframes down to responsive TFs
# Set threshold to RNAi with at least 3 significant DAMs
threshold <- 3
rnai_to_keep <- regulation_summary$RNAi[regulation_summary$sig >= threshold]
# Add controls to RNAi to keep
controls_keep <- c('w1118_F', 'attP2_F', 'attP40_F')
rnai_to_keep <- c(rnai_to_keep, controls_keep)

metab_data_resp <- metab_data %>%
  filter(metab_data$RNAi %in% rnai_to_keep)
write.csv(metab_data_resp, 'processed_data/filtered_metab_data.csv', row.names = FALSE)
metadata_resp <- metadata %>%
  filter(metadata$RNAi %in% rnai_to_keep)
write.csv(metadata_resp, 'processed_data/filtered_metab_metadata.csv', row.names = FALSE)
regulation_summary <- regulation_summary %>%
  filter(regulation_summary$RNAi %in% rnai_to_keep)

# Add fractions of upregulated and downregulated genes
regulation_summary$frac_up <- regulation_summary$Upregulated / 
  (regulation_summary$sig + regulation_summary$`Not significant`)
regulation_summary$frac_down <- regulation_summary$Downregulated / 
  (regulation_summary$sig + regulation_summary$`Not significant`)

# Plot
p1 <- ggplot(regulation_summary, aes(x=frac_up)) +
  geom_histogram(bins=30, fill = 'gray', color='black') +
  xlab('Fraction of Upregulated DAMs') +
  ylab('Number of Perturbations') +
  ggtitle('Histogram of Fractions of Upregulated DAMs') +
  theme_minimal()
p2 <- ggplot(regulation_summary, aes(x=frac_down)) +
  geom_histogram(bins=15, fill = 'gray', color='black') +
  xlab('Fraction of Downregulated DAMs') +
  ylab('Number of Perturbations') +
  coord_cartesian(xlim = c(0, 0.4)) +
  ggtitle('Histogram of Fractions of Downregulated DAMs') +
  theme_minimal()

fractions <- p1 + p2
fractions
ggsave('results/diff_analysis/reg_fractions_metab.png', width = 12, height = 6)

# For each metabolite, how many TFs affect it?
# Number of up/down/ns by metabolite
results_df_filtered <- results_df %>%
  filter(results_df$RNAi %in% rnai_to_keep) %>%
  arrange(p.adj)
results_df_filtered_disc <- results_df_disc %>%
  filter(results_df_disc$RNAi %in% rnai_to_keep) %>%
  arrange(p.adj)
write.csv(results_df_filtered, 'results/diff_analysis/DAMs_all.csv', row.names = FALSE)
write.csv(results_df_filtered_disc, 'results/diff_analysis/DAMs_all_discovery.csv', row.names = FALSE)
results_df_sig <- results_df_filtered %>%
  filter(abs(log2FC) > 0.6) %>%
  filter(p.adj < 0.05)%>%
  arrange(p.adj)
write.csv(results_df_sig, 'results/diff_analysis/DAMs_sig.csv', row.names = FALSE)

reg_summary_metab <- results_df_filtered %>%
  count(metabolite, Reg) %>%
  tidyr::pivot_wider(
    names_from = Reg,
    values_from = n,
    values_fill = 0)
reg_summary_metab$sig <- reg_summary_metab$Downregulated + reg_summary_metab$Upregulated
write.csv(reg_summary_metab, 'results/diff_analysis/reg_summary_metab.csv', row.names = FALSE)
regulation_summary_disc <- results_df_filtered_disc %>%
  group_by(metabolite) %>%
  summarise(
    Downregulated = sum(Reg == "Downregulated"),
    Upregulated = sum(Reg == "Upregulated"),
    Not.significant = sum(Reg == "Not significant"),
    sig = Downregulated + Upregulated,
    sig_RNAi = paste(RNAi[Reg %in% c("Upregulated","Downregulated")], collapse=";"),
    .groups = 'drop')
write.csv(regulation_summary_disc, 'results/diff_analysis/reg_summary_metab_disc.csv', row.names = FALSE)

regulation_summary_tf <- results_df_filtered %>%
  count(Target, Reg) %>%
  tidyr::pivot_wider(
    names_from = Reg,
    values_from = n,
    values_fill = 0) 
regulation_summary_tf$sig <- regulation_summary_tf$Downregulated + regulation_summary_tf$Upregulated
write.csv(regulation_summary_tf, 'results/diff_analysis/reg_summary_metab_tf.csv', row.names = FALSE)

# Plot
p3 <- ggplot(reg_summary_metab, aes(x=sig)) +
  geom_histogram(bins=30, fill = 'gray', color='black') +
  xlab('Number of TFs Affecting a Metabolite') +
  ylab('Frequency') +
  ggtitle('Histogram of Number of TFs Affecting Metabolites') +
  theme_minimal()
p3
ggsave('results/diff_analysis/tfs_affecting_metab.png', width = 10, height = 6)
ggsave('results/for_figures/diff_analysis/tfs_affecting_metab.png', width = 6, height = 4)

# Volcano plot
# Get top 20 DAM labels
df_ordered <- results_df %>% arrange(p.adj)
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
  ggtitle('Differentially Abundant Metabolites') 
volcano_plot
ggsave('results/diff_analysis/volcano_metab.png', width = 10, height = 10)
ggsave('results/for_figures/diff_analysis/volcano_metab.png', width = 6, height = 4)

# Add target information to regulation summary
regulation_summary <- regulation_summary %>%
  left_join(metadata_unique, by = "RNAi")
write.csv(regulation_summary, 'results/diff_analysis/reg_summary_metab_rnai.csv', row.names = FALSE)


### RNAi-RNAi correlation within targets V2
expr_mat <- metab_data_resp %>% column_to_rownames('Sample') %>% select(-RNAi) %>% as.matrix()
expr_mat <- expr_mat[metadata_resp$Sample, ]

# Get control mean
ctrl_samps <- metadata_resp$Sample[metadata_resp$Target == 'w1118']
ctrl_means <- colMeans(expr_mat[ctrl_samps, ])

# Center matrix to controls
expr_mat_centered <- sweep(expr_mat, 2, ctrl_means, '-')

# Add RNAi and target information back
expr_df <- as.data.frame(expr_mat_centered) %>% 
  tibble::rownames_to_column('Sample') %>%
  left_join(metadata_resp %>% select(Sample, RNAi, Target), by = 'Sample')

# Group by target and calculate correlation
tf_cor_list <- expr_df %>%
  filter(Target != "w1118") %>%
  group_split(Target)

tf_cor_results <- map_df(tf_cor_list, function(df_target) {
  reps <- df_target$Sample
  target <- unique(df_target$Target)
  
  mat_subset <- expr_mat_centered[reps, , drop = FALSE]
  cor_mat <- cor(t(mat_subset))   # correlation across metabolites
  
  cor_vals <- cor_mat[upper.tri(cor_mat)]
  
  data.frame(Target = target,
             Correlation = mean(cor_vals))
})

# Full correlation matrix
full_cor_mat <- cor(t(expr_mat_centered))

# Get all pair combinations
rnai_info <- expr_df[, c("RNAi", "Target")]

bg_values <- c()

# Get the background mean and standard deviation from all other pairs
for (i in 1:(nrow(full_cor_mat)-1)) {
  for (j in (i+1):nrow(full_cor_mat)) {
    target_i <- rnai_info$Target[i]
    target_j <- rnai_info$Target[j]
    
    # Calculate background corr for every pair where the target is not the same
    if (target_i != target_j) {
      bg_values <- c(bg_values, full_cor_mat[i, j])
    }
  }
}

bg_mean <- mean(bg_values) # 0.39
bg_sd <- sd(bg_values) # 0.14

# Calculating z-score and p-value
tf_cor_results <- tf_cor_results %>%
  mutate(z.score = (Correlation - bg_mean) / bg_sd,
         p.value = 1 - pnorm(z.score),
         p.adj = p.adjust(p.value, method = 'BH')) %>%
  arrange(p.adj)
write.csv(tf_cor_results, 'results/diff_analysis/null_dist/tf_null_dist_metab.csv', row.names = FALSE)
