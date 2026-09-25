### Filtering and preprocessing of proteomics data

# Remove unwanted samples, filter by missingness, impute NAs, 
# normalize, scale, remove bad replicates
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggrepel)
library(readxl)
library(purrr)
library(patchwork)
library(umap)

################## READ IN RAW DATA ##################

prot_data_raw <- read_excel('raw_data/TF_screen_proteomics.xlsx')
metadata <- prot_data_raw %>% select(c(1:19))
prot_data <- prot_data_raw %>% select(c(1, 20:1229)) %>%
  column_to_rownames('Sample')

################## CLEAN UP METADATA ##################

metadata <- metadata %>%
  dplyr::rename(RNAi = Unique_Identifier, target = Gene, control = `Control to use2`, batch = `Fly B-day (month-year)`, fly_count = `fly number`) %>%
  mutate(plate = as.factor(plate), batch = as.factor(batch), fly_count = as.factor(fly_count)) %>%
  as.data.frame()
rownames(metadata) <- metadata$Sample
metadata$plate <- substr(metadata$plate, 6, 6)

### CHECK FOR DUPLICATE PROTEINS (none found)
# prot_cols <- names(prot_data)[!names(prot_data) %in% c("Sample", "plate")]
# duplicated_names <- unique(prot_cols[duplicated(prot_cols)]) 

################## FILTERING PROTEINS ##################

# Get number of zero values per protein
col_zeroes <- colSums(prot_data == 0) %>% 
  as.data.frame() %>% 
  rownames_to_column('Protein') %>% 
  dplyr::rename(Zeroes = '.')

# Explore data for batch effects by plate
prot_data_plate <- prot_data %>% 
  rownames_to_column('Sample') %>%
  left_join(y = metadata %>% select(c('Sample', 'plate')), by = 'Sample')
prot_data_plate <- prot_data_plate %>% group_by(plate)

# Compute min, mean, sd, and max for each plate
qc_plate <- prot_data_plate %>%
  summarise(mean = mean(unlist(across(where(is.numeric)))),
            sd = sd(unlist(across(where(is.numeric)))),
            min = min(unlist(across(where(is.numeric)))),
            max = max(unlist(across(where(is.numeric)))))

# And for batch
prot_data_batch <- prot_data %>% 
  rownames_to_column('Sample') %>%
  left_join(y = metadata %>% select(c('Sample', 'batch')), by = 'Sample')
prot_data_batch <- prot_data_batch %>% group_by(batch)

# Compute min, mean, sd, and max for each batch
qc_batch <- prot_data_batch %>%
  summarise(mean = mean(unlist(across(where(is.numeric)))),
            sd = sd(unlist(across(where(is.numeric)))),
            min = min(unlist(across(where(is.numeric)))),
            max = max(unlist(across(where(is.numeric)))))

### FILTER for missingness

threshold <- .5 # Change when deciding actual threshold; currently 50% missing
# 506 and not 497 because even if we're not using those lines they still had missing or nonmissing data
prot_filter <- col_zeroes$Protein[(col_zeroes$Zeroes / 506) > threshold] 
prot_data_filtered <- prot_data[, !(colnames(prot_data) %in% prot_filter)]

### IMPUTE with plate half-minimum per protein 
prot_filtered_plate <- prot_data_filtered %>%
  rownames_to_column('Sample') %>%
  left_join(y = metadata %>% select(c('Sample', 'plate')), by = 'Sample')
prot_filtered_plate <- prot_filtered_plate %>% group_by(plate)

prot_imputed <- prot_filtered_plate %>% 
  pivot_longer(cols = where(is.numeric),
               names_to = 'protein',
               values_to = 'intensity') %>%
  group_by(plate, protein) %>%
  mutate(nonzero_min = if(any(intensity > 1e4, na.rm = TRUE)) {
    min(intensity[intensity > 1e4], na.rm = TRUE)} else {0},
    impute_val = nonzero_min / 2,
    intensity = ifelse(intensity < 1e4, impute_val, intensity)) %>%
  ungroup() 

# Compute mean per protein per plate
protein_plate_means <- prot_imputed %>%
  group_by(protein, plate) %>%
  summarise(plate_mean = mean(intensity, na.rm = TRUE), .groups = "drop")

plot_df <- prot_imputed %>%
  left_join(protein_plate_means %>% select(protein, plate), by = c("protein", "plate")) %>%
  select(Sample, plate, protein, intensity)

prot_imputed <- prot_imputed %>%
  left_join(protein_plate_means %>% select(protein, plate), by = c("protein", "plate")) %>%
  select(Sample, plate, protein, intensity) %>%
  pivot_wider(names_from = protein, values_from = intensity)

# Check intensities
sample_order <- prot_imputed %>%
  distinct(Sample, plate) %>%
  arrange(plate, Sample) %>%
  pull(Sample)
plot_df <- plot_df %>%
  mutate(Sample = factor(Sample, levels = sample_order))
plate_annot <- prot_imputed %>%
  distinct(Sample, plate) %>%
  mutate(y = -500000) %>%
  mutate(Sample = factor(Sample, levels = sample_order))

ggplot(plot_df, aes(x = Sample, y = intensity, fill = plate)) +
  geom_boxplot(outlier.size = 0, width = 0.7) +
  geom_tile(data = plate_annot, mapping = aes(x = Sample, y = y, fill = plate),
            height = 500000, inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Distribution of Protein Intensities by Sample",
       x = "Sample", y = "Abundance", fill = "Plate") +
  ylim(c(-1000000, max(plot_df$intensity))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), # Hide x labels (too many)
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

#ggsave('results/for_figures/preprocessing/intensities_preproc_prot.png', height = 4, width = 12)


# Log2 transform
prot_imputed[c(3:ncol(prot_imputed))] <- log2(prot_imputed[c(3:ncol(prot_imputed))])

# Scale per plate using log2 values
prot_imputed <- prot_imputed %>% 
  pivot_longer(cols = where(is.numeric),
               names_to = 'protein',
               values_to = 'intensity') %>%
  group_by(plate, protein) %>%
  mutate(plate_mean = mean(intensity), plate_sd = sd(intensity),
         intensity_scaled = ifelse(plate_sd > 0,
                                   (intensity - plate_mean) / plate_sd, 0)) %>%
  ungroup() %>%
  select(Sample, plate, protein, intensity_scaled) %>%
  pivot_wider(names_from = protein, values_from = intensity_scaled)

### FILTER out bad replicates 

prot_cols <- setdiff(colnames(prot_imputed), c("Sample", "plate", "RNAi"))

# Function to check correlations between RNAi, find those over 3 with the least correlation
corr_check <- function(rnai) {
  rnai_name <- unique(rnai$RNAi)
  
  # Skip if control
  if (rnai_name == "w1118_F") {
    return(invisible(NULL))
  }
  
  # Skip if less than 3 replicates
  n_reps <- nrow(rnai)
  if (n_reps < 3) {
    return(invisible(NULL))
  }
  
  # Get correlation matrix
  prot_mat <- as.matrix(rnai[, prot_cols])
  rownames(prot_mat) <- rnai$Sample
  corr_mat <- cor(t(prot_mat), use = 'pairwise.complete.obs')
  
  # Calculate each replicate's average correlation to others
  mean_cor <- apply(corr_mat, 1, function(x) mean(x[-which.max(x==1)], na.rm=TRUE))
  
  # Names of low correlation replicates
  low_cor_reps <- names(mean_cor[mean_cor < 0.5])
  
  # If there are bad replicates, return name of RNAi
  if (length(low_cor_reps > 0)) {
    return(rnai_name)
  }
}

# Get names of RNAIs with bad reps
prot_imputed <- prot_imputed %>% left_join(metadata %>% select(Sample, RNAi))

rnais_bad_reps <- prot_imputed %>%
  group_by(RNAi) %>%
  group_split() %>%
  map(corr_check) %>%
  compact() %>%
  unlist

# Function to eliminate low correlation reps with backward elimination
remove_bad_reps <- function(norm_data, cutoff = 0.5) {
  reps_to_remove <- list() 
  
  # Separate RNAi into groups
  for (rnai_name in rnais_bad_reps) {
    rnai <- subset(norm_data, RNAi == rnai_name)
    removed <- character(0)
    
    repeat {
      # Build correlation matrix
      prot_mat <- as.matrix(rnai[, prot_cols])
      rownames(prot_mat) <- rnai$Sample
      corr_mat <- cor(t(prot_mat), use = "pairwise.complete.obs")
      
      # Stop if there are only 3 replicates left
      if (nrow(corr_mat) <= 3) break
      
      # Find the lowest correlation; if it is above the cutoff then stop
      min_corr <- min(corr_mat[lower.tri(corr_mat)], na.rm = TRUE)
      if (min_corr >= cutoff) break
      
      # Find pair w lowest correlation
      low_pair <- which(corr_mat == min_corr, arr.ind = TRUE)[1, ]
      pair_samples <- rownames(corr_mat)[low_pair]
      
      # Calculate average correlation for each sample in the pair, remove one with lowest
      avg_corrs <- sapply(pair_samples, function(s) {
        mean(corr_mat[s, setdiff(rownames(corr_mat), s)], na.rm = TRUE)
      })
      to_remove <- names(which.min(avg_corrs))
      removed <- c(removed, to_remove)
      rnai <- rnai[rnai$Sample != to_remove, ]
    }
    
    reps_to_remove[[rnai_name]] <- removed
  }
  return(reps_to_remove)
}

# Find replicates to remove
bad_reps <- remove_bad_reps(prot_imputed)

# Remove replicates
samples_to_remove <- unlist(bad_reps, use.names = FALSE)
prot_imputed_nbr <- prot_imputed %>%
  filter(!Sample %in% samples_to_remove)
metadata_nbr <- metadata %>%
  filter(!Sample %in% samples_to_remove)
metadata_nbr$Sample <- as.character(metadata_nbr$Sample)

# Filter out RNAi that only have 1 replicate
prot_imputed_nbr <- prot_imputed_nbr %>%
  add_count(RNAi) %>%
  filter(n > 1) %>%
  select(-n)
metadata_nbr <- metadata_nbr %>%
  add_count(RNAi) %>%
  filter(n > 1) %>%
  select(-n)

# Filter out male RNAi lines (_M in the name)
male_rnai_samps <- metadata_nbr$Sample[grepl('_M', metadata_nbr$RNAi)]
prot_imputed_nbr <- prot_imputed_nbr %>%
  filter(!Sample %in% male_rnai_samps)
metadata_nbr <- metadata_nbr %>%
  filter(!Sample %in% male_rnai_samps)

# Fix metadata
metadata_nbr <- metadata_nbr %>%
  mutate(target = ifelse(is.na(target), 'control', target),
         target = as.character(target))

# Remove all data with cuff and daw targets, added by accident
cuff_samp <- metadata_nbr$Sample[metadata_nbr$target == 'cuff']
daw_samp <- metadata_nbr$Sample[metadata_nbr$target == 'daw']
jra_samp <- metadata_nbr$Sample[metadata_nbr$RNAi == 'Jra_7216_F']
samp_rm <- c(cuff_samp, daw_samp, jra_samp)
prot_imputed_nbr <- prot_imputed_nbr %>%
  filter(!Sample %in% samp_rm) 
metadata_nbr <- metadata_nbr %>%
  filter(!Sample %in% samp_rm)

# Save filtered data and metadata
#write.csv(prot_imputed_nbr, 'processed_data/preprocessed_prot_data.csv', row.names = FALSE)
#write.csv(metadata_nbr, 'processed_data/preprocessed_prot_metadata.csv', row.names = FALSE)

# Plot intensities again
long_df <- prot_imputed_nbr %>%
  select(-RNAi) %>%
  pivot_longer(-c(Sample, plate), names_to = "Protein", values_to = "Value")

sample_order <- long_df %>%
  distinct(Sample, plate) %>%
  arrange(plate, Sample) %>%
  pull(Sample)
long_df <- long_df %>%
  mutate(Sample = factor(Sample, levels = sample_order))
plate_annot <- long_df %>%
  distinct(Sample, plate) %>%
  mutate(y = -11) %>%
  mutate(Sample = factor(Sample, levels = sample_order))

ggplot(long_df, aes(x = Sample, y = Value, fill = plate)) +
  geom_boxplot(outlier.size = 0, width = 0.7) +
  geom_tile(data = plate_annot, mapping = aes(x = Sample, y = y, fill = plate),
            height = 0.5, inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Distribution of Protein Intensities (log2) by Sample",
       x = "Sample", y = "log2(Abundance)", fill = "Plate") +
  ylim(c(-11.5, max(long_df$Value))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), # Hide x labels (too many)
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

ggsave('results/for_figures/preprocessing/intensities_postproc_prot.png', height = 4, width = 12)

### PCA
dat <- prot_imputed_nbr %>%
  select(-c(Sample, plate, RNAi)) %>% 
  as.matrix()
rownames(dat) <- prot_imputed_nbr$Sample 

# Run PCA
controls <- "w1118_F"

pca <- prcomp(dat, scale. = TRUE, center = TRUE)  # already median-normalized

# Proportion of variance explained
expl_var <- pca$sdev^2 / sum(pca$sdev^2)
pc1_var <- round(100 * expl_var[1], 1)
pc2_var <- round(100 * expl_var[2], 1)

# Extract PCA scores (first 2 PCs)
pca_df <- as.data.frame(pca$x[,1:2]) %>%
  mutate(Sample = rownames(pca$x)) %>%
  left_join(metadata_nbr, by = "Sample")

pca_df_temp <- pca_df %>%
  mutate(RNAi_grouped = ifelse(RNAi %in% controls, RNAi, "Other"),
         label = ifelse(RNAi %in% controls, 'Control', NA))
pca_df_temp$RNAi_grouped <- factor(pca_df_temp$RNAi_grouped,
                                   levels = c('w1118_F', 'Other'))

# Plot, colored by Plate
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = plate)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(
    title = "PCA Colored by Plate",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

# Plot, colored by Plate
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(
    title = "PCA Colored by Batch",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

p3 <- ggplot(pca_df_temp, aes(x = PC1, y = PC2, color = RNAi_grouped)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text_repel(
    aes(label = label), 
    na.rm = TRUE, size = 3, show.legend = FALSE
  ) +
  theme_classic() +
  labs(
    title = "PCA Colored by RNAi (Controls Labeled)",
    color = "RNAi",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

p <- p1 + p3
p
ggsave('results/preprocessing/PCAplots_prot.png', p, width=13, height=5)
ggsave('results/for_figures/preprocessing/PCAplots_prot.png', p, width=8, height=4)

# Get PCA for old data
prot_data <- prot_data_raw %>% select(c(1, 20:1229)) %>%
  column_to_rownames('Sample')
prot_to_remove <- c("A53E_DROME",  "CH16_DROME",  "T2FA_DROME",  "MS57A_DROME", "CP304_DROME")
prot_data <- prot_data %>% select(-prot_to_remove)
dat <- prot_data %>% as.matrix()
pca <- prcomp(dat, scale. = TRUE, center = TRUE)  # already median-normalized

# Proportion of variance explained
expl_var <- pca$sdev^2 / sum(pca$sdev^2)
pc1_var <- round(100 * expl_var[1], 1)
pc2_var <- round(100 * expl_var[2], 1)

metadata <- prot_data_raw %>% select(c(1:19))
metadata <- metadata %>% 
  dplyr::rename(batch = `Fly B-day (month-year)`, RNAi = Unique_Identifier) %>%
  mutate(plate = gsub('Plate', '', plate)) %>%
  mutate(plate = gsub('-1', '', plate))
# Extract PCA scores (first 2 PCs)
pca_df <- as.data.frame(pca$x[,1:2]) %>%
  mutate(Sample = rownames(pca$x)) %>%
  left_join(metadata, by = "Sample")

pca_df_temp <- pca_df %>%
  mutate(RNAi_grouped = ifelse(RNAi %in% controls, RNAi, "Other"),
         label = ifelse(RNAi %in% controls, 'Control', NA))
pca_df_temp$RNAi_grouped <- factor(pca_df_temp$RNAi_grouped,
                                   levels = c('w1118_F', 'Other'))

# Plot, colored by Plate
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = plate)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(
    title = "PCA Colored by Plate",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

# Plot, colored by Plate
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(
    title = "PCA Colored by Batch",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

p3 <- ggplot(pca_df_temp, aes(x = PC1, y = PC2, color = RNAi_grouped)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text_repel(
    aes(label = label), 
    na.rm = TRUE, size = 3, show.legend = FALSE
  ) +
  theme_classic() +
  labs(title = "PCA Colored by RNAi (Controls Labeled)", color = "RNAi", 
       x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)"))

p <- p1 + p3
p
ggsave('results/for_figures/preprocessing/PCAplots_prot_old.png', p, width=8, height=4)


### UMAP
dat <- prot_imputed_nbr %>%
  select(-c(Sample, plate, RNAi)) %>% 
  as.matrix()
rownames(dat) <- prot_imputed_nbr$Sample 

umap_res <- umap(dat, n_neighbors=5, metric="euclidean", random_state=123)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- rownames(dat)
umap_df <- left_join(umap_df, metadata %>% select(Sample, plate, batch, RNAi), 
                     by = "Sample")

# Plot colored by plate
p1 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = plate)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(title = "UMAP - colored by Plate")

# Plot colored by batch
p2 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(title = "UMAP - colored by Batch")

# specific UMAP; RNAi_grouped for plotting
umap_df_temp <- umap_df %>%
  mutate(RNAi_grouped = ifelse(RNAi %in% controls, RNAi, "Other"))
umap_df_temp <- umap_df_temp %>%
  mutate(Plate_label = ifelse(RNAi %in% controls, as.character(plate), NA))
p3 = ggplot(umap_df_temp, aes(x = UMAP1, y = UMAP2, color = RNAi_grouped)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text_repel(
    aes(label = Plate_label), 
    na.rm = TRUE,           # Only labels non-NA, i.e. controls
    size = 3, 
    show.legend = FALSE
  ) +
  theme_classic() +
  labs(title = "UMAP - controls labeled, colored by RNAi", color = "RNAi")

umap_plots <- p1 + p2 + p3
ggsave('results/preprocessing/UMAPplots_prot.png', umap_plots, width=12, height=8)

# Old UMAP
prot_data <- prot_data_raw %>% select(c(1, 20:1229)) %>%
  column_to_rownames('Sample')
prot_to_remove <- c("A53E_DROME",  "CH16_DROME",  "T2FA_DROME",  "MS57A_DROME", "CP304_DROME")
prot_data <- prot_data %>% select(-prot_to_remove)
dat <- prot_data %>% as.matrix()

umap_res <- umap(dat, n_neighbors=5, metric="euclidean", random_state=123)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- rownames(dat)
umap_df <- left_join(umap_df, metadata %>% select(Sample, plate, batch, RNAi), 
                     by = "Sample")

# Plot colored by plate
p1 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = plate)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(title = "UMAP - colored by Plate")

# Plot colored by batch
p2 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(title = "UMAP - colored by Batch")

# specific UMAP; RNAi_grouped for plotting
umap_df_temp <- umap_df %>%
  mutate(RNAi_grouped = ifelse(RNAi %in% controls, RNAi, "Other"))
umap_df_temp <- umap_df_temp %>%
  mutate(Plate_label = ifelse(RNAi %in% controls, as.character(plate), NA))
p3 = ggplot(umap_df_temp, aes(x = UMAP1, y = UMAP2, color = RNAi_grouped)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text_repel(aes(label = Plate_label), 
                  na.rm = TRUE, size = 3, show.legend = FALSE) +
  theme_classic() +
  labs(title = "UMAP - controls labeled, colored by RNAi", color = "RNAi")

umap_plots <- p1 + p2 + p3
ggsave('results/for_figures/UMAPplots_prot_OLD.png', umap_plots, width=13, height=5)

