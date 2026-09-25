### Filtering and preprocessing of metabolomics data

# Remove duplicates/unwanted samples, filter by missingness, impute NAs, 
# normalize, center by plate, remove bad replicates
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggrepel)
library(purrr)

################## READ IN SPECTRAL INTENSITY RAW DATA ##################

metab_metadata <- readxl::read_xlsx("raw_data/Metabolomics_metadata.xlsx", sheet = 1)
metab_metadata <- metab_metadata %>%
  mutate(
    Plate_number = str_extract(Plate, "(?<=Plate)\\d+"),
    Gene = ifelse(is.na(Gene), Genotype, Gene),
    control = ifelse(is.na(Genotype), "RNAi", "Control")
  )
metab_data <- read.csv("raw_data/TF-3-metadata_analysis.csv", na.strings = 0, check.names = F)

# remove outlier control
outlier_control <- c("5-G8","5-G9","5-G10")

################## CLEAN UP METADATA ##################

metab_metadata <- metab_metadata %>%
  select(Sample, Unique_Identifier, Gene, `Control to use2`, Plate_number, `Fly B-day (month-year)`, `Number of flies`, Sex...8) %>%
  rename(Plate = Plate_number, Sample = Sample, Sex = Sex...8, RNAi = Unique_Identifier, Target = Gene, Control = `Control to use2`, Batch = `Fly B-day (month-year)`, Fly_Count = `Number of flies`) %>%
  as.data.frame()
metab_metadata <- metab_metadata %>% mutate(Plate = as.factor(Plate), Batch = as.factor(Batch), Fly_Count = as.factor(Fly_Count))
rownames(metab_metadata) <- metab_metadata$Sample

################## HANDLE DUPLICATE METABOLITES ##################

# Function to resolve duplicate metabolites
resolve_duplicates <- function(data, metadata) {
  # Get column names (excluding Sample and Plate number)
  metabolite_cols <- names(data)[!names(data) %in% c("Sample", "Plate number")]
  
  # Find duplicated names
  duplicated_names <- unique(metabolite_cols[duplicated(metabolite_cols)])
  
  if(length(duplicated_names) > 0) {
    cat("Found duplicated metabolites:", paste(duplicated_names, collapse = ", "), "\n")
    
    # Identify control samples (w1118_F)
    # Adjust this condition based on your actual control identification
    control_samples <- metadata$Sample[metadata$RNAi == "w1118_F" | 
                                         metadata$Target == "w1118" | 
                                         grepl("w1118", metadata$RNAi, ignore.case = TRUE)]
    
    cat("Found", length(control_samples), "control samples\n")
    
    columns_to_remove <- c()
    
    for(dup_name in duplicated_names) {
      # Get all columns with this name
      dup_indices <- which(metabolite_cols == dup_name)
      dup_col_indices <- which(names(data) == dup_name)
      
      cat("\nProcessing duplicates for:", dup_name, "\n")
      cat("Found at positions:", dup_col_indices, "\n")
      
      # Calculate metrics for each duplicate
      metrics <- data.frame(
        col_index = dup_col_indices,
        na_count = numeric(length(dup_col_indices)),
        median_control = numeric(length(dup_col_indices))
      )
      
      for(i in 1:length(dup_col_indices)) {
        col_idx <- dup_col_indices[i]
        col_data <- data[, col_idx]
        
        # Count NAs
        metrics$na_count[i] <- sum(is.na(col_data))
        
        # Calculate median in control samples
        if(length(control_samples) > 0) {
          control_indices <- which(data$Sample %in% control_samples)
          control_values <- col_data[control_indices]
          metrics$median_control[i] <- median(control_values, na.rm = TRUE)
        } else {
          metrics$median_control[i] <- median(col_data, na.rm = TRUE)
        }
      }
      
      # Sort by NA count (ascending), then by median (descending)
      metrics <- metrics[order(metrics$na_count, -metrics$median_control), ]
      
      cat("Metrics for duplicates:\n")
      print(metrics)
      
      # Keep the first one (fewest NAs, highest median in case of tie)
      keep_col <- metrics$col_index[1]
      remove_cols <- metrics$col_index[-1]
      
      cat("Keeping column at position:", keep_col, "\n")
      cat("Removing columns at positions:", remove_cols, "\n")
      
      columns_to_remove <- c(columns_to_remove, remove_cols)
    }
    
    # Remove duplicate columns
    if(length(columns_to_remove) > 0) {
      data <- data[, -columns_to_remove]
      cat("\nRemoved", length(columns_to_remove), "duplicate columns\n")
    }
  } else {
    cat("No duplicate metabolite names found\n")
  }
  
  return(data)
}

# Apply the function to resolve duplicates
metab_data <- resolve_duplicates(metab_data, metab_metadata)

# Continue with the rest of your analysis
# Filter out outlier control
metab_metadata <- metab_metadata %>% filter(! Sample %in% outlier_control)
metab_data <- metab_data %>% filter(! Sample %in% outlier_control)

metab_data <- metab_data %>% select(-`Plate number`)
rownames(metab_data) <- metab_data$Sample

arranging_order <- metab_data$Sample

metab_data <- metab_data[arranging_order,]
metab_metadata <- metab_metadata[arranging_order,]

# Get metabolite columns
metab_columns_unfiltered <- colnames(metab_data[2:ncol(metab_data)])

# Verify no more duplicates
cat("\nFinal check - Duplicate column names:", any(duplicated(names(metab_data))), "\n")
cat("Dataset dimensions:", nrow(metab_data), "samples x", ncol(metab_data)-1, "metabolites\n")

###### FILTERING METABOLITES ######

# Step 1: Filter metabolites by missingness and create summary table
filter_by_missingness <- function(metab_data, metab_columns) {
  
  # Calculate missingness for each metabolite across all samples
  metabolite_missingness <- metab_data %>%
    select(all_of(metab_columns)) %>%
    summarise(across(everything(), ~ sum(is.na(.x)) / length(.x))) %>%
    pivot_longer(everything(), names_to = "Metabolite", values_to = "Missingness_Fraction") %>%
    mutate(
      Missingness_Percent = Missingness_Fraction * 100,
      Missingness_Category = case_when(
        Missingness_Percent < 20 ~ "Low (0-<20%)",
        Missingness_Percent < 50 ~ "Moderate (20-<50%)", 
        Missingness_Percent < 75 ~ "High (50-<75%)",
        Missingness_Percent >= 75 ~ "Filtered (>=75%)"
      )
    ) %>%
    arrange(desc(Missingness_Percent))
  
  # Summary counts
  missingness_summary <- metabolite_missingness %>%
    count(Missingness_Category, name = "N_Metabolites") %>%
    mutate(Percentage = round(N_Metabolites / sum(N_Metabolites) * 100, 1))
  
  # Filter out highly missing metabolites (>=75%)
  kept_metabolites <- metabolite_missingness %>%
    filter(Missingness_Percent < 75) %>%
    pull(Metabolite)
  
  cat("=== MISSINGNESS FILTERING ===\n")
  print(missingness_summary)
  cat("Metabolites kept:", length(kept_metabolites), "out of", nrow(metabolite_missingness), "\n")
  
  return(list(
    metabolite_info = metabolite_missingness,
    summary = missingness_summary,
    kept_metabolites = kept_metabolites
  ))
}

# Step 1: Filter by missingness (same as before)
missingness_results <- filter_by_missingness(metab_data, metab_columns_unfiltered)
filtered_metab <- metab_data %>% select(Sample, missingness_results$kept_metabolites)

# Step 2: Impute with half-minimum per plate
# ---- add Plate information to your intensity data ----
filtered_metab_with_plate <- filtered_metab %>%
  left_join(metab_metadata %>% select(Sample, Plate), by = "Sample")

# Gather to long (Sample, Plate, Metabolite, Value)
metab_long <- filtered_metab_with_plate %>%
  pivot_longer(
    -c(Sample, Plate),
    names_to = "Metabolite",
    values_to = "Value"
  )

# Impute NAs and intensities lower than 1e4 by plate/metabolite with half-minimum
metab_long_imputed <- metab_long %>%
  group_by(Plate, Metabolite) %>%
  mutate(
    # Compute the half-minimum for this metabolite/plate
    HalfMin = 0.5 * min(Value, na.rm = TRUE),
    Value = ifelse(Value < 1e4, HalfMin, Value),
    Value = ifelse(is.na(Value), HalfMin, Value)
  ) %>%
  ungroup() %>%
  select(-HalfMin)

# Check intensities
plate_annot <- metab_long_imputed %>%
  distinct(Sample, Plate) %>%
  mutate(y = -5000000)

ggplot(metab_long_imputed, aes(x = Sample, y = Value, fill = Plate)) +
  geom_boxplot(outlier.size = 0, width = 0.7) +
  geom_tile(data = plate_annot, mapping = aes(x = Sample, y = y, fill = Plate),
            height = 5000000, inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Distribution of Metabolite Intensities by Sample",
       x = "Sample", y = "Abundance", fill = "Plate") +
  ylim(c(-10000000, max(metab_long_imputed$Value))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), # Hide x labels (too many)
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

# ggsave('results/for_figures/preprocessing/intensities_preproc.png', height = 4, width = 12)

# Pivot back to wide format, keeping Plate as needed (or drop Plate column)
filtered_metab_imputed <- metab_long_imputed %>%
  select(-Plate) %>% # or keep Plate if you want
  pivot_wider(
    id_cols = Sample,
    names_from = Metabolite,
    values_from = Value
  ) %>%
  arrange(match(Sample, filtered_metab$Sample)) # keep original order

d4_succinate <- read.csv("raw_data/d4-succinate.csv") %>%
  mutate(Sample = paste0(Plate ,"-", Position)) %>% filter(Type == "Sample")

# normalize by d4-succinate
global_median <- median(d4_succinate$d4.Succinate, na.rm = TRUE)
# Left join d4_succinate to your metabolites by Sample
data_norm <- filtered_metab_imputed %>%
  left_join(select(d4_succinate, Sample, d4.Succinate), by = "Sample")
data_norm <- data_norm %>%
  mutate(norm_factor = global_median / d4.Succinate)
# Identify metabolite columns: all except Sample, d4.Succinate and norm_factor
metab_cols <- setdiff(names(data_norm), c("Sample", "d4.Succinate", "norm_factor"))

# Apply normalization
data_norm <- data_norm %>%
  mutate(across(all_of(metab_cols), ~ .x * norm_factor))
# Optional: drop norm_factor and d4.Succinate columns if no longer needed
final_norm <- data_norm %>%
  select(-d4.Succinate, -norm_factor)

# Step 2.5: Show intensity by plate/batch colored by quartile
log2_filtered_metab_imputed <- log2(final_norm[,colnames(final_norm)[2:ncol(final_norm)]])

# long_df <- log2_filtered_metab_imputed %>%
#   mutate(Sample = metab_metadata$Sample) %>%     # ensure 'Sample' is column
#   left_join(metab_metadata %>% select(Sample, Plate), by = "Sample") %>%
#   pivot_longer(-c(Sample, Plate), names_to = "Metabolite", values_to = "Value")
# 
# ggplot(long_df, aes(x = Sample, y = Value, fill = Plate)) +
#   geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, width = 0.7) +
#   scale_fill_brewer(palette = "Set1") +
#   labs(title = "Distribution of Metabolite Intensities (log2) by Sample",
#        x = "Sample", y = "log2 Abundance", fill = "Plate") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(), # Hide x labels (too many)
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())


# Step 3: Normalize samples by median metabolite intensity
log2_filtered_metab_imputed_norm <- log2_filtered_metab_imputed %>%
  mutate(Sample = metab_metadata$Sample)

metab_cols <- setdiff(colnames(log2_filtered_metab_imputed_norm), "Sample")

log2_filtered_metab_imputed_norm[metab_cols] <-
  t(apply(log2_filtered_metab_imputed_norm[metab_cols], 1,
          function(x) x - median(x, na.rm=TRUE)))
log2_filtered_metab_imputed_norm <- as.data.frame(log2_filtered_metab_imputed_norm)
log2_filtered_metab_imputed_norm$Sample <- metab_metadata$Sample

# Step 3.5: check intensities
# long_df <- log2_filtered_metab_imputed_norm %>%
#   mutate(Sample = metab_metadata$Sample) %>%     # ensure 'Sample' is column
#   left_join(metab_metadata %>% select(Sample, Plate), by = "Sample") %>%
#   pivot_longer(-c(Sample, Plate), names_to = "Metabolite", values_to = "Value")
# 
# ggplot(long_df, aes(x = Sample, y = Value, fill = Plate)) +
#   geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, width = 0.7) +
#   scale_fill_brewer(palette = "Set1") +
#   labs(title = "Distribution of Metabolite Intensities (log2) by Sample",
#        x = "Sample", y = "log2 Abundance", fill = "Plate") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(), # Hide x labels (too many)
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())

# Step 4: Center by plate
log2_centered_by_plate <- log2_filtered_metab_imputed %>%
  mutate(Sample = metab_metadata$Sample) %>%
  left_join(metab_metadata %>% select(Sample, Plate), by = "Sample") %>%
  pivot_longer(-c(Sample, Plate), names_to = "Metabolite", values_to = "Value") %>%
  group_by(Plate, Metabolite) %>%
  mutate(Centered = Value - median(Value, na.rm = TRUE)) %>%
  ungroup() %>%
  select(Sample, Metabolite, Centered) %>%
  pivot_wider(names_from = Metabolite, values_from = Centered)

long_df <- log2_centered_by_plate %>%
  mutate(Sample = metab_metadata$Sample) %>%     # ensure 'Sample' is column
  left_join(metab_metadata %>% select(Sample, Plate), by = "Sample") %>%
  pivot_longer(-c(Sample, Plate), names_to = "Metabolite", values_to = "Value")

ggplot(long_df, aes(x = Sample, y = Value, fill = Plate)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, width = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Distribution of Metabolite Intensities (log2) by Sample",
       x = "Sample", y = "log2 Abundance", fill = "Plate") +
  theme_bw() +
  theme(axis.text.x = element_blank(), # Hide x labels (too many)
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

# Step 4.25: Filter out bad replicates + male RNAi + one replicate + unwanted targets
# Add RNAi data back
norm_data <- log2_centered_by_plate %>%
  left_join(metab_metadata %>% select(Sample, RNAi), by = "Sample")

# Get metabolite columns (all but sample and RNAi) 
metab_cols <- setdiff(colnames(norm_data), c("Sample", "RNAi"))

# Function to calculate correlation within an RNAi group and find bad reps
corr_check <- function(rnai) {
  rnai_name <- unique(rnai$RNAi)
  
  # Skip if less than 3 replicates
  n_reps <- nrow(rnai)
  if (n_reps < 3) {
    return(invisible(NULL))
  }
  
  # Get correlation matrix
  metab_mat <- as.matrix(rnai[, metab_cols])
  rownames(metab_mat) <- rnai$Sample
  corr_mat <- cor(t(metab_mat), use = 'pairwise.complete.obs')
  
  # Calculate each replicate's average correlation to others
  mean_cor <- apply(corr_mat, 1, function(x) mean(x[-which.max(x)], na.rm = TRUE))
  
  # Names of low correlation replicates
  low_cor_reps <- names(mean_cor[mean_cor < 0.5])
  
  # # Print results
  # if (length(low_cor_reps > 0)) {
  #   cat('\nRNAi:', rnai_name, '\n')
  #   cat("Replicates:", paste(rnai$Sample, collapse = ", "), "\n")
  #   cat("Average correlations:\n")
  #   print(round(mean_cor, 3))
  # } else {
  #   cat('\nAll replicates in', rnai_name, 'OK\n')
  # }
  
  # If there are bad replicates, return name of RNAi
  if (length(low_cor_reps > 0)) {
    return(rnai_name)
  }
}

# # Check correlations (to run, uncomment lines in function above)
# norm_data %>%
#   group_split(RNAi) %>%
#   walk(corr_check)

# Get names of RNAIs with bad reps
rnais_bad_reps <- norm_data %>%
  group_by(RNAi) %>%
  group_split() %>%
  map(corr_check) %>%
  compact() %>%
  unlist

# Function to eliminate low correlation reps with backward elimination
remove_bad_reps <- function(norm_data, cutoff = 0.5) {
  reps_to_remove <- list() # List to store replicates to remove
  
  # Separate each RNAi into groups
  for (rnai_name in rnais_bad_reps) {
    rnai <- subset(norm_data, RNAi == rnai_name)
    removed <- character(0)
    
    repeat {
      # Build correlation matrix
      metab_mat <- as.matrix(rnai[, metab_cols])
      rownames(metab_mat) <- rnai$Sample
      corr_mat <- cor(t(metab_mat), use = "pairwise.complete.obs")
      
      # Stop if there are only 3 replicates left
      if (nrow(corr_mat) <= 3) break
      
      # Find the lowest correlation; if it is above the cutoff, stop
      min_corr <- min(corr_mat[lower.tri(corr_mat)], na.rm = TRUE)
      if (min_corr >= cutoff) break
      
      # Find the pair with lowest correlation
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
bad_reps <- remove_bad_reps(norm_data)

# Remove replicates
samples_to_remove <- unlist(bad_reps, use.names = FALSE)
samples_to_remove <- samples_to_remove[!samples_to_remove %in% '3-C2'] # Add back
norm_data_nobadreps <- norm_data %>%
  filter(!Sample %in% samples_to_remove)

# Remove all data with cuff and daw targets and outlier controls + Jra
cuff_samp <- metab_metadata$Sample[metab_metadata$Target == 'cuff']
daw_samp <- metab_metadata$Sample[metab_metadata$Target == 'daw']
jra_samp <- metab_metadata$Sample[metab_metadata$RNAi == 'Jra_7216_F']
samp_rm <- c(cuff_samp, daw_samp, jra_samp)
norm_data_nobadreps <- norm_data_nobadreps %>%
  filter(!Sample %in% samp_rm)
metab_metadata <- metab_metadata %>%
  filter(!Sample %in% samp_rm)

# Filter out RNAi that only have 1 replicate
# names(which(table(metab_data$RNAi) == 1)); abd-A_35644_F, aop_26759_F, ftz-f1_27659_F, Gug_51414_F, Hr39_27086_F, lbl_60001_F, pan_26743_F
norm_data_nobadreps <- norm_data_nobadreps %>%
  add_count(RNAi) %>%
  filter(n > 1) %>%
  select(-n)
metab_metadata <- metab_metadata %>%
  add_count(RNAi) %>%
  filter(n > 1) %>%
  select(-n)

# Filter out male RNAi lines (_M in the name)
male_rnai_samps <- metab_metadata$Sample[grepl('_M', metab_metadata$RNAi)]
norm_data_nobadreps <- norm_data_nobadreps %>%
  filter(!Sample %in% male_rnai_samps)
metab_metadata <- metab_metadata %>%
  filter(!Sample %in% male_rnai_samps)

# write.csv(norm_data_nobadreps, 'processed_data/preprocessed_metab_data.csv', row.names = FALSE)

# Step 4.5: Show intensity by plate/batch colored by quartile
# Ensure Sample column exists and is unique row identifier:
dat <- norm_data_nobadreps %>%
  select(-c(Sample, RNAi)) %>% 
  as.matrix()

rownames(dat) <- norm_data_nobadreps$Sample  # Just in case
meta_nobadreps <- metab_metadata %>%
  filter(!Sample %in% samples_to_remove)
# write.csv(meta_nobadreps, 'processed_data/preprocessed_metab_metadata.csv', row.names = FALSE)
meta_nobadreps$Sample <- as.character(meta_nobadreps$Sample)

# Plot intensities again
long_df <- norm_data_nobadreps %>%
  select(-RNAi) %>%
  left_join(metab_metadata %>% select(Sample, Plate), by = "Sample") %>%
  pivot_longer(-c(Sample, Plate), names_to = "Metabolite", values_to = "Value")

plate_annot <- long_df %>%
  distinct(Sample, Plate) %>%
  mutate(y = -11)
  
ggplot(long_df, aes(x = Sample, y = Value, fill = Plate)) +
  geom_boxplot(outlier.size = 0, width = 0.7) +
  geom_tile(data = plate_annot, mapping = aes(x = Sample, y = y, fill = Plate),
            height = 0.5, inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Distribution of Metabolite Intensities (log2) by Sample",
       x = "Sample", y = "log2(Abundance)", fill = "Plate") +
  ylim(c(-11.5, max(long_df$Value))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), # Hide x labels (too many)
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

# ggsave('results/for_figures/preprocessing/intensities_postproc.png', height = 4, width = 12)


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
  left_join(meta_nobadreps, by = "Sample")

pca_df_temp <- pca_df %>%
  mutate(RNAi_grouped = ifelse(RNAi %in% controls, RNAi, "Other"),
         label = ifelse(RNAi %in% controls, 'Control', NA))
pca_df_temp$RNAi_grouped <- factor(pca_df_temp$RNAi_grouped,
                                   levels = c('w1118_F', 'Other'))

# Plot, colored by Plate
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Plate)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(
    title = "PCA Colored by Plate",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

# Plot, colored by batch
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(
    title = "PCA Colored by Batch",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

p3 <- ggplot(pca_df_temp, aes(x = PC1, y = PC2, color = RNAi_grouped)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text_repel(aes(label = label), na.rm = TRUE, 
                  size = 3, show.legend = FALSE, alpha = 1) +
  theme_classic() +
  labs(
    title = "PCA Colored by RNAi (Controls Labeled)",
    color = "RNAi",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

library(ggrepel)
library(patchwork)

# p <- p1 + p2 + p3
p <- p1 + p3
p

ggsave('results/preprocessing/PCAplots_metab.png', plot=p, height = 8, width = 20)
ggsave('results/for_figures/preprocessing/PCAplots_metab.png', plot=p, height = 4, width = 8)

# Old PCA
metab_data <- read.csv("raw_data/TF-3-metadata_analysis.csv", na.strings = 0, check.names = F)
metab_data <- resolve_duplicates(metab_data, metab_metadata)
metab_data <- metab_data %>%
  select(-`Plate number`) %>%
  tibble::column_to_rownames('Sample')
metab_data[is.na(metab_data)] <- 0
metab_metadata <- readxl::read_xlsx("raw_data/Metabolomics_metadata.xlsx", sheet = 1)
metab_metadata <- metab_metadata %>%
  mutate(Plate_number = str_extract(Plate, "(?<=Plate)\\d+"),
         Gene = ifelse(is.na(Gene), Genotype, Gene),
         control = ifelse(is.na(Genotype), "RNAi", "Control"))
metab_metadata <- metab_metadata %>%
  select(Sample, Unique_Identifier, Gene, `Control to use2`, Plate_number, `Fly B-day (month-year)`, `Number of flies`, Sex...8) %>%
  rename(Plate = Plate_number, Sample = Sample, Sex = Sex...8, RNAi = Unique_Identifier, Target = Gene, Control = `Control to use2`, Batch = `Fly B-day (month-year)`, Fly_Count = `Number of flies`) %>%
  as.data.frame()

pca <- prcomp(metab_data, scale. = TRUE, center = TRUE)  # already median-normalized

# Proportion of variance explained
expl_var <- pca$sdev^2 / sum(pca$sdev^2)
pc1_var <- round(100 * expl_var[1], 1)
pc2_var <- round(100 * expl_var[2], 1)

# Extract PCA scores (first 2 PCs)
pca_df <- as.data.frame(pca$x[,1:2]) %>%
  mutate(Sample = rownames(pca$x)) %>%
  left_join(metab_metadata, by = "Sample")

pca_df_temp <- pca_df %>%
  mutate(RNAi_grouped = ifelse(RNAi %in% controls, RNAi, "Other"),
         label = ifelse(RNAi %in% controls, 'Control', NA))
pca_df_temp$RNAi_grouped <- factor(pca_df_temp$RNAi_grouped,
                                   levels = c('w1118_F', 'Other'))

# Plot, colored by Plate
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Plate)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(
    title = "PCA Colored by Plate",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

# Plot, colored by Plate
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
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

ggsave('results/for_figures/preprocessing/PCAplots_metab_old.png', plot=p, height = 4, width = 8)


# Run UMAP (use default settings or tune as needed)
library(umap)
library(factoextra)  # For PCA visualization if desired

umap_res <- umap(dat, n_neighbors=5, metric="euclidean", random_state=123)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- rownames(dat)

umap_df <- left_join(umap_df, metab_metadata, by = "Sample")

# Plot, colored by Plate
p1 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Plate)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(title = "UMAP - colored by Plate")


# Plot, colored by Plate
p2 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Batch)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(title = "UMAP - colored by Batch")

# specific UMAP:
# 1. Make "RNAi_grouped" for plotting

umap_df_temp <- umap_df %>%
  mutate(RNAi_grouped = ifelse(RNAi %in% controls, RNAi, "Other"))

# 2. Label Plate ONLY for the specified controls
umap_df_temp <- umap_df_temp %>%
  mutate(Plate_label = ifelse(RNAi %in% controls, as.character(Plate), NA))

# 3. Plot -- label only controls
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

p1 + p2 + p3

ggsave('results/preprocessing/UMAPplots_metab.png', plot=p, height = 8, width = 20)

# Old UMAP
umap_res <- umap(metab_data, n_neighbors=5, metric="euclidean", random_state=123)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- rownames(metab_data)

umap_df <- left_join(umap_df, metab_metadata, by = "Sample")

# Plot, colored by Plate
p1 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Plate)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(title = "UMAP - colored by Plate")


# Plot, colored by Plate
p2 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Batch)) +
  geom_point(size = 2, alpha=0.7) + 
  theme_classic() + 
  labs(title = "UMAP - colored by Batch")

# specific UMAP:
# 1. Make "RNAi_grouped" for plotting

umap_df_temp <- umap_df %>%
  mutate(RNAi_grouped = ifelse(RNAi %in% controls, RNAi, "Other"))

# 2. Label Plate ONLY for the specified controls
umap_df_temp <- umap_df_temp %>%
  mutate(Plate_label = ifelse(RNAi %in% controls, as.character(Plate), NA))

# 3. Plot -- label only controls
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

p1 + p2 + p3

ggsave('results/for_figures/UMAPplots_metab_OLD.png', plot=p, height = 5, width = 12)


# Step 7: Heatmap with correlation of replicates (RNAi)

# Average metabolite data by RNAi
averaged_data <- norm_data_nobadreps %>%
  group_by(RNAi) %>%
  summarise(across(-Sample, mean, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("RNAi") %>%
  as.matrix()

# Calculate RNAi-to-RNAi correlation
rnai_cor <- cor(t(averaged_data), use = "complete.obs")

# Create RNAi-level metadata (take first occurrence of each RNAi)
rnai_metadata <- metab_metadata %>%
  group_by(RNAi) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(RNAi, Target, Plate, Batch) %>%
  column_to_rownames("RNAi")

# Ensure row order matches correlation matrix
rnai_annotation <- rnai_metadata[rownames(rnai_cor), ]

library(pheatmap)
library(RColorBrewer)

# Define annotation colors for ALL RNAi (keep full colors)
annotation_colors <- list(
  Target = rainbow(length(unique(rnai_annotation$Target))),
  Plate = RColorBrewer::brewer.pal(n = length(unique(rnai_annotation$Plate)), name = "Set2"),
  Batch = RColorBrewer::brewer.pal(n = length(unique(rnai_annotation$Batch)), name = "Set3")
)

# Name the colors
names(annotation_colors$Target) <- unique(rnai_annotation$Target)
names(annotation_colors$Plate) <- unique(rnai_annotation$Plate)
names(annotation_colors$Batch) <- unique(rnai_annotation$Batch)

# CREATE CUSTOM LABELS: Only show labels for controls
controls <- c("w1118_F", "attP2_F", "attP40_F")

# Row labels - only controls get their names, others get empty strings
custom_row_labels <- ifelse(rownames(rnai_cor) %in% controls, 
                            rownames(rnai_cor), 
                            "")

# Column labels - only controls get their names, others get empty strings  
custom_col_labels <- ifelse(colnames(rnai_cor) %in% controls, 
                            colnames(rnai_cor), 
                            "")

# Create heatmap with custom labels but full annotations
pheatmap::pheatmap(
  rnai_cor,
  annotation_row = rnai_annotation,        # Full annotation (colors for all)
  annotation_col = rnai_annotation,        # Full annotation (colors for all)
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  labels_row = custom_row_labels,          # Custom row labels (only controls)
  labels_col = custom_col_labels,          # Custom column labels (only controls)
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "RNAi-to-RNAi Correlation Heatmap (Controls Labeled)"
)
