### Linear regression to determine which proteins affect metabolites
library(dplyr)
library(tibble)
library(tidyverse)
library(glmnet)
library(forcats)
library(ggpubr)
library(ggridges)

# Read in data
metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE)
metab_metadata <- read.csv('processed_data/filtered_metab_metadata.csv')
prot_data <- read.csv('processed_data/filtered_prot_data.csv', check.names = FALSE)
prot_metadata <- read.csv('processed_data/filtered_prot_metadata.csv')
metab_paths <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')

# Uncomment to filter down to metabolic proteins
# prot_paths <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv') 
# prots_keep <- prot_paths %>% filter(grepl('Metabolism', classes)) %>% pull(protein)

# Combine dfs
intersect <- intersect(prot_metadata$Sample.info1, metab_metadata$Sample)

protein <- prot_data %>% 
  left_join(prot_metadata %>% select(Sample, Sample.info1), by = 'Sample') %>%
  filter(Sample.info1 %in% intersect) %>%
  arrange(Sample.info1) %>%
  column_to_rownames('Sample.info1') %>%
  #select(any_of(prots_keep)) # Uncomment for metabolic prots
  select(-c(Sample, RNAi)) 
metab <- metab_data %>%
  filter(Sample %in% intersect) %>%
  arrange(Sample) %>%
  column_to_rownames('Sample') %>%
  select(-RNAi) 

# Uncomment to change prots to FBgns
# mapf <- read.csv('raw_data/uniprot_to_fbgn.csv')
# protein <- protein %>% select(-c(LAM0_DROME, LAMC_DROME)) # Only if using all prots
# prots <- as.data.frame(colnames(protein)) %>% 
#   rename(protein = `colnames(protein)`) %>%
#   left_join(mapf, by = 'protein')
# colnames(protein) <- prots$fbgn

# Make sure tables are sorted
all(rownames(metab) == rownames(protein))
X <- protein
Y <- metab

# Test with subsets of metabs
# top_ten_met <- c('AC(10:0)', 'AC(12:1)', 'AC(14:1-OH)', 'AC(10:1)', 'AC(8:0)', 
#                  'FA(14:0)', '5-Oxoproline', 'dTMP', 'AC(12-OH)', 'N6-N6-N6-Trimethyl-L-lysine')
# test3 <- c('AC(10:0)', '5-Oxoproline', 'dTMP', 'N6-N6-N6-Trimethyl-L-lysine',
#            'Malate', 'D-Glucose', 'NAD+')
# cysmet_metabs <- metab_paths %>%
#   filter(grepl('Cysteine and methionine metabolism', pathways)) %>%
#   pull(metabolite)
# argpro_metabs <- metab_paths %>%
#   filter(grepl('Arginine and proline metabolism', pathways)) %>%
#   pull(metabolite)
# Y <- Y %>% select(all_of(argpro_metabs))

hits_plots <- list()
prediction_plots <- list()
coef_dfs <- list()

# Loop to run lasso and plot
for (i in colnames(Y)) {
  
  y <- as.matrix(Y[[i]])
  
  if(any(is.na(y))) {
    next
  }
  
  set.seed(12345)
  
  lasso_model <- cv.glmnet(
    as.matrix(X),
    y,
    nfolds = 10
  )
  
  coef_lasso <- as.matrix(coef(lasso_model, s = "lambda.min"))[-1, , drop = FALSE]
  coef_matrix <- as.data.frame(coef_lasso)
  colnames(coef_matrix) <- "coef"
  coef_matrix$protein <- rownames(coef_matrix)
  
  coef_selected <- coef_matrix[coef_matrix$coef != 0, ]
  
  if(nrow(coef_selected) == 0) {
    next
  }
  
  coef_selected <- coef_selected %>%
    mutate(protein = fct_reorder(protein, coef))
  
  if (nrow(coef_selected) > 40) {
    coef_selected <- bind_rows(
      slice_min(coef_selected, coef, n = 20),
      slice_max(coef_selected, coef, n = 20)
    ) %>%
      distinct(protein, .keep_all = TRUE)
  }
  
  hits_plot <- ggplot(coef_selected, aes(x = protein, y = coef)) +
    geom_bar(stat = "identity", fill = "#f68060", alpha = 0.6, width = 0.4) +
    coord_flip() +
    theme_bw() +
    labs(
      title = paste("LASSO-selected proteins for", i),
      x = "",
      y = "Coefficient"
    )
  
  hits_plots[[i]] <- hits_plot
  
  predictions <- predict(lasso_model, as.matrix(X), s = "lambda.min")
  
  plot_df <- data.frame(
    predicted = as.numeric(predictions),
    actual = as.numeric(y)
  )
  
  prediction_plot <- ggscatter(
    plot_df,
    x = "predicted",
    y = "actual",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    xlab = "LASSO prediction",
    ylab = paste(i, "abundance")
  )
  
  prediction_plots[[i]] <- prediction_plot
  
  coef_selected <- coef_selected %>%
    mutate(metabolite = i)
  
  coef_dfs[[i]] <- coef_selected

}

pdf("results/linear_regression/v1_filtered/all_hitplots.pdf", width = 8, height = 8)
for (p in hits_plots) {
  print(p)}
dev.off()

pdf("results/linear_regression/v1_filtered/all_predplots.pdf", width = 10, height = 8)
for (p in prediction_plots) {
  print(p)}
dev.off()

coef_df <- bind_rows(coef_dfs)
coef_df <- coef_df %>% relocate(metabolite)
write.csv(coef_df, 'results/linear_regression/v1_filtered/coefficient_data.csv', row.names = FALSE)



## Get lists of predictors
predictor_list <- list()

for (i in colnames(Y)) {
  y <- as.matrix(Y[[i]])
  
  if(any(is.na(y))) {next}
  
  set.seed(12345)
  
  # Build model
  lasso_model <- cv.glmnet(as.matrix(X),
                           y,
                           nfolds = 10)
  
  # Extract coefficients
  coef_lasso <- as.matrix(coef(lasso_model, s = "lambda.min"))[-1, , drop = FALSE]
  coef_matrix <- data.frame(protein = rownames(coef_lasso),
                            coef = coef_lasso[,1])
  coef_selected <- coef_matrix %>% filter(coef != 0)
  
  # Skip if no predictors
  if(nrow(coef_selected) == 0) {next}
  
  # Collapse to string and write out
  pred_string <- paste0(coef_selected$protein, collapse = '; ')
  predictor_list[[i]] <- data.frame(metabolite = i,
                                    predictors = pred_string)
}  

predictor_df <- bind_rows(predictor_list)  
write.csv(predictor_df, 'results/linear_regression/predictor_lists.csv', row.names = FALSE)  
  
  

