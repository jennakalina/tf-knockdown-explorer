### Concordance analysis on filtered proteomics and metabolomics data 
library(readxl)
library(tidyverse)
library(magrittr)
library(purrr)
library(parallel)
library(survival)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Read in data (change file paths to run at different thresholds)
metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE)
metab_metadata <- read.csv('processed_data/filtered_metab_metadata.csv')
prot_data <- read.csv('processed_data/filtered_prot_data.csv', check.names = FALSE)
prot_metadata <- read.csv('processed_data/filtered_prot_metadata.csv')

# Get list of samples in each
intersect <- intersect(prot_metadata$Sample.info1, metab_metadata$Sample)
prot_not_metab <- setdiff(prot_metadata$Sample.info1, metab_metadata$Sample)
metab_not_prot <- setdiff(metab_metadata$Sample, prot_metadata$Sample.info1)
all_samps <- c(intersect, prot_not_metab, metab_not_prot)

# Filter down data to only samples that are in the intersection
metab_data <- metab_data %>% filter(Sample %in% intersect)
metab_metadata <- metab_metadata %>% filter(Sample %in% intersect)

prot_data <- prot_data %>% 
  left_join(prot_metadata %>% select(Sample, Sample.info1), by = 'Sample') %>%
  rename(Sample1 = Sample, Sample = Sample.info1) %>%
  relocate(Sample) %>%
  filter(Sample %in% intersect) %>%
  select(-Sample1)
prot_metadata <- prot_metadata %>%
  rename(Sample1 = Sample, Sample = Sample.info1) %>%
  relocate(Sample) %>%
  filter(Sample %in% intersect)

prot_data <- prot_data %>% tibble::column_to_rownames('Sample') %>% select (-RNAi)
metab_data <- metab_data %>% tibble::column_to_rownames('Sample') %>% select (-RNAi)

# Center data to control
controls <- c("2-E2", "2-E3", "2-E4")

# compute control means
prot_ctrl_mean <- colMeans(prot_data[controls, ], na.rm = TRUE)
met_ctrl_mean  <- colMeans(metab_data[controls, ], na.rm = TRUE)

# center all samples to control mean
prot_centered <- sweep(prot_data, 2, prot_ctrl_mean, "-")
met_centered  <- sweep(metab_data, 2, met_ctrl_mean, "-")

# Scale data
prot_scaled <- scale(prot_centered)
met_scaled  <- scale(met_centered)

# Calculate concordance
proteins <- colnames(prot_scaled)
metabolites <- colnames(met_scaled)

# Make sure the same samples are present
intersect_samples <- intersect(rownames(prot_scaled), rownames(met_scaled))

prot_scaled <- prot_scaled[intersect_samples, , drop = FALSE]
met_scaled <- met_scaled[intersect_samples, , drop = FALSE]

# Ensure identical ordering
met_scaled <- met_scaled[rownames(prot_scaled), , drop = FALSE]

# Check
stopifnot(identical(rownames(prot_scaled), rownames(met_scaled)))

conc_results <- mclapply(proteins, function(g){
  x <- prot_scaled[, g]  # protein vector
  res <- lapply(metabolites, function(m){
    y <- met_scaled[, m]  # metabolite vector
    # remove NA
    keep <- !is.na(x) & !is.na(y)
    x_full <- x[keep]
    y_full <- y[keep]
    
    # concordance
    conc <- survival::concordance(y_full ~ x_full)
    data.frame(
      protein = g,
      metabolite = m,
      n = length(x_full),
      concordance = conc$concordance,
      variance = conc$var,
      zscore = (conc$concordance - 0.5) / sqrt(conc$var)
    )
  })
  
  do.call(rbind, res)
  
}, mc.cores = detectCores() - 1) %>% do.call(rbind, .)

conc_results <- conc_results %>%
  mutate(p.value = 2 * pnorm(-abs(zscore)),
         padj = p.adjust(p.value, method = "BH")) %>%
  arrange(padj, -concordance)
write.csv(conc_results, 'results/concordance/conc_results.csv', row.names = FALSE)

### Show which prot-metab interactions show strong interactions
# How many are significant
sum(conc_results$padj < 0.05) 
volcano_df <- conc_results %>% 
  mutate(neglog10p = -log10(padj), 
         sig = case_when(padj < 0.05 & concordance > 0.5 ~ "Concordant (significant)", 
                         padj < 0.05 & concordance < 0.5 ~ "Discordant (significant)", 
                         padj > 0.05 ~ 'Not significant',
                         TRUE ~ "Not significant"),
         label = paste(protein, metabolite, sep = "-"))
# Set labels for only top 5 concordant/discordant
label_df <- volcano_df %>%
  filter(sig != "Not significant") %>%
  group_by(sig) %>%
  slice_max(order_by = neglog10p, n = 5, with_ties = FALSE) %>%
  ungroup()

ggplot(volcano_df, aes(x = concordance, y = neglog10p, color = sig)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Concordant (significant)" = "firebrick3", 
                                "Discordant (significant)" = "royalblue", 
                                "Not significant" = "gray50" )) +
  labs(title = "Volcano Plot of Protein–Metabolite Concordance",
       x = "Concordance Index",
       y = "-log10(p-value)",
       color = "Category") +
  geom_text_repel(data = label_df, aes(label = label), size = 3.5) +
  theme_minimal(base_size = 14)
#ggsave('results/concordance/all_prots/volcano.png', height=8, width=10)

# Of significant pairs, what are the top 20 and bottom 20 pairs
conc_results_sig <- conc_results %>% filter(padj < 0.05)
top_bottom_20_pmi <- conc_results_sig %>% arrange(desc(concordance))
top_bottom_20_pmi <- top_bottom_20_pmi %>% slice(1:20, (n() - 19):n())
write.csv(top_bottom_20_pmi, 'results/concordance/all_prots/top_bottom_20_pmi.csv', row.names = FALSE)

# For the top 10 PMI pairs, make the concordance plot
# Function to plot
plot_quadrant <- function(protein, metabolite, conc_value, padj_value, prot_scaled, met_scaled) {
  x <- met_scaled[, metabolite]
  y <- prot_scaled[, protein]
  
  df_plot <- data.frame(metabolite = x, protein = y) %>%
    mutate(quadrant = case_when(
        metabolite >= 0 & protein >= 0 ~ "Concordant (Up)",
        metabolite <= 0 & protein <= 0 ~ "Concordant (Down)",
        metabolite >= 0 & protein <= 0 ~ "Discordant",
        metabolite <= 0 & protein >= 0 ~ "Discordant"))
  
  ggplot(df_plot, aes(x = metabolite, y = protein, color = quadrant)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(size = 3, alpha = 0.8) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, 
             label = paste0("Concordance = ", round(conc_value, 3), 
                            "\nAdjusted p-value = ", signif(padj_value, 3) ), size = 4 ) +
    labs(
      title = paste("Concordance Quadrant Plot:", protein, "vs", metabolite),
      x = paste("Metabolite abundance:", metabolite),
      y = paste("Protein abundance:", protein),
      color = "Quadrant") +
    theme_minimal(base_size = 14)
}

# Get only top 10
top10 <- top_bottom_20_pmi[1:10, ]
plots <- lapply(1:10, function(i) {
  plot_quadrant(
    protein = top10$protein[i],
    metabolite = top10$metabolite[i],
    conc_value = top10$concordance[i],
    padj_value = top10$padj[i],
    prot_scaled = prot_scaled,
    met_scaled = met_scaled
  )
})

pdf("results/concordance/all_prots/top10_ci_plots.pdf", width = 10, height = 8)
for (p in plots) {
  print(p)}
dev.off()


## How many sig proteins per metab and vice versa
# Count significant correlations per metabolite
met_counts <- conc_results %>%
  group_by(metabolite) %>%
  summarise(n_sig = sum(padj < 0.05),
            n_total = n()) %>%
  arrange(desc(n_sig)) %>%
  mutate(rank = row_number())
top10_met <- met_counts %>%
  slice_head(n = 10)
# Count significant correlations per protein
prot_counts <- conc_results %>%
  group_by(protein) %>%
  summarise(n_sig = sum(padj < 0.05),
            n_total = n()) %>%
  arrange(desc(n_sig)) %>%
  mutate(rank = row_number())
top10_prot <- prot_counts %>%
  slice_head(n = 10)

# Plot
ggplot(met_counts, aes(x = rank, y = n_sig)) +
  geom_point(color = "grey40", size = 2) +
  geom_point(data = top10_met, color = "red", size = 2.5) +
  # Labels only for top 10
  geom_text_repel(
    data = top10_met,
    aes(label = paste0(metabolite, " (", n_sig, ")")),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    max.overlaps = 20) +
  theme_bw() +
  labs(x = "Sorted metabolites",
       y = "Number of significant correlations") +
  theme(panel.grid = element_blank())
ggsave('results/concordance/all_prots/cor_per_metab.png', height=6, width=10)

ggplot(prot_counts, aes(x = rank, y = n_sig)) +
  geom_point(color = "grey40", size = 2) +
  geom_point(data = top10_prot, color = "red", size = 2.5) +
  # Labels only for top 10
  geom_text_repel(
    data = top10_prot,
    aes(label = paste0(protein, " (", n_sig, ")")),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    max.overlaps = 30) +
  theme_bw() +
  labs(x = "Sorted proteins",
       y = "Number of significant correlations") +
  theme(panel.grid = element_blank())
ggsave('results/concordance/all_prots/cor_per_prot.png', height=6, width=10)



