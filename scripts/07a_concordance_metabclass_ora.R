### Concordance analysis downstream; filter down to metabolism class, 
### check associations of each top metab/prot
### ORA on top pairs to validate novel pairs
library(tidyverse)
library(magrittr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

# Read in data 
metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE)
metab_metadata <- read.csv('processed_data/filtered_metab_metadata.csv')
prot_data <- read.csv('processed_data/filtered_prot_data.csv', check.names = FALSE)
prot_metadata <- read.csv('processed_data/filtered_prot_metadata.csv')
prot_path <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
proteins_metab_class <- prot_path %>% filter(grepl('Metabolism', classes))

# Get list of samples in each
intersect <- intersect(prot_metadata$Sample.info1, metab_metadata$Sample)
prot_not_metab <- setdiff(prot_metadata$Sample.info1, metab_metadata$Sample)
metab_not_prot <- setdiff(metab_metadata$Sample, prot_metadata$Sample.info1)
all_samps <- c(intersect, prot_not_metab, metab_not_prot)
sample_locs <- data.frame(Sample = all_samps)
sample_locs$Location <- ifelse(sample_locs$Sample %in% intersect, 'Both',
                               ifelse(sample_locs$Sample %in% prot_not_metab, 'Protein', 'Metab'))
# Get sample to RNAi information
rnai_to_samp_prot <- prot_metadata %>% select(Sample.info1, RNAi) %>% rename(Sample = Sample.info1)
rnai_to_samp_metab <- metab_metadata %>% filter(Sample %in% metab_not_prot) %>% select(Sample, RNAi)
rnai_to_samp <- rbind(rnai_to_samp_prot, rnai_to_samp_metab)
sample_locs <- sample_locs %>%
  left_join(rnai_to_samp, by = 'Sample')

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
met_ctrl_mean <- colMeans(metab_data[controls, ], na.rm = TRUE)

# center all samples to control mean
prot_centered <- sweep(prot_data, 2, prot_ctrl_mean, "-")
met_centered <- sweep(metab_data, 2, met_ctrl_mean, "-")

# Scale data
prot_scaled <- scale(prot_centered)
met_scaled <- scale(met_centered)

# Read in condordance results, filter down to proteins in metabolism class
conc_results <- read.csv('results/concordance/conc_results.csv')
conc_results <- conc_results %>% filter(protein %in% proteins_metab_class$protein)

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
  geom_text_repel(data = label_df, aes(label = label), 
                  size = 3.5, box.padding = 0.5, show.legend = FALSE) +
  theme_minimal(base_size = 14)
ggsave('results/concordance/met_class_prots/volcano.png', height=8, width=10)
ggsave('results/for_figures/concordance/volcano_metclass.png', height=8, width=10)

# Of significant pairs, what are the top 20 and bottom 20 pairs
conc_results_sig <- conc_results %>% filter(padj < 0.05)
top_bottom_20_pmi <- conc_results %>% arrange(desc(concordance))
top_bottom_20_pmi <- top_bottom_20_pmi %>%  slice(1:20, (n() - 19):n())
write.csv(top_bottom_20_pmi, 'results/concordance/met_class_prots/top_bottom_20_pmi.csv', row.names = FALSE)

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
      x = paste("Metabolite Log2FC:", metabolite),
      y = paste("Protein Log2FC:", protein),
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

pdf("results/concordance/met_class_prots/top10_ci_plots.pdf", width = 10, height = 8)
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
    size = 4.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black") +
  theme_bw() +
  labs(title = 'Significant Correlations per Metabolite (Top 10 Labeled)',
       x = "Sorted metabolites",
       y = "Number of significant correlations") +
  theme(panel.grid = element_blank())
ggsave('results/concordance/met_class_prots/cor_per_metab.png', height=6, width=10)

ggplot(prot_counts, aes(x = rank, y = n_sig)) +
  geom_point(color = "grey40", size = 2) +
  geom_point(data = top10_prot, color = "red", size = 2.5) +
  # Labels only for top 10
  geom_text_repel(
    data = top10_prot,
    aes(label = paste0(protein, " (", n_sig, ")")),
    size = 4.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black") +
  theme_bw() +
  labs(title = 'Significant Correlations per Protein (Top 10 Labeled)',
       x = "Sorted proteins",
       y = "Number of significant correlations") +
  theme(panel.grid = element_blank())
ggsave('results/concordance/met_class_prots/cor_per_prot.png', dpi=800, height=6, width=10)

# Get all prots per top metab and vice versa
topmets <- conc_results %>%
  group_by(metabolite) %>%
  summarise(n_sig = sum(padj < 0.05),
            n_total = n()) %>%
  arrange(desc(n_sig)) %>%
  mutate(rank = row_number()) %>%
  slice_head(n = 15) %>%
  pull(metabolite)
topprots <- conc_results %>%
  group_by(protein) %>%
  summarise(n_sig = sum(padj < 0.05),
            n_total = n()) %>%
  arrange(desc(n_sig)) %>%
  mutate(rank = row_number()) %>%
  slice_head(n = 15) %>%
  pull(protein)

topmets_assoc <- conc_results %>%
  filter(metabolite %in% topmets) %>%
  group_by(metabolite) %>%
  summarise(associated = paste0(protein, collapse = '; '),
            .groups = 'drop') %>%
  mutate(type = 'metabolite') %>%
  rename(analyte = metabolite) %>%
  select(c(analyte, type, associated))
topprots_assoc <- conc_results %>%
  filter(protein %in% topprots) %>%
  group_by(protein) %>%
  summarise(associated = paste0(metabolite, collapse = '; '),
            .groups = 'drop') %>%
  mutate(type = 'protein') %>%
  rename(analyte = protein) %>%
  select(c(analyte, type, associated))
associations <- rbind(topmets_assoc, topprots_assoc)
write.csv(associations, 'results/concordance/met_class_prots/top15met_prot_associations.csv', row.names = FALSE)

### For top 10 proteins with most significant pairs, get list of metabs associated with it, do ORA on them
### And vice versa for top 10 metabs

## Top 10 proteins - ORA on their metabolites
top10prot <- top10_prot$protein
conc_top10prot <- conc_results %>% filter(protein %in% top10prot) %>% filter(padj < 0.05)

# Get list of sig metabs per protein
sig_metabs_per_prot <- conc_top10prot %>%
  dplyr::select(protein, metabolite) %>%
  arrange(protein) 

# Get pathway data
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')
path_classes <- read.csv("raw_data/kegg_mapping/kegg_pathways_manual.csv")

# Get valid path and subclass names for only metabolism class 
valid_subcl <- unique(path_classes$subclass[path_classes$class == 'Metabolism'])
valid_subcl <- valid_subcl[!valid_subcl %in% 'Global and overview maps']

# Filter metab path to valid subclasses
metab_path_all <- metab_path_all %>% 
  separate_rows(subclasses, sep = ";\\s*") %>%
  filter(subclasses %in% valid_subcl) %>%
  select(metabolite, subclasses) %>%
  distinct()

# get all metabolites tested as background
all_metabs <- unique(conc_results$metabolite)

# Function to perform ORA analysis per protein
ora_per_prot <- function(prot) {
  # Get list of metabs per protein
  metabs_sig <- sig_metabs_per_prot %>%
    filter(protein == prot) %>%
    pull(metabolite) %>%
    unique()
  
  # Get subclasses to test
  subclasses <- unique(metab_path_all$subclasses)
  
  map_df(subclasses, function(subcl) {
    metabs_in_path <- metab_path_all %>%
      filter(subclasses == subcl) %>%
      pull(metabolite) %>%
      unique()
    
    ### Fisher's test setup
    # metabs in path in protein     | metabs not in path in protein 
    # metabs in path not in protein | metabs not in path not in protein
    a <- length(intersect(metabs_sig, metabs_in_path))
    b <- length(setdiff(metabs_sig, metabs_in_path))
    c <- length(setdiff(metabs_in_path, metabs_sig))
    d <- length(setdiff(all_metabs, union(metabs_sig, metabs_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow=2)
    
    pval <- fisher.test(mat)$p.value
    
    tibble(protein = prot,
           subclass = subcl,
           overlap = a,
           p.value = pval)
    
  }) %>%
    mutate(padj = p.adjust(p.value, method = 'BH'),
           enrichment = -log10(padj))
}

# Run for top 10 prots
ora_results <- map_df(top10prot, ora_per_prot)

# Plot pathway enrichment per protein
plot_cluster <- function(prot){
  ora_results %>%
    filter(protein == prot,
           enrichment > 0,
           subclass != "Global and overview maps") %>%
    arrange(enrichment) %>%
    tail(20) %>%
    ggplot(aes(enrichment, reorder(subclass, enrichment))) +
    geom_point(aes(size = overlap, fill = enrichment), shape=21) +
    scale_fill_continuous() +
    labs(
      title = paste(prot, "Pathway Enrichment"),
      x = "-log10(adjusted p-value)",
      y = "Pathway"
    ) +
    theme_bw()
}

plots <- map(top10prot, plot_cluster) 
# Some proteins have no significant subclasses
p <- (plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] +
        plots[[5]] + plots[[6]] + plots[[9]] + plots[[10]] + plot_spacer()) +
  plot_layout(ncol = 3)
ggsave('results/concordance/met_class_prots/subclass_per_top10prot.png', p, width=20, height=12)


## Protein - pathway level
valid_paths <- unique(path_classes$pathway_name[path_classes$class == 'Metabolism']) 

# Filter metab path to valid subclasses
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv') %>% 
  separate_rows(pathways, sep = ";\\s*") %>% 
  filter(pathways %in% valid_paths) %>% 
  select(metabolite, pathways) %>% 
  distinct()

# get all metabolites tested as background
all_metabs <- unique(conc_results$metabolite)

# Function to perform ORA analysis per protein
ora_per_prot <- function(prot) {
  # Get list of metabs per protein
  metabs_sig <- sig_metabs_per_prot %>%
    filter(protein == prot) %>%
    pull(metabolite) %>%
    unique()
  
  # Get subclasses to test
  pathways <- unique(metab_path_all$pathways) 
  
  map_df(pathways, function(path) {
    metabs_in_path <- metab_path_all %>% 
      filter(pathways == path) %>% 
      pull(metabolite) %>%
      unique()
    
    ### Fisher's test setup
    # metabs in path in protein     | metabs not in path in protein 
    # metabs in path not in protein | metabs not in path not in protein
    a <- length(intersect(metabs_sig, metabs_in_path))
    b <- length(setdiff(metabs_sig, metabs_in_path))
    c <- length(setdiff(metabs_in_path, metabs_sig))
    d <- length(setdiff(all_metabs, union(metabs_sig, metabs_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow=2)
    
    pval <- fisher.test(mat)$p.value
    
    tibble(protein = prot,
           pathway = path, 
           overlap = a,
           p.value = pval)
    
  }) %>%
    mutate(padj = p.adjust(p.value, method='BH'),
           enrichment = -log10(padj))
}

# Run for top 10 prots
ora_results <- map_df(top10prot, ora_per_prot)

# Plot pathway enrichment per protein
plot_cluster <- function(prot){
  ora_results %>%
    filter(protein == prot,
           enrichment > 0.25) %>%
    arrange(enrichment) %>%
    tail(20) %>%
    ggplot(aes(enrichment, reorder(pathway, enrichment))) + 
    geom_point(aes(size = overlap, fill = enrichment), shape=21) +
    scale_fill_continuous() +
    labs(
      title = paste(prot, "Pathway Enrichment"),
      x = "-log10(adjusted p-value)",
      y = "Pathway"
    ) +
    theme_bw()
}

plots <- map(top10prot, plot_cluster) 
# Save plots
p <- (plots[[1]] + plots[[2]] + plots[[3]] + plots[[6]] + plots[[9]] + plots[[10]]) +
  plot_layout(ncol = 3)
ggsave('results/concordance/met_class_prots/pathway_per_top10prot.png', p, width=20, height=8)


######### Same but for metabolites top 10
## Top 10 metabolites - ORA on their proteins
#conc_results <- read.csv('processed_data/concordance/conc_results.csv') # Uncomment to run this on ALL proteins not just metabolism-associated
top10metab <- top10_met$metabolite
conc_top10metab <- conc_results %>% filter(metabolite %in% top10metab) %>% filter(padj < 0.05)

# Get list of sig proteins per metab
sig_prots_per_metab <- conc_top10metab %>%
  dplyr::select(metabolite, protein) %>%
  arrange(metabolite) 

# Get pathway data
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')

# Filter prot path to valid subclasses
prot_path_all <- prot_path_all %>% 
  separate_rows(subclasses, sep = ";\\s*") %>%
  filter(subclasses %in% valid_subcl) %>%
  select(protein, subclasses) %>%
  distinct()

# get all proteins as background
all_prots <- unique(conc_results$protein)

# Function to perform ORA analysis per metabolite
ora_per_metab <- function(metab) {
  # Get list of prots per metab
  prots_sig <- sig_prots_per_metab %>%
    filter(metabolite == metab) %>%
    pull(protein) %>%
    unique()
  
  # Get subclasses to test
  subclasses <- unique(prot_path_all$subclasses)
  
  map_df(subclasses, function(subcl) {
    prots_in_path <- prot_path_all %>%
      filter(subclasses == subcl) %>%
      pull(protein) %>%
      unique()
    
    ### Fisher's test setup
    # prots in path in metabolite     | prots not in path in metabolite 
    # prots in path not in metabolite | prots not in path not in metabolite
    a <- length(intersect(prots_sig, prots_in_path))
    b <- length(setdiff(prots_sig, prots_in_path))
    c <- length(setdiff(prots_in_path, prots_sig))
    d <- length(setdiff(all_prots, union(prots_sig, prots_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow=2)
    
    pval <- fisher.test(mat)$p.value
    
    tibble(metabolite = metab,
           subclass = subcl,
           overlap = a,
           p.value = pval)
    
  }) %>%
    mutate(padj = p.adjust(p.value, method='BH'),
           enrichment = -log10(padj))
}

# Run for top 10 metabs
ora_results_metab <- map_df(top10metab, ora_per_metab)

# Plot pathway enrichment per protein
plot_cluster <- function(metab){
  ora_results_metab %>%
    filter(metabolite == metab,
           enrichment > 0,
           subclass != "Global and overview maps") %>%
    arrange(enrichment) %>%
    tail(20) %>%
    ggplot(aes(enrichment, reorder(subclass, enrichment))) +
    geom_point(aes(size = overlap, fill = enrichment), shape=21) +
    scale_fill_continuous() +
    labs(
      title = paste(metab, "Pathway Enrichment"),
      x = "-log10(adjusted p-value)",
      y = "Pathway"
    ) +
    theme_bw()
}

plots <- map(top10metab, plot_cluster) 
p <- (plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] +
        plots[[8]] + plots[[9]] + plot_spacer() + plot_spacer()) +
  plot_layout(ncol = 3)
ggsave('results/concordance/met_class_prots/subclass_per_top10met.png', p, width=20, height=12)


# Pathway level
# Filter prot path to valid subclasses
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv') %>% 
  separate_rows(pathways, sep = ";\\s*") %>%
  filter(pathways %in% valid_paths) %>%
  select(protein, pathways) %>%
  distinct()

# get all proteins as background
all_prots <- unique(conc_results$protein)

# Function to perform ORA analysis per metabolite
ora_per_metab <- function(metab) {
  # Get list of prots per metab
  prots_sig <- sig_prots_per_metab %>%
    filter(metabolite == metab) %>%
    pull(protein) %>%
    unique()
  
  # Get pathways to test
  pathways <- unique(prot_path_all$pathways)
  
  map_df(pathways, function(path) {
    prots_in_path <- prot_path_all %>%
      filter(pathways == path) %>%
      pull(protein) %>%
      unique()
    
    ### Fisher's test setup
    # prots in path in metabolite     | prots not in path in metabolite 
    # prots in path not in metabolite | prots not in path not in metabolite
    a <- length(intersect(prots_sig, prots_in_path))
    b <- length(setdiff(prots_sig, prots_in_path))
    c <- length(setdiff(prots_in_path, prots_sig))
    d <- length(setdiff(all_prots, union(prots_sig, prots_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow=2)
    
    pval <- fisher.test(mat)$p.value
    
    tibble(metabolite = metab,
           pathway = path,
           overlap = a,
           p.value = pval)
    
  }) %>%
    mutate(padj = p.adjust(p.value, method='BH'),
           enrichment = -log10(padj))
}

# Run for top 10 metabs
ora_results_metab <- map_df(top10metab, ora_per_metab)

# Plot pathway enrichment per protein
plot_cluster <- function(metab){
  ora_results_metab %>%
    filter(metabolite == metab,
           enrichment > 0) %>%
    arrange(enrichment) %>%
    tail(20) %>%
    ggplot(aes(enrichment, reorder(pathway, enrichment))) +
    geom_point(aes(size = overlap, fill = enrichment), shape=21) +
    scale_fill_continuous() +
    labs(
      title = paste(metab, "Pathway Enrichment"),
      x = "-log10(adjusted p-value)",
      y = "Pathway"
    ) +
    theme_bw()
}

plots <- map(top10metab, plot_cluster) 
p <- (plots[[3]] + plots[[5]] + plots[[6]] + plots[[10]]) +
  plot_layout(ncol = 2)
ggsave('results/concordance/met_class_prots/pathway_per_top10met.png', p, width=12, height=8)


### Enhanced visualization just for TML-KAD1
# Read in condordance results, filter down to proteins in metabolism class
conc_results <- read.csv('results/concordance/conc_results.csv')
conc_results <- conc_results %>% filter(protein %in% proteins_metab_class$protein)

### Plot concordance for TML-KAD1
conc_res_filt <- conc_results %>% filter(metabolite == 'N6-N6-N6-Trimethyl-L-lysine', protein == 'KAD1_DROME')

df_plot <- data.frame(metabolite = met_scaled[, 'N6-N6-N6-Trimethyl-L-lysine'], 
                      protein = prot_scaled[, 'KAD1_DROME']) %>%
  mutate(quadrant = case_when(
    metabolite >= 0 & protein >= 0 ~ "Concordant (Up)",
    metabolite <= 0 & protein <= 0 ~ "Concordant (Down)",
    metabolite >= 0 & protein <= 0 ~ "Discordant",
    metabolite <= 0 & protein >= 0 ~ "Discordant")) %>%
  tibble::rownames_to_column('Sample') %>%
  left_join(metab_metadata %>% select(Sample, Target), by = 'Sample') %>% 
  mutate(Outline = ifelse(grepl('Attp', Target), 'Control Line', 'Perturbed'),
         Outline = ifelse(Target == 'Mondo', 'Mondo', Outline),
         Outline = ifelse(Target == 'bigmax', 'Bigmax', Outline))

ggplot(df_plot, aes(x = metabolite, y = protein, fill = quadrant, color = Outline)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, alpha = 0.8, shape = 21) +
  scale_color_manual(values = c('Control Line' = 'brown', 'Mondo' = 'black', 'Bigmax' = 'navyblue')) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, 
           label = paste0("Concordance = ", round(conc_res_filt$concordance, 3), 
                          "\nAdjusted p-value = ", signif(conc_res_filt$padj, 3) ), size = 4 ) +
  geom_text_repel(aes(label = Target)) +
  labs(
    title = 'Concordance Quadrant Plot: KAD1 vs N6-N6-N6-Trimethyl-L-lysine',
    x = "Metabolite Log2FC: TML",
    y = "Protein Log2FC: KAD1",
    fill = "Quadrant",
    color = 'Targets of Interest') +
  theme_minimal(base_size = 14)

ggsave('results/concordance/tml_kad1/conc_plot.png', height = 10, width = 12)


# TF per quadrant
quads <- data.frame(metabolite = met_scaled[, 'N6-N6-N6-Trimethyl-L-lysine'], 
                    protein = prot_scaled[, 'KAD1_DROME']) %>%
  mutate(quadrant = case_when(
    metabolite >= 0 & protein >= 0 ~ "Concordant (Up)",
    metabolite <= 0 & protein <= 0 ~ "Concordant (Down)",
    metabolite >= 0 & protein <= 0 ~ "Discordant (Met Up Prot Down)",
    metabolite <= 0 & protein >= 0 ~ "Discordant (Prot Up Met Down)")) %>%
  tibble::rownames_to_column('Sample') %>%
  left_join(metab_metadata %>% select(Sample, Target), by = 'Sample') %>%
  select(Target, quadrant) %>%
  distinct()
write.csv(quads, 'results/concordance/tml_kad1/quadrants.csv', row.names = FALSE)



### Making spoke graphs (top metabs -> associated proteins and vice versa) for top concordance results
rm(list=ls())
library(tidygraph)
library(ggraph)
library(tibble)
library(dplyr)
library(patchwork)

# Read in data
conc_results <- read.csv('results/concordance/conc_results.csv')
prot_path <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
proteins_metab_class <- prot_path %>% filter(grepl('Metabolism', classes))

# Filter down to metab-associated proteins
conc_results <- conc_results %>% filter(protein %in% proteins_metab_class$protein)

metabs_keep <- c('Citrate', 'Coumarate', 'N6-N6-N6-Trimethyl-L-lysine', 'AC(14:1-OH)', 
                 '6-Phospho-D-gluconate', 'D-Ribose 5-diphosphate', '2-Hydroxyglutarate')
top_mets <- conc_results %>%
  filter(padj < 0.05) %>% 
  distinct(metabolite, protein) %>%
  count(metabolite, name = "n_prot") %>% 
  arrange(desc(n_prot)) %>%
  #slice_head(n = 3)
  filter(metabolite %in% metabs_keep)

prot_lists <- conc_results %>%
  filter(padj < 0.05) %>%
  semi_join(top_mets, by = 'metabolite') %>%
  distinct(metabolite, protein)

# Reg summaries
regsum_prot <- read.csv('results/diff_analysis/reg_summary_prot.csv') %>% select(protein, sig)
regsum_met <- read.csv('results/diff_analysis/reg_summary_metab.csv') %>% select(metabolite, sig)

# Loop through each metab to get graph
metabs_top3 <- top_mets$metabolite
plots <- list()

for (i in 1:length(metabs_top3)) {
  metab <- metabs_top3[[i]]
  proteins <- unique(prot_lists$protein[prot_lists$metabolite == metab])
  
  # Set up data for graph
  nodes <- tibble(name = c(metab, proteins),
                  x = c(0, cos(seq(0, 2*pi, length.out = length(proteins)+1)[-1])),
                  y = c(0, sin(seq(0, 2*pi, length.out = length(proteins)+1)[-1])))
  
  # Add regulation direction and number of TF data
  conc_res_filt <- conc_results %>% 
    filter(metabolite == metab, padj < 0.05) %>%
    select(protein, concordance) %>%
    mutate(direction = ifelse(concordance > 0.5, 'Concordant', 'Discordant'),
           n_tfs = ifelse(direction == 'Concordant', concordance * 159, (1 - concordance) * 159)) %>%
    rename(name = protein) %>%
    select(-concordance) %>%
    left_join(regsum_prot %>% rename(name = protein), by = 'name')
  
  edges <- tibble(from = 1, to = 2:(length(proteins)+1)) %>%
    left_join(conc_res_filt %>%
                select(name, n_tfs) %>%
                mutate(to = match(name, nodes$name)), by = "to")
  
  nodes <- nodes %>% 
    left_join(conc_res_filt, by = 'name') %>%
    mutate(direction = ifelse(name == metab, 'NA', direction),
           n_tfs = ifelse(name == metab, 100, n_tfs)) %>%
    left_join(regsum_met %>% rename(name = metabolite), by = 'name') %>%
    mutate(sig.x = ifelse(is.na(sig.x), sig.y, sig.x)) %>%
    rename(sig = sig.x) %>% select(-sig.y)
  
  graph <- tbl_graph(nodes = nodes, edges = edges)
  
  # Make graph
  g <- ggraph(graph, layout = "manual", x = x, y = y) +
    geom_edge_link(aes(width = n_tfs)) +
    scale_edge_width(range = c(0.5, 1.75)) +
    geom_node_point(aes(color = direction, size = sig)) +
    geom_node_label(aes(label = name), repel = TRUE, force = 10,
                    segment.color = NA, box.padding = 1, point.padding = 0.5) +
    scale_color_manual(values = c("Concordant" = "tomato", "Discordant" = "steelblue")) +
    scale_size_continuous(range = c(5,10)) +
    theme_void() +
    labs(color = 'Direction', size = 'Number of TFs Regulating', edge_width = 'Number of Concordant Pairs')
  plots[[i]] <- g
}


p <- wrap_plots(plots, ncol = 2) & theme(plot.margin = margin(25, 25, 25, 25))
ggsave('results/concordance/met_class_prots/prots_per_metab_spokes_volcano_top.png', p, width = 16, height = 16)

### Proteins
prots_keep <- c('VATE_DROME', 'VATD1_DROME', 'IMDH_DROME', 'DCUP_DROME', 
                'PROD_DROME', 'KAD1_DROME', 'GSTD3_DROME', 'SDHA_DROME')
top_prots <- conc_results %>%
  filter(padj < 0.05) %>% 
  distinct(metabolite, protein) %>%
  count(protein, name = "n_met") %>% 
  arrange(desc(n_met)) %>%
  #slice_head(n = 3)
  filter(protein %in% prots_keep)

met_lists <- conc_results %>%
  filter(padj < 0.05) %>%
  semi_join(top_prots, by = 'protein') %>%
  distinct(metabolite, protein)


# Loop through each prot to get graph
prots_top3 <- top_prots$protein
plots <- list()

for (i in 1:length(prots_top3)) {
  prot <- prots_top3[[i]]
  metabs <- unique(met_lists$metabolite[met_lists$protein == prot])
  
  # Set up data for graph
  nodes <- tibble(name = c(prot, metabs),
                  x = c(0, cos(seq(0, 2*pi, length.out = length(metabs)+1)[-1])),
                  y = c(0, sin(seq(0, 2*pi, length.out = length(metabs)+1)[-1]))) 
  
  # Add regulation direction and number of TF data
  conc_res_filt <- conc_results %>% 
    filter(protein == prot, padj < 0.05) %>%
    select(metabolite, concordance) %>%
    mutate(direction = ifelse(concordance > 0.5, 'Concordant', 'Discordant'),
           n_tfs = ifelse(direction == 'Concordant', concordance * 159, (1 - concordance) * 159)) %>%
    rename(name = metabolite) %>%
    select(-concordance) %>%
    left_join(regsum_met %>% rename(name = metabolite), by = 'name')
  
  edges <- tibble(from = 1, to = 2:(length(metabs)+1)) %>%
    left_join(conc_res_filt %>%
                select(name, n_tfs) %>%
                mutate(to = match(name, nodes$name)), by = "to")
  
  nodes <- nodes %>% 
    left_join(conc_res_filt, by = 'name') %>%
    mutate(direction = ifelse(name == prot, 'NA', direction),
           n_tfs = ifelse(name == prot, 100, n_tfs)) %>%
    left_join(regsum_prot %>% rename(name = protein), by = 'name') %>%
    mutate(sig.x = ifelse(is.na(sig.x), sig.y, sig.x)) %>%
    rename(sig = sig.x) %>% select(-sig.y)
  
  graph <- tbl_graph(nodes = nodes, edges = edges)
  
  # Make graph
  g <- ggraph(graph, layout = "manual", x = x, y = y) +
    geom_edge_link(aes(width = n_tfs)) +
    scale_edge_width(range = c(0.5, 1.75)) +
    geom_node_point(aes(color = direction, size = sig)) +
    geom_node_label(aes(label = name), repel = TRUE, force = 10,
                    segment.color = NA, box.padding = 1, point.padding = 0.25) +
    scale_color_manual(values = c("Concordant" = "tomato", "Discordant" = "steelblue")) +
    scale_size_continuous(range = c(5,10)) +
    theme_void() +
    labs(color = 'Direction', size = 'Number of TFs Regulating', edge_width = 'Number of Concordant Pairs')
  plots[[i]] <- g
}


p <- wrap_plots(plots, ncol = 2) & theme(plot.margin = margin(25, 25, 25, 25))
ggsave('results/concordance/met_class_prots/network_graphs/metabs_per_prot_spokes_volcanotop.png', p, width = 16, height = 20)


