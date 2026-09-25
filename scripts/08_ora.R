### Pathway over-representation analysis on the TF clusters
# Author: Jenna Kalina
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

# Load data (TF FOR PROTEOMICS FIRST)
# clust_mem_prot <- read.csv('results/corr_analysis/clust_mem/cluster_mem_tf_prot_metclassonly_k9.csv') %>% rename(tf = rnai)
clust_mem_prot <- read.csv('results/rCCA/clustering/cluster_membership.csv') %>% rename(tf = TF)
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv("raw_data/kegg_mapping/kegg_pathways_manual.csv")
#daps <- read.csv('results/diff_analysis/DAPs_sig.csv')
daps <- read.csv('results/diff_analysis/DAPs_all_discovery.csv')
daps <- daps %>% filter(p.adj < 0.1)
  
# Get valid path and subclass names for only metabolism class 
valid_subcl <- unique(path_classes$subclass[path_classes$class == 'Metabolism'])

# Get list of sig proteins per TF
prots_per_tf <- daps %>%
  filter(target %in% clust_mem_prot$tf) %>%
  group_by(target) %>%
  summarise(prots = paste0(protein, collapse = '; '),
            .groups = 'drop') %>%
  rename(tf = target)

# Add pathway data
prot_path_all <- prot_path_all %>% 
  select(protein, subclasses)

# Expand lists of analytes per tf
tf_prot_long <- clust_mem_prot %>%
  left_join(prots_per_tf, by = 'tf') %>%
  separate_rows(prots, sep = ";\\s*") %>%
  rename(protein = prots) %>%
  distinct()
tf_prot_paths <- tf_prot_long %>%
  left_join(prot_path_all, by = 'protein') %>%
  separate_rows(subclasses, sep = ";\\s*") %>%
  filter(subclasses %in% valid_subcl)
  #mutate(subclasses = ifelse(grepl('metab', subclasses, ignore.case = TRUE), subclasses, 'Unannotated'))


# Set background
all_prot_data <- read.csv('processed_data/filtered_prot_data.csv') %>% select(-c(Sample, RNAi))
all_prots <- colnames(all_prot_data)
#all_prots <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv') %>% 
#  filter(grepl('Metabolism', classes)) %>% pull(protein) %>% unique()

prot_path_exp <- prot_path_all %>%
  separate_rows(subclasses, sep = ";\\s*") %>%
  filter(subclasses %in% valid_subcl) %>%
  #mutate(subclasses = ifelse(grepl('metab', subclasses, ignore.case = TRUE), subclasses, 'Unannotated')) %>%
  distinct()

# Function to run cluster enrichment
clust_enrichment <- function(clust) {
  tfs_in_clust <- tf_prot_paths %>%
    filter(cluster == clust) %>%
    pull(tf) %>%
    unique()
  
  prots_in_clust <- tf_prot_paths %>%
    filter(tf %in% tfs_in_clust) %>%
    pull(protein) %>%
    unique()
  
  pathways <- unique(tf_prot_paths$subclasses) 
  
  res <- map_df(pathways, function(path) {
    prots_in_path <- prot_path_exp %>%
      filter(subclasses == path) %>% 
      pull(protein) %>%
      unique()
    
    ### Fisher's test setup
    # prots in path in cluster     | prots not in path in cluster 
    # prots in path not in cluster | prots not in path not in cluster
    
    a <- length(intersect(prots_in_clust, prots_in_path))
    b <- length(setdiff(prots_in_path, prots_in_clust))
    c <- length(setdiff(prots_in_clust, prots_in_path))
    d <- length(setdiff(all_prots, union(prots_in_clust, prots_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    
    pval <- fisher.test(mat, alternative = 'greater')$p.value
    
    tibble(cluster = clust,
           pathway = path,
           overlap = a,
           p.value = pval)
  })
  res <- res %>% mutate(p.adj = p.adjust(p.value, method = 'BH'),
                        enrichment = -log10(p.adj))
  return(res)
}

clusters <- sort(unique(clust_mem_prot$cluster))
enrich_results_prot <- map_df(clusters, clust_enrichment)

# Plot all in one
all_plot <- enrich_results_prot %>%
  filter(enrichment > 1.3,
         pathway != "Unannotated") %>%
  group_by(cluster) %>%
  slice_max(enrichment, n = 10) %>%
  ungroup() %>%
  ggplot(aes(x = factor(cluster), y = pathway)) +
  geom_point(aes(fill = enrichment, size = overlap), shape = 21) +
  scale_fill_viridis_c(name = '-log10(adj p-value)') +
  labs(x = 'Cluster', y = 'Pathway', title = 'Pathway Enrichment per TF Cluster (Protein Data)') +
  theme_bw()
all_plot
ggsave('results/ora/tf_tf_annotation/subclass_all_tf_prot_k9_metclassonly.png', all_plot, height=5, width=7)
ggsave('results/rCCA/clustering/ora/subclass_all_tf_prot_k6.png', all_plot, height=5, width=8)

# DAPs discovery, filter down to metab only
enrich_clust1 <- enrich_results_prot %>% 
  filter(cluster == 1) %>%
  mutate(sig = ifelse(enrichment > 5, 'Adjusted p-value < 1e-5', 'Not significant'))

ggplot(enrich_clust1, aes(x = reorder(pathway, p.adj), y = enrichment)) + 
  geom_point(aes(size = overlap, fill = sig), shape = 21) +
  scale_fill_manual(values = c('salmon', 'grey')) +
  geom_hline(yintercept = 5, color = 'gray50') +
  theme_bw() +
  labs(title = 'Over-represented Pathways in Cluster 1 Protein Data',
       x = 'Pathway Subclass', y = 'Enrichment',
       size = '# of Significant Proteins in Pathway',
       fill = 'Significance') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('results/rCCA/clustering/ora/prot_clust1_subclasses.png', height=7, width=9)

# DAPs normal significance, don't filter down to metab only
all_plot <- enrich_results_prot %>%
  filter(p.value < 0.1,
         pathway != "Unannotated") %>%
  group_by(cluster) %>%
  slice_max(enrichment, n = 10) %>%
  ungroup() %>%
  ggplot(aes(x = factor(cluster), y = pathway)) +
  geom_point(aes(fill = enrichment, size = overlap), shape = 21) +
  scale_fill_viridis_c(name = '-log10(adj p-value)') +
  labs(x = 'Cluster', y = 'Pathway', title = 'Pathway Enrichment per TF Cluster (Protein Data)') +
  theme_bw()
all_plot
ggsave('results/rCCA/clustering/ora/subclass_all_tf_prot_k6_ALLPROT.png', height=5, width=8)


## PATHWAY LEVEL
rm(list=ls())
# clust_mem_prot <- read.csv('results/corr_analysis/clust_mem/cluster_mem_tf_prot_metclassonly_k9.csv') %>% rename(tf = rnai)
clust_mem_prot <- read.csv('results/rCCA/clustering/cluster_membership.csv') %>% rename(tf = TF)
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv("raw_data/kegg_mapping/kegg_pathways_manual.csv")
daps <- read.csv('results/diff_analysis/DAPs_sig.csv')

# Get valid path and subclass names for only metabolism class 
valid_path <- unique(path_classes$pathway_name[path_classes$class == 'Metabolism'])

# Get list of sig proteins per TF
prots_per_tf <- daps %>%
  filter(target %in% clust_mem_prot$tf) %>%
  group_by(target) %>%
  summarise(prots = paste0(protein, collapse = '; '),
            .groups = 'drop') %>%
  rename(tf = target)

# Add pathway data
prot_path_all <- prot_path_all %>% 
  select(protein, pathways)

# Expand lists of analytes per tf
tf_prot_long <- clust_mem_prot %>%
  left_join(prots_per_tf, by = 'tf') %>%
  separate_rows(prots, sep = ";\\s*") %>%
  rename(protein = prots) %>%
  distinct()
tf_prot_paths <- tf_prot_long %>%
  left_join(prot_path_all, by = 'protein') %>%
  separate_rows(pathways, sep = ";\\s*") %>%
  filter(pathways %in% valid_path)


# Set background
all_prot_data <- read.csv('processed_data/filtered_prot_data.csv') %>% select(-c(Sample, RNAi))
all_prots <- colnames(all_prot_data)

prot_path_exp <- prot_path_all %>%
  separate_rows(pathways, sep = ";\\s*") %>%
  filter(pathways %in% valid_path) %>%
  distinct()

# Function to run cluster enrichment
clust_enrichment <- function(clust) {
  tfs_in_clust <- tf_prot_paths %>%
    filter(cluster == clust) %>%
    pull(tf) %>%
    unique()
  
  prots_in_clust <- tf_prot_paths %>%
    filter(tf %in% tfs_in_clust) %>%
    pull(protein) %>%
    unique()
  
  pathways <- unique(tf_prot_paths$pathways) 
  
  res <- map_df(pathways, function(path) {
    prots_in_path <- prot_path_exp %>%
      filter(pathways == path) %>% 
      pull(protein) %>%
      unique()
    
    ### Fisher's test setup
    # prots in path in cluster     | prots not in path in cluster 
    # prots in path not in cluster | prots not in path not in cluster
    
    a <- length(intersect(prots_in_clust, prots_in_path))
    b <- length(setdiff(prots_in_path, prots_in_clust))
    c <- length(setdiff(prots_in_clust, prots_in_path))
    d <- length(setdiff(all_prots, union(prots_in_clust, prots_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    
    pval <- fisher.test(mat, alternative = 'greater')$p.value
    
    tibble(cluster = clust,
           pathway = path,
           overlap = a,
           p.value = pval)
  })
  res <- res %>% mutate(p.adj = p.adjust(p.value, method = 'BH'),
                        enrichment = -log10(p.adj))
  return(res)
}

clusters <- sort(unique(clust_mem_prot$cluster))
enrich_results_prot <- map_df(clusters, clust_enrichment)

# Plot all in one
variance_pathways <- enrich_results_prot %>%
  filter(enrichment > 1.3,
         pathway != "Unannotated") %>%
  complete(pathway, cluster, fill = list(overlap = 0)) %>%
  group_by(pathway) %>%
  summarise(variance = var(overlap),
            .groups = 'drop') %>%
  filter(variance > 2) %>%
  pull(pathway)
all_plot_df <- enrich_results_prot %>%
  filter(enrichment > 1.3,
         pathway %in% variance_pathways) %>%
  group_by(cluster) %>%
  ungroup() 
all_plot <- ggplot(all_plot_df, aes(x = factor(cluster), y = pathway)) +
  geom_point(aes(fill = enrichment, size = overlap), shape = 21) +
  scale_fill_viridis_c(name = '-log10(adj p-value)') +
  labs(x = 'Cluster', y = 'Pathway', title = 'Pathway Enrichment per TF Cluster (Protein Data)') +
  theme_bw()
all_plot
ggsave('results/ora/tf_tf_annotation/pathway_all_prot_k9_metclassonly.png', all_plot, height=5, width=9)
ggsave('results/rCCA/clustering/ora/path_all_tf_prot_k6.png', all_plot, height=5, width=8)


### Metabolite now
# Load data
rm(list=ls())
# clust_mem_metab <- read.csv('results/corr_analysis/clust_mem/cluster_mem_tf_metab_k7.csv') %>% rename(tf = rnai)
clust_mem_metab <- read.csv('results/rCCA/clustering/cluster_membership.csv') %>% rename(tf = TF)
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')
path_classes <- read.csv("raw_data/kegg_mapping/kegg_pathways_manual.csv")
#dams <- read.csv('results/diff_analysis/DAMs_sig.csv')
dams <- read.csv('results/diff_analysis/DAMs_all_discovery.csv')
dams <- dams %>% filter(p.adj < 0.1)

# Get valid path and subclass names for only metabolism class 
valid_subcl <- unique(path_classes$subclass[path_classes$class == 'Metabolism'])

# Get list of sig proteins per TF
metabs_per_tf <- dams %>%
  filter(Target %in% clust_mem_metab$tf) %>%
  group_by(Target) %>%
  summarise(metabs = paste0(metabolite, collapse = '; '),
            .groups = 'drop') %>%
  rename(tf = Target)

# Add pathway data
metab_path_all <- metab_path_all %>% 
  select(metabolite, subclasses)

# Expand lists of analytes per tf
tf_metab_long <- clust_mem_metab %>%
  left_join(metabs_per_tf, by = 'tf') %>%
  separate_rows(metabs, sep = ";\\s*") %>%
  rename(metabolite = metabs) %>%
  distinct()
tf_metab_paths <- tf_metab_long %>%
  left_join(metab_path_all, by = 'metabolite') %>%
  separate_rows(subclasses, sep = ";\\s*") %>%
  filter(subclasses %in% valid_subcl)


# Set background
all_metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE) %>% select(-c(Sample, RNAi))
all_metabs <- colnames(all_metab_data)

metab_path_exp <- metab_path_all %>%
  separate_rows(subclasses, sep = ";\\s*") %>%
  filter(subclasses %in% valid_subcl) %>%
  distinct()

# Function to run cluster enrichment
clust_enrichment <- function(clust) {
  tfs_in_clust <- tf_metab_paths %>%
    filter(cluster == clust) %>%
    pull(tf) %>%
    unique()
  
  metabs_in_clust <- tf_metab_paths %>%
    filter(tf %in% tfs_in_clust) %>%
    pull(metabolite) %>%
    unique()
  
  pathways <- unique(tf_metab_paths$subclasses) 
  
  res <- map_df(pathways, function(path) {
    metabs_in_path <- metab_path_exp %>%
      filter(subclasses == path) %>% 
      pull(metabolite) %>%
      unique()
    
    ### Fisher's test setup
    # prots in path in cluster     | prots not in path in cluster 
    # prots in path not in cluster | prots not in path not in cluster
    
    a <- length(intersect(metabs_in_clust, metabs_in_path))
    b <- length(setdiff(metabs_in_clust, metabs_in_path))
    c <- length(setdiff(metabs_in_path, metabs_in_clust))
    d <- length(setdiff(all_metabs, union(metabs_in_clust, metabs_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    
    pval <- fisher.test(mat, alternative = 'greater')$p.value
    
    tibble(cluster = clust,
           pathway = path,
           overlap = a,
           p.value = pval)
  })
  #res <- res %>% mutate(p.adj = p.adjust(p.value, method = 'BH'),
  #                      enrichment = -log10(p.adj))
  res <- res %>% mutate(enrichment = -log10(p.value))
  return(res)
}

clusters <- sort(unique(clust_mem_metab$cluster))
enrich_results_metab <- map_df(clusters, clust_enrichment)

# Plot all in one
all_plot <- enrich_results_metab %>%
  filter(enrichment > 0.5,
         !pathway %in% c("Unannotated", 'Global and overview maps')) %>%
  add_row(cluster = 2, pathway = 'Energy metabolism', overlap = NA, p.value = NA, enrichment = NA) %>% ### ONLY TO INCLUDE CLUSTER 2
  group_by(cluster) %>%
  slice_max(enrichment, n = 10) %>%
  ungroup() %>%
  ggplot(aes(x = factor(cluster), y = pathway)) +
  geom_point(aes(fill = enrichment, size = overlap), shape = 21) +
  scale_fill_viridis_c(name = '-log10(adj p-value)') +
  labs(x = 'Cluster', y = 'Pathway', title = 'Pathway Enrichment per TF Cluster (Metabolite Data)') +
  theme_bw()
all_plot
ggsave('results/ora/tf_tf_annotation/subclass_all_tf_metab_k7.png', all_plot, height=5, width=8)
ggsave('results/rCCA/clustering/ora/subclass_all_tf_metab_k6.png', all_plot, height=5, width=8)


## PATHWAY LEVEL
rm(list=ls())
#clust_mem_metab <- read.csv('results/corr_analysis/clust_mem/cluster_mem_tf_metab_k7.csv') %>% rename(tf = rnai)
clust_mem_metab <- read.csv('results/rCCA/clustering/cluster_membership.csv') %>% rename(tf = TF)
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')
path_classes <- read.csv("raw_data/kegg_mapping/kegg_pathways_manual.csv")
dams <- read.csv('results/diff_analysis/DAMs_sig.csv')

# Get valid path and subclass names for only metabolism class 
valid_path <- unique(path_classes$pathway_name[path_classes$class == 'Metabolism'])

# Get list of sig proteins per TF
metabs_per_tf <- dams %>%
  filter(Target %in% clust_mem_metab$tf) %>%
  group_by(Target) %>%
  summarise(metabs = paste0(metabolite, collapse = '; '),
            .groups = 'drop') %>%
  rename(tf = Target)

# Add pathway data
metab_path_all <- metab_path_all %>% 
  select(metabolite, pathways)

# Expand lists of analytes per tf
tf_met_long <- clust_mem_metab %>%
  left_join(metabs_per_tf, by = 'tf') %>%
  separate_rows(metabs, sep = ";\\s*") %>%
  rename(metabolite = metabs) %>%
  distinct()
tf_met_paths <- tf_met_long %>%
  left_join(metab_path_all, by = 'metabolite') %>%
  separate_rows(pathways, sep = ";\\s*") %>%
  filter(pathways %in% valid_path)


# Set background
all_metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE) %>% select(-c(Sample, RNAi))
all_mets <- colnames(all_metab_data)

met_path_exp <- metab_path_all %>%
  separate_rows(pathways, sep = ";\\s*") %>%
  filter(pathways %in% valid_path) %>%
  distinct()

# Function to run cluster enrichment
clust_enrichment <- function(clust) {
  tfs_in_clust <- tf_met_paths %>%
    filter(cluster == clust) %>%
    pull(tf) %>%
    unique()
  
  mets_in_clust <- tf_met_paths %>%
    filter(tf %in% tfs_in_clust) %>%
    pull(metabolite) %>%
    unique()
  
  pathways <- unique(tf_met_paths$pathways) 
  
  res <- map_df(pathways, function(path) {
    mets_in_path <- met_path_exp %>%
      filter(pathways == path) %>% 
      pull(metabolite) %>%
      unique()
    
    ### Fisher's test setup
    # prots in path in cluster     | prots not in path in cluster 
    # prots in path not in cluster | prots not in path not in cluster
    
    a <- length(intersect(mets_in_clust, mets_in_path))
    b <- length(setdiff(mets_in_path, mets_in_clust))
    c <- length(setdiff(mets_in_clust, mets_in_path))
    d <- length(setdiff(all_mets, union(mets_in_clust, mets_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    
    pval <- fisher.test(mat, alternative = 'greater')$p.value
    
    tibble(cluster = clust,
           pathway = path,
           overlap = a,
           p.value = pval)
  })
  res <- res %>% mutate(p.adj = p.adjust(p.value, method = 'BH'),
                        enrichment = -log10(p.value))
  return(res)
}

clusters <- sort(unique(clust_mem_metab$cluster))
enrich_results_metab <- map_df(clusters, clust_enrichment)

# Plot all in one
all_plot_df <- enrich_results_metab %>%
  filter(enrichment > 1.2)
all_plot <- ggplot(all_plot_df, aes(x = factor(cluster), y = pathway)) +
  geom_point(aes(fill = enrichment, size = overlap), shape = 21) +
  scale_fill_viridis_c(name = '-log10(p-value)') +
  labs(x = 'Cluster', y = 'Pathway', title = 'Pathway Enrichment per TF Cluster (Metabolite Data)') +
  theme_bw()
all_plot
ggsave('results/ora/tf_tf_annotation/pathway_all_met_k7.png', all_plot, height=5, width=9)
ggsave('results/rCCA/clustering/ora/path_all_tf_metab_k6.png', all_plot, height=5, width=8)


### Now on met-met and prot-prot clusters for identification
rm(list=ls())
# Load data
clust_mem_metab <- read.csv('results/corr_analysis/clust_mem/cluster_mem_metab_k4.csv') 
clust_mem_prot <- read.csv('results/corr_analysis/clust_mem/cluster_mem_prot_k4.csv') 
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv("raw_data/kegg_mapping/kegg_pathways_manual.csv")

# Get valid path and subclass names for only metabolism class 
valid_paths <- unique(path_classes$pathway_name[path_classes$class == 'Metabolism'])
valid_subcl <- unique(path_classes$subclass[path_classes$class == 'Metabolism'])

# Filter back down to just metab/prot and cluster
clust_metab <- clust_mem_metab %>% select(metabolite, cluster) %>% distinct() %>% left_join(metab_path_all %>% select(-kegg_id), by = 'metabolite')
clust_prot <- clust_mem_prot %>% select(protein, cluster) %>% distinct() %>% left_join(prot_path_all %>% select(-kegg_id), by = 'protein')

# Expand pathway
clust_metab_exp <- clust_metab %>% 
  separate_rows(subclasses, sep = ";\\s*") %>% #### CHANGE for pathway/subclass
  filter(subclasses %in% valid_subcl) #### CHANGE
clust_prot_exp <- clust_prot %>% 
  separate_rows(subclasses, sep = ";\\s*") %>% #### CHANGE
  filter(subclasses %in% valid_subcl) #### CHANGE

# For each cluster, are the metabolites within it correlated with any certain pathway
distinct_metab <- clust_metab_exp %>%
  select(metabolite, cluster, subclasses) %>% #### CHANGE
  distinct()
all_metabs <- unique(distinct_metab$metabolite)

# Function to run cluster enrichment
clust_enrichment <- function(clust) {
  metabs_in_clust <- distinct_metab %>%
    filter(cluster == clust) %>%
    pull(metabolite) %>%
    unique()
  
  pathways <- unique(distinct_metab$subclasses) ### CHANGE
  
  res <- map_df(pathways, function(path) {
    metabs_in_path <- distinct_metab %>%
      filter(subclasses == path) %>% ### CHANGE
      pull(metabolite) %>%
      unique()
    
    ### Fisher's test setup
    # metabs in path in cluster     | metabs not in path in cluster 
    # metabs in path not in cluster | metabs not in path not in cluster
    
    a <- length(intersect(metabs_in_clust, metabs_in_path))
    b <- length(setdiff(metabs_in_clust, metabs_in_path))
    c <- length(setdiff(metabs_in_path, metabs_in_clust))
    d <- length(setdiff(all_metabs, union(metabs_in_clust, metabs_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    
    pval <- fisher.test(mat)$p.value
    
    tibble(cluster = clust,
           pathway = path,
           overlap = a,
           p.value = pval)
  })
  #res %>% mutate(p.adj = p.adjust(p.value, method = 'BH'),
  #               enrichment = -log10(p.adj))
  res %>% mutate(enrichment = -log10(p.value))
}

clusters <- sort(unique(distinct_metab$cluster))
enrich_results_metab <- map_df(clusters, clust_enrichment)

# Plot pathway enrichment per cluster
plot_cluster <- function(clust){
  enrich_results_metab %>%
    filter(cluster == clust,
           enrichment > 0.25,
           pathway != "Unannotated") %>%
    arrange(enrichment) %>%
    tail(10) %>%
    ggplot(aes(enrichment, reorder(pathway, enrichment))) +
    geom_point(aes(size = overlap, fill = enrichment), shape=21) +
    scale_fill_continuous() +
    labs(
      title = paste("Cluster", clust, "Pathway Enrichment"),
      x = "-log10(p-value)",
      y = "Pathway"
    ) +
    theme_bw()
}

plots <- map(clusters, plot_cluster) 
p <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])
p
ggsave('results/ora/subclass_per_clust_metab_k4.png', p, height=8, width=14)

# Plot all in one
all_plot <- enrich_results_metab %>%
  filter(enrichment > 0.25,
         pathway != "Unannotated") %>%
  group_by(cluster) %>%
  slice_max(enrichment, n = 10) %>%
  ungroup() %>%
  ggplot(aes(x = factor(cluster), y = pathway)) +
  geom_point(aes(fill = enrichment, size = overlap), shape = 21) +
  scale_fill_viridis_c(name = '-log10(p-value)') +
  labs(x = 'Cluster', y = 'Pathway', title = 'Pathway Enrichment per Cluster (Metabolite)') +
  theme_bw()
all_plot
ggsave('results/for_figures/ora/subclass_all_metab_k4.png', all_plot, height=5, width=7)


##### Proteins
# For each cluster, are the proteins within it correlated with any certain pathway
distinct_prot <- clust_prot_exp %>%
  select(protein, cluster, subclasses) %>% ### CHANGE
  distinct()
all_prots <- unique(distinct_prot$protein)

# Function to run cluster enrichment
clust_enrichment <- function(clust) {
  prots_in_clust <- distinct_prot %>%
    filter(cluster == clust) %>%
    pull(protein) %>%
    unique()
  
  pathways <- unique(distinct_prot$subclasses) ### CHANGE
  
  res <- map_df(pathways, function(path) {
    prots_in_path <- distinct_prot %>%
      filter(subclasses == path) %>% ### CHANGE
      pull(protein) %>%
      unique()
    
    ### Fisher's test setup
    # prots in path in cluster     | prots not in path in cluster 
    # prots in path not in cluster | prots not in path not in cluster
    
    a <- length(intersect(prots_in_clust, prots_in_path))
    b <- length(setdiff(prots_in_clust, prots_in_path))
    c <- length(setdiff(prots_in_path, prots_in_clust))
    d <- length(setdiff(all_prots, union(prots_in_clust, prots_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    
    pval <- fisher.test(mat)$p.value
    
    tibble(cluster = clust,
           pathway = path,
           overlap = a,
           p.value = pval)
  })
  #res %>% mutate(p.adj = p.adjust(p.value, method = 'BH'),
  #               enrichment = -log10(p.adj))
  res %>% mutate(enrichment = -log10(p.value))
}

clusters <- sort(unique(distinct_prot$cluster))
enrich_results_prot <- map_df(clusters, clust_enrichment)

# Plot pathway enrichment per cluster
plot_cluster <- function(clust){
  enrich_results_prot %>%
    filter(cluster == clust,
           enrichment > 0.25,
           pathway != "Unannotated") %>%
    arrange(enrichment) %>%
    tail(10) %>%
    ggplot(aes(enrichment, reorder(pathway, enrichment))) +
    geom_point(aes(size = overlap, fill = enrichment), shape=21) +
    scale_fill_continuous() +
    labs(
      title = paste("Cluster", clust, "Pathway Enrichment"),
      x = "-log10(p-value)",
      y = "Pathway"
    ) +
    theme_bw()
}

plots <- map(clusters, plot_cluster)
p <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])
p
ggsave('results/ora/subclass_per_clust_prot_k4.png', p, height=8, width=14)

# Plot all in one
all_plot <- enrich_results_prot %>%
  filter(enrichment > 0.25,
         pathway != "Unannotated") %>%
  group_by(cluster) %>%
  slice_max(enrichment, n = 10) %>%
  ungroup() %>%
  ggplot(aes(x = factor(cluster), y = pathway)) +
  geom_point(aes(fill = enrichment, size = overlap), shape = 21) +
  scale_fill_viridis_c(name = '-log10(p-value)') +
  labs(x = 'Cluster', y = 'Pathway', title = 'Pathway Enrichment per Cluster (Protein)') +
  theme_bw()
all_plot
ggsave('results/ora/subclass_all_prot_k4.png', all_plot, height=5, width=7)
ggsave('results/ora/subclass_all_prot_k4.png', all_plot, height=5, width=7)




