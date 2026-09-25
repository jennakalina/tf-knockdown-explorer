library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(purrr)

dams <- read.csv('results/diff_analysis/DAMs_all_discovery.csv')
daps <- read.csv('results/diff_analysis/DAPs_all_discovery.csv')
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv', check.names = FALSE)
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv('raw_data/kegg_mapping/kegg_pathways_manual.csv')

# Take top RNAi per df (don't run to do RNAi level)
dams_top_rnai <- dams %>%
  filter(Reg != 'Not significant') %>%
  filter(!Target %in% c('Attp40', 'KH017-Attp2')) %>%
  count(RNAi, Target) %>%
  group_by(Target) %>%
  slice_max(n, n=1, with_ties = FALSE) %>%
  pull(RNAi)
dams <- dams %>% filter(RNAi %in% dams_top_rnai)

daps_top_rnai <- daps %>%
  filter(Reg != 'Not significant') %>%
  rename(Target = target) %>%
  count(RNAi, Target) %>%
  group_by(Target) %>%
  slice_max(n, n=1, with_ties = FALSE) %>%
  pull(RNAi)
daps <- daps %>% filter(RNAi %in% daps_top_rnai) %>% rename(Target = target)

# Unnest pathways
valid_sub <- path_classes %>% 
  filter(class == 'Metabolism',
         subclass != 'Global and overview maps') %>% 
  pull(subclass) %>%
  unique()

metab_path_all_long <- metab_path_all %>%
  mutate(subclasses = strsplit(subclasses, ";\\s*")) %>% 
  unnest(subclasses) %>%
  rename(subclass = subclasses) %>%
  select(metabolite, subclass) %>%
  filter(subclass %in% valid_sub)
prot_path_all_long <- prot_path_all %>%
  mutate(subclasses = strsplit(subclasses, ";\\s*")) %>% 
  unnest(subclasses) %>%
  rename(subclass = subclasses) %>%
  select(protein, subclass) %>%
  filter(subclass %in% valid_sub)

# Filter down to just TF/DAM(P) information
dams_per_tf <- dams %>% filter(!Reg == 'Not significant') %>% select(metabolite, Target) %>% unique() 
daps_per_tf <- daps %>% filter(!Reg == 'Not significant') %>% select(protein, Target) %>% unique() 

# Set backgrounds
all_tfs_met <- dams_per_tf %>% pull(Target) %>% unique()
all_tfs_prot <- daps_per_tf %>% pull(Target) %>% unique()
all_subclass <- valid_sub
all_mets <- unique(dams$metabolite)
all_prots <- unique(daps$protein)

# Metabolite dataset - Fisher test
enrichment_met <- function(tf) {

  mets_in_tf <- dams_per_tf %>% filter(Target == tf) %>% pull(metabolite) %>% unique()
  
  res <- map_df(all_subclass, function(subcl) {
    mets_in_path <- metab_path_all_long %>%
      filter(subclass == subcl) %>%
      pull(metabolite) %>%
      unique()
    
    ### Fisher's test setup
    # metabs in path in TF     | metabs not in path in TF 
    # metabs in path not in TF | metabs not in path not in TF
    
    a <- length(intersect(mets_in_tf, mets_in_path))
    b <- length(setdiff(mets_in_tf, mets_in_path))
    c <- length(setdiff(mets_in_path, mets_in_tf))
    d <- length(setdiff(all_mets, union(mets_in_tf, mets_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    
    pval <- fisher.test(mat, alternative = 'greater')$p.value
    
    tibble(tf = tf,
           subclass = subcl,
           overlap = a,
           p.value = pval)
  })
  res %>% mutate(p.adj = p.adjust(p.value, method = 'BH'),
                 enrichment = -log10(p.value))
}

enrich_res_met <- map_df(all_tfs_met, enrichment_met)

# Order by cluster
clust_mem <- read.csv('results/corr_analysis/clust_mem/cluster_mem_tf_metab_k7.csv')
clust_mem <- clust_mem %>% 
  filter(rnai %in% unique(enrich_res_met$tf)) %>%
  mutate(order = seq(1, 59),
         label = paste0(rnai, ' (Cluster ', cluster, ')'))
enrich_res_met <- enrich_res_met %>% left_join(clust_mem %>% rename(tf = rnai), by = 'tf')

ggplot(enrich_res_met %>% filter(p.value < 0.1), aes(x = subclass, y = reorder(label, desc(order)))) +
  geom_point(aes(fill = enrichment, size = overlap), shape = 21) +
  scale_fill_viridis_c(name = '-log10(p-value)') +
  labs(x = 'Pathway Subclass', y = 'TF', title = 'Pathway Enrichment per TF (Metabolite Data)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
ggsave('results/tf_systems/fisher_tf_subclass/metabolite_res_relaxed_ClusterLabel.png', width = 6, height = 7)

# Now for protein
enrichment_prot <- function(tf) {
  
  prots_in_tf <- daps_per_tf %>% filter(Target == tf) %>% pull(protein) %>% unique()
  
  res <- map_df(all_subclass, function(subcl) {
    prots_in_path <- prot_path_all_long %>%
      filter(subclass == subcl) %>%
      pull(protein) %>%
      unique()
    
    ### Fisher's test setup
    # prots in path in TF     | prots not in path in TF 
    # prots in path not in TF | prots not in path not in TF
    
    a <- length(intersect(prots_in_tf, prots_in_path))
    b <- length(setdiff(prots_in_tf, prots_in_path))
    c <- length(setdiff(prots_in_path, prots_in_tf))
    d <- length(setdiff(all_prots, union(prots_in_tf, prots_in_path)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    
    pval <- fisher.test(mat, alternative = 'greater')$p.value
    
    tibble(tf = tf,
           subclass = subcl,
           overlap = a,
           p.value = pval)
  })
  res %>% mutate(p.adj = p.adjust(p.value, method = 'BH'),
                 enrichment = -log10(p.value))
}

enrich_res_prot <- map_df(all_tfs_prot, enrichment_prot)

ggplot(enrich_res_prot %>% filter(p.value < 0.05), aes(x = subclass, y = tf)) +
  geom_point(aes(fill = enrichment, size = overlap), shape = 21) +
  scale_fill_viridis_c(name = '-log10(p-value)') +
  labs(x = 'Pathway Subclass', y = 'TF', title = 'Pathway Enrichment per TF (Protein Data)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
ggsave('results/tf_systems/fisher_tf_subclass/protein_res.png', width = 6, height = 9)




