### Examining systems-level TF effects
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(stringr)
library(pheatmap)

# Read in data
dams <- read.csv('results/diff_analysis/DAMs_sig.csv')
daps <- read.csv('results/diff_analysis/DAPs_sig.csv')
metab_paths <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')
prot_paths <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')

# Filter down pathways to just metabolism subclasses
metab_paths_met <- metab_paths %>%
  select(metabolite, subclasses) %>%
  separate_rows(subclasses, sep = ';\\s*') %>%
  group_by(metabolite) %>%
  mutate(has_metab = any(str_detect(subclasses, regex("metab", ignore_case = TRUE)))) %>%
  filter((has_metab & str_detect(subclasses, regex('metab', ignore_case = TRUE))) | !has_metab) %>%
  mutate(subclasses = ifelse(!has_metab, 'Unannotated', subclasses)) %>%
  summarise(subclasses = paste0(subclasses, collapse = '; '),
            .groups = 'drop')
prot_paths_met <- prot_paths %>%
  select(protein, subclasses) %>%
  separate_rows(subclasses, sep = ';\\s*') %>%
  group_by(protein) %>%
  mutate(has_metab = any(str_detect(subclasses, regex("metab", ignore_case = TRUE)))) %>%
  filter((has_metab & str_detect(subclasses, regex('metab', ignore_case = TRUE))) | !has_metab) %>%
  mutate(subclasses = ifelse(!has_metab, 'Unannotated', subclasses)) %>%
  distinct() %>%
  summarise(subclasses = paste0(subclasses, collapse = '; '),
            .groups = 'drop') %>%
  filter(!protein == '')

# Keep RNAi line per TF that has the most significant analytes
dams_top_rnai <- dams %>%
  filter(!Target %in% c('Attp40', 'KH017-Attp2')) %>%
  count(RNAi, Target) %>%
  group_by(Target) %>%
  slice_max(n, n=1, with_ties = FALSE) %>%
  pull(RNAi)
dams <- dams %>% filter(RNAi %in% dams_top_rnai)

daps_top_rnai <- daps %>%
  rename(Target = target) %>%
  count(RNAi, Target) %>%
  group_by(Target) %>%
  slice_max(n, n=1, with_ties = FALSE) %>%
  pull(RNAi)
daps <- daps %>% filter(RNAi %in% daps_top_rnai) %>% rename(Target = target)

dams_pathways <- dams %>%
  select(Target, metabolite, Reg) %>%
  left_join(metab_paths_met, by = 'metabolite') %>%
  separate_rows(subclasses, sep = ';\\s*') %>%
  rename(pathway = subclasses)
daps_pathways <- daps %>%
  select(Target, protein, Reg) %>%
  left_join(prot_paths_met, by = 'protein') %>%
  separate_rows(subclasses, sep = ';\\s*') %>%
  rename(pathway = subclasses)

## Add analyte type before merging
dam_path_ct <- dams_pathways %>%
  group_by(Target, pathway, Reg) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(type = 'Metabolite')
dap_path_ct <- daps_pathways %>%
  group_by(Target, pathway, Reg) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(type = 'Protein')

## Combine
merged_path_ct <- bind_rows(dam_path_ct, dap_path_ct) %>%
  filter(Target != 'control',
         pathway != 'Unannotated')

## Make downregulated negative
merged_path_ct <- merged_path_ct %>%
  mutate(plot_n = ifelse(Reg == 'Downregulated', -n, n))

## Order TFs by total analytes across both types
order_df <- merged_path_ct %>%
  group_by(Target) %>%
  summarise(total_all = sum(abs(plot_n)),
            .groups = 'drop')
merged_path_ct <- merged_path_ct %>%
  left_join(order_df, by = 'Target')

# Set colors
path_levels <- unique(merged_path_ct$pathway)
cols <- RColorBrewer::brewer.pal(n = length(path_levels), name = "Set3")
names(cols) <- path_levels

# Plot
ggplot(merged_path_ct, aes(x = plot_n, y = reorder(Target, total_all), fill = pathway)) +
  geom_col(color = 'black', linewidth = 0.1) +
  geom_vline(xintercept = 0, color = 'gray25') +
  facet_wrap(~type, ncol=2, scales='free_x') +
  scale_fill_manual(values = cols) +
  labs(x = 'Number of Analytes',
       y = 'TF',
       fill = 'Pathway Subclass',
       title = 'Pathway Enrichment per TF') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = 'bold'))
ggsave('results/tf_systems/bar_count_TFlevel.png', width = 12, height = 8)

### Bar percent upregulated TF level
merged_path_pct_up <- merged_path_ct %>%
  filter(Reg == 'Upregulated') %>%
  select(-c(plot_n, total_all, Reg))
total_cts <- merged_path_pct_up %>%
  group_by(Target, type) %>%
  summarise(n_total = sum(n),
            .groups = 'drop')
merged_path_pct_up <- merged_path_pct_up %>%
  left_join(total_cts, by = c('Target', 'type')) %>%
  mutate(pct = n / n_total * 100) %>%
  left_join(order_df, by = 'Target')
ann_df <- merged_path_pct_up %>%
  select(Target, type, n_total, total_all) %>%
  distinct()

# Plot
ggplot(merged_path_pct_up, aes(x = pct, y = reorder(Target, total_all), fill = pathway)) +
  geom_col(color = 'black', linewidth = 0.1) +
  geom_text(data = ann_df, aes(x = 100, y = reorder(Target, total_all), label = n_total),
            hjust = -0.2, size = 3, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, color = 'gray25') +
  facet_wrap(~type, ncol=2, scales='free_x') +
  scale_fill_manual(values = cols) +
  labs(x = 'Pathway Makeup (%)',
       y = 'TF',
       fill = 'Pathway Subclass',
       title = 'Upregulated Pathway Enrichment per TF') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = 'bold'))
ggsave('results/tf_systems/bar_pct_upreg_TFlevel.png', width = 14, height = 10)

### Bar percent downregulated RNAi level
merged_path_pct_dn <- merged_path_ct %>%
  filter(Reg == 'Downregulated') %>%
  select(-c(plot_n, total_all, Reg))
total_cts <- merged_path_pct_dn %>%
  group_by(Target, type) %>%
  summarise(n_total = sum(n),
            .groups = 'drop')
merged_path_pct_dn <- merged_path_pct_dn %>%
  left_join(total_cts, by = c('Target', 'type')) %>%
  mutate(pct = n / n_total * 100) %>%
  left_join(order_df, by = 'Target')
ann_df <- merged_path_pct_dn %>%
  select(Target, type, n_total, total_all) %>%
  distinct()

# Plot
ggplot(merged_path_pct_dn, aes(x = pct, y = reorder(Target, total_all), fill = pathway)) +
  geom_col(color = 'black', linewidth = 0.1) +
  geom_text(data = ann_df, aes(x = 100, y = reorder(Target, total_all), label = n_total),
            hjust = -0.2, size = 3, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, color = 'gray25') +
  facet_wrap(~type, ncol=2, scales='free_x') +
  scale_fill_manual(values = cols) +
  labs(x = 'Pathway Makeup (%)',
       y = 'TF',
       fill = 'Pathway Subclass',
       title = 'Downregulated Pathway Enrichment per TF') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = 'bold'))
ggsave('results/tf_systems/bar_pct_downreg_TFlevel.png', width = 14, height = 10)

# Same but switching around regulation direction and type
### Bar percent metabolite level
merged_path_pct_met <- merged_path_ct %>%
  filter(type == 'Metabolite') %>%
  select(-c(plot_n, total_all, type))
total_cts <- merged_path_pct_met %>%
  group_by(Target, Reg) %>%
  summarise(n_total = sum(n),
            .groups = 'drop')
merged_path_pct_met <- merged_path_pct_met %>%
  left_join(total_cts, by = c('Target', 'Reg')) %>%
  mutate(pct = n / n_total * 100) %>%
  left_join(order_df, by = 'Target')
ann_df <- merged_path_pct_met %>%
  select(Target, Reg, n_total, total_all) %>%
  distinct()

# Plot
ggplot(merged_path_pct_met, aes(x = pct, y = reorder(Target, total_all), fill = pathway)) +
  geom_col(color = 'black', linewidth = 0.1) +
  geom_text(data = ann_df, aes(x = 100, y = reorder(Target, total_all), label = n_total),
            hjust = -0.2, size = 3, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, color = 'gray25') +
  facet_wrap(~Reg, ncol=2, scales='free_x') +
  scale_fill_manual(values = cols) +
  labs(x = 'Pathway Makeup (%)',
       y = 'TF',
       fill = 'Pathway Subclass',
       title = 'Metabolite Pathway Enrichment per TF') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = 'bold'))
ggsave('results/tf_systems/bar_pct_met_TFlevel.png', width = 12, height = 8)

### Bar percent downregulated RNAi level
merged_path_pct_prot <- merged_path_ct %>%
  filter(type == 'Protein') %>%
  select(-c(plot_n, total_all, type))
total_cts <- merged_path_pct_prot %>%
  group_by(Target, Reg) %>%
  summarise(n_total = sum(n),
            .groups = 'drop')
merged_path_pct_prot <- merged_path_pct_prot %>%
  left_join(total_cts, by = c('Target', 'Reg')) %>%
  mutate(pct = n / n_total * 100) %>%
  left_join(order_df, by = 'Target')
ann_df <- merged_path_pct_prot %>%
  select(Target, Reg, n_total, total_all) %>%
  distinct()

# Plot
ggplot(merged_path_pct_prot, aes(x = pct, y = reorder(Target, total_all), fill = pathway)) +
  geom_col(color = 'black', linewidth = 0.1) +
  geom_text(data = ann_df, aes(x = 100, y = reorder(Target, total_all), label = n_total),
            hjust = -0.2, size = 3, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, color = 'gray25') +
  facet_wrap(~Reg, ncol=2, scales='free_x') +
  scale_fill_manual(values = cols) +
  labs(x = 'Pathway Makeup (%)',
       y = 'TF',
       fill = 'Pathway Subclass',
       title = 'Protein Pathway Enrichment per TF') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = 'bold'))
ggsave('results/tf_systems/bar_pct_prot_TFlevel.png', width = 12, height = 10)



########### Same, but RNAi level
rm(list=ls())
dams <- read.csv('results/diff_analysis/DAMs_sig.csv')
daps <- read.csv('results/diff_analysis/DAPs_sig.csv')
metab_paths <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')
prot_paths <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')

# Filter down pathways to just metabolism subclasses
metab_paths_met <- metab_paths %>%
  select(metabolite, subclasses) %>%
  separate_rows(subclasses, sep = ';\\s*') %>%
  group_by(metabolite) %>%
  mutate(has_metab = any(str_detect(subclasses, regex("metab", ignore_case = TRUE)))) %>%
  filter((has_metab & str_detect(subclasses, regex('metab', ignore_case = TRUE))) | !has_metab) %>%
  mutate(subclasses = ifelse(!has_metab, 'Unannotated', subclasses)) %>%
  summarise(subclasses = paste0(subclasses, collapse = '; '),
            .groups = 'drop')
prot_paths_met <- prot_paths %>%
  select(protein, subclasses) %>%
  separate_rows(subclasses, sep = ';\\s*') %>%
  group_by(protein) %>%
  mutate(has_metab = any(str_detect(subclasses, regex("metab", ignore_case = TRUE)))) %>%
  filter((has_metab & str_detect(subclasses, regex('metab', ignore_case = TRUE))) | !has_metab) %>%
  mutate(subclasses = ifelse(!has_metab, 'Unannotated', subclasses)) %>%
  distinct() %>%
  summarise(subclasses = paste0(subclasses, collapse = '; '),
            .groups = 'drop') %>%
  filter(!protein == '')

# Filter out controls
dams <- dams %>%
  filter(!Target %in% c('Attp40', 'KH017-Attp2')) 

dams_pathways <- dams %>%
  select(RNAi, metabolite, Reg) %>%
  left_join(metab_paths_met, by = 'metabolite') %>%
  separate_rows(subclasses, sep = ';\\s*') %>%
  rename(pathway = subclasses)
daps_pathways <- daps %>%
  select(RNAi, protein, Reg) %>%
  left_join(prot_paths_met, by = 'protein') %>%
  separate_rows(subclasses, sep = ';\\s*') %>%
  rename(pathway = subclasses)

## Add analyte type before merging
dam_path_ct <- dams_pathways %>%
  group_by(RNAi, pathway, Reg) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(type = 'Metabolite')
dap_path_ct <- daps_pathways %>%
  group_by(RNAi, pathway, Reg) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(type = 'Protein')

## Combine
merged_path_ct <- bind_rows(dam_path_ct, dap_path_ct) %>%
  filter(pathway != 'Unannotated')

## Make downregulated negative
merged_path_ct <- merged_path_ct %>%
  mutate(plot_n = ifelse(Reg == 'Downregulated', -n, n))

## Order TFs by total analytes across both types
order_df <- merged_path_ct %>%
  group_by(RNAi) %>%
  summarise(total_all = sum(abs(plot_n)),
            .groups = 'drop')
merged_path_ct <- merged_path_ct %>%
  left_join(order_df, by = 'RNAi')

# Set colors
path_levels <- unique(merged_path_ct$pathway)
cols <- RColorBrewer::brewer.pal(n = length(path_levels), name = "Set3")
names(cols) <- path_levels

# Plot
ggplot(merged_path_ct, aes(x = plot_n, y = reorder(RNAi, total_all), fill = pathway)) +
  geom_col(color = 'black', linewidth = 0.1) +
  geom_vline(xintercept = 0, color = 'gray25') +
  facet_wrap(~type, ncol=2, scales='free_x') +
  scale_fill_manual(values = cols) +
  labs(x = 'Number of Analytes',
       y = 'TF',
       fill = 'Pathway Subclass',
       title = 'Pathway Enrichment per RNAi') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = 'bold'))
ggsave('results/tf_systems/bar_count_RNAiLevel.png', width = 14, height = 14)

### Bar percent upregulated RNAi level
merged_path_pct_up <- merged_path_ct %>%
  filter(Reg == 'Upregulated') %>%
  select(-c(plot_n, total_all, Reg))
total_cts <- merged_path_pct_up %>%
  group_by(RNAi, type) %>%
  summarise(n_total = sum(n),
            .groups = 'drop')
merged_path_pct_up <- merged_path_pct_up %>%
  left_join(total_cts, by = c('RNAi', 'type')) %>%
  mutate(pct = n / n_total * 100) %>%
  left_join(order_df, by = 'RNAi')
ann_df <- merged_path_pct_up %>%
  select(RNAi, type, n_total, total_all) %>%
  distinct()

# Plot
ggplot(merged_path_pct_up, aes(x = pct, y = reorder(RNAi, total_all), fill = pathway)) +
  geom_col(color = 'black', linewidth = 0.1) +
  geom_text(data = ann_df, aes(x = 100, y = reorder(RNAi, total_all), label = n_total),
            hjust = -0.2, size = 3, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, color = 'gray25') +
  facet_wrap(~type, ncol=2, scales='free_x') +
  scale_fill_manual(values = cols) +
  labs(x = 'Pathway Makeup (%)',
       y = 'RNAi',
       fill = 'Pathway Subclass',
       title = 'Upregulated Pathway Enrichment per RNAi') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = 'bold'))
ggsave('results/tf_systems/bar_pct_upreg_RNAiLevel.png', width = 14, height = 14)

### Bar percent downregulated RNAi level
merged_path_pct_dn <- merged_path_ct %>%
  filter(Reg == 'Downregulated') %>%
  select(-c(plot_n, total_all, Reg))
total_cts <- merged_path_pct_dn %>%
  group_by(RNAi, type) %>%
  summarise(n_total = sum(n),
            .groups = 'drop')
merged_path_pct_dn <- merged_path_pct_dn %>%
  left_join(total_cts, by = c('RNAi', 'type')) %>%
  mutate(pct = n / n_total * 100) %>%
  left_join(order_df, by = 'RNAi')
ann_df <- merged_path_pct_dn %>%
  select(RNAi, type, n_total, total_all) %>%
  distinct()

# Plot
ggplot(merged_path_pct_dn, aes(x = pct, y = reorder(RNAi, total_all), fill = pathway)) +
  geom_col(color = 'black', linewidth = 0.1) +
  geom_text(data = ann_df, aes(x = 100, y = reorder(RNAi, total_all), label = n_total),
            hjust = -0.2, size = 3, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, color = 'gray25') +
  facet_wrap(~type, ncol=2, scales='free_x') +
  scale_fill_manual(values = cols) +
  labs(x = 'Pathway Makeup (%)',
       y = 'RNAi',
       fill = 'Pathway Subclass',
       title = 'Downregulated Pathway Enrichment per RNAi') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = 'bold'))
ggsave('results/tf_systems/bar_pct_downreg_RNAiLevel.png', width = 14, height = 14)


### HEATMAPS - TF Level
rm(list=ls())
dams <- read.csv('results/diff_analysis/DAMs_all_discovery.csv')
daps <- read.csv('results/diff_analysis/DAPs_all_discovery.csv')
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv', check.names = FALSE)
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv('raw_data/kegg_mapping/kegg_pathways_manual.csv')

# Take top RNAi per df
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

# Add pathway data to DAMs/DAPs
dams <- dams %>% 
  select(metabolite, Target, Reg) %>%
  left_join(metab_path_all_long, by = 'metabolite')
da_scores_met <- dams %>%
  group_by(Target, subclass) %>%
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass)) %>%
  mutate(da_score = (n_up - n_down) / n_total)
da_scores_met_wider <- da_scores_met %>%
  select(Target, subclass, da_score) %>%
  pivot_wider(names_from = subclass, values_from = da_score) %>%
  tibble::column_to_rownames('Target')

daps <- daps %>% 
  select(protein, Target, Reg) %>%
  left_join(prot_path_all_long, by = 'protein')
da_scores_prot <- daps %>%
  group_by(Target, subclass) %>%
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass)) %>%
  mutate(da_score = (n_up - n_down) / n_total)
da_scores_prot_wider <- da_scores_prot %>%
  select(Target, subclass, da_score) %>%
  pivot_wider(names_from = subclass, values_from = da_score) %>%
  tibble::column_to_rownames('Target') %>%
  mutate(across(where(is.numeric), coalesce, 0))

# Plot
p <- pheatmap(da_scores_met_wider,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_rows = TRUE,
         clustering_distance_rows = 'correlation',
         cluster_cols = FALSE,
         angle_col = 45,
         main = "Pathway per TF DA Scores (Metabolite Data)")
ggsave('results/tf_systems/heatmap_metab_TFlevel.png', p, width = 8, height = 10)

p2 <- pheatmap(da_scores_prot_wider,
              color = colorRampPalette(c("blue", "white", "red"))(100),
              breaks = seq(-1, 1, length.out = 101),
              show_rownames = TRUE,
              show_colnames = TRUE,
              fontsize_row = 8,
              fontsize_col = 8,
              cluster_rows = TRUE,
              clustering_distance_rows = 'correlation',
              cluster_cols = FALSE,
              angle_col = 45,
              main = "Pathway per TF DA Scores (Protein Data)")
ggsave('results/tf_systems/heatmap_prot_TFlevel.png', p2, width = 8, height = 14)

### SAME but RNAi level
dams <- read.csv('results/diff_analysis/DAMs_all_discovery.csv')
daps <- read.csv('results/diff_analysis/DAPs_all_discovery.csv') %>% rename(Target = target)

# Add pathway data to DAMs/DAPs
dams <- dams %>% 
  select(metabolite, RNAi, Reg) %>%
  left_join(metab_path_all_long, by = 'metabolite')
da_scores_met <- dams %>%
  group_by(RNAi, subclass) %>%
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass)) %>%
  mutate(da_score = (n_up - n_down) / n_total)
da_scores_met_wider <- da_scores_met %>%
  select(RNAi, subclass, da_score) %>%
  pivot_wider(names_from = subclass, values_from = da_score) %>%
  tibble::column_to_rownames('RNAi')

daps <- daps %>% 
  select(protein, RNAi, Reg) %>%
  left_join(prot_path_all_long, by = 'protein')
da_scores_prot <- daps %>%
  group_by(RNAi, subclass) %>%
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass)) %>%
  mutate(da_score = (n_up - n_down) / n_total)
da_scores_prot_wider <- da_scores_prot %>%
  select(RNAi, subclass, da_score) %>%
  pivot_wider(names_from = subclass, values_from = da_score) %>%
  tibble::column_to_rownames('RNAi') %>%
  mutate(across(where(is.numeric), coalesce, 0))

# Plot
p3 <- pheatmap(da_scores_met_wider,
              color = colorRampPalette(c("blue", "white", "red"))(100),
              breaks = seq(-1, 1, length.out = 101),
              show_rownames = TRUE,
              show_colnames = TRUE,
              fontsize_row = 8,
              fontsize_col = 8,
              cluster_rows = TRUE,
              clustering_distance_rows = 'correlation',
              cluster_cols = FALSE,
              angle_col = 45,
              main = "Pathway per RNAi DA Scores (Metabolite Data)")
ggsave('results/tf_systems/heatmap_metab_RNAiLevel.png', p3, width = 8, height = 12)

p4 <- pheatmap(da_scores_prot_wider,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-1, 1, length.out = 101),
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize_row = 8,
               fontsize_col = 8,
               cluster_rows = TRUE,
               clustering_distance_rows = 'euclidean',
               cluster_cols = FALSE,
               angle_col = 45,
               main = "Pathway per RNAi DA Scores (Protein Data)")
ggsave('results/tf_systems/heatmap_prot_RNAiLevel.png', p4, width = 8, height = 16)

