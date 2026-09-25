### DA Heatmaps (bubble plots of pathways/subclasses with DA scores)
library(ggplot2)
library(dplyr)

# Load data
reg_sum_metab <- read.csv('results/diff_analysis/reg_summary_metab_disc.csv')
reg_sum_prot <- read.csv('results/diff_analysis/reg_summary_prot_disc.csv') 
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv', check.names = FALSE)
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv('raw_data/kegg_mapping/kegg_pathways_manual.csv')
########## PATHWAY LEVEL ########## 
# Unnest pathways
metab_path_all_long <- metab_path_all %>%
  mutate(pathways = strsplit(pathways, ";\\s*")) %>% 
  unnest(pathways) %>%
  rename(Pathway = pathways)
reg_sum_metab_path <- reg_sum_metab %>%
  left_join(metab_path_all_long, by = "metabolite")
# Calculate pathway scores
pathway_scores_metab <- reg_sum_metab_path %>%
  mutate(sig_RNAi = strsplit(sig_RNAi, ";")) %>% 
  unnest(sig_RNAi) %>%
  group_by(Pathway) %>%
  summarise(
    n_tot = sum(sig) + sum(Not.significant),
    total_up = sum(Upregulated),
    total_down = sum(Downregulated),
    unique_rnai = n_distinct(sig_RNAi),
    DA_score = (total_up - total_down) / n_tot,
    DF_score = (total_up + total_down) / n_tot,
    .groups = "drop"
  ) %>%
  arrange(desc(abs(DA_score)))
# Now for protein
prot_path_all_long <- prot_path_all %>%
  mutate(pathways = strsplit(pathways, ";\\s*")) %>% 
  unnest(pathways) %>%
  rename(Pathway = pathways)
reg_sum_prot_path <- reg_sum_prot %>%
  left_join(prot_path_all_long, by = "protein")
pathway_scores_prot <- reg_sum_prot_path %>%
  mutate(sig_RNAi = strsplit(sig_RNAi, ";")) %>% 
  unnest(sig_RNAi) %>%
  group_by(Pathway) %>%
  summarise(n_tot = sum(sig) + sum(Not.significant) + sum(NAs),
            total_up = sum(Upregulated),
            total_down = sum(Downregulated),
            unique_rnai = n_distinct(sig_RNAi),
            DA_score = (total_up - total_down) / n_tot,
            DF_score = (total_up + total_down) / n_tot,
            .groups = 'drop') %>%
  arrange(desc(abs(DA_score)))

# Get data for heatmap
das_per_path <- data.frame(pathway = union(pathway_scores_metab$Pathway, pathway_scores_prot$Pathway))
pathway_scores_metab <- pathway_scores_metab %>% rename(pathway = Pathway)
das_per_path <- das_per_path %>%
  left_join(pathway_scores_metab %>% select(DA_score, pathway, unique_rnai), by = 'pathway') %>%
  rename(DA_score_metabolite = DA_score, unique_rnai_metabolite = unique_rnai)
das_per_path$DA_score_metabolite[is.na(das_per_path$DA_score_metabolite)] <- 0
das_per_path$unique_rnai_metabolite[is.na(das_per_path$unique_rnai_metabolite)] <- 0
pathway_scores_prot <- pathway_scores_prot %>% rename(pathway = Pathway)
das_per_path <- das_per_path %>%
  left_join(pathway_scores_prot %>% select(DA_score, pathway, unique_rnai), by = 'pathway') %>%
  rename(DA_score_protein = DA_score, unique_rnai_protein = unique_rnai)
das_per_path$DA_score_protein[is.na(das_per_path$DA_score_protein)] <- 0
das_per_path$unique_rnai_protein[is.na(das_per_path$unique_rnai_protein)] <- 0

# Convert to long
das_long <- das_per_path %>%
  pivot_longer(cols = c(DA_score_metabolite, unique_rnai_metabolite, DA_score_protein, unique_rnai_protein),
               names_to = c('.value', 'type'),
               names_pattern = '(DA_score|unique_rnai)_(.*)')
# Order by highest absolute DA score
das_long <- das_long %>%
  group_by(pathway) %>%
  mutate(order_val = max(abs(DA_score), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!pathway %in% 'Unannotated')

### OPTIONAL: Comment/uncomment to filter down to just in metabolism class
metab_paths <- path_classes$pathway_name[path_classes$class == 'Metabolism']
das_long <- das_long %>%
  filter(pathway %in% metab_paths)

# Bubble heatmap
ggplot(das_long, aes(x = type, y = reorder(pathway, order_val))) +
  geom_point(aes(size = unique_rnai, fill = DA_score),
             shape = 21, color = 'black', stroke = 0.5) +
  scale_fill_gradient2(low = 'royalblue', mid = 'white', high = 'tomato3', midpoint = 0) +
  scale_size_continuous(range = c(1, 5)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = 'gray90'),
        panel.grid.minor = element_blank()) +
  labs(x = '', y = 'Pathway', title = 'Heatmap of DA Scores per Pathway',
       fill = 'DA score', size = 'Number of Significant Perturbations')
ggsave('results/diff_analysis/da_heatmaps/all_pathways.png', height=12, width=8)

# Filter down; only 68 pathways are present in both
paths_keep <- intersect(pathway_scores_metab$pathway, pathway_scores_prot$pathway)
paths_keep <- paths_keep[c(1:34,36:67)] # Cut out Butanoate metabolism; has a zero value 
das_long_intersect <- das_long %>% 
  filter(pathway %in% paths_keep) %>%
  group_by(pathway) %>%
  mutate(Synchronicity = ifelse(sign(DA_score[type == 'metabolite']) == sign(DA_score[type == 'protein']),
                                'Synchronous', 'Asynchronous')) %>%
  ungroup()

### Comment/uncomment to filter down to just in metabolism class
das_long_intersect <- das_long_intersect %>%
  filter(pathway %in% metab_paths)

ggplot(das_long_intersect, aes(x = type, y = reorder(pathway, order_val), group = pathway)) +
  geom_line(aes(linetype = Synchronicity), linewidth = 0.3, alpha = 0.6) +
  scale_linetype_manual(values = c('Synchronous' = 'solid',
                                   'Asynchronous' = 'dashed')) +
  geom_point(aes(size = unique_rnai, fill = DA_score),
             shape = 21, color = 'black', stroke = 0.5) +
  scale_fill_gradient2(low = 'royalblue', mid = 'white', high = 'tomato3', midpoint = 0) +
  scale_size_continuous(range = c(1, 6)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = 'gray90'),
        panel.grid.minor = element_blank()) +
  labs(x = '', y = 'Pathway', title = 'Heatmap of DA Scores per Pathway',
       fill = 'DA score', size = 'Number of Significant Perturbations')
ggsave('results/diff_analysis/da_heatmaps/intersect_pathways.png', height=10, width=8)
ggsave('results/for_figures/diff_analysis/intersect_pathways.png', height=10, width=8)


########## SUBCLASS LEVEL ########## 
# Unnest pathways
metab_path_all_long <- metab_path_all %>%
  mutate(pathways = strsplit(subclasses, ";\\s*")) %>% 
  unnest(pathways) %>%
  rename(Pathway = pathways)
reg_sum_metab_path <- reg_sum_metab %>%
  left_join(metab_path_all_long, by = "metabolite")
# Calculate pathway scores
pathway_scores_metab <- reg_sum_metab_path %>%
  mutate(sig_RNAi = strsplit(sig_RNAi, ";")) %>% 
  unnest(sig_RNAi) %>%
  group_by(Pathway) %>%
  summarise(
    n_tot = sum(sig) + sum(Not.significant),
    total_up = sum(Upregulated),
    total_down = sum(Downregulated),
    unique_rnai = n_distinct(sig_RNAi),
    DA_score = (total_up - total_down) / n_tot,
    DF_score = (total_up + total_down) / n_tot,
    .groups = "drop"
  ) %>%
  arrange(desc(abs(DA_score)))
# Now for protein
prot_path_all_long <- prot_path_all %>%
  mutate(pathways = strsplit(subclasses, ";\\s*")) %>% 
  unnest(pathways) %>%
  rename(Pathway = pathways)
reg_sum_prot_path <- reg_sum_prot %>%
  left_join(prot_path_all_long, by = "protein")
pathway_scores_prot <- reg_sum_prot_path %>%
  mutate(sig_RNAi = strsplit(sig_RNAi, ";")) %>% 
  unnest(sig_RNAi) %>%
  group_by(Pathway) %>%
  summarise(n_tot = sum(sig) + sum(Not.significant) + sum(NAs),
            total_up = sum(Upregulated),
            total_down = sum(Downregulated),
            unique_rnai = n_distinct(sig_RNAi),
            DA_score = (total_up - total_down) / n_tot,
            DF_score = (total_up + total_down) / n_tot,
            .groups = 'drop') %>%
  arrange(desc(abs(DA_score)))

# Get data for heatmap
das_per_path <- data.frame(pathway = union(pathway_scores_metab$Pathway, pathway_scores_prot$Pathway))
pathway_scores_metab <- pathway_scores_metab %>% rename(pathway = Pathway)
das_per_path <- das_per_path %>%
  left_join(pathway_scores_metab %>% select(DA_score, pathway, unique_rnai), by = 'pathway') %>%
  rename(DA_score_metabolite = DA_score, unique_rnai_metabolite = unique_rnai)
das_per_path$DA_score_metabolite[is.na(das_per_path$DA_score_metabolite)] <- 0
das_per_path$unique_rnai_metabolite[is.na(das_per_path$unique_rnai_metabolite)] <- 0
pathway_scores_prot <- pathway_scores_prot %>% rename(pathway = Pathway)
das_per_path <- das_per_path %>%
  left_join(pathway_scores_prot %>% select(DA_score, pathway, unique_rnai), by = 'pathway') %>%
  rename(DA_score_protein = DA_score, unique_rnai_protein = unique_rnai)
das_per_path$DA_score_protein[is.na(das_per_path$DA_score_protein)] <- 0
das_per_path$unique_rnai_protein[is.na(das_per_path$unique_rnai_protein)] <- 0

# Convert to long
das_long <- das_per_path %>%
  pivot_longer(cols = c(DA_score_metabolite, unique_rnai_metabolite, DA_score_protein, unique_rnai_protein),
               names_to = c('.value', 'type'),
               names_pattern = '(DA_score|unique_rnai)_(.*)')
# Order by highest absolute DA score
das_long <- das_long %>%
  group_by(pathway) %>%
  mutate(order_val = max(abs(DA_score), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!pathway %in% 'Unannotated')

### Comment/uncomment to filter down to metabolism subclasses
metab_subcl <- unique(path_classes$subclass[path_classes$class == 'Metabolism'])
das_long <- das_long %>%
  filter(pathway %in% metab_subcl) %>%
  filter(!pathway == 'Global and overview maps')

# Bubble heatmap
ggplot(das_long, aes(x = type, y = reorder(pathway, order_val))) +
  geom_point(aes(size = unique_rnai, fill = DA_score),
             shape = 21, color = 'black', stroke = 0.5) +
  scale_fill_gradient2(low = 'royalblue', mid = 'white', high = 'tomato3', midpoint = 0) +
  scale_size_continuous(range = c(1, 5)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = 'gray90'),
        panel.grid.minor = element_blank()) +
  labs(x = '', y = 'Pathway Subclass', title = 'Heatmap of DA Scores per Subclass',
       fill = 'DA score', size = 'Number of Significant Perturbations')
ggsave('results/diff_analysis/da_heatmaps/all_subclasses.png', height=5, width=8)

# Filter down to subclasses in both
paths_keep <- intersect(pathway_scores_metab$pathway, pathway_scores_prot$pathway)
das_long_intersect <- das_long %>% 
  filter(pathway %in% paths_keep) %>%
  group_by(pathway) %>%
  mutate(Synchronicity = ifelse(sign(DA_score[type == 'metabolite']) == sign(DA_score[type == 'protein']),
                                'Synchronous', 'Asynchronous')) %>%
  ungroup()

### Comment/uncomment to filter down to metabolism subclasses
das_long_intersect <- das_long_intersect %>%
  filter(pathway %in% metab_subcl) %>%
  filter(!pathway == 'Global and overview maps')

ggplot(das_long_intersect, aes(x = type, y = reorder(pathway, order_val), group = pathway)) +
  geom_line(aes(linetype = Synchronicity), linewidth = 0.3, alpha = 0.6) +
  scale_linetype_manual(values = c('Synchronous' = 'solid',
                                   'Asynchronous' = 'dashed')) +
  geom_point(aes(size = unique_rnai, fill = DA_score),
             shape = 21, color = 'black', stroke = 0.5) +
  scale_fill_gradient2(low = 'royalblue', mid = 'white', high = 'tomato3', midpoint = 0) +
  scale_size_continuous(range = c(1, 6)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = 'gray90'),
        panel.grid.minor = element_blank()) +
  labs(x = '', y = 'Pathway Subclass', title = 'Heatmap of DA Scores per Subclass',
       fill = 'DA score', size = 'Number of Significant Perturbations')
ggsave('results/diff_analysis/da_heatmaps/intersect_subclasses.png', height=5, width=8)
