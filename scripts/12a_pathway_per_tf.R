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

# Add pathway data to DAMs
dams <- dams %>% 
  select(metabolite, Target, Reg) %>% # Switch for RNAi level
  left_join(metab_path_all_long, by = 'metabolite')
da_scores_met <- dams %>%
  group_by(Target, subclass) %>% # Switch for RNAi level
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
full_to_nick <- data.frame(subclass = unique(da_scores_met$subclass),
                           sub_short = c('AA', 'Carbohydrate', 'Energy', 'Glycan', 'Lipid', 
                                         'Cofactors/Vitamins', 'Nucleotide', 'Xenobiotics'))
da_scores_met <- da_scores_met %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

### For pairwise, skip to ###PAIRWISE
da_scores_met <- da_scores_met %>%
  mutate(pct_up = n_up / n_total,
         pct_down = -n_down / n_total)

plot_df <- da_scores_met %>%
  pivot_longer(cols = c(pct_up, pct_down),
               names_to = "direction",
               values_to = "pct") %>%
  mutate(direction = recode(direction,
                            pct_up = "Upregulated",
                            pct_down = "Downregulated"),
         # point size values
         n = ifelse(direction == "Upregulated", n_up, n_down),
         
         # keep ordering stable
         subclass = factor(subclass,
                           levels = unique(subclass)))

# Set factor labels to show up correctly
plot_df$direction <- factor(plot_df$direction,
                            levels = c("Upregulated", "Downregulated"))

p_met <- ggplot(plot_df, aes(x = direction, y = Target, color = pct)) + # Switch for RNAi level
  geom_point(aes(size=log2(n), fill = pct), shape = 21, color = 'black') +
  scale_size_continuous(range = c(0.1, 4)) +
  scale_x_discrete(guide = guide_axis(angle=90)) +
  scale_fill_gradient2(low = 'steelblue4', mid = 'white', 
                       high = 'tomato3', midpoint = 0) +
  labs(x = '', y = '', fill = 'DA Score', size = 'log2(count)', title = 'DA Scores per TF - Metabolite Data') + # Switch for RNAi level
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~subclass, nrow=1) +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text=element_text(size=7), 
        axis.title=element_text(size=7),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color="black"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        strip.text.x = element_text(angle = 90, size = 7))
p_met
ggsave('results/tf_systems/subclass_metab_TFlevel.png', width = 8, height = 9)

# Add pathway data to DAPs
daps <- daps %>% 
  mutate(Reg = ifelse(is.na(Reg), 'Not significant', Reg)) %>%
  select(protein, Target, Reg) %>% # Switch for RNAi
  left_join(prot_path_all_long, by = 'protein')
da_scores_prot <- daps %>%
  group_by(Target, subclass) %>% # Switch for RNAi
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
da_scores_prot <- da_scores_prot %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

da_scores_prot <- da_scores_prot %>%
  mutate(pct_up = n_up / n_total,
         pct_down = -n_down / n_total)

plot_df <- da_scores_prot %>%
  pivot_longer(cols = c(pct_up, pct_down),
               names_to = "direction",
               values_to = "pct") %>%
  mutate(direction = recode(direction,
                            pct_up = "Upregulated",
                            pct_down = "Downregulated"),
         # point size values
         n = ifelse(direction == "Upregulated", n_up, n_down),
         
         # keep ordering stable
         subclass = factor(subclass,
                           levels = unique(subclass)))

# Set factor labels to show up correctly
plot_df$direction <- factor(plot_df$direction,
                            levels = c("Upregulated", "Downregulated"))

p_prot <- ggplot(plot_df, aes(x = direction, y = Target, color = pct)) + # Switch for RNAi
  geom_point(aes(size=log2(n), fill = pct), shape = 21, color = 'black') +
  scale_size_continuous(range = c(0.1, 3.5)) +
  scale_x_discrete(guide = guide_axis(angle=90)) +
  scale_fill_gradient2(low = 'steelblue4', mid = 'white', 
                       high = 'tomato3', midpoint = 0) +
  labs(x = '', y = '', fill = 'DA Score', size = 'log2(count)', title = 'DA Scores per TF - Protein Data') + # Switch for RNAi
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~subclass, nrow=1) +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text=element_text(size=7), 
        axis.title=element_text(size=7),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color="black"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        strip.text.x = element_text(angle = 90, size = 7))
p_prot
ggsave('results/tf_systems/subclass_prot_RNAiLevel.png', width = 8, height = 14)


### Similar, but met v prot instead of up v down
rm(list=ls())
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

dams <- dams %>% 
  select(metabolite, Target, Reg) %>% # Switch for RNAi level
  left_join(metab_path_all_long, by = 'metabolite')
da_scores_met <- dams %>%
  group_by(Target, subclass) %>% # Switch for RNAi level
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
full_to_nick <- data.frame(subclass = unique(da_scores_met$subclass),
                           sub_short = c('AA', 'Carbohydrate', 'Energy', 'Glycan', 'Lipid', 
                                         'Cofactors/Vitamins', 'Nucleotide', 'Xenobiotics'))
da_scores_met <- da_scores_met %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

# Add pathway data to DAPs
daps <- daps %>% 
  mutate(Reg = ifelse(is.na(Reg), 'Not significant', Reg)) %>%
  select(protein, Target, Reg) %>% # Switch for RNAi
  left_join(prot_path_all_long, by = 'protein')
da_scores_prot <- daps %>%
  group_by(Target, subclass) %>% # Switch for RNAi
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
da_scores_prot <- da_scores_prot %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

# Calculate DA scores
da_scores_met <- da_scores_met %>%
  mutate(da_score_met = (n_up - n_down) / n_total) %>%
  rename(n_total_met = n_total) %>%
  select(-c(n_up, n_down))

da_scores_prot <- da_scores_prot %>%
  mutate(da_score_prot = (n_up - n_down) / n_total) %>%
  rename(n_total_prot = n_total) %>%
  select(-c(n_up, n_down))

plot_df <- da_scores_prot %>% left_join(da_scores_met, by = c('Target', 'subclass')) # Switcgh for RNAi

plot_df <- plot_df %>%
  pivot_longer(cols = c(da_score_met, da_score_prot),
               names_to = "type",
               values_to = "da_score") %>%
  mutate(type = recode(type,
                       da_score_met = "Metabolite",
                       da_score_prot = "Protein"),
         # point size values
         n = ifelse(type == "Metabolite", n_total_met, n_total_prot),
         
         # keep ordering stable
         subclass = factor(subclass,
                           levels = unique(subclass)))

ggplot(plot_df, aes(x = type, y = Target, color = pct)) + # Switch for RNAi
  geom_point(aes(size=log2(n), fill = da_score), shape = 21, color = 'black') +
  scale_size_continuous(range = c(0.1, 4)) +
  scale_x_discrete(guide = guide_axis(angle=90)) +
  scale_fill_gradient2(low = 'steelblue4', mid = 'white', 
                       high = 'tomato3', midpoint = 0) +
  labs(x = '', y = '', fill = 'DA Score', size = 'log2(count)', title = 'DA Scores per TF') + # Switch for RNAi
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~subclass, nrow=1) +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text=element_text(size=7), 
        axis.title=element_text(size=7),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color="black"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        strip.text.x = element_text(angle = 90, size = 7))
ggsave('results/tf_systems/subclass_TFlevel.png', width = 8, height = 12)




###PAIRWISE
# Make pairs
pairwise_df <- da_scores_met %>%
  group_by(Target) %>%
  reframe(pair = combn(subclass, 2, simplify = FALSE)) %>%
  mutate(subclass1 = map_chr(pair, 1),
         subclass2 = map_chr(pair, 2)) %>%
  select(-pair) %>%
  left_join(da_scores_met %>% 
              select(Target, subclass, n_up, n_down, n_total) %>%
              rename(subclass1 = subclass, n_up_1 = n_up, n_down_1 = n_down, n_total_1 = n_total),
            by = c('Target', 'subclass1')) %>%
  left_join(da_scores_met %>% 
              select(Target, subclass, n_up, n_down, n_total) %>%
              rename(subclass2 = subclass, n_up_2 = n_up, n_down_2 = n_down, n_total_2 = n_total),
            by = c('Target', 'subclass2'))

pairwise_plot_df <- pairwise_df %>%
  mutate(subclass_label = paste0(subclass1, '-', subclass2),
         n_up = n_up_1 + n_up_2,
         n_down = n_down_1 + n_down_2,
         n_total = n_total_1 + n_total_2,
         pct_up = n_up / n_total,
         pct_down = -n_down / n_total) %>%
  select(Target, subclass_label, n_up, n_down, pct_up, pct_down)

plot_df <- pairwise_plot_df %>%
  pivot_longer(cols = c(pct_up, pct_down),
               names_to = "direction",
               values_to = "pct") %>%
  mutate(direction = recode(direction,
                            pct_up = "Upregulated",
                            pct_down = "Downregulated"),
    # point size values
    n = ifelse(direction == "Upregulated", n_up, n_down),
    
    # keep ordering stable
    subclass_label = factor(subclass_label,
                            levels = unique(subclass_label)))

# Set factor labels to show up correctly
plot_df$direction <- factor(plot_df$direction,
                            levels = c("Upregulated", "Downregulated"))

ggplot(plot_df, aes(x = direction, y = Target, color = pct)) +
  geom_point(aes(size=log2(n), fill = pct), shape = 21, color = 'black') +
  scale_size_continuous(range = c(0.1, 4)) +
  scale_x_discrete(guide = guide_axis(angle=90)) +
  scale_fill_gradient2(low = 'steelblue4', mid = 'white', 
                       high = 'tomato3', midpoint = 0) +
  labs(x = '', y = '', fill = 'DA Score', size = 'log2(count)') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~subclass_label, nrow=1) +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        text=element_text(size=7), 
        axis.text=element_text(size=7), 
        axis.title=element_text(size=7),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color="black"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        strip.text.x = element_text(angle = 90, size = 7))
ggsave('results/tf_systems/path_pairs_metab_TFlevel.png', width = 18, height = 11)

### RNAi level
dams <- read.csv('results/diff_analysis/DAMs_all_discovery.csv')

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
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
da_scores_met <- da_scores_met %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

# Make pairs
pairwise_df <- da_scores_met %>%
  group_by(RNAi) %>%
  reframe(pair = combn(subclass, 2, simplify = FALSE)) %>%
  mutate(subclass1 = map_chr(pair, 1),
         subclass2 = map_chr(pair, 2)) %>%
  select(-pair) %>%
  left_join(da_scores_met %>% 
              select(RNAi, subclass, n_up, n_down, n_total) %>%
              rename(subclass1 = subclass, n_up_1 = n_up, n_down_1 = n_down, n_total_1 = n_total),
            by = c('RNAi', 'subclass1')) %>%
  left_join(da_scores_met %>% 
              select(RNAi, subclass, n_up, n_down, n_total) %>%
              rename(subclass2 = subclass, n_up_2 = n_up, n_down_2 = n_down, n_total_2 = n_total),
            by = c('RNAi', 'subclass2'))

pairwise_plot_df <- pairwise_df %>%
  mutate(subclass_label = paste0(subclass1, '-', subclass2),
         n_up = n_up_1 + n_up_2,
         n_down = n_down_1 + n_down_2,
         n_total = n_total_1 + n_total_2,
         pct_up = n_up / n_total,
         pct_down = -n_down / n_total) %>%
  select(RNAi, subclass_label, n_up, n_down, pct_up, pct_down)

plot_df <- pairwise_plot_df %>%
  pivot_longer(cols = c(pct_up, pct_down),
               names_to = "direction",
               values_to = "pct") %>%
  mutate(direction = recode(direction,
                            pct_up = "Upregulated",
                            pct_down = "Downregulated"),
         # point size values
         n = ifelse(direction == "Upregulated", n_up, n_down),
         
         # keep ordering stable
         subclass_label = factor(subclass_label,
                                 levels = unique(subclass_label)))

# Set factor labels to show up correctly
plot_df$direction <- factor(plot_df$direction,
                            levels = c("Upregulated", "Downregulated"))

ggplot(plot_df, aes(x = direction, y = RNAi, color = pct)) +
  geom_point(aes(size=log2(n), fill = pct), shape = 21, color = 'black') +
  scale_size_continuous(range = c(0.1, 4)) +
  scale_x_discrete(guide = guide_axis(angle=90)) +
  scale_fill_gradient2(low = 'steelblue4', mid = 'white', 
                       high = 'tomato3', midpoint = 0) +
  labs(x = '', y = '', fill = 'DA Score', size = 'log2(count)') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~subclass_label, nrow=1) +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        text=element_text(size=7), 
        axis.text=element_text(size=7), 
        axis.title=element_text(size=7),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color="black"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        strip.text.x = element_text(angle = 90, size = 7))
ggsave('results/tf_systems/path_pairs_metab_RNAiLevel.png', width = 18, height = 14)


### PROTEIN TF level
# Add pathway data to DAPs
daps <- daps %>% 
  mutate(Reg = ifelse(is.na(Reg), 'Not significant', Reg)) %>%
  select(protein, Target, Reg) %>%
  left_join(prot_path_all_long, by = 'protein')
da_scores_prot <- daps %>%
  group_by(Target, subclass) %>%
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
da_scores_prot <- da_scores_prot %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

# Make pairs
pairwise_df <- da_scores_prot %>%
  group_by(Target) %>%
  reframe(pair = combn(subclass, 2, simplify = FALSE)) %>%
  mutate(subclass1 = map_chr(pair, 1),
         subclass2 = map_chr(pair, 2)) %>%
  select(-pair) %>%
  left_join(da_scores_prot %>% 
              select(Target, subclass, n_up, n_down, n_total) %>%
              rename(subclass1 = subclass, n_up_1 = n_up, n_down_1 = n_down, n_total_1 = n_total),
            by = c('Target', 'subclass1')) %>%
  left_join(da_scores_prot %>% 
              select(Target, subclass, n_up, n_down, n_total) %>%
              rename(subclass2 = subclass, n_up_2 = n_up, n_down_2 = n_down, n_total_2 = n_total),
            by = c('Target', 'subclass2'))

pairwise_plot_df <- pairwise_df %>%
  mutate(subclass_label = paste0(subclass1, '-', subclass2),
         n_up = n_up_1 + n_up_2,
         n_down = n_down_1 + n_down_2,
         n_total = n_total_1 + n_total_2,
         pct_up = n_up / n_total,
         pct_down = -n_down / n_total) %>%
  select(Target, subclass_label, n_up, n_down, pct_up, pct_down)

plot_df <- pairwise_plot_df %>%
  pivot_longer(cols = c(pct_up, pct_down),
               names_to = "direction",
               values_to = "pct") %>%
  mutate(direction = recode(direction,
                            pct_up = "Upregulated",
                            pct_down = "Downregulated"),
         # point size values
         n = ifelse(direction == "Upregulated", n_up, n_down),
         
         # keep ordering stable
         subclass_label = factor(subclass_label,
                                 levels = unique(subclass_label)))

# Set factor labels to show up correctly
plot_df$direction <- factor(plot_df$direction,
                            levels = c("Upregulated", "Downregulated"))

ggplot(plot_df, aes(x = direction, y = Target, color = pct)) +
  geom_point(aes(size=log2(n), fill = pct), shape = 21, color = 'black') +
  scale_size_continuous(range = c(0.1, 3.5)) +
  scale_x_discrete(guide = guide_axis(angle=90)) +
  scale_fill_gradient2(low = 'steelblue4', mid = 'white', 
                       high = 'tomato3', midpoint = 0) +
  labs(x = '', y = '', fill = 'DA Score', size = 'log2(count)') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~subclass_label, nrow=1) +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        text=element_text(size=7), 
        axis.text=element_text(size=7), 
        axis.title=element_text(size=7),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color="black"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        strip.text.x = element_text(angle = 90, size = 7))
ggsave('results/tf_systems/path_pairs_prot_TFlevel.png', width = 18, height = 12)

# RNAi level
daps <- read.csv('results/diff_analysis/DAPs_all_discovery.csv') %>% rename(Target = target)

daps <- daps %>% 
  mutate(Reg = ifelse(is.na(Reg), 'Not significant', Reg)) %>%
  select(protein, RNAi, Reg) %>%
  left_join(prot_path_all_long, by = 'protein')
da_scores_prot <- daps %>%
  group_by(RNAi, subclass) %>%
  summarise(n_up = sum(Reg == 'Upregulated'),
            n_down = sum(Reg == 'Downregulated'),
            n_total = n(),
            .groups = 'drop') %>%
  filter(!is.na(subclass),
         !subclass %in% c('Metabolism of terpenoids and polyketides',
                          'Biosynthesis of other secondary metabolites')) 

# Add nicknames for subclasses
da_scores_prot <- da_scores_prot %>% 
  left_join(full_to_nick, by = 'subclass') %>% 
  select(-subclass) %>% 
  rename(subclass = sub_short)

# Make pairs
pairwise_df <- da_scores_prot %>%
  group_by(RNAi) %>%
  reframe(pair = combn(subclass, 2, simplify = FALSE)) %>%
  mutate(subclass1 = map_chr(pair, 1),
         subclass2 = map_chr(pair, 2)) %>%
  select(-pair) %>%
  left_join(da_scores_prot %>% 
              select(RNAi, subclass, n_up, n_down, n_total) %>%
              rename(subclass1 = subclass, n_up_1 = n_up, n_down_1 = n_down, n_total_1 = n_total),
            by = c('RNAi', 'subclass1')) %>%
  left_join(da_scores_prot %>% 
              select(RNAi, subclass, n_up, n_down, n_total) %>%
              rename(subclass2 = subclass, n_up_2 = n_up, n_down_2 = n_down, n_total_2 = n_total),
            by = c('RNAi', 'subclass2'))

pairwise_plot_df <- pairwise_df %>%
  mutate(subclass_label = paste0(subclass1, '-', subclass2),
         n_up = n_up_1 + n_up_2,
         n_down = n_down_1 + n_down_2,
         n_total = n_total_1 + n_total_2,
         pct_up = n_up / n_total,
         pct_down = -n_down / n_total) %>%
  select(RNAi, subclass_label, n_up, n_down, pct_up, pct_down)

plot_df <- pairwise_plot_df %>%
  pivot_longer(cols = c(pct_up, pct_down),
               names_to = "direction",
               values_to = "pct") %>%
  mutate(direction = recode(direction,
                            pct_up = "Upregulated",
                            pct_down = "Downregulated"),
         # point size values
         n = ifelse(direction == "Upregulated", n_up, n_down),
         
         # keep ordering stable
         subclass_label = factor(subclass_label,
                                 levels = unique(subclass_label)))

# Set factor labels to show up correctly
plot_df$direction <- factor(plot_df$direction,
                            levels = c("Upregulated", "Downregulated"))

ggplot(plot_df, aes(x = direction, y = RNAi, color = pct)) +
  geom_point(aes(size=log2(n), fill = pct), shape = 21, color = 'black') +
  scale_size_continuous(range = c(0.1, 3.5)) +
  scale_x_discrete(guide = guide_axis(angle=90)) +
  scale_fill_gradient2(low = 'steelblue4', mid = 'white', 
                       high = 'tomato3', midpoint = 0) +
  labs(x = '', y = '', fill = 'DA Score', size = 'log2(count)') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~subclass_label, nrow=1) +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        text=element_text(size=7), 
        axis.text=element_text(size=7), 
        axis.title=element_text(size=7),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color="black"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        strip.text.x = element_text(angle = 90, size = 7))
ggsave('results/tf_systems/path_pairs_prot_RNAiLevel.png', width = 18, height = 16)


