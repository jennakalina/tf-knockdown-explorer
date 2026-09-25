### Additional figures
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readxl)

# Load data
metab_data <- read.csv('processed_data/filtered_metab_data.csv', check.names = FALSE)
prot_data <- read.csv('processed_data/filtered_prot_data.csv', check.names = FALSE)
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv("raw_data/kegg_mapping/kegg_pathways_manual.csv")

# Expand to get one class per row
metabolites <- colnames(metab_data %>% select(-c(Sample, RNAi)))
class_metab <- metab_path_all %>%
  filter(metabolite %in% metabolites) %>%
  separate_rows(classes, sep = ";\\s*") %>%
  mutate(classes = ifelse(classes == 'Genetic_Information_Processing', 'Genetic Information Processing', classes))
class_cts_metab <- as.data.frame(table(class_metab$classes))
class_cts_metab <- class_cts_metab %>%
  dplyr::rename(Class = Var1) %>%
  mutate(pct = round((Freq / sum(Freq) * 100), 2)) %>%
  arrange(desc(Freq)) %>%
  mutate(Class = paste0(Class, ' (', pct, '%)'),
         Label = c('M', 'EIP', 'OS', 'GIP', 'CP', 'U'))

ggplot(class_cts_metab, aes(x = '', y = Freq, fill = Class)) +
  geom_bar(stat = 'identity', width = 1, color = 'white') +
  coord_polar('y') +
  theme_void() +
  geom_label_repel(aes(label = Label), 
                   position = position_stack(vjust = 0.5), 
                   show.legend = FALSE, 
                   color = 'white',
                   size = 3,
                   force = 0.25,
                   box.padding = 0.05) +
  labs(title = ' 245 Metabolites by Class')

ggsave('results/for_figures/piecharts/metab_filt_class.png', height=4, width = 6)

# Do the same for subclasses, filtering out to only subclasses in metabolism class
subclass_metab <- metab_path_all %>%
  filter(metabolite %in% metabolites) %>%
  separate_rows(subclasses, sep = ";\\s*") %>%
  filter(grepl('Metabolism', classes))
length(unique(subclass_metab$metabolite)) # 243
# Get counts of each pathway
subclass_cts_metab <- as.data.frame(table(subclass_metab$subclasses))

# Make some dfs for filtering purposes; filter non-metabolism pathways, filter low-expression pathways to "other"
subcl_in_metab <- unique(path_classes %>% filter(class == 'Metabolism') %>% select(subclass) %>% filter(subclass != 'Global and overview maps'))
other_subclass_cat <- subclass_cts_metab %>% filter(Var1 %in% subcl_in_metab$subclass) %>% filter(Freq < 5)

subclass_cts_metab <- subclass_cts_metab %>%
  dplyr::rename(Subclass = Var1) %>%
  arrange(desc(Freq)) %>%
  filter(Subclass %in% subcl_in_metab$subclass) %>%
  add_row(Subclass = 'Other - pathways with < 5 proteins', Freq = sum(other_subclass_cat$Freq)) %>%
  filter(Freq >= 5) %>%
  mutate(pct = round((Freq / sum(Freq) * 100), 2),
         Subclass = paste0(Subclass, ' (', pct, '%)'),
         Label = c('AA', 'Lipid', 'Cofac/Vit', 'Carb', 'Nuc',  
                   'Energy', 'Glycan', 'Other'))

ggplot(subclass_cts_metab, aes(x = '', y = Freq, fill = Subclass)) +
  geom_bar(stat = 'identity', width = 1, color = 'white') +
  coord_polar('y') +
  theme_void() +
  geom_label_repel(aes(label = Label), 
                   position = position_stack(vjust = 0.5), 
                   show.legend = FALSE, 
                   color = 'white',
                   size = 3,
                   force = 0.25,
                   box.padding = 0.05) +
  labs(title = ' 240 Metabolites by Subclass (Within Metabolism Class)')  
ggsave('results/for_figures/piecharts/metab_filt_subclass.png', height=4, width=6)


##### Now Protein #####
# Expand to get one class per row
proteins <- colnames(prot_data %>% select(-c(Sample, RNAi)))
class_prot <- prot_path_all %>%
  filter(protein %in% proteins) %>%
  separate_rows(classes, sep = ";\\s*") %>%
  mutate(classes = ifelse(classes == 'Genetic_Information_Processing', 'Genetic Information Processing', classes))
class_cts_prot <- as.data.frame(table(class_prot$classes))
class_cts_prot <- class_cts_prot %>%
  dplyr::rename(Class = Var1) %>%
  arrange(desc(Freq)) %>%
  mutate(pct = round((Freq / sum(Freq) * 100), 2),
         Class = paste0(Class, ' (', pct, '%)'),
         Label = c('U', 'GIP', 'M', 'CP', 'EIP', 'OS'))

ggplot(class_cts_prot, aes(x = '', y = Freq, fill = Class)) +
  geom_bar(stat = 'identity', width = 1, color = 'white') +
  coord_polar('y') +
  theme_void() +
  geom_label(aes(label = Label), 
             position = position_stack(vjust = 0.5), 
             show.legend = FALSE, 
             color = 'white',
             linewidth = 0,
             size = 3) +
  labs(title = ' 736 Proteins by Class')

# ggsave('results/for_figures/piecharts/prot_filt_class.png', height=4, width=6)

# Do the same for subclasses, filtering out to only subclasses in GIP + metabolism class
subclass_prot <- prot_path_all %>%
  filter(protein %in% proteins) %>%
  separate_rows(subclasses, sep = ";\\s*") %>%
  filter(grepl('Metabolism', classes))
# Get counts of each pathway
subclass_cts_prot <- as.data.frame(table(subclass_prot$subclasses))

# Make some dfs for filtering purposes; filter non-metabolism pathways, filter low-expression pathways to "other"
subcl_in_prot <- unique(path_classes %>% filter(class == 'Metabolism') %>% select(subclass) %>% filter(subclass != 'Global and overview maps'))
other_subclass_cat <- subclass_cts_prot %>% filter(Var1 %in% subcl_in_prot$subclass) %>% filter(Freq < 5)

subclass_cts_prot <- subclass_cts_prot %>%
  dplyr::rename(Subclass = Var1) %>%
  arrange(desc(Freq)) %>%
  filter(Subclass %in% subcl_in_prot$subclass) %>%
  add_row(Subclass = 'Other - pathways with < 5 proteins', Freq = sum(other_subclass_cat$Freq)) %>%
  filter(Freq >= 5) %>%
  mutate(pct = round((Freq / sum(Freq) * 100), 2),
         Subclass = paste0(Subclass, ' (', pct, '%)'),
         Label = c('Energy', 'Carb', 'AA', 'Cofac/Vit', 'Lipid', 
                   'Nuc', 'Xenobio', 'Other'))

ggplot(subclass_cts_prot, aes(x = '', y = Freq, fill = Subclass)) +
  geom_bar(stat = 'identity', width = 1, color = 'white') +
  coord_polar('y') +
  theme_void() +
  geom_label_repel(aes(label = Label), 
                   position = position_stack(vjust = 0.5), 
                   show.legend = FALSE, 
                   color = 'white',
                   size = 3,
                   force = 0.25,
                   box.padding = 0.05) +
  labs(title = ' 148 Proteins by Subclass (Within Metabolism Class)')  
# ggsave('results/for_figures/piecharts/prot_filt_subclass.png', height=4, width=6)

########## UNFILTERED
rm(list=ls())
metab_data <- read.csv("raw_data/TF-3-metadata_analysis.csv", na.strings = 0, check.names = F)
prot_data_raw <- read_excel('raw_data/TF_screen_proteomics.xlsx') %>% 
  select(c(1, 20:1229)) %>%
  tibble::column_to_rownames('Sample')
metab_path_all <- read.csv('raw_data/kegg_mapping/metab_paths_manual.csv')
prot_path_all <- read.csv('raw_data/kegg_mapping/prot_paths_manual.csv')
path_classes <- read.csv("raw_data/kegg_mapping/kegg_pathways_manual.csv")

# Get analyte names
metabolites <- unique(colnames(metab_data))
metabolites <- metabolites[!metabolites %in% c('Sample', 'Plate number')]
proteins <- colnames(prot_data_raw)

# Metabolites unfiltered
setdiff(metabolites, metab_path_all$metabolite)
metab_path_all <- metab_path_all %>%
  add_row(metabolite = 'N1-N12-Diacetylspermine', kegg_id = 'C03413',
          pathways = 'Unannotated', subclasses = 'Unannotated', classes = 'Unannotated') %>%
  add_row(metabolite = 'Cytosine', kegg_id = 'C00380', pathways = 'Pyrimidine metabolism',
          subclasses = 'Nucleotide metabolism', classes = 'Metabolism')

class_metab <- metab_path_all %>%
  filter(metabolite %in% metabolites) %>%
  separate_rows(classes, sep = ";\\s*") %>%
  mutate(classes = ifelse(classes == 'Genetic_Information_Processing', 'Genetic Information Processing', classes))
class_cts_metab <- as.data.frame(table(class_metab$classes))
class_cts_metab <- class_cts_metab %>%
  dplyr::rename(Class = Var1) %>%
  mutate(pct = round((Freq / sum(Freq) * 100), 2)) %>%
  arrange(desc(Freq)) %>%
  mutate(Class = paste0(Class, ' (', pct, '%)'),
         Label = c('M', 'EIP', 'OS', 'GIP', 'CP', 'U'))

ggplot(class_cts_metab, aes(x = '', y = Freq, fill = Class)) +
  geom_bar(stat = 'identity', width = 1, color = 'white') +
  coord_polar('y') +
  theme_void() +
  geom_label(aes(label = Label), 
             position = position_stack(vjust = 0.5), 
             show.legend = FALSE, 
             color = 'white',
             linewidth = 0,
             size = 3) +
  labs(title = ' 251 Metabolites (Unfiltered) by Class')
ggsave('results/for_figures/piecharts/prot_unfilt_class.png', height=4, width=6)

# Protein
prot_to_cg <- read.delim('raw_data/uniprot_to_fbgn.tsv')
prot_to_cg_clean <- prot_to_cg %>%
  group_by(uniprot) %>%
  summarise(cg = {
    x <- cg[!is.na(cg)] 
    if (length(x) == 0) 'Unmapped' else x[1]
  }, .groups = 'drop') %>%
  add_row(uniprot = 'LAM0_DROME', cg = 'Unmapped') %>%
  add_row(uniprot = 'LAMC_DROME', cg = 'Unmapped') %>%
  mutate(cg = ifelse(cg != 'Unmapped', paste0('dme:Dmel_', cg), cg))

compound_mapping <- data.frame(
  original_name = proteins,
  kegg_id = NA_character_,
  stringsAsFactors = FALSE
)
# Map compound names to KEGG IDs
for (i in 1:nrow(compound_mapping)) {
  compound_name <- compound_mapping$original_name[i]
  match_idx <- which(prot_to_cg_clean$uniprot == compound_name)
  if (length(match_idx) > 0) {
    compound_mapping$kegg_id[i] <- prot_to_cg_clean$cg[match_idx[1]]  # Take first match if multiple
  } else {
    compound_mapping$kegg_id[i] <- compound_name  # Keep original name if no match
  }}
compound_mapping <- tibble(protein = compound_mapping$original_name,
                           kegg_id = gsub('dme:Dmel_', '', compound_mapping$kegg_id))
# Gene to pathway mapping
gene_to_path <- KEGGREST::keggLink("pathway", "dme")
all_pathways <- KEGGREST::keggList("pathway", "dme")

pathway_lookup <- tibble(
  pathway_id = gsub("path:", "", names(all_pathways)),
  pathway_name = as.character(all_pathways)
)
pathway_lookup <- pathway_lookup  %>% 
  mutate(pathway_name = gsub(" - Drosophila.*$", "", pathway_name))

gene_path_df <- tibble(
  kegg_id = gsub("dme:Dmel_", "", names(gene_to_path)),
  pathway_id = gsub("path:", "", as.character(gene_to_path))
)

# Add pathway id data to compound ids
prot_to_path <- compound_mapping %>%
  left_join(gene_path_df, by = "kegg_id") %>%
  left_join(pathway_lookup, by = "pathway_id") %>%
  dplyr::select(protein, pathway_name) %>%
  mutate(pathway_name = ifelse(is.na(pathway_name), 'Unannotated', pathway_name))

prot_with_class <- prot_to_path %>%
  left_join(path_classes, by = 'pathway_name') %>%
  mutate(subclass = ifelse(is.na(subclass), 'Unannotated', subclass),
         class = ifelse(is.na(class), 'Unannotated', class)) %>%
  filter(!pathway_name %in% 'Metabolic pathways',
         !subclass %in% 'Global and overview maps')

# Filter down to manually selected pathways
valid_paths <- unique(path_classes$pathway_name)
valid_subcl <- unique(path_classes$subclass)
valid_class <- unique(path_classes$class)
# Filter down to manually selected pathways
prot_with_class_filt <- prot_with_class %>%
  mutate(pathway_name = ifelse(pathway_name %in% valid_paths, pathway_name, 'Unannotated'),
         subclass = ifelse(subclass %in% valid_subcl, subclass, 'Unannotated'),
         class = ifelse(class %in% valid_class, class, 'Unannotated')) %>%
  distinct() %>%
  group_by(protein) %>%  # Filter out Unannotated as long as there's another pathway there
  filter(!(pathway_name == "Unannotated" & any(pathway_name != "Unannotated"))) %>%
  ungroup()

# Get each pathway per individual protein
prot_path_summ <- prot_with_class_filt %>%
  group_by(protein) %>%
  summarise(pathways = paste(unique(pathway_name), collapse = '; '), 
            subclasses = paste(unique(subclass), collapse = '; '), 
            classes = paste(unique(class), collapse = '; '),
            .groups = 'drop') %>%
  arrange(pathways) %>%
  left_join(compound_mapping, by = 'protein') %>%
  dplyr::select(protein, kegg_id, pathways, subclasses, classes)

class_prot <- prot_path_summ %>%
  separate_rows(classes, sep = ";\\s*") %>%
  mutate(classes = ifelse(classes == 'Genetic_Information_Processing', 'Genetic Information Processing', classes))
class_cts_prot <- as.data.frame(table(class_prot$classes))
class_cts_prot <- class_cts_prot %>%
  dplyr::rename(Class = Var1) %>%
  arrange(desc(Freq)) %>%
  mutate(pct = round((Freq / sum(Freq) * 100), 2),
         Class = paste0(Class, ' (', pct, '%)'),
         Label = c('U', 'GIP', 'M', 'CP', 'EIP', 'OS'))

ggplot(class_cts_prot, aes(x = '', y = Freq, fill = Class)) +
  geom_bar(stat = 'identity', width = 1, color = 'white') +
  coord_polar('y') +
  theme_void() +
  geom_label(aes(label = Label), 
             position = position_stack(vjust = 0.5), 
             show.legend = FALSE, 
             color = 'white',
             linewidth = 0,
             size = 3) +
  labs(title = ' 1210 Proteins (Unfiltered) by Class')
ggsave('results/for_figures/piecharts/prot_unfilt_class.png', width = 6, height = 4)




