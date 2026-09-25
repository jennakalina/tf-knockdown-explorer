library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(patchwork)

# Read in protein and counts data
prot_data <- read.csv('processed_data/filtered_noscale_prot.csv')
prot_metadata <- read.csv('processed_data/filtered_noscale_prot_meta.csv')
cts <- read.csv('../../bulkRNA/Kerui/TF_knockdown/results/counts/normalized_deseq_counts.csv')

# Filter down proteins to just the TF perturbations of interest
tfs <- list('Blimp-1', 'kay', 'Pdp1', 'pros', 'shn', 'sim')
samps_keep <- prot_metadata$Sample[prot_metadata$target %in% tfs]

prot_data <- prot_data[prot_data$Sample %in% samps_keep,]
prot_data <- prot_data %>% 
  select(-RNAi) %>% 
  left_join(prot_metadata %>% select(Sample, target), by = 'Sample') %>%
  select(-c(Sample, plate)) 

# Average protein data per target
avg_prot_data <- prot_data %>%
  group_by(target) %>%
  summarise(across(where(is.numeric), ~mean(.x)),
            .groups = 'drop') %>%
  column_to_rownames('target') %>%
  select(-c(LAM0_DROME, LAMC_DROME)) # Two proteins without an FBgn

# Filter counts to genes that correspond to proteins and vice versa
mapf <- read.delim('raw_data/uniprot_to_fbgn_filt.tsv')
colnames(mapf) <- c('protein', 'fbgn')

# List of overlapping FBgns
proteins <- colnames(avg_prot_data)
fbgns_keep <- intersect(cts$FBGN, mapf$fbgn) 
prots_keep <- mapf$protein[mapf$fbgn %in% fbgns_keep] %>% unique()
prots_keep <- prots_keep[!prots_keep %in% c('MYSP2_DROME', 'TPM4_DROME')] # Duplicate mapping

mapf <- mapf %>% filter(fbgn %in% fbgns_keep, protein %in% prots_keep)

# Counts filtered down
avg_cts <- cts %>% 
  select(-GENEID) %>% column_to_rownames('FBGN') %>% 
  t() %>%
  as.data.frame() %>%
  select(any_of(fbgns_keep)) %>%
  rownames_to_column('sample') %>%
  mutate(condition = sub("_[0-9]+$", "", sample)) %>%
  group_by(condition) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("condition")
rownames(avg_cts) <- c('kay', 'Pdp1', 'pros', 'shn', 'sim', 'Blimp-1', 'HSD', 'WT')

# Prot data filtered down, map names to FBgn
avg_prot_data <- avg_prot_data %>% 
  select(any_of(prots_keep)) %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column('protein') %>%
  left_join(mapf, by = 'protein') %>%
  select(-protein) %>%
  column_to_rownames('fbgn') %>%
  t() %>%
  as.data.frame()

# Get DEG lists
degs <- read.csv('../../bulkRNA/Kerui/TF_knockdown/results/degs/all_tf_vs_wt_degs.csv')
cond_map <- data.frame(condition = c("kay_vs_wt","blimp_vs_wt","hsd_vs_wt","pdp_vs_wt","pro_vs_wt","shn_vs_wt","sim_vs_wt"),
                       TF = c('kay','Blimp-1','hsd','Pdp1','pros','shn','sim'))
degs <- degs %>% left_join(cond_map, by = 'condition') %>% dplyr::select(-condition)

# Plot
plots <- list()
for (i in 1:length(tfs)) {
  tf <- tfs[[i]]
  
  prot_data_sub <- avg_prot_data[tf,]
  cts_sub <- avg_cts[tf,]
  
  ### Filter down to DEGs
  fbgn_list <- degs %>% filter(TF == tf) %>% pull(FBGN) %>% unique()
  fbgn_keep <- intersect(fbgn_list, colnames(prot_data_sub))
  
  prot_data_sub <- prot_data_sub %>% select(any_of(fbgn_keep))
  cts_sub <- cts_sub %>% select(any_of(fbgn_keep))
  
  plot_df <- prot_data_sub %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column('fbgn') %>% 
    rename(prot_ct = as.character(tf)) %>%
    left_join(cts_sub %>% t() %>% as.data.frame() %>% 
                rownames_to_column('fbgn') %>%
                rename(RNA_ct = as.character(tf)), by = 'fbgn') %>%
    mutate(RNA_ct_log = log10(RNA_ct)) 

  p <- ggplot(plot_df, aes(x = RNA_ct_log, y = prot_ct)) + 
    geom_point() + 
    geom_smooth(method = 'lm', formula = 'y ~ x', fullrange = TRUE) +
    theme_bw() +
    xlim(c(0, 6)) + ylim(c(0, 9)) +
    labs(title = paste0('mRNA-Protein Correlation of ', tf, ' Perturbation'),
         x = expression('log10(mRNA Count)'), y = 'log10(Protein Intensity)')
  
  m <- lm(prot_ct ~ RNA_ct_log, plot_df)
  r2 <- format(summary(m)$r.squared, digits = 3)
  r2_lab <- paste0('R^2 == ', r2)
  
  r <- cor(plot_df$RNA_ct_log, plot_df$prot_ct,
           method = "spearman",
           use = "complete.obs")
  
  r_lab <- paste0("italic(r) == ", format(r, digits = 3))
  
  p1 <- p + geom_text(x = 0.5, y = 8, label = r2_lab, parse = TRUE)
  
  plots[[i]] <- p1
}

combined_plot <- wrap_plots(plots, ncol = 3)
ggsave('results/gene_protein_corr/correlations_per_tf_log10_deseq_r2.png', combined_plot, width = 16, height = 10)


### Add annotations
library(AnnotationDbi)
library(org.Dm.eg.db)
library(KEGGREST)

mapf <- read.csv('results/gene_protein_corr/path_maps.csv')

ids <- clusterProfiler::bitr(fbgns_keep, fromType = 'FLYBASE', toType = 'FLYBASECG', OrgDb = org.Dm.eg.db)
tab <- limma::getGeneKEGGLinks(species = "dme")
tab$GeneID <- gsub("Dmel_", "", tab$GeneID)

paths <- keggList("pathway", "dme")
path_df <- data.frame(kegg = sub("path:", "", names(paths)),
                      pathway = unname(paths))

fbgn_to_path <- ids %>% 
  dplyr::rename(GeneID = FLYBASECG) %>% 
  left_join(tab) %>%
  left_join(path_df %>% dplyr::rename(PathwayID = kegg)) %>%
  mutate(pathway = sub(' - Drosophila melanogaster \\(fruit fly\\)', '', pathway),
         pathway = ifelse(is.na(pathway), 'Unannotated', pathway)) %>%
  left_join(mapf %>% dplyr::rename(pathway = Pathway)) %>%
  dplyr::select(FLYBASE, Subclass) %>% distinct() %>%
  dplyr::rename(pathway = Subclass)

### BY PATHWAY
pathways_keep <- unique(fbgn_to_path$pathway)
pathways_keep <- pathways_keep[!pathways_keep %in% c('Unannotated', 'Glycan biosynthesis and metabolism',
                                                     'Biosynthesis of other secondary metabolites')]
pathway_cols <- c("Amino acid metabolism" = "steelblue",
                  "Metabolism of cofactors and vitamins" = "orange2",
                  "Carbohydrate metabolism" = "salmon",
                  'Energy metabolism' = 'springgreen4',
                  'Lipid metabolism' = 'orchid3',
                  'Nucleotide metabolism' = 'royalblue3',
                  'Xenobiotics metabolism' = 'red3')

plots <- list()
for (tf in tfs) {
  prot_data_sub <- avg_prot_data[tf,]
  cts_sub <- avg_cts[tf,]
  
  # ### Filter down to DEGs
  # fbgn_list <- degs %>% filter(TF == tf) %>% pull(FBGN) %>% unique()
  # fbgn_keep <- intersect(fbgn_list, colnames(prot_data_sub))
  # 
  # prot_data_sub <- prot_data_sub %>% dplyr::select(any_of(fbgn_keep))
  # cts_sub <- cts_sub %>% dplyr::select(any_of(fbgn_keep))
  
  plot_df <- prot_data_sub %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("fbgn") %>%
    dplyr::rename(prot_ct = all_of(tf)) %>%
    left_join(cts_sub %>%
                t() %>%
                as.data.frame() %>%
                rownames_to_column("fbgn") %>%
                dplyr::rename(RNA_ct = all_of(tf)), by = "fbgn") %>%
    mutate(RNA_ct_log = log10(RNA_ct)) #%>% filter(fbgn != 'FBgn0000055')
  
  for (pw in pathways_keep) {
    # Get genes just in path
    pathway_genes <- fbgn_to_path %>%
      filter(pathway == pw) %>%
      pull(FLYBASE) %>%
      unique()
    
    plot_df2 <- plot_df %>%
      mutate(in_pathway = fbgn %in% pathway_genes)
    
    # Calculate R squared
    pathway_df <- plot_df2 %>% filter(in_pathway)
    
    m <- lm(prot_ct ~ RNA_ct_log, data = pathway_df)
    r2 <- summary(m)$r.squared
    r2_label <- paste0("R² = ", round(r2, 3))
    
    # Plot
    p <- ggplot(plot_df2, aes(RNA_ct_log, prot_ct)) +
      geom_point(color = "grey80", size = 2) +
      geom_point(data = subset(plot_df2, in_pathway),
                 color = pathway_cols[pw],
                 size = 2.5) +
      geom_smooth(data = pathway_df, aes(RNA_ct_log, prot_ct),
                  method = "lm",
                  formula = y ~ x,
                  color = "black",
                  alpha = 0.3) +
      annotate("text", x = 1, y = 8, label = r2_label) +
      theme_bw() +
      coord_cartesian(xlim = c(0,7), ylim = c(0,9)) +
      labs(title = paste(tf, "-", pw),
        subtitle = paste(sum(plot_df2$in_pathway), "genes highlighted"),
        x = expression(log[10]("mRNA Count")),
        y = expression(log[10]("Protein Intensity")))
    
    plots[[paste(tf, pw, sep = "_")]] <- p
  }
}

combined_plot <- wrap_plots(plots, ncol = 7)
ggsave('results/gene_protein_corr/correlations_annotated_subclass.png', combined_plot, width = 30, height = 24)


# All genes vs DEGs, filtered down to metabolic genes
met_genes <- fbgn_to_path %>% filter(!pathway == 'Unannotated') %>% pull(FLYBASE)

plots <- list()
for (i in 1:length(tfs)) {
  tf <- tfs[[i]]
  
  prot_data_sub <- avg_prot_data[tf,]
  cts_sub <- avg_cts[tf,]
  
  ### Filter down to metabolic genes
  prot_data_sub <- prot_data_sub %>% dplyr::select(any_of(met_genes))
  cts_sub <- cts_sub %>% dplyr::select(any_of(met_genes))

  plot_df <- prot_data_sub %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column('fbgn') %>% 
    dplyr::rename(prot_ct = as.character(tf)) %>%
    left_join(cts_sub %>% t() %>% as.data.frame() %>% 
                rownames_to_column('fbgn') %>%
                dplyr::rename(RNA_ct = as.character(tf)), by = 'fbgn') %>%
    mutate(RNA_ct_log = log10(RNA_ct)) 
  
  p <- ggplot(plot_df, aes(x = RNA_ct_log, y = prot_ct)) + 
    geom_point() + 
    geom_smooth(method = 'lm', formula = 'y ~ x', fullrange = TRUE) +
    theme_bw() +
    xlim(c(0, 7)) + ylim(c(0, 9)) +
    labs(title = paste0('mRNA-Protein Correlation of ', tf, ' - All'),
         x = expression('log10(mRNA Count)'), y = 'log10(Protein Intensity)')
  
  m <- lm(prot_ct ~ RNA_ct_log, plot_df)
  r2 <- format(summary(m)$r.squared, digits = 3)
  r2_lab <- paste0('R^2 == ', r2)
  
  r <- cor(plot_df$RNA_ct_log, plot_df$prot_ct,
           method = "spearman",
           use = "complete.obs")
  
  r_lab <- paste0("italic(r) == ", format(r, digits = 3))
  
  p1 <- p + geom_text(x = 0.5, y = 8, label = r2_lab, parse = TRUE)
  
  plots[[i]] <- p1
}

plots_degs <- list()
for (i in 1:length(tfs)) {
  tf <- tfs[[i]]
  
  prot_data_sub <- avg_prot_data[tf,]
  cts_sub <- avg_cts[tf,]
  
  ### Filter down to metabolic genes
  prot_data_sub <- prot_data_sub %>% dplyr::select(any_of(met_genes))
  cts_sub <- cts_sub %>% dplyr::select(any_of(met_genes))
  
  ### Filter down to DEGs
  fbgn_list <- degs %>% filter(TF == tf) %>% pull(FBGN) %>% unique()
  fbgn_keep <- intersect(fbgn_list, colnames(prot_data_sub))
  
  prot_data_sub <- prot_data_sub %>% dplyr::select(any_of(fbgn_keep))
  cts_sub <- cts_sub %>% dplyr::select(any_of(fbgn_keep))
  
  plot_df <- prot_data_sub %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column('fbgn') %>% 
    dplyr::rename(prot_ct = as.character(tf)) %>%
    left_join(cts_sub %>% t() %>% as.data.frame() %>% 
                rownames_to_column('fbgn') %>%
                dplyr::rename(RNA_ct = as.character(tf)), by = 'fbgn') %>%
    mutate(RNA_ct_log = log10(RNA_ct)) 
  
  p <- ggplot(plot_df, aes(x = RNA_ct_log, y = prot_ct)) + 
    geom_point() + 
    geom_smooth(method = 'lm', formula = 'y ~ x', fullrange = TRUE) +
    theme_bw() +
    xlim(c(0, 7)) + ylim(c(0, 9)) +
    labs(title = paste0('mRNA-Protein Correlation of ', tf, ' - DEGs'),
         x = expression('log10(mRNA Count)'), y = 'log10(Protein Intensity)')
  
  m <- lm(prot_ct ~ RNA_ct_log, plot_df)
  r2 <- format(summary(m)$r.squared, digits = 3)
  r2_lab <- paste0('R^2 == ', r2)
  
  r <- cor(plot_df$RNA_ct_log, plot_df$prot_ct,
           method = "spearman",
           use = "complete.obs")
  
  r_lab <- paste0("italic(r) == ", format(r, digits = 3))
  
  p1 <- p + geom_text(x = 0.5, y = 8, label = r2_lab, parse = TRUE)
  
  plots_degs[[i]] <- p1
}

all_plots <- c(plots[[1]], plots_degs[[1]], plots[[2]], plots_degs[[2]], 
               plots[[3]], plots_degs[[3]], plots[[4]], plots_degs[[4]], 
               plots[[5]], plots_degs[[5]], plots[[6]], plots_degs[[6]])
combined_plot <- wrap_plots(all_plots, ncol = 2)


ggsave('results/gene_protein_corr/correlations_per_tf_log10_metabolic_comparison.png', combined_plot, width = 10, height = 20)


