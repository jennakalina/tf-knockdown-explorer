# TF Knockdown Data Explorer

Repository for a web interface browsing a global TF knockdown dataset and the associated scripts used to complete the analysis of raw data. See paper methods for additional information on all analyses performed.

https://jennakalina.github.io/tf-knockdown-explorer/index.html

## Scripts

The directory "scripts" contains all of the R scripts used to process the raw data and generate the data used for the web interface. Methods are described in the paper.

## Web Interface Tabs

### Differential Analysis

Differential analysis was completed by calculating the log2 fold change and *p*-value of each metabolite or protein (analyte) at an RNAi level compared to the control. 

Search by TF (target) name, metabolite, or protein to see all results for that query sorted by adjusted *p*-value.

### Concordance 

Concordance meta-analysis was used to identify protein-metabolite interaction (PMI) pairs that were consistently regulated in the same direction across all samples. 

Search by a single analyte to see results for all possible pairs, browse the volcano plot to see significant results, and search by specific PMI pair to produce a quadrant plot showing which knocked-down TFs are concordant vs discordant for any given pair.

### Perturbation Clusters

Perturbation clusters were determined by *k*-means clustering based on the sample projection plot output by rCCA analysis at a TF level. 

View the clustered projection plot and click on any TF (or select from the drop-down menu) to view a heatmap of the DAMs and DAPs associated with its cluster.
