# Scripts pertaining to downstream analyses of RNAseq data


## Usage
```
git clone https://github.com/vaguiarpulido/covid19-research.git
cd covid19-research
```

### Exploratory analyses (Principal Component Analysis (PCA)) and Differential gene expression analysis

```
Rscript scripts/downstream_analyses/dataExploration_and_DEAnalysis_primary_dataset.R --data-dir data/GSE147507 --output-dir output/diff_expr_analysis/GSE147507
Rscript scripts/downstream_analyses/dataExploration_and_DEAnalysis_additional_dataset.R --data-dir data/GSE150316 --output-dir output/diff_expr_analysis/GSE150316
```

**Output**

For each comparison is created a directory containing:
- a PCA plot with the first 2 Principal Components;
- a subdirectory containing:
  - the Venn diagram for the results of DESeq2-edgeR-limmaVoom;
  - the lists of the differentially expressed genes detected in at least 2 of the 3 tools used (DESeq2-edgeR-limmaVoom), applying a statistical cut-off of 0.05 for the adjusted p-value;
  - a subdirectory called *single_tool* with the lists of the differentially expressed genes detected by each tool used (DESeq2-edgeR-limmaVoom), applying a statistical cut-off of 0.05 for the adjusted p-value;
  - a subdirectory called *unfiltered* with the lists of the differentially expressed genes detected by each tool used (DESeq2-edgeR-limmaVoom), without applying any filtering;

### Enrichment analyses

#### Gene Ontology overrepresentation analysis
    GO_COVID.R is the main script for this analysis. Ensembl_to_entrez.R, GOstats.human.R, and Load.GO.R are loaded directly inside the main script.
    Packages which should be installed: biomaRt, org.Hs.eg.db, GO.db, and GOstats.

The main script is used directly on the folders created from the differencial expression analysis as input. We used DESEq2 lists of differentially expressed genes with a cut-off of 0.05 for the adjusted p-value.

The script generates a subdirectory for each study containing one table (csv format) for each comparison tested.
  
#### Pathway analysis
This script is called with the following command:
```
Rscript --vanilla SPIA_Code.Rmd <name_of_DEG_output_from_DESeq2>
```
