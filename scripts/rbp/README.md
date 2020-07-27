# Analysis of human RBP binding motifs in SARS-CoV-2 genomes

`rbp_analysis.R` is the main script for this analysis. `rbp_functions.R` contains custom functions that are imported by `rbp_analysis.R`.

## Requirements: 
- data.table
- TFBSTools
- seqinr
- Biostrings
- rtracklayer
- foreach
- doParallel
- readxl

## Versions for original run:
This script was run in the following environment:

R version 3.5.1 (2018-07-02)

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] readxl_1.3.1         doParallel_1.0.15    iterators_1.0.12     foreach_1.5.0        rtracklayer_1.42.2   GenomicRanges_1.34.0 GenomeInfoDb_1.18.2 
 [8] Biostrings_2.50.2    XVector_0.22.0       IRanges_2.16.0       S4Vectors_0.20.1     BiocGenerics_0.28.0  seqinr_3.6-1         TFBSTools_1.20.0    
[15] data.table_1.12.2   

## Workflow
1. Filtering the list of human RNA Binding Proteins (RBPs) and their Position Weight Matrices (PWMs) from the ATtRACT database.
2. Scanning the genome of SARS-CoV-2 and related viruses with the resulting PWMs
3. Calculating enrichment of RBP binding motifs in viral genomes and genomic regions, using randomly scrambled sequences to calculate a background distribution.
4. Calculating conservation of RBP binding motifs across SARS-CoV-2 isolates.
5. Intersecting the list of potential SARS-CoV-2 interacting RBPs with other datasets.
