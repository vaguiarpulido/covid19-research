# Genome-wide bioinformatic analyses predict key host and viral factors in SARS-CoV-2 pathogenesis

## Introduction

The novel betacoronavirus named Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) caused a worldwide pandemic (COVID-19) after initially emerging in Wuhan, China. To better understand how SARS-CoV-2 interacts with human cells, we devised a novel, comprehensive bioinformatic strategy to analyze public RNA sequencing and viral genome sequencing data, which included:
1. Detection of differentially expressed human genes along with functional and signaling pathway enrichment analyses.
2. Analysis of changes in transcript isoform expression and usage
3. Identification of transposable elements (TEs) deregulated during infection.
4. Prediction of conserved interactions between viral RNA and human RNA-binding proteins (RBPs).
5. Prediction of viral sequence variants potentially affecting host response.

<figure>
  <p align="center">
  <img src="figures/Fig1.png" width="700" align="center">
  </p>
</figure>

To our knowledge, this is the *first meta-analysis to predict host factors that play a specific role in SARS-CoV-2 pathogenesis, distinct from other respiratory viruses*. 

## Contents

This GitHub repository contains the code to reproduce the analyses included in our manuscript submission ([preprint available on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.07.28.225581v1)).

Each of the following sections contains instructions and code to replicate our analyses into specific aspects of SARS-CoV-2 host-virus interactions.

### [RNA-Seq Data Processing](https://github.com/vaguiarpulido/covid19-research/tree/master/scripts/data_processing)
Instructions and code to process RNA-seq datasets from virus-infected and control samples, including read alignment and transcript quantification.

### [Exploratory Analyses of RNA-Seq Data](https://github.com/vaguiarpulido/covid19-research/tree/master/scripts/downstream_analyses)
Instructions and code to perform PCA analysis, differential expression analysis, and GO term enrichment based on RNA-seq data.

### [Pathway Analysis](https://github.com/vaguiarpulido/covid19-research/tree/master/scripts/pathways)
Instructions and code to identify statistically significant alterations in pathway expression based on the differentially expressed genes detected by DESeq2.

### [Metabolic Projection](https://github.com/vaguiarpulido/covid19-research/tree/master/scripts/Metabolic-fluxes)
Instructions and code to project gene expression data onto the human metabolic network, and predict increased or decreased metabolic fluxes.

### [Isoform Usage Analysis](https://github.com/vaguiarpulido/covid19-research/tree/master/scripts/isoform_analysis)
Instructions and code to analyze isoform usage and predict functional consequences of isoform-level changes in virus-infected samples.

### [Transposable Element Analysis](https://github.com/vaguiarpulido/covid19-research/tree/master/scripts/TE-analysis)
Instructions and code to analyze transposable element activity in virus-infected samples.

### [RBP analysis](https://github.com/vaguiarpulido/covid19-research/tree/master/scripts/rbp)
Instructions and code to predict putative binding sites for human RBPs on viral genome sequences, and to analyze binding site conservation across multiple SARS-CoV-2 isolates.
