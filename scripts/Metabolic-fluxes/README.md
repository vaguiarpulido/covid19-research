# Projection of transcriptomic data on the human metabolic network

[Moomin](https://github.com/htpusa/moomin) runs in matlab, the script **Moomim.m** is the main script needed to run the analysis.
**EBSeq.R** and the function **run_ebseq.R** should be run in R prior to Moomin in order to create the input files.
Several versions of the human metabolic network recon can be downloaded [here](https://www.vmh.life/#downloadview).

## Requirements: 
Moomin: Matlab, CPLEX, [Cobra Toolbox](https://opencobra.github.io/cobratoolbox/stable/)
R: [EBSeq](https://bioconductor.org/packages/release/bioc/html/EBSeq.html), data.table

## Workflow
1. Differential expression analysis with EBSeq in R to obtain posterior probability of a gene being differentially expressed;
2. Create tables for Moomin input in R;
3. Load human metabolic network in Cobra Toolbox in Matlab;
4. Change the model gene IDs (in entrez_gene) to ensembl_gene_id dictionary;
5. Run Moomin on the modified model.
