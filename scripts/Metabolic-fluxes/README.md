# Projection of transcriptomic data on the human metabolic network

Moomin runs in matlab, the script **Moomim.m** is the main script needed to run the analysis.
**EBSeq.R** and the function **run_ebseq.R** should be run in R prior to Moomin in order to create the input files.
Several versions of the human metabolic network recon can be downloaded [here](https://www.vmh.life/#downloadview).

## Requirements: 
Moomin: Matlab, CPLEX, Cobra Toolbox

R: EBSeq, data.table

## Workflow
1. Differential expression analysis with EBSeq to obtain posterior probability of a gene being differentially expressed.
2. Create tables for Moomin input.
3. Load human metabolic network in Cobra Toolbox (Matlab)
4. Run Moomin script.
