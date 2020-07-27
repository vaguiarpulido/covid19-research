# Isoform Usage and Predictive Functional Consequence Analyses of COVID19 Samples
IsoformSwitchAnalyzeR R package: https://www.bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html

Version used in this analysis: v1.11.0

# Inputs
  - StringTie transcript quantification data
  - [CPAT](http://lilab.research.bcm.edu/cpat), [IUPred2A](https://iupred2a.elte.hu), [SignalIP-5.0](http://www.cbs.dtu.dk/services/SignalP), [PFAM Hmmscan](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan) external webtool analyses results (performed after completion of isoformSwitchAnalysisPart1())

# Outputs
- Biological consequence and splicing analyses plots
- Isoform usage vs gene expression volcano plots
- Isoform features data (i.e. gene & isoform expression; isoform usage (dIF), consequence analysis results)

# Steps
1) Import RNA-seq StringTie transcript data
2) Integrate transcript data into aSwitchList
3) Perform isoformSwitchAnalysisPart1() and isoformSwitchAnalysisPart2()
4) Import differential expression analyses from DESeq2 R package results
5) Run biological and splicing enrichment analyses and plot graphs
6) Visualize individual isoforms via volcano plots and switchPlot()
7) Export isoform features data (cross-reference gene names and their product functions using [GeneCards GeneALaCart](https://genealacart.genecards.org))
