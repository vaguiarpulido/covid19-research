# IsoformSwitchAnalyzeR R Package for COVID-19 Data Analysis
# Bioconductor: https://www.bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html

# Install packages needed
# if (!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
#   BiocManager::install()
#   }
# BiocManager::install("IsoformSwitchAnalyzeR")

# Load libraries
library(IsoformSwitchAnalyzeR)
library(DESeq2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggrepel)
library(grid)
library(gridExtra)
library(stringr)
library(readr)
library(plyr)
library(dplyr)
library(purrr)

# Set working directory
setwd("/")

# Import RNA-seq StringTie sample transcript data
# Contains gene ID, gene name, gene length, # exons, coverage, FPKM
# Make sure sample StringTie data separated into subdirectories
  # Verify subdirectories contain files in "t_data.ctab" format
# Need to specify read length (can obtain from external Qualimap analysis)
COVID_data <- importIsoformExpression(
  parentDir = "./STRINGTIE_DIR", 
  readLength = 150
)

# View abundance estimates for 1st 3 isoforms across samples
head(COVID_data$abundance,3)

# View read counts for 1st 3 isoforms across samples
head(COVID_data$counts, 3)

# Create design matrix containing which biological replicates/samples belong to which condition
# Extract condition names before 3nd "_"
names_split <- strsplit(colnames(COVID_data$abundance), "_")[-1] # don't include "isoform_id" col
COVID_conds <- as.character(lapply(names_split, function(x) {
  paste(x[1:3], sep="_", collapse = "_")
}))
# Create design matrix
myDesign <- data.frame(
  sampleID = colnames(COVID_data$abundance)[-1], 
  condition = COVID_conds)
myDesign

# Identify comparisons want to make b/w samples
comparisonsToMake <- as.data.frame(cbind(condition_1 = c("Series1_NHBE_Mock", "Series2_A549_Mock", "Series3_A549_Mock", "Series4_A549_Mock", "Series5_A549_Mock", "Series7_Calu3_Mock", "Series8_A549_Mock", "Series8_A549_Mock", "Series9_NHBE_Mock"),
                                         condition_2 = c("Series1_NHBE_SARS-CoV-2", "Series2_A549_SARS-CoV-2", "Series3_A549_RSV", "Series4_A549_IAV", "Series5_A549_SARS-CoV-2", "Series7_Calu3_SARS-CoV-2", "Series8_A549_HPIV3", "Series8_A549_RSV" ,"Series9_NHBE_IAV")),
                                   stringsAsFactors =F)
comparisonsToMake

# Integrate data into aSwitchList
# isoformExonAnnotation= human genome reference annotation (hg19) downloaded from GENCODE
aSwitchList <- importRdata(
  isoformCountMatrix   = COVID_data$counts,
  isoformRepExpression = COVID_data$abundance,
  designMatrix         = myDesign,
  comparisonsToMake= comparisonsToMake,
  isoformExonAnnoation = "./gencode19/gencode.v19.annotation.gtf",
  showProgress = F
)
aSwitchList



# ---------------------------------------------------------------------------



# ISOFORM USAGE & FUNCTIONAL CONSEQUENCE ANALYSES



# PART 1 of ISOFORM SWITCH ANALYSIS WORKFLOW
# Filters for non-expressed genes/isoforms, identifies isoform switches (dIFcutoff), annotates ORFs
aSwitchList <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = aSwitchList,
  dIFcutoff            = 0.3, # default
  genomeObject = BSgenome.Hsapiens.UCSC.hg19,
  pathToOutput =  '.', # output nucleotide & AA sequence FASTA files for external consequence analysis for PT2
)
aSwitchList

# View summary of isoform switching features & comparisons across samples & export
ISAP1 <- data.frame(extractSwitchSummary(aSwitchList))
ISAP1
write.csv(ISAP1, file = "./COVID19_IsoSwitchAnalysisPt1_Table.csv")

# PART 2 of ISOFORM SWITCH ANALYSIS WORKFLOW
# Import/incorporate external sequence annotation, analyze alternative splicing, predict functional consequences,
# visualize general & individual isoform switching activities
# Run CPAT, Pfam, SignalP, IUPred2 separately via webserver or locally & reimport
aSwitchList <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = aSwitchList,
  n = 10,
  dIFcutoff                 = 0.3, # same as Part1
  removeNoncodinORFs        = F,  # also want to look at non-coding regions for intron retention
  pathToCPATresultFile      = './CPAT_FILE.txt',
  pathToPFAMresultFile      = './PFAM_SCAN_FILE.txt',
  pathToSignalPresultFile   = './SIGNALIP_FILE.txt',
  pathToIUPred2AresultFile = './IUPRED2A_FILE.txt',
  codingCutoff              = 0.725, # coding potential cutoff suggested for humans 
  outputPlots               = T,
  pathToOutput              = '.'
)
aSwitchList

# Save aSwitchList
save(aSwitchList,file="./aSwitchList.rda")



# ---------------------------------------------------------------------------



# DIFFERENTIAL GENE EXPRESSION EXTERNAL DATA IMPORT



# Import DESeq2 differential gene expression results
  # Performed by collaborators & downloaded from Google Drive
DE <- list.files(path="/", pattern="*.tsv", full.names=T)
DE

# Combine results into a list and rename
DE.list <- lapply(myfiles, 
                  read.csv,
                  sep="\t", 
                  header=T)
names(DE.list) <- c("Series1_Series1_NHBE_SARS_CoV_2",
                    "Series2_A549_SARS_CoV_2",
                    "Series3_A549_RSV",
                    "Series4_A549_IAV",
                    "Series5_A549_SARS_CoV_2",
                    "Series7_Calu3_SARS_CoV_2",
                    "Series8_A549_HPIV3",
                    "Series8_A549_RSV", 
                    "Series9_NHBE_IAV")

# Subset DESeq2 adjusted pvalues & gene name
DE.list <- lapply(DE.list, function(x) {
  x[c("external_gene_name", "DESeq2_padj")]
})

# Add "condition" col w/ sample name for each DF
DE.list <- purrr::imap(my.data.2, ~mutate(.x, condition = .y))

# Reduce DE.list & add DESeq2 gene pvalues to isoformFeatures gene_q_val col
DE.DF <- Reduce(rbind, my.data.2)
aSwitchList$isoformFeatures$gene_q_value <- ifelse(match(aSwitchList$isoformFeatures$gene_name, DE.DF$external_gene_name) & 
                                                     match(aSwitchList$isoformFeatures$condition_2, DE.DF$condition), 
                                                   DE.DF$DESeq2_padj[match(aSwitchList$isoformFeatures$gene_name, DE.DF$external_gene_name)],
                                                   NA)

# Renumber rows in isoformFeatures
rownames(aSwitchList$isoformFeatures) <- 1:nrow(aSwitchList$isoformFeatures)



# ---------------------------------------------------------------------------



# FUNCTIONAL BIOLOGICAL CONSEQUENCE AND SPLICING ANALYSES



# Plot biological consequence summary histogram
extractConsequenceSummary(
  aSwitchList,
  consequencesToAnalyze='all',
  plotGenes = F,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = F      # enables analysis of fraction of significant features
)

# Plot biological consequence summary enrichment
extractConsequenceEnrichment(
  aSwitchList, 
  returnResult = F)

# Plot splicing summary histogram
extractSplicingSummary(
  aSwitchList,
  splicingToAnalyze='all',
  plotGenes = F,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = F      # enables analysis of fraction of significant features
)

# Plot splicing enrichment
extractSplicingEnrichment(
  aSwitchList, 
  returnResult = F)



# ---------------------------------------------------------------------



# INDIVIDUAL ISOFORM VISUALIZATION: volcano plots & splice graphs



# Extract isoform features from aSwitchList
aif <- aSwitchList$isoformFeatures

# Plot top 30 isos experiencing highest gene switch qval for each comparison
# Split data into comparisons (cond1 vs cond2) into a list
aif_list <- split(aif, f=aif["condition_2"])

# Extract top 30 isos per comparison
top30_isos <- lapply(aif_list, function(z) {
  # Subset sig data from gene_log2FC sig criteria
  gene_subset_sig = z[which(z$color == 'FDR < 0.05 + \nLog2FC + dIF'), ]
  # Find top 30 isos for each comparison
  top_30 = head(gene_subset_sig[order(gene_subset_sig$gene_switch_q_value),], 30)
  return(top_30)
})
top30_isos

# Re-combine top 30 isos from all comparisons into one DF
aif_top30_isos <- data.frame(do.call(rbind, top30_isos))
rownames(aif_top30_isos) <- 1:nrow(aif_top30_isos)
top_30_isos_gene_names <- as.vector(aif_top30_isos$gene_name) # make vectors for easier name search
top_30_isos_gene_names

# Create col labeling isoforms that meet sig critria
aif$top30 <- ifelse(aif$gene_name %in% top_30_isos_gene_names & aif$color == "FDR < 0.05 + \nLog2FC + dIF",
                    "Top 30","Other")
# Verify isos
which(aif$top30=="Top 30")

# Set significance criteria
  # 1) Gene log2FC < -1.3 or Gene log2FC > 1.3
  # 2) dIF > |0.3|
  # 3) Isoform switch qvalue (FDR) < 0.05
aif$color <- ifelse(aif$isoform_switch_q_value < 0.05 & (aif$gene_log2_fold_change > 1.3 | aif$gene_log2_fold_change < -1.3) & (aif$dIF > 0.3 | aif$dIF < -0.3), "FDR < 0.05 + \nLog2FC + dIF", 
                    ifelse(aif$isoform_switch_q_value < 0.05 & (aif$dIF > 0.3 | aif$dIF < -0.3), "FDR < 0.05 + dIF","Not Sig"))

# Plot top 30 isos w/ sig criteria + gene names
p <- ggplot(data=aif, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(aes(color=color),size=2) + 
  coord_cartesian(xlim = c(-10, 10)) +
  geom_hline(yintercept = -0.3, linetype='dashed') +
  geom_hline(yintercept = 0.3, linetype='dashed') +
  geom_vline(xintercept = -1.3, linetype='dashed') +
  geom_vline(xintercept = 1.3, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switching', values = c('blue','red','gray')) +
  labs(x='Gene log2 fold change', y='dIF') +
  ggtitle("Top 30 Isoforms with Significant Gene Switch Usage") +
  theme_bw(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right", legend.text = element_text(size = 20)) +
  geom_text_repel(data=subset(aif, gene_name %in% top_30_isos_gene_names & color == "FDR < 0.05 + \nLog2FC + dIF"),
                  aes(label=gene_name),
                  size=5,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))

# Separate data into multiple graphs based on condition_2
# Save graph as 25" x 25"
p + facet_wrap(~ condition_2)


# SPLICE GRAPH PLOTTING (ex: IL6 gene)
# Plot IL6 gene transcript, gene expression, iso expression, iso usage
switchPlot(
  aSwitchList,
  gene='IL6',
  condition1 = 'Series1_NHBE_Mock',
  condition2 = 'Series1_NHBE_SARS_CoV_2'
)

# Plot transcripts (i.e. splice graphs) only
switchPlotTranscript(
  aSwitchList,
  gene = 'IL6',
  condition1 = 'Series1_NHBE_Mock',
  condition2 = 'Series1_NHBE_SARS_CoV_2'
)

# Plot gene expression data only
switchPlotGeneExp(
  aSwitchList,
  gene = 'IL6',
  condition1 = 'Series1_NHBE_Mock',
  condition2 = 'Series1_NHBE_SARS_CoV_2',
  logYaxis = T
)

# Plot isoform usage data only
switchPlotIsoUsage(aSwitchList,
                   gene = 'IL6',
                   condition1 = 'Series1_NHBE_Mock',
                   condition2 = 'Series1_NHBE_SARS_CoV_2',
                   addErrorbars = T,
                   confidenceIntervalErrorbars = T,
)



# ---------------------------------------------------------------------



# EXPORT DATA



# Subset and save all genes meeting sig criteria
subset_sig_genes <- aif[which(aif$color == 'FDR < 0.05 + \nLog2FC + dIF'), ] # sig criteria
names(subset_sig_genes)[names(subset_sig_genes)=="color"] <- "Significance Criteria" # rename "color" col
write.csv(subset_sig_genes, file = "./COVID_ISA_isoFeatures_results.csv")

# Subset and save top 30 sig isos isoform data
top30_save <- subset(aif_top30_isos)
names(top30_save)[names(top30_save)=="color"] <- "Significance Criteria" # rename "color" col
write.csv(top30_save, file = "./COVID_ISA_top30_sig_isos_results.csv")
