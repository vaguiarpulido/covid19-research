##################################################################
#### Perform GO enrichment analysis from DESeq2 tables of DEG ####
##################################################################

rm(list=ls())     

#####################################
### Load functions and libraries  ###
#####################################
library(readr)
source("load/Load.GO.R")
source("load/GOstats.human.R")
source("load/Ensembl_to_entrez.R")

#######################################
### Analysis for GSE147507 Dataset  ###
#######################################

GSE147507<-'DEA_results/GSE147507_MainDataset'
filenames_GSE147507<-list.dirs(path = 'DEA_results/GSE147507_MainDataset', full.names = FALSE, recursive = FALSE )

for (series in filenames_GSE147507){
  
  ### Load DESeq2 results  ###
  filename<-paste('DEA_results/GSE147507_MainDataset/',series,'/des_Condition/ControlVsInfected_genes.DESeq2.tsv', sep='')
  DESeq <- read_delim(filename, "\t", escape_double = FALSE, trim_ws = TRUE)
  
  ### Run GO analysis function ###
  GO.human.list(DESeq$gene_id)
  
  ### Correct p-values for multitesting and create log2OddsRatio  ###
  out$padj<-p.adjust(out$Pvalue,method = "bonferroni")
  out$log2OddsRatio<-log(out$OddsRatio,base=2)

  ### Write csv file with results ###   
  write.csv(out,file=paste("output/GSE147507",series,"_GO.csv",sep=""))
}


#######################################
### Analysis for GSE150316 Dataset  ###
#######################################

GSE150316<-'DEA_results/GSE150316_ValidationDataset/WithoutThreeOutliers'
filenames_GSE150316<-list.dirs(path = 'DEA_results/GSE150316_ValidationDataset/WithoutThreeOutliers', full.names = FALSE, recursive = FALSE )

for (series in filenames_GSE150316){
  ### Load DESeq2 results  ###
  filename<-paste('DEA_results/GSE150316_ValidationDataset/WithoutThreeOutliers/',series,'/des_Case/ControlVs',series,'_genes.DESeq2.tsv', sep='')
  DESeq <- read_delim(filename, "\t", escape_double = FALSE, trim_ws = TRUE)
  
  ### Validation dataset IDs are not in Ensembl format as the previous dataset  ###
  Ensembl_genes<-subset(ensembl_to_entrez, external_gene_name %in% DESeq$gene_id)
  
  ### Run GO analysis function ###
  GO.human.list(Ensembl_genes$ensembl_gene_id)
  
  ### Correct p-values for multitesting and create log2OddsRatio  ###
  out$padj<-p.adjust(out$Pvalue,method = "bonferroni")
  out$log2OddsRatio<-log(out$OddsRatio,base=2)
  
  ### Write csv file with results ###   
  write.csv(out,file=paste("output/GSE150316/",series,"_GO.csv",sep=""))
}