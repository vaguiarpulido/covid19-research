### Get entrez_gene_id for all ensembl_gene_ids

library(biomaRt)
mart<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
              path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

### If error: use instead 
#mart<-useEnsembl(biomart="ensembl", host="grch37.ensembl.org", 
#                dataset="hsapiens_gene_ensembl",mirror="useast")

ensembl_to_entrez<-getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'),mart=mart)