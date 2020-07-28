### Load annotation from biomart prior to running GO.human.list.R

library(biomaRt)
library(org.Hs.eg.db)

#######################
#######################

### Load GO Annotations

#######################
#######################

### org.Hs.egGO is an R object that provides
### mappings between Entrez gene identifers and the GO
### identifers that they are directly associated with
entrez_object <- org.Hs.egGO

### Get the Entrez gene IDs that are mapped to a GO ID
mapped_genes <- mappedkeys(entrez_object)

### map GO IDs to Entrez IDs and convert to a list
entrez_to_go <- as.list(entrez_object[mapped_genes])

#### perform the reverse and map GO terms to Entrez gene ids
go_object <- as.list(org.Hs.egGO2EG)


