GO.human.list<-function(list,outfile){
  ### Gene Ontology Analysis with Human Datasets
  
  out<<-as.data.frame(matrix(ncol = 8, nrow = 0))
  
  if (length(list)>0){
    library(biomaRt)
    library(org.Hs.eg.db)
    library(GO.db)
    library(GOstats)
    
    #######################
    #######################
    
    ### Load ensembl_to_entrez script
    ### ensembl_2_entrez_gene.R
    
    ### Load GO Annotations
    ### Load.GO.R
    
    #######################
    #######################
    
    ### GO testing
    
    #######################
    #######################
 
    df<-subset(ensembl_to_entrez, ensembl_gene_id %in% list)
    df_entrez<-as.character(na.omit(df[,2]))
       
    # Use all the genes with GO terms as the universe
    universe <- mapped_genes
    paramsBP <- new('GOHyperGParams',
                    geneIds = df_entrez,
                    universeGeneIds = universe,
                    ontology = 'BP',
                    pvalueCutoff = 0.05,
                    conditional = FALSE,
                    testDirection = 'over',
                    annotation = "org.Hs.eg.db")
    
    paramsMF <- new('GOHyperGParams',
                    geneIds = df_entrez,
                    universeGeneIds = universe,
                    ontology = 'MF',
                    pvalueCutoff = 0.05,
                    conditional = FALSE,
                    testDirection = 'over',
                    annotation = "org.Hs.eg.db")
    
    paramsCC <- new('GOHyperGParams',
                    geneIds = df_entrez,
                    universeGeneIds = universe,
                    ontology = 'CC',
                    pvalueCutoff = 0.05,
                    conditional = FALSE,
                    testDirection = 'over',
                    annotation = "org.Hs.eg.db")
    
    BP <- hyperGTest(paramsBP)
    MF <- hyperGTest(paramsMF)
    CC <- hyperGTest(paramsCC)
    
    if (nrow(summary(BP))>0){
      outBP<-as.data.frame(summary(BP))
      colnames(outBP)[1] <- "GO"
      outBP$Type<-'biological_process'
    }else{
      outBP<-as.data.frame(matrix(ncol = 8, nrow = 0))
      colnames(outBP) <- c("GO","Pvalue","OddsRatio","ExpCount","Count","Size","Term","Type")
    }
    
    if (nrow(summary(MF))>0){
      outMF<-as.data.frame(summary(MF))
      colnames(outMF)[1] <- "GO"
      outMF$Type<-'molecular_function'
    }else{
      outMF<-as.data.frame(matrix(ncol = 8, nrow = 0))
      colnames(outMF) <- c("GO","Pvalue","OddsRatio","ExpCount","Count","Size","Term","Type")
    }
    
    if (nrow(summary(CC))>0){
      outCC<-as.data.frame(summary(CC))
      colnames(outCC)[1] <- "GO"
      outCC$Type<-'cellular_compartment'
    }else{
      outCC<-as.data.frame(matrix(ncol = 8, nrow = 0))
      colnames(outCC) <- c("GO","Pvalue","OddsRatio","ExpCount","Count","Size","Term","Type")
    }
    
    out<<-do.call("rbind", list(outBP, outMF, outCC))
    
    if(missing(outfile)){
      return(out)
    }else{
      write.csv(out,file=outfile)
    }
  }
}  