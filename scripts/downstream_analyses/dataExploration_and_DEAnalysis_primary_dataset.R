# For Linux users:
# - sudo apt-get install libcurl4-openssl-dev
# - sudo apt-get install libxml2-dev

if(!requireNamespace('optparse', quietly = TRUE)){
  install.packages('optparse', dependencies = T)
}
library("optparse")

option_list <- list(
  make_option(c("-d", "--data-dir"), type="character", default=NULL, 
              help="input dataset directory", metavar="character"),
  make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$d)){
  print_help(opt_parser)
  stop("Specify the input dataset directory", call.=FALSE)
}
if (is.null(opt$o)){
  print_help(opt_parser)
  stop("Specify the output directory", call.=FALSE)
}

for (package in c('RCurl', 'XML', 'tidyverse', 'BiocManager', 'gplots', 'ggfortify', 'hexbin')){
  if(!requireNamespace(package, quietly = TRUE)){
    install.packages(package, dependencies = T)
  }
}

for (bioconductor_package in c('DESeq2', 'apeglm', 'limma', 'edgeR')){
  print(bioconductor_package)
  print(!requireNamespace(bioconductor_package, quietly = TRUE))
  
  if(!requireNamespace(bioconductor_package, quietly = TRUE)){
    BiocManager::install(bioconductor_package, dependencies = T)
  }
}

library(tidyverse)
library(DESeq2)
library(edgeR)
library(limma)
library(gplots)
library(ggfortify)

genome_annotation <- 'GENCODE_v19'

#====================================================================================
# Paths and directories
#======================
path_metadata <- file.path('/home/guarracino/git/covid19-research/data/GSE147507/', 'metadata.tsv')
path_read_count <- file.path('/home/guarracino/git/covid19-research/data/GSE147507/', 'gene_count_matrix.csv')
path_ensembl_annotation_all <- file.path('/home/guarracino/git/covid19-research/data/GSE147507/', paste0('ensembl_annotation.all.', genome_annotation, '.txt'))
#====================================================================================

#====================================================================================
# Prepare data and metadata
#==========================
read_counts_df <- read_csv(path_read_count) %>%
  mutate(gene_id = sub(pattern = "(.*)\\..*$", replacement = "\\1", gene_id)) %>% # It removes the version
  column_to_rownames(var = "gene_id")

colnames(read_counts_df) <- gsub("^X", "",  colnames(read_counts_df)) # Fix column names starting with a number (R adds an X)

ensembl_annotation_df <- read.csv(path_ensembl_annotation_all, header = T, sep = '\t') %>%
  column_to_rownames(var = "ensembl_gene_id")

metadata_df <- read.csv(path_metadata, sep = '\t') %>% column_to_rownames(var = 'SampleCode')
metadata_df$Condition <- factor(metadata_df$Condition, levels = c("Control", "Infected"))

normalized_counts <- counts(DESeq.ds, normalized=TRUE)

ensg_to_plot <- 'ENSG00000136244'
plotCounts(
  dds = DESeq.ds,
  intgroup = c("Condition"),
  # gene = gsub('...$', '', ensg_to_plot),
  gene = ensg_to_plot,
  normalized = TRUE, transform = TRUE,
  returnData = F, # To save the output of plotCounts() to a variable
)
d <- plotCounts(
  dds = DESeq.ds,
  intgroup = c("Condition"),
  # gene = gsub('...$', '', ensg_to_plot),
  gene = ensg_to_plot,
  normalized = TRUE, transform = TRUE,
  returnData = TRUE, # To save the output of plotCounts() to a variable
)


diff_expr_analysis <- function(
  read_counts_df,
  metadata_df,
  ensembl_annotation_df,
  
  group_str,
  numerator_level_str,
  denominator_level_str,
  coef,
  
  design_formula,
  pv_adj_threshold,
  dir_output
) {
  # edgeR
  edgeR.DGElist <- DGEList(
    counts = read_counts_df,
    samples = metadata_df,
    group = metadata_df[, group_str]
  )
  
  design <- model.matrix(
    data = edgeR.DGElist$samples,
    design_formula
  )
  
  # Filtering
  to_keep_list  <- rowSums(read_counts_df >= 5) >= (min(table(metadata_df[, group_str])) / 2)
  edgeR.DGElist <- edgeR.DGElist[to_keep_list, , keep.lib.sizes=FALSE]
  
  edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
  edgeR.DGElist <- estimateDisp(edgeR.DGElist, design)
  
  edger_ql_fit <- glmQLFit(edgeR.DGElist, design)
  edger_qlf <- glmQLFTest(edger_ql_fit, coef = coef)
  
  edgeR.ql.res.sorted <- topTags(
    edger_qlf,
    n = Inf , # To retrieve all genes
    p.value = 1,
    sort.by = "PValue", adjust.method = "BH"
  )
  
  # DESeq2
  DESeq.ds <- DESeqDataSetFromMatrix(
    countData = read_counts_df,
    colData = metadata_df,
    design = design_formula
  )
  
  DESeq.ds <- DESeq.ds[to_keep_list, ]
  DESeq.ds <- DESeq(DESeq.ds)
  
  coef_deseq2 <- paste(group_str, numerator_level_str, 'vs', denominator_level_str, sep = '_')
  
  DESeq.res <- results(
    object = DESeq.ds,
    contrast = c(group_str, numerator_level_str, denominator_level_str),
    independentFiltering = TRUE,
    alpha = pv_adj_threshold,
    pAdjustMethod = "BH",
  )
  DESeq.res_lfcShrink <- lfcShrink(
    res = DESeq.res,
    dds = DESeq.ds,
    coef = coef_deseq2,
    type = "apeglm",
  )
  
  # limma-voom
  voomTransformed <- voom(edgeR.DGElist, design, plot = FALSE)
  voomed.fitted <- lmFit(voomTransformed, design = design)
  
  voomed.fitted <- eBayes(voomed.fitted)
  
  limmaV.res.sorted <- topTable(
    voomed.fitted, coef = coef,
    number = Inf,
    p.value = 1,
    sort.by = "P", adjust.method = "BH"
  )
  
  # Results
  dir_xxx <- file.path(
    dir_output,
    paste0('design_', gsub("\\ \\+\\ ", "", as.character(design_formula)[2]))
  )
  
  
  DEG_All_unfiltered <- rownames(DESeq.res_lfcShrink)
  
  DE_list <- list(
    LimmaV = rownames(subset(limmaV.res.sorted, adj.P.Val <= pv_adj_threshold)),
    edgeR = rownames(subset(edgeR.ql.res.sorted$table, FDR <= pv_adj_threshold)),
    DESeq2 = rownames(subset(DESeq.res_lfcShrink, padj <= pv_adj_threshold))
  )

  DEG_All.expFCpadj <- data.frame(
    LimmaV_AveExpr = limmaV.res.sorted[DEG_All_unfiltered, ]$AveExpr,
    edgeR_logCPM = edgeR.ql.res.sorted[DEG_All_unfiltered, ]$table$logCPM,
    DESeq2_baseMean = DESeq.res_lfcShrink[DEG_All_unfiltered, ]$baseMean,
    
    LimmaV_log2FC = limmaV.res.sorted[DEG_All_unfiltered, ]$logFC,
    edgeR_log2FC = edgeR.ql.res.sorted[DEG_All_unfiltered, ]$table$logFC,
    DESeq2_log2FC = DESeq.res_lfcShrink[DEG_All_unfiltered, ]$log2FoldChange,
    
    LimmaV_padj = limmaV.res.sorted[DEG_All_unfiltered, ]$adj.P.Val,
    edgeR_padj = edgeR.ql.res.sorted[DEG_All_unfiltered, ]$table$FDR,
    DESeq2_padj = DESeq.res_lfcShrink[DEG_All_unfiltered, ]$padj,
    
    row.names = DEG_All_unfiltered
  )
  
  den_sample_list <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(!!sym(group_str) == denominator_level_str) %>% pull(sample)
  num_sample_list <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(!!sym(group_str) == numerator_level_str) %>% pull(sample)
  
  DEG_All.expFCpadj <- merge(
    DEG_All.expFCpadj %>% rownames_to_column(var = 'gene_id'),
    read_counts_df %>% rownames_to_column(var = 'gene_id') %>% select(all_of(c('gene_id', den_sample_list, num_sample_list))),
    by = 'gene_id'
  ) %>% column_to_rownames(var = 'gene_id')
  
  
  dir.create(dir_xxx, recursive = TRUE, showWarnings = FALSE)
  
  png(file.path(dir_xxx, paste0(denominator_level_str, 'Vs', numerator_level_str, '_edgeR-DESeq2-limmaVoom_venn.png')))
  gplots::venn(DE_list)
  dev.off()
  
  # Get the values that appear in two or more lists
  rl <- rle(sort(unlist(DE_list)))
  DEG_atLeast_i_samples <- rl$values[rl$lengths >= 2]
  
  for(gene_list_and_name in list(
    list(l = DEG_atLeast_i_samples, name = 'AtLeast2', col_value = 'DESeq2_log2FC'),
    
    list(l = rownames(subset(limmaV.res.sorted, adj.P.Val <= pv_adj_threshold)), name = 'limmaVoom', col_value = 'LimmaV_log2FC', dir = '/single_tool'),
    list(l = rownames(subset(edgeR.ql.res.sorted$table, FDR <= pv_adj_threshold)), name = 'edgeR', col_value = 'edgeR_log2FC', dir = '/single_tool'),
    list(l = rownames(subset(DESeq.res_lfcShrink, padj <= pv_adj_threshold)), name = 'DESeq2', col_value = 'DESeq2_log2FC', dir = '/single_tool'),
    
    list(l = rownames(limmaV.res.sorted), name = 'limmaVoom_unfiltered', col_value = 'LimmaV_log2FC', dir = '/unfiltered'),
    list(l = rownames(edgeR.ql.res.sorted$table), name = 'edgeR_unfiltered', col_value = 'edgeR_log2FC', dir = '/unfiltered'),
    list(l = rownames(DESeq.res_lfcShrink), name = 'DESeq2_unfiltered', col_value = 'DESeq2_log2FC', dir = '/unfiltered')
  )
  ){
    gene_list <- gene_list_and_name$l
    print(length(gene_list))
    name_tools <- gene_list_and_name$name
    col_value <- gene_list_and_name$col_value
    
    x_df <- DEG_All.expFCpadj[gene_list,]
    gene_list.under <- rownames(x_df[x_df[col_value] < 0, ])
    gene_list.over <- rownames(x_df[x_df[col_value] > 0, ])
    
    print(paste0(name_tools, ': ', length(gene_list), ' = ', length(gene_list.under), ' (under) + ', length(gene_list.over), ' (over)'))
    
    dir_xxx_yyy <- paste0(dir_xxx, gene_list_and_name$dir)
    print(dir_xxx_yyy)
    dir.create(dir_xxx_yyy, recursive = TRUE, showWarnings = FALSE)
    if (name_tools == 'AtLeast2'){
      name_tools <- ''
    }else{
      name_tools <- paste0(name_tools, '.')
    }
    
    if (!is.null(ensembl_annotation_df)){
      merge(
        ensembl_annotation_df[gene_list, ] %>% rownames_to_column(var = 'gene_id'),
        DEG_All.expFCpadj[gene_list, ] %>% rownames_to_column(var = 'gene_id'),
        by = 'gene_id'
      ) %>% arrange(desc(!!sym(col_value))) %>%
        write_tsv(
          file.path(
            dir_xxx_yyy,
            paste0(denominator_level_str, 'Vs', numerator_level_str, '_genes.', name_tools, 'tsv')
          ),
          quote = FALSE
        )
      
      merge(
        ensembl_annotation_df[gene_list.under, ] %>% rownames_to_column(var = 'gene_id'),
        DEG_All.expFCpadj[gene_list.under, ] %>% rownames_to_column(var = 'gene_id'),
        by = 'gene_id'
      ) %>% arrange(!!sym(col_value)) %>%
        write_tsv(
          file.path(
            dir_xxx_yyy,
            paste0(denominator_level_str, 'Vs', numerator_level_str, '_genes.under.', name_tools, 'tsv')
          ),
          quote = FALSE
        )
      
      merge(
        ensembl_annotation_df[gene_list.over, ] %>% rownames_to_column(var = 'gene_id'),
        DEG_All.expFCpadj[gene_list.over, ] %>% rownames_to_column(var = 'gene_id'),
        by = 'gene_id'
      ) %>% arrange(desc(!!sym(col_value))) %>%
        write_tsv(
          file.path(
            dir_xxx_yyy,
            paste0(denominator_level_str, 'Vs', numerator_level_str, '_genes.over.', name_tools, 'tsv')
          ),
          quote = FALSE
        )
    }else{
      DEG_All.expFCpadj[gene_list, ] %>% rownames_to_column(var = 'gene_id') %>% arrange(desc(!!sym(col_value))) %>%
        write_tsv(
          file.path(
            dir_xxx_yyy,
            paste0(denominator_level_str, 'Vs', numerator_level_str, '_genes.', name_tools, 'tsv')
          ),
          quote = FALSE
        )
      
      DEG_All.expFCpadj[gene_list.under, ] %>% rownames_to_column(var = 'gene_id') %>% arrange(desc(!!sym(col_value))) %>%
        write_tsv(
          file.path(
            dir_xxx_yyy,
            paste0(denominator_level_str, 'Vs', numerator_level_str, '_genes.under.', name_tools, 'tsv')
          ),
          quote = FALSE
        )
      
      DEG_All.expFCpadj[gene_list.over, ] %>% rownames_to_column(var = 'gene_id') %>% arrange(!!sym(col_value)) %>%
        write_tsv(
          file.path(
            dir_xxx_yyy,
            paste0(denominator_level_str, 'Vs', numerator_level_str, '_genes.over.', name_tools, 'tsv')
          ),
          quote = FALSE
        )
    }
  }
}

group_str <- 'Condition'
numerator_level_str <- 'Infected'
denominator_level_str <- 'Control'
coef <- 'ConditionInfected'
pv_adj_threshold <- 0.05
design_formula <- ~ Condition

for (batch in unique(metadata_df$Batch)){
  print(batch)
  
  virus_list <- metadata_df %>% filter(Batch %in% batch) %>% filter(Virus != 'Mock') %>% pull(Virus) %>% unique()
  
  for(virus in virus_list){
    print(paste0('    ', virus))
    
    dir_output = file.path(opt$o, paste0(paste0(batch, collapse = '.'), '_', virus))
    
    sample_to_take <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(Batch %in% batch & Virus %in% c('Mock', virus)) %>% pull(sample)
    
    metadata_sub_df <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(sample %in% sample_to_take) %>% column_to_rownames(var = 'sample')
    metadata_sub_df$Batch <- factor(metadata_sub_df$Batch)
    
    diff_expr_analysis(
      read_counts_df = read_counts_df %>% select(all_of(sample_to_take)),
      metadata_df = metadata_sub_df,
      
      ensembl_annotation_df = ensembl_annotation_df,
      group_str = group_str,
      numerator_level_str = numerator_level_str,
      denominator_level_str = denominator_level_str,
      coef = coef,
      
      design_formula = design_formula,
      
      pv_adj_threshold = pv_adj_threshold,
      
      dir_output = dir_output
    )
    
    DESeq.ds <- DESeqDataSetFromMatrix(
      countData = read_counts_df %>% select(all_of(sample_to_take)),
      colData = metadata_sub_df,
      design = design_formula
    )
    to_keep_list  <- rowSums(counts(DESeq.ds) >= 5) >= (min(table(metadata_sub_df$Condition))/2)
    sum(to_keep_list)
    DESeq.ds <- DESeq.ds[to_keep_list, ]
    
    DESeq.vst <- vst(DESeq.ds, blind = TRUE)
    vst.norm.counts <- assay(DESeq.vst)

    pca <- prcomp(t(vst.norm.counts))
    percent <- round(summary(pca)$importance[2,]*100)
    
    pca_and_meta_df <- merge(
      metadata_sub_df %>% rownames_to_column(var = 'sample') %>% dplyr::select(sample, Line, Virus, Treatment, Condition, BiolRep, Batch),
      pca$x %>% as.data.frame() %>% select(PC1, PC2) %>% rownames_to_column(var = 'sample'),
      by = 'sample'
    ) %>% column_to_rownames('sample')
    
    if(length(unique(metadata_sub_df$Line)) > 1){
      autoplot(pca,
               data = pca_and_meta_df, 
               fill="Condition", 
               shape="Line",
               size=4) + 
        ggtitle(paste0("vst transformed counts")) + theme(plot.title = element_text(hjust = 0.5)) +
        scale_shape_manual(values=c(21, 22, 24))+ scale_fill_brewer(palette="Dark2") + 
        guides(fill = guide_legend(override.aes=list(shape=22)))
      ggsave(filename = file.path(dir_output, 'PCA_ConditionLine.png'), width = 6, height = 4)
    }else{
      autoplot(pca,
               data = pca_and_meta_df, 
               fill="Condition",
               shape="Treatment",
               size=4) + 
        ggtitle(paste0("vst transformed counts")) + theme(plot.title = element_text(hjust = 0.5)) +
        scale_shape_manual(values=c(21, 22, 24))+ scale_fill_brewer(palette="Dark2") + 
        guides(fill = guide_legend(override.aes=list(shape=22)))
      ggsave(filename = file.path(dir_output, 'PCA_ConditionTreatment.png'), width = 6, height = 4)
    }
    
    if(length(unique(metadata_sub_df$Batch)) > 1){
      autoplot(pca,
               data = pca_and_meta_df, 
               fill="Condition", 
               shape="Batch",
               size=4) + 
        ggtitle(paste0("vst transformed counts")) + theme(plot.title = element_text(hjust = 0.5)) +
        scale_shape_manual(values=c(21, 22, 24))+ scale_fill_brewer(palette="Dark2") + 
        guides(fill = guide_legend(override.aes=list(shape=22)))
      ggsave(filename = file.path(dir_output, 'PCA_ConditionBatch.png'), width = 6, height = 4)
    }
  }

}
