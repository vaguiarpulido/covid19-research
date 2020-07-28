#### Metabolism COVID-19
rm(list = ls())
library(data.table)
countData<-read.csv("gene_count_matrix.csv",row.names = "gene_id")

#### Run EBSeq2 by Series
source("run_ebseq.R")

### Series1
count_table<-countData[,1:6]
DE=as.factor(c("Mock","Mock","Mock","SARS-CoV-2","SARS-CoV-2","SARS-CoV-2"))
outfile=paste("Series1.csv",sep="")
run_ebseq(count_table,DE,outfile)

### Series2
count_table<-countData[,7:12]
DE=as.factor(c("Mock","Mock","Mock","SARS-CoV-2","SARS-CoV-2","SARS-CoV-2"))
outfile=paste("Series2.csv",sep="")
run_ebseq(count_table,DE,outfile)

### Series3
count_table<-countData[,13:16]
DE=as.factor(c("Mock","Mock","RSV","RSV"))
outfile=paste("Series3.csv",sep="")
run_ebseq(count_table,DE,outfile)

### Series4
count_table<-countData[,c(19,20,17,18)]
DE=as.factor(c("Mock","Mock","IAV","IAV"))
outfile=paste("Series4.csv",sep="")
run_ebseq(count_table,DE,outfile)

### Series5
count_table<-countData[,21:26]
DE=as.factor(c("Mock","Mock","Mock","SARS-CoV-2","SARS-CoV-2","SARS-CoV-2"))
outfile=paste("Series5.csv",sep="")
run_ebseq(count_table,DE,outfile)

### Series7
count_table<-countData[,27:32]
DE=as.factor(c("Mock","Mock","Mock","SARS-CoV-2","SARS-CoV-2","SARS-CoV-2"))
outfile=paste("Series7.csv",sep="")
run_ebseq(count_table,DE,outfile)

### Series8 HPIV3
count_table<-countData[,c(36,37,38,33,34,35)]
DE=as.factor(c("Mock","Mock","Mock","HPIV3","HPIV3","HPIV3"))
outfile=paste("Series8H.csv",sep="")
run_ebseq(count_table,DE,outfile)

### Series8 RSV
count_table<-countData[,36:41]
DE=as.factor(c("Mock","Mock","Mock","RSV","RSV","RSV"))
outfile=paste("Series8R.csv",sep="")
run_ebseq(count_table,DE,outfile)

### Series9 IAV
count_table<-countData_entrez[,c(46:49,42:45)]
DE=as.factor(c("Mock","Mock","Mock","Mock","IAV","IAV","IAV","IAV"))
outfile=paste("Series9.csv",sep="")
run_ebseq(count_table,DE,outfile)

