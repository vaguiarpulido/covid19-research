run_ebseq<-function(counts,conditions,outfile,threshold_counts,maxround){
  
  require(EBSeq)
  if(missing(threshold_counts)){
    threshold_counts=10}
  if(missing(maxround)){
    maxround=50}
  
  low_counts<-count_table[rowSums(count_table)<threshold_counts,]
  counts<-count_table[rowSums(count_table)>=threshold_counts,]  
  
  sizes <- MedianNorm(counts)
  GeneNormData <- GetNormalizedMat(counts, sizes)
  EBOut=EBTest(Data=GeneNormData,Conditions=conditions,sizeFactors=sizes, maxround=maxround)
  GeneFC=PostFC(EBOut)

  df<-as.data.frame(cbind(GeneNormData,EBOut$PPMat,GeneFC$PostFC))
  write.csv(df,file=paste("EBSeq_",outfile,sep=""))
  write.csv(low_counts,file=paste("low_counts_",outfile,sep=""))

}