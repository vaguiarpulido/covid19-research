library(data.table)
library(seqinr)
library(Biostrings)
library(TFBSTools)
library(foreach)
library(doParallel)
library(rtracklayer)

###################################################################################

# Data I/O functions
###################################################################################

# Function to read GFFs
readGFF=function(gff_file, skipfirst=T){
  gff = import(gff_file)
  if(skipfirst){
    gff = gff[2:length(gff)]
  }
  gff = as.data.table(gff)
  gff = gff[, .(seqnames, start, end, type, gene)]
  return(gff)
}

# Function to read fasta > string
readFasta=function(fasta_file, toRNA=F){
  fa = read.fasta(fasta_file, as.string = T, forceDNAtolower = F)
  seq = fa[[1]][[1]]
  if(toRNA){
    seqString = RNAString(gsub("T", "U", seq))
  } else{
    seqString = DNAString(seq)
  }
  seqString = list(seqString)
  names(seqString) = names(fa)
  return(seqString)
}


# Function to read PWMs
readPWMsFromFasta = function (pwm_file) {

  # Read all lines from PWM file
  lines = readLines(pwm_file)

  # Find header lines, start and end of each PWM
  ind = which(substr(lines, 1L, 1L) == ">")
  nseq = length(ind)
  start = ind + 1
  end = ind - 1
  end = c(end[-1], length(lines))

  # Get PWM IDs
  ids = lapply(seq_len(nseq), function(i) {
    firstword <- strsplit(lines[ind[i]], "\t")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })

  # Split PWMs
  pwms = lapply(seq_len(nseq), function(i) strsplit(lines[start[i]:end[i]], "\t"))

  # Format as numeric matrix
  pwms = lapply(pwms, function(x) matrix(as.numeric(unlist(x)), ncol=4, byrow=T))

  # Convert to PWMatrix class
  pwms = lapply(seq_len(nseq), function(i){
    PWMatrix(profileMatrix=matrix(c(pwms[[i]]), byrow=TRUE, nrow=4, dimnames=list(c("A", "C", "G", "T"))), ID=ids[[i]])
  })

  # Name with PWM ID
  names(pwms) = ids
  return(pwms)
}


###################################################################################

# Function to calculate PWM information content
GetPWMEntropy = function(pwm){
  m=as.matrix(pwm)
  ic = -sum(m*log2(m))
  return(ic)
}

addGapsGFF=function(gff){
  gaps = as.data.table(gaps(IRanges(gff$start, gff$end)))
  gaps[, seqnames:=gff$seqnames[1]]
  gaps[, type:="intergenic"][, gene:=NA]
  gapped_gff = rbind(gff, gaps[, .(seqnames, start, end, type, gene)])
  gapped_gff = gapped_gff[order(start),]
  return(gapped_gff)
}

# Function to split genome sequence into regions
split_genome_gff=function(ref, gff){
  gff_seqs=list()
  for(i in 1:nrow(gff)){
    gff_seqs[[i]] = ref[[1]][gff$start[i]:gff$end[i]]
  }
  gff_seqs = DNAStringSet(gff_seqs)
  gff_seqs_by_type = list()
  for(typ in gff$type){
    gff_seqs_by_type[[typ]] = gff_seqs[gff$type==typ]
  }
  gff_seqs_by_type[["genome"]] = DNAStringSet(ref)
  #gff_seqs_by_type[["neg_genome"]] = reverseComplement(DNAStringSet(ref))
  return(gff_seqs_by_type)
}

###################################################################################

# Binding site discovery functions
###################################################################################

# Multithreaded function to scan sequence(s) with multiple PWMs
ScanSeqWithPWMs = function(seqString, pwmList, rbp_to_pwm, seqName, strand="*", print_num=TRUE){

  sites = foreach(i = 1:length(pwmList), .combine=rbind) %dopar% {
    # Read PWM ID
    id = as.character(names(pwmList)[i])
    # Scan genome
    curr_sites = searchSeq(pwmList[[i]], seqString, min.score="90%", strand=strand, seqname = seqName)
    if(length(curr_sites) > 0){
      # Convert to data table
      curr_sites = as.data.table(writeGFF3(curr_sites))
      if(length(curr_sites) > 0){
        # format
        curr_sites[, seq:= curr_sites[, tstrsplit(attributes, split=";|=", perl=T)][, V6]]
        curr_sites[, attributes:=NULL]
        curr_sites[, Matrix_id:= id]
      }
    }
  }
  # Match binding sites with protein name
  sites = merge(sites, rbp_to_pwm, by="Matrix_id")
  
  # Find duplicate binding sites for the same protein and choose the one with the highest score
  sites = sites[sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name, seqname)]$V1]
  sites = sites[!duplicated(sites[, .(start, end, strand, Gene_name, seqname)]),]

  # sort by position
  sites = sites[order(seqname, start),]
  
  # Filter columns
  sites[, source:=NULL]
  sites[, feature:=NULL]
  
  # Add site length
  sites[, len:=nchar(seq)]
  
  if(print_num){
    print(paste0("Found ", nrow(sites), " sites."))
  }
  
  return(sites)
}


# Annotate sites with genomic features
annotateSites=function(sites, gff){

  sites_ranges = GRanges(sites)
  gff_ranges = GRanges(gff)
  site_overlaps = findOverlaps(sites_ranges, gff_ranges, ignore.strand=F, type="any")
  
  annotated_sites = sites_ranges[queryHits(site_overlaps)]
  annotated_sites$type = gff_ranges$type[subjectHits(site_overlaps)]
  annotated_sites$gene = gff_ranges$gene[subjectHits(site_overlaps)]
  
  # Convert to data table
  annotated_sites = as.data.table(annotated_sites)
  return(annotated_sites)
}
###################################################################################

# Scrambling functions
###################################################################################

# Function to scramble a single sequence
scrambleSingleSeq = function(sequence, N){

  # Calculate nucleotide frequency in real sequence
  freqs = alphabetFrequency(sequence)[1:4]
  
  # List all bases to scramble
  bases=c()
  for(base_type in names(freqs)){
    bases = c(bases, rep(base_type, freqs[base_type]))
  }
  
  # Generate scrambled sequences
  sim_seqs = list()
  for(i in 1:N){
    sim_seqs[[i]] = DNAString(paste(sample(bases), collapse=""))
  }
  sim_seqs = DNAStringSet(sim_seqs)
  return(sim_seqs)
}

# Function to scramble a set of sequences
scrambleMultiSeq = function(seqs, N){
  
  # Calculate nucleotide frequency in real sequences
  widths = width(seqs)
  combined_seqs = DNAString(paste0(seqs, collapse=""))
  freqs = alphabetFrequency(combined_seqs)[1:4]
  
  # List all bases to scramble
  bases=c()
  for(base_type in names(freqs)){
    bases = c(bases, rep(base_type, freqs[base_type]))
  }
  
  # Get positions at which to split combined sequence
  positions = c(1)
  for(w in widths){
    positions = c(positions, positions[length(positions)] + w)
  }
  
  # Generate scrambled sequences
  sim_seqs = list()
  for(i in 1:N){
    sim_seqs[[i]] = paste(sample(bases), collapse="")
  }
  
  # Split scrambled sequences
  for(i in 1:N){
    stringSet = list()
    for(spos in 2:length(positions)){
      stringSet[[spos-1]] = DNAString(substr(sim_seqs[[i]], positions[spos-1], positions[spos]-1))
    }
    sim_seqs[[i]] = DNAStringSet(stringSet)
    names(sim_seqs[[i]]) = rep(i, length(sim_seqs[[i]]))
  }
  
  # Combine all sequences
  sim_seqs = do.call(c, sim_seqs)
  
  return(sim_seqs)
}



###################################################################################
# Enrichment functions
###################################################################################
# Enrichment test for RBP
enrich_rbps_real=function(real_sites, sim_sites){
  
  print("Counting real binding sites per protein per strand")
  site_count = real_sites[, .N, by=.(Gene_name, strand)]
  
  print("Counting binding sites per protein per strand, on the simulated sequence")
  sim_site_count = sim_sites[, .N, by=.(seqname, Gene_name, strand)]
  sim_site_count = dcast(sim_site_count, seqname+strand~Gene_name, value.var = "N", fill=0)
  sim_site_count = melt(sim_site_count, id.vars=1:2, variable.name = "Gene_name", value.name = "N")
  sim_site_count = sim_site_count[, .(mean_count=mean(N), sd_count=sd(N)), by=.(Gene_name, strand)]
  
  print("Comparing binding site counts")
  site_count = merge(site_count, sim_site_count, by = c("Gene_name", "strand"), all=T)
  site_count[is.na(N), N:=0]
  
  # Calculate z-score
  print("Calculating z-scores")
  site_count[, z:=(N-mean_count)/sd_count]
  site_count[, pval1:=pnorm(-z)]
  site_count[, pval2:=2*pnorm(-abs(z))]
  
  # Multiple hypothesis correction
  print("FDR correction")
  site_count[, padj1:=p.adjust(pval1, "fdr"), by=strand]
  site_count[, padj2:=p.adjust(pval2, "fdr"), by=strand]
  
  return(site_count)
}


######################################################################

# Apply functions to each annotated region of a genome

SimulateSeqsByRegion=function(seqsets, N, regions=c("genome", "three_prime_UTR", "five_prime_UTR", "intergenic")){
  result=list()
  for(region in regions){
    print(paste0("Simulating ", N, " copies of ", region))
    if(region=="intergenic"){
      result[[region]] = scrambleMultiSeq(seqsets[[region]], N)
    }
    else{
      result[[region]] = scrambleSingleSeq(seqsets[[region]][[1]], N)
    }
  }
  return(result)  
}


EnrichSeqSetsByRegion=function(seqsets, sim_seqsets, pwmList, rbp_to_pwm, prefix, real_sites, regions=c("three_prime_UTR", "five_prime_UTR", "intergenic", "genome")){
  
  for(region in regions){
    sim = sim_seqsets[[region]]
    if(region!="genome"){
      print(paste0("Finding binding sites in simulated ", region))
      sim_sites = ScanSeqWithPWMs(sim, pwmList, rbp_to_pwm, strand="+")
      print("saving")
      save(sim_sites, file=paste0(prefix, "_sim_sites_", region, ".RData"))
      print("Enrichment")
      enr=enrich_rbps_real(real_sites[type==region][strand=="+"], sim_sites)
      print("saving")
      save(enr, file=paste0(prefix, "_sim_enrich_", region, ".RData"))
      rm(sim_sites)
      rm(enr)
    }
    else{
      # Positive strand
      print("Finding binding sites")
      sim_sites = ScanSeqWithPWMs(sim, pwmList, rbp_to_pwm, strand="+")
      print("saving")
      save(sim_sites, file=paste0(prefix, "_sim_sites_genome.RData"))
      print("Enrichment")
      enr=enrich_rbps_real(real_sites[strand=="+"], sim_sites)
      save(enr, file=paste0(prefix, "_sim_enrich_genome.RData"))
      rm(sim_sites)
      rm(enr)
      # Negative strand
      print("Finding binding sites - negative")
      sim_sites = ScanSeqWithPWMs(sim, pwmList, rbp_to_pwm, strand="-")
      print("saving")
      save(sim_sites, file=paste0(prefix, "_sim_sites_neg_genome.RData"))
      print("Enrichment - negative")
      enr = enrich_rbps_real(real_sites[strand=="-"], sim_sites)
      print("saving - negative")
      save(enr, file=paste0(prefix, "_sim_enrich_neg_genome.RData"))
      rm(sim_sites)
      rm(enr)
    }
  }
}


ConstructEnrichmentList=function(prefix, regions=c("genome", "neg_genome", "three_prime_UTR", "five_prime_UTR", "intergenic")){

  result=list()

  for(region in regions){
    load(paste0(prefix, "_sim_enrich_", region, ".RData"))
    result[[region]] = enr
  }

  return(result)
}


FilterToEnriched=function(sites, enrichment, minq, regions=c("genome", "three_prime_UTR", "five_prime_UTR", "intergenic", "neg_genome")){
  
  results=list()
  
  for(region in regions){
    if(region=="genome"){
      results[[region]] = sites[strand=="+"][Gene_name %in% enrichment[[region]][padj1<minq]$Gene_name]
    }
    else if(region=="neg_genome"){
      results[[region]] = sites[strand=="-"][Gene_name %in% enrichment[[region]][padj1<minq]$Gene_name]
    }
    else {
      results[[region]] = sites[strand=="+"][type==region][Gene_name %in% enrichment[[region]][padj1<minq]$Gene_name]
    }
  }
  
  return(results)
}


SiteConservation=function(candidate_sites, mat, ref_base_count){
for(i in 1:nrow(candidate_sites)){
  if(i %% 10==0){print(i)}
  start_col=ref_base_count[base==candidate_sites[i, start], col]
  end_col=ref_base_count[base==candidate_sites[i, end], col]
  if(length(start_col)>1){start_col = start_col[1]}
  if(length(end_col)>1){end_col = end_col[1]}
  site_mat=mat[, start_col:end_col]
  site_mat=as.vector(apply(site_mat, 1, function(x){paste0(x, collapse="")}))
  ref_seq=site_mat[[1]]
  chars = sapply(site_mat, function(x){unique(strsplit(x, "")[[1]])})
  candidate_sites$n_N[i]=length(grep("N|Y|K|R|W|V", chars))
  candidate_sites$n_match[i]=sum(site_mat==ref_seq)
}
candidate_sites[,N_cand:=nrow(mat)-n_N]
return(candidate_sites)
}