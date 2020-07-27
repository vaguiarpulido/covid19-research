######################################################################################

# Import requirements
library(data.table)
library(TFBSTools)
library(seqinr)
library(Biostrings)
library(rtracklayer)
library(foreach)
library(doParallel)
library(readxl)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

registerDoParallel(detectCores())

######################################################################################

# Parameters

n_sims=1000

######################################################################################
# List input data files

## ATtRACT # Downloaded from https://attract.cnic.es/download
rbp_file = "ATtRACT/ATtRACT_db.txt"
pwm_file = "ATtRACT/pwm.txt"

# SARS-CoV-2 genome files # Downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.
ref_file="reference/sequence.fasta.txt"
gff_file = "reference/GCF_009858895.2_ASM985889v3_genomic.gff"

# SARS1 reference genome
ref_file_sars1 = "reference/SARS1/sequence.fasta.txt"
gff_file_sars1 = "reference/SARS1/sequence.gff3"

# RaTG13 reference genome
ref_file_ratg13 = "reference/RaTG13/MN996532.1_sequence.fasta.txt"
gff_file_ratg13 = "reference/RaTG13/sequence.gff3"


######################################################################################
# # Read input data

print("Reading RBPs")
rbp = fread(rbp_file, select=c(1:4, 7:9, 11:13))
rbp = unique(rbp)

print("Reading PWMs")
pwm = readPWMsFromFasta(pwm_file)

print("Reading reference genome")
ref = readFasta(ref_file)

print("Reading GFF")
gff = readGFF(gff_file)
gff = gff[type!="CDS", ]
gff = addGapsGFF(gff)
  
# Reading SARS1 reference genome
ref_sars1 = readFasta(ref_file_sars1)
gff_sars1 = readGFF(gff_file_sars1)
gff_sars1 = gff_sars1[type!="CDS", ]

# Reading RaTG13 reference genome
ref_ratg13 = readFasta(ref_file_ratg13)
gff_ratg13 = readGFF(gff_file_ratg13)
gff_ratg13 = gff_ratg13[type!="CDS", ]

######################################################################################

# 1. Filter RBPs and PWMs

# Count number of RBPs in database
initial_n = rbp[, length(unique(Gene_name))]

# Select only human RBPs
print("Filtering human RBPs")
rbp = rbp[Organism == "Homo_sapiens",]
print(paste0("Reduced number of RBPs from ", initial_n, " to ", rbp[, length(unique(Gene_name))]))

# Remove PWMs which consist of only a single motif
initial_n = rbp[, length(unique(Gene_name))]
rbp = rbp[Score!="1.000000**",]
rbp = unique(rbp[, Score:=NULL])
print(paste0("Reduced number of RBPs from ", initial_n, " to ", rbp[, length(unique(Gene_name))]))

# # Select PWMs that match these RBPs
print("Selecting PWMs that match to human RBPs")
initial_n = length(pwm)
pwm = pwm[names(pwm) %in% rbp[, Matrix_id]]
print(paste0("Reduced number of PWMs from ", initial_n, " to ", length(pwm)))

# # Filter out high-entropy PWMs
print("Removing high-entropy PWMs")
entropy = sapply(pwm, GetPWMEntropy)
#plot(density(entropy))
initial_n = length(pwm)
pwm = pwm[entropy < 8]
print(paste0("Reduced number of PWMs from ", initial_n, " to ", length(pwm)))

# # Filter the RBPs corresponding to these PWMs
print("Filtering RBP-PWM matches after removing high-entropy PWMs")
initial_n = rbp[, length(unique(Gene_name))]
rbp = rbp[Matrix_id %in% names(pwm),]
print(paste0("Reduced number of RBPs in the search from ", initial_n, " to ", rbp[, length(unique(Gene_name))]))

# # Match PWMs to RBPs
rbp_to_pwm = unique(rbp[, .(Gene_name, Matrix_id)])

# #####################################################################################
# # Save filtered data
# print("Saving")
# save(rbp, file="output/filtered_rbp.RData")
# save(pwm, file="output/filtered_pwm.RData")
# save(rbp_to_pwm, file="output/rbp_to_pwm.RData")

load("output/filtered_rbp.RData")
load("output/filtered_pwm.RData")
load("output/rbp_to_pwm.RData")
#####################################################################################

# 2. Find binding sites in the SARS-CoV-2 reference genome and other genomes

# SARS-CoV-2

# # Scan the genome with these PWMs
print("Scanning the reference genome for PWMs")
sites = ScanSeqWithPWMs(ref[[1]], pwm, rbp_to_pwm, seqName=names(ref))
sites = annotateSites(sites, gff)

# Save sites
save(sites, file="output/sites.RData")

#SARS1
print("Scanning the SARS1 reference genome for PWMs")
sites_sars1 = ScanSeqWithPWMs(ref_sars1[[1]], pwm, rbp_to_pwm, seqName=names(ref_sars1))
sites_sars1 = annotateSites(sites_sars1, gff_sars1)
save(sites_sars1, file="output/sites_sars1.RData")

#RaTG13
gff_ratg13=rbind(gff_ratg13, data.table(seqnames="MN996532.1", start=1, end=250, type="five_prime_UTR", gene=NA))
gff_ratg13=rbind(gff_ratg13, data.table(seqnames="MN996532.1", start=29500, end=29855, type="three_prime_UTR", gene=NA))

print("Scanning the RaTG13 reference genome for PWMs")
sites_ratg13 = ScanSeqWithPWMs(ref_ratg13[[1]], pwm, rbp_to_pwm, seqName=names(ref_ratg13))
sites_ratg13 = annotateSites(sites_ratg13, gff_ratg13)
save(sites_ratg13, file="output/sites_ratg13.RData")


############################################################

#2. Simulate genomes

# SARS2
seqsets = split_genome_gff(ref, gff)
save(seqsets, file="output/seqsets.RData")

sim_seqsets = SimulateSeqsByRegion(seqsets, n_sims)
save(sim_seqsets, file="output/sim_seqsets.RData")

#SARS1
seqsets_sars1 = split_genome_gff(ref_sars1, gff_sars1)
save(seqsets_sars1, file="output/seqsets_sars1.RData")

sim_seqsets_sars1 = SimulateSeqsByRegion(seqsets_sars1, n_sims, regions=c("three_prime_UTR", "five_prime_UTR"))
save(sim_seqsets_sars1, file="output/sim_seqsets_sars1.RData")

#RaTG13
seqsets_ratg13 = split_genome_gff(ref_ratg13, gff_ratg13)
save(seqsets_ratg13, file="output/seqsets_ratg13.RData")
sim_seqsets_ratg13 = SimulateSeqsByRegion(seqsets_ratg13, n_sims, regions=c("three_prime_UTR", "five_prime_UTR"))
save(sim_seqsets_ratg13, file="output/sim_seqsets_ratg13.RData")

###################################################

# Enrichment

# SARS2
load("output/sites.RData")
load("output/seqsets.RData")
load("output/sim_seqsets.RData")

# scan each seqSet
EnrichSeqSetsByRegion(seqsets, sim_seqsets, pwm, rbp_to_pwm, "output/sars2", sites)

# Combine
enrichment=ConstructEnrichmentList("output/sars2")
save(enrichment, file="output/enrichment.RData")

# Filter sites of enriched proteins by region
sites_by_enriched = FilterToEnriched(sites, enrichment, 0.01)
save(sites_by_enriched, file="output/sites_by_enriched.RData")

rm(sites)
rm(seqsets)
rm(sites_by_enriched)
rm(enrichment)

# SARS1

load("output/sites_sars1.RData")
load("output/seqsets_sars1.RData")
load("output/sim_seqsets_sars1.RData")

EnrichSeqSetsByRegion(seqsets_sars1, sim_seqsets_sars1, pwm, rbp_to_pwm, "output/sars1", sites_sars1, regions=c("three_prime_UTR", "five_prime_UTR", "genome", "neg_genome"))
enrichment_sars1=ConstructEnrichmentList("output/sars1", regions=c("three_prime_UTR", "five_prime_UTR", "genome", "neg_genome"))
save(enrichment_sars1, file="output/enrichment_sars1.RData")

sites_by_enriched_sars1 = FilterToEnriched(sites_sars1, enrichment_sars1, 0.01, regions=c("three_prime_UTR", "five_prime_UTR", "genome", "neg_genome"))
save(sites_by_enriched_sars1, file="output/sars1_sites_by_enriched.RData")

# RaTG13

load("output/sites_ratg13.RData")
load("output/seqsets_ratg13.RData")
load("output/sim_seqsets_ratg13.RData")

EnrichSeqSetsByRegion(seqsets_ratg13, sim_seqsets_ratg13, pwm, rbp_to_pwm, "output/ratg13", sites_ratg13, regions=c("three_prime_UTR", "five_prime_UTR", "genome", "neg_genome"))
enrichment_ratg13=ConstructEnrichmentList("output/ratg13", regions=c("three_prime_UTR", "five_prime_UTR", "genome", "neg_genome"))
save(enrichment_ratg13, file="output/enrichment_ratg13.RData")

sites_by_enriched_ratg13 = FilterToEnriched(sites_ratg13, enrichment_ratg13, 0.01, regions=c("three_prime_UTR", "five_prime_UTR", "genome", "neg_genome"))
save(sites_by_enriched_ratg13, file="output/ratg13_sites_by_enriched.RData")

######################################################################

# # Load full gisaid MSA
aln=read.alignment("reference/GISAID_data/msa_0520.fasta", format="fasta", forceToLower = FALSE)

# # Move reference to row 1
grep("Wuhan-Hu-1", aln$nam)
sq1 = aln$seq[[18361]]
nm1 = aln$nam[[18361]]
aln$nam[[18361]] = aln$nam[[1]]
aln$seq[[18361]] = aln$seq[[1]]
aln$nam[[1]] = nm1
aln$seq[[1]] = sq1

# # Save
save(aln, file="gisaid_msa_0520.RData")

# # Convert to matrix
mat = matrix("-", nrow = aln[[1]], ncol = nchar(aln$seq[[1]]))
rownames(mat) = aln$nam
for(i in 1:nrow(mat)){
  mat[i, ] = strsplit(toupper(aln$seq[[i]]), "")[[1]]
}
save(mat, file="gisaid_msa_0520.mat.RData")

# # Count bases of reference in MSA
ref_base_count=data.table(col = 1:ncol(mat), base=0)
base_reached = 0
for(i in 1:ncol(mat)){
  aln_col = mat[1, i]
  if(aln_col != "-"){
    base_reached = base_reached + 1
  }
  ref_base_count[i, "base"] = base_reached
}
save(ref_base_count, file="gisaid_msa_0520.ref_base_count.RData")

##########################################################

load("output/sites_by_enriched.RData")

# # Calculate binding site conservation number
enriched_site_conservation=list()
for(region in c("five_prime_UTR", "three_prime_UTR", "intergenic")){
	print(paste0("Calculating conservation for enriched sites in ", region))
	enriched_site_conservation[[region]] = SiteConservation(sites_by_enriched[[region]], mat, ref_base_count)
}

# Save
save(enriched_site_conservation, file="enriched_site_conservation.RData")

# ######################################################################

# Check for RBPs in PPI
t1 = as.data.table(read_xlsx("gordon_supplement/media-6.xlsx", skip = 1))
t1_rbps = t1[PreyGene %in% rbp$Gene_name]

# Check GTex
gtex=fread("gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")

# Check single-cell data
lungsc = fread("../scrna_demo/data/GSE122960/Donor01-08_count_matrix.csv")
dp=intersect(which(lungsc[gene=="ACE2"]>0), which(lungsc[gene=="TMPRSS2"]>0))
lungsc_dp=lungsc[, ..dp]
save(lungsc_dp, file="output/lungsc_dp.RData")
rbps_dp=melt(lungsc_dp[gene %in% rbp$Gene_name], id.vars=1)[value>0, unique(gene)]
intersect(rbps_dp, as.character(enrichment$genome[padj1<0.01, Gene_name]))
intersect(rbps_dp, as.character(enrichment$neg_genome[padj1<0.01, Gene_name]))
intersect(rbps_dp, as.character(enrichment$five_prime_UTR[padj1<0.01, Gene_name]))
intersect(rbps_dp, as.character(enrichment$three_prime_UTR[padj1<0.01, Gene_name]))
intersect(rbps_dp, as.character(enrichment$intergenic[padj1<0.01, Gene_name]))


##########################################################
