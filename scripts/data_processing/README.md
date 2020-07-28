# RNAseq data processing

## Software
- [STAR](https://github.com/alexdobin/STAR) (v2.3.7a)
- [Samtools](http://www.htslib.org/) (v1.10)
- [StringTie](http://ccb.jhu.edu/software/stringtie/index.shtml) (v2.1.1)
- [prepDE.py](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py) (provided by the creators of StringTie)

## Reference genome
- [GRCh37.p13 (GENCODE Release 19)](https://www.gencodegenes.org/human/release_19.html)

## Input
- Fastq.gz files

## Output
- Count matrices

# Steps
1. Map reads to the reference genome using `STAR`
2. Sort and index alignment file (BAM) using `Samtools`
3. Perform transcript quantification using `StringTie`
4. Generate count matrices to be employed by DESeq2 using `prepDE.py`
