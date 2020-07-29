# TE analysis

## Requirements
[TETools](https://github.com/l-modolo/TEtools)

## Input
Fastq files

## Workflow

Fastq files were used as input to the TEcount module from TEtools:

```
python3 TEcount.py -rosette XXX -column 3 -TE XXX -count outfile.txt -RNA reads.fq
```
