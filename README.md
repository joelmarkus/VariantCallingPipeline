# Variant Calling Pipeline

**A wrapper tool to generate variant calling .vcf and .bed files from raw reads.**

The pipeline will have the following steps:

1. FASTQ reads alignment to a reference genome to create an alignment file.
2. Processing the alignment file (file format conversion, sorting, alignment improvement).
3. Calling the variants.

The following tools are used for the pipeline:

1. bwa for the alignment
2. samtools/HTS package for processing and calling variants
3. GATK v3.7.0 for improving the alignment.

Before running: 
- Please specify the entire path of your files.
- The jar location of GATK should be ~/bin/GenomeAnalysisTK.jar.

**Syntax:**

```
./ NWAlign.py -a <input reads 1> -b <input reads 2> -r <reference genome>
```

The code also uses a realignment option that uses Mills file for improved alignment. The base file can be downloaded from [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false&pli=1).

The resulting BED file from VCF will be in the following format:

```
Chromosome  Start Stop  Length
```
