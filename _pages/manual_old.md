---
layout: page
title: Manual Pre v1.1
hide: True
permalink: /manual_old/
order: 3
share: false
---

{: .no_toc}

- TOC
{:toc}

# Introduction

This is the old manual for versions earlier than 1.1. Please use this information for those versions only.

---

# Scripts

## hisatgenotype.py -- Analysis of a whole human genome
The hisatgenotype.py python script will analyze a whole human genome using whole genome sequencing reads. It will align reads to genotype genome, extract reads belonging to each locus of interest, and perform typing and assembly. 

Argument | Description | Default | Example
-------- | ----------- | ------- | -------
**\--base** | Base file name for index, variants, haplotypes, etc. | genotype_genome | `--base genotype_genome`
**\--region-list** | Comma separated list of regions | (empty) meaning every region available | `--region-list hla.A,hla.B,hla.C,codis.FGA,cyp.2D6`
**-p / \--threads** | Number of threads to be used  | 1 | `-p 4`
**-U** | Single-end read file name | *None* | `-U read.fq.gz`
**-1** | Paired-end read file name 1 | *None* | `-1 read.1.fq.gz`
**-2** | Paired-end read file name 2 | *None* | `-2 read.2.fq.gz`
**-f / \--fasta** | Reads are provided in FASTA format | *None* | `-f`
**\--alignment-file** | Sorted bam alignment file converted from HISAT2's sam alignment output | *None* | `--alignment-file NA12878.bam`
**\--num-editdist** | Number of maximum edit distance or mismatches allowed in the alignment, typing, and assembly of a read (not a pair) | 2 | `--num-editdist 0`
**\--assembly** | Perform assembly of each locus of interest | disabled | `--assembly`
**\--verbose** | Provide more information | disabled | `--verbose`
**-h / \--help** | Output help message | disabled | `--help `

Examples:
```bash
$ hisatgenotype.py --base genotype_genome -p 4 -1 read.1.fq.gz -2 read.2.fq.gz
$ hisatgenotype.py --verbose --base genotype_genome -p 4 --region-list hla.A -1 read.1.fq.gz -2 read.2.fq.gz
$ hisatgenotype.py --base genotype_genome --alignment-file NA12878.bam  
$ hisatgenotype.py --base genotype_genome --region-list hla.A,hla.B,hla.C,hla.DQA1,hla.DQB1,hla.DRB1,codis,cyp.2D6 --alignment-file NA12878.bam
```

## hisatgenotype_locus.py -- Analysis of a specific gene, genomic region, or gene family 
The hisatgenotype_locus.py python script will analyze a particular locus or a set of genes or genomic regions using extracted reads. As the name of the script (locus) implies, it will align reads to a part of genotype genome and perform typing and assembly.

Argument | Description | Default | Example
-------- | ----------- | ------- | -------
**\--base** | Base file name for index, variants, haplotypes, etc. | hla | `--base codis`
**\--locus-list** | Comma separated list of loci | (empty) meaning every locus available | `--locus-list A,B,C,DQA1,DQB1,DRB1`
**-p / \--threads** | Number of threads to be used | 1 | `-p 4`
**-U** | Single-end read file name | *None* | `-U read.fq.gz`
**-1** | Paired-end read file name 1 | *None* | `-1 read.1.fq.gz`
**-2** | Paired-end read file name 2 | *None* | `-2 read.2.fq.gz`
**-f / \--fasta** | Reads are provided in FASTA format | false | `-f`
**\--num-editdist** | Number of maximum edit distance or mismatches allowed in the alignment, typing, and assembly of a read (not a pair) | 2 | `--num-editdist 0`
**\--assembly** | Perform assembly of each locus of interest | disabled | `--assembly`
**\--genotype-genome** | Base name for genotype genome, which the program will use instead of region-based small indexes | *None* | `--genotype-genome genotype_genome`
**\--verbose** | Provide more information | disabled | `--verbose`
**-h / \--help** | Output help message | disabled | `--help` 

Examples:
```bash
 $ hisatgenotype_locus.py --base hla --locus-list A --num-editdist 2 --assembly -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz
 $ hisatgenotype_locus.py --base hla --locus-list A,B,C,DQA1,DQB1,DRB1 --num-editdist 2 -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz
```

### For advanced users
In addition to analyzing real sequencing reads, the script comes with additional testing capabilities, which might come in handy when determining whether a problematic result is due to noise in sequencing reads or algorithms. To do so, the script can generate simulation reads with sequencing errors and SNPs, and test and evaluate the accuracy of its typing and assembly algorithms, with additional options described below.

Argument | Description | Default | Example
-------- | ----------- | ------- | -------
**\--read-len** | Length of a simulated read | 100 | `--read-len 150`
**\--fragment-len** | Length of a fragment | 350 | `--fragment-len 500`
**\--perbase-errorrate** | Sequencing error rate per base in percentage | 0.0 (%) | `--perbase-errorrate 0.2`
**\--perbase-snprate** | SNP rate per base in percentage | 0.0 (%) | `--perbase-snprate 0.1`
**\--skip-fragment-regions** | Comma separated list of regions from which no reads originate | *None* | `--skip-fragment-regions 500-600,1200-1400`
**\--random-seed** | Seeding number for random number generation | 1 | `--random-seed 10`
**\--debug** | Comma separated list of debugging parameters | *None* | `--debug "pair,test_list:A*03:01:01:05,A*24:02:01:01"`

Examples:
```bash
 $ hisatgenotype_locus.py --base hla --locus-list A 
 $ hisatgenotype_locus.py --base hla --locus-list A --debug "pair"
 $ hisatgenotype_locus.py --base hla --locus-list B --assembly --debug "pair,test_list:A*03:01:01:05,A*24:02:01:01"
 $ hisatgenotype_locus.py --base hla --locus-list A --assembly --simulate-interval 10 --perbase-errorrate 0.3 --perbase-snprate 0.1 --num-editdist 2 --debug "pair,full"
 $ hisatgenotype_locus.py --base hla --locus-list A --assembly --simulate-interval 10 --perbase-errorrate 0.3 --perbase-snprate 0.1 --num-editdist 2 --debug "pair,full" --skip-fragment-regions 1000-1020 --random-seed 0
```

## hisatgenotype_build_genome.py -- Building a genotype genome 
The hisatgenotype_build_genome.py python script will create a genotype genome using the human reference genome (GRCh38) and several genomic databases. A tutorial for the script can be found here.

Argument | Description | Default | Example
-------- | ----------- | ------- | -------
**\--base** | Base file name for index, variants, haplotypes, etc. | genotype_genome | `--base genotype_genome`
**\--database-list** | Comma separated list of databases | hla | `--database-list hla,codis,cyp`
**\--commonvar** | Include common variants from dbSNP (about 13 million) | false | `--commonvar`
**\--clinvar** | Include clinically relevant variants from NCBI ClinVar database | false | `--clinvar`
**-p / \--threads** | Number of threads to be used | 1 | `-p 4`
**\--verbose** | Provide more information | disabled | `--verbose`
**-h / \--help** | Output help message | disabled | `--help`

Examples:
```bash
 $ hisatgenotype_build_genome.py -p 4 --base genotype_genome --database-list hla --commonvar
 $ hisatgenotype_build_genome.py -p 4 --base genotype_genome --database-list hla,codis,cyp --commonvar
```

## hisatgenotype_extract_vars.py -- Extraction of variants, haplotypes, etc.
This hisatgenotype_extract_vars.py python script is typically run by other scripts, though instructions are provided here for users who wish to directly run the script. The script will extract allele sequences, variants, haplotypes, etc. from multiple sequence alignments available in hisatgenotype_db. Tutorials for the script can be found at HLA typing and DNA fingerprinting anlaysis.

Argument | Description | Default | Example
-------- | ----------- | ------- | -------
**\--base** | Base file name for index, variants, haplotypes, etc. | hla | `--base codis`
**\--locus-list** | Comma separated list of loci | (empty) meaning every locus available | `--locus-list A,B,C,DQA1,DQB1,DRB1`
**\--min-var-freq** | Exclude variants whose freq is below than this value in percentage | 0.0 (%) | `--min-var-freq 0.1`
**\--leftshift** | Shift deletions to the leftmost | disabled | `--leftmost`
**\--inter-gap** | Maximum distance for variants to be in the same haplotype | 30 | `--inter-gap 30`
**\--intra-gap** | Break a haplotype into several haplotypes | 50 | `--intra-gap 50`
**\--verbose** | Provide more information | disabled | `--verbose`
**-h / \--help** | Output help message | disabled | `--help`

Examples:
```bash
 $ hisatgenotype_extract_vars.py --base hla --locus-list A,B,C,DQA1,DQB1,DRB1 --min-var-freq 0.1
 $ hisatgenotype_extract_vars.py --base codis
```

## hisatgenotype_extract_reads.py -- Extraction of reads that belong to genomic regions of interest
The hisatgenotype_extract_reads.py python script will extract reads that belong to loci of interest from whole sequencing reads. The script can also be used to extract reads from many samples. A tutorial for the script can be found here.

Argument | Description | Default | Example
-------- | ----------- | ------- | -------
**\--base** | Base file name for index, variants, haplotypes, etc. | genotype_genome | `--base genotype_genome`
**\--database-list** | Comma separated list of databases | (empty) meaning every database available | `--database-list hla,codis,cyp`
**-U** | Single-end read file name | *None* | `-U read.fq.gz`
**-1** | Paired-end read file name 1 | *None* | `-1 read.1.fq.gz`
**-2** | Paired-end read file name 2 | *None* | `-2 read.2.fq.gz`
**\--read-dir** | Name of directory where read files are placed | *None* | `--read-dir test_input`
**\--out-dir** | Name of directory where extracted read files will be placed | *None* | `--out-dir test_output`
**\--suffix** | Read file suffix name | fq.gz | `--suffix fastq.gz`
**-f / \--fasta** | Reads are provided in FASTA format | false | `-f`
**\--num-editdist** | Number of maximum edit distance or mismatches allowed in the alignment, typing, and assembly of a read (not a pair) | 2 | `--num-editdist 0`
**-p / \--threads** | Number of threads to be used | 1 | `-p 4`
**\--max-sample** | Number of samples to be extracted | sys.maxint | `--max-sample 917`
**\--job-range** | Comma separated two numbers representing a range | 0,1 | `--job-range 1,4`
**\--verbose** | Provide more information | disabled | `--verbose`
**-h / \--help** | Output help message | disabled | `--help`

Examples:
```bash
 $ hisatgenotype_extract_reads.py --base genotype_genome --database-list hla,codis,cyp -1 test.1.fq.gz -2 test.2.fq.gz
 $ hisatgenotype_extract_reads.py --base genotype_genome --database-list hla,codis,cyp --read-dir test_input --out-dir test_output
```

## hisatgenotype_locus_samples.py -- Analysis of many samples 
The hisatgenotype_locus_samples.py python script will analyze genomic regions of interest using reads that are extracted from using hisatgenotype_extract_reads.py. Please avoid directly running the script on whole genomic sequencing reads as doing so will lead to inaccurate results due to reasons partly described here. A tutorial for the script can be found at HLA typing. The script is designed to analyze many samples.

Argument | Description | Default | Example
-------- | ----------- | ------- | -------
**\--base** | Base file name for index, variants, haplotypes, etc. | genotype_genome | `--base genotype_genome`
**\--region-list** | Comma separated list of regions | (empty) meaning every region available | `--region-list hla.A,hla.B,hla.C,codis.FGA,cyp.2D6`
**\--read-dir** | Name of directory where read files are placed | *None* | `--read-dir test_input`
**\--out-dir** | Output directory name | *None* | `--out-dir test_output`
**-f / \--fasta** | Reads are provided in FASTA format | false | `-f`
**\--num-editdist** | Number of maximum edit distance or mismatches allowed in the alignment, typing, and assembly of a read (not a pair) | 2 | `--num-editdist 0`
**-p / \--threads** | Number of threads to be used | 1 | `-p 4`
**\--assembly** | Perform assembly of each locus of interest | disabled | `--assembly`
**\--verbose** | Provide more information | disabled | `--verbose`
**-h / \--help** | Output help message | disabled | `--help`

Examples:
```bash
 $ hisatgenotype_locus_samples.py -p 4 --region-list hla --read-dir test_input --out-dir test_output
 $ hisatgenotype_locus_samples.py -p 4 --region-list hla.A,hla.B --read-dir test_input --out-dir test_output
 $ hisatgenotype_locus_samples.py -p 4 --region-list codis --read-dir test_input --out-dir test_output
```