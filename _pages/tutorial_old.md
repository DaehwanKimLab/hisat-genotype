---
layout: page
title: Tutorials Pre v1.1
hide: True
permalink: /tutorials_old/
order: 5
share: false
---

{: .no_toc}

- TOC
{:toc}

# HLA Typing
## General Instructions
Using the human reference genome, referred to as the linear reference (e.g. GRCh38 and hg38), for genomic analysis would be rather straightforward if our variants were uniformly distributed with only one nucleotide difference every 1,000 nucleotides, which most of the currently used alignment programs could handle with great accuracy. However, the distribution of variants is not uniform. Some genomic regions such as HLA genes and DNA fingerprinting loci are highly polymorphic. So using the reference genome for analyzing such highly polymorphic regions may not be the most effective approach, and this is where our graph reference comes into play. Here we describe one case, HLA-typing, where our graph reference/alignment method outperforms currently used approaches.

The tutorial requires a 64-bit computer running either Linux or Mac OS X and 8 GB of RAM. All the commands used should be run from the Unix shell prompt within a terminal window and are prefixed with a '$' character.

This tutorial is under active development and subject to change at any time.

## Initial Setup
We use HISAT2 for graph representation and alignment, which is currently the most practical and quickest program available. We refer to hisat-genotype-top as our top directory where all of our programs are located. hisat-genotype-top is a place holder that you can change to whatever name you’d like to use.

In order to install HISAT2, please run the following commands.

```bash
 $ git clone https://github.com/DaehwanKimLab/hisat2 hisat-genotype
 $ cd hisat-genotype
 hisat-genotype-top$ git checkout hisatgenotype_v1.0.1_beta
 $ make hisat2-align-s hisat2-build-s hisat2-inspect-s
```

Add the above directory (hisat-genotype) to your PATH environment variable (e.g. ~/.bashrc) to make the binaries we just built above and other python scripts available everywhere:

```bash
 export PATH=hisat-genotype:hisat-genotype/hisatgenotype_scripts:$PATH
 export PYTHONPATH=hisat-genotype/hisatgenotype_modules:$PYTHONPATH
```

After that, you may have to run the following command to reflect the change,

```bash
 $ source ~/.bashrc
```

Create a directory where we perform our analysis for HLA typing and assembly, which we will refer to as hla-analysis. hla-analysis is a place holder that you can change to whatever name you’d like to use.

```bash
 $ mkdir hla-analysis
```

Change the current directory to hla-analysis.

```bash
 $ cd hla-analysis
```

Additional program requirements:

```bash
 SAMtools (version 1.3 or later)
```

## Downloading or Building a Graph Reference and Index
The graph reference we are going to build incorporates variants of numerous HLA alleles into the linear reference using a graph. The graph reference also includes some known variants of other regions of the genome (e.g. common small variants).

We provide a pre-built graph reference and index here, which you can download as follows.

```bash
 hla-analysis$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/genotype_genome_20180128.tar.gz
 hla-analysis$ tar xvzf genotype_genome_20180128.tar.gz
```

Or alternatively, if you want to build a graph reference, please refer to the Building a graph reference.

## Read Extraction
Since whole genome sequencing (WGS) data includes reads that are from the whole genome, the first step is to extract the reads that belong to the HLA genes by aligning them to the graph reference with HISAT2. The graph reference enables substantially more sensitive alignment compared to the linear reference (e.g. based on our experiment, for read extraction, we got twice as many reads than from using the reference only.)

```bash
 hla-analysis$ hisatgenotype_extract_reads.py --base genotype_genome --read-dir ILMN_reads --out-dir ILMN
```

, where ILMN_reads is a placeholder for a directory containing exemplary read files as follows. Please make sure that your files end with 1.fq.gz or 2.fq.gz.

```bash
 NA12892.extracted.1.fq.gz
 NA12892.extracted.2.fq.gz
 NA12878.extracted.1.fq.gz
 NA12878.extracted.2.fq.gz
```

Or alternatively, we provide a set of reads that we extracted for HLA genes available at ILMN.tar.gz. The reads are from Illumina Platinum Genomes, which include 17 individuals with the CEPH pedigree 1463.

```bash
 hla-analysis$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/hla/ILMN.tar.gz
 hla-analysis$ tar xvzf ILMN.tar.gz
```

## Typing and Assembly
HISAT-genotype performs both HLA typing and assembly as follows.

You can perform HLA typing and assembly for HLA-A gene on sequencing reads from a genome, NA12892 (Illumina's HiSeq 2000 platform).

```bash
 hla-analysis$ hisatgenotype_locus.py --base hla --locus-list A -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz
```

## Interpreting Output
### Typing Output

```bash
Number of reads aligned: 1507
  1 A*02:01:01:02L (count: 571)
  2 A*02:01:31 (count: 557)
  3 A*02:20:02 (count: 557)
  4 A*02:29 (count: 557)
  5 A*02:321N (count: 556)
  6 A*02:372 (count: 556)
  7 A*02:610:02 (count: 556)
  8 A*02:249 (count: 555)
  9 A*02:479 (count: 555)
 10 A*02:11:01 (count: 554)
```

The above lines show the top ten alleles that the most number of reads are mapped to or compatible with. For example, the allele first ranked, A*02:01:01:02L, is compatible with 571 reads. This raw estimate based on the number of reads should not be used to determine the two true alleles because the alleles that resemble both but are not true alleles often tend to be compatible with more reads than either of the true alleles. Thus, we apply a statistical model to identify the two true alleles as described here.

```bash
Abundance of alleles
  1 ranked A*02:01:01:01 (abundance: 54.32%)
  2 ranked A*11:01:01:01 (abundance: 45.20%)
  3 ranked A*24:33 (abundance: 0.48%)
``` 

The above rankings show the top three alleles that are most abundant in the sample. Normally, the top two alleles in this estimate (e.g. A*02:01:01:01 and A*11:01:01:01) are considered as the two alleles that best match a given sequencing data.

### Assembly Output

These are screenshots of the assembly output as a result of the command line from Typing and Assembly. The actual output is available here in HTML format. Please note that currently only Google Chrome browser supports the HTML file.

HLA Assembly NA12892 1.png

The first two bands are two alleles predicted by HISAT-genotype, in this case A*02:01:01:01 in green and A*11:01:01:01 in yellow. Below are shorter bands indicating read alignments whose color is determined according to their compatibility with either allele. If reads are compatible with both alleles, they are shown in white.

HLA Assembly NA12892 2.png

As above, the first two bands are two alleles predicted by HISAT-genotype, and the next two bands are two alleles assembled by HISAT-genotype. In most of the cases, the predicted alleles are the same as the assembled alleles.

## Analyzing Many Samples
Here we will explain how to analyze hundreds or more samples on a compute cluster, in particular, SLURM.

Since whole genome sequencing (WGS) data includes reads that are from the whole genome, the first step is to extract the reads that belong to the HLA genes by aligning them to the graph reference with HISAT2. The read extraction step consumes the most time during the whole analysis.

You can create a script named run_extract.sh that you will run with sbatch as follows. This command assumes that we use one high-performance machine with 48 CPU cores and 512 GB or more of RAM to extract reads from all the samples. We use 48 CPU cores within the machine to extract reads from 48 samples in parallel.

```bash
 #!/bin/bash -l
 #SBATCH --job-name=HLA.extract
 #SBATCH --nodes=1
 #SBATCH --cpus-per-task=48
 #SBATCH --mem=480G
 #SBATCH --partition=lrgmem
 #SBATCH --time=100:0:0 
 #SBATCH --workdir=hla-analysis
 hisatgenotype_extract_reads.py --base genotype_genome --database-list hla --read-dir genomes --out-dir genomes.HLA -p 48 
```

Then run the following command,

```bash
 hla-analysis$ sbatch run_extract.sh  
```

Or alternatively, if you have access to multiple high performance machines, for example, four systems, then use the following run_extract.sh, and change the --job-range option for each system.

```bash
 #!/bin/bash -l
 #SBATCH --job-name=HLA.extract
 #SBATCH --nodes=1
 #SBATCH --cpus-per-task=48
 #SBATCH --mem=480G
 #SBATCH --partition=lrgmem
 #SBATCH --time=100:0:0 
 #SBATCH --workdir=hla-analysis
 hisatgenotype_extract_reads.py --base genotype_genome --database-list hla --read-dir genomes --out-dir genomes.HLA -p 48 --job-range 0,4
```

Then run the following command,

```bash
 hla-analysis$ sbatch run_extract.sh
 hla-analysis$ sbatch run_extract.sh (after changing the script to use --job-range 1,4)
 hla-analysis$ sbatch run_extract.sh (after changing the script to use --job-range 2,4)
 hla-analysis$ sbatch run_extract.sh (after changing the script to use --job-range 3,4)
```

These commands will enable four high-performance computers to extract reads from 192 (= 48 x 4) samples in parallel.

Perform typing for HLA-A, HLA-B, HLA-C, HLA-DQA1, HLA-DQB1, and HLA-DRB1 using the following command,

```bash
 hla-analysis$ hisatgenotype_locus_samples.py -p 4 --region-list hla.A,hla.B,hla.C,hla.DQA1,hla.DQB1,hla.DRB1 --read-dir genomes.HLA > CP_HLA.txt
```

Perform assembly for HLA-A gene or other genes (this currently runs for one gene at a time) as follows,

```bash
 hla-analysis$ hisatgenotype_locus_samples.py -p 4 --assembly --region-list hla.A --read-dir genomes.HLA --out-dir CP_A > CP_HLA_A.txt
```


# DNA Fingerprinting
## Typing and Assembly
HISAT-genotype performs both typing and assembly for 13 DNA fingerprinting loci as follows.

Change the current directory to codis-analysis, where codis stands for Combined DNA Index System.

```bash
 hisat-genotype-top$ cd codis-analysis
```

After that you can perform CODIS typing and assembly for D13S317 locus on sequencing reads from a genome, NA12892 (Illumina's HiSeq 2000 platform).

```bash
 codis-analysis$ hisatgenotype_locus.py --base codis --locus-list D13S317 --num-editdist 2 -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz
```

## Internal Usage

```bash
 codis-analysis$ git clone https://github.com/infphilo/hisatgenotype_db
 codis-analysis$ hisat-genotype-top/hisatgenotype_modules/hisatgenotype_extract_codis_data.py
 (Manually fix some repeat structures in codis.dat)
 codis-analysis$ mv codis.dat hisatgenotype_db/CODIS/
 codis-analysis$ hisat-genotype-top/hisatgenotype_modules/hisatgenotype_convert_codis.py
 codis-analysis$ mv *_gen.fasta hisatgenotype_db/CODIS/fasta/;  mv *_gen.msf hisatgenotype_db/CODIS/msf/
 codis-analysis$ cd hisatgenotype_db; git commit -a -m "."; git push origin master; cd ..
```

# CYP Typing
## Typing and Assembly

HISAT-genotype performs both typing and assembly for CYP gene family as follows.

Change the current directory to cyp-analysis.

```bash
 hisat-genotype-top$ cd cyp-analysis
```

After that you can perform CYP typing and assembly for CYP2D6 on sequencing reads from a genome, NA12892 (Illumina's HiSeq 2000 platform).

```bash
 cyp-analysis$ hisatgenotype_locus.py --base cyp --locus-list 2D6 -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz
```

## Internal Usage

```bash
 cyp-analysis$ git clone https://github.com/infphilo/hisatgenotype_db
 cyp-analysis$ hisat-genotype-top/hisatgenotype_modules/hisatgenotype_extract_cyp.py
 cyp-analysis$ mv *_gen.fasta hisatgenotype_db/CYP/fasta/;  mv *_gen.msf hisatgenotype_db/CYP/msf/
 cyp-analysis$ cd hisatgenotype_db; git commit -a -m "."; git push origin master; cd ..
```