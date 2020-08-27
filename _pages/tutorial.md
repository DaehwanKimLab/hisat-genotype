---
layout: page
title: Tutorials
permalink: /tutorials/
order: 5
share: false
---

{: .no_toc}

- TOC
{:toc}

# Past Tutorials

[Tutorial Pre v1.1]({{ site.baseurl }}{% link _pages/tutorial_old.md %})

# HLA Typing
## General Instructions
Using the human reference genome, referred to as the linear reference (e.g. GRCh38 and hg38), for genomic analysis would be rather straightforward if our variants were uniformly distributed with only one nucleotide difference every 1,000 nucleotides, which most of the currently used alignment programs could handle with great accuracy. However, the distribution of variants is not uniform. Some genomic regions such as HLA genes and DNA fingerprinting loci are highly polymorphic. So using the reference genome for analyzing such highly polymorphic regions may not be the most effective approach, and this is where our graph reference comes into play. Here we describe one case, HLA-typing, where our graph reference/alignment method outperforms currently used approaches.

The tutorial requires a 64-bit computer running either Linux or Mac OS X and 8 GB of RAM. All the commands used should be run from the Unix shell prompt within a terminal window and are prefixed with a '$' character.

This tutorial is under active development and subject to change at any time.

## Initial Setup
We use HISAT2 for graph representation and alignment, which is currently the most practical and quickest program available. We refer to hisatgenotype as our top directory where all of our programs are located. hisatgenotype is a place holder that you can change to whatever name you’d like to use.

In order to install HISAT2, please run the following commands for automated installation (Mac/Linux using Bash). For manual installation please see the Manual.

```bash
git clone https://github.com/DaehwanKimLab/hisat-genotype.git ~/hisatgenotype
cd hisatgenotype
bash setup.sh -r
```

Create a directory where we perform our analysis for HLA typing and assembly, which we will refer to as hla-analysis. hla-analysis is a place holder that you can change to whatever name you’d like to use.

```bash
mkdir hla-analysis
```

Change the current directory to hla-analysis.

```bash
cd hla-analysis
```

Additional program requirements:

```bash
SAMtools (version 1.3 or later)
```

## Downloading or Building a Graph Reference and Index
The graph reference we are going to build incorporates variants of numerous HLA alleles into the linear reference using a graph. The graph reference also includes some known variants of other regions of the genome (e.g. common small variants).

We provide a pre-built graph reference and index which was automatically downloaded in the last step.

Or alternatively, if you want to build a graph reference, please refer to the Building a graph reference.

## Typing and Assembly
HISAT-genotype performs both HLA typing and assembly and now had read extraction and database management built in. The following steps are a simple methodology for typing HLA gene using a test set of reads we provide. The reads are from Illumina Platinum Genomes, which include 17 individuals with the CEPH pedigree 1463.

The following scrip downloads the test dataset. Alternatively, you can use your own data in place of the ILMN files.

```bash
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/hla/ILMN.tar.gz
tar xvzf ILMN.tar.gz
```

You can perform HLA typing and assembly for HLA-A gene on sequencing reads from a genome, NA12892 (Illumina's HiSeq 2000 platform).

```bash
hisatgenotype --base hla --locus-list A -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz
```

Even though the ILMN data is already preextracted for HLA, HISAT-genotype will attempt to extract the reads and place them in a new folder. Note that you can add more loci to the `--locus-list` option above. EX `--locus-list A,B,C,DRB1,DQA1`

You can add the `--assembly` option to get HISAT-genotype to assemble the reads into alleles.

## Interpreting Output
Output will be found by default in a folder called *hisatgenotype_out*.
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

The above lines show the top ten alleles that the most number of reads are mapped to or compatible with. For example, the allele first ranked, A\*02:01:01:02L, is compatible with 571 reads. This raw estimate based on the number of reads should not be used to determine the two true alleles because the alleles that resemble both but are not true alleles often tend to be compatible with more reads than either of the true alleles. Thus, we apply a statistical model to identify the two true alleles as described here.

```bash
Abundance of alleles
  1 ranked A*02:01:01:01 (abundance: 54.32%)
  2 ranked A*11:01:01:01 (abundance: 45.20%)
  3 ranked A*24:33 (abundance: 0.48%)
```

The above rankings show the top three alleles that are most abundant in the sample. Normally, the top two alleles in this estimate (e.g. A\*02:01:01:01 and A\*11:01:01:01) are considered as the two alleles that best match a given sequencing data.

By running the following code the output can be reformated into a csv file for easier use:

```bash
hisatgenotype_toolkit parse-results --csv --in-dir hisatgenotype_out
```

### Assembly Output
When using the *--assembly* option in HISATgenotype two additional files are produced. A png file with a graphical representation of the alleles typed and assembled, and a fasta file with the sequence of the contigs reported. The contigs will be annotated with the allele they derive from if the algorithm was able to match it to a seqence. 

If used with the example files, the first two bands in the png are two alleles predicted by HISAT-genotype, in this case A\*02:01:01:01 in green and A\*11:01:01:01 in yellow. Below are shorter bands indicating read alignments whose color is determined according to their compatibility with either allele. If reads are compatible with both alleles, they are shown in white.

The middle two bands are two alleles predicted by HISAT-genotype, and the next two bands are two alleles assembled by HISAT-genotype. In most of the cases, the predicted alleles are the same as the assembled alleles.
