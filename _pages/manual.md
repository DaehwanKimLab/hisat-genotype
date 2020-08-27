---
layout: page
title: Manual
permalink: /manual/
order: 3
share: false
---

{: .no_toc}

- TOC
{:toc}

# Introduction

Using the human reference genome, referred to as the linear reference (e.g. GRCh38 and hg38), for genomic analysis would be rather straightforward if our variants were uniformly distributed with only one nucleotide difference every 1,000 nucleotides, which most of the currently used alignment programs could handle with great accuracy. However, the distribution of variants is not uniform. Some genomic regions such as HLA genes and DNA fingerprinting loci are highly polymorphic. So using the reference genome for analyzing such highly polymorphic regions may not be the most effective approach, and this is where our graph reference comes into play.

---

# HISAT-genotype Set-up

We use [HISAT2] for graph representation and alignment, which is currently the most practical and quickest program available. We refer to hisat-genotype as our top directory where all of our programs are located. hisat-genotype is a place holder that you can change to whatever name you’d like to use.

[HISAT2]:          http://ccb.jhu.edu/software/hisat2

## Requirements
### List
+ Python 3<
+ Samtools 1.3<
+ C++ Compiler

### Description
HISAT-genotype is encoded in __PYTHON 3__ and contains standard python libraries. A
python 3.7 release is recommended. Additional software required is [Samtools]
version 1.3 or later.

Building HISAT2 from source requires a GNU-like environment with GCC, GNU Make
and other basics.  It should be possible to build HISAT2 on most vanilla Linux
installations or on a Mac installation with [Xcode] installed.  HISAT2 can
also be built on Windows using [Cygwin] or [MinGW] (MinGW recommended). For a
MinGW build the choice of what compiler is to be used is important since this
will determine if a 32 or 64 bit code can be successfully compiled using it. If
there is a need to generate both 32 and 64 bit on the same machine then a multilib
MinGW has to be properly installed. [MSYS], the [zlib] library, and depending on
architecture [pthreads] library are also required. We are recommending a 64 bit
build since it has some clear advantages in real life research problems. In order
to simplify the MinGW setup it might be worth investigating popular MinGW personal
builds since these are coming already prepared with most of the toolchains needed.

[Cygwin]:   http://www.cygwin.com/
[MinGW]:    http://www.mingw.org/
[MSYS]:     http://www.mingw.org/wiki/msys
[zlib]:     http://cygwin.com/packages/mingw-zlib/
[pthreads]: http://sourceware.org/pthreads-win32/
[GnuWin32]: http://gnuwin32.sf.net/packages/coreutils.htm
[Download]: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/
[Xcode]: https://developer.apple.com/Xcode
[Samtools]: https://www.htslib.org

## Automated Install - (Mac/Linux)
This is a simple automated method for installing and getting HISAT-genotype setup on a Linux/Mac system using Bash.
Replace the ~ which whichever directory you'd like to store HISAT-genotype in. The -r option in the setup script will pre-download all of the basic requited indicies into the HISAT-genotype source directory. If you want to manually direct where these indicies are downloaded use the -x option followed by the absolute path to the desired location.

```bash
git clone https://github.com/DaehwanKimLab/hisat-genotype.git ~/hisatgenotype
cd hisatgenotype
bash setup.sh -r
```

Check if setup was successful using the following commands:
```bash
hisatgenotype --help
hisat2 --help
```

If there is an error then something did not set-up properly and you'll need to run a manual install

### setup.sh options
* **\-h** | *Default* : *none*
> Show help screen

* **\-b** | *Default* : False
> Do not try to automatically add HISAT-genotype and HISAT2 to your Path environment

* **\-r** | *Default* : False
> Pre-download the base indicies for HISAT-genotype

* **\-x** | *Default* : *[PATH_TO_HISATGENOTYPE]*/indicies
> If -r option is set, set desired location for indicies if different than default

***NOTE:*** If you want to predownload all indicies before running HISAT-genotype but after HISAT-genotype install and not automatically during HISAT-genotype run, you can use `bash setup.sh -brx PATH_TO_DIR` while in HISAT-genotype install directory.

## Manual Install - (Mac/Linux/Windows)
### Downloading HISAT-genotype and Building HISAT2 from source
This download example will place HISAT-genotype in your home (~) directory if you are using a linux system. 
Change the ~ to whichever directory you desire if this is not the behavior you want.

```bash
git clone --recurse-submodules https://github.com/DaehwanKimLab/hisat-genotype ~/hisatgenotype
cd ~/hisatgenotype/hisat2

$ make
```

### Adding HISAT-genotype to PATH
Add the above directory (hisat-genotype) to your PATH environment variable
(e.g. ~/.bashrc) to make the binaries built above and other python scripts
available everywhere:

```bash
$ export PATH=~/hisatgenotype:~/hisatgenotype/hisat2:$PATH
$ export PYTHONPATH=~/hisatgenotype/hisatgenotype_modules:$PYTHONPATH
```

---

# Running HISAT-genotype 

## Past Manuals

[Manual Pre v1.1]({{ site.baseurl }}{% link _pages/manual_old.md %})

## hisatgenotype - Analysis of a whole human genome
The hisatgenotype.py python script will analyze a whole human genome using whole genome sequencing reads. It will align reads to genotype genome, extract reads belonging to each locus of interest, and perform typing and assembly. 

Usage:
```bash
$ hisatgenotype -x [GENOME] --base [GENE_GROUP] -z [INDEX_DIR] [OPTIONS] -1 [FASTQ_PAIR1] -2 [FASTQ_PAIR2]
```

### Standard Options

* **-x / \--ref-genome** | *Default* : *None* 
> Base name for genome index if not genotype_genome. Generally reserved for custom graph genomes end user may want to use  
> Example: `-x custome_genome`

* **\--base / \--base-fname** | *Default* : (empty) all databases 
> Base file name for index, variants, haplotypes, etc. (e.g. hla, rbg, codis). This will be anything in the hisatgenotype_db folder  
> Example: `--base hla`

* **\--locus-list** | *Default* : (empty) all genes 
> A comma-separated list of gene names
> Example: `--locus-list A,B,C,DRB1,DQA1,DQB1`

* **\-z / \--index_dir** | *Default* : pre-downloaded directory or link file (hg_ix.link) location
> Set location for the indecies HISATgenotype requires
> Example: `-z ~/hisatgenotype/indicies`

* **-f / \--fasta** | *Default* : `False`  
> Bool to indicate if reads are provided in FASTA format  
> Example: `-f`

* **-U** | *Default* : *None*  
> Single-end read file name  
> Example: `-U read.fq.gz`

* **-1** | *Default* : *None*  
> Paired-end read file name 1  
> Example: `-1 read.1.fq.gz`

* **-2** | *Default* : *None*  
> Paired-end read file name 2  
> Example: `-2 read.2.fq.gz`

* **\--keep-alignemnt** | *Default* : `False`  
> Bool to keep the alignment BAM file if typing from FASTQ(A). If typing from a BAM file this is irrelevent.  
> Example: `--keep-alignment`

* **\--in-dir** | *Default* : *None*  
> Directory HISATgenotpye will search for FASTQ(A) files and batch process. HISAT-genotype will attempt to automatically pair the files if `--single-end` isn't set.
> Try to have the names be similar between the pairs with a single difference to make it easier for HISAT-genotype to pair the files (e.g. *hg_granulocyte_samp1_L.fastq* and *hg_granulocyte_samp1_R.fastq*)  
> Example: `--in-dir input_fastq_dir`

* **\--out-dir** | *Default* : `/hisatgenotype_out`  
> Directory where all resulting files will be placed. This includes the typing results, BAM files, and assembly files.  
> Example: `--out-dir results`

* **\--bamfile** | *Default* : *None*  
> BAM file name if using already aligned reads  
> Example: `--bamfile hg_granulocyte_samp1.bam`

* **\--single-end** | *Default* : `False`  
> Bool to indicate if file(s) is/are single ended. Only needed with `--in-dir` or `--bamfile`  
> Example: `--single-end`

* **\--assembly** | *Default* : disabled  
> Perform assembly of each locus of interest  
> Example: `--assembly`

* **\--assembly-name** | *Default* : `assemply_graph`  
> Assembly base file name to use  
> Example: `--assembly-name hg_granulocyte`

* **\--assembly-verbose** | *Default* : `False` 
> Bool to output additional assembly information  
> Example: `--assembly-verbose`

* **-p / \--threads** | *Default* : 1   
> Number of threads to be used  
> Example: `-p 4`

* **\--verbose** | *Default* : `False`  
> Provide more information
> Example: `--verbose`

* **-h / \--help** | *Default* : `False` 
> Output help message  
> Example: `--help `

* **\--advanced-help** | *Default* : `False` 
> Output help message for advanced options  
> Example: `--advanced-help`

### Advanced Options

* **\--keep-extract** | *Default* : `False`    
> Bool to keep extracted read fastq files  
> Example: `--keep-extract`

* **\--build-base** | *Default* : `False`
> Build the indexes listed in \--base-fname

* **\--aligner** | *Default* : `hisat2`
> Set aligner to use (ex. hisat2, bowtie2)

* **\--linear-index** | *Default* : `False`
> Use linear index

* **\--num-mismatch** | *Default* : `0`
> Maximum number of mismatches per read alignment to be considered

* **\--inter-gap**
> Maximum distance for variants to be in the same haplotype

* **\--intra-gap**
> Break a haplotype into several haplotypes

* **\--whole-haplotype**
> Include partial alleles (e.g. A_nuc.fasta)

* **\--min-var-freq** | *Default* : `0.0`
> Exclude variants whose freq is below than this value in percentage

* **\--ext-seq** | *Default* : `0`
> Length of extra sequences flanking backbone sequences

* **\--leftshift** | *Default* : `False`
> Shift deletions to the leftmost

* **\--suffix** | *Default* : `fq.gz`
> Read file suffix

* **\--simulation** | *Default* : `False`
> Simulated reads (Default: False)

* **\--pp, \--threads-aprocess** | *Default* : `1`
> Number of threads a process

* **\--max-sample** | *Default* : `sys.maxint`
> Number of samples to be extracted

* **\--job-range** | *Default* : `0,1`
> two numbers (e.g. 1,3)

* **\--extract-whole** | *Default* : `False`
> Extract all reads

* **\--no-partial** | *Default* : `False`
> Include partial alleles (e.g. A_nuc.fasta)

* **\--simulate-interval** | *Default* : `10`
> Reads simulated at every these base pairs

* **\--read-len** | *Default* : `100`
> Length of simulated reads

* **\--fragment-len** | *Default* : `350
> Length of fragments

* **\--best-alleles** | *Default* : `False`
> *Placeholder*

* **\--random-seed** | *Default* : `1`
> A seeding number for randomness

* **\--num-editdist** | *Default* : `2`
> Maximum number of mismatches per read alignment to be considered

* **\--perbase-errorrate** | *Default* : `0.0`
> Per basepair error rate in percentage when simulating

* **\--perbase-snprate** | *Default* : `0.0`
> Per basepair SNP rate in percentage when simulating

* **\--skip-fragment-regions** | *Default* : *None*
> A comma-separated list of regions from which no reads originate, e.g., 500-600,1200-1400

* **\--verbose-level** | *Default* : `0`
> also print some statistics to stderr

* **\--no-error-correction** | *Default* : `False`
> Correct sequencing errors

* **\--only-locus-list** | *Default* : (empty) all genes
> A comma-separated list of genes

* **\--discordant** | *Default* : `False`
> Allow discordantly mapped pairs or singletons

* **\--type-primary-exons** | *Default* : `False`
> Look at primary exons first

* **\--keep-low-abundance-alleles** | *Default* : `False`
> Do not remove alleles with low abundance while performing typing

* **\--display-alleles** | *Default* : *None*
> A comma-separated list of alleles to display in HTML

* **\--debug** | *Default* : *None*
> Test database or code  
> (options: basic, pair, full, single-end, test_list, test_id)  
> Example: `--debug test_id:10,basic`

## hisatgenotype_toolkit - Individual Scripts and Tools for Custom Pipelines

*Work In Progress to document*

Usage:
```bash
$ hisatgenotype_toolkit <BASE_TOOL> [TOOL_OPTIONS]
```

### build-genome        
*(hisatgenotype_build_genome.py)*

### call-variants       
*(hisatgenotype_call_variants.py)*

### convert-codis       
*(hisatgenotype_convert_codis.py)*

### extract-RBG         
*(hisatgenotype_extract_RBG.py)*

### extract-codis-data  
*(hisatgenotype_extract_codis_data.py)*

### extract-cyp-data    
*(hisatgenotype_extract_cyp_data.py)*

### extract-reads       
*(hisatgenotype_extract_reads.py)*

### extract-vars        
*(hisatgenotype_extract_vars.py)*

### legacy              
*(hisatgenotype_legacy.py)*

### locus               
*(hisatgenotype_locus.py)*

### locus-samples       
*(hisatgenotype_locus_samples.py)*

### parse-results
*(hisatgenotype_parse_results.py)*

* **\--in-dir** | *Default* : Current Directory
> Input directory where HISAT-genotype *.report* files can be found

* **\-t / \--trim**  | *Default* : 4/All
> Trim the reported alleles to 1 (A\*01), 2 (A\*01:01), 3 (A\*01:01:01), or all (A\*01:01:01:01) fields

* **\--csv** | *Default* : False
> Return the formated output in a tab deliminated csv file (tsv) for use in spreadsheets

* **\--output-file** | *Default* : HG_report_results.csv
> Name of csv file


