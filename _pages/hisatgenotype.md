---
layout: page
title: Main
permalink: /main/
order: 1
share: false
---

**HISAT-genotype** is a next-generation genomic analysis software platform capable of assembling and genotyping human genes and genomic regions. Thie software leverages **HISAT2**s graph FM index and graph alignemnt algorithm to align reads to a specially constructed graph genome. An Expectation-Maximization (EM) algorithm finds the maximum likelihood estimates for each gene allele and a guided de Bruijn graph is used to construct the allele sequences.  
[__Get HISAT-genotype__](https://github.com/DaehwanKimLab/hisat-genotype)

### News and Releases:
+ __*08/27/2020*: HISAT-genotype 1.3.1 release__
    - setup.sh script now works properly and has been given option to predownload required indicies
    - parse-results script has been added to hisatgenotype\_toolkit and can be used to parse the report files into something more managable in downstream applications and can limit the reported allele digits/files
    - *--assembly* option now will also output a fasta file with the constructed allele sequences
    - Added -z/--index\_dir option to hisatgenotype that will tell hisatgenotype where to look for/download indicies
    - General bug fixes and improvements to stability
<br>

+ __*05/20/2020*: HISAT-genotype 1.3.0 release__
    - Python 3.7 is now the base python HISATgenotype is built upon
    - Core scripts have been further consolidated and collapsed with redundant code being removed
    - hisatgenotype\_script folder has been renamed etc and now contains mostly accessory scripts
    - debug argument is now fully random when testing allele combinations
    - Each python module now has a flag for degugging and running sanity checks
    - Results are now sorted alphabetically
    - Added tool to concolidate and simplify results from HISAT-genotype
    - Added Setup script to aid HISAT-genotype initial set-up
    - HISAT-genotype now uses the latest version of HISAT2 and tracks versions
    - Results are output to a result file by default
    - Reworked Github repo structure for easy access to different version
<br>

+ __*01/23/2020*: HISAT-genotype 1.1.3 release__
    - Incremental change to address bugs with using -U option
<br>

+ __*07/11/2019*: HISAT-genotype 1.1.2-beta release__
    - Wrapper has been added to HISAT-genotype. `hisatgenotype` and `hisatgenotype_toolkit` runs the entirety of HISAT-genotype
    - Argument names for scripts have been made more consistent
    - HISAT-genotype code has been consolidated into modules. All core functionality is now in modules
    - Scripts for HISAT-genotype have been consolidated into single directory and are runnable through `hisatgenotype_toolkit`
    - Changes made to debug mode to prevent overwritting of debug results through each iteration
    - Advanced arguments are hidden by default in `hisatgentoype`  
<br>

+ __*02/01/2018*: HISAT-genotype 1.0.1-beta release__
    - HLA assembly results are reported in the PDF format, instead of the HTML format.  
<br>

+ __*06/08/2017*: HISAT-genotype 1.0.0-beta release__ (first release)
    - HISAT-genotype currently supports HLA typing, the discovery of novel HLA alleles, DNA fingerprinting analysis, and CYP2D6 typing.
    - The platform is currently designed for processing whole genome sequencing reads produced by Illumina platforms (e.g. 100 to 300 bps) and works great with paired-end reads and at least 20x coverage.
    - We plan to create consortia to work together in enabling the platform to analyze the whole human genome with its >20,000 genes.  
<br>

+ __*06/08/2017*: The HISAT-genotype source code is available in a [public GitHub repository](https://github.com/DaehwanKimLab/hisat-genotype).__  
<br>

### Other links:
[HISAT2](http://ccb.jhu.edu/software/hisat2/index.shtml), [github](https://github.com/DaehwanKimLab/hisat2)  
[HISAT-genotype FTP](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data)
