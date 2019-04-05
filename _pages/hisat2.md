---
layout: page
title: Main
permalink: /main/
order: 1
share: false
---

**HISAT-genotype** is a next-generation genomic analysis software platform capable of assembling and genotyping human genes and genomic regions. Thie software leverages **HISAT2**s graph FM index and graph alignemnt algorithm to align reads to a specially constructed graph genome. An Expectation-Maximization (EM) algorithm finds the maximum likelihood estimates for each gene allele and a guided de Bruijn graph is used to construct the allele sequences.

### HISAT-genotype 1.0.1-beta release 2/1/2018
 * HLA assembly results are reported in the PDF format, instead of the HTML format.

### HISAT-genotype 1.0.0-beta release 6/8/2017 (first release)
 * HISAT-genotype currently supports HLA typing, the discovery of novel HLA alleles, DNA fingerprinting analysis, and CYP2D6 typing.
 * The platform is currently designed for processing whole genome sequencing reads produced by Illumina platforms (e.g. 100 to 300 bps) and works great with paired-end reads and at least 20x coverage.
 * We plan to create consortia to work together in enabling the platform to analyze the whole human genome with its >20,000 genes.

### The HISAT-genotype source code is available in a [public GitHub repository](https://github.com/DaehwanKimLab/hisat-genotype) (6/8/2017).
