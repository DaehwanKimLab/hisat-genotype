#!/usr/bin/env python
# --------------------------------------------------------------------------- #
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>                            #
#                                                                             #
# This file is part of HISAT-genotype. It contains the core algorithms for    #
# pre-processing and processing reads for the core.                           #
#                                                                             #
# HISAT-genotype is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# HISAT-genotype is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with HISAT-genotype.  If not, see <http://www.gnu.org/licenses/>.     #
# --------------------------------------------------------------------------- #

import sys
import os
import subprocess
import re
import resource
import inspect
import random
import glob
import multiprocessing
import hisatgenotype_typing_common as typing_common
import hisatgenotype_validation_check as validation_check

""" Flag to turn on file debugging to run sanity checks """
SANITY_CHECK = True
# --------------------------------------------------------------------------- #
# Scripts for use in extract_vars.                                            #
#                                                                             #
# Builds a HISAT-genotype database using fasta, MSF, and data file from       #
# an allele database                                                          #
# --------------------------------------------------------------------------- #
"""
Mapping from base pair to a location in MSF format given an MSF sequence
"""
def create_map(seq):
    seq_map = {}
    count = 0
    for i in range(len(seq)):
        bp = seq[i]
        if bp in '.EN~':
            continue
        assert bp in "ACGT", '%s not a valid basepair' % bp
        seq_map[count] = i
        count += 1
    return seq_map

"""
Given list of MSF sequences, build a "consesus sequence based on most 
  common nucleotide sequences present
"""
def create_consensus_seq(seqs,
                         seq_len,
                         min_var_freq,
                         remove_empty = True):
    consensus_freq = [[0, 0, 0, 0, 0] for i in range(seq_len)]
    seq_coverage = [0 for i in range(seq_len)]
    for i in range(len(seqs)):                
        seq = seqs[i]
        if len(seq) != seq_len:
            continue                    
        for j in range(seq_len):
            nt = seq[j]
            if nt == "~":
                continue
            
            seq_coverage[j] += 1 # coverage of sequence
            assert nt in "ACGT.EN", "Nucleotide %s not supported" % nt
            if nt == 'A':
                consensus_freq[j][0] += 1
            elif nt == 'C':
                consensus_freq[j][1] += 1
            elif nt == 'G':
                consensus_freq[j][2] += 1
            elif nt == 'T':
                consensus_freq[j][3] += 1
            else:
                assert nt in ".EN"
                consensus_freq[j][4] += 1

    for j in range(len(consensus_freq)):
        for k in range(len(consensus_freq[j])):
            consensus_freq[j][k] /= float(seq_coverage[j])
            consensus_freq[j][k] *= 100.0

    consensus_seq = ""
    has_empty = False
    for c in range(len(consensus_freq)):
        freq = consensus_freq[c]
        A, C, G, T, E = freq
        # No alleles have bases at this particular location
        if E >= 100.0:
            has_empty = True
            consensus_seq += 'E'
            continue
        if E >= 100.0 - min_var_freq: ## TODO: Revisit necessity of this
            idx = 4                   ## Goal is to include only viable seqs
        else:
            idx = freq.index(max([A, C, G, T]))
        assert idx < 5
        consensus_seq += "ACGT."[idx]
    consensus_seq = ''.join(consensus_seq)

    # Remove dots (deletions) (Empty spaces only at this case)
    skip_pos = set()
    if has_empty and remove_empty:
        for seq_i in range(len(seqs)):
            seqs[seq_i] = list(seqs[seq_i])
        for i in range(len(consensus_seq)):
            if consensus_seq[i] != 'E':
                continue
            skip_pos.add(i)
            for seq_i in range(len(seqs)):
                if i >= len(seqs[seq_i]):
                    continue
                seqs[seq_i][i] = 'E'
        for seq_i in range(len(seqs)):
            seqs[seq_i] = ''.join(seqs[seq_i])
            seqs[seq_i] = seqs[seq_i].replace('E', '')
        consensus_seq = consensus_seq.replace('E', '')

    # Convert a list form of consensus_freq to a dictionary form
    temp_freq = []
    for j in range(len(consensus_freq)):
        if j in skip_pos:
            continue
        freq_dic = {}
        for k in range(len(consensus_freq[j])):
            freq = consensus_freq[j][k]
            if freq <= 0.0:
                continue
            nt = "ACGT."[k]                    
            freq_dic[nt] = freq
        temp_freq.append(freq_dic)
    consensus_freq = temp_freq

    assert len(consensus_seq) == len(consensus_freq)                
    return consensus_seq, consensus_freq

"""
Left-shift deletions if poissble
"""
def leftshift_deletions(backbone_seq, 
                        seq, 
                        debug = False):
    if len(seq) != len(backbone_seq):
        return seq
    seq = list(seq)
    seq_len = len(seq)
    bp_i = 0
    # Skip the first deletion
    while bp_i < seq_len:
        if seq[bp_i] in "ACGT":
            break
        bp_i += 1

    while bp_i < seq_len:
        bp = seq[bp_i]
        if bp != '.':
            bp_i += 1
            continue
        bp_j = bp_i + 1
        while bp_j < seq_len:
            bp2 = seq[bp_j]
            if bp2 != '.':
                break
            else:
                bp_j += 1

        if bp_j >= seq_len:
            bp_i = bp_j
            break

        if debug:
            # Print location of pointers in string and string seq
            print((bp_i, 
                   bp_j, 
                   backbone_seq[bp_i-10:bp_i], 
                   backbone_seq[bp_i:bp_j], 
                   backbone_seq[bp_j:bp_j+10]), 
                  file=sys.stderr)
            print((bp_i, 
                   bp_j, 
                   ''.join(seq[bp_i-10:bp_i]), 
                   ''.join(seq[bp_i:bp_j]), 
                   ''.join(seq[bp_j:bp_j+10])), 
                  file=sys.stderr)
        prev_i, prev_j = bp_i, bp_j

        while bp_i > 0 \
                and seq[bp_i-1] in "ACGT" \
                and backbone_seq[bp_j-1] in "ACGT":
            if seq[bp_i-1] != backbone_seq[bp_j-1]:
                break
            seq[bp_j-1] = seq[bp_i-1]
            seq[bp_i-1] = '.'
            bp_i -= 1
            bp_j -= 1
        bp_i = bp_j
        while bp_i < seq_len:
            if seq[bp_i] in "ACGT":
                break
            bp_i += 1

        if debug:
            # Print last pointer position in sequence
            print((prev_i, 
                   prev_j, 
                   ''.join(seq[prev_i-10:prev_i]),
                   ''.join(seq[prev_i:prev_j]), 
                   ''.join(seq[prev_j:prev_j+10])),
                  file=sys.stderr)

    return ''.join(seq)

""" Split haplotypes that include large gaps inside """
def split_haplotypes(haplotypes, intra_gap):
    split_haps = set()
    for haplotype in haplotypes:
        haplotype = haplotype.split('#')
        assert len(haplotype) > 0
        if len(haplotype) == 1:
            split_haps.add(haplotype[0])
            continue
        prev_s, s = 0, 1
        while s < len(haplotype):
            prev_locus, prev_type, prev_data = haplotype[s-1].split('-')
            locus, type, data = haplotype[s].split('-')
            prev_locus = int(prev_locus)
            locus      = int(locus)
            if prev_type == 'D':
                prev_locus += (int(prev_data) - 1)
            if prev_locus + intra_gap < locus:
                split_haps.add('#'.join(haplotype[prev_s:s]))
                prev_s = s
            s += 1
            if s == len(haplotype):
                split_haps.add('#'.join(haplotype[prev_s:s]))
    return split_haps

""" Check sequences are of equal length """
def find_seq_len(seqs):
    seq_lens = {}
    for s in range(len(seqs)):
        seq_len = len(seqs[s])
        if seq_len not in seq_lens:
            seq_lens[seq_len] = 1
        else:
            seq_lens[seq_len] += 1

    max_seq_count = 0
    for tmp_seq_len, tmp_seq_count in seq_lens.items():
        if tmp_seq_count > max_seq_count:
            seq_len       = tmp_seq_len
            max_seq_count = tmp_seq_count
    return seq_len

""" Key for sorting Variants in POS-TYPE-INFO format """
def key_varKey(x):
    type_ord = {"I" : 0, "M" : 1, "D" : 2}
    nt_order = {"A" : 0, "C" : 1, "G" : 2, "T" : 3}
    def pat2num(nt):
        if len(nt) == 0:
            return 0
        num = len(nt_order) * pat2num(nt[:-1]) + nt_order[nt[-1]]
        return num

    locus, type_, data = x.split("-")
    locus = int(locus)
    next_ = type_ord[type_]
    if type_ == "D":
        last_ = int(data)
    else:
        if len(data) > 1:
            last_ = pat2num(data)
        else:
            last_ = nt_order[data]

    return(locus, next_, last_)

""" Key for sorting haplotypes in a String""" 
def hapKey(x):
    x = x.split('#')
    x_srt, _, _ = x[0].split('-')
    x_end, x_type, x_data = x[-1].split('-')
    x_srt = int(x_srt)
    x_end = int(x_end)
    if x_type == 'D':
        x_end += (int(x_data) - 1)
    return(x_srt, x_end)

"""
Main script to extract variants:
Inputs: Fasta, MSF, dat file with allele sequences
Outputs: HISAT-genotype formatted database
"""
def extract_vars(base_fname,
                 ix_dir,
                 locus_list,
                 inter_gap,
                 intra_gap,
                 whole_haplotype,
                 min_var_freq,
                 ext_seq_len,
                 leftshift,
                 partial,
                 verbose):

    base_fullpath_name = ix_dir + base_fname

    # Download human genome and HISAT2 index
    typing_common.download_genome_and_index(ix_dir)

    # CB TODO: Make this read from a file or read the gen in the file name
    spliced_gene   = ['hla', 'rbg']
    unspliced_gene = ['codis', 'cyp' ,'rrna']
    
    # Corresponding genomic loci found by HISAT2 (reference is GRCh38)
    #   e.g. hisat2 --no-unal --score-min C,0 -x grch38/genome 
    #               -f hisatgenotype_db/HLA/fasta/A_gen.fasta
    locus_file = open(base_fullpath_name + ".locus", 'w')
    left_ext_seq_dic  = {}
    right_ext_seq_dic = {}
    genes             = {}
    gene_strand       = {}
    fasta_dname = "hisatgenotype_db/%s/fasta" % base_fname.upper()

    # Check HLA genes
    gene_names = []
    if base_fname in spliced_gene:
        fasta_fnames = glob.glob("%s/*_gen.fasta" % fasta_dname)
    else:
        assert base_fname in unspliced_gene
        fasta_fnames = glob.glob("%s/*_gen.fasta" % fasta_dname)
    for gen_fname in fasta_fnames:
        gene_name = gen_fname.split('/')[-1].split('_')[0]
        if gene_name == "hla":
            continue
        gene_names.append(gene_name)

    gene_names = typing_common.sort_genall(gene_names)
    if locus_list == []:
        locus_list = gene_names

    cigar_re = re.compile('\d+\w')
    remove_locus_list = []
    for gene in locus_list:
        aligner_cmd = ["hisat2"]
        # only include genes with no exact match to reference genome
        if base_fname not in ["cyp", "rbg"]: 
            aligner_cmd += ["--score-min", "C,-12"]
        aligner_cmd += ["--no-unal",
                        "-x", "grch38/genome",
                        "-f", "%s/%s_gen.fasta" % (fasta_dname, gene)]
        align_proc = subprocess.Popen(aligner_cmd,
                                      universal_newlines = True,
                                      stdout = subprocess.PIPE,
                                      stderr = open("/dev/null", 'w'))
        allele_id   = ""
        best_chr    = ""
        best_left   = -1
        best_right  = -1
        best_AS     = -sys.maxsize
        best_strand = ""
        for line in align_proc.stdout:
            if line.startswith('@'):
                continue
            line = line.strip()
            cols = line.split()
            
            temp_allele_id, flag, chr, left, _, cigar_str = cols[:6]
            left   = int(left) - 1
            right  = left
            cigars = cigar_re.findall(cigar_str)
            cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
            if len(cigars) > 1 or cigars[0][0] != 'M':
                continue
            for i in range(len(cigars)):
                cigar_op, length = cigars[i]
                if cigar_op in "MND":
                    right += length

            flag   = int(flag)
            strand = '-' if flag & 0x10 else '+'
            AS     = ""
            for i in range(11, len(cols)):
                col = cols[i]
                if col.startswith("AS"):
                    AS = col[5:]
            assert AS != ""
            AS = int(AS)
            if AS > best_AS:
                allele_id   = temp_allele_id
                best_chr    = chr
                best_left   = left
                best_right  = right
                best_AS     = AS
                best_strand = strand

        chr    = best_chr
        left   = best_left
        right  = best_right
        strand = best_strand
        align_proc.communicate()
        if allele_id == "":
            remove_locus_list.append(gene)
            continue
        if base_fname in spliced_gene:
            allele_name = ""
            for line in open("%s/%s_gen.fasta" % (fasta_dname, gene)):
                line = line.strip()
                if not line.startswith('>'):
                    continue
                if base_fname == 'hla':
                    tmp_allele_id, tmp_allele_name = line[1:].split()[:2]
                else:
                    tmp_allele_id   = line[1:].split()[0]
                    tmp_allele_name = line[1:].split()[0]
                if allele_id == tmp_allele_id:
                    allele_name = tmp_allele_name
                    break
        else:
            allele_name = allele_id
        assert allele_name != "" and strand != ''
        genes[gene]       = allele_name
        gene_strand[gene] = strand
        print("%s-%s's reference allele is %s on '%s' strand of chromosome %s" % \
               (base_fname.upper(), gene, allele_name, strand, chr), 
              file=sys.stderr)

        assert chr != "" and left >= 0 and right > left
        if ext_seq_len > 0:
            left_ext_seq  = ""
            right_ext_seq = ""
            left1 = max(1, left - ext_seq_len)
            left2 = max(1, left - 1)
            if left2 > 0:
                extract_seq_cmd = ["samtools", 
                                   "faidx", 
                                   "genome.fa", 
                                   "%s:%d-%d" % (chr, left1, left2)]
                extract_seq_proc = subprocess.Popen(extract_seq_cmd,
                                                    universal_newlines = True,
                                                    stdout = subprocess.PIPE,
                                                    stderr = open("/dev/null", 'w'))
                for line in extract_seq_proc.stdout:
                    if line.startswith('>'):
                        continue
                    line = line.strip()
                    left_ext_seq += line
            extract_seq_cmd = ["samtools", 
                               "faidx", 
                               "genome.fa", 
                               "%s:%d-%d" % (chr, right, right + ext_seq_len - 1)]
            extract_seq_proc = subprocess.Popen(extract_seq_cmd,
                                                universal_newlines = True,
                                                stdout = subprocess.PIPE,
                                                stderr = open("/dev/null", 'w'))
            for line in extract_seq_proc.stdout:
                if line.startswith('>'):
                    continue
                line = line.strip()
                right_ext_seq += line

            if strand == '-':
                left_ext_seq  = typing_common.reverse_complement(right_ext_seq)
                right_ext_seq = typing_common.reverse_complement(left_ext_seq)
            left_ext_seq_dic[gene]  = left_ext_seq
            right_ext_seq_dic[gene] = right_ext_seq
            
    # Extract exon information from hla.data
    gene_exons       = {}
    gene_exon_counts = {}
    
    if os.path.exists("hisatgenotype_db/VERSION"):
        dbversion = open("hisatgenotype_db/VERSION", 'r').read()
    else:
        dbversion = "NONE"

    if base_fname in spliced_gene:        
        skip = False
        look_exon_num = False
        for line in open("hisatgenotype_db/%s/%s.dat" % (base_fname.upper(), 
                                                         base_fname)):
            if line.startswith("DE"):
                if not line.split()[1][-1].isdigit():
                    allele_name = line.split()[1][:-1]
                else:
                    allele_name = line.split()[1]
                
                if allele_name.startswith("%s-" % base_fname.upper()):
                    allele_name = allele_name[len("%s-" % base_fname.upper()):]
                gene = allele_name.split('*')[0]
                
                if not gene in genes:
                    skip = True
                else:
                    skip = False
            if skip:
                continue
            if not line.startswith("FT"):
                continue
            
            if line.find("exon") != -1:
                look_exon_num = True
                if allele_name == genes[gene]:
                    exon_range = line.split()[2].split("..")
                    exon_left  = int(exon_range[0]) - 1
                    exon_right = int(exon_range[1]) - 1
                    assert exon_left >= 0
                    assert exon_left < exon_right
                    if not gene in gene_exons:
                        gene_exons[gene] = []
                    if gene in left_ext_seq_dic:
                        left_ext_seq_len = len(left_ext_seq_dic[gene])
                    else:
                        left_ext_seq_len = 0
                    gene_exons[gene].append([exon_left + left_ext_seq_len, 
                                             exon_right + left_ext_seq_len])
            elif look_exon_num:
                assert line.find("number")
                look_exon_num = False
                digits = [x for x in filter(str.isdigit, line)]
                num    = int("".join(digits))-1
                if gene not in gene_exon_counts:
                    gene_exon_counts[gene] = {}
                if num not in gene_exon_counts[gene]:
                    gene_exon_counts[gene][num] = 1
                else:
                    gene_exon_counts[gene][num] += 1
                
        for gene, exon_counts in gene_exon_counts.items():
            print(("%s exon counts:" % gene, exon_counts), 
                  file=sys.stderr)

    itr             = 0
    tmp_genes       = {}
    tmp_gene_strand = {}
    while True:
        if itr == len(locus_list):
            break
        gene = locus_list[itr]
        if gene in remove_locus_list \
                or (base_fname in spliced_gene \
                    and gene not in gene_exons):
            locus_list.pop(itr)
            continue

        if gene in genes:
            tmp_genes[gene]       = genes[gene]
            tmp_gene_strand[gene] = gene_strand[gene]
        
        itr += 1
    genes       = tmp_genes
    gene_strand = tmp_gene_strand

    # Write the backbone sequences into a fasta file
    backbone_file  = open(base_fullpath_name + "_backbone.fa", 'w')        
    # variants w.r.t the backbone sequences into a SNP file
    var_file       = open(base_fullpath_name + ".snp", 'w')
    var_index_file = open(base_fullpath_name + ".index.snp", 'w')
    # variant frequence
    var_freq_file  = open(base_fullpath_name + ".snp.freq", 'w')
    # haplotypes
    haplotype_file = open(base_fullpath_name + ".haplotype", 'w')
    # pairs of a variant and the corresponding HLA allels into a LINK file    
    link_file      = open(base_fullpath_name + ".link", 'w')
    # Write all the sequences with dots removed into a file
    input_file     = open(base_fullpath_name + "_sequences.fa", 'w')
    # Write allele names into a file
    allele_file    = open("%s.allele" % base_fullpath_name, 'w')
    # Read partial alleles from hla.data, and write them into a file
    partial_file   = open("%s.partial" % base_fullpath_name, 'w')
    # Write HISATgenotype Database version into file
    with open("%s.version" % base_fullpath_name, 'w') as fnversion:
        fnversion.write('Database %s derived from HISATgenotype DB version: %s' \
                            % (base_fullpath_name, dbversion))

    num_vars       = 0
    num_haplotypes = 0
    full_alleles   = {}
    for gene in gene_names:
        if gene not in genes:
            continue # Break if the gene was filtered above

        ref_gene      = genes[gene]
        strand        = gene_strand[gene]
        left_ext_seq  = ""
        right_ext_seq = ""
        if gene in left_ext_seq_dic:
            left_ext_seq  = left_ext_seq_dic[gene]
            right_ext_seq = right_ext_seq_dic[gene]

        if base_fname in spliced_gene:
            MSA_fname = "hisatgenotype_db/%s/msf/%s_gen.msf" \
                            % (base_fname.upper(), gene)
        else:
            MSA_fname = "hisatgenotype_db/%s/msf/%s_gen.msf" \
                            % (base_fname.upper(), gene)
            
        if not os.path.exists(MSA_fname):
            print("Warning: %s does not exist" % MSA_fname, 
                  file=sys.stderr)
            continue

        names, seqs = typing_common.read_MSF_file(MSA_fname,
                                                  full_alleles, 
                                                  left_ext_seq, 
                                                  right_ext_seq)
        full_allele_names   = set(names.keys())

        # Identify a consensus sequence
        assert len(seqs) > 0

        seq_len       = find_seq_len(seqs)        
        backbone_name = "%s*BACKBONE" % gene
        backbone_seq, backbone_freq = create_consensus_seq(seqs,
                                                           seq_len,
                                                           min_var_freq,
                                                           not partial) # Remove empty sequences?
        
        # Allele sequences can shrink, so readjust the sequence length
        if not partial:
            seq_len = find_seq_len(seqs)

        if partial and base_fname in spliced_gene:
            partial_MSA_fname = "hisatgenotype_db/%s/msf/%s_nuc.msf" \
                                    % (base_fname.upper(), gene)
            if not os.path.exists(partial_MSA_fname):
                print("Warning: %s does not exist" % partial_MSA_fname, 
                      file=sys.stderr)
                continue
            partial_names, partial_seqs = typing_common.read_MSF_file(partial_MSA_fname,
                                                                      full_alleles)
            
            ref_seq             = seqs[names[ref_gene]]
            ref_seq_map         = create_map(ref_seq)
            ref_partial_seq     = partial_seqs[partial_names[ref_gene]]
            ref_partial_seq_map = create_map(ref_partial_seq)
            exons               = gene_exons[gene]
            exon_len            = 0
            ref_exons           = [] # converted exons to MSF file (e.g. A_gen.msf)
            ref_partial_exons   = [] # converted exons to MSF file (e.g. A_nuc.msf)

            complete = True
            for exon in exons:
                left, right = exon
                ref_exons.append([ref_seq_map[left], ref_seq_map[right]])
                next_exon_len = right - left + exon_len
                if next_exon_len >= len(ref_partial_seq_map):
                    print("Warning: partial sequences (%s) seem to be incomplete" \
                            % gene, 
                          file=sys.stderr)
                    complete = False
                    seq_len  = find_seq_len(seqs)
                    break
                ref_partial_exons.append([ref_partial_seq_map[exon_len], 
                                          ref_partial_seq_map[next_exon_len]])
                exon_len += (right - left + 1)
                # Make sure two MSF files (e.g. A_gen.msf and A_nuc.msf) 
                # share the same MSF lengths in the exonic sequences
                ref_exon_len         = ref_exons[-1][1] \
                                        - ref_exons[-1][0] \
                                        + 1
                ref_partial_exon_len = ref_partial_exons[-1][1] \
                                        - ref_partial_exons[-1][0] \
                                        + 1
                assert ref_exon_len == ref_partial_exon_len

            if complete:
                partial_seq_len = find_seq_len(partial_seqs)
                partial_backbone_seq, \
                    partial_backbone_freq = create_consensus_seq(partial_seqs,
                                                                 partial_seq_len,
                                                                 min_var_freq,
                                                                 False) # Remove
                                                                        # empty 
                                                                        # sequences?
                for name in partial_names:
                    if name in names:
                        continue
                    seq_id  = partial_names[name]
                    seq     = partial_seqs[seq_id]
                    new_seq = ""
                    right   = 0
                    for e in range(len(exons)):
                        ref_exon         = ref_exons[e]
                        ref_partial_exon = ref_partial_exons[e]
                        new_seq         += backbone_seq[right:ref_exon[0]]
                        exon_seq         = seq[ref_partial_exon[0]\
                                               :ref_partial_exon[1] + 1]
                        nt_exon_seq      = exon_seq\
                                                .replace('.', '')\
                                                .replace('~', '')
                        if len(nt_exon_seq) == 0:
                            exon_seq = partial_backbone_seq[ref_partial_exon[0]\
                                                            :ref_partial_exon[1] \
                                                                + 1]
                        new_seq += exon_seq
                        right    = ref_exon[1] + 1
                    new_seq    += backbone_seq[right:]
                    names[name] = len(seqs)
                    seqs.append(new_seq)

                backbone_seq, backbone_freq = create_consensus_seq(seqs,
                                                                   seq_len,
                                                                   min_var_freq,
                                                                   remove_empty=True)
                seq_len = find_seq_len(seqs)
                
        # CB - Experimental code: Fill in internal ~ 
        # with consensus background sequence
        # Required if there are ~ regions in the center of alleles
        missing_seq = False
        for itr in range(len(seqs)):
            if "~" not in seqs[itr]:
                continue

            missing_seq = True
            seq = seqs[itr]    
            assert len(seq) == len(backbone_seq)
            seq_ = ''
            for s in range(len(seq)):
                if seq[s] == "~":
                    seq_ += backbone_seq[s]
                else:
                    assert seq[s] in "ACGT.N"
                    seq_ += seq[s]
            seqs[itr] = seq_
        
        if missing_seq:
            print("Warning %s contains missing sequence "\
                    "in the data. Filling in with consensus" % gene)

        names, seqs, collapsed = typing_common.collapse_alleles(names, 
                                                                seqs, 
                                                                list_collapse=True, 
                                                                verbose=True) 

        if ref_gene in collapsed:
            ref_gene    = collapsed[ref_gene]
            genes[gene] = ref_gene

        # Check for empty sequences or omitted nucleotides. 
        # This can be empty if there is a problem with the database.
        # Solution is to omit the offending gene from the 
        # database to prevent premature termination of database build.
        if min_var_freq <= 0.0:
            omits = ['.', 'E', '~']
            breakout = False

            for omit in omits:
                if omit in backbone_seq:
                    if verbose:
                        print('%s in backbone of %s with no "\
                                  "minimum variation set' % (omit, gene))
                    print('Error in database: Omitting %s!!' % (gene),
                          file=sys.stderr)
                    breakout = True
                    break

            if breakout:
                continue

        # Reverse complement MSF if this gene is on '-' strand
        if strand == '-':
            # Reverse exons
            ref_seq     = seqs[names[ref_gene]].replace('.', '').replace('~', '')
            ref_seq_len = len(ref_seq)
            if base_fname in spliced_gene:
                exons = []
                for left, right in reversed(gene_exons[gene]):
                    rleft  = ref_seq_len - right - 1
                    rright = ref_seq_len - left - 1
                    exons.append([rleft, rright])
                gene_exons[gene] = exons
                exon_counts      = {}
                for exon_i, count in gene_exon_counts[gene].items():
                    exon_counts[len(gene_exons[gene]) - exon_i - 1] = count
                gene_exon_counts[gene] = exon_counts

            for i in range(len(seqs)):
                seqs[i] = typing_common.reverse_complement(seqs[i])
            backbone_seq, backbone_freq = create_consensus_seq(seqs, 
                                                               seq_len, 
                                                               min_var_freq, 
                                                               remove_empty = True)
            seq_len = find_seq_len(seqs)

        if leftshift:
            for seq_i in range(len(seqs)):
                seqs[seq_i] = leftshift_deletions(backbone_seq, seqs[seq_i])
            backbone_seq, backbone_freq = create_consensus_seq(seqs, 
                                                               seq_len, 
                                                               min_var_freq, 
                                                               remove_empty = True)
            seq_len = find_seq_len(seqs)

        print("%s: number of alleles is %d." % (gene, len(names)), 
              file=sys.stderr)

        Vars = {}
        for cmp_name, id in names.items():
            if cmp_name == backbone_name:
                continue
            assert id < len(seqs)
            cmp_seq = seqs[id]
            if len(cmp_seq) != seq_len:
                print("Warning: the length of %s (%d) is different from %d" \
                        % (cmp_name, len(cmp_seq), seq_len), 
                      file=sys.stderr)
                continue

            # Insert Variants into Var dictionary variable
            # Inherits 4 variables from above and shouldn't be moved
            def insertVar(type, info):
                pos, backbone_pos, data = info
                if type in "MI":
                    varKey = "%d-%s-%s" % (pos, type, data)
                else:
                    varKey = "%d-%s-%d" % (pos, type, data)

                if varKey not in Vars:
                    if type == 'M':
                        assert data in backbone_freq[backbone_pos], \
                            "Data %s not in backbone %s of type %s" % (
                                data, 
                                backbone_freq[backbone_pos], 
                                type
                            )
                        freq = backbone_freq[backbone_pos][data]
                    elif type == 'D':
                        del_len = int(data)
                        freq    = 100.0
                        for d in range(del_len):
                            assert '.' in backbone_freq[backbone_pos + d]
                            freq2 = backbone_freq[backbone_pos + d]['.']
                            if freq2 < freq:
                                freq = freq2
                    else:
                        assert type == 'I'
                        ins_len = len(data)
                        freq    = 100.0
                        for i in range(ins_len):
                            nt = data[i]
                            assert nt in backbone_freq[backbone_pos + i]
                            freq2 = backbone_freq[backbone_pos + i][nt]
                            if freq2 < freq:
                                freq = freq2
                        assert freq <= min_var_freq
                    
                    Vars[varKey] = [freq, [cmp_name]]
                else:
                    Vars[varKey][1].append(cmp_name)

            insertion = []
            deletion  = []
            ndots     = 0
            for s in range(seq_len):
                assert not (insertion and deletion)
                bc = backbone_seq[s]
                cc = cmp_seq[s]
                if bc not in '.~' and cc not in '.~':
                    if insertion:
                        insertVar('I', insertion)
                        insertion = []
                    elif deletion:
                        insertVar('D', deletion)
                        deletion = []
                    if bc != cc:
                        mismatch = [s - ndots, s, cc]
                        insertVar('M', mismatch)
                elif bc == '.' and cc not in '.~':
                    if deletion:
                        insertVar('D', deletion)
                        deletion = []
                    if insertion:
                        insertion[2] += cc
                    else:
                        insertion = [s - ndots, s, cc]
                elif bc not in '.~' and cc == '.':
                    if insertion:
                        insertVar('I', insertion)
                        insertion = []
                    if deletion:
                        deletion[2] += 1
                    else:
                        deletion = [s - ndots, s, 1]

                if bc == '.':
                    ndots += 1

            if insertion:
                insertVar('I', insertion)
            elif deletion:
                insertVar('D', deletion)
        print("Number of variants is %d." % (len(Vars.keys())), 
              file=sys.stderr)       

        Vars_ = {}
        for key, values in Vars.items():
            freq, names_ = values
            for name in names_:
                if not name in Vars_:
                    Vars_[name] = [key]
                else:
                    Vars_[name].append(key)
        for name, vars in Vars_.items():
            Vars_[name] = sorted(vars, key=key_varKey)

            if SANITY_CHECK: # SANITY CHECK for sorting
                validation_check.validate_variants(Vars_[name])

        if SANITY_CHECK: # SANITY CHECK for rebuilding seqs
            validation_check.validate_constructs(names,
                                                 backbone_name,
                                                 backbone_seq,
                                                 Vars_,
                                                 seqs)

        # Remap the backbone allele, which is sometimes slighly different from
        #   fasta version
        ref_backbone_id  = names[ref_gene]
        ref_backbone_seq = seqs[ref_backbone_id]
        aligner_cmd = ["hisat2"]
        if base_fname in spliced_gene:
            # Guarantee two missmatchs then one per 5000bp
            aligner_cmd += ["--score-min", "L,-12,-0.0012"] 
        aligner_cmd += ["--no-unal",
                        "-x", "grch38/genome",
                        "-f", 
                        "-c", "%s" % ref_backbone_seq.replace('.', '')]
        align_proc = subprocess.Popen(aligner_cmd,
                                      universal_newlines = True,
                                      stdout = subprocess.PIPE,
                                      stderr = open("/dev/null", 'w'))
        best_chr   = ""
        best_left  = 0
        best_right = 0
        best_AS    = -sys.maxsize
        for line in align_proc.stdout:
            if line.startswith('@'):
                continue
            line = line.strip()
            cols = line.split()
            allele_id, flag, chr, left, mapQ, cigar_str = cols[:6]
            flag  = int(flag)
            left  = int(left) - 1
            right = left
            AS    = ""

            assert flag & 0x10 == 0, 'Allele %s with flag %d' % (line, flag)
            for i in range(11, len(cols)):
                col = cols[i]
                if col.startswith("AS"):
                    AS = col[5:]
            AS     = int(AS)
            cigars = cigar_re.findall(cigar_str)
            cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
            for i in range(len(cigars)):
                cigar_op, length = cigars[i]
                if cigar_op in "MND":
                    right += length
            if AS > best_AS:
                best_chr   = chr
                best_left  = left
                best_right = right
                best_AS    = AS
        chr   = best_chr
        left  = best_left
        right = best_right
        align_proc.communicate()
        if left == right:
            print("Warning: %s (%s) is not remapped: Removing" \
                    % (gene, ref_gene), 
                  file=sys.stderr)
            continue # TODO: This is causing no var printout in .snp file
        assert left < right

        # Write the backbone sequences into a fasta file
        print(">%s" % (backbone_name), 
              file=backbone_file)
        backbone_seq_ = backbone_seq.replace('.', '')
        assert "~" not in backbone_seq_
        for s in range(0, len(backbone_seq_), 60):
            print(backbone_seq_[s:s+60], 
                  file=backbone_file)

        base_locus  = 0                
        ref_seq     = seqs[names[ref_gene]]
        ref_seq_map = create_map(ref_seq)
        del_count   = []
        for nt in backbone_seq:
            assert nt in "ACGT."
            add = 1 if nt == '.' else 0
            if len(del_count) == 0:
                del_count.append(add)
            else:
                del_count.append(del_count[-1] + add)
        
        if base_fname in spliced_gene:
            exon_str = ""
            for exon_i in range(len(gene_exons[gene])):
                exon_left, exon_right = gene_exons[gene][exon_i]
                exon_left   = ref_seq_map[exon_left]
                exon_right  = ref_seq_map[exon_right]
                exon_left  -= del_count[exon_left]
                exon_right -= del_count[exon_right]
                if exon_str != "":
                    exon_str += ','
                primary = gene_exon_counts[gene][exon_i] \
                            == max(gene_exon_counts[gene].values())
                exon_str += ("%d-%d%s" 
                                % (exon_left, exon_right, 'p' if primary else ''))

            # Sanity check for exonic sequence
            if SANITY_CHECK \
                    and os.path.exists("hisatgenotype_db/%s/fasta/%s_nuc.fasta" \
                                            % (base_fname.upper(), gene)):
                validation_check.validate_exons(exon_str,
                                                backbone_seq,
                                                Vars_,
                                                ref_gene,
                                                ref_seq,
                                                gene_strand,
                                                gene,
                                                base_fname)
        else:
            exon_str = "%d-%d" % (left - left + 1, right - left + 1)

        print("%s\t%s\t%d\t%d\t%d\t%s\t%s" \
                % (backbone_name, 
                   chr, 
                   left, 
                   right - 1, 
                   len(backbone_seq.replace('.', '')), 
                   exon_str, 
                   gene_strand[gene]), 
              file=locus_file)

        # Write
        #       (1) variants w.r.t the backbone sequences into a SNP file
        #       (2) pairs of a variant and the corresponding HLA allels 
        #               into a LINK file    
        keys = sorted(Vars.keys(), key=key_varKey)
        if SANITY_CHECK:
            validation_check.validate_variants(keys)
        var2ID = {}
        for k in range(len(keys)):
            locus, type, data = keys[k].split('-')
            locus = int(locus)
            if type == 'M':
                type_str = "single"
            elif type == 'I':
                type_str = "insertion"
            else:
                assert type == 'D'
                type_str = "deletion"

            freq, names_      = Vars[keys[k]]
            names_            = sorted(names_)            
            varID             = "hv%d" % (num_vars)
            tmp_backbone_name = backbone_name
            print("%s\t%s\t%s\t%d\t%s" \
                    % (varID, 
                       type_str, 
                       tmp_backbone_name, 
                       base_locus + locus, 
                       data),
                  file=var_file)
            if freq >= min_var_freq:
                print("%s\t%s\t%s\t%d\t%s" \
                        % (varID, 
                           type_str, 
                           tmp_backbone_name, 
                           base_locus + locus, 
                           data), 
                      file=var_index_file)
            print("%s\t%.2f" % (varID, freq), 
                  file=var_freq_file)
            print("%s\t%s" % (varID, ' '.join(names_)), 
                  file=link_file)
            var2ID[keys[k]] = num_vars
            num_vars += 1

        # Write haplotypes
        add_seq_len   = 0
        excluded_vars = set()
        var_leftmost  = sys.maxsize
        var_rightmost = -1
        for k in range(len(keys)):
            key = keys[k]
            if Vars[key][0] < min_var_freq:
                excluded_vars.add(key)

            # Update leftmost and rightmost of Vars
            locus, type, data = key.split('-')
            left = right = int(locus)
            if type == 'D':
                right = left + int(data) - 1
            if k == 0:
                var_leftmost = left
            if var_rightmost < right:
                var_rightmost = right

        i = 0
        while i < len(keys):
            key_i = keys[i]
            locus, type, data = key_i.split('-')
            locus = int(locus)
            if type == 'D':
                locus += (int(data) - 1)
            prev_locus = locus
            if whole_haplotype:
                j = len(keys)
            else:
                j = i + 1
                while j < len(keys):
                    key_j = keys[j]
                    locus2, type2, data2 = key_j.split('-')
                    locus2 = int(locus2)
                    if prev_locus + inter_gap < locus2:
                        break
                    prev_locus = locus2
                    if type == 'D':
                        prev_locus += (int(data) - 1)
                    j += 1

            alleles = set()
            for k in range(i, j):
                key_k = keys[k]
                freq, names_ = Vars[key_k]
                if freq < min_var_freq:
                    continue
                add_alleles = set(names_)
                alleles    |= add_alleles

            haplotypes = set()
            cur_vars   = set(keys[i:j]) - excluded_vars
            for allele in alleles:
                allele_vars = set(Vars_[allele]) - excluded_vars
                allele_var_sort = sorted(list(cur_vars & allele_vars), 
                                         key=key_varKey)
                if SANITY_CHECK:
                    validation_check.validate_variants(allele_var_sort)

                allele_cur_vars = '#'.join(allele_var_sort)
                haplotypes.add(allele_cur_vars)

            if not whole_haplotype:
                haplotypes = split_haplotypes(haplotypes, intra_gap)

            haplotypes = sorted(list(haplotypes), key = hapKey)
            if SANITY_CHECK: # SANITY CHECK for sorting
                validation_check.validate_haplotype(haplotypes)

            # Write haplotypes
            sanity_vars = set()
            for h_i in range(len(haplotypes)):
                h = haplotypes[h_i].split('#')
                varIDs = []
                for var in h:
                    varIDs.append("hv%s" % var2ID[var])
                    sanity_vars.add(var2ID[var])
                if whole_haplotype:
                    h_begin = var_leftmost
                    h_end   = var_rightmost
                else:
                    h1_locus, _, _ = h[0].split('-')
                    h2_locus, h2_type, h2_data = h[-1].split('-')
                    h_begin = int(h1_locus)
                    h_end   = int(h2_locus)
                    if h2_type == 'D':
                        h_end += (int(h2_data) - 1)
                    assert h_begin <= h_end
                    h_new_begin = h_begin
                    for h_j in reversed(range(0, h_i)):
                        hc = haplotypes[h_j].split('#')
                        hc_begin, hc_type, hc_data = hc[-1].split('-')
                        hc_begin = int(hc_begin)
                        hc_end   = hc_begin
                        if hc_type == 'D':
                            hc_end += (int(hc_data) - 1)
                        if hc_end + inter_gap < h_begin:
                            break
                        if h_new_begin > hc_end:
                            h_new_begin = hc_end
                    assert h_new_begin <= h_begin
                    h_begin = h_new_begin
                tmp_backbone_name = backbone_name
                print("ht%d\t%s\t%d\t%d\t%s" \
                        % (num_haplotypes, 
                           tmp_backbone_name, 
                           base_locus + h_begin, base_locus + h_end, 
                           ','.join(varIDs)),
                      file=haplotype_file)
                num_haplotypes += 1
                add_seq_len    += (h_end - h_begin + 1)
            assert len(sanity_vars) == len(cur_vars)
            i = j

        print(("Length of additional sequences for haplotypes:", add_seq_len), 
              file=sys.stderr)

        sorted_all_alleles = typing_common.sort_genall(list(names.keys()), 
                                                       alleles = True)

        # Write all the sequences with dots removed into a file
        for name in sorted_all_alleles:
            if name not in names:
                continue
            ID = names[name]
            print(">%s" % (name), file=input_file)
            assert ID < len(seqs)
            seq = seqs[ID].replace('.', '').replace('~', '')
            for s in range(0, len(seq), 60):
                print(seq[s:s+60], file=input_file)
            print(name, file=allele_file)
            if name not in full_allele_names:
                print(name, file=partial_file)

    backbone_file.close()
    locus_file.close()
    var_file.close()
    var_index_file.close()
    var_freq_file.close()
    haplotype_file.close()
    link_file.close()
    input_file.close()
    allele_file.close()
    partial_file.close()
   
# --------------------------------------------------------------------------- #
# Scripts for extract_reads                                                   #
#                                                                             #
# Given input fastq(a) extract reads from a given location chosen by user     #
# --------------------------------------------------------------------------- #
"""
Process for parallelizing extractions
TODO: consider shifting to subprocess
"""
def parallel_work(pids, 
                  work, 
                  fq_fname_base, 
                  fq_fname, 
                  fq_fname2, 
                  ranges,
                  simulation,
                  verbose):
    child = -1
    for i in range(len(pids)):
        if pids[i] == 0:
            child = i
            break

    while child == -1:
        status = os.waitpid(0, 0)
        for i in range(len(pids)):
            if status[0] == pids[i]:
                child = i
                pids[i] = 0
                break

    child_id = os.fork()
    if child_id == 0:
        work(fq_fname_base, 
             fq_fname, 
             fq_fname2, 
             ranges,
             simulation,
             verbose)
        os._exit(os.EX_OK)
    else:
        pids[child] = child_id

"""
Waiting for all of the processes to end
"""
def wait_pids(pids):
    for pid in pids:
        if pid > 0:
            os.waitpid(pid, 0)   

""" Write reads to fasta or fastq file """
def write_read(gzip_proc, 
                read_name, 
                seq, 
                qual,
                prev_read_name,
                fastq):
    if fastq:
        gzip_proc.stdin.write("@%s\n" % read_name)
        gzip_proc.stdin.write("%s\n" % seq)
        gzip_proc.stdin.write("+\n")
        gzip_proc.stdin.write("%s\n" % qual)
    else:
        gzip_proc.stdin.write(">%s\n" % prev_read_name)
        gzip_proc.stdin.write("%s\n" % seq) 

"""
Extract reads from given fastq(a) file

General Input: HISAT-genotype genome, base information, fastq(a) files
Output: Reads aligning to base region from genotype genome
"""
def extract_reads(base_fname,    # Base file name of genome to use
                  ix_dir,
                  database_list,
                  read_dir,
                  out_dir,
                  suffix,
                  read_fname,
                  fastq,
                  paired,
                  simulation,
                  threads,
                  threads_aprocess,
                  max_sample,
                  job_range,
                  aligner,
                  block_size,
                  verbose):
    if block_size > 0:
        resource.setrlimit(resource.RLIMIT_NOFILE, (1000, 1000))
        resource.setrlimit(resource.RLIMIT_NPROC, (1000, 1000))

    if not typing_common.check_base(base_fname, aligner, ix_dir):
        if base_fname == "genotype_genome":
            print("Downloading genotype_genome into %s" % ix_dir)
            typing_common.download_genotype_genome(ix_dir)
        else:
            exit(1)

    fname_list    = {} # For use in Hisatgenotype script
    filter_region = len(database_list) > 0
    ranges        = []
    regions       = {}
    region_loci   = {}
    for line in open("%s.locus" % base_fname):
        family, allele_name, chr, left, right = line.strip().split()[:5]
        if filter_region and family.lower() not in database_list:
            continue
        region_name = "%s-%s" % (family, allele_name.split('*')[0])
        assert region_name not in regions
        regions[region_name] = allele_name
        left  = int(left)
        right = int(right)

        if chr not in region_loci:
            region_loci[chr] = {}
        region_loci[chr][region_name] = [allele_name, chr, left, right]
        if family.lower() not in database_list:
            database_list.append(family.lower())

    if out_dir != "":
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        out_dir = out_dir if out_dir.endswith("/") else out_dir + "/"    
    elif out_dir == "":
        out_dir = "./"

    # Extract reads
    if len(read_fname) > 0:
        if paired:
            fq_fnames  = [read_fname[0]]
            fq_fnames2 = [read_fname[1]]
        else:
            fq_fnames = read_fname
    else:
        fq_fnames = glob.glob("%s/*.%s" % (read_dir, suffix)) 
        fq_fnames = sorted(fq_fnames)
        if paired:
            fq_fnames, fq_fnames2, paired_fq_basen \
                = typing_common.get_filename_match(fq_fnames)
        
        if len(fq_fnames) == 0:
            print("Error: no files identified in %s directory with suffix .%s" \
                    % (read_dir, suffix), 
                  file=sys.stderr)
            exit(1)

    count = 0
    pids = [0 for i in range(threads)]
    for file_i in range(len(fq_fnames)):
        if file_i >= max_sample:
            break
        fq_fname = fq_fnames[file_i]
        if job_range[1] > 1:
            if job_range[0] != (file_i % job_range[1]):
                continue

        if paired:
            if len(read_fname) > 0:
                fq_fname_base = fq_fname.split('/')[-1]
            else:
                fq_fname_base = paired_fq_basen[file_i].split('/')[-1]
        else:
            fq_fname_base = fq_fname.split('/')[-1].split('.')[0]
            
        if paired:
            fq_fname2 = fq_fnames2[file_i]
            if not os.path.exists(fq_fname2):
                print("%s does not exist." % (fq_fname2),
                      file=sys.stderr)
                continue
        else:
            fq_fname2 = ""

        # Process database and if files exist do not reextract
        omit_extract = True
        for database in database_list:
            if database not in fname_list:
                fname_list.update({ database : [] })
            fname_list[database].append('%s-%s' % (fq_fname_base, database))

            if paired \
                    and os.path.exists("%s%s-%s-extracted-1.fq.gz" \
                                           % (out_dir, fq_fname_base, database)):
                continue
            elif os.path.exists("%s%s-%s-extracted.fq.gz" \
                                    % (out_dir, fq_fname_base, database)):
                continue
            omit_extract = False

        if omit_extract:
            print("\tFiles found: Omitted extracting reads from %s" \
                      % (fq_fname_base), 
                  file=sys.stderr)
            continue

        count += 1        
        print("\t%d: Extracting reads from %s" \
                  % (count, fq_fname_base), 
              file=sys.stderr)
        
        def work(fq_fname_base,
                 fq_fname, 
                 fq_fname2, 
                 ranges,
                 simulation,
                 verbose):
            aligner_cmd = [aligner]
            if threads_aprocess > 1:
                aligner_cmd += ["-p", "%d" % threads_aprocess]
            if not fastq:
                aligner_cmd += ["-f"]
            aligner_cmd += ["-x", base_fname]
            if aligner == "hisat2":
                aligner_cmd += ["--no-spliced-alignment"]
                # aligner_cmd += ["--max-altstried", "64"]
            aligner_cmd += ["-X", "1000"]
            if paired:
                aligner_cmd += ["-1", fq_fname,
                                "-2", fq_fname2]
            else:
                aligner_cmd += ["-U", fq_fname]
            if verbose:
                print("\t\trunning", ' '.join(aligner_cmd), 
                      file=sys.stderr)

            align_proc = subprocess.Popen(aligner_cmd,
                                          universal_newlines = True,
                                          stdout = subprocess.PIPE,
                                          stderr = open("/dev/null", 'w'))

            gzip_dic = {}
            for database in database_list:
                if paired:
                    # LP6005041-DNA_A01-extracted-1.fq.gz
                    gzip1_proc = subprocess.Popen(
                        ["gzip"],
                        universal_newlines = True,
                        stdin=subprocess.PIPE,
                        stdout=open("%s%s-%s-extracted-1.fq.gz" \
                                        % (out_dir, 
                                           fq_fname_base, 
                                           database),
                                    'w'),
                        stderr=open("/dev/null", 'w')
                    )

                    # LP6005041-DNA_A01-extracted-2.fq.gz
                    gzip2_proc = subprocess.Popen(
                        ["gzip"],
                        universal_newlines = True,
                        stdin=subprocess.PIPE,
                        stdout=open("%s%s-%s-extracted-2.fq.gz" \
                                        % (out_dir, 
                                           fq_fname_base, 
                                           database),
                                    'w'),
                        stderr=open("/dev/null", 'w')
                    )

                else:
                    # LP6005041-DNA_A01-extracted-fq.gz
                    gzip1_proc = subprocess.Popen(
                        ["gzip"],
                        universal_newlines = True,
                        stdin=subprocess.PIPE,
                        stdout=open("%s%s-%s-extracted.fq.gz" \
                                        % (out_dir, 
                                           fq_fname_base, 
                                           database), 
                                    'w'),
                        stderr=open("/dev/null", 'w'))
                gzip_dic[database] = [gzip1_proc, gzip2_proc if paired else None]

            whole_gzip_dic = {}
            if block_size > 0:
                mult = block_size / 1000000
                for chr_line in open("%s.fa.fai" % base_fname):
                    chr, length = chr_line.strip().split('\t')[:2]
                    length = int(length)
                    if chr not in [str(i+1) for i in range(22)] \
                                      + ['X', 'Y', 'MT']:
                        continue
                    length = (length + block_size - 1) / block_size
                    assert chr not in whole_gzip_dic
                    whole_gzip_dic[chr] = []
                    for region_i in range(length):
                        if paired:
                            # LP6005041-DNA_A01.extracted.1.fq.gz
                            gzip1_proc = subprocess.Popen(
                                ["gzip"],
                                universal_newlines = True,
                                stdin=subprocess.PIPE,
                                stdout=open("%s%s-%s-%d_%dM-extracted-1.fq.gz" 
                                                % (out_dir, 
                                                   fq_fname_base, 
                                                   chr, 
                                                   region_i * mult, 
                                                   (region_i + 1) * mult), 
                                            'w'),
                                stderr=open("/dev/null", 'w')
                            )

                            # LP6005041-DNA_A01.extracted.2.fq.gz
                            gzip2_proc = subprocess.Popen(
                                ["gzip"],
                                universal_newlines = True,
                                stdin=subprocess.PIPE,
                                stdout=open("%s%s-%s-%d_%dM-extracted-2.fq.gz" 
                                                % (out_dir, 
                                                   fq_fname_base, 
                                                   chr, 
                                                   region_i * mult, 
                                                   (region_i + 1) * mult), 
                                            'w'),
                                stderr=open("/dev/null", 'w')
                            )
                        else:
                            # LP6005041-DNA_A01.extracted.fq.gz
                            gzip1_proc = subprocess.Popen(
                                ["gzip"],
                                universal_newlines = True,
                                stdin=subprocess.PIPE,
                                stdout=open("%s%s-%s-%d_%dM-extracted.fq.gz" 
                                                % (out_dir, 
                                                   fq_fname_base, 
                                                   chr, 
                                                   region_i * mult, 
                                                   (region_i + 1) * mult), 
                                            'w'),
                                stderr=open("/dev/null", 'w')
                            )
                        whole_gzip_dic[chr].append([gzip1_proc, 
                                                    gzip2_proc if paired \
                                                                else None])                 

            prev_read_name     = ""
            extract_read       = set()
            whole_extract_read = set()
            read1              = []
            read2              = []
            read1_first        = True
            read2_first        = True
            chk_line           = True
            for line in align_proc.stdout:
                if line.startswith('@'):
                    continue
                line = line.strip()
                cols = line.split()
                read_name, flag, chr, pos, mapQ, cigar, _, _, _, read, qual = cols[:11]
                
                flag = int(flag)
                pos  = int(pos) - 1                 
                AS   = ""
                XS   = ""
                NH   = ""
                for i in range(11, len(cols)):
                    col = cols[i]
                    if col.startswith("AS"):
                        AS = int(col[5:])
                    elif col.startswith("XS"):
                        XS = int(col[5:])
                    elif col.startswith("NH"):
                        NH = int(col[5:])

                if (chk_line and prev_read_name != ""):
                    chk_line = False
                    if (read_name != prev_read_name \
                            and not simulation \
                            and paired):
                        print("Error: Paired read names are not the same", 
                              file=sys.stderr)
                        sys.exit(1)

                if (not simulation and read_name != prev_read_name) \
                        or (simulation \
                              and read_name.split('|')[0] \
                                    != prev_read_name.split('|')[0]):
                    for region in extract_read:
                        write_read(gzip_dic[region][0], 
                                   prev_read_name, 
                                   read1[0], 
                                   read1[1],
                                   prev_read_name,
                                   fastq)
                        if paired:
                            write_read(gzip_dic[region][1], 
                                       prev_read_name, 
                                       read2[0], 
                                       read2[1],
                                       prev_read_name,
                                       fastq)
                            
                    for chr_region_num in whole_extract_read:
                        region_chr, region_num = chr_region_num.split('-')
                        region_num = int(region_num)
                        if region_chr not in whole_gzip_dic:
                            continue

                        assert region_num < len(whole_gzip_dic[region_chr])
                        write_read(whole_gzip_dic[region_chr][region_num][0], 
                                   prev_read_name, 
                                   read1[0], 
                                   read1[1],
                                   prev_read_name,
                                   fastq)
                        if paired:
                            write_read(whole_gzip_dic[region_chr][region_num][1], 
                                       prev_read_name, 
                                       read2[0], 
                                       read2[1],
                                       prev_read_name,
                                       fastq)

                    prev_read_name     = read_name
                    extract_read       = set()
                    whole_extract_read = set()
                    read1              = []
                    read2              = []
                    read1_first        = True
                    read2_first        = True

                if flag & 0x4 == 0 \
                        and ((aligner == "hisat2" and NH == 1) \
                                or (aligner == "bowtie2" \
                                        and AS > XS \
                                        and read1_first if flag & 0x40 \
                                        or not paired else read2_first)):
                    if chr in region_loci:
                        for region, loci in region_loci[chr].items():
                            region = region.split('-')[0].lower()
                            _, _, loci_left, loci_right = loci
                            # there might be a different candidate 
                            # region for each of left and right reads
                            if pos >= loci_left and pos < loci_right:
                                extract_read.add(region)
                                break
                    if block_size > 0:
                        chr_region_num = "%s-%d" % (chr, pos / block_size)
                        whole_extract_read.add(chr_region_num)

                if flag & 0x40 or not paired: # left read
                    read1_first = False
                    if not read1:
                        if flag & 0x10: # reverse complement
                            read1 = [typing_common.reverse_complement(read), 
                                     qual[::-1]]
                        else:
                            read1 = [read, qual]
                else:
                    assert flag & 0x80 # right read
                    read2_first = False
                    if flag & 0x10: # reverse complement
                        read2 = [typing_common.reverse_complement(read), 
                                 qual[::-1]]
                    else:
                        read2 = [read, qual]

            for region in extract_read:
                write_read(gzip_dic[region][0], 
                           prev_read_name, 
                           read1[0], 
                           read1[1],
                           prev_read_name,
                           fastq)
                if paired:
                    write_read(gzip_dic[region][1], 
                               prev_read_name, 
                               read2[0], 
                               read2[1],
                               prev_read_name,
                               fastq)

            for chr_region_num in whole_extract_read:
                region_chr, region_num = chr_region_num.split('-')
                region_num = int(region_num)
                if region_chr not in whole_gzip_dic:
                    continue
                assert region_num < len(whole_gzip_dic[region_chr])
                write_read(whole_gzip_dic[region_chr][region_num][0], 
                           prev_read_name, 
                           read1[0], 
                           read1[1],
                           prev_read_name,
                           fastq)
                if paired:
                    write_read(whole_gzip_dic[region_chr][region_num][1], 
                               prev_read_name, 
                               read2[0], 
                               read2[1],
                               prev_read_name,
                               fastq)

            for gzip1_proc, gzip2_proc in gzip_dic.values():
                gzip1_proc.stdin.close()
                if paired:
                    gzip2_proc.stdin.close()

            for gzip_list in whole_gzip_dic.values():
                for gzip1_proc, gzip2_proc in gzip_list:
                    gzip1_proc.stdin.close()
                    if paired:
                        gzip2_proc.stdin.close()         

        if threads <= 1:
            work(fq_fname_base, 
                 fq_fname, 
                 fq_fname2,
                 ranges,
                 simulation,
                 verbose)
        else:
            parallel_work(pids, 
                          work, 
                          fq_fname_base, 
                          fq_fname, 
                          fq_fname2, 
                          ranges,
                          simulation,
                          verbose)

    if threads > 1:
        wait_pids(pids)

    return fname_list
