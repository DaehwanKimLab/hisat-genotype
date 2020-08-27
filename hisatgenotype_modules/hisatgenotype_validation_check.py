#!/usr/bin/env python
# --------------------------------------------------------------------------- #
# Copyright 2020, Christopher Bennett <christopher@bennett-tech.dev>          #
#                                                                             #
# This file is part of HISAT-genotype. It contains basic unit tests and       #
# validation tests.                                                           #
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
import hisatgenotype_typing_common as typing_common

# --------------------------------------------------------------------------- #
# Sorting Unit tests                                                          #
# --------------------------------------------------------------------------- #
""" Validate the haplotype sorting in typing_process ~line 1300 """
def validate_haplotype(haplotypes):
    valid = True
    msg   = []
    if not haplotypes:
        valid = False
        msg.append("Empty data")

    def cmp_haplotype(a, b):
        a = a.split('#')
        a1_locus, _, _ = a[0].split('-')
        a2_locus, a2_type, a2_data = a[-1].split('-')
        a_begin = int(a1_locus)
        a_end   = int(a2_locus)
        if a2_type == 'D':
            a_end += (int(a2_data) - 1)
        b = b.split('#')
        b1_locus, _, _ = b[0].split('-')
        b2_locus, b2_type, b2_data = b[-1].split('-')
        b_begin = int(b1_locus)
        b_end   = int(b2_locus)
        if b2_type == 'D':
            b_end += (int(b2_data) - 1)
        if a_begin != b_begin:
            return a_begin - b_begin
        return a_end - b_end

    for i in range(len(haplotypes)-1):
        a = haplotypes[i]
        b = haplotypes[i+1]

        val = cmp_haplotype(a, b)
        if val > 0:
            valid = False
            msg.append(a + "_" + b)

    if not valid:
        print("Error: Failed Heplotype Sorting!!",
              file=sys.stderr)
        print("\n".join(msg),
              file=sys.stderr)                
        exit(1)

""" Validate the variant sorting in typing_process ~line 950 """
def validate_variants(vars):
    valid = True
    msg   = []
    if not vars:
        valid = False
        msg.append("Empty data")

    def cmp_varKey(a, b):
        a_locus, a_type, a_data = a.split('-')
        b_locus, b_type, b_data = b.split('-')
        a_locus = int(a_locus)
        b_locus = int(b_locus)
        if a_locus != b_locus:
            return a_locus - b_locus
        if a_type != b_type:
            if a_type == 'I':
                return -1
            elif b_type == 'I':
                return 1
            elif a_type == 'M':
                return -1
            else:
                assert b_type == 'M'
                return 1
        assert a_data != b_data
        if a_type in "MI":
            if a_data < b_data:
                return -1
            else:
                return 1
        else:
            assert a_type == 'D'
            return int(a_data) - int(b_data)

    for i in range(len(vars)-1):
        a = vars[i]
        b = vars[i+1]

        val = cmp_varKey(a, b)
        if val > 0:
            valid = False
            msg.append(a + "_" + b)

    if not valid:
        print("Error: Failed Variant Sorting!!",
              file=sys.stderr)
        print("\n".join(msg),
              file=sys.stderr)
        exit(1)

""" Validate node sorting in assembly_graph script line ~600 """
def validate_node_sorting(nodes):
    def node_cmp(a, b):
        if a[2] != b[2]:
            return a[2] - b[2]
        else:
            return a[1] - b[1]
    
    valid = True
    msg   = []
    for i in range(len(nodes)-1):
        a = nodes[i]
        b = nodes[i+1]
        val = node_cmp(a, b)

        if val > 0:
            valid = False
            msg.append("\t%s\n\t%s" % (str(a), str(b)))

    if not valid:
        print("Error: Failed Node Sorting!!",
              file=sys.stderr)
        print("\n".join(msg),
              file=sys.stderr)
        exit(1)

# --------------------------------------------------------------------------- #
# Data Validity Unit tests                                                    #
# --------------------------------------------------------------------------- #
""" Sanity check -
    #    (1) Reconstruct the other sequences from the backbone 
    #           sequence and variants and
    #    (2) Confirm these constructed sequences are the same 
    #           as those input sequences.
"""
def validate_constructs(names, 
                        backbone_name, 
                        backbone_seq, 
                        Vars_,
                        seqs):
    for cmp_name, id in names.items():
        if cmp_name == backbone_name:
            continue

        constr_seq = backbone_seq.replace('.', '')
        assert "~" not in constr_seq
        constr_seq = list(constr_seq)
        locus_diff = 0
        if cmp_name not in Vars_:
            continue
        
        for var in Vars_[cmp_name]:
            try:
                locus, type, data = var.split('-')
                locus = int(locus)
            except ValueError:
                continue

            if type == 'M':
                assert len(data) == 1
                constr_seq[locus + locus_diff] = data[0]
            elif type == 'I':
                assert locus + locus_diff >= 0
                assert locus + locus_diff <= len(constr_seq)
                constr_seq = constr_seq[:locus + locus_diff] \
                                + list(data) \
                                + constr_seq[locus + locus_diff:]
                locus_diff += len(data)
            else:
                assert type == 'D'
                assert locus + locus_diff + len(data) <= len(constr_seq)
                assert locus + locus_diff >= 0
                del_len = int(data)
                constr_seq = constr_seq[:locus + locus_diff] \
                                + constr_seq[locus + locus_diff + del_len:]
                locus_diff -= del_len

        assert id < len(seqs)
        cmp_seq = seqs[id].replace('.', '')
        if len(constr_seq) != len(cmp_seq):
            print("Error: reconstruction fails (%s)! \
                        Lengths different: %d vs. %d" \
                            % (cmp_name, len(constr_seq), len(cmp_seq)), 
                    file=sys.stderr)
            exit(1)

        # Add missing sequence markers
        if "~" in cmp_seq:
            for s in range(len(constr_seq)):
                if cmp_seq[s] == "~":
                    constr_seq[s] = "~"
        constr_seq = "".join(constr_seq)          

        # Sanity check
        for s in range(len(constr_seq)):
            if constr_seq[s] != cmp_seq[s]:
                print("Differ at %d: %s vs. %s (reconstruction vs. original)" \
                        % (s, constr_seq[s], cmp_seq[s]), 
                        file=sys.stderr)
                print("%s:%s vs. %s:%s" \
                        % (constr_seq[s-10:s], 
                            constr_seq[s:s+10], 
                            cmp_seq[s-10:s], 
                            cmp_seq[s:s+10]),
                        file=sys.stderr)

        if constr_seq != cmp_seq.replace('.', ''):
            print("Error: reconstruction fails for %s" % (cmp_name), 
                    file=sys.stderr)
            exit(1)

""" Sanity Check for exonic sequences """
def validate_exons(exon_str,
                   backbone_seq,
                   Vars_,
                   ref_gene,
                   ref_seq,
                   gene_strand,
                   gene,
                   base_fname,
                   hisatgenotype_db):
    exons_ = []
    for exon in exon_str.split(','):
        if exon.endswith('p'):
            exon = exon[:-1]
        exon_left, exon_right = exon.split('-')
        exon_left  = int(exon_left)
        exon_right = int(exon_right)
        exons_.append([exon_left, exon_right])

    backbone_seq_ = backbone_seq.replace('.', '')
    if ref_gene in Vars_:
        vars_ = Vars_[ref_gene]
    else:
        vars_ = []
    seq_ = list(backbone_seq_)
    has_insertion = False
    for var_ in vars_:
        var_pos, var_type, var_data = var_.split('-')
        var_pos = int(var_pos)
        assert var_pos >= 0 and var_pos < len(backbone_seq_)
        if var_type == 'M':
            seq_[var_pos] = var_data
        elif var_type == 'D':
            del_len = int(var_data)
            assert var_pos + del_len <= len(ref_seq)
            seq_[var_pos:var_pos + del_len] = ['.'] * del_len
        else:
            assert var_type == 'I'
            has_insertion = True

    seq_ = ''.join(seq_)
    exon_seq_ = ""
    for exon_left, exon_right in exons_:
        exon_seq_ += seq_[exon_left:exon_right+1]
    exon_seq_ = exon_seq_.replace('.', '').replace('~', '')
    if gene_strand[gene] == '-':
        exon_seq_ = typing_common.reverse_complement(exon_seq_)

    cmp_exon_seq_ = ""
    allele_name_  = ""
    for line in open("%s/%s/fasta/%s_nuc.fasta" % (hisatgenotype_db,
                                                   base_fname.upper(), 
                                                   gene)):
        if line.startswith(">"):
            if allele_name_ == ref_gene:
                break
            if base_fname == "hla":
                allele_name_ = line.strip().split()[1] 
            else:
                allele_name_ = line.strip().split()[0].replace('>','')
            cmp_exon_seq_ = ""
        else:
            cmp_exon_seq_ += line.strip()

    if exon_seq_ != cmp_exon_seq_:
        print("Warning: exonic sequences do not match (%s)" 
                % gene, 
               file=sys.stderr) 
        if len(exon_seq_) > len(cmp_exon_seq_) and cmp_exon_seq_ in exon_seq_:
            print("Warning: exonic sequences from file are smaller than rebuild",
                file=sys.stderr)
        elif len(cmp_exon_seq_) > len(exon_seq_) and exon_seq_ in cmp_exon_seq_:
            print("Warning: exonic sequences from file are larger than rebuild",
                file=sys.stderr)
        else:
            print("Warning: exonic sequences are not built properly",
                file=sys.stderr)
            exit(1)

""" Check uniqueness of sequences in typing_common line ~ 1850"""
def check_amb_uniqueness(cmp_list,
                         cmp_left,
                         cmp_right,
                         left_alt_set,
                         right_alt_set):
    ht_set_ = set()
    for ht in left_alt_set:
        ht = '-'.join(ht.split('-')[1:])
        if ht == "":
            continue
        if ht in ht_set_:
            print(("Error: %s should not be in" % ht, ht_set_), 
                  file=sys.stderr)

            # DK - debugging purposes
            print("DK: cmp_list_range: [%d, %d]" % (cmp_left, cmp_right))
            print("DK: cmp_list:", cmp_list)
            print("DK: left_alt_set:", left_alt_set, "right_alt_set:", right_alt_set)
            
            exit(1)
        ht_set_.add(ht)
    for ht in right_alt_set:
        ht = '-'.join(ht.split('-')[:-1])
        if ht == "":
            continue
        if ht in ht_set_:
            print(("Error: %s should not be in" % ht, ht_set_), file=sys.stderr)
            exit(1)
        ht_set_.add(ht)

""" Check if all exons alleles are in allele_rep_set on typing_core line ~550"""
def check_repset_inclusion(allele_rep_set,
                           allele_reps,
                           primary_exon_allele_reps):
    for exon_allele in primary_exon_allele_reps.keys():
        # DK - debugging purposes
        if exon_allele not in allele_rep_set:
            print("Error: %s not in Rep set!" % exon_allele,
                  file=sys.stderr)
            print((allele_reps[exon_allele], 
                   exon_allele in primary_exon_allele_reps.keys()))
            exit(1)

""" Check validity of genes in Genes variable in typing_core line ~2400 """
def check_allele_validity(base_fname, Genes): 
    Genes2 = {}
    typing_common.read_allele_seq(base_fname + "_backbone.fa", Genes2, True)
    typing_common.read_allele_seq(base_fname + "_sequences.fa", Genes2, True)
    for gene_name, alleles in Genes.items():
        assert gene_name in Genes2
        for allele_name, allele_seq in alleles.items():
            assert allele_name in Genes2[gene_name]
            allele_seq2 = Genes2[gene_name][allele_name]
            assert allele_seq == allele_seq2, \
                'Problems with: %s - length %d vs %d' \
                    % (allele_name, len(allele_seq), len(allele_seq2))
