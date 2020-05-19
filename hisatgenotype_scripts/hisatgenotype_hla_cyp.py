#!/usr/bin/env python

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#


import sys, os, subprocess, re
import inspect, random
import math
from argparse import ArgumentParser, FileType
import hisatgenotype_args as arguments

"""
"""
def simulate_reads(HLAs,
                   test_HLA_list,
                   simulate_interval):
    HLA_reads_1, HLA_reads_2 = [], []
    for test_HLA_names in test_HLA_list:
        gene = test_HLA_names[0].split('*')[0]
        # ref_allele = refHLAs[gene]
        # ref_seq = HLAs[gene][ref_allele]

        # Simulate reads from two HLA alleles
        def simulate_reads_impl(seq, simulate_interval = 1, frag_len = 250, read_len = 100):
            comp_table = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
            reads_1, reads_2 = [], []
            for i in range(0, len(seq) - frag_len + 1, simulate_interval):
                reads_1.append(seq[i:i+read_len])
                tmp_read_2 = reversed(seq[i+frag_len-read_len:i+frag_len])
                read_2 = ""
                for s in tmp_read_2:
                    if s in comp_table:
                        read_2 += comp_table[s]
                    else:
                        read_2 += s
                reads_2.append(read_2)
            return reads_1, reads_2

        for test_HLA_name in test_HLA_names:
            HLA_seq = HLAs[gene][test_HLA_name]
            tmp_reads_1, tmp_reads_2 = simulate_reads_impl(HLA_seq, simulate_interval)
            HLA_reads_1 += tmp_reads_1
            HLA_reads_2 += tmp_reads_2

    # Write reads into a fasta read file
    def write_reads(reads, idx):
        read_file = open('hla_input_%d.fa' % idx, 'w')
        for read_i in range(len(reads)):
            print(">%d" % (read_i + 1), file=read_file)
            print(reads[read_i], file=read_file)
        read_file.close()
    write_reads(HLA_reads_1, 1)
    write_reads(HLA_reads_2, 2)


"""
Align reads, and sort the alignments into a BAM file
"""
def align_reads(ex_path,
                base_fname,
                aligner,
                index_type,
                read_fname,
                fastq,
                threads,
                verbose):
    if aligner == "hisat2":
        hisat2 = os.path.join(ex_path, "hisat2")
        aligner_cmd = [hisat2,
                       "--no-unal",
                       "--mm"]
        if index_type == "linear":
            aligner_cmd += ["-k", "10"]
        aligner_cmd += ["-x", "%s.%s" % (base_fname, index_type)]
    elif aligner == "bowtie2":
        aligner_cmd = [aligner,
                       "--no-unal",
                       "-k", "10",
                       "-x", base_fname]
    else:
        assert False
    assert len(read_fname) in [1,2]
    aligner_cmd += ["-p", str(threads)]
    if not fastq:
        aligner_cmd += ["-f"]
    if len(read_fname) == 1:
        aligner_cmd += ["-U", read_fname[0]]
    else:
        aligner_cmd += ["-1", "%s" % read_fname[0],
                        "-2", "%s" % read_fname[1]]

    if verbose:
        print(' '.join(aligner_cmd), file=sys.stderr)
    align_proc = subprocess.Popen(aligner_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=open("/dev/null", 'w'))

    sambam_cmd = ["samtools",
                  "view",
                  "-bS",
                  "-"]
    sambam_proc = subprocess.Popen(sambam_cmd,
                                   stdin=align_proc.stdout,
                                   stdout=open("hla_input_unsorted.bam", 'w'),
                                   stderr=open("/dev/null", 'w'))
    sambam_proc.communicate()
    if index_type == "graph":
        bamsort_cmd = ["samtools",
                       "sort",
                       "hla_input_unsorted.bam",
                       "-o", "hla_input.bam"]
        bamsort_proc = subprocess.Popen(bamsort_cmd,
                                        stderr=open("/dev/null", 'w'))
        bamsort_proc.communicate()

        bamindex_cmd = ["samtools",
                        "index",
                        "hla_input.bam"]
        bamindex_proc = subprocess.Popen(bamindex_cmd,
                                         stderr=open("/dev/null", 'w'))
        bamindex_proc.communicate()

        os.system("rm hla_input_unsorted.bam")            
    else:
        os.system("mv hla_input_unsorted.bam hla_input.bam")


"""
""" 
def normalize(prob):
    total = sum(prob.values())
    for allele, mass in prob.items():
        prob[allele] = mass / total

        
"""
"""
def prob_diff(prob1, prob2):
    diff = 0.0
    for allele in prob1.keys():
        if allele in prob2:
            diff += abs(prob1[allele] - prob2[allele])
        else:
            diff += prob1[allele]
    return diff


"""
"""
def HLA_prob_cmp(a, b):
    if a[1] != b[1]:
        if a[1] < b[1]:
            return 1
        else:
            return -1
    assert a[0] != b[0]
    if a[0] < b[0]:
        return -1
    else:
        return 1


"""
"""
def single_abundance(HLA_cmpt,
                     HLA_length):
    def normalize2(prob, length):
        total = 0
        for allele, mass in prob.items():
            assert allele in length
            total += (mass / length[allele])
        for allele, mass in prob.items():
            assert allele in length
            prob[allele] = mass / length[allele] / total

    HLA_prob, HLA_prob_next = {}, {}
    for cmpt, count in HLA_cmpt.items():
        alleles = cmpt.split('-')
        for allele in alleles:
            if allele not in HLA_prob:
                HLA_prob[allele] = 0.0
            HLA_prob[allele] += (float(count) / len(alleles))

    # normalize2(HLA_prob, HLA_length)
    normalize(HLA_prob)
    def next_prob(HLA_cmpt, HLA_prob, HLA_length):
        HLA_prob_next = {}
        for cmpt, count in HLA_cmpt.items():
            alleles = cmpt.split('-')
            alleles_prob = 0.0
            for allele in alleles:
                assert allele in HLA_prob
                alleles_prob += HLA_prob[allele]
            for allele in alleles:
                if allele not in HLA_prob_next:
                    HLA_prob_next[allele] = 0.0
                HLA_prob_next[allele] += (float(count) * HLA_prob[allele] / alleles_prob)
        # normalize2(HLA_prob_next, HLA_length)
        normalize(HLA_prob_next)
        return HLA_prob_next

    diff, iter = 1.0, 0
    while diff > 0.0001 and iter < 1000:
        HLA_prob_next = next_prob(HLA_cmpt, HLA_prob, HLA_length)
        diff = prob_diff(HLA_prob, HLA_prob_next)
        HLA_prob = HLA_prob_next
        iter += 1
    for allele, prob in HLA_prob.items():
        allele_len = HLA_length[allele]
        HLA_prob[allele] /= float(allele_len)
    normalize(HLA_prob)
    HLA_prob = [[allele, prob] for allele, prob in HLA_prob.items()]
    HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)
    return HLA_prob

    
"""
"""
def joint_abundance(HLA_cmpt,
                    HLA_length):
    allele_names = set()
    for cmpt in HLA_cmpt.keys():
        allele_names |= set(cmpt.split('-'))
    
    HLA_prob, HLA_prob_next = {}, {}
    for cmpt, count in HLA_cmpt.items():
        alleles = cmpt.split('-')
        for allele1 in alleles:
            for allele2 in allele_names:
                if allele1 < allele2:
                    allele_pair = "%s-%s" % (allele1, allele2)
                else:
                    allele_pair = "%s-%s" % (allele2, allele1)
                if not allele_pair in HLA_prob:
                    HLA_prob[allele_pair] = 0.0
                HLA_prob[allele_pair] += (float(count) / len(alleles))

    if len(HLA_prob) <= 0:
        return HLA_prob

    # Choose top allele pairs
    def choose_top_alleles(HLA_prob):
        HLA_prob_list = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
        HLA_prob_list = sorted(HLA_prob_list, cmp=HLA_prob_cmp)
        HLA_prob = {}
        best_prob = HLA_prob_list[0][1]
        for i in range(len(HLA_prob_list)):
            allele_pair, prob = HLA_prob_list[i]
            if prob * 2 <= best_prob:
                break                        
            HLA_prob[allele_pair] = prob
        normalize(HLA_prob)
        return HLA_prob
    HLA_prob = choose_top_alleles(HLA_prob)

    def next_prob(HLA_cmpt, HLA_prob):
        HLA_prob_next = {}
        for cmpt, count in HLA_cmpt.items():
            alleles = cmpt.split('-')
            prob = 0.0
            for allele in alleles:
                for allele_pair in HLA_prob.keys():
                    if allele in allele_pair:
                        prob += HLA_prob[allele_pair]
            for allele in alleles:
                for allele_pair in HLA_prob.keys():
                    if not allele in allele_pair:
                        continue
                    if allele_pair not in HLA_prob_next:
                        HLA_prob_next[allele_pair] = 0.0
                    HLA_prob_next[allele_pair] += (float(count) * HLA_prob[allele_pair] / prob)
        normalize(HLA_prob_next)
        return HLA_prob_next

    diff, iter = 1.0, 0
    while diff > 0.0001 and iter < 1000:
        HLA_prob_next = next_prob(HLA_cmpt, HLA_prob)
        diff = prob_diff(HLA_prob, HLA_prob_next)
        HLA_prob = HLA_prob_next
        HLA_prob = choose_top_alleles(HLA_prob)
        iter += 1

    HLA_prob = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
    HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)
    return HLA_prob


"""
"""
def HLA_typing(ex_path,
               base_fname,
               simulation,
               reference_type,
               hla_list,
               partial,
               refHLAs,
               HLAs,
               HLA_names,
               HLA_lengths,
               refHLA_loci,
               Vars,
               Var_list,
               Links,
               exclude_allele_list,
               aligners,
               num_mismatch,
               fastq,
               read_fname,
               alignment_fname,
               threads,
               enable_coverage,
               best_alleles,
               verbose):

    def lower_bound(Var_list, pos):
        low, high = 0, len(Var_list)
        while low < high:
            m = (low + high) / 2
            m_pos = Var_list[m][0]
            if m_pos < pos:
                low = m + 1
            elif m_pos > pos:
                high = m
            else:
                assert m_pos == pos
                while m > 0:
                    if Var_list[m-1][0] < pos:
                        break
                    m -= 1
                return m
        return low        
            
    if simulation:
        test_passed = {}
    for aligner, index_type in aligners:
        if index_type == "graph":
            print("\n\t\t%s %s on %s" % (aligner, index_type, reference_type), file=sys.stderr)
        else:
            print("\n\t\t%s %s" % (aligner, index_type), file=sys.stderr)

        if alignment_fname == "":
            # Align reads, and sort the alignments into a BAM file
            align_reads(ex_path,
                        base_fname,
                        aligner,
                        index_type,
                        read_fname,
                        fastq,
                        threads,
                        verbose)
            
        for test_HLA_names in hla_list:
            if simulation:
                gene = test_HLA_names[0].split('*')[0]
            else:
                gene = test_HLA_names
            
            ref_allele = refHLAs[gene]
            ref_seq = HLAs[gene][ref_allele]
            ref_exons = refHLA_loci[gene][-1]

            # Read alignments
            alignview_cmd = ["samtools",
                             "view"]
            if alignment_fname == "":
                alignview_cmd += ["hla_input.bam"]
            else:
                if not os.path.exists(alignment_fname + ".bai"):
                    os.system("samtools index %s" % alignment_fname)
                alignview_cmd += [alignment_fname]
            base_locus = 0
            if index_type == "graph":
                if reference_type == "gene":
                    alignview_cmd += ["%s" % ref_allele]
                else:
                    assert reference_type in ["chromosome", "genome"]
                    _, chr, left, right, _ = refHLA_loci[gene]
                    base_locus = left
                    alignview_cmd += ["%s:%d-%d" % (chr, left + 1, right + 1)]

                bamview_proc = subprocess.Popen(alignview_cmd,
                                                stdout=subprocess.PIPE,
                                                stderr=open("/dev/null", 'w'))

                sort_read_cmd = ["sort", "-k", "1", "-n"]
                alignview_proc = subprocess.Popen(sort_read_cmd,
                                                  stdin=bamview_proc.stdout,
                                                  stdout=subprocess.PIPE,
                                                  stderr=open("/dev/null", 'w'))
            else:
                alignview_proc = subprocess.Popen(alignview_cmd,
                                             stdout=subprocess.PIPE,
                                             stderr=open("/dev/null", 'w'))

            # Count alleles
            HLA_counts, HLA_cmpt = {}, {}
            coverage = [0 for i in range(len(ref_seq) + 1)]
            num_reads, total_read_len = 0, 0
            prev_read_id = None
            prev_exon = False
            if index_type == "graph":
                # Cigar regular expression
                cigar_re = re.compile('\d+\w')
                for line in alignview_proc.stdout:
                    cols = line.strip().split()
                    read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
                    read_seq, qual = cols[9], cols[10]
                    num_reads += 1
                    total_read_len += len(read_seq)
                    flag, pos = int(flag), int(pos)
                    pos -= (base_locus + 1)
                    if pos < 0:
                        continue

                    if flag & 0x4 != 0:
                        continue

                    NM, Zs, MD = "", "", ""
                    for i in range(11, len(cols)):
                        col = cols[i]
                        if col.startswith("Zs"):
                            Zs = col[5:]
                        elif col.startswith("MD"):
                            MD = col[5:]
                        elif col.startswith("NM"):
                            NM = int(col[5:])

                    if NM > num_mismatch:
                        continue

                    # daehwan - for debugging purposes
                    debug = False
                    if read_id in ["2339"] and False:
                        debug = True
                        print("read_id: %s)" % read_id, pos, cigar_str, "NM:", NM, MD, Zs)
                        print("            ", read_seq)

                    vars = []
                    if Zs:
                        vars = Zs.split(',')

                    assert MD != ""
                    MD_str_pos, MD_len = 0, 0
                    read_pos, left_pos = 0, pos
                    right_pos = left_pos
                    cigars = cigar_re.findall(cigar_str)
                    cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
                    cmp_list = []
                    for i in range(len(cigars)):
                        cigar_op, length = cigars[i]
                        if cigar_op == 'M':
                            # Update coverage
                            if enable_coverage:
                                if right_pos + length < len(coverage):
                                    coverage[right_pos] += 1
                                    coverage[right_pos + length] -= 1
                                elif right_pos < len(coverage):
                                    coverage[right_pos] += 1
                                    coverage[-1] -= 1

                            first = True
                            MD_len_used = 0
                            while True:
                                if not first or MD_len == 0:
                                    if MD[MD_str_pos].isdigit():
                                        num = int(MD[MD_str_pos])
                                        MD_str_pos += 1
                                        while MD_str_pos < len(MD):
                                            if MD[MD_str_pos].isdigit():
                                                num = num * 10 + int(MD[MD_str_pos])
                                                MD_str_pos += 1
                                            else:
                                                break
                                        MD_len += num
                                # Insertion or full match followed
                                if MD_len >= length:
                                    MD_len -= length
                                    cmp_list.append(["match", right_pos + MD_len_used, length - MD_len_used])
                                    break
                                first = False
                                read_base = read_seq[read_pos + MD_len]
                                MD_ref_base = MD[MD_str_pos]
                                MD_str_pos += 1
                                assert MD_ref_base in "ACGT"
                                cmp_list.append(["match", right_pos + MD_len_used, MD_len - MD_len_used])
                                cmp_list.append(["mismatch", right_pos + MD_len, 1])
                                MD_len_used = MD_len + 1
                                MD_len += 1
                                # Full match
                                if MD_len == length:
                                    MD_len = 0
                                    break
                        elif cigar_op == 'I':
                            cmp_list.append(["insertion", right_pos, length])
                        elif cigar_op == 'D':
                            if MD[MD_str_pos] == '0':
                                MD_str_pos += 1
                            assert MD[MD_str_pos] == '^'
                            MD_str_pos += 1
                            while MD_str_pos < len(MD):
                                if not MD[MD_str_pos] in "ACGT":
                                    break
                                MD_str_pos += 1
                            cmp_list.append(["deletion", right_pos, length])
                        elif cigar_op == 'S':
                            cmp_list.append(["soft", right_pos, length])
                        else:                    
                            assert cigar_op == 'N'
                            cmp_list.append(["intron", right_pos, length])

                        if cigar_op in "MND":
                            right_pos += length

                        if cigar_op in "MIS":
                            read_pos += length

                    exon = False
                    for exon in ref_exons:
                        exon_left, exon_right = exon
                        if right_pos <= exon_left or pos > exon_right:
                            continue
                        else:
                            exon = True
                            break

                    if right_pos > len(ref_seq):
                        continue

                    def add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, exon = True):
                        max_count = max(HLA_count_per_read.values())
                        cur_cmpt = set()
                        for allele, count in HLA_count_per_read.items():
                            if count < max_count:
                                continue
                            if allele in exclude_allele_list:
                                continue                                
                            cur_cmpt.add(allele)                    
                            if not allele in HLA_counts:
                                HLA_counts[allele] = 1
                            else:
                                HLA_counts[allele] += 1

                        if len(cur_cmpt) == 0:
                            return

                        # daehwan - for debugging purposes                            
                        alleles = ["", ""]
                        # alleles = ["B*40:304", "B*40:02:01"]
                        allele1_found, allele2_found = False, False
                        for allele, count in HLA_count_per_read.items():
                            if count < max_count:
                                continue
                            if allele == alleles[0]:
                                allele1_found = True
                            elif allele == alleles[1]:
                                allele2_found = True
                        if allele1_found != allele2_found:
                            print(alleles[0], HLA_count_per_read[alleles[0]])
                            print(alleles[1], HLA_count_per_read[alleles[1]])
                            if allele1_found:
                                print("%s\tread_id %s - %d vs. %d]" % (alleles[0], prev_read_id, max_count, HLA_count_per_read[alleles[1]]))
                            else:
                                print("%s\tread_id %s - %d vs. %d]" % (alleles[1], prev_read_id, max_count, HLA_count_per_read[alleles[0]]))
                            print(read_seq)

                        cur_cmpt = sorted(list(cur_cmpt))
                        cur_cmpt = '-'.join(cur_cmpt)
                        add = 1
                        if partial and not exon:
                            add *= 0.2
                        if not cur_cmpt in HLA_cmpt:
                            HLA_cmpt[cur_cmpt] = add
                        else:
                            HLA_cmpt[cur_cmpt] += add

                    if read_id != prev_read_id:
                        if prev_read_id != None:
                            add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, prev_exon)

                        HLA_count_per_read = {}
                        for HLA_name in HLA_names[gene]:
                            if HLA_name.find("BACKBONE") != -1:
                                continue
                            HLA_count_per_read[HLA_name] = 0

                    def add_count(var_id, add):
                        assert var_id in Links
                        alleles = Links[var_id]
                        for allele in alleles:
                            if allele.find("BACKBONE") != -1:
                                continue
                            HLA_count_per_read[allele] += add
                            # daehwan - for debugging purposes
                            if debug:
                                if allele in ["DQA1*05:05:01:01", "DQA1*05:05:01:02"]:
                                    print(allele, add, var_id)

                    # Decide which allele(s) a read most likely came from
                    # also sanity check - read length, cigar string, and MD string
                    for var_id, data in Vars[gene].items():
                        var_type, var_pos, var_data = data
                        if var_type != "deletion":
                            continue
                        if left_pos >= var_pos and right_pos <= var_pos + int(var_data):
                            add_count(var_id, -1)                            
                    ref_pos, read_pos, cmp_cigar_str, cmp_MD = left_pos, 0, "", ""
                    cigar_match_len, MD_match_len = 0, 0            
                    for cmp in cmp_list:
                        type = cmp[0]
                        length = cmp[2]
                        if type == "match":
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos + length <= var_pos:
                                    break
                                if ref_pos <= var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "insertion":
                                        if ref_pos < var_pos and ref_pos + length > var_pos + len(var_data):
                                            add_count(var_id, -1)
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print(cmp, var_id, Links[var_id])
                                    elif var_type == "deletion":
                                        del_len = int(var_data)
                                        if ref_pos < var_pos and ref_pos + length > var_pos + del_len:
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print(cmp, var_id, Links[var_id], -1, Vars[gene][var_id])
                                            # Check if this might be one of the two tandem repeats (the same left coordinate)
                                            cmp_left, cmp_right = cmp[1], cmp[1] + cmp[2]
                                            test1_seq1 = ref_seq[cmp_left:cmp_right]
                                            test1_seq2 = ref_seq[cmp_left:var_pos] + ref_seq[var_pos + del_len:cmp_right + del_len]
                                            # Check if this happens due to small repeats (the same right coordinate - e.g. 19 times of TTTC in DQA1*05:05:01:02)
                                            cmp_left -= read_pos
                                            cmp_right += (len(read_seq) - read_pos - cmp[2])
                                            test2_seq1 = ref_seq[cmp_left+int(var_data):cmp_right]
                                            test2_seq2 = ref_seq[cmp_left:var_pos] + ref_seq[var_pos+int(var_data):cmp_right]
                                            if test1_seq1 != test1_seq2 and test2_seq1 != test2_seq2:
                                                add_count(var_id, -1)
                                    else:
                                        if debug:
                                            print(cmp, var_id, Links[var_id], -1)
                                        add_count(var_id, -1)
                                var_idx += 1

                            read_pos += length
                            ref_pos += length
                            cigar_match_len += length
                            MD_match_len += length
                        elif type == "mismatch":
                            read_base = read_seq[read_pos]
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "single":
                                        if var_data == read_base:
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print(cmp, var_id, 1, var_data, read_base, Links[var_id])

                                            # daehwan - for debugging purposes
                                            if False:
                                                read_qual = ord(qual[read_pos])
                                                add_count(var_id, (read_qual - 60) / 60.0)
                                            else:
                                                add_count(var_id, 1)
                                        # daehwan - check out if this routine is appropriate
                                        # else:
                                        #    add_count(var_id, -1)
                                var_idx += 1

                            cmp_MD += ("%d%s" % (MD_match_len, ref_seq[ref_pos]))
                            MD_match_len = 0
                            cigar_match_len += 1
                            read_pos += 1
                            ref_pos += 1
                        elif type == "insertion":
                            ins_seq = read_seq[read_pos:read_pos+length]
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            # daehwan - for debugging purposes
                            if debug:
                                print(left_pos, cigar_str, MD, vars)
                                print(ref_pos, ins_seq, Var_list[gene][var_idx], Vars[gene][Var_list[gene][var_idx][1]])
                                # sys.exit(1)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "insertion":                                
                                        if var_data == ins_seq:
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print(cmp, var_id, 1, Links[var_id])
                                            add_count(var_id, 1)
                                var_idx += 1

                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            read_pos += length
                            cmp_cigar_str += ("%dI" % length)
                        elif type == "deletion":
                            del_len = length
                            # Deletions can be shifted bidirectionally
                            temp_ref_pos = ref_pos
                            while temp_ref_pos > 0:
                                last_bp = ref_seq[temp_ref_pos + del_len - 1]
                                prev_bp = ref_seq[temp_ref_pos - 1]
                                if last_bp != prev_bp:
                                    break
                                temp_ref_pos -= 1
                            var_idx = lower_bound(Var_list[gene], temp_ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if temp_ref_pos < var_pos:
                                    first_bp = ref_seq[temp_ref_pos]
                                    next_bp = ref_seq[temp_ref_pos + del_len]
                                    if first_bp == next_bp:
                                        temp_ref_pos += 1
                                        continue
                                    else:
                                        break
                                if temp_ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "deletion":
                                        var_len = int(var_data)
                                        if var_len == length:
                                            if debug:
                                                print(cmp, var_id, 1, Links[var_id])
                                                print(ref_seq[var_pos - 10:var_pos], ref_seq[var_pos:var_pos+int(var_data)], ref_seq[var_pos+int(var_data):var_pos+int(var_data)+10])
                                            add_count(var_id, 1)
                                var_idx += 1

                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_MD += ("%d" % MD_match_len)
                            MD_match_len = 0
                            cmp_cigar_str += ("%dD" % length)
                            cmp_MD += ("^%s" % ref_seq[ref_pos:ref_pos+length])
                            ref_pos += length
                        elif type == "soft":
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            read_pos += length
                            cmp_cigar_str += ("%dS" % length)
                        else:
                            assert type == "intron"
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_cigar_str += ("%dN" % length)
                            ref_pos += length                    
                    if cigar_match_len > 0:
                        cmp_cigar_str += ("%dM" % cigar_match_len)
                    cmp_MD += ("%d" % MD_match_len)
                    if read_pos != len(read_seq) or \
                            cmp_cigar_str != cigar_str or \
                            cmp_MD != MD:
                        print(("Error:", cigar_str, MD), file=sys.stderr)
                        print(("\tcomputed:", cmp_cigar_str, cmp_MD), file=sys.stderr)
                        print(("\tcmp list:", cmp_list), file=sys.stderr)
                        assert False            

                    prev_read_id = read_id
                    prev_exon = exon

                if num_reads <= 0:
                    continue

                if prev_read_id != None:
                    add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read)

                # Coverage
                # it is not used by the default
                if enable_coverage:
                    assert num_reads > 0
                    read_len = int(total_read_len / float(num_reads))
                    coverage_sum = 0
                    for i in range(len(coverage)):
                        if i > 0:
                            coverage[i] += coverage[i-1]
                        coverage_sum += coverage[i]
                    coverage_avg = coverage_sum / float(len(coverage))
                    assert len(ref_seq) < len(coverage)
                    for i in range(len(ref_seq)):
                        coverage_threshold = 1.0 * coverage_avg
                        if i < read_len:
                            coverage_threshold *= ((i+1) / float(read_len))
                        elif i + read_len > len(ref_seq):
                            coverage_threshold *= ((len(ref_seq) - i) / float(read_len))
                        if coverage[i] >= coverage_threshold:
                            continue
                        pseudo_num_reads = (coverage_threshold - coverage[i]) / read_len
                        var_idx = lower_bound(Var_list[gene], i + 1)
                        if var_idx >= len(Var_list[gene]):
                            var_idx = len(Var_list[gene]) - 1
                        cur_cmpt = set()
                        while var_idx >= 0:
                            var_pos, var_id = Var_list[gene][var_idx]
                            var_type, _, var_data = Vars[gene][var_id]
                            if var_type == "deletion":
                                del_len = int(var_data)
                                if i < var_pos:
                                    break
                                if i + read_len < var_pos + int(var_data):
                                    assert var_id in Links
                                    cur_cmpt = cur_cmpt.union(set(Links[var_id]))
                            var_idx -= 1
                        if cur_cmpt:
                            cur_cmpt = '-'.join(list(cur_cmpt))
                            if not cur_cmpt in HLA_cmpt:
                                HLA_cmpt[cur_cmpt] = 0
                            HLA_cmpt[cur_cmpt] += pseudo_num_reads
            else:
                assert index_type == "linear"
                def add_alleles(alleles):
                    if not allele in HLA_counts:
                        HLA_counts[allele] = 1
                    else:
                        HLA_counts[allele] += 1

                    cur_cmpt = sorted(list(alleles))
                    cur_cmpt = '-'.join(cur_cmpt)
                    if not cur_cmpt in HLA_cmpt:
                        HLA_cmpt[cur_cmpt] = 1
                    else:
                        HLA_cmpt[cur_cmpt] += 1

                prev_read_id, prev_AS = None, None
                alleles = set()
                for line in alignview_proc.stdout:
                    cols = line[:-1].split()
                    read_id, flag, allele = cols[:3]
                    flag = int(flag)
                    if flag & 0x4 != 0:
                        continue
                    if not allele.startswith(gene):
                        continue
                    if allele.find("BACKBONE") != -1:
                        continue

                    AS = None
                    for i in range(11, len(cols)):
                        col = cols[i]
                        if col.startswith("AS"):
                            AS = int(col[5:])
                    assert AS != None
                    if read_id != prev_read_id:
                        if alleles:
                            if aligner == "hisat2" or \
                                    (aligner == "bowtie2" and len(alleles) < 10):
                                add_alleles(alleles)
                            alleles = set()
                        prev_AS = None
                    if prev_AS != None and AS < prev_AS:
                        continue
                    prev_read_id = read_id
                    prev_AS = AS
                    alleles.add(allele)

                if alleles:
                    add_alleles(alleles)

            HLA_counts = [[allele, count] for allele, count in HLA_counts.items()]
            def HLA_count_cmp(a, b):
                if a[1] != b[1]:
                    return b[1] - a[1]
                assert a[0] != b[0]
                if a[0] < b[0]:
                    return -1
                else:
                    return 1
            HLA_counts = sorted(HLA_counts, cmp=HLA_count_cmp)
            for count_i in range(len(HLA_counts)):
                count = HLA_counts[count_i]
                if simulation:
                    found = False
                    for test_HLA_name in test_HLA_names:
                        if count[0] == test_HLA_name:
                            print("\t\t\t*** %d ranked %s (count: %d)" % (count_i + 1, test_HLA_name, count[1]))
                            found = True
                            """
                            if count_i > 0 and HLA_counts[0][1] > count[1]:
                                print >> sys.stderr, "Warning: %s ranked first (count: %d)" % (HLA_counts[0][0], HLA_counts[0][1])
                                assert False
                            else:
                                test_passed += 1
                            """
                    if count_i < 5 and not found:
                        print("\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1]))
                else:
                    print("\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1]))
                    if count_i >= 9:
                        break
            print("\n")

            HLA_prob = single_abundance(HLA_cmpt, HLA_lengths[gene])

            success = [False for i in range(len(test_HLA_names))]
            found_list = [False for i in range(len(test_HLA_names))]
            for prob_i in range(len(HLA_prob)):
                prob = HLA_prob[prob_i]
                found = False
                if simulation:
                    for name_i in range(len(test_HLA_names)):
                        test_HLA_name = test_HLA_names[name_i]
                        if prob[0] == test_HLA_name:
                            rank_i = prob_i
                            while rank_i > 0:
                                if prob == HLA_prob[rank_i - 1][1]:
                                    rank_i -= 1
                                else:
                                    break
                            print("\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (rank_i + 1, test_HLA_name, prob[1] * 100.0))
                            if rank_i < len(success):
                                success[rank_i] = True
                            found_list[name_i] = True
                            found = True                        
                    if not False in found_list:
                        break
                if not found:
                    print("\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, prob[0], prob[1] * 100.0))
                    if best_alleles and prob_i < 2:
                        print("SingleModel %s (abundance: %.2f%%)" % (prob[0], prob[1] * 100.0))
                if not simulation and prob_i >= 9:
                    break
            print("\n")

            if len(test_HLA_names) == 2 or not simulation:
                HLA_prob = joint_abundance(HLA_cmpt, HLA_lengths[gene])
                if len(HLA_prob) <= 0:
                    continue
                success = [False]
                for prob_i in range(len(HLA_prob)):
                    allele_pair, prob = HLA_prob[prob_i]
                    allele1, allele2 = allele_pair.split('-')
                    if best_alleles and prob_i < 1:
                        print("PairModel %s (abundance: %.2f%%)" % (allele_pair, prob * 100.0))
                    if simulation:
                        if allele1 in test_HLA_names and allele2 in test_HLA_names:
                            rank_i = prob_i
                            while rank_i > 0:
                                if HLA_prob[rank_i-1][1] == prob:                                        
                                    rank_i -= 1
                                else:
                                    break
                            print("\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (rank_i + 1, allele_pair, prob * 100.0))
                            if rank_i == 0:
                                success[0] = True
                            break
                    print("\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, allele_pair, prob * 100.0))
                    if not simulation and prob_i >= 9:
                        break
                print("\n")

                # Li's method
                """
                li_hla = os.path.join(ex_path, "li_hla/hla")
                if os.path.exists(li_hla):
                    li_hla_cmd = [li_hla,
                                  "hla",
                                  "hla_input.bam",
                                  "-b", "%s*BACKBONE" % gene]
                    li_hla_proc = subprocess.Popen(li_hla_cmd,
                                                   stdout=subprocess.PIPE,
                                                   stderr=open("/dev/null", 'w'))

                    # read in the result of Li's hla
                    for line in li_hla_proc.stdout:
                        allele1, allele2, score = line.strip().split()
                        score = float(score)
                        if simulation:
                            if allele1 in test_HLA_names and allele2 in test_HLA_names:
                                print >> sys.stderr, "\t\t\t*** 1 ranked %s-%s (score: %.2f)" % (allele1, allele2, score)
                                success[0] = True
                            else:
                                print >> sys.stderr, "\t\t\tLiModel fails"
                        if best_alleles:
                            print >> sys.stdout, "LiModel %s-%s (score: %.2f)" % (allele1, allele2, score)
                    li_hla_proc.communicate()
                """

            if simulation and not False in success:
                aligner_type = "%s %s" % (aligner, index_type)
                if not aligner_type in test_passed:
                    test_passed[aligner_type] = 1
                else:
                    test_passed[aligner_type] += 1

    if simulation:
        return test_passed


def read_HLA_alleles(fname, HLAs):
    for line in open(fname):
        if line.startswith(">"):
            HLA_name = line.strip().split()[0][1:]
            HLA_gene = HLA_name.split('*')[0]
            if not HLA_gene in HLAs:
                HLAs[HLA_gene] = {}
            if not HLA_name in HLAs[HLA_gene]:
                HLAs[HLA_gene][HLA_name] = ""
        else:
            HLAs[HLA_gene][HLA_name] += line.strip()
    return HLAs

"""
"""
def genotyping(base_fname,
               reference_type,
               hla_list,
               partial,
               aligners,
               read_fname,
               alignment_fname,
               threads,
               simulate_interval,
               enable_coverage,
               best_alleles,
               exclude_allele_list,
               default_allele_list,
               num_mismatch,
               verbose,
               daehwan_debug):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(genotyping))
    ex_path = os.path.dirname(curr_script)

    # Clone a git repository, IMGTHLA
    if not os.path.exists("IMGTHLA"):
        os.system("git clone https://github.com/jrob119/IMGTHLA.git")

    # Clone hisat2 genotype database, hisat_genotype_db
    """
    if not os.path.exists("hisat_genotype_db"):
        os.system("git clone https://github.com/infphilo/hisat_genotype_db.git")
    """

    simulation = (read_fname == [] and alignment_fname == "")

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        hisat2_inspect = os.path.join(ex_path, "hisat2-inspect")
        os.system("%s grch38/genome > genome.fa" % hisat2_inspect)
        os.system("samtools faidx genome.fa")

    # Check if the pre-existing files (hla*) are compatible with the current parameter setting
    if os.path.exists("%s.ref" % base_fname):
        left = 0
        HLA_genes = set()
        BACKBONE = False
        for line in open("%s.ref" % base_fname):
            HLA_name = line.strip().split()[0]
            if HLA_name.find("BACKBONE") != -1:
                BACKBONE = True
            HLA_gene = HLA_name.split('*')[0]
            HLA_genes.add(HLA_gene)
        delete_hla_files = False
        if reference_type == "gene":
            if not BACKBONE:
                delete_hla_files = True
        elif reference_type in ["chromosome", "genome"]:
            if BACKBONE:
                delete_hla_files = True
        else:
            assert False
        if not set(hla_list).issubset(HLA_genes):
            delete_hla_files = True
        if base_fname == "hla":
            if delete_hla_files:
                os.system("rm %s*" % base_fname)
    
    # Extract HLA variants, backbone sequence, and other sequeces  
    HLA_fnames = [base_fname+"_backbone.fa",
                  base_fname+"_sequences.fa",
                  base_fname+".ref",
                  base_fname+".snp",
                  base_fname+".haplotype",
                  base_fname+".link",
                  base_fname+"_alleles_excluded.txt"]
    
    # Check if excluded alleles in current files match
    excluded_alleles_match = False
    if(os.path.exists(HLA_fnames[6])):
        afile = open(HLA_fnames[6],'r')
        afile.readline()
        lines = afile.read().split()
        excluded_alleles_match = (set(exclude_allele_list) == set(lines))
        afile.close()
    elif len(exclude_allele_list) == 0:
        excluded_alleles_match = True
        try:
            temp_name = HLA_fnames[6]
            HLA_fnames.remove(HLA_fnames[6])
            os.remove(temp_name)
        except OSError:
            pass
        
    if not excluded_alleles_match:
        print("Creating Allele Exclusion File.\n")
        afile = open(HLA_fnames[6],'w')
        afile.write("Alleles excluded:\n")
        afile.write("\n".join(exclude_allele_list))
        afile.close()
        
    if (not check_files(HLA_fnames)) or (not excluded_alleles_match) :
        extract_hla_script = os.path.join(ex_path, "hisatgenotype_extract_vars.py")
        extract_cmd = [extract_hla_script,
                       "--base", base_fname,
                       "--reference-type", reference_type]

        if base_fname == "hla":
            extract_cmd += ["--hla-list", ','.join(hla_list)]

        if len(exclude_allele_list) > 0:
            print(exclude_allele_list)
            extract_cmd += ["--exclude-allele-list", ",".join(exclude_allele_list)]

        if len(base_fname) > 3:
            extract_cmd += ["--base", base_fname]

        if partial:
            extract_cmd += ["--partial"]
        extract_cmd += ["--inter-gap", "30",
                        "--intra-gap", "50"]
        if verbose:
            print("\tRunning:", ' '.join(extract_cmd), file=sys.stderr)
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        
        if not check_files(HLA_fnames):
            print("Error: extract_HLA_vars failed!", file=sys.stderr)
            sys.exit(1)
            
    for aligner, index_type in aligners:
        # Build HISAT2 graph indexes based on the above information
        if aligner == "hisat2" and index_type == "graph":
            HLA_hisat2_graph_index_fnames = ["%s.graph.%d.ht2" % (base_fname, i+1) for i in range(8)]
            if not check_files(HLA_hisat2_graph_index_fnames) or (not excluded_alleles_match):
                hisat2_build = os.path.join(ex_path, "hisat2-build")
                build_cmd = [hisat2_build,
                             "-p", str(threads),
                             "--snp", HLA_fnames[3],
                             "--haplotype", HLA_fnames[4] ,
                             HLA_fnames[0],
                             "%s.graph" % base_fname]
                if verbose:
                    print("\tRunning:", ' '.join(build_cmd), file=sys.stderr)
                proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                proc.communicate()        
                if not check_files(HLA_hisat2_graph_index_fnames):
                    print("Error: indexing HLA failed!  Perhaps, you may have forgotten to build hisat2 executables?", file=sys.stderr)
                    sys.exit(1)

        # Build HISAT2 linear indexes based on the above information
        elif aligner == "hisat2" and index_type == "linear":
            HLA_hisat2_linear_index_fnames = ["%s.linear.%d.ht2" % (base_fname, i+1) for i in range(8)]
            if reference_type == "gene" and (not check_files(HLA_hisat2_linear_index_fnames) or (not excluded_alleles_match)):
                hisat2_build = os.path.join(ex_path, "hisat2-build")
                build_cmd = [hisat2_build,
                             "%s,%s"%(HLA_fnames[0],HLA_fnames[1]),
                             "%s.linear" % base_fname]
                proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                proc.communicate()        
                if not check_files(HLA_hisat2_linear_index_fnames):
                    print("Error: indexing HLA failed!", file=sys.stderr)
                    sys.exit(1)

        # Build Bowtie2 indexes based on the above information
        else:
            assert aligner == "bowtie2" and index_type == "linear"
            HLA_bowtie2_index_fnames = ["%s.%d.bt2" % (base_fname, i+1) for i in range(4)]
            HLA_bowtie2_index_fnames += ["%s.rev.%d.bt2" % (base_fname, i+1) for i in range(2)]
            if reference_type == "gene" and (not check_files(HLA_bowtie2_index_fnames) or (not excluded_alleles_match)):
                build_cmd = ["bowtie2-build",
                             "%s,%s"%(HLA_fnames[0],HLA_fnames[1]),
                             base_fname]
                proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'))
                proc.communicate()        
                if not check_files(HLA_bowtie2_index_fnames):
                    print("Error: indexing HLA failed!", file=sys.stderr)
                    sys.exit(1)
        
    # Read partial alleles from hla.data (temporary)
    partial_alleles = set()
    if base_fname == "hla":
        for line in open("IMGTHLA/hla.dat"):
            if not line.startswith("DE"):
                continue
            allele_name = line.split()[1][4:-1]
            gene = allele_name.split('*')[0]
            if line.find("partial") != -1:
                partial_alleles.add(allele_name)

    if len(default_allele_list) != 0:
        #print os.getcwd()
        if not os.path.exists("./Default-HLA/hla_backbone.fa"):
            #current_path = os.getcwd()
            try:
                os.mkdir("./Default-HLA")
            except:
                pass
            #os.chdir(current_path + "/Default-HLA")
            
            extract_hla_script = os.path.join(ex_path, "hisat2_extract_HLA_vars.py")
            extract_cmd = [extract_hla_script,
                           "--reference-type", reference_type,
                           "--hla-list", ','.join(hla_list),
                           "--base", "./Default-HLA/hla"]

            if partial:
                extract_cmd += ["--partial"]
            extract_cmd += ["--inter-gap", "30",
                            "--intra-gap", "50"]
            if verbose:
                print("\tRunning:", ' '.join(extract_cmd), file=sys.stderr)
            proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            proc.communicate()
            
            if not os.path.exists("./Default-HLA/hla_backbone.fa"):
                print("Error: extract_HLA_vars (Default) failed!", file=sys.stderr)
                sys.exit(1)
    
    # Read HLA alleles (names and sequences)
    refHLAs, refHLA_loci = {}, {}
    for line in open("%s.ref" % base_fname):
        HLA_name, chr, left, right, length, exon_str = line.strip().split()
        HLA_gene = HLA_name.split('*')[0]
        assert not HLA_gene in refHLAs
        refHLAs[HLA_gene] = HLA_name
        left, right = int(left), int(right)
        exons = []
        for exon in exon_str.split(','):
            exon_left, exon_right = exon.split('-')
            exons.append([int(exon_left), int(exon_right)])
        refHLA_loci[HLA_gene] = [HLA_name, chr, left, right, exons]
        
    HLAs = {}
    if reference_type == "gene":
        read_HLA_alleles(HLA_fnames[0], HLAs)
    read_HLA_alleles(HLA_fnames[1], HLAs)
    
    # HLA gene alleles
    HLA_names = {}
    for HLA_gene, data in HLAs.items():
        HLA_names[HLA_gene] = list(data.keys())

    # HLA gene allele lengths
    HLA_lengths = {}
    for HLA_gene, HLA_alleles in HLAs.items():
        HLA_lengths[HLA_gene] = {}
        for allele_name, seq in HLA_alleles.items():
            HLA_lengths[HLA_gene][allele_name] = len(seq)

    # Construct excluded alleles (Via default backbone data)
    custom_allele_check = False
    if len(default_allele_list) > 0:
        custom_allele_check = True
        HLAs_default = {}
        read_HLA_alleles("./Default-HLA/hla_backbone.fa",HLAs_default)
        read_HLA_alleles("./Default-HLA/hla_sequences.fa",HLAs_default)
        
        for HLA_gene, HLA_alleles in HLAs_default.items():
            for allele_name, seq in HLA_alleles.items():
                if allele_name in default_allele_list:
                    HLA_lengths[HLA_gene][allele_name] = len(seq)

    # Read HLA variants, and link information
    Vars, Var_list = {}, {}
    for line in open("%s.snp" % base_fname):
        var_id, var_type, allele, pos, data = line.strip().split('\t')
        pos = int(pos)
        if reference_type != "gene":
            allele, dist = None, 0
            for tmp_gene, values in refHLA_loci.items():
                allele_name, chr, left, right, exons = values
                if allele == None:
                    allele = allele_name
                    dist = abs(pos - left)
                else:
                    if dist > abs(pos - left):
                        allele = allele_name
                        dist = abs(pos - left)
            
        gene = allele.split('*')[0]
        if not gene in Vars:
            Vars[gene] = {}
            assert not gene in Var_list
            Var_list[gene] = []
            
        assert not var_id in Vars[gene]
        left = 0
        if reference_type != "gene":
            _, _, left, _, _ = refHLA_loci[gene]
        Vars[gene][var_id] = [var_type, pos - left, data]
        Var_list[gene].append([pos - left, var_id])
        
    for gene, in_var_list in Var_list.items():
        Var_list[gene] = sorted(in_var_list)
        
    Links = {}
    for line in open("%s.link" % base_fname):
        var_id, alleles = line.strip().split('\t')
        alleles = alleles.split()
        assert not var_id in Links
        Links[var_id] = alleles

    # Test HLA typing
    test_list = []
    if simulation:
        basic_test, pair_test = True, False
        if daehwan_debug:
            if "basic_test" in daehwan_debug:
                basic_test, pair_test = True, False
            else:
                basic_test, pair_test = False, True

        test_passed = {}
        test_list = []
        if base_fname == "hla":
            genes = list(set(hla_list) & set(HLA_names.keys()))
        else:
            genes = HLA_names.keys()
            
        if basic_test:
            for gene in genes:
                HLA_gene_alleles = HLA_names[gene]
                for HLA_name in HLA_gene_alleles:
                    if HLA_name.find("BACKBONE") != -1:
                        continue
                    test_list.append([[HLA_name]])
        if pair_test:
            test_size = 500
            allele_count = 2
            for test_i in range(test_size):
                test_pairs = []
                for gene in genes:
                    HLA_gene_alleles = []                    
                    for allele in HLA_names[gene]:
                        if allele.find("BACKBONE") != -1:
                            continue
                        HLA_gene_alleles.append(allele)

                    # DK - temporary
                    if len(HLA_gene_alleles) < 2:
                        continue
                        
                    nums = [i for i in range(len(HLA_gene_alleles))]
                    random.shuffle(nums)
                    test_pairs.append(sorted([HLA_gene_alleles[nums[i]] for i in range(allele_count)]))
                test_list.append(test_pairs)

        for test_i in range(len(test_list)):
            if "test_id" in daehwan_debug:
                daehwan_test_ids = daehwan_debug["test_id"].split('-')
                if str(test_i + 1) not in daehwan_test_ids:
                    continue

            print("Test %d" % (test_i + 1), file=sys.stderr)
            test_HLA_list = test_list[test_i]
           
            # daehwan - for debugging purposes
            # test_HLA_list = [["A*11:50Q", "A*11:01:01:01", "A*01:01:01:01"]]
            for test_HLA_names in test_HLA_list:
                for test_HLA_name in test_HLA_names:
                    if custom_allele_check:
                        gene = test_HLA_name.split('*')[0]
                        test_HLA_seq = HLAs_default[gene][test_HLA_name]
                        seq_type = "partial" if test_HLA_name in partial_alleles else "full"
                        print("\t%s - %d bp (%s sequence)" % (test_HLA_name, len(test_HLA_seq), seq_type), file=sys.stderr)
                        continue
                    gene = test_HLA_name.split('*')[0]
                    test_HLA_seq = HLAs[gene][test_HLA_name]
                    seq_type = "partial" if test_HLA_name in partial_alleles else "full"
                    print("\t%s - %d bp (%s sequence)" % (test_HLA_name, len(test_HLA_seq), seq_type), file=sys.stderr)
            if custom_allele_check:
                simulate_reads(HLAs_default, test_HLA_list, simulate_interval)
            else:
                simulate_reads(HLAs, test_HLA_list, simulate_interval)

            if "test_id" in daehwan_debug:
                read_fname = ["hla_input_1.fa"]
            else:
                read_fname = ["hla_input_1.fa", "hla_input_2.fa"]

            fastq = False
            
            tmp_test_passed = HLA_typing(ex_path,
                                         base_fname,
                                         simulation,
                                         reference_type,
                                         test_HLA_list,
                                         partial,
                                         refHLAs,
                                         HLAs,                       
                                         HLA_names,
                                         HLA_lengths,
                                         refHLA_loci,
                                         Vars,
                                         Var_list,
                                         Links,
                                         exclude_allele_list,
                                         aligners,
                                         num_mismatch,
                                         fastq,
                                         read_fname,
                                         alignment_fname,
                                         threads,
                                         enable_coverage,
                                         best_alleles,
                                         verbose)

            for aligner_type, passed in tmp_test_passed.items():
                if aligner_type in test_passed:
                    test_passed[aligner_type] += passed
                else:
                    test_passed[aligner_type] = passed

                print("\t\tPassed so far: %d/%d (%.2f%%)" % (test_passed[aligner_type], test_i + 1, (test_passed[aligner_type] * 100.0 / (test_i + 1))), file=sys.stderr)


        for aligner_type, passed in test_passed.items():
            print("%s:\t%d/%d passed (%.2f%%)" % (aligner_type, passed, len(test_list), passed * 100.0 / len(test_list)), file=sys.stderr)
    
    else: # With real reads or BAMs
        if base_fname == "hla":
            gene_list = hla_list
        else:
            gene_list = Vars.keys()
        print("\t", ' '.join(gene_list), file=sys.stderr)

        fastq = True
        HLA_typing(ex_path,
                   base_fname,
                   simulation,
                   reference_type,
                   gene_list,
                   partial,
                   refHLAs,
                   HLAs,                       
                   HLA_names,
                   HLA_lengths,
                   refHLA_loci,
                   Vars,
                   Var_list,
                   Links,
                   exclude_allele_list,
                   aligners,
                   num_mismatch,
                   fastq,
                   read_fname,
                   alignment_fname,
                   threads,
                   enable_coverage,
                   best_alleles,
                   verbose)

        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='genotyping')

    # Add Arguments
    arguments.args_databases(parser)
    arguments.args_set_aligner(parser)
    arguments.args_bamfile(parser)
    arguments.args_reference_type(parser)
    arguments.args_hla_cyp(parser)
    arguments.args_common(parser, 
                          debug = True) # Need Debug option

    args = parser.parse_args()
    if not args.base_fname:
        args.base_fname = 'hla'
    if not args.locus_list:
        args.locus_list = 'A,B,C,DQA1,DQB1,DRB1'

    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print("Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type), file=sys.stderr)
        sys.exit(1)

    args.hla_list = args.locus_list.split(',')
    if args.aligners == "":
        print("Error: --aligners must be non-empty.", file=sys.stderr)
        sys.exit(1)    

    if args.aligners:
        args.aligners = args.aligners.split(',')
        for i in range(len(args.aligners)):
            args.aligners[i] = args.aligners[i].split('.')
    else:
        args.aligners = [[args.aligner, 'graph' if args.graph_index else 'linear' ]]

    if args.read_fname:
        args.read_fname = args.read_fname.split(',')
    else:
        args.read_fname = []
    if args.alignment_fname != "" and \
            not os.path.exists(args.alignment_fname):
        print("Error: %s doesn't exist." % args.alignment_fname, file=sys.stderr)
        sys.exit(1)
    
    if len(args.default_allele_list) > 0:
        args.default_allele_list = args.default_allele_list.split(',')
        
    if len(args.exclude_allele_list) > 0:
        if args.exclude_allele_list.strip().isdigit():
            num_alleles = int(args.exclude_allele_list)
            if not os.path.exists("./Default-HLA/hla_backbone.fa"):
                try:
                    os.mkdir("./Default-HLA")
                except:
                    pass
                
                extract_hla_script = os.path.join(args.ex_path, "hisat2_extract_HLA_vars.py")
                extract_cmd = [extract_hla_script,
                               "--reference-type", args.reference_type,
                               "--hla-list", ','.join(args.hla_list),
                               "--base", "./Default-HLA/hla"]
                if args.partial:
                    extract_cmd += ["--partial"]
                extract_cmd += ["--inter-gap", "30",
                                "--intra-gap", "50"]
                if args.verbose:
                    print("\tRunning:", ' '.join(extract_cmd), file=sys.stderr)
                proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                proc.communicate()
                if not os.path.exists("./Default-HLA/hla_backbone.fa"):
                    print("Error: extract_HLA_vars (Default) failed!", file=sys.stderr)
                    sys.exit(1)
       
            HLAs_default = {}
            #read_HLA_alleles("./Default-HLA/hla_backbone.fa",HLAs_default)
            read_HLA_alleles("./Default-HLA/hla_sequences.fa",HLAs_default)
    
            allele_names = list(HLAs_default['A'].keys())
            random.shuffle(allele_names)
            args.exclude_allele_list = allele_names[0:num_alleles]
            args.default_allele_list = allele_names[num_alleles:2*num_alleles]
            
            args.default_allele_list = args.default_allele_list + args.exclude_allele_list
        else:
            args.exclude_allele_list = args.exclude_allele_list.split(',')
        
    debug = {}
    if args.debug != "":
        for item in args.debug.split(','):
            if ':' in item:
                key, value = item.split(':')
                debug[key] = value
            else:
                debug[item] = 1

    random.seed(1)
    genotyping(args.base_fname,
               args.reference_type,
               args.hla_list,
               args.partial,
               args.aligners,
               args.read_fname,
               args.alignment_fname,
               args.threads,
               args.simulate_interval,
               args.coverage,
               args.best_alleles,
               args.exclude_allele_list,
               args.default_allele_list,
               args.num_mismatch,
               args.verbose,
               debug)

    
