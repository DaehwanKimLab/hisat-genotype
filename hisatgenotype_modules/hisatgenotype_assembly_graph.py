#!/usr/bin/env python
# --------------------------------------------------------------------------- #
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>                            #
#                                                                             #
# This file is part of HISAT-genotype. It contains the core algorithms for    #
# assembly and phasing.                                                       #
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
import math
import random
from datetime import datetime, date, time
from collections import deque
from copy import deepcopy

# --------------------------------------------------------------------------- #
# Basic Functions that are used in main classes to build de Bruijn graph      #
# --------------------------------------------------------------------------- #
"""
Find the most prevelent nucleotide in a nucleotide dictionary
"""
def get_major_nt(nt_dic):
    nt = ''
    max_count = 0
    for tmp_nt, tmp_value in nt_dic.items():
        tmp_count, tmp_var_id = tmp_value
        if len(tmp_nt) == 1:
            assert tmp_nt in "ACGTDN"
        else:
            assert len(tmp_nt) == 2 \
                        and tmp_nt[0] == 'I' \
                        and tmp_nt[1] in "ACGT"
        if tmp_count > max_count:
            max_count = tmp_count
            nt = tmp_nt
    if len(nt) == 1:
        assert nt in "ACGTDN"
    else:
        assert len(nt) == 2 and nt[0] == 'I' and nt[1] in "ACGT"
    return nt                

"""
Match scores between two dictionaries of nucleotides containing ACGT
"""
def match_score(nt_dic1, nt_dic2):
    sum_1  = sum([count for count, _ in nt_dic1.values()])
    sum_2  = sum([count for count, _ in nt_dic2.values()])
    total1 = sum_1 * 2.0
    total2 = sum_2 * 2.0
    best  = 0.0
    for nt in "ACGT":
        if nt not in nt_dic1 or nt not in nt_dic2:
            continue
        tmp_best = nt_dic1[nt][0] / total1 + nt_dic2[nt][0] / total2
        if tmp_best > best:
            best = tmp_best
    return best

"""
Remove any deletions in the sequence. May need to change the D to another char
"""
def get_ungapped_seq(seq):
    ungapped_seq = []
    for i in range(len(seq)):
        nt_dic = seq[i]
        nt = get_major_nt(nt_dic)
        if nt == 'D': # TODO consider changing this D to . and/or N
            continue
        ungapped_seq.append(nt_dic)
    return ungapped_seq

"""
Convert coordinates of the sequnece to a position by removing indels
"""
def get_ungapped_seq_pos(seq, pos):
    tot_del_len = 0 
    tot_ins_len = 0
    for i in range(len(seq)):
        nt_dic = seq[i]
        nt     = get_major_nt(nt_dic)
        if nt == 'D':
            tot_del_len += 1
        elif nt[0] == 'I':
            tot_ins_len += 1
        if i - tot_ins_len == pos:
            return pos - tot_del_len
    return -1

"""
Get mate node id
    HSQ1008:141:D0CC8ACXX:3:2304:4780:36964|L to 
    HSQ1008:141:D0CC8ACXX:3:2304:4780:36964|R or vice versa
"""
def get_mate_node_id(node_id):
    node_id2, end = node_id.split('|')
    if end == 'L':
        end = 'R'
    else:
        end = 'L'
    node_id2 = '|'.join([node_id2, end])
    return node_id2

"""
Viterbi Algorithm for longest path through contig graph
"""
def viterbi_path(trellis, states, verbose = False):
    vit     = [[]]
    endpath = [-1,None]

    #initialize path
    node_score = -sys.maxsize
    for i in range(len(trellis[0])):
        weight = trellis[0][i]
        if weight > node_score:
            endpath    = [0,i]
            node_score = weight
        vit[0].append({"weight" : trellis[0][i], "prev" : None})
    
    #extend path
    for t in range(1, len(trellis)):
        vit.append([])
        node_score = -sys.maxsize
        for j in range(len(trellis[t])):
            (weight, state) = max([(vit[t-1][n]['weight'] \
                                        + trellis[t][j], n) \
                                    for n in range(len(vit[t-1]))], 
                                  key = lambda x : x[0])
            if weight > node_score:
                endpath    = [t,j]
                node_score = weight
            vit[t].append({"weight" : weight, "prev" : state})

    assert endpath[1] is not None

    # Backtrace Path
    path = []
    while endpath[1] is not None:
        t, node = endpath
        path.append(states[t][node])
        prev    = vit[t][node]["prev"]
        t      -= 1
        endpath = [t,prev]

    if verbose:
        print("States: %s" % str(states))
        print("Trellis: %s" % str(trellis))
        print("Viterbi: %s" % str(vit))

    return node_score, path[::-1]


class Node:
    # Initialize
    def __init__(self,
                 id,
                 left,
                 seq,
                 qual,
                 var,
                 ref_seq,
                 ref_vars,
                 mpileup,
                 simulation):
        self.next = [] # list of next nodes

        if simulation:
            id = id.split('_')[0]
        self.id   = id   # Node ID
        self.left = left # starting position

        # sequence that node represents
        #   with information about how the sequence is related to backbone
        assert len(seq) == len(var)
        assert len(seq) == len(qual)
        self.seq     = []
        self.ins_len = 0
        for s in range(len(seq)):
            nt = seq[s]
            if len(nt) == 1:
                assert nt in "ACGTDN"
            else:
                assert len(nt) == 2 and nt[0] == 'I' and nt[1] in "ACGT"
                self.ins_len += 1                
            var_id = var[s]
            self.seq.append({nt : [1, var_id]})
        self.qual = []
        for q in qual:
            if q != '':
                self.qual.append(max(0, ord(q) / 10 - 3))
            else:
                self.qual.append(0)

        self.right = self.left + len(seq) - 1 - self.ins_len

        self.read_ids = set([id])
        self.mate_ids = set([id.split('|')[0]])
        self.calculate_avg_cov()

        self.ref_seq  = ref_seq
        self.ref_vars = ref_vars
        self.mpileup  = mpileup

    # Check how compatible allele is in regard to read or pair
    def compatible_with_rnode(self, rnode):
        assert rnode.left + len(rnode.seq) <= len(self.seq)
        score = 0
        for i in range(len(rnode.seq)):
            allele_bp = self.seq[rnode.left + i]
            read_bp   = rnode.seq[i]
            if allele_bp == read_bp:
                score += 1

        return float(score) / len(rnode.seq)


    # Check how nodes overlap with each other without considering deletions
    def overlap_with(self, 
                     other, 
                     vars, 
                     skipN = False, 
                     debug = False):
        assert self.left <= other.left
        if self.right < other.left:
            return -1, -1
        seq       = get_ungapped_seq(self.seq)
        other_seq = get_ungapped_seq(other.seq)
        add_mm    = len(self.mate_ids & other.mate_ids)
        i_left    = get_ungapped_seq_pos(self.seq, other.left - self.left)
        for i in range(i_left - 5, i_left + 6):
            max_mm = 0.012 * (len(seq) - i) # 1 mismatch per 83 bases
            tmp_mm = 0.0
            for j in range(len(other_seq)):
                if i + j >= len(seq):
                    break
                nt_dic       = seq[i+j]
                other_nt_dic = other_seq[j]
                nt           = get_major_nt(nt_dic)
                other_nt     = get_major_nt(other_nt_dic)
                mismatch     = 0.0
                if skipN and (nt == 'N' or other_nt == 'N'):
                    mismatch = 0.0
                elif nt != other_nt:
                    mismatch = 1.0 - match_score(seq[i+j], other_seq[j])
                    
                    # Higher penalty for mismatches in variants
                    nt_var       = nt_dic[nt][1], 
                    other_nt_var = other_nt_dic[other_nt][1]
                    if nt_var != other_nt_var:
                        mismatch = 5.0
                        adjust = min(1.0, nt_dic[nt][0] \
                                            / self.get_avg_cov()) \
                                    * min(1.0, other_nt_dic[other_nt][0] \
                                                / other.get_avg_cov())
                        mismatch *= adjust
                        if mismatch < 1.0:
                            mismatch = 1.0

                assert mismatch >= 0.0
                tmp_mm += mismatch
                if tmp_mm > max_mm:
                    break

            if debug:
                print("at %d (%d) with overlap of %d and mismatch of %.2f" \
                        % (i, self.left + i, j, tmp_mm),
                      file=sys.stderr)

            if tmp_mm <= max_mm:
                return i, min(len(seq) - i, len(other_seq)), tmp_mm
                
        return -1, -1, sys.maxsize

    # Combine two nodes with considering deletions
    def combine_with(self, other):
        # DK - debugging purposes
        if self.left > other.left:
            self.print_info()
            other.print_info()
            return
        
        assert self.left <= other.left

        # Merge two sequences
        assert len(other.seq) > 0 and 'D' not in other.seq[0].keys()
        j = 0        
        # Merge the overlapped parts
        if self.right >= other.left:
            overlap = False 
            ins_len = 0
            for i in range(len(self.seq)):
                nt_dic = self.seq[i]
                nt     = get_major_nt(nt_dic)
                if nt.startswith('I'):
                    ins_len += 1
                if i == other.left - self.left + ins_len:
                    overlap = True
                    break
            assert overlap
            new_seq = self.seq[:i]
            while i < len(self.seq) and j < len(other.seq):
                nt_dic  = self.seq[i]
                nt_dic2 = other.seq[j]
                for nt, value in nt_dic2.items():
                    count, var_id = value
                    if nt in nt_dic:
                        nt_dic[nt][0] += count
                    else:
                        nt_dic[nt] = [count, var_id]
                new_seq.append(nt_dic)
                i += 1
                j += 1
            # this node contains the other node
            if i < len(self.seq):
                new_seq += self.seq[i:]
        # Fill in the gap between the two nodes if exists
        else:
            new_seq = self.seq[:]
            sum_1 = sum([count for count, _ in self.seq[-1].values()])
            sum_2 = sum([count for count, _ in other.seq[0].values()])
            flank_cov = (sum_1 + sum_2) / 2.0
            for k in range(other.left - self.right - 1):
                ref_nt_dic = self.mpileup[k + 1 + self.right][1]
                nt_dic = {}
                # Fill in the gap with Ns for now
                if len(ref_nt_dic) == 0 or True:
                    nt_dic = {'N' : [1, ""]}
                else:
                    weight = flank_cov \
                                / max(1.0, 
                                      sum([count \
                                            for count, _ in ref_nt_dic.values()]))
                    for nt, value in ref_nt_dic.items():
                        count, var_id = value
                        nt_dic[nt] = [count * weight, var_id]
                new_seq.append(nt_dic)

        # Append the rest of the other sequence to it
        if j < len(other.seq):
            new_seq += deepcopy(other.seq[j:])
        self.read_ids |= other.read_ids
        self.mate_ids |= other.mate_ids

        self.seq     = new_seq
        self.ins_len = 0
        for i in range(len(self.seq)):
            nt_dic = self.seq[i]
            nt     = get_major_nt(nt_dic)
            if nt[0] == 'I':
                self.ins_len += 1
        self.right = self.left + len(self.seq) - 1 - self.ins_len
        
        # Update coverage
        self.calculate_avg_cov()

    # Return the length of the ungapped sequence
    def ungapped_length(self):
        return len(get_ungapped_seq(self.seq))

    # Contains Ns?
    def contain_Ns(self):
        for i in range(len(self.seq)):
            nt_dic = self.seq[i]
            nt     = get_major_nt(nt_dic)
            if nt == 'N':
                return True
        return False

    # Get variant ids
    def get_var_ids(self, 
                    left = 0, 
                    right = sys.maxsize):
        vars    = []
        left    = max(left, self.left)
        right   = min(right, self.right)
        ins_len = 0
        for pos in range(left, right + 1):
            var_i = pos - self.left + ins_len
            while var_i < len(self.seq):
                nt_dic = self.seq[var_i]
                nt     = get_major_nt(nt_dic)
                if nt.startswith('I'):
                    var_i   += 1
                    ins_len += 1
                else:
                    break            
            for _, var in nt_dic.values():
                if var == "" or var == "unknown":
                    continue
                assert var in self.ref_vars
                if len(vars) > 0 and var == vars[-1]:
                    continue
                type, pos, data = self.ref_vars[var]
                if (type == "single" and data == nt) \
                        or (type == "deletion" and nt == 'D') \
                        or (type == "insertion" and len(nt) == 2 and nt[1] == data):
                    vars.append(var)
        return vars

    # Get variant ids
    #   left and right are gene-level coordinates
    def get_vars(self, 
                 left = 0, 
                 right = sys.maxsize):
        vars     = []
        left     = max(left, self.left)
        right    = min(right, self.right)
        skip_pos = -1
        ins_len  = 0
        for pos in range(left, right + 1):
            if pos <= skip_pos:
                continue
            var_i = pos - self.left + ins_len
            while var_i < len(self.seq):
                nt_dic = self.seq[var_i]
                nt     = get_major_nt(nt_dic)
                if nt.startswith('I'):
                    var_i   += 1
                    ins_len += 1
                    var      = nt_dic[nt][1]
                    if len(vars) > 0 and var != vars[-1][0]:
                        vars.append([var, pos])
                else:
                    break
            if nt == self.ref_seq[pos]:
                continue
            if nt == 'N':
                vars.append(["gap", pos])
                continue            
            added = False
            for _, var in nt_dic.values():
                if var == "" or var == "unknown":
                    continue
                if len(vars) > 0 and var == vars[-1][0]:
                    continue
                assert var in self.ref_vars
                type, var_pos, data = self.ref_vars[var]                    
                if data == nt or (type == "deletion" and nt == 'D'):
                    assert pos + ins_len >= var_pos
                    if type == "deletion" and pos > var_pos:
                        continue                    
                    if type == "deletion":
                        skip_pos = pos + int(data) - 1
                    added = True
                    vars.append([var, pos])
            if not added \
                    and "unknown" in [var_id for _, var_id in nt_dic.values()]:
                vars.append(["unknown", pos])

        return vars

    # Get average coverage
    def get_avg_cov(self):
        return self.avg

    # Calculate average coverage
    def calculate_avg_cov(self):
        self.avg = 0.0
        for nt_dic in self.seq:
            for count, _ in nt_dic.values():
                self.avg += count
        self.avg /= len(self.seq)
        return self.avg

    # Display node information
    def print_info(self, output=sys.stderr):
        seq      = "" 
        var_str  = ""
        prev_var = ""
        ins_len  = 0
        for i in range(len(self.seq)):
            if (self.left + i - ins_len) % 100 == 0:
                seq += ("|%d|" % (self.left + i - ins_len))
            elif (self.left + i - ins_len) % 20 == 0:
                seq += '|'
            nt_dic = self.seq[i]
            nt     = get_major_nt(nt_dic)
            if nt[0] == 'I':
                seq += "\033[93m"
            elif nt != self.ref_seq[self.left + i - ins_len]:
                var_id = nt_dic[nt][1]
                if var_id == "unknown" or var_id.startswith("nv"):
                    seq += "\033[91m" # red
                else:
                    seq += "\033[94m" # blue
            if nt[0] == 'I':
                seq += nt[1]
            else:
                seq += nt
            if nt[0] == 'I' or nt != self.ref_seq[self.left + i - ins_len]:
                seq += "\033[00m"

            var = []
            for _, var_id in nt_dic.values():
                if var_id == "":
                    continue
                var.append(var_id)
            var = '-'.join(var)
            if var != "" and var != prev_var:
                var_str += "\t%d: %s %s" \
                                % (self.left + i - ins_len, var, str(nt_dic))
            prev_var = var
            if nt[0] == 'I':
                ins_len += 1
        
        outlist = [
            "Node ID: %s" % self.id,
            "Pos: [%d, %d], Avg. coverage: %.1f" \
                % (self.left, self.right, self.get_avg_cov()),
            "\t%s" % seq,
            "\t%s" % var_str,
            "mates: %d" % len(self.mate_ids), # sorted(self.mate_ids)
            "reads: %d" % len(self.read_ids), # sorted(self.read_ids)
            "\n"         
        ]

        for out in outlist:
            print(out, file=output)

class Graph:
    def __init__(self,
                 backbone,
                 gene_vars,
                 exons,
                 primary_exons,
                 partial_allele_ids,
                 true_allele_nodes = {},
                 predicted_allele_nodes = {},
                 display_allele_nodes = {},
                 simulation = False):
        self.backbone               = backbone # backbone sequence
        self.gene_vars              = gene_vars
        self.exons                  = exons
        self.primary_exons          = primary_exons
        self.partial_allele_ids     = partial_allele_ids
        self.true_allele_nodes      = true_allele_nodes
        self.predicted_allele_nodes = predicted_allele_nodes
        self.allele_node_order      = []
        self.display_allele_nodes   = display_allele_nodes
        self.simulation             = simulation

        self.read_nodes = self.nodes = {}
        self.nodes2                  = None
        self.other_nodes             = {}
        self.edges                   = {}
        self.to_node                 = {}
        self.from_node               = {}

        self.left_margin     = 350
        self.right_margin    = 20
        self.top_margin      = 20
        self.bottom_margin   = 20
        self.scalex          = 5
        self.scaley          = 2
        self.width           = len(self.backbone) \
                                * self.scalex \
                                + self.left_margin \
                                + self.right_margin
        self.unscaled_height = 6000
        self.height          = self.unscaled_height \
                                * self.scaley
        self.coverage        = {}

    # Add node, which is an alignment w.r.t. the reference
    def add_node(self, 
                 id, 
                 id_i, 
                 node, 
                 simulation = False):
        if simulation:
            id = id.split('_')[0]
            
        if id_i == 0:
            if id in self.nodes:
                print("Warning multi-mapped read: %s" % id, 
                      file=sys.stderr)
                return
            assert id not in self.nodes
            self.nodes[id] = node
        else:
            if id not in self.other_nodes:
                self.other_nodes[id] = []
            self.other_nodes[id].append(node)

    # Remove nodes that are inside other nodes or with low coverage
    def remove_nodes(self, nodes):
        delete_ids = set()
        node_list  = [[id, node.left, node.right] for id, node in nodes.items()]
        def node_cmp(a, b):
            if a[2] != b[2]:
                return a[2] - b[2]
            else:
                return a[1] - b[1]
        node_list = sorted(node_list, cmp=node_cmp)
        for n in range(len(node_list)):
            id, left, right = node_list[n]
            node = nodes[id]
            i    = n - 1
            while i >= 0:
                id2, left2, right2 = node_list[i]
                if right2 < left:
                    break
                node2 = nodes[id2]
                if left <= left2 and right2 <= right:
                    at, overlap, mm = node.overlap_with(node2, self.gene_vars)
                    if mm < 1.0:
                        mult = overlap / float(max(right - left, right2 - left2))
                        if node2.get_avg_cov() * mult * 10 < node.get_avg_cov():
                            delete_ids.add(id2)
                        elif left == left2 and right == right2:
                            delete_ids.add(id)
                    elif overlap > 0:
                        if node2.get_avg_cov() * 10 < node.get_avg_cov():
                            delete_ids.add(id2)
                        elif node.get_avg_cov() * 10 < node2.get_avg_cov():
                            delete_ids.add(id)
                i -= 1

        for delete_id in delete_ids:
            del nodes[delete_id]
       
    """
    Full De Bruijn assembly
    """
    def guided_DeBruijn(self,
                        print_msg = False):
        assert len(self.nodes) > 0
        k          = 60 # k-mer
        DRB1_debug = False
        node_seq   = {}
        def add_node_seq(node_seq, id):
            nodes = [self.nodes[id]]
            if id in self.other_nodes:
                nodes += self.other_nodes[id]
            for node_i in range(len(nodes)):
                node = nodes[node_i]
                s    = 0 
                seq  = []
                while s < len(node.seq):
                    nt_dic = node.seq[s] # {'C': [1, '']}
                    nt     = get_major_nt(nt_dic)
                    if nt in "ACGTND":
                        seq.append(nt)
                    else:
                        assert len(nt) == 2 \
                                   and nt[0] == 'I' \
                                   and nt[1] in "ACGT"
                    s += 1

                if len(seq) < k:
                    continue

                def leftshift(seq, ref_seq):
                    seq_len = len(seq)
                    assert seq_len > 0 and seq[0] != 'D'

                    bp_i = 0
                    while bp_i < seq_len:
                        bp = seq[bp_i]
                        if bp != 'D':
                            bp_i += 1
                            continue
                        bp_j = bp_i + 1
                        while bp_j < seq_len:
                            bp2 = seq[bp_j]
                            if bp2 != 'D':
                                break
                            else:
                                bp_j += 1

                        if bp_j >= seq_len:
                            bp_i = bp_j
                            break

                        prev_i = bp_i 
                        prev_j = bp_j
                        while bp_i > 0 \
                                and seq[bp_i-1] in "ACGT" \
                                and ref_seq[bp_j-1] in "ACGT":
                            if seq[bp_i-1] != ref_seq[bp_j-1]:
                                break
                            seq[bp_j-1] = seq[bp_i-1]
                            seq[bp_i-1] = 'D'
                            bp_i       -= 1
                            bp_j       -= 1
                        bp_i = bp_j
                        while bp_i < seq_len:
                            if seq[bp_i] in "ACGT":
                                break
                            bp_i += 1

                if DRB1_debug:
                    leftshift(seq, self.backbone[node.left:node.left + len(seq)])
                node_seq["%s.%d" % (id, node_i)] = seq
            
        for id in self.nodes.keys():
            # Rough dictionary of node sequences indexed by node id
            add_node_seq(node_seq, id) 
            
        # AAA.1 => AAA, 1
        def get_id_and_sub(id):
            id_split = id.split('.')
            return '.'.join(id_split[:-1]), int(id_split[-1])

        # Build and clean de Bruijn Graph
        try_hard = False
        while True:
            delete_ids = set()
            nodes      = [] # list of nodes to contain kmers
            for id, node in self.nodes.items():
                nodes_ = [node]
                if id in self.other_nodes:
                    nodes_ += self.other_nodes[id]
                for node_i in range(len(nodes_)):
                    node = nodes_[node_i]
                    id_  = "%s.%d" % (id, node_i)
                    if id_ not in node_seq:
                        continue
                    seq = node_seq[id_]

                    if len(seq) < k or 'N' in seq:
                        continue
                    kmer = seq[:k] 
                    seq  = seq[k:]
                    nodes.append([id_, node.left, node.right, kmer, seq])
                
            def node_cmp(a, b):
                if a[1] != b[1]:
                    return a[1] - b[1]
                else:
                    return a[2] - b[2]
            nodes = sorted(nodes, cmp=node_cmp) # Sort by left position

            # Generate numerical read IDs (Positional Indecies)
            id_to_num = {}
            num_to_id = []
            for id in [node[0] for node in nodes]:
                id_to_num[id] = len(id_to_num)
                num_to_id.append(id)

            # Construct De Bruijn graph with 60-mer
            self.debruijn = debruijn = [[] \
                                for i in range(len(self.backbone) - k + 1)]
            min_n = 0
            for pos in range(len(debruijn)):
                for n in range(min_n, len(nodes)):
                    id, node_pos, node_right, kmer, seq = nodes[n]
                    if node_pos < pos:
                        min_n = n + 1
                        continue
                    elif node_pos > pos:
                        break

                    assert len(kmer) == k

                    # Add a new node or update the De Bruijn graph
                    curr_vertices = debruijn[pos]
                    found         = False
                    kmer_seq      = ''.join(kmer)
                    for v in range(len(curr_vertices)):
                        cmp_nt, cmp_k_m1_mer = curr_vertices[v][:2]
                        if kmer_seq == cmp_k_m1_mer + cmp_nt:                        
                            curr_vertices[v][3].append(n)
                            found = True
                            break

                    if not found:
                        predecessors = []
                        if pos > 0:
                            prev_vertices = debruijn[pos - 1]
                            for v in range(len(prev_vertices)):
                                cmp_nt, cmp_k_m1_mer = prev_vertices[v][:2]
                                if kmer_seq[:-1] == cmp_k_m1_mer[1:] + cmp_nt:
                                    predecessors.append(v)
                        debruijn[pos].append([kmer_seq[-1],           # base
                                              ''.join(kmer_seq[:-1]), # (k-1)-mer
                                              predecessors,           # predecessors
                                              [n]])                   # int read IDs

                    # Update k-mer
                    if len(seq) > 0:
                        kmer = kmer[1:] + seq[:1]
                        seq  = seq[1:]
                        nodes[n] = [id, node_pos + 1, node_right, kmer, seq]

            # Average number of kmers
            total_kmers = 0
            for pos in range(len(debruijn)):
                vertices = debruijn[pos]
                for _, _, _, num_ids in vertices:
                    total_kmers += len(num_ids)
            avg_kmers = float(total_kmers) / len(debruijn)

            # Filter out reads
            for pos in range(len(debruijn)):
                vertices     = debruijn[pos]
                num_vertices = 0
                num_kmers    = 0
                # TODO Check to see if this DRB1 handling is needed
                for v in range(len(vertices)):
                    _, _, predecessors, num_ids = vertices[v]
                    if not (set(num_ids) <= delete_ids):
                        num_vertices += 1
                        if DRB1_debug:
                            num_kmers = len(set(num_ids) - delete_ids)
                if num_vertices <= 1:
                    if DRB1_debug:
                        if pos > 300 and pos + 300 < len(debruijn):
                            if num_vertices == 1 and num_kmers * 8 < avg_kmers:
                                for _, _, _, num_ids in vertices:
                                    delete_ids |= set(num_ids)
                    continue
                
                vertice_count = [0] * len(vertices)
                for v in range(len(vertices)): # Iter through each position vertex
                    _, _, predecessors, num_ids = vertices[v]
                    for num_id in num_ids: # Iter through each node id
                        if num_id in delete_ids:
                            continue
                        read_id = get_id_and_sub(num_to_id[num_id])[0]
                        if read_id in self.other_nodes:
                            continue
                        mate_read_id = get_mate_node_id(read_id)
                        if mate_read_id in self.nodes:
                            vertice_count[v] += 1

                # First look at and remove reads that are multi-aligned locally
                first_pair = None
                for v in range(len(vertices)):
                    read_ids = set([get_id_and_sub(num_to_id[num_id])[0] \
                                        for num_id in vertices[v][3]])
                    for v2 in range(v + 1, len(vertices)):
                        read_ids2 = set([get_id_and_sub(num_to_id[num_id])[0] \
                                        for num_id in vertices[v2][3]])
                        if read_ids & read_ids2:
                            first_pair = [v, v2, read_ids & read_ids2]
                            break

                debug_msg = False
                if debug_msg:
                    print(("at", pos, vertices), 
                           file=sys.stderr)
                    print(("count:", vertice_count), 
                            file=sys.stderr)

                if try_hard:
                    vertice_with_id = [[vertice_count[v], v] \
                                            for v in range(len(vertice_count))]
                    vertice_with_id = sorted(vertice_with_id, key=lambda a: a[0])
                    for v in range(len(vertice_count) - 2):
                        v           = vertice_with_id[v][1]
                        num_ids     = vertices[v][3]
                        delete_ids |= set(num_ids)
                        if debug_msg:
                            print((v, "is removed with", num_ids), 
                                  file=sys.stderr)
                else:
                    if first_pair:
                        v, v2, multi_read_ids = first_pair
                        v_ = v if vertice_count[v] < vertice_count[v2] else v2
                        for num_id in vertices[v_][3]:
                            id = get_id_and_sub(num_to_id[num_id])[0]
                            if id in multi_read_ids:
                                delete_ids.add(num_id)
                    else:
                        assert len(vertices) >= 2
                        relative_avg = (sum(vertice_count) - vertice_count[v]) \
                                            / float(len(vertice_count) - 1)
                        if len(vertices) == 2:
                            for v in range(len(vertices)):
                                # Eliminate reads that have conflicts with other reads due to a deletion
                                if vertice_count[v] * 2 < relative_avg:
                                    nt, kmer, _, num_ids = vertices[1-v]
                                    if nt == 'D':
                                        num_id  = num_ids[0]
                                        id_sub  = num_to_id[num_id]
                                        id, sub = get_id_and_sub(id_sub)
                                        if sub == 0:
                                            left = pos - self.nodes[id].left
                                        else:
                                            left = pos - self.other_nodes[id][sub - 1].left
                                        seq       = node_seq[id_sub]
                                        seq_right = ''.join(seq[left+k:])
                                        seq_right = seq_right.replace('D', '')
                                        success   = True
                                        for num_id2 in vertices[v][3]:
                                            id_sub2 = num_to_id[num_id2]
                                            id2, sub2 = get_id_and_sub(id_sub2)
                                            if sub2 == 0:
                                                left2 \
                                                    = pos \
                                                      - self.nodes[id2].left
                                            else:
                                                left2 \
                                                    = pos \
                                                      - self.other_nodes\
                                                          [id2][sub2 - 1].left
                                            seq2       = node_seq[id_sub2]
                                            seq2_right = ''.join(seq2[left2+k:])
                                            if seq_right.find(seq2_right) != 0:
                                                success = False
                                                break
                                        if success:
                                            delete_ids |= set(vertices[v][3])

                                # DK - working on ...
                                if DRB1_debug:
                                    if vertice_count[v] * 8 < relative_avg:
                                        num_ids     = vertices[v][3]
                                        delete_ids |= set(num_ids)
                                        if debug_msg:
                                            print((v, "is removed with", num_ids), 
                                                  file=sys.stderr)
                                    elif vertice_count[v] * 8 < avg_kmers:
                                        num_ids     = vertices[v][3]
                                        delete_ids |= set(num_ids)
                        else:
                            second2last = sorted(vertice_count)[1]
                            for v in range(len(vertices)):
                                if vertice_count[v] < second2last:
                                    num_ids     = vertices[v][3]
                                    delete_ids |= set(num_ids)
                                    if debug_msg:
                                        print((v, "is removed with", num_ids), 
                                              file=sys.stderr)
                if debug_msg:
                    print("\n\n", file=sys.stderr)      
                
            # delete nodes based on delete_ids
            ids_to_be_updated = set()
            for num_id in delete_ids:
                id_sub  = num_to_id[num_id]
                id, sub = get_id_and_sub(id_sub)
                ids_to_be_updated.add(id)
                if sub == 0:
                    self.nodes[id] = None
                else:
                    self.other_nodes[id][sub-1] = None

            # Cleaning up other nodes
            for id in self.nodes.keys():
                other_nodes = []
                if id in self.other_nodes:
                    for other_node in self.other_nodes[id]:
                        if other_node != None:
                            other_nodes.append(other_node)
                if self.nodes[id] == None:
                    if len(other_nodes) == 0:
                        del self.nodes[id]
                    else:
                        self.nodes[id] = other_nodes[0]
                        del other_nodes[0]
                if id in self.other_nodes:
                    if len(other_nodes) == 0:
                        del self.other_nodes[id]
                    else:
                        self.other_nodes[id] = other_nodes

            for id in ids_to_be_updated:
                if id in self.nodes:
                    add_node_seq(node_seq, id) # Update node sequences

            if len(delete_ids) == 0:
                if try_hard:
                    break
                else:
                    try_hard = True

        # Print De Bruijn graph
        if print_msg:
            for i in range(len(debruijn)):
                curr_vertices = debruijn[i]
                if len(curr_vertices) == 0:
                    continue
                consensus_seq = [{} for j in range(k)]
                for v in range(len(curr_vertices)):
                    nt, k_m1_mer = curr_vertices[v][:2]
                    kmer = k_m1_mer + nt
                    assert len(kmer) == k
                    for j in range(k):
                        nt = kmer[j]
                        if nt not in consensus_seq[j]:
                            consensus_seq[j][nt] = 1
                        else:
                            consensus_seq[j][nt] += 1

                print(i, file=sys.stderr)
                for v in range(len(curr_vertices)):
                    nt, k_m1_mer, predecessors, num_ids = curr_vertices[v]
                    kmer     = k_m1_mer + nt
                    kmer_seq = ""
                    for j in range(k):
                        nt = kmer[j]
                        if len(consensus_seq[j]) >= 2:
                            kmer_seq += "\033[94m"
                        kmer_seq += nt
                        if len(consensus_seq[j]) >= 2:
                            kmer_seq += "\033[00m"
                        
                    print(("\t%d:" % v, 
                           kmer_seq, 
                           len(num_ids), 
                           predecessors, 
                           num_ids),
                          file=sys.stderr)

        id_to_num = {} # Regenerate new ID_to_num
        for num in range(len(num_to_id)):
            id_sub         = num_to_id[num]
            id             = get_id_and_sub(id_sub)[0]
            num_to_id[num] = id
            if id not in id_to_num:
                id_to_num[id] = set()
            id_to_num[id].add(num)          
                    
        ### Compress read nodes
        paths      = []
        path_queue = deque() 
        done       = set()
        for i in range(len(debruijn)):
            if len(debruijn[i]) == 0:
                continue
            for i2 in range(len(debruijn[i])):
                path_queue.append("%d-%d" % (i, i2))
            break

        while len(path_queue) > 0:
            i_str = path_queue.popleft()
            if i_str in done:
                continue

            i, i2   = i_str.split('-')
            i, i2   = int(i), int(i2)
            num_ids = debruijn[i][i2][3]
            j       = i + 1
            while j < len(debruijn):
                merge       = len(debruijn[j-1]) > len(debruijn[j])
                branch      = len(debruijn[j-1]) < len(debruijn[j])
                new_i2      = -1
                tmp_num_ids = []
                found       = False
                for j2 in range(len(debruijn[j])):
                    _, _, predecessors, add_read_ids = debruijn[j][j2]
                    if len(predecessors) == 0:
                        branch = True
                        path_queue.append("%d-%d" % (j, j2))
                    elif i2 in predecessors:
                        found = True
                        # merge into one node
                        if len(predecessors) > 1:
                            merge = True
                        if new_i2 >= 0:
                            branch = True
                        new_i2 = j2
                        tmp_num_ids += add_read_ids

                if merge or branch:
                    for j2 in range(len(debruijn[j])):
                        _, _, predecessors, add_num_ids = debruijn[j][j2]
                        if i2 in predecessors:
                            path_queue.append("%d-%d" % (j, j2))
                    break
                if not found:
                    break
                
                num_ids += tmp_num_ids
                i2       = new_i2
                j       += 1

            done.add(i_str)
            num_ids = set(num_ids)
            paths.append([i, j, num_ids])
            if j < len(debruijn) and len(debruijn[j]) == 0:
                j += 1
                while j < len(debruijn) and len(debruijn[j]) == 0:
                    j += 1
                if j < len(debruijn):
                    for j2 in range(len(debruijn[j])):
                        path_queue.append("%d-%d" % (j, j2))
                        
        def get_mate_num_ids(num_ids):
            mate_num_ids = set()
            for num_id in num_ids:
                read_id = num_to_id[num_id]
                mate_read_id = get_mate_node_id(read_id)
                if mate_read_id in id_to_num:
                    mate_num_id   = id_to_num[mate_read_id]
                    mate_num_ids |= mate_num_id
                    
            return mate_num_ids
        
        ### Generate a compressed assembly graph
        def path_cmp(a, b):
            if a[0] != b[0]:
                return a[0] - b[0]
            else:
                return a[1] - b[1]
        paths = sorted(paths, cmp=path_cmp)

        for p in range(len(paths)):
            if print_msg: 
                print(("path:", p, paths[p]), 
                      file=sys.stderr)

        # Build list of equivalent nodes
        excl_num_ids = set() # exclusive num ids
        equiv_list   = []
        p            = 0
        while p < len(paths):
            left, right, num_ids = paths[p]
            p2 = p + 1
            while p2 < len(paths):
                next_left, next_right, next_num_ids = paths[p2]
                if next_left >= right:
                    break
                p2 += 1

            equiv_list.append([])
            for i in range(p, p2):
                left, right, num_ids = paths[i]
                equiv_list[-1].append([[i], 
                                       num_ids, 
                                       num_ids | get_mate_num_ids(num_ids), 
                                       set()])
                if p + 1 < p2:
                    assert p + 2 == p2
                    excl_num_ids |= num_ids
            p = p2

        new_equiv_list = []
        for classes in equiv_list:
            if len(classes) > 1:
                new_equiv_list.append(classes)
                continue
            assert len(classes) == 1
            num_ids = classes[0][1] - excl_num_ids
            if len(num_ids) <= 0:
                continue
            classes[0][1] = num_ids
            classes[0][2] = num_ids | get_mate_num_ids(num_ids)
            new_equiv_list.append(classes)
        equiv_list = new_equiv_list

        # Compress and phase de Bruijn to allele
        known_alleles = False
        while True:
            if print_msg:
                for i in range(len(equiv_list)):
                    classes = equiv_list[i]
                    for j in range(len(classes)):
                        ids, num_ids, all_ids, alleles = classes[j]
                        print((i, 
                               j, 
                               ids, 
                               len(num_ids), 
                               sorted(list(num_ids))[:20], 
                               alleles), 
                              file=sys.stderr)
                    print("\n", file=sys.stderr)

            # Collapse nodes to reference (annotation)
            def annotate_contig(viterbi = False):
                # Phasing Algorithm based on graph coloring and the principles behind viturbi algorithm
                if viterbi:
                    def jaccard(setA, setB):
                        setA  = set(setA) 
                        setB  = set(setB)
                        inter = len(setA & setB)+1
                        union = len(setA | setB)+1
                        return math.log10(float(inter)/float(union))

                    # TODO: Work in progress building maximixing path
                    alleles = self.predicted_allele_nodes.keys()
                    vitres  = {'key' : [], 'value' : [], 'path' : []}
                    anodes  = [None, None]
                    # For all allele pairs ...
                    for i in range(len(alleles)):
                        anodes[0] = self.predicted_allele_nodes[alleles[i]]
                        for j in range(i,len(alleles)):
                            vitres['key'].append([alleles[i], alleles[j]])                            
                            anodes[1] = self.predicted_allele_nodes[alleles[j]]
                            trellis   = [] 
                            states    = []

                            # ... go through length of all 
                            # contigs across allele ...
                            node_IDs = []
                            for k in range(len(equiv_list)):
                                classes = equiv_list[k]
                                mx = []
                                node_IDs.append([])

                                # ... and for each contig at a position ...
                                for l in range(len(classes)):
                                    mx.append([])
                                    num_id    = sorted(list(classes[l][1]))[0]
                                    node_id   = "(%d-%d)%s" \
                                                    % (k, l, num_to_id[num_id])
                                    node      = self.nodes2[node_id]
                                    node_vars = node.get_var_ids()
                                    node_IDs[-1].append(node_id)
                                    # ... compair the contig to the allele region the contig aligns to
                                    for m in range(len(anodes)):
                                        allele_vars = anodes[m].get_var_ids(
                                                        node.left, 
                                                        node.right
                                                      )
                                        mx[-1].append(jaccard(node_vars, 
                                                              allele_vars))

                                # The add the comparison to a trellis graph
                                assert mx
                                if len(mx) > 1:
                                    state = [[0,1],[1,0]]
                                    mx[1] = mx[1][::-1]
                                    # Find best value for pair
                                    mx    = [sum(k) for k in zip(*mx)]
                                else:
                                    state = [[0,0], [0,0]]
                                    mx = mx[0]
                                states.append(state)
                                trellis.append(mx)

                            score, path = viterbi_path(trellis, states)
                            vitres['path'].append(path)
                            vitres['value'].append(score)

                    print(vitres)
                    ix = max(range(len(vitres['value'])), 
                             key=vitres['value'].__getitem__)
                    best_alleles = vitres['key'][ix]
                    best_path    = vitres['path'][ix]
                    best_score   = vitres['value'][ix]

                    assert len(best_path) == len(equiv_list)
                    for i in range(len(equiv_list)):
                        classes = equiv_list[i]
                        for j in range(len(best_path[i])):
                            if best_alleles[j] not in classes[best_path[i][j]][3]:
                                classes[best_path[i][j]][3].add(best_alleles[j])

                    return [best_alleles, best_score]

                # Alternative phasing algorithm
                else:
                    for i in range(len(equiv_list)):
                        classes = equiv_list[i]
                        for j in range(len(classes)):
                            num_ids     = sorted(list(classes[j][1]))
                            node_id     = "(%d-%d)%s" \
                                             % (i, j, num_to_id[num_ids[0]])
                            node        = self.nodes2[node_id]
                            node_vars   = node.get_var_ids()
                            max_alleles = set() 
                            max_common  = -sys.maxsize
                            for anode in self.predicted_allele_nodes.values():
                                allele_vars = anode.get_var_ids(node.left, 
                                                                node.right)
                                tmp_common = len(set(node_vars) \
                                                    & set(allele_vars)) \
                                                - len(set(node_vars) \
                                                    | set(allele_vars))
                                if tmp_common > max_common:
                                    max_common  = tmp_common
                                    max_alleles = set([anode.id])
                                elif tmp_common == max_common:
                                    max_alleles.add(anode.id)
                            classes[j][3] = max_alleles

            if known_alleles:
                v_coloring = annotate_contig(True) # Use Viturbi coloring

            # Identify identical stretches of nodes to merge
            best_common_mat = [], 
            best_stat       = -sys.maxsize, 
            best_i          = -1 
            best_i2         = -1
            for i in range(len(equiv_list) - 1):
                classes = equiv_list[i]
                for i2 in range(i + 1, len(equiv_list)):
                    classes2   = equiv_list[i2]
                    common_mat = []
                    for j in range(len(classes)):
                        common_mat.append([])
                        if known_alleles:
                            ids = classes[j][3]
                        else:
                            ids = classes[j][2]
                        for j2 in range(len(classes2)):
                            if known_alleles:
                                ids2 = classes2[j2][3]
                            else:
                                ids2 = classes2[j2][2]
                            common_mat[-1].append(len(ids & ids2))

                    # Calculate stat
                    common_stat = 0
                    if len(classes) == 1 or len(classes2) == 1:
                        for row in common_mat:
                            common_stat += sum(row)
                    else:
                        for row in common_mat:
                            sorted_row   = sorted(row, reverse=True)
                            common_stat += (sorted_row[0] - sorted_row[1])
                        if common_mat[0][0] + common_mat[1][1] \
                                == common_mat[1][0] + common_mat[0][1]:
                            common_stat = -1

                    if common_stat > best_stat:
                        best_common_mat = common_mat
                        best_stat       = common_stat
                        best_i          = i
                        best_i2         = i2

            if print_msg:
                print(("best:", best_i, best_i2, best_stat, best_common_mat), 
                      file=sys.stderr)
                print("\n\n", file=sys.stderr)

            if known_alleles and best_stat < 0:
                self.remove_nodes(self.nodes2)
                break

            # Compress nodes: Collapse nodes to contigs 
            if best_stat < 0: 
                known_alleles = True
                new_nodes     = {}
                for i in range(len(equiv_list)):
                    classes = equiv_list[i]
                    for j in range(len(classes)):
                        ids, num_ids, all_ids, alleles = classes[j]
                        num_ids = sorted(list(num_ids))

                        if print_msg: 
                            print((i, j, num_ids), 
                                  file=sys.stderr)

                        assert (num_ids) > 0
                        read_id = num_to_id[num_ids[0]]
                        node    = deepcopy(self.nodes[read_id])
                        for num_id2 in num_ids[1:]:
                            read_id2 = num_to_id[num_id2]
                            node2    = self.nodes[read_id2]
                            node.combine_with(node2)

                        new_read_id = "(%d-%d)%s" % (i, j, read_id)
                        node.id     = new_read_id
                        assert new_read_id not in new_nodes
                        new_nodes[new_read_id] = node
                        
                self.nodes  = new_nodes                
                self.nodes2 = deepcopy(self.nodes)
                self.remove_nodes(self.nodes)
                continue

            mat      = best_common_mat
            classes  = equiv_list[best_i]
            classes2 = equiv_list[best_i2]

            # Filter vertices further if necessary
            def del_row(classes, mat, r):
                return classes[:r] + classes[r+1:], mat[:r] + mat[r+1:]
            
            def del_col(classes, mat, c):                    
                new_mat = []
                for row in mat:
                    row = row[:c] + row[c+1:]
                    new_mat.append(row)
                return classes[:c] + classes[c+1:], new_mat
                
            ## Forcing Collapse to at most two alleles
            assert len(classes) <= 2 and len(classes2) <= 2
            if len(classes) == 2 and len(classes2) == 2:
                # Check row
                num_ids1 = len(classes[0][1])
                num_ids2 = len(classes[1][1])
                if num_ids1 * 6 < num_ids2 or num_ids2 * 6 < num_ids1:
                    row_sum1 = sum(mat[0])
                    row_sum2 = sum(mat[1])
                    if row_sum1 > max(2, row_sum2 * 6):
                        classes, mat   = del_row(classes, mat, 1)
                        classes[0][1] -= excl_num_ids
                    elif row_sum2 > max(2, row_sum1 * 6):
                        classes, mat   = del_row(classes, mat, 0)
                        classes[0][1] -= excl_num_ids
                # Check column
                if len(classes) == 2:
                    num_ids1 = len(classes2[0][1])
                    num_ids2 = len(classes2[1][1])
                    if num_ids1 * 6 < num_ids2 or num_ids2 * 6 < num_ids1:
                        col_sum1 = mat[0][0] + mat[1][0] 
                        col_sum2 = mat[0][1] + mat[1][1]
                        if col_sum1 > max(2, col_sum2 * 6):
                            classes2, mat   = del_col(classes2, mat, 1)
                            classes2[0][1] -= excl_num_ids
                        elif col_sum2 > max(2, col_sum1 * 6):
                            classes2, mat   = del_col(classes2, mat, 0)
                            classes2[0][1] -= excl_num_ids

            merge_list = []
            def add_merge(classes, 
                          classes2, 
                          i, 
                          j, 
                          k):
                if known_alleles:
                    num_id1 = sorted(list(classes[i][1]))[0]
                    num_id2 = sorted(list(classes2[j][1]))[0]

                    node_id1 = "(%d-%d)%s" % (best_i, i, num_to_id[num_id1])
                    node_id2 = "(%d-%d)%s" % (best_i2, j, num_to_id[num_id2])
                    node_id3 = "(%d-%d)%s" % (best_i, k, num_to_id[min(num_id1, 
                                                                       num_id2)])
                    merge_list.append([node_id1, node_id2, node_id3])

                classes[i][0]  = sorted(classes[i][0] + classes2[j][0])
                classes[i][1] |= classes2[j][1]

            copy_list = []
            def add_copy(classes, 
                         classes2, 
                         i, 
                         j, 
                         k):
                if known_alleles:
                    num_ids  = classes2[j][1]
                    num_ids  = sorted(list(num_ids))
                    num_id   = num_ids[0]
                    node_id  = "(%d-%d)%s" % (best_i2, j, num_to_id[num_id])
                    node_id2 = "(%d-%d)%s" % (best_i, k, num_to_id[num_id])
                    copy_list.append([node_id, node_id2])

                classes[i] = classes2[j]

            remove_list = []
            def add_remove(classes, i):
                if known_alleles:
                    num_ids = classes[i][1]
                    num_ids = sorted(list(num_ids))
                    num_id  = num_ids[0]
                    node_id = "(%d-%d)%s" % (best_i, i, num_to_id[num_id])
                    remove_list.append([node_id])
                classes = [classes[1-i]]
                         
            if len(classes) == 1 and len(classes2) == 1:
                add_merge(classes, classes2, 0, 0, 0)
                
            elif len(classes) == 1:
                if 0 not in classes[0][0] \
                        and mat[0][0] > max(2, mat[0][1] * 6) \
                        and len(classes2[0][1]) > len(classes2[1][1]) * 2:
                    add_merge(classes, classes2, 0, 0, 0)
                elif 0 not in classes[0][0] \
                        and mat[0][1] > max(2, mat[0][0] * 6) \
                        and len(classes2[1][1]) > len(classes2[0][1]) * 2:
                    add_merge(classes, classes2, 0, 1, 0)
                else:
                    classes.append(deepcopy(classes[0]))

                    # Handle a special case at 5' end
                    if 0 in classes[0][0] \
                            and len(classes[0][0]) == 1 \
                            and (mat[0][0] > mat[0][1] * 2 \
                                or mat[0][1] > mat[0][0] * 2):
                        if mat[0][0] > mat[0][1]:
                            add_merge(classes, classes2, 0, 0, 0)
                            add_copy(classes, classes2, 1, 1, 1)
                        else:
                            assert mat[0][1] > mat[0][0]
                            add_copy(classes, classes2, 0, 0, 0)
                            add_merge(classes, classes2, 1, 1, 1)
                    else:
                        add_merge(classes, classes2, 0, 0, 0)
                        add_merge(classes, classes2, 1, 1, 1)
                        
            elif len(classes2) == 1:
                if mat[0][0] > max(2, mat[1][0] * 6):
                    add_merge(classes, classes2, 0, 0, 0)
                    if len(classes[0][1]) > len(classes[1][1]) * 6:
                        add_remove(classes, 1)
                elif mat[1][0] > max(2, mat[0][0] * 6):
                    add_merge(classes, classes2, 1, 0, 0)
                    if len(classes[1][1]) > len(classes[0][1]) * 6:
                        add_remove(classes, 0)
                else:
                    add_merge(classes, classes2, 0, 0, 0)
                    add_merge(classes, classes2, 1, 0, 1)
                    
            else:                
                score00 = mat[0][0] + mat[1][1]
                score01 = mat[0][1] + mat[1][0]
                if score00 > score01:
                    add_merge(classes, classes2, 0, 0, 0)
                    add_merge(classes, classes2, 1, 1, 1)
                elif score00 < score01:
                    add_merge(classes, classes2, 0, 1, 0)
                    add_merge(classes, classes2, 1, 0, 1)
                else:
                    break

            for c in range(len(classes)):
                classes[c][2] = classes[c][1] | get_mate_num_ids(classes[c][1])

            equiv_list[best_i] = classes            
            equiv_list         = equiv_list[:best_i2] + equiv_list[best_i2+1:]
            
            # Phasing to allele
            if known_alleles:
                exclude_ids = set()
                new_nodes   = {}
                for node_id1, node_id2, node_id3 in merge_list:
                    if self.nodes2[node_id1].left <= self.nodes2[node_id2].left:
                        node  = deepcopy(self.nodes2[node_id1])
                        node2 = self.nodes2[node_id2]
                    else:                        
                        node  = deepcopy(self.nodes2[node_id2])
                        node2 = self.nodes2[node_id1]
                    node.combine_with(node2)
                    node.id = node_id3
                    new_nodes[node_id3] = node
                    exclude_ids.add(node_id1)
                    exclude_ids.add(node_id2)

                for node_id1, node_id2 in copy_list:
                    node    = self.nodes2[node_id1]
                    node.id = node_id2
                    new_nodes[node_id2] = node
                    exclude_ids.add(node_id1)

                exclude_ids |= set(remove_list)

                for node_id, node in self.nodes2.items():
                    if node_id in exclude_ids:
                        continue
                    num, id = node_id.split(')')
                    i, i2   = num[1:].split('-')
                    i, i2   = int(i), int(i2)
                    if i > best_i2:
                        i -= 1
                    node_id = "(%d-%d)%s" % (i, i2, id)
                    node.id = node_id
                    new_nodes[node_id] = node
                        
                self.nodes2 = new_nodes
        
        if known_alleles:
            return v_coloring
        else:
            return [['No Known alleles to match'], -1]

    # Display graph information
    def print_info(self): 
        print("Backbone len: %d" % len(self.backbone), 
              file=sys.stderr)
        print("\t%s" % self.backbone, 
              file=sys.stderr)   

    # Compare nodes and get information
    def get_node_comparison_info(self, node_dic):
        assert len(node_dic) > 0
        nodes = [[id, node.left, node.right] for id, node in node_dic.items()]
        def node_cmp(a, b):
            if a[1] != b[1]:
                return a[1] - b[1]
            else:
                return a[2] - b[2]
        nodes  = sorted(nodes, cmp=node_cmp)
        seqs   = [] 
        colors = []
        for p in range(len(self.backbone)):
            nts = set()
            for n in range(len(nodes)):
                id, left, right = nodes[n]
                node = node_dic[id]
                if p >= left and p <= right:
                    nt_dic = node.seq[p - left]
                    nt     = get_major_nt(nt_dic)
                    nts.add(nt)

            for n in range(len(nodes)):
                if p == 0:
                    seqs.append([])
                    colors.append([])
                id, left, right = nodes[n]
                node = node_dic[id]
                if p >= left and p <= right:
                    nt_dic = node.seq[p - left]
                    nt     = get_major_nt(nt_dic)
                    seqs[n].append(nt)
                    if nt != self.backbone[p]:
                        if len(nts) > 1:
                            colors[n].append('R')
                        else:
                            colors[n].append('B')
                    else:
                        colors[n].append('N')
                else:
                    seqs[n].append(' ')

        assert len(nodes) == len(seqs)
        for n in range(len(nodes)):
            node      = nodes[n] 
            seq       = seqs[n]
            color     = colors[n]
            new_left  = 0 
            new_right = len(seq) - 1
            while seq[new_left] == 'D':
                new_left += 1
            while seq[new_right] == 'D':
                new_right -= 1

            node[1]   = new_left
            node[2]   = new_right
            seqs[n]   = seq[new_left:new_right+1]
            colors[n] = color[new_left:new_right+1]

        return nodes, seqs, colors

    # Compare nodes
    def print_node_comparison(self, node_dic):
        nodes, seqs, colors = self.get_node_comparison_info(node_dic)
        interval = 100
        for p in range(0, (len(self.backbone) + interval - 1) \
                                / interval * interval, interval):
            cur_seqs = []
            for n in range(len(nodes)):
                id, left, right = nodes[n] # inclusive coordinate
                right    += 1
                seq       = []
                seq_left  = max(p, left)
                seq_right = min(p+interval, right)
                if seq_left >= seq_right:
                    continue
                if p < left:
                    seq += ([' '] * (left - p))
                for s in range(seq_left, seq_right):
                    nt    = seqs[n][s-left]
                    color = colors[n][s-left]
                    if color in "RB":
                        if color == 'R':
                            nt = "\033[91m" + nt
                        else:
                            nt = "\033[94m" + nt
                        nt += "\033[00m"        
                    seq.append(nt)
                if right < p + interval:
                    seq += ([' '] * (p + interval - right))
                seq = ''.join(seq)
                cur_seqs.append([seq, id])

            if len(cur_seqs) <= 0:
                continue
                
            print(p, file=sys.stderr)
            for seq, id in cur_seqs:
                print(("\t", seq, id), 
                      file=sys.stderr)

    # Calculate coverage
    def calculate_coverage(self):
        if self.simulation:
            allele_nodes = self.true_allele_nodes
        else:
            allele_nodes = self.predicted_allele_nodes
        allele_nodes = [[id, node.left, node.right] \
                            for id, node in allele_nodes.items()]
        coverage = {}
        for allele_id, _, _ in allele_nodes:
            coverage[allele_id] = [0.0 for _ in range(len(self.backbone))]

        nodes = [[id, node.left, node.right] for id, node in self.nodes.items()]
        for id, left, right in nodes:
            node   = self.nodes[id]
            nodes2 = [[node, left, right]]
            if id in self.other_nodes:
                for node in self.other_nodes[id]:
                    nodes2.append([node, node.left, node.right])

            for node, left, right in nodes2:
                node_vars           = node.get_vars()
                node_var_ids        = node.get_var_ids()
                max_common          = -sys.maxsize
                max_allele_node_ids = []
                for allele_node_id, allele_left, allele_right in allele_nodes:
                    if right - left <= 500 and (left < allele_left or right > allele_right):
                        continue
                    if self.simulation:
                        allele_node = self.true_allele_nodes[allele_node_id]
                    else:
                        allele_node = self.predicted_allele_nodes[allele_node_id]
                    allele_vars = allele_node.get_var_ids(left, right)
                    common_vars = set(node_var_ids) & set(allele_vars)
                    tmp_common  = len(common_vars) \
                                    - len(set(node_var_ids) | set(allele_vars))
                    if max_common < tmp_common:
                        max_common          = tmp_common
                        max_allele_node_ids = [allele_node_id]
                    elif max_common == tmp_common:
                        max_allele_node_ids.append(allele_node_id)
                if len(max_allele_node_ids) <= 0:
                    continue
                add_cov = 1.0 / len(nodes2) / len(max_allele_node_ids)
                assert add_cov > 0.0
                for allele_node_id in max_allele_node_ids:
                    for p in range(left, right + 1):
                        coverage[allele_node_id][p] += add_cov

        max_cov = 0.0
        for allele_id, cov in coverage.items():
            max_cov = max(max_cov, max(cov))
        for allele_id, cov in coverage.items():
            cov2 = [c / max_cov for c in cov]
            coverage[allele_id] = cov2
        self.coverage = coverage
                                
    # Begin drawing graph
    def begin_draw(self, fname_base):
        pdfDraw = self.pdfDraw = open(fname_base + '.pdf', 'w')
        print(r'%PDF-1.7', file=pdfDraw)
        self.objects    = []
        self.stream     = []
        self.draw_items = []
        
    # End drawing graph
    def end_draw(self):
        self.unscaled_height += 50
        self.height = self.unscaled_height * self.scaley
        
        def get_x(x):
            return self.left_margin + x * self.scalex

        def get_y(y):
            return self.height - self.top_margin - y * self.scaley

        # Get scalar
        def get_sx(x):
            return x * self.scalex

        def get_sy(y):
            return y * self.scaley
        
        pdfDraw = self.pdfDraw
        self.add_pdf_object('<</Type /Catalog /Pages 2 0 R>>')
        self.add_pdf_object('<</Type /Pages /Kids [3 0 R] /Count 1>>')
        self.add_pdf_object('<</Type /Page /Parent 2 0 R /Resources 4 0 R '\
                            '/MediaBox [0 0 %d %d] /Contents 6 0 R>>' \
                                % (self.width, self.height))
        self.add_pdf_object('<</Font <</F1 5 0 R>>>>')
        self.add_pdf_object('<</Type /Font /Subtype /Type1 /BaseFont /Helvetica>>')

        # Draw vertical dotted lines at every 100nt and thick lines at every 500nt
        pre_items = []
        for pos in range(0, len(self.backbone), 100):
            main_line = (pos != 0 and pos % 500 == 0)
            dic = {"coord": [pos, 2, pos, self.unscaled_height - 2],
                   "stroke" : "0.5 0.5 0.5",
                   "line_width" : 1 if main_line else 0.2}
            if not main_line:
                dic["line_dash"] = "[3] 0"
            pre_items.append(["line", dic])
        self.draw_items = pre_items + self.draw_items

        fill       = "0 0 0"
        stroke     = "0 0 0"
        line_width = 2.0
        line_dash  = ""
        for type, dic in self.draw_items:
            commands = []
            if type != "state":
                assert "coord" in dic

            if "fill" in dic and dic["fill"] != fill:
                fill = dic["fill"]
                commands.append("%s rg" % fill)
            if "stroke" in dic and dic["stroke"] != stroke:
                stroke = dic["stroke"]
                commands.append("%s RG" % stroke)
            if "line_width" in dic and dic["line_width"] != line_width:
                line_width = dic["line_width"]
                commands.append("%.1f w" % line_width)
            if "line_dash" in dic:
                if dic["line_dash"] != line_dash:
                    line_dash = dic["line_dash"]
                    commands.append("%s d" % line_dash)
            elif line_dash != "":
                line_dash = ""
                commands.append("[] 0 d")
                    
            if type == "rect":
                x, y, sx, sy = dic["coord"]
                re_str = "%d %d %d %d" % (get_x(x), 
                                          get_y(y), 
                                          get_sx(sx), 
                                          get_sy(sy))
                if "fill" in dic:
                    commands.append("%s re f" % re_str)
                if "stroke" in dic:
                    commands.append("%s re S" % re_str)
                    
            elif type == "line":
                x, y, x2, y2 = dic["coord"]
                commands.append("%d %d m %d %d l h S" \
                                    % (get_x(x), get_y(y), get_x(x2), get_y(y2)))
            elif type == "text":
                assert "text" in dic and "font_size" in dic
                x, y = dic["coord"]
                commands.append("BT /F1 %d Tf %d %d Td (%s) Tj ET" \
                                    % (dic["font_size"], 
                                       get_x(x), 
                                       get_y(y), 
                                       dic["text"]))
            else:
                assert type == "state"
                
            self.stream.append(' '.join(commands))

        # Write stream
        self.add_pdf_stream('\n'.join(self.stream))

        # Write xref and trailer
        to_xref = pdfDraw.tell()
        print('xref', file=pdfDraw)
        print("0 %d" % (len(self.objects) + 1), file=pdfDraw)
        print(r'0000000000 65535 f', file=pdfDraw)
        for object in self.objects:
            print("%s 00000 n" % "{:010}".format(object), file=pdfDraw)
        print('trailer <</Size %d /Root 1 0 R>>' \
                % (len(self.objects) + 1), file=pdfDraw)
        print('startxref', file=pdfDraw)
        print(str(to_xref), file=pdfDraw)
        print(r'%%EOF', file=pdfDraw)
        
        self.pdfDraw.close()

    def add_pdf_object(self, obj):
        self.objects.append(self.pdfDraw.tell())
        print("%d 0 obj %s" % (len(self.objects), obj), file=self.pdfDraw)
        print('endobj', file=self.pdfDraw)

    def add_pdf_stream(self, stream):
        self.add_pdf_object("<</Length %d>>\nstream\n%s\nendstream" \
                                % (len(stream), stream))
        
    # Draw graph
    #   Top left as (0, 0) and Bottom right as (width, height)
    def draw(self,
             begin_y,
             title = ""):
        assert len(self.nodes) > 0
        nodes = [[id, node.left, node.right] for id, node in self.nodes.items()]
        def node_cmp(a, b):
            return a[1] - b[1]
        nodes     = sorted(nodes, cmp=node_cmp)
        max_right = len(self.backbone)

        # display space
        end_y  = begin_y + 10000
        dspace = [[[begin_y, end_y]]] * (max_right + 1)
        def get_dspace(left, right, height):
            assert left < len(dspace) and right < len(dspace)
            range1 = dspace[left]
            for range2 in dspace[left + 1:right + 1]:
                new_range = []
                # sub range
                for t1, b1 in range1:
                    for t2, b2 in range2:
                        if b1 < t2:
                            break
                        if b2 < t1:
                            continue
                        t = max(t1, t2)
                        b = min(b1, b2)
                        if b - t >= height:
                            new_range.append([t, b])

                range1 = new_range
            if len(range1) <= 0:
                return -1

            t, b = range1[0]
            assert b - t >= height
            b = t + height
            for i in range(left, right+1):
                range1 = dspace[i]
                range2 = []
                found  = False
                for j in range(len(range1)):
                    t2, b2 = range1[j]
                    if t2 <= t and b <= b2:
                        found = True
                        if t2 < t:
                            range2.append([t2, t])
                        if b < b2:
                            range2.append([b, b2])
                    else:
                        range2.append([t2, b2])
                dspace[i] = range2
                assert found
            return t

        def get_x(x):
            return self.left_margin + x * self.scalex

        def get_y(y):
            return self.height - self.top_margin - y * self.scaley

        # Get scalar
        def get_sx(x):
            return x * self.scalex

        def get_sy(y):
            return y * self.scaley

        # Draw exons
        y = get_dspace(0, max_right, 14)
        for e in range(len(self.exons)):
            left, right = self.exons[e]
            right += 1

            # Draw exon
            self.draw_items.append(["rect",
                                    {"coord" : [left, y + 10, right - left, 10],
                                     "fill" : "1 1 1",
                                     "stroke" : "0 0 0",
                                     "line_width" : 2}])

            primary = False
            for left_, _ in self.primary_exons:
                if left == left_:
                    primary = True
                    break                

            # Draw label
            self.draw_items.append(["text",
                                    {"coord" : [left + 2, y + 7],
                                     "text" : "Exon %d%s" % (e+1, 
                                                            " (primary)" \
                                                                if primary else ""),
                                     "fill" : "0 0 0",
                                     "font_size" : 12}])
            if e > 0:
                prev_right = self.exons[e-1][1] + 1
                self.draw_items.append(["line",
                                        {"coord": [prev_right, 
                                                   y + 5, 
                                                   left, 
                                                   y + 5],
                                         "line_width" : 2}])

        # Draw backbone sequence
        y = get_dspace(0, max_right, 4)
        for pos in range(len(self.backbone)):
            base = self.backbone[pos]
            self.draw_items.append(["text",
                                    {"coord" : [pos, y + 2],
                                     "text" : base,
                                     "fill" : "0.5 0 0.5",
                                     "font_size" : 8}])

        # Draw true or predicted alleles
        node_colors = ["1 1 0", "0 1 0", "1 0.8 0.64", "0.76 0.27 0.5"]
        allele_node_colors = ["0.87 0.87 0", 
                              "0 0.53 0", 
                              "0.87 0.66 0.5", 
                              "0.63 0.14 0.38"]
        def draw_alleles(allele_node_dic, 
                         allele_node_colors, 
                         display = False):
            if len(allele_node_dic) <= 0:
                return
            allele_nodes, \
              seqs, \
              colors \
                = self.get_node_comparison_info(allele_node_dic)

            def draw_coverage(allele_node, 
                              allele_id, 
                              left, 
                              right, 
                              allele_node_color):
                if allele_id not in self.coverage:
                    return
                y = get_dspace(0, max_right, 14)
                for p in range(left, right):
                    cov = math.ceil(self.coverage[allele_id][p] * 12)
                    self.draw_items.append(["rect",
                                            {"coord" : [p, y + 13, 1, cov],
                                             "fill" : allele_node_color}])

            for n_ in range(len(allele_nodes)):
                n    = -1
                prob = ""
                if not display \
                        and not self.simulation \
                        and len(self.allele_node_order) == len(allele_node_dic):
                    allele_id, prob = self.allele_node_order[n_]
                    for n2_ in range(len(allele_nodes)):
                        if allele_id == allele_nodes[n2_][0]:
                            n = n2_
                            break
                    prob = ": %.2f" % prob
                else:
                    n = n_
                assert n >= 0 and n < len(allele_nodes)
                allele_id, left, right = allele_nodes[n]
                right += 1
                allele_node = allele_node_dic[allele_id]
                allele_node_color = allele_node_colors[n % len(allele_node_colors)]

                draw_coverage(allele_node, 
                              allele_id, 
                              left, 
                              right, 
                              allele_node_color)
                
                y = get_dspace(0, max_right, 14)

                # Draw allele name
                if display:
                    allele_type = "display"
                else:
                    if self.simulation:
                        allele_type = "true"
                    else:
                        allele_type = "predicted"
                text_ = "%s (%s, %s)" % (allele_id, 
                                         "partial" \
                                             if allele_id in self.partial_allele_ids \
                                             else "full", 
                                         allele_type)
                self.draw_items.append(["text",
                                    {"coord" : [-55, y + 7],
                                     "text" : text_,
                                     "fill" : "0 0 1",
                                     "font_size" : 18}])
                # Draw node
                self.draw_items.append(["rect",
                                        {"coord" : [left, y + 10, right - left, 10],
                                         "fill" : allele_node_color,
                                         "stroke" : "0 0 0",
                                         "line_width" : 2}])


                color_boxes = []
                c = 0
                while c < len(colors[n]):
                    color = colors[n][c]
                    c2    = c + 1
                    if color != 'N':                        
                        while c2 < len(colors[n]):
                            color2 = colors[n][c2]
                            if color != color2:
                                break
                            c2 += 1
                        color_boxes.append([c, c2, color])
                    c = c2

                # Draw variants
                for color_box in color_boxes:
                    cleft, cright, color = color_box
                    cleft  += left
                    cright += left
                    if color == 'B':
                        color = "0 0 1" # blue 
                    else:
                        color = "0.12 0.56 1"
                    # DK - debugging purposes
                    color = "0 0 1"
                    self.draw_items.append(["rect",
                                            {"coord" : [cleft, 
                                                        y + 9, 
                                                        cright - cleft, 
                                                        8],
                                             "fill" : color}])

            return allele_nodes, seqs, colors

        if self.simulation:
            allele_nodes, seqs, colors = draw_alleles(self.true_allele_nodes,
                                                      allele_node_colors)
        else:
            allele_nodes, seqs, colors = draw_alleles(self.predicted_allele_nodes,
                                                      allele_node_colors)

        draw_alleles(self.display_allele_nodes,
                     ["1 0.96 0.95"],
                     True) # display alleles?

        # Draw location at every 100bp
        y = get_dspace(0, nodes[-1][2], 14)
        for pos in range(0, nodes[-1][2], 100):
            # Draw label
            self.draw_items.append(["text",
                                    {"coord" : [pos + 1, y + 2],
                                     "text" : "%d" % (pos + 1),
                                     "fill" : "0 0 0",
                                     "font_size" : 10}])
                
        # Draw nodes
        node_to_y  = {}
        draw_title = False
        for id, left, right in nodes:
            node   = self.nodes[id]
            nodes2 = [[node, left, right]]
            if id in self.other_nodes:
                for node in self.other_nodes[id]:
                    nodes2.append([node, node.left, node.right])
                    if left > node.left:
                        left = node.left
                    if right < node.right:
                        right = node.right

            # Get y position
            y = get_dspace(left, right, 14 * len(nodes2))
            for node, left, right in nodes2:
                if y < 0:
                    continue
                node_to_y[id] = y

                node_vars    = node.get_vars()
                node_var_ids = node.get_var_ids()
                if len(nodes2) > 1:
                    color = "0.85 0.85 0.85"
                elif len(allele_nodes) > 0:
                    color = "1 1 1"
                    max_common = -sys.maxsize
                    for a in range(len(allele_nodes)):
                        allele_node_id, \
                          allele_left, \
                          allele_right \
                            = allele_nodes[a]
                        if right - left <= 500 \
                                and (left < allele_left or right > allele_right):
                            continue
                        if self.simulation:
                            allele_node = self.true_allele_nodes[allele_node_id]
                        else:
                            allele_node = self.predicted_allele_nodes[allele_node_id]
                        allele_vars = allele_node.get_var_ids(left, right)
                        common_vars = set(node_var_ids) & set(allele_vars)
                        tmp_common = len(common_vars) \
                                        - len(set(node_var_ids) | set(allele_vars))
                        if max_common < tmp_common:
                            max_common = tmp_common
                            color = node_colors[a % len(node_colors)]
                        elif max_common == tmp_common:
                            color = "1 1 1"
                else:
                    color = "1 1 0" # yellow

                # Draw node
                right += 1
                self.draw_items.append(["rect",
                                        {"coord" : [left, y + 10, right - left, 10],
                                         "fill" : color,
                                         "stroke" : "0 0 0",
                                         "line_width" : 2}])
                
                # Draw variants
                for var_id, pos in node_vars:
                    if var_id == "gap":
                        var_type, var_left = "single", pos
                        color = "0 0 0"
                    elif var_id == "unknown" or var_id.startswith("nv"):
                        var_type, var_left = "single", pos
                        color = "1 0 0"
                    else:
                        var_type, var_left, var_data = self.gene_vars[var_id]
                        color = "0 0 1"
                    if var_type == "single":
                        var_right = var_left + 1
                    elif var_type == "insertion":
                        var_right = var_left + len(var_data)
                    else:
                        assert var_type == "deletion"
                        var_right = var_left + int(var_data)
                    self.draw_items.append(["rect",
                                            {"coord" : [var_left, 
                                                        y + 9, 
                                                        var_right - var_left, 
                                                        8],
                                             "fill" : color}])

                # Draw label
                if get_sx(right - left) >= 300:
                    self.draw_items.append(["text",
                                            {"coord" : [left + 2, y + 7],
                                             "text" : node.id,
                                             "fill" : "0 0 1",
                                             "font_size" : 12}])
            
                if not draw_title:
                    draw_title = True
                    self.draw_items.append(["text",
                                            {"coord" : [-68, y + 7],
                                             "text" : title,
                                             "fill" : "0 0 0",
                                             "font_size" : 24}])        
                y += 14

        curr_y = get_dspace(0, nodes[-1][2], 1)
        self.unscaled_height = curr_y if curr_y > 0 else end_y
        return self.unscaled_height

