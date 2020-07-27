#!/usr/bin/env python
# --------------------------------------------------------------------------- #
# Copyright 2017, Daehwan Kim <infphilo@gmail.com>                            #
#                                                                             #
# This file is part of HISAT-genotype. Base wrapper version of hisatgenotype  #
# to run typing only                                                          #
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
import random
from argparse import ArgumentParser, FileType
from hisatgenotype_typing_core import genotyping_locus
import hisatgenotype_args as arguments

# --------------------------------------------------------------------------- #
# This is the Wrapper script that runs the core code found in                 #
# hisatgenotype_modules/hisatgenotype_typing_core                             #
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    parser = ArgumentParser(
        description='hisatgenotype_locus')

    # Add Arguments
    arguments.args_databases(parser,
                             True) # Add option to rename genotype_genome
    arguments.args_bamfile(parser)
    arguments.args_input_output(parser, 
                                indir=False) # Remove --read-dir option
    arguments.args_aligner_inputs(parser,
                                  True) # Add option to keep alignments
    arguments.args_set_aligner(parser,
                               False) # No option to set missmatch
    arguments.args_locus(parser)
    arguments.args_no_partial(parser)
    arguments.args_assembly(parser)
    arguments.args_common(parser,
                     debug = True) # Add option for debug

    args = parser.parse_args()
    if args.locus_list == "":
        locus_list = []
    else:
        locus_list = args.locus_list.split(',')
        if args.base_fname == "genome":
            assert ':' in args.locus_list
            for i in range(len(locus_list)):
                assert ':' in locus_list[i] and '-' in locus_list[i]
                chr, coord    = locus_list[i].split(':')
                left, right   = coord.split('-')
                locus_list[i] = [chr, int(left), int(right)]

    if args.only_locus_list == "":
        only_locus_list = []
    else:
        locus_list = only_locus_list = args.only_locus_list.split(',')

    if args.aligner == "":
        print("Error: --aligners must be non-empty.", file=sys.stderr)
        sys.exit(1)
    
    if ',' not in args.aligner and '.' not in args.aligner:
        if args.graph_index:
            args.aligner += '.graph'
        else:
            args.aligner += '.linear'
    args.aligner = args.aligner.split(',')
    for i in range(len(args.aligner)):
        args.aligner[i] = args.aligner[i].split('.')

    if args.read_fname_U != "":
        args.read_fname = [args.read_fname_U]
    elif args.read_fname_1 != "" or args.read_fname_2 != "":
        if args.read_fname_1 == "" or args.read_fname_2 == "":
            print("Error: please specify both -1 and -2.", 
                  file=sys.stderr)
            sys.exit(1)
        args.read_fname = [args.read_fname_1, args.read_fname_2]
    else:
        args.read_fname = []
    if args.alignment_fname != "" and \
            not os.path.exists(args.alignment_fname):
        print("Error: %s doesn't exist." % args.alignment_fname, 
              file=sys.stderr)
        sys.exit(1)

    if args.verbose and args.verbose_level == 0:
        args.verbose_level = 1
        
    debug = {}
    debug_opts = ["basic", "full", "pair", "test_list", "test_id", "single-end"]
    if args.debug != "":
        for item in args.debug.split(','):
            if ':' in item:
                fields = item.split(':')
                assert len(fields) >= 2
                key, value = fields[0], ':'.join(fields[1:])
                debug[key] = value
            else:
                debug[item] = 1
        for item in debug:
            if item not in debug_opts:
                print("Warning: %s not valid option for debug" % item, 
                      file=sys.stderr)
                exit(1)
            else:
                continue

    if not args.partial:
        print("Warning: --no-partial should be used for debugging purpose only.", 
              file=sys.stderr)

    if args.read_len * 2 > args.fragment_len:
        print("Warning: fragment might be too short (%d)" % (args.fragment_len), 
              file=sys.stderr)

    skip_fragment_regions = []
    if args.skip_fragment_regions != "":
        prev_left, prev_right = -1, -1
        for region in args.skip_fragment_regions.split(','):
            left, right = region.split('-')
            left  = int(left),
            right = int(right)
            assert left < right
            assert prev_right < left
            prev_left  = left
            prev_right = right
            skip_fragment_regions.append([left, right])

    if args.display_alleles == "":
        display_alleles = []
    else:
        display_alleles = args.display_alleles.split(',')

    if args.out_dir == "":
        args.out_dir = os.getcwd()
    else:
        if not os.path.exists(args.out_dir):
            print("Out directory doesn't exist", 
                  file=sys.stderr)
        exit(1)

    random.seed(args.random_seed)
    genotyping_locus(args.base_fname.lower(),
                     locus_list,
                     args.genotype_genome,
                     only_locus_list,
                     args.partial,
                     args.aligner,
                     args.read_fname,
                     args.fastq,
                     args.alignment_fname,
                     args.threads,
                     args.simulate_interval,
                     args.read_len,
                     args.fragment_len,
                     args.best_alleles,
                     args.num_editdist,
                     args.perbase_errorrate,
                     args.perbase_snprate,
                     skip_fragment_regions,
                     args.assembly,
                     args.output_base,
                     args.error_correction,
                     args.keep_alignment,
                     args.discordant,
                     args.type_primary_exons,
                     args.remove_low_abundance_alleles,
                     display_alleles,
                     args.verbose_level,
                     args.assembly_verbose,
                     args.out_dir,
                     debug)

