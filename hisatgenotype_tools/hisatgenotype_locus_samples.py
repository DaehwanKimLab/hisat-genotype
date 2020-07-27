#!/usr/bin/env python
# --------------------------------------------------------------------------- #
# Copyright 2017, Daehwan Kim <infphilo@gmail.com>                            #
#                                                                             #
# This file is part of HISAT-genotype. This was the preliminary method for    #
# running samples multithreaded on a cloud                                    #
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
import threading
import inspect
import random
import glob
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common
import hisatgenotype_args as arguments


# Platinum genomes - CEPH pedigree (17 family members)
CEPH_pedigree = {
    "NA12889" : {"gender" : "M", 
                 "spouse" : "NA12890", 
                 "children" : ["NA12877"]},
    "NA12890" : {"gender" : "F", 
                 "spouse" : "NA12889", 
                 "children" : ["NA12877"]},
    "NA12877" : {"gender" : "M", 
                 "father" : "NA12889", 
                 "mother" : "NA12890", 
                 "spouse" : "NA12878", 
                 "children" : ["NA12879", 
                               "NA12880", 
                               "NA12881", 
                               "NA12882", 
                               "NA12883", 
                               "NA12884", 
                               "NA12885", 
                               "NA12886", 
                               "NA12887", 
                               "NA12888", 
                               "NA12893"]},
    "NA12891" : {"gender" : "M", 
                 "spouse" : "NA12892", 
                 "children" : ["NA12878"]},
    "NA12892" : {"gender" : "F", 
                 "spouse" : "NA12891", 
                 "children" : ["NA12878"]},
    "NA12878" : {"gender" : "F", 
                 "father" : "NA12892", 
                 "mother" : "NA12891", 
                 "spouse" : "NA12877", 
                 "children" : ["NA12879", 
                               "NA12880", 
                               "NA12881", 
                               "NA12882", 
                               "NA12883",
                               "NA12884", 
                               "NA12885", 
                               "NA12886", 
                               "NA12887", 
                               "NA12888", 
                               "NA12893"]},
    "NA12879" : {"gender" : "F", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12880" : {"gender" : "F", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12881" : {"gender" : "F", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12882" : {"gender" : "M", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12883" : {"gender" : "M", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12884" : {"gender" : "M", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12885" : {"gender" : "F", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12886" : {"gender" : "M", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12887" : {"gender" : "F", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12888" : {"gender" : "M", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    "NA12893" : {"gender" : "M", 
                 "father" : "NA12877", 
                 "mother" : "NA12878"},
    }

# --------------------------------------------------------------------------- #
# Threadding Class to multiprocess HISAT-genotype                             #
# --------------------------------------------------------------------------- #
class myThread(threading.Thread):
    def __init__(self,
                 lock, 
                 paths,
                 reference_type,
                 region_list,
                 num_editdist,
                 max_sample,
                 assembly,
                 out_dir,
                 genotype_results,
                 verbose):
        threading.Thread.__init__(self)
        self.lock             = lock
        self.paths            = paths
        self.reference_type   = reference_type
        self.region_list      = region_list
        self.num_editdist     = num_editdist
        self.max_sample       = max_sample
        self.assembly         = assembly
        self.out_dir          = out_dir
        self.genotype_results = genotype_results
        self.verbose          = verbose

    def run(self):
        global work_idx
        while True:
            self.lock.acquire()
            my_work_idx = work_idx
            work_idx   += 1
            self.lock.release()
            if my_work_idx >= len(self.paths) \
                    or my_work_idx >= self.max_sample:
                return
            worker(self.lock,
                   self.paths[my_work_idx],
                   self.reference_type,
                   self.region_list,
                   self.num_editdist,
                   self.assembly,
                   self.out_dir,
                   self.genotype_results,
                   self.verbose)

            
# --------------------------------------------------------------------------- #
# Set up workers to manage threads and paths                                  #
# --------------------------------------------------------------------------- #
work_idx = 0
def worker(lock,
           path,
           reference_type,
           region_list,
           num_editdist,
           assembly,
           out_dir,
           genotype_results,
           verbose):
    fq_name  = path.split('/')[-1]
    read_dir = '/'.join(path.split('/')[:-1])
    genome   = fq_name.split('.')[0]
    if not fq_name.endswith("extracted.1.fq.gz"):
        return
    read_basename = fq_name[:fq_name.find("extracted.1.fq.gz")]
    read_fname_1  = "%s/%sextracted.1.fq.gz" % (read_dir, read_basename)
    read_fname_2  = "%s/%sextracted.2.fq.gz" % (read_dir, read_basename)

    if not os.path.exists(read_fname_1) or not os.path.exists(read_fname_2):
        return
    lock.acquire()
    print(genome, file=sys.stderr)
    lock.release()

    for family, loci in region_list.items():
        test_hla_cmd = ["hisatgenotype_toolkit locus",
                        "--base", family]
        if len(loci) > 0:
            test_hla_cmd += ["--locus", ','.join(loci)]
        test_hla_cmd += ["--num-editdist", str(num_editdist)]
        test_hla_cmd += ["-1", read_fname_1, "-2", read_fname_2]
        if assembly:
            test_hla_cmd += ["--assembly"]
            test_hla_cmd += ["--assembly-base"]
            if out_dir != "":
                test_hla_cmd += ["%s/%s" % (out_dir, genome)]
            else:
                test_hla_cmd += [genome]        

        if verbose:
            lock.acquire()
            print(' '.join(test_hla_cmd), file=sys.stderr)
            lock.release()

        proc = subprocess.Popen(test_hla_cmd, 
                                universal_newlines=True,
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        test_alleles = set()
        output_list  = []
        for line in proc.stdout:
            line = line.strip()
            if line.find("abundance") == -1:
                continue

            rank, _, allele, _, abundance = line.split()        
            output_list.append([allele, float(abundance[:-2])])

    lock.acquire()
    for allele, abundance in output_list:
        print("%s\t%s\t%.2f" % (genome, allele, abundance))
        genotype_results.append([genome, allele, abundance])
    sys.stdout.flush()
    lock.release()


# --------------------------------------------------------------------------- #
# Main Genotyping script to run database building and genotyping              #
# --------------------------------------------------------------------------- #
def genotyping(read_dir,
               reference_type,
               region_list,
               num_editdist,
               nthreads,
               max_sample,
               assembly,
               out_dir,
               verbose,
               platinum_check):
    for database_name in region_list:
        # Extract variants, backbone sequence, and other sequeces
        typing_common.extract_database_if_not_exists(database_name,
                                                     [])            # locus_list
        # Build HISAT2's graph index
        typing_common.build_index_if_not_exists(database_name,
                                                "hisat2",
                                                "graph",
                                                1,            # threads
                                                verbose)
    
    if not os.path.exists(read_dir):
        print("Error: %s does not exist." % read_dir, file=sys.stderr)
        sys.exit(1)

    if out_dir != "" and not os.path.exists(out_dir):
        os.mkdir(out_dir)        

    # fastq files
    fq_fnames = glob.glob("%s/*.extracted.1.fq.gz" % read_dir)

    genotype_results = []

    lock = threading.Lock()
    threads = []
    for t in range(nthreads):
        thread = myThread(lock,
                          fq_fnames,
                          reference_type,
                          region_list,
                          num_editdist,
                          max_sample,
                          assembly,
                          out_dir,
                          genotype_results,
                          verbose)
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()


    if platinum_check:
        genotype_dic = {}
        for genome, allele, abundance in genotype_results:
            region, _ = allele.split('*')
            if region not in genotype_dic:
                genotype_dic[region] = {}
            if genome not in genotype_dic[region]:
                genotype_dic[region][genome] = []
            if len(genotype_dic[region][genome]) >= 2:
                continue
            # DK - debugging purposes
            # if abundance < 0.15 * 100:
            #    continue
            genotype_dic[region][genome].append([allele, abundance])

        for region, region_genotype in genotype_dic.items():
            print(region, file=sys.stderr)
            included = 0
            total    = 0
            for genome, genome_alleles in region_genotype.items():
                genome_alleles = set([allele for allele, _ in genome_alleles])
                if "father" in CEPH_pedigree[genome]:
                    assert "mother" in CEPH_pedigree[genome]
                    parents = [CEPH_pedigree[genome]["father"], 
                               CEPH_pedigree[genome]["mother"]]
                else:
                    parents = []
                parent_allele_sets = []
                assert len(parents) in [0, 2]
                if len(parents) == 2 \
                        and parents[0] in region_genotype \
                        and parents[1] in region_genotype:
                    for parent_allele, _ in region_genotype[parents[0]]:
                        for parent_allele2, _ in region_genotype[parents[1]]:
                            parent_allele_sets.append(set([parent_allele, 
                                                           parent_allele2]))
                print(("\t", genome, genome_alleles, parent_allele_sets), 
                      file=sys.stderr)
                if len(parent_allele_sets) > 0:
                    total += 1
                    if genome_alleles in parent_allele_sets:
                        included += 1
            print("\t%d / %d" % (included, total), file=sys.stderr)

# --------------------------------------------------------------------------- #
# Main function fro script                                                    #
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    parser = ArgumentParser(
        description='genotyping on many samples')

    # Add Arguments
    arguments.args_databases(parser)
    arguments.args_input_output(parser)
    arguments.args_reference_type(parser)
    arguments.args_assembly(parser)
    arguments.args_locus_samples(parser)
    arguments.args_common(parser)

    args = parser.parse_args()

    if args.read_dir == "":
        print("Error: please specify --read-dir.", file=sys.stderr)
        sys.exit(1)

    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print("Error: --reference-type (%s) must be one of gene, "\
                "chromosome, and genome." % (args.reference_type), 
              file=sys.stderr)
        sys.exit(1)

    region_list = {}
    if args.region_list != "":
        for region in args.region_list.split(','):
            region = region.split('.')
            if len(region) < 1 or len(region) > 2:
                print("Error: --region-list is incorrectly formatted.", 
                      file=sys.stderr)
                sys.exit(1)
                
            family = region[0].lower()
            if len(region) == 2:
                locus_name = region[1].upper()
            if family not in region_list:
                region_list[family] = set()
            if len(region) == 2:
                region_list[family].add(locus_name)
    else:
        args.base_fname = args.base_fname.split(',')
        args.locus_list = args.locus_list.split(';')
        if len(args.base_fname) != len(args.locus_list):
            print("Error: --base and --locus-list not correct", file=sys.stderr)
        for itr in range(len(args.base_fname.split)):
            family     = args.base_fname[itr]
            loci_names = args.locus_list[itr].split(',')
            region_list.update({ family : loci_names })

    genotyping(args.read_dir,
               args.reference_type,
               region_list,
               args.num_editdist,
               args.threads,
               args.max_sample,
               args.assembly,
               args.out_dir,
               args.verbose,
               args.platinum_check)

