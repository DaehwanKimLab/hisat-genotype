#!/usr/bin/env python
# --------------------------------------------------------------------------- #
# Copyright 2016, Daehwan Kim <infphilo@gmail.com>                            #
#                                                                             #
# This file is part of HISAT 2. Builds the genotype genome for use with       #
# HISAT-genotype                                                              #
#                                                                             #
# HISAT 2 is free software: you can redistribute it and/or modify             #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# HISAT 2 is distributed in the hope that it will be useful,                  #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.            #
# --------------------------------------------------------------------------- #

import os
import sys
import subprocess
import re
import shutil
import inspect
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common
import hisatgenotype_args as arguments

# --------------------------------------------------------------------------- #
#  General Functions                                                          #
# --------------------------------------------------------------------------- #
def read_clnsig(fname):
    clnsig_dic = {}
    for line in open(fname):
        var_id, gene, clnsig = line.strip().split('\t')
        clnsig_dic[var_id] = [gene, clnsig]
    return clnsig_dic


# --------------------------------------------------------------------------- #
#  Main function for building genotype_genome
# --------------------------------------------------------------------------- #
def build_genotype_genome(base_fname,                          
                          inter_gap,
                          intra_gap,
                          threads,
                          database_list,
                          use_clinvar,
                          use_commonvar,
                          aligner,
                          graph_index,
                          verbose):    
    # Download HISAT2 index
    typing_common.download_genome_and_index()

    # Load genomic sequences
    chr_dic, chr_names, chr_full_names = typing_common.read_genome("genome.fa")

    genotype_vars       = {}
    genotype_haplotypes = {} 
    genotype_clnsig     = {}
    if use_clinvar:
        # Extract variants from the ClinVar database
        CLINVAR_fnames = ["clinvar.vcf.gz",
                          "clinvar.snp",
                          "clinvar.haplotype",
                          "clinvar.clnsig"]

        if not typing_common.check_files(CLINVAR_fnames):
            if not os.path.exists("clinvar.vcf.gz"):
                os.system("wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/"\
                            "vcf_GRCh38/archive/2017/clinvar_20170404.vcf.gz")
            assert os.path.exists("clinvar.vcf.gz")

            extract_cmd = ["hisat2_extract_snps_haplotypes_VCF.py"]
            extract_cmd += ["--inter-gap", str(inter_gap),
                            "--intra-gap", str(intra_gap),
                            "--genotype-vcf", "clinvar.vcf.gz",
                            "genome.fa", "/dev/null", "clinvar"]
            if verbose:
                print("\tRunning:", ' '.join(extract_cmd), 
                      file=sys.stderr)
            proc = subprocess.Popen(extract_cmd, 
                                    stdout=open("/dev/null", 'w'), 
                                    stderr=open("/dev/null", 'w'))
            proc.communicate()
            if not typing_common.check_files(CLINVAR_fnames):
                print("Error: extract variants from clinvar failed!", 
                      file=sys.stderr)
                sys.exit(1)

        # Read variants to be genotyped
        genotype_vars = typing_common.read_variants("clinvar.snp")

        # Read haplotypes
        genotype_haplotypes = typing_common.read_haplotypes("clinvar.haplotype")

        # Read information about clinical significance
        genotype_clnsig = read_clnsig("clinvar.clnsig")

    if use_commonvar:
        # Extract variants from dbSNP database
        # TODO: CB Write script to make local uptodate SNP database from dbSNP 
        # ftp://ftp.ncbi.nlm.nih.gov/snp/database/README.create_local_dbSNP.txt 
        commonvar_fbase = "snp144Common"
        commonvar_fnames = ["%s.snp" % commonvar_fbase,
                            "%s.haplotype" % commonvar_fbase]
        if not typing_common.check_files(commonvar_fnames):
            if not os.path.exists("%s.txt.gz" % commonvar_fbase):
                os.system("wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/"\
                               "database/%s.txt.gz" % commonvar_fbase)
            assert os.path.exists("%s.txt.gz" % commonvar_fbase)
            os.system("gzip -cd %s.txt.gz "\
                         "| awk 'BEGIN{OFS=\"\t\"} "\
                             "{if($2 ~ /^chr/) {$2 = substr($2, 4)}; "\
                              "if($2 == \"M\") {$2 = \"MT\"} print}' > %s.txt" \
                                  % (commonvar_fbase, commonvar_fbase))
            extract_cmd = ["hisat2_extract_snps_haplotypes_UCSC.py",
                           "--inter-gap", str(inter_gap),
                           "--intra-gap", str(intra_gap),
                           "genome.fa", 
                           "%s.txt" % commonvar_fbase, 
                           commonvar_fbase]
            if verbose:
                print("\tRunning:", ' '.join(extract_cmd), file=sys.stderr)
            proc = subprocess.Popen(extract_cmd, 
                                    stdout=open("/dev/null", 'w'), 
                                    stderr=open("/dev/null", 'w'))
            proc.communicate()
            if not typing_common.check_files(commonvar_fnames):
                print("Error: extract variants from clinvar failed!", 
                      file=sys.stderr)
                sys.exit(1)

        # Read variants to be genotyped
        genotype_vars = typing_common.read_variants(commonvar_fnames[0])

        # Read haplotypes
        genotype_haplotypes = typing_common.read_haplotypes(commonvar_fnames[1])

    # Genes to be genotyped
    genotype_genes = {}

    # Read genes or genomics regions
    for database_name in database_list:
        # Extract HLA variants, backbone sequence, and other sequeces
        typing_common.extract_database_if_not_exists(database_name,
                                                     [],            # locus_list
                                                     inter_gap,
                                                     intra_gap,
                                                     True,          # partial?
                                                     verbose)
        locus_fname = "%s.locus" % database_name
        assert os.path.exists(locus_fname)
        for line in open(locus_fname):
            locus_name, \
              chr, \
              left, \
              right, \
              length, \
              exon_str, \
              strand \
                   = line.strip().split()
            left   = int(left)
            right  = int(right)
            length = int(length)
            if chr not in chr_names:
                continue
            if chr not in genotype_genes:
                genotype_genes[chr] = []
            genotype_genes[chr].append([left, 
                                        right, 
                                        length, 
                                        locus_name, 
                                        database_name, 
                                        exon_str, 
                                        strand])

    # Write genotype genome
    var_num       = 0
    haplotype_num = 0
    genome_out_file    = open("%s.fa" % base_fname, 'w')
    locus_out_file     = open("%s.locus" % base_fname, 'w')
    var_out_file       = open("%s.snp" % base_fname, 'w')
    index_var_out_file = open("%s.index.snp" % base_fname, 'w')
    haplotype_out_file = open("%s.haplotype" % base_fname, 'w')
    link_out_file      = open("%s.link" % base_fname, 'w')
    coord_out_file     = open("%s.coord" % base_fname, 'w')
    clnsig_out_file    = open("%s.clnsig" % base_fname, 'w')
    for c in range(len(chr_names)):
        chr           = chr_names[c]
        chr_full_name = chr_full_names[c]
        assert chr in chr_dic
        chr_seq = chr_dic[chr]
        chr_len = len(chr_seq)
        if chr in genotype_genes:
            chr_genes = genotype_genes[chr]
            chr_genes = sorted(chr_genes, key = lambda x: (x[1], x[2], x[3]))
        else:
            chr_genes = []

        chr_genotype_vars = []
        chr_genotype_vari = 0
        if graph_index:
            if chr in genotype_vars:
                chr_genotype_vars = genotype_vars[chr]
            chr_genotype_haplotypes = []
            chr_genotype_hti        = 0
            if chr in genotype_haplotypes:
                chr_genotype_haplotypes = genotype_haplotypes[chr]

        def add_vars(left, 
                     right, 
                     chr_genotype_vari, 
                     chr_genotype_hti, 
                     haplotype_num):
            # Output variants with clinical significance
            while chr_genotype_vari < len(chr_genotype_vars):
                var_left, \
                  var_type, \
                  var_data, \
                  var_id \
                    = chr_genotype_vars[chr_genotype_vari]
                var_right = var_left
                if var_type == "deletion":
                    var_right += var_data
                if var_right > right:
                    break
                if var_right >= left:
                    chr_genotype_vari += 1
                    continue

                out_str = "%s\t%s\t%s\t%d\t%s" % (var_id, 
                                                  var_type, 
                                                  chr, 
                                                  var_left + off, 
                                                  var_data)
                print(out_str, file=var_out_file)
                print(out_str, file=index_var_out_file)

                if var_id in genotype_clnsig:
                    var_gene, clnsig = genotype_clnsig[var_id]
                    print("%s\t%s\t%s" \
                             % (var_id, var_gene, clnsig), file=clnsig_out_file)
                
                chr_genotype_vari += 1

            # Output haplotypes
            while chr_genotype_hti < len(chr_genotype_haplotypes):
                ht_left, ht_right, ht_vars = chr_genotype_haplotypes[chr_genotype_hti]
                if ht_right > right:
                    break
                if ht_right >= left:
                    chr_genotype_hti += 1
                    continue

                print("ht%d\t%s\t%d\t%d\t%s" \
                        % (haplotype_num, 
                           chr, 
                           ht_left + off, 
                           ht_right + off, 
                           ','.join(ht_vars)), 
                      file=haplotype_out_file)
                chr_genotype_hti += 1
                haplotype_num    += 1

            return chr_genotype_vari, chr_genotype_hti, haplotype_num

        out_chr_seq = ""
        off         = 0
        prev_right  = 0
        for gene in chr_genes:
            left, right, length, name, family, exon_str, strand = gene

            if not graph_index:
                # Output gene (genotype_genome.gene)
                print("%s\t%s\t%s\t%d\t%d\t%s\t%s" \
                        % (family.upper(), 
                           name, 
                           chr, 
                           left, 
                           right, 
                           exon_str, 
                           strand), 
                      file=locus_out_file)
                continue            

            chr_genotype_vari, \
              chr_genotype_hti, \
              haplotype_num \
                = add_vars(left, 
                           right, 
                           chr_genotype_vari, 
                           chr_genotype_hti, 
                           haplotype_num)

            # Read gene family sequences and information
            allele_seqs       = typing_common.read_allele_seq("%s_backbone.fa" % family)
            allele_vars       = typing_common.read_variants("%s.snp" % family)
            allele_index_vars = typing_common.read_variants("%s.index.snp" % family)
            allele_haplotypes = typing_common.read_haplotypes("%s.haplotype" % family)
            links             = typing_common.read_links("%s.link" % family, True)

            if name not in allele_seqs:
                continue
            if name not in allele_vars or name not in allele_index_vars:
                vars       = [] 
                index_vars = []
            else:
                vars       = allele_vars[name]
                index_vars = allele_index_vars[name]
                
            allele_seq    = allele_seqs[name]
            index_var_ids = set()
            for _, _, _, var_id in index_vars:
                index_var_ids.add(var_id)

            if name not in allele_haplotypes:
                haplotypes = []
            else:
                haplotypes = allele_haplotypes[name]
            assert length == len(allele_seq)
            assert left < chr_len and right < chr_len
            # Skipping overlapping genes
            if left < prev_right:
                print("Warning: skipping %s ..." % (name), 
                      file=sys.stderr)
                continue

            varID2htID = {}
            assert left < right
            prev_length = right - left + 1
            assert prev_length <= length

            if prev_right < left:
                out_chr_seq += chr_seq[prev_right:left]

            # Output gene (genotype_genome.locus)
            print("%s\t%s\t%s\t%d\t%d\t%s\t%s" \
                    % (family.upper(), 
                       name, 
                       chr, 
                       len(out_chr_seq), 
                       len(out_chr_seq) + length - 1, 
                       exon_str, 
                       strand), 
                  file=locus_out_file)

            # Output coord (genotype_genome.coord)
            print("%s\t%d\t%d\t%d" \
                    % (chr, 
                       len(out_chr_seq), 
                       left, 
                       right - left + 1), 
                  file=coord_out_file)
            out_chr_seq += allele_seq

            # Output variants (genotype_genome.snp and genotype_genome.index.snp)
            for var in vars:
                var_left, var_type, var_data, var_id = var
                new_var_id         = "hv%d" % var_num
                varID2htID[var_id] = new_var_id
                new_var_left       = var_left + left + off
                assert var_type in ["single", "deletion", "insertion"]
                assert new_var_left < len(out_chr_seq)
                if var_type == "single":                    
                    assert out_chr_seq[new_var_left] != var_data
                elif var_type == "deletion":
                    assert new_var_left + var_data <= len(out_chr_seq)
                else:
                    assert var_type == "insertion"

                out_str = "%s\t%s\t%s\t%d\t%s" \
                            % (new_var_id, var_type, chr, new_var_left, var_data)
                print(out_str, 
                      file=var_out_file)
                if var_id in index_var_ids:
                    print(out_str, 
                      file=index_var_out_file)
                var_num += 1
                
            # Output haplotypes (genotype_genome.haplotype)
            for haplotype in haplotypes:
                ht_left, ht_right, ht_vars = haplotype
                new_ht_left  = ht_left + left + off
                assert new_ht_left < len(out_chr_seq)
                new_ht_right = ht_right + left + off
                assert new_ht_left <= new_ht_right
                assert new_ht_right <= len(out_chr_seq)
                new_ht_vars = []
                for var_id in ht_vars:
                    assert var_id in varID2htID
                    new_ht_vars.append(varID2htID[var_id])
                print("ht%d\t%s\t%d\t%d\t%s" \
                        % (haplotype_num, 
                           chr, 
                           new_ht_left, 
                           new_ht_right, 
                           ','.join(new_ht_vars)), 
                      file=haplotype_out_file)
                haplotype_num += 1

            # Output link information between alleles and variants (genotype_genome.link)
            for link in links:
                var_id, allele_names = link
                if var_id not in varID2htID:
                    continue
                new_var_id = varID2htID[var_id]
                print("%s\t%s" % (new_var_id, " ".join(allele_names)), 
                      file=link_out_file)
                
            off       += (length - prev_length)
            prev_right = right + 1

        if not graph_index:
            continue

        # Write the rest of the Vars
        chr_genotype_vari, \
          chr_genotype_hti, \
          haplotype_num \
            = add_vars(sys.maxsize, 
                       sys.maxsize, 
                       chr_genotype_vari, 
                       chr_genotype_hti, 
                       haplotype_num)            
            
        print("%s\t%d\t%d\t%d" \
                % (chr, 
                   len(out_chr_seq), 
                   prev_right, 
                   len(chr_seq) - prev_right), 
              file=coord_out_file)
        out_chr_seq += chr_seq[prev_right:]

        assert len(out_chr_seq) == len(chr_seq) + off

        # Output chromosome sequence
        print(">%s" % (chr_full_name), 
              file=genome_out_file)
        line_width = 60
        for s in range(0, len(out_chr_seq), line_width):
            print(out_chr_seq[s:s+line_width], 
                  file=genome_out_file)

    genome_out_file.close()
    locus_out_file.close()
    var_out_file.close()
    index_var_out_file.close()
    haplotype_out_file.close()
    link_out_file.close()
    coord_out_file.close()
    clnsig_out_file.close()

    allele_out_file = open("%s.allele" % base_fname, 'w')
    if graph_index:
        for database in database_list:
            for line in open("%s.allele" % database):
                allele_name = line.strip()
                print("%s\t%s" % (database.upper(), allele_name), 
                      file=allele_out_file)
    allele_out_file.close()

    partial_out_file = open("%s.partial" % base_fname, 'w')
    if graph_index:
        for database in database_list:
            for line in open("%s.partial" % database):
                allele_name = line.strip()
                print("%s\t%s" % (database.upper(), allele_name), 
                      file=partial_out_file)
    partial_out_file.close()

    if not graph_index:
        shutil.copyfile("genome.fa", "%s.fa" % base_fname)

    # Index genotype_genome.fa
    index_cmd = ["samtools", "faidx", "%s.fa" % base_fname]
    subprocess.call(index_cmd)

    # Build indexes based on the above information
    if graph_index:
        assert aligner == "hisat2"
        build_cmd = ["hisat2-build",
                     "-p", str(threads),
                     "--snp", "%s.index.snp" % base_fname,
                     "--haplotype", "%s.haplotype" % base_fname,
                     "%s.fa" % base_fname,
                     "%s" % base_fname]
    else:        
        assert aligner in ["hisat2", "bowtie2"]
        build_cmd = ["%s-build" % aligner,
                     "-p" if aligner == "hisat2" else "--threads", str(threads),
                     "%s.fa" % base_fname,
                     "%s" % base_fname]
    if verbose:
        print("\tRunning:", ' '.join(build_cmd), 
              file=sys.stderr)
        
    subprocess.call(build_cmd, 
                    stdout=open("/dev/null", 'w'), 
                    stderr=open("/dev/null", 'w'))

    if aligner == "hisat2":
        index_fnames = ["%s.%d.ht2" % (base_fname, i+1) for i in range(8)]
    else:
        index_fnames = ["%s.%d.bt2" % (base_fname, i+1) for i in range(4)]
        index_fnames += ["%s.rev.%d.bt2" % (base_fname, i+1) for i in range(2)]
    if not typing_common.check_files(index_fnames):
        print("Error: indexing failed! "\
               "Perhaps, you may have forgotten to build %s executables?" \
                    % aligner, 
              file=sys.stderr)
        sys.exit(1)


# --------------------------------------------------------------------------- #
# Main function to build index and run script                                 #
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Build genotype genome")

    # Add Arguments
    arguments.args_databases(parser)
    arguments.args_var_gaps(parser)
    arguments.args_set_aligner(parser, 
                               False) # no missmatch option
    arguments.args_build_genome(parser)
    arguments.args_common(parser)

    args = parser.parse_args()
    if args.inter_gap > args.intra_gap:
        print("Error: --inter-gap (%d) must be smaller than --intra-gap (%d)" \
                % (args.inter_gap, args.intra_gap), 
              file=sys.stderr)
        sys.exit(1)
    
    if not args.base_fname:
        args.base_fname = 'genotype_genome'    

    if not args.locus_list:
        database_list = []
        typing_common.clone_hisatgenotype_database()
        for database in os.listdir("hisatgenotype_db"):
            if database in ['.git', 'README.md']:
                continue
            database_list.append(database.lower())
    else:
        database_list = args.locus_list.split(',')

    if args.use_clinvar and args.use_commonvar:
        print("Error: both --clinvar and --commonvar cannot be used together.", 
              file=sys.stderr)
        sys.exit(1)

    if args.aligner not in ["hisat2", "bowtie2"]:
        print("Error: --aligner should be either hisat2 or bowtie2.", 
              file=sys.stderr)
        sys.exit(1)        

    build_genotype_genome(args.base_fname,
                          args.inter_gap,
                          args.intra_gap,
                          args.threads,
                          database_list,
                          args.use_clinvar,
                          args.use_commonvar,
                          args.aligner,
                          args.graph_index,
                          args.verbose)
