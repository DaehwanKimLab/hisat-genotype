#!/usr/bin/env python

import sys, os, subprocess, re
import itertools
from multiprocessing import Pool
from argparse import ArgumentParser
import hisatgenotype_typing_common as typing_common
import hisatgenotype_args as arguments


""" WORK IN PROGRESS
def get_positions (alignemnt_fname):
    pull_align = ['samtools', 'view', alignment_fname]
    proc = subprocess.Popen(pull_align,
                            stdout=subprocess.PIPE,
                            stderr=open("/dev/null", "w"))

    snp_ranges = {}
    cigar_re = re.compile('\d+\w')
    for line in proc.stdout:
        line = line.strip()
        cols = line.split()
        _ , flag, chrom , pos, mapq, cigar_str, _ , _ , _ , seq, qual = cols[:11] # All standard fields from SAM
        flag, pos, mapq, read_len = int(flag), int(pos), int(mapq), len(seq)
        
        if flag & 0x4 != 0:
            continue
        if chrom not in snp_ranges:
            snp_ranges.update({ chrom : set() })

        md = ''
        for itr in range(11, len(cols)):
            opt, _ , value = cols[itr].split(':')
            if opt == "MD":
                md = value # MD field is string of missmatches

        assert md

        cigars = cigar_re.findall(cigar_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]

        # Skip all reads that don't have a variant
        if all(cigar[0] in "MSH" for cigar in cigars) and md.isdigit():
            md = int(md)
            for cigar in cigars:
                if cigar[0] == 'M' and cigar[1] == md:
                    continue
        else:
            for i in range(pos-read_len, pos+(2*read_len)):
                snp_ranges[chrom].add(i)

    def make_ranges(plist):
        plist = sorted(plist)
        for key, group in itertools.groupby(enumerate(plist), lambda t: t[1] - t[0]):
            group = list(group)
            yield group[0][1], group[-1][1]

    for chrom, pos_list in snp_ranges.items():
        pos_list = make_ranges(pos_list)
        snp_ranges.update({ chrom : pos_list })

    return snp_ranges

def get_snps(alignment_fname,
             chromosome,
             positions):
    left, right = positions
    snp_pos, vcf_entry = {}, []
    for itr in range(left,right+1):

    pull_align = ['samtools', 'view', alignment_fname, '%s:%s-%s' % (chromosome, left, right)]
    proc = subprocess.Popen(pull_align,
                            stdout=subprocess.PIPE,
                            stderr=open("/dev/null", "w"))

    snp_pos = {}
    cigar_re = re.compile('\d+\w')
    for line in proc.stdout:
        line = line.strip()
        cols = line.split()
        _ , flag, chrom , pos, mapq, cigar_str, _ , _ , _ , seq, qual = cols[:11] # All standard fields from SAM
        flag, pos, mapq, read_len= int(flag), int(pos), int(mapq), len(seq)
        
        assert chrom == chromosome
        if flag & 0x4 != 0:
            continue

        md, zs, nh = '', '', ''
        for itr in range(11, len(cols)):
            opt, _ , value = cols[itr].split(':')
            if opt == "MD":
                md = value # MD field is string of missmatches
            elif opt == "Zs":
                zs = value # Zs field from Hisat2 contains known snps the read aligned to
            elif opt == "NH" or opt == "X0":
                nh = int(value) # NH in hisat2 and bowtie is number of mapped locations X0 in BWA is number of best locations

        assert md and nh

        cigars = cigar_re.findall(cigar_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]

        # Skip all reads that don't have a variant
        if all(cigar[0] in "MSH" for cigar in cigars) and md.isdigit():
            md = int(md)
            for cigar in cigars:
                if cigar[0] == 'M' and cigar[1] == md:
                    continue

        if zs:
            zs = zs.split(',')
            zs = [[zs_entry.split('|')] for zs_entry in zs]

        # Generate lists of nucleotides and scores
        seq_cstr, weight_cstr,  = [], []
        pval_read = 1-(10**(-mapq/10))

        # for Substitution
        for itr in range(len(seq)):
            qual_score = float(ord(qual[itr])) - 33
            pval_base = 1-(10**(-qual_score/10))
            base_weight = (pval_base*pval_read)/nh

            seq_cstr.append(seq[itr])
            weight_cstr.append(base_weight)

        # for indel
        seq_cstr.append('indel')
        weight_cstr.append(pval_read/nh)

        # Parse md and zs field
        md = re.split('(\d+)', md)
        md_list, zs_list, ptr = [], [], 0 
        for field in md:
            if field.isdigit():
                for itr in range(int(field)):
                    md_list.append("M")
                    zs_list.append("-")
                    ptr += 1
            else:
                md_list.append(field)
                if field == "^":
                    zs_list.append('-')

                if zs:
                    for jtr in range(len(zs)):
                        if ptr == zs[jtr][0]+1:
                            zs_list.append("%s_%s" % (zs[jtr][1], zs[jtr][2]))
                else:
                    zs_list.append('-')

                if field != '^':
                    ptr += 1

        assert len(md_list) == len(zs_list)

        # CIGAR string expansion and replacement
        snp_cstr, known_cstr = [], []
        md_ptr, seq_ptr = 0, 0
        for itr in range(len(cigars)):
            cigar_op, lenth = cigars[itr]

            if cigar_op in "SH":
                seq_ptr += length

            if cigar_op == 'N':
                pos += length

            if cigar_op in "MX=":
                for jtr in range(length):
                    if md_list[md_ptr] != 'M': 
                        snp_cstr.append([pos, md_list[md_ptr]])

                    md_ptr += 1
                    seq_ptr += 1
                    

                    cigar_cstr.append(md_list[md_pos])
                    if zs_list[md_pos] == '-':
                        if md_list[md_pos] == "M":
                            known_cstr.append('-')
                        else:
                            known_cstr.append('U')
                    else:
                        known_cstr.append(zs_list[md_pos])
                    
                    md_pos += 1

            if cigar_op == "D":
                assert md_list[md_pos] == '^'
                md_pos += 1 # Add one to point to correct position in string
                for jtr in range(length):
                    cigar_cstr.append(cigar_op)
                    if zs_list[md_pos] == '-':
                        known_cstr.append('U')
                    else:
                        known_cstr.append(zs_list[md_pos])

                    md_pos += 1

            if cigar_op in "IN":
                for jtr in range(length):
                    cigar_cstr.append(cigar_op)

""" 

def samtools_caller(bam_fname,
                    genome_fname,
                    threads,
                    opts):
    # file names
    bam_sort_fname = 'sorted_%s' % bam_fname
    bcf_fname = bam_fname.replace('.bam', '.bcf')
    
    # commands
    sort_cmd = ['samtools', 'sort'] + opts[0] + ['-@', str(threads), '-o', bam_sort_fname, bam_fname]
    genome_index_cmd = ['samtools', 'faidx', genome_fname]
    bcf_call_cmd_in = ['bcftools', 'mpileup'] + opts[1] + ['-Ou', '--threads', str(threads), '-f', genome_fname, bam_sort_fname]
    bcf_call_cmd_out = ['bcftools', 'call'] + opts[2] + ['--threads', str(threads), '-Ob', '-o', bcf_fname]

    # execute commands
    if not os.path.exists(bam_sort_fname):
        subprocess.call(sort_cmd)
    if not os.path.exists('%s.fai' % genome_fname):
        subprocess.call(genome_index_cmd)
    bcf_in = subprocess.Popen(bcf_call_cmd_in, stdout=subprocess.PIPE)
    bcf_out = subprocess.Popen(bcf_call_cmd_out, stdin=bcf_in.stdout)

def record_variants(aligner,
                    is_fastq,
                    reference,
                    use_graph,
                    read_fnames,
                    alignment_fname,
                    first_align,
                    keep_align,
                    sam_bcf_opts,
                    threads):
    if first_align:
        if reference == "":
            print("Error: --ref-genome or -x option not set", file=sys.stderr)
            exit(1)

        fname_prefix = read_fnames[0].split('.')[0]
        bam_fname = fname_prefix + '.bam'

        cmd = [aligner]
        if aligner == "hisat2" or aligner == "bowtie2":
            if aligner == "hisat2":
                cmd += ['--no-spliced-alignment']
            cmd += ['-x', reference, '-p', threads]
            if len(read_fnames) > 1:
                cmd += ['-1', read_fnames[0], '-2', read_fnames[1]]
            else:
                cmd += ['-U', read_fnames[0]]
        elif aligner == "bwa":
            cmd += ["mem", reference]
            cmd += read_fnames
            
        cmd += ['|', 'samtools', 'view', '-Sb', '-', '>', bam_fname] 
        os.system(" ".join(cmd))
        alignment_fname.append(bam_fname)

    if sam_bcf_opts:
        samtools_caller(alignment_fname[0], reference, threads, sam_bcf_opts)          
 
if __name__ == '__main__':
    parser = ArgumentParser(
        description='HISAT-Genotype Call Variants')
    parser.add_argument('-x', '--ref-genome',
                        dest="reference",
                        type=str,
                        default = "",
                        required = True,
                        help="Name of reference to use with aligner of choice (Required; ex genome.fa)")        
    arguments.args_aligner_inputs(parser,
                                  keep=True) # Add option to keep alignments
    arguments.args_bamfile(parser)
    arguments.args_set_aligner(parser,
                               missmatch = False) # Turn off option for setting missmatch
    arguments.args_common(parser)

    subparser = parser.add_subparsers(title="Method",
                        dest="method",
                        help="Choose SNP calling method")    
    btools_parser = subparser.add_parser("bcftools",
                         help = 'Set BCFtools as variant caller')
    btools_parser.add_argument('--sam-opts',
                               type=str,
                               dest='sam_opts',
                               default='',
                               help='String of options for samtools sort (-@ option already included)')
    btools_parser.add_argument('--mplp-opts',
                               type=str,
                               dest='mplp_opts',
                               default='',
                               help='String of additional options for bcftools mpileup (-Ou already used)')
    btools_parser.add_argument('--vcall-opts',
                               type=str,
                               dest='vcall_opts',
                               default='-mv',
                               help='String of additional options for bcftools call (-mv are defaults -Ob already used)')

    args = parser.parse_args()

    args.aligner = args.aligner.lower()
    if args.aligner not in ['hisat2', 'bowtie2', 'bwa']:
        print("Error: --alinger supports hisat2, bowtie2, or bwa", file=sys.stderr)
        exit(1)

    read_fnames, bam_fname = [], []
    if args.read_fname_U:
        read_fnames = [args.read_fname_U]
    elif args.read_fname_1 or args.read_fname_2:
        if not args.read_fname_1 or not args.read_fname_2:
            print("Error: Please specify both -1 and -2 options", file=sys.stderr)
            exit(1)
        read_fnames = [args.read_fname_1,
                       args.read_fname_2]
    elif args.alignment_fname:
        bam_fname = [args.alignment_fname]
    else:
        print("Error: Please provide file options", file=sys.stderr)
        exit(1)

    if not ((read_fnames or bam_fname) and not (read_fnames and bam_fname)):
        print("Error: Please don't use both --bamfile and -U / -1,-2 options", file=sys.stderr)
        exit(1)

    sopts = []
    if args.method == 'bcftools':
        try:
            sopts = [args.sam_opts.split(),
                     args.mplp_opts.split(),
                     args.vcall_opts.split()]
        except:
            pass

    record_variants(args.aligner,
                    args.fastq,
                    args.reference,
                    args.graph_index,
                    read_fnames,
                    bam_fname,
                    False if bam_fname else True,
                    args.keep_alignment,
                    sopts,
                    args.threads)

