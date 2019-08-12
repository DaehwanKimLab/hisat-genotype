#!/usr/bin/env python
#
# Copyright 2019, Christopher Bennett <christopher@bennett-tech.dev>
#
# This file is part of HISAT-genotype.
#
# HISAT-genotype is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT-genotype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT-genotype.  If not, see <http://www.gnu.org/licenses/>.
#

import os, sys, glob
from argparse import ArgumentParser
import hisatgenotype_typing_common as typing_common
import hisatgenotype_args as hg_args

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Script for simplifying HISAT-genotype results')

    hg_args.args_input_output(parser)
    args = parser.parse_args()

    if args.read_dir:
        indir = args.read_dir
    else:
        indir = '.'

    reports = glob.glob('%s/*.report' % indir)

    report_results = {}
    for report in reports:
        report_results[report] = typing_common.call_nuance_results(report)

    for file_ in report_results:
        print report_results[file_]['tree']