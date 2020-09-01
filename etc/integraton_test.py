#!/usr/bin/env python
# --------------------------------------------------------------------------- #
# Copyright 2020, Christopher Bennett <christopher@bennett-tech.dev>          #
#                                                                             # 
# This file is part of HISAT-genotype. It's purpose is to test any new        #
# release of HISATgenotype for proper function.                               #
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

# --------------------------------------------------------------------------- #
# Basic internal tests for HISATgenotypte proper function and data            #
# processing                                                                  #
# --------------------------------------------------------------------------- #
def test_basic():
    script = ['hisatgenotype', 
                '--base', 'hla', 
                '--locus-list', 'A', 
                '--debug', '"basic,test_size:5,set_seed:101"']
    try:
        proc = subprocess.Popen(script,
                                universal_newlines = True,
                                stdout = subprocess.PIPE,
                                stderr = open("/dev/null", 'w'))

        tests = {
                    "Using HISAT2"   : False,   # Test if properly using HISAT2
                    "Aligning"       : False,   # Test if getting alignments back
                    "Counting Reads" : False,   # Test if alleles are getting reads
                    "EM Results"     : False,   # Test if EM results are returned
                    "Passing"        : False,   # Test if passing results are returned
                    "Test number"    : False    # Test if 5 tests are run
                }

        for line in proc.stdout:
            line = line.strip()
            if "hisat2 graph" in line:
                tests["Using HISAT2"] = True
            elif line.endswith("aligned"):
                tests["Aligning"] = True
            elif "count:" in line:
                tests["Counting Reads"] = True
            elif "abundance:" in line:
                tests["EM Results"] = True
            elif "Passed so far: 5/5 (100.00%)" in line:
                tests["Passing"] = True
            elif "Test 5" in line:
                tests["Test number"] = True
            else:
                continue

    except Exception as err:
        print(err, file = sys.stderr)
        print("Error on line {}".format(sys.exc_info()[-1].tb_lineno))
        exit(1)

def test_assembly():
    script = ['hisatgenotype', 
                '--base', 'hla', 
                '--locus-list', 'A',
                '--assembly',
                '--debug', '"basic,test_size:1,set_seed:101"']
    try:
        proc = subprocess.Popen(script,
                                universal_newlines = True,
                                stdout = subprocess.PIPE,
                                stderr = open("/dev/null", 'w'))

        tests = {
                    "Proper Run"    : False,    # Test if basic run starts
                    "Assembly"      : True,    # Test if assembly error is active
                    "Viturbi"       : False,    # Test fidelity of Viturbi Phasing
                    "Output Filled" : False     # Test if outpit PDF is filled
                }

        for line in proc.stdout:
            line = line.strip()
            if "abundance:" in line:
                tests["Proper Run"] = True
            elif line == "Error in building and calling viterbi":
                tests["Assembly"] = False
            elif line == "A: A*11:29 : A*11:29 (Group score: 1.00000)":
                tests["Viturbi"] = True
            else:
                continue

        with open("hisatgenotype_out/assembly_graph-1.hla.A.pdf", "r") as ifi:
            for line in ifi:
                line = line.strip()
                if line == "%\%\EOF":
                    tests["Output Filled"] = True

    except Exception as err:
        print(err, file = sys.stderr)
        print("Error on line {}".format(sys.exc_info()[-1].tb_lineno))
        exit(1)

        