#!/usr/bin/env bash
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

hisatgenotype --base hla --locus-list A --debug "basic,test_size:5,set_seed:101" --out-dir hg_test1_basic

hisatgenotype --base hla --locus-list A,B --debug "pair,test_size:5,set_seed:100" --out-dir hg_test2_paired

hisatgenotype --base hla --locus-list A --assembly --debug "basic,test_size:1,set_seed:101" --out-dir hg_test3_assembly

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/hla/ILMN.tar.gz
tar xvzf ILMN.tar.gz

hisatgenotype --base hla --locus-list A -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz --out-dir hg_test4_realbasic

hisatgenotype --base hla --locus-list A --assembly -v -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz --out-dir hg_test5_realassembly
