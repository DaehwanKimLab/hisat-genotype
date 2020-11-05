#!/usr/bin/env bash
#
# Copyright 2020, Christopher Bennett <christopher@bennett-tech.dev>
#
# This file is part of HISAT-genotype. It is designed to set-up HISAT-genotype after 
# downloading from github and will add all appropriate links to your path
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

# Set working directory and argument types using getopts
HG_DIR=$(pwd)
OPTIND=1
GET_REF=NO
OMMIT_BASH=NO
OUTFILE=""
while getopts "hbrx:" opt; do
    case "$opt" in
        h) 
            echo "USAGE: setup.sh -hrb"
            echo " -h : Show this help screen"
            echo " -b : Ommit adding PATH to .bashrc or .bash_profile file (Defualt - False)"
            echo " -r : Download the base references for HISAT-genotype to default folder (Default - False)"
            echo " -x : Download base references for HISAT-genotype to folder of choice"
            exit 0
            ;;
        b)  OMMIT_BASH=YES
            ;;
        r)  GET_REF=YES
            OUTFILE=indicies
            ;;           
        x)  GET_REF=YES
            OUTFILE=$OPTARG
            ;;
    esac
done

shift $((OPTIND-1))
[ "${1:-}" = "--" ] && shift

add_to_bash(){
    echo "Adding to $1 and sourcing"
    echo "export PATH=$HG_DIR:$HG_DIR/hisat2:\$PATH" >> $1
    echo "export PYTHONPATH=$HG_DIR/hisatgenotype_modules:\$PYTHONPATH" >> $1
    source $1
}

# Files to check for
DOWNLOADED="hisat2.cpp"
BUILT="hisat2-align-l"
BASHRC=~/.bashrc
BASH_PROFILE=~/.bash_profile

# Exit code for `type -P` is 1 if samtools cannot be found in PATH
type -P samtools &> /dev/null
if [[ $? == 1 ]]; then
    echo "Could not find samtools in PATH"
    exit 1;
fi

## This section downloads and sets-up hisat2 submodule
echo "Setting up HISAT2"

if ! command -v hisat2 &> /dev/null; then
    if test ! -f "$BUILT"; then
        echo "> No HISAT2 found on system"
        if test ! -f "$DOWNLOADED"; then
            echo "> Gathering Module"
            git submodule init
            git submodule update
        fi
        echo "> Initiating Build"
        # Move to hisat2 submodule directory
        cd hisat2
        make
        cd ../
    fi
fi

# Add PATH lines to BASH
if [ "$OMMIT_BASH" == "NO" ]; then
    if test -f "$BASHRC"; then
        add_to_bash "$BASHRC"
    elif test -f "$BASH_PROFILE"; then
        add_to_bash "$BASH_PROFILE"
    else
        add_to_bash "$BASHRC"
    fi
fi

# Download all references for HISAT-genotype
if [  "$GET_REF" == "YES" ]; then
    if ! command -v hisat2 &> /dev/null; then
        echo "Cannot Build Indicies without HISAT2"
        exit 1
    fi

    if [ "$OUTFILE" != "indicies" ]; then
        echo $OUTFILE > hg_ix.link
    fi

    mkdir $OUTFILE
    cd $OUTFILE

    # genotype_genome
    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/genotype_genome_20180128.tar.gz
    tar xvzf genotype_genome_20180128.tar.gz
    rm genotype_genome_20180128.tar.gz

    #grch38
    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
    tar xvzf grch38.tar.gz
    rm grch38.tar.gz
    hisat2-inspect grch38/genome > genome.fa
    samtools faidx genome.fa

    #HISATgenotpye Database
    git clone https://github.com/DaehwanKimLab/hisatgenotype_db.git
    cd -
fi

