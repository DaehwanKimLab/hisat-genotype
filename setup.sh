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

# Set working directory
HG_DIR=$(pwd)

add_to_bash(){
    echo "Adding to $1 and sourcing"
    echo "export PATH=$HG_DIR:\$PATH" >> $1
    echo "export PYTHONPATH=$HG_DIR/hisatgenotype_modules:\$PYTHONPATH" >> $1
    source $1
}

# Files to check for
DOWNLOADED="hisat2.cpp"
BUILT="hisat2-align-s"
BASHRC=~/.bashrc_tmp
BASH_PROFILE=~/.bash_profile_tmp

### This section downloads and sets-up hisat2 submodule
echo "Setting up HISAT2"
# Move to hisat2 submodule directory
cd hisat2

if test ! -f "$BUILT"; then
    if test ! -f "$DOWNLOADED"; then
        echo "> Gathering Module"
        git sumbodule init
        git submodule update

    fi

    echo "> Initiating Build"
    make hisat2-align-s hisat2-build-s hisat2-inspect-section
fi

# Return to hisatgenotype directory
cd ../

# Add PATH lines to BASH
if test -f "$BASHRC"; then
    add_to_bash "$BASHRC"
elif test -f "$BASH_PROFILE"; then
    add_to_bash "$BASH_PROFILE"
else
    add_to_bash "$BASHRC"
fi