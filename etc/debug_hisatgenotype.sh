#!/bin/bash

# This bash script is used to test some common errors in hisatgenotype and 
#  should be used when there is no output from hisatgenotype command.

####### Set up script
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>HISATgenoptye_debug.log

TITLE="Error Information"
RIGHT_NOW=$(date +"%x %r %Z")
TIME_STAMP="Event Occured on $RIGHT_NOW"

# Functions
function system_info() {
    echo "System release info: "
    lsb_release -a
    uname -r
}

function assertInstalled() {
    for var in "$@"; do 
        if ! which $var &> /dev/null; then
            echo >&2 "$var cannot be found on your system. Please install before running Debug"
            exit 1
        fi
    done
}

assertInstalled curl samtools

######## Arguments
if [ "$1" == "" ]; then
    echo >&2 "Please provide file to test"
    exit 1
else
    if test -f "$1"; then
        echo "Running test with $1"
    else
        echo >&2 "$1 does not exsist"
        exit 1
    fi
fi

###### Test environment
echo >&1 $(system_info)

samtools --version >&1

echo >&1 "Testing builds"
if type hisat2 > /dev/null 2>&1; then
    hisat2 --version >&1
else
    echo >&2 "No HISAT2 system command."
    exit 1
fi

if type hisatgenotype > /dev/null 2>&1; then
    echo >&1 "hisatgenotype found"
else
    echo >&1 "No HISATgenotype system command."
fi

echo >&1 "Downloading and unpacking:"
if curl https://cloud.biohpc.swmed.edu/index.php/s/grch38_snp/download > grch38_snp.tar.gz; then
    echo >&1 "Downloaded GRCh38 HISAT2 file"
    tar -xzvf grch38_snp.tar.gz
else
    echo >&1 "Cannot download HISAT2 index"
    exit 1
fi

if wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/hla/ILMN.tar.gz; then
    echo >&1 "Running preliminary test"
    tar xvzf ILMN.tar.gz

    hisat2 -x grch38_snp/genome_snp --no-spliced-alignment -X 1000 -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz 2>&1 > /dev/null

fi

echo >&1 "Running File Test"
hisat2 -x grch38_snp/genome_snp --no-spliced-alignment -X 1000 -U $1 2>&1 > /dev/null

##### Cleanup
rm -r grch38_snp/ ILMN/
rm grch38_snp.tar.gz ILMN.tar.gz
