#!/bin/bash

BLUE='\033[0;34m'
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# terminate script
die() {
	echo -e "${RED}$1${NC}" >&2
	echo
	exit 1
}

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 "$@"
    else
        "$@"
    fi
}

mkdir -p test/tmp || die "Creating the tmp directory failed"

if [ ! -f test/tmp/genome_chr22.fa ]; then
    wget  -N -O test/tmp/genome_chr22.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" || die "Downloading the genome chr22 failed"
    gzip -d test/tmp/genome_chr22.fa.gz || die "Unzipping the genome chr22 failed"
fi

echo -e "${BLUE}Test 1: view hifi${NC}"
ex  ./minimod view -r test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test1.tsv  || die "Test 1: Running the tool failed"
diff -q test/expected/test1.tsv test/tmp/test1.tsv || die "Test 1: diff failed"

echo -e "${BLUE}Test 2: view ont${NC}"
ex  ./minimod view -r test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2.tsv || die "Test 2: Running the tool failed"
diff -q test/expected/test2.tsv test/tmp/test2.tsv || die "Test 2: diff failed"

echo -e "${BLUE}Test 3: meth-freq hifi${NC}"
ex  ./minimod meth-freq -r test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test3.tsv  || die "Test 3: Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test3.tsv > test/expected/test3.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test3.tsv > test/tmp/test3.tsv.sorted
diff -q test/expected/test3.tsv.sorted test/tmp/test3.tsv.sorted || die "Test 3: tsv diff failed"

echo -e "${BLUE}Test 4: meth-freq hifi bedmethyl output${NC}"
ex  ./minimod meth-freq -b -r test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test3.bedmethyl  || die "Test 3: Running the tool failed"
sort -k1,1 -k2,2n -k6,6 test/expected/test3.bedmethyl > test/expected/test3.bedmethyl.sorted
sort -k1,1 -k2,2n -k6,6 test/tmp/test3.bedmethyl > test/tmp/test3.bedmethyl.sorted
diff -q test/expected/test3.bedmethyl.sorted test/tmp/test3.bedmethyl.sorted || die "Test 3: bedmethyl diff failed"

echo -e "${BLUE}Test 5: meth-freq ont${NC}"
ex  ./minimod meth-freq -r test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test4.tsv || die "Test 4: Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test4.tsv > test/expected/test4.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test4.tsv > test/tmp/test4.tsv.sorted
diff -q test/expected/test4.tsv.sorted test/tmp/test4.tsv.sorted || die "Test 4: tsv diff failed"

echo -e "${BLUE}Test 6: meth-freq ont bedmethyl output${NC}"
ex  ./minimod meth-freq -b -r test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test4.bedmethyl || die "Test 4: Running the tool failed"
sort -k1,1 -k2,2n -k6,6 test/expected/test4.bedmethyl > test/expected/test4.bedmethyl.sorted
sort -k1,1 -k2,2n -k6,6 test/tmp/test4.bedmethyl > test/tmp/test4.bedmethyl.sorted
diff -q test/expected/test4.bedmethyl.sorted test/tmp/test4.bedmethyl.sorted || die "Test 4: bedmethyl diff failed"

echo -e "${GREEN}ALL TESTS PASSED !${NC}"