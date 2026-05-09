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
    wget  -N -O test/tmp/genome_chr22.fa "https://raw.githubusercontent.com/imsuneth/shared-files/main/genome_chr22.fa" || die "Downloading the genome chr22 failed"
fi

if [ ! -f test/tmp/genome_chr1.fa ]; then
    wget  -N -O test/tmp/genome_chr1.fa "https://github.com/imsuneth/shared-files/raw/refs/heads/main/genome_chr1.fa" || die "Downloading the genome chr1 failed"
fi

if [ ! -f test/tmp/truth.tsv ]; then
    wget  -N -O test/tmp/truth.tsv "https://raw.githubusercontent.com/imsuneth/shared-files/main/truth.tsv" || die "Downloading the truthset failed"
fi


testname="varview example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -c "m,h" /genome/hg38noAlt.fa test/data/example-ont.bam test/data/example-ont-clair.vcf > test/tmp/example-ont.mm.varview.tsv || die "${testname} failed"
diff -q test/tmp/example-ont.mm.varview.tsv test/expected/example-ont.mm.varview.tsv || die "${testname} failed: output does not match expected output"

testname="varview -b example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" /genome/hg38noAlt.fa test/data/example-ont.bam test/data/example-ont-clair.vcf > test/tmp/example-ont.mm.varview.bed || die "${testname} failed"
diff -q test/tmp/example-ont.mm.varview.bed test/expected/example-ont.mm.varview.bed || die "${testname} failed: output does not match expected output"

testname="varview dna_4mC_5mC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" /genome/hg38noAlt.fa test/data/dna_4mC_5mC_mm_chr22.bam test/data/dna_4mC_5mC_mm_chr22.vcf > test/tmp/dna_4mC_5mC_mm_chr22.mm.varview.bed || die "${testname} failed"
diff -q test/tmp/dna_4mC_5mC_mm_chr22.mm.varview.bed test/expected/dna_4mC_5mC_mm_chr22.mm.varview.bed || die "${testname} failed: output does not match expected output"

testname="varview dna_5mCG_5hmCG_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" /genome/hg38noAlt.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/data/dna_5mCG_5hmCG_mm_chr22.vcf > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.varview.bed || die "${testname} failed"
diff -q test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.varview.bed test/expected/dna_5mCG_5hmCG_mm_chr22.mm.varview.bed || die "${testname} failed: output does not match expected output"

testname="varview dna_5mC_5hmC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" /genome/hg38noAlt.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/data/dna_5mC_5hmC_mm_chr22.vcf > test/tmp/dna_5mC_5hmC_mm_chr22.mm.varview.bed || die "${testname} failed"
diff -q test/tmp/dna_5mC_5hmC_mm_chr22.mm.varview.bed test/expected/dna_5mC_5hmC_mm_chr22.mm.varview.bed || die "${testname} failed: output does not match expected output"

# varfreq tests

testname="varfreq example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" /genome/hg38noAlt.fa test/data/example-ont.bam test/data/example-ont-clair.vcf > test/tmp/example-ont.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/example-ont.mm.varfreq.tsv test/expected/example-ont.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq dna_4mC_5mC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" /genome/hg38noAlt.fa test/data/dna_4mC_5mC_mm_chr22.bam test/data/dna_4mC_5mC_mm_chr22.vcf > test/tmp/dna_4mC_5mC_mm_chr22.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/dna_4mC_5mC_mm_chr22.mm.varfreq.tsv test/expected/dna_4mC_5mC_mm_chr22.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq dna_5mCG_5hmCG_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" /genome/hg38noAlt.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/data/dna_5mCG_5hmCG_mm_chr22.vcf > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.varfreq.tsv test/expected/dna_5mCG_5hmCG_mm_chr22.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq dna_5mC_5hmC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" /genome/hg38noAlt.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/data/dna_5mC_5hmC_mm_chr22.vcf > test/tmp/dna_5mC_5hmC_mm_chr22.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/dna_5mC_5hmC_mm_chr22.mm.varfreq.tsv test/expected/dna_5mC_5hmC_mm_chr22.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq -b example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" /genome/hg38noAlt.fa test/data/example-ont.bam test/data/example-ont-clair.vcf > test/tmp/example-ont.mm.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/example-ont.mm.varfreq.bedmethyl test/expected/example-ont.mm.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

