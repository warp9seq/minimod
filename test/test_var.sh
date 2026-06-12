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
ex ./minimod varview -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont.bam test/data/example-ont-clair.vcf > test/tmp/example-ont.mm.varview.tsv || die "${testname} failed"
diff -q test/tmp/example-ont.mm.varview.tsv test/expected/example-ont.mm.varview.tsv || die "${testname} failed: output does not match expected output"

testname="varview -b example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont.bam test/data/example-ont-clair.vcf > test/tmp/example-ont.mm.varview.bed || die "${testname} failed"
diff -q test/tmp/example-ont.mm.varview.bed test/expected/example-ont.mm.varview.bed || die "${testname} failed: output does not match expected output"

testname="varview dna_4mC_5mC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/data/dna_4mC_5mC_mm_chr22.vcf > test/tmp/dna_4mC_5mC_mm_chr22.mm.varview.bed || die "${testname} failed"
diff -q test/tmp/dna_4mC_5mC_mm_chr22.mm.varview.bed test/expected/dna_4mC_5mC_mm_chr22.mm.varview.bed || die "${testname} failed: output does not match expected output"

testname="varview dna_5mCG_5hmCG_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/data/dna_5mCG_5hmCG_mm_chr22.vcf > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.varview.bed || die "${testname} failed"
diff -q test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.varview.bed test/expected/dna_5mCG_5hmCG_mm_chr22.mm.varview.bed || die "${testname} failed: output does not match expected output"

testname="varview dna_5mC_5hmC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/data/dna_5mC_5hmC_mm_chr22.vcf > test/tmp/dna_5mC_5hmC_mm_chr22.mm.varview.bed || die "${testname} failed"
diff -q test/tmp/dna_5mC_5hmC_mm_chr22.mm.varview.bed test/expected/dna_5mC_5hmC_mm_chr22.mm.varview.bed || die "${testname} failed: output does not match expected output"

# varfreq tests

testname="varfreq example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont.bam test/data/example-ont-clair.vcf > test/tmp/example-ont.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/example-ont.mm.varfreq.tsv test/expected/example-ont.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq -b example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont.bam test/data/example-ont-clair.vcf > test/tmp/example-ont.mm.varfreq.bed || die "${testname} failed"
diff -q test/tmp/example-ont.mm.varfreq.bed test/expected/example-ont.mm.varfreq.bed || die "${testname} failed: output does not match expected output"

testname="varfreq -b dna_4mC_5mC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/data/dna_4mC_5mC_mm_chr22.vcf > test/tmp/dna_4mC_5mC_mm_chr22.mm.varfreq.bed || die "${testname} failed"
diff -q test/tmp/dna_4mC_5mC_mm_chr22.mm.varfreq.bed test/expected/dna_4mC_5mC_mm_chr22.mm.varfreq.bed || die "${testname} failed: output does not match expected output"

testname="varfreq -b dna_5mCG_5hmCG_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/data/dna_5mCG_5hmCG_mm_chr22.vcf > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.varfreq.bed || die "${testname} failed"
diff -q test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.varfreq.bed test/expected/dna_5mCG_5hmCG_mm_chr22.mm.varfreq.bed || die "${testname} failed: output does not match expected output"

testname="varfreq -b dna_5mC_5hmC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/data/dna_5mC_5hmC_mm_chr22.vcf > test/tmp/dna_5mC_5hmC_mm_chr22.mm.varfreq.bed || die "${testname} failed"
diff -q test/tmp/dna_5mC_5hmC_mm_chr22.mm.varfreq.bed test/expected/dna_5mC_5hmC_mm_chr22.mm.varfreq.bed || die "${testname} failed: output does not match expected output"

# hg002 PGXXSX240041 promethion chr22 - covers all variant types: SNP, INS, DEL, multiallelic, RefCall (chr22:17280000-17380000)
testname="varfreq hg002_prom_chr22_snp_ins_del_multiallele_refcall"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" test/tmp/genome_chr22.fa test/data/hg002_prom_chr22_snp_ins_del_multiallele_refcall.bam test/data/hg002_prom_chr22_snp_ins_del_multiallele_refcall.vcf > test/tmp/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.tsv test/expected/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq -b hg002_prom_chr22_snp_ins_del_multiallele_refcall"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/hg002_prom_chr22_snp_ins_del_multiallele_refcall.bam test/data/hg002_prom_chr22_snp_ins_del_multiallele_refcall.vcf > test/tmp/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.bed || die "${testname} failed"
diff -q test/tmp/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.bed test/expected/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.bed || die "${testname} failed: output does not match expected output"


# HG002 PromethION - insertion with multiple CpG sites (T->TGCCGCGCGCGCAC at chr22:15689408), validates offset per CpG
testname="varfreq hg002_prom_chr22_ins_multi_cg"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" test/tmp/genome_chr22.fa test/data/hg002_prom_chr22_ins_multi_cg.bam test/data/hg002_prom_chr22_ins_multi_cg.vcf > test/tmp/hg002_prom_chr22_ins_multi_cg.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/hg002_prom_chr22_ins_multi_cg.mm.varfreq.tsv test/expected/hg002_prom_chr22_ins_multi_cg.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq -b hg002_prom_chr22_ins_multi_cg"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/hg002_prom_chr22_ins_multi_cg.bam test/data/hg002_prom_chr22_ins_multi_cg.vcf > test/tmp/hg002_prom_chr22_ins_multi_cg.mm.varfreq.bed || die "${testname} failed"
diff -q test/tmp/hg002_prom_chr22_ins_multi_cg.mm.varfreq.bed test/expected/hg002_prom_chr22_ins_multi_cg.mm.varfreq.bed || die "${testname} failed: output does not match expected output"

# phased VCF: example-ont reads with alternating HP=1/HP=2 + VCF mixing 1|1, 0|1, 1|0, 0|0, 0/1.
# Exercises hom-alt (any hap), per-haplotype filtering, hom-ref skip, unphased fallthrough.
testname="varview example-ont-phased"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont-hpmix.bam test/data/example-ont-phased.vcf > test/tmp/example-ont-phased.mm.varview.tsv || die "${testname} failed"
diff -q test/tmp/example-ont-phased.mm.varview.tsv test/expected/example-ont-phased.mm.varview.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq --haplotypes example-ont-phased"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq --haplotypes -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont-hpmix.bam test/data/example-ont-phased.vcf > test/tmp/example-ont-phased.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/example-ont-phased.mm.varfreq.tsv test/expected/example-ont-phased.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

# Synthetic compound-variant test: two adjacent phased SNPs together create a new CpG that
# neither produces alone. Verifies the per-hap cluster pass + is_compound scan-time gating.
testname="varview --haplotypes example-ont-compound"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview --haplotypes -c "m" test/tmp/genome_chr22.fa test/data/example-ont-compound.bam test/data/example-ont-compound.vcf > test/tmp/example-ont-compound.m.varview.tsv || die "${testname} failed"
diff -q test/tmp/example-ont-compound.m.varview.tsv test/expected/example-ont-compound.m.varview.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq --haplotypes example-ont-compound"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq --haplotypes -b -c "m" test/tmp/genome_chr22.fa test/data/example-ont-compound.bam test/data/example-ont-compound.vcf > test/tmp/example-ont-compound.m.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/example-ont-compound.m.varfreq.bedmethyl test/expected/example-ont-compound.m.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

# Indel compound test: an insertion (A->AC) followed by a SNP (T->G) on the same hap. Neither
# variant alone produces a CpG; together the inserted C and the SNP G form one. Verifies
# indel-aware cluster patching and the insertion-style scan-time match.
testname="varview --haplotypes example-ont-compound-indel"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview --haplotypes -c "m" test/tmp/genome_chr22.fa test/data/example-ont-compound-indel.bam test/data/example-ont-compound-indel.vcf > test/tmp/example-ont-compound-indel.m.varview.tsv || die "${testname} failed"
diff -q test/tmp/example-ont-compound-indel.m.varview.tsv test/expected/example-ont-compound-indel.m.varview.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq --haplotypes example-ont-compound-indel"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq --haplotypes -b -c "m" test/tmp/genome_chr22.fa test/data/example-ont-compound-indel.bam test/data/example-ont-compound-indel.vcf > test/tmp/example-ont-compound-indel.m.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/example-ont-compound-indel.m.varfreq.bedmethyl test/expected/example-ont-compound-indel.m.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

echo -e "${GREEN}All tests passed!${NC}"
