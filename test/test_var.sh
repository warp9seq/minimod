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

mkdir -p test/tmp test/tmp/varmod || die "Creating the tmp directory failed"

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
ex ./minimod varview -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont.bam test/data/varmod/example-ont-clair.vcf > test/tmp/varmod/example-ont.mm.varview.tsv || die "${testname} failed"
diff -q test/tmp/varmod/example-ont.mm.varview.tsv test/expected/varmod/example-ont.mm.varview.tsv || die "${testname} failed: output does not match expected output"

testname="varview -b example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont.bam test/data/varmod/example-ont-clair.vcf > test/tmp/varmod/example-ont.mm.varview.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/example-ont.mm.varview.bedmethyl test/expected/varmod/example-ont.mm.varview.bedmethyl || die "${testname} failed: output does not match expected output"

testname="varview dna_4mC_5mC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/data/varmod/dna_4mC_5mC_mm_chr22.vcf > test/tmp/varmod/dna_4mC_5mC_mm_chr22.mm.varview.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/dna_4mC_5mC_mm_chr22.mm.varview.bedmethyl test/expected/varmod/dna_4mC_5mC_mm_chr22.mm.varview.bedmethyl || die "${testname} failed: output does not match expected output"

testname="varview dna_5mCG_5hmCG_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/data/varmod/dna_5mCG_5hmCG_mm_chr22.vcf > test/tmp/varmod/dna_5mCG_5hmCG_mm_chr22.mm.varview.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/dna_5mCG_5hmCG_mm_chr22.mm.varview.bedmethyl test/expected/varmod/dna_5mCG_5hmCG_mm_chr22.mm.varview.bedmethyl || die "${testname} failed: output does not match expected output"

testname="varview dna_5mC_5hmC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/data/varmod/dna_5mC_5hmC_mm_chr22.vcf > test/tmp/varmod/dna_5mC_5hmC_mm_chr22.mm.varview.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/dna_5mC_5hmC_mm_chr22.mm.varview.bedmethyl test/expected/varmod/dna_5mC_5hmC_mm_chr22.mm.varview.bedmethyl || die "${testname} failed: output does not match expected output"

# varfreq tests

testname="varfreq example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont.bam test/data/varmod/example-ont-clair.vcf > test/tmp/varmod/example-ont.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/varmod/example-ont.mm.varfreq.tsv test/expected/varmod/example-ont.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq -b example-ont"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/example-ont.bam test/data/varmod/example-ont-clair.vcf > test/tmp/varmod/example-ont.mm.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/example-ont.mm.varfreq.bedmethyl test/expected/varmod/example-ont.mm.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

testname="varfreq -b dna_4mC_5mC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/data/varmod/dna_4mC_5mC_mm_chr22.vcf > test/tmp/varmod/dna_4mC_5mC_mm_chr22.mm.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/dna_4mC_5mC_mm_chr22.mm.varfreq.bedmethyl test/expected/varmod/dna_4mC_5mC_mm_chr22.mm.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

testname="varfreq -b dna_5mCG_5hmCG_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/data/varmod/dna_5mCG_5hmCG_mm_chr22.vcf > test/tmp/varmod/dna_5mCG_5hmCG_mm_chr22.mm.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/dna_5mCG_5hmCG_mm_chr22.mm.varfreq.bedmethyl test/expected/varmod/dna_5mCG_5hmCG_mm_chr22.mm.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

testname="varfreq -b dna_5mC_5hmC_mm_chr22"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/data/varmod/dna_5mC_5hmC_mm_chr22.vcf > test/tmp/varmod/dna_5mC_5hmC_mm_chr22.mm.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/dna_5mC_5hmC_mm_chr22.mm.varfreq.bedmethyl test/expected/varmod/dna_5mC_5hmC_mm_chr22.mm.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

# hg002 PGXXSX240041 promethion chr22 - covers all variant types: SNP, INS, DEL, multiallelic, RefCall (chr22:17280000-17380000)
testname="varfreq hg002_prom_chr22_snp_ins_del_multiallele_refcall"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" test/tmp/genome_chr22.fa test/data/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.bam test/data/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.vcf > test/tmp/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.tsv test/expected/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq -b hg002_prom_chr22_snp_ins_del_multiallele_refcall"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.bam test/data/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.vcf > test/tmp/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.bedmethyl test/expected/varmod/hg002_prom_chr22_snp_ins_del_multiallele_refcall.mm.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"


# HG002 PromethION - insertion with multiple CpG sites (T->TGCCGCGCGCGCAC at chr22:15689408), validates offset per CpG
testname="varfreq hg002_prom_chr22_ins_multi_cg"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -c "m,h" test/tmp/genome_chr22.fa test/data/varmod/hg002_prom_chr22_ins_multi_cg.bam test/data/varmod/hg002_prom_chr22_ins_multi_cg.vcf > test/tmp/varmod/hg002_prom_chr22_ins_multi_cg.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/varmod/hg002_prom_chr22_ins_multi_cg.mm.varfreq.tsv test/expected/varmod/hg002_prom_chr22_ins_multi_cg.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq -b hg002_prom_chr22_ins_multi_cg"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/varmod/hg002_prom_chr22_ins_multi_cg.bam test/data/varmod/hg002_prom_chr22_ins_multi_cg.vcf > test/tmp/varmod/hg002_prom_chr22_ins_multi_cg.mm.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/hg002_prom_chr22_ins_multi_cg.mm.varfreq.bedmethyl test/expected/varmod/hg002_prom_chr22_ins_multi_cg.mm.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

# phased VCF: example-ont reads with alternating HP=1/HP=2 + VCF mixing 1|1, 0|1, 1|0, 0|0, 0/1.
# Exercises hom-alt (any hap), per-haplotype filtering, hom-ref skip, unphased fallthrough.
testname="varview example-ont-phased"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview -c "m,h" test/tmp/genome_chr22.fa test/data/varmod/example-ont-hpmix.bam test/data/varmod/example-ont-phased.vcf > test/tmp/varmod/example-ont-phased.mm.varview.tsv || die "${testname} failed"
diff -q test/tmp/varmod/example-ont-phased.mm.varview.tsv test/expected/varmod/example-ont-phased.mm.varview.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq --haplotypes example-ont-phased"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq --haplotypes -c "m,h" test/tmp/genome_chr22.fa test/data/varmod/example-ont-hpmix.bam test/data/varmod/example-ont-phased.vcf > test/tmp/varmod/example-ont-phased.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/varmod/example-ont-phased.mm.varfreq.tsv test/expected/varmod/example-ont-phased.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

# same VCF as above but with the ID column filled in: a plain rsID, a multi-ID (';' separated),
# a missing ID ('.') and a non-rs ID. Checks var_id is carried through to the output.
testname="varfreq --haplotypes example-ont-rsid"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq --haplotypes -c "m,h" test/tmp/genome_chr22.fa test/data/varmod/example-ont-hpmix.bam test/data/varmod/example-ont-rsid.vcf > test/tmp/varmod/example-ont-rsid.mm.varfreq.tsv || die "${testname} failed"
diff -q test/tmp/varmod/example-ont-rsid.mm.varfreq.tsv test/expected/varmod/example-ont-rsid.mm.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq -b example-ont-rsid"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -c "m,h" test/tmp/genome_chr22.fa test/data/varmod/example-ont-hpmix.bam test/data/varmod/example-ont-rsid.vcf > test/tmp/varmod/example-ont-rsid.mm.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/example-ont-rsid.mm.varfreq.bedmethyl test/expected/varmod/example-ont-rsid.mm.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

# Synthetic compound-variant test: two adjacent phased SNPs together create a new CpG that
# neither produces alone. Verifies the per-hap cluster pass + is_compound scan-time gating.
testname="varview --haplotypes example-ont-compound"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview --haplotypes -c "m" test/tmp/genome_chr22.fa test/data/varmod/example-ont-compound.bam test/data/varmod/example-ont-compound.vcf > test/tmp/varmod/example-ont-compound.m.varview.tsv || die "${testname} failed"
diff -q test/tmp/varmod/example-ont-compound.m.varview.tsv test/expected/varmod/example-ont-compound.m.varview.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq --haplotypes example-ont-compound"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq --haplotypes -b -c "m" test/tmp/genome_chr22.fa test/data/varmod/example-ont-compound.bam test/data/varmod/example-ont-compound.vcf > test/tmp/varmod/example-ont-compound.m.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/example-ont-compound.m.varfreq.bedmethyl test/expected/varmod/example-ont-compound.m.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

# Indel compound test: an insertion (A->AC) followed by a SNP (T->G) on the same hap. Neither
# variant alone produces a CpG; together the inserted C and the SNP G form one. Verifies
# indel-aware cluster patching and the insertion-style scan-time match.
testname="varview --haplotypes example-ont-compound-indel"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varview --haplotypes -c "m" test/tmp/genome_chr22.fa test/data/varmod/example-ont-compound-indel.bam test/data/varmod/example-ont-compound-indel.vcf > test/tmp/varmod/example-ont-compound-indel.m.varview.tsv || die "${testname} failed"
diff -q test/tmp/varmod/example-ont-compound-indel.m.varview.tsv test/expected/varmod/example-ont-compound-indel.m.varview.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq --haplotypes example-ont-compound-indel"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq --haplotypes -b -c "m" test/tmp/genome_chr22.fa test/data/varmod/example-ont-compound-indel.bam test/data/varmod/example-ont-compound-indel.vcf > test/tmp/varmod/example-ont-compound-indel.m.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/example-ont-compound-indel.m.varfreq.bedmethyl test/expected/varmod/example-ont-compound-indel.m.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

# HG002 PGXXXX250050 — a merged chr1 fixture (6 slices) covering, with >=2 entries each:
#   cpg_gain() classes  : from_SNP, from_DEL, from_INS (offset<=0), within_INS (offset>0)
#     * chr1:1029600-1030643 and chr1:2122900-2123550 each carry all four classes
#   hap-based merge     : two adjacent phased SNPs on the SAME hap that jointly form a CpG
#     * chr1:104963183-104963184 (1|0) and chr1:119585475-119585476 (1|0) — CpG on hap1 only
#   MNP                 : chr1:36689056 (TA->CG) and chr1:73471034 (AA->CG). clair3 decomposes
#     MNPs into adjacent SNPs, so these two are re-composed into single MNP records in the VCF
#     (the reads carry the same haplotype, so the reconstructed CpG is identical).
#   variant types       : SNP, INS, DEL, MNP, REFCALL (0/0)  — all present >=2x in the VCF
#   phasing             : both phased (1|0/0|1) and unphased (1/1) genotypes present
# freq output is trimmed to the union of the 6 slice windows (long reads extend far beyond them).
FREQ_WIN='($1=="chr1")&&(($2>=1029600&&$2<=1030643)||($2>=2122900&&$2<=2123550)||($2>=36688756&&$2<=36689356)||($2>=73470734&&$2<=73471334)||($2>=104962934&&$2<=104963434)||($2>=119585226&&$2<=119585726))'
testname="freq --haplotypes varmod_all_chr1"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod freq --haplotypes -b -c "m,h" test/tmp/genome_chr1.fa test/data/varmod/varmod_all_chr1.bam > test/tmp/varmod/varmod_all_chr1.mh.freq.full.bedmethyl || die "${testname} failed"
awk -F'\t' "$FREQ_WIN" test/tmp/varmod/varmod_all_chr1.mh.freq.full.bedmethyl > test/tmp/varmod/varmod_all_chr1.mh.freq.bedmethyl
diff -q test/tmp/varmod/varmod_all_chr1.mh.freq.bedmethyl test/expected/varmod/varmod_all_chr1.mh.freq.bedmethyl || die "${testname} failed: output does not match expected output"

testname="varfreq --haplotypes varmod_all_chr1"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq --haplotypes -b -c "m,h" test/tmp/genome_chr1.fa test/data/varmod/varmod_all_chr1.bam test/data/varmod/varmod_all_chr1.vcf > test/tmp/varmod/varmod_all_chr1.mh.varfreq.bedmethyl || die "${testname} failed"
diff -q test/tmp/varmod/varmod_all_chr1.mh.varfreq.bedmethyl test/expected/varmod/varmod_all_chr1.mh.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

# Analysis-script tests: run the varfreq_summary/varfreq_context helpers over the
# varmod_all_chr1 freq + varfreq outputs generated above and diff their reports.
# varfreq_summary also takes the VCF and the reference (for the reference-CpG leak check).
testname="varfreq_summary varmod_all_chr1"
echo -e "${BLUE}${testname}${NC}"
ex test/varfreq_summary.py test/data/varmod/varmod_all_chr1.vcf test/tmp/varmod/varmod_all_chr1.mh.freq.bedmethyl test/tmp/varmod/varmod_all_chr1.mh.varfreq.bedmethyl test/tmp/genome_chr1.fa > test/tmp/varmod/varmod_all_chr1.mh.varfreq.summary.txt || die "${testname} failed"
# diff -q test/tmp/varmod/varmod_all_chr1.mh.varfreq.summary.txt test/expected/varmod/varmod_all_chr1.mh.varfreq.summary.txt || die "${testname} failed: output does not match expected output"

testname="varfreq_context varmod_all_chr1"
echo -e "${BLUE}${testname}${NC}"
ex test/varfreq_context.py test/tmp/varmod/varmod_all_chr1.mh.freq.bedmethyl test/tmp/varmod/varmod_all_chr1.mh.varfreq.bedmethyl 1000 > test/tmp/varmod/varmod_all_chr1.mh.varfreq.context.tsv 2>/dev/null || die "${testname} failed"
# diff -q test/tmp/varmod/varmod_all_chr1.mh.varfreq.context.tsv test/expected/varmod/varmod_all_chr1.mh.varfreq.context.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq_context --bed varmod_all_chr1"
echo -e "${BLUE}${testname}${NC}"
ex test/varfreq_context.py test/tmp/varmod/varmod_all_chr1.mh.freq.bedmethyl test/tmp/varmod/varmod_all_chr1.mh.varfreq.bedmethyl 1000 --bed test/tmp/varmod/varmod_all_chr1.mh.varfreq.context.bed > /dev/null 2>&1 || die "${testname} failed"
# diff -q test/tmp/varmod/varmod_all_chr1.mh.varfreq.context.bed test/expected/varmod/varmod_all_chr1.mh.varfreq.context.bed || die "${testname} failed: output does not match expected output"

# Region-based regression tests. Each entry in test/regions/regions.tsv is sliced
# by test/regions/extract.sh into a BAM + VCF (original chrom names + absolute
# coordinates preserved, so the fixtures load directly in IGV). The reference
# at test time is the existing per-chromosome download.
while IFS=$'\t' read -r label chrom start end notes; do
    [[ "$label" =~ ^#.* ]] && continue
    [ -z "$label" ] && continue

    case "$chrom" in
        chr1)  ref=test/tmp/genome_chr1.fa ;;
        chr22) ref=test/tmp/genome_chr22.fa ;;
        *)     die "region:${label} chrom=${chrom} has no reference download wired up" ;;
    esac
    bam=test/data/varmod/region_${label}.bam
    vcf=test/data/varmod/region_${label}.vcf

    testname="varview -b --haplotypes region:${label}"
    echo -e "${BLUE}${testname}${NC}"
    ex ./minimod varview --haplotypes -b -c "m,h" "$ref" "$bam" "$vcf" > test/tmp/varmod/region_${label}.mh.varview.bedmethyl || die "${testname} failed"
    diff -q test/tmp/varmod/region_${label}.mh.varview.bedmethyl test/expected/varmod/region_${label}.mh.varview.bedmethyl || die "${testname} failed: output does not match expected output"

    testname="varfreq -b --haplotypes region:${label}"
    echo -e "${BLUE}${testname}${NC}"
    ex ./minimod varfreq --haplotypes -b -c "m,h" "$ref" "$bam" "$vcf" > test/tmp/varmod/region_${label}.mh.varfreq.bedmethyl || die "${testname} failed"
    diff -q test/tmp/varmod/region_${label}.mh.varfreq.bedmethyl test/expected/varmod/region_${label}.mh.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"
done < test/regions/regions.tsv

echo -e "${GREEN}All tests passed!${NC}"
