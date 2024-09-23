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

# ======= Extensive tests (prerequisits: buttery-eel, minimap2) - tested on gtgpu =======
# blow5=/home/hasindu/scratch/hg2_prom_lsk114_5khz/chr22/PGXXXX230339_reads_chr22.blow5
# guppybin=/data/suneth/tools/ont-dorado-server/bin
# model=dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac.cfg
# ref=/genome/hg38noAlt.fa
# refidx=/genome/hg38noAlt.idx

# buttery-eel -g $guppybin --port 5000 --use_tcp --device cuda:all --call_mods --config $model -i $blow5 -o "test/tmp/reads.unaligned.sam" || die "${testname} Running buttery-eel failed"
# samtools fastq -@ 32 -TMM,ML "test/tmp/reads.unaligned.sam" | minimap2 -t 32 -x map-ont --sam-hit-only -Y -a -y --secondary=no $refidx - | samtools sort -@ 32 - > "reads.bam" || die "${testname} Running samtools, minimap2 failed"
# samtools index -@ 32 "test/tmp/reads.bam" || die "${testname} Running samtools index failed"

# testname="Ext. Test 1: view extensive test using PGXXXX230339_reads_chr22"
# expected=/data/suneth/do_not_delete/ext_test_minimod/ext.test1.tsv
# echo -e "${BLUE}${testname}${NC}"
# ex  ./minimod view -t 32 $ref /data/suneth/tool-validation/minimod/PGXXXX230339_reads_chr22/PGXXXX230339_reads_chr22.remora.bam > test/tmp/ext.test1.tsv  || die "${testname} Running the tool failed"
# diff -q $expected test/tmp/ext.test1.tsv || die "${testname} diff failed"

testname="Ext. Test 1: mod-freq extensive test using PGXXXX230339_reads_chr22"
echo -e "${BLUE}${testname}${NC}"
exp_corr=0.90
ex  ./minimod mod-freq -b -t 32 /genome/hg38noAlt.fa /data/suneth/tool-validation/minimod/PGXXXX230339_reads_chr22/PGXXXX230339_reads_chr22.remora.bam > test/tmp/ext.test1.bedmethyl || die "${testname} Running the tool failed"
corr=`./test/compare.py /data/suneth/tool-validation/hg2_chr22_bi.tsv test/tmp/ext.test1.bedmethyl`
if (( $(echo "$corr >= $exp_corr" | bc -l) )); then
    echo -e "${GREEN}Corr: $corr\tExpected: $exp_corr\tPassed${NC}\n"
elif (( $(echo "$exp_corr > $corr" | bc -l) )); then
    echo -e "${RED}Corr: $corr\tExpected: $exp_corr\tDecreased${NC}\n"
    die "${testname} Correlation decreased"
fi

testname="Ext. Test 2: mod-freq extensive test using na12878_prom_lsk114"
echo -e "${BLUE}${testname}${NC}"
exp_corr=0.90
ex  ./minimod mod-freq -b -t 32 /genome/hg38noAlt.fa /home/hasindu/scratch/na12878_prom_lsk114/compare-bisulfite/remora_with_supple/remora_mapped.bam > test/tmp/ext.test2.bedmethyl || die "${testname} Running the tool failed"
corr=`./test/compare.py /data/suneth/tool-validation/hg2_chr22_bi.tsv test/tmp/ext.test2.bedmethyl`
if (( $(echo "$corr >= $exp_corr" | bc -l) )); then
    echo -e "${GREEN}Corr: $corr\tExpected: $exp_corr\tPassed${NC}\n"
elif (( $(echo "$exp_corr > $corr" | bc -l) )); then
    echo -e "${RED}Corr: $corr\tExpected: $exp_corr\tDecreased${NC}\n"
    die "${testname} Correlation decreased"
fi
# ============================================================================================

echo -e "${GREEN}ALL TESTS PASSED !${NC}"