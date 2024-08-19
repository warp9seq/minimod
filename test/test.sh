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

if [ ! -f test/tmp/truth.tsv ]; then
    wget  -N -O test/tmp/truth.tsv "https://raw.githubusercontent.com/imsuneth/shared-files/main/truth.tsv" || die "Downloading the truthset failed"
fi

exp_corr=0.87 # update this if the expected correlation changes
testname="Accuracy Test: mod-freq results correlation with truthset"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -t 8 -b test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/accu.bedmethyl  || die "${testname} Running the tool failed"
corr=`./test/compare.py test/tmp/truth.tsv test/tmp/accu.bedmethyl`
if (( $(echo "$corr >= $exp_corr" | bc -l) )); then
    echo -e "${GREEN}Corr: $corr\tExpected: $exp_corr\tPassed${NC}\n"
elif (( $(echo "$exp_corr > $corr" | bc -l) )); then
    echo -e "${RED}Corr: $corr\tExpected: $exp_corr\tDecreased${NC}\n"
    die "${testname} Correlation decreased"
fi

testname="Test 1: view hifi"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -m 0.2 -t 8 test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test1.tsv  || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test1.tsv > test/tmp/test1.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test1.tsv > test/tmp/test1.tsv.sorted
diff -q test/tmp/test1.exp.tsv.sorted test/tmp/test1.tsv.sorted || die "${testname} diff failed"

testname="Test 2: view ont"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -m 0.0 -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2.tsv > test/tmp/test2.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2.tsv > test/tmp/test2.tsv.sorted
diff -q test/tmp/test2.exp.tsv.sorted test/tmp/test2.tsv.sorted || die "${testname} diff failed"

testname="Test 3: mod-freq hifi"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -t 8 test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test3.tsv  || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test3.tsv > test/tmp/test3.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test3.tsv > test/tmp/test3.tsv.sorted
diff -q test/tmp/test3.exp.tsv.sorted test/tmp/test3.tsv.sorted || die "${testname} diff failed"

testname="Test 4: mod-freq hifi bedmethyl output batch size 1"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -b -t 8 -K 1 test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test4.bedmethyl  || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k6,6 test/expected/test4.bedmethyl > test/tmp/test4.bedmethyl.exp.sorted
sort -k1,1 -k2,2n -k6,6 test/tmp/test4.bedmethyl > test/tmp/test4.bedmethyl.sorted
diff -q test/tmp/test4.bedmethyl.exp.sorted test/tmp/test4.bedmethyl.sorted || die "${testname} diff failed"

testname="Test 5: mod-freq ont"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test5.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test5.tsv > test/tmp/test5.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test5.tsv > test/tmp/test5.tsv.sorted
diff -q test/tmp/test5.exp.tsv.sorted test/tmp/test5.tsv.sorted || die "${testname} diff failed"

testname="Test 6: mod-freq ont bedmethyl output"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -b -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test6.bedmethyl || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k6,6 test/expected/test6.bedmethyl > test/tmp/test6.bedmethyl.exp.sorted
sort -k1,1 -k2,2n -k6,6 test/tmp/test6.bedmethyl > test/tmp/test6.bedmethyl.sorted
diff -q test/tmp/test6.bedmethyl.exp.sorted test/tmp/test6.bedmethyl.sorted || die "${testname} diff failed"

testname="Test 7: mod-freq ont with mod threshold"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -m 0.2 -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test7.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test7.tsv > test/tmp/test7.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test7.tsv > test/tmp/test7.tsv.sorted
diff -q test/tmp/test7.exp.tsv.sorted test/tmp/test7.tsv.sorted || die "${testname} diff failed"

testname="Test 8: mod-freq ont with mod codes m and h"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -c "mh" -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test8.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test8.tsv > test/tmp/test8.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test8.tsv > test/tmp/test8.tsv.sorted
diff -q test/tmp/test8.exp.tsv.sorted test/tmp/test8.tsv.sorted || die "${testname} diff failed"

testname="Test 9: mod-freq ont with mod codes h only"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -c "h" -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test9.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test9.tsv > test/tmp/test9.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test9.tsv > test/tmp/test9.tsv.sorted
diff -q test/tmp/test9.exp.tsv.sorted test/tmp/test9.tsv.sorted || die "${testname} diff failed"

testname="Test 10: view ont with mod threshold"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -m 0.2 -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test10.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test10.tsv > test/tmp/test10.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test10.tsv > test/tmp/test10.tsv.sorted
diff -q test/tmp/test10.exp.tsv.sorted test/tmp/test10.tsv.sorted || die "${testname} diff failed"

testname="Test 11: view ont with mod codes m and h with different thresholds"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c "mh" -m "0.2,0.5"  -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test11.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test11.tsv > test/tmp/test11.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test11.tsv > test/tmp/test11.tsv.sorted
diff -q test/tmp/test11.exp.tsv.sorted test/tmp/test11.tsv.sorted || die "${testname} diff failed"

testname="Test 12: mod-freq ont with mod codes m and h with different thresholds"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -c "mh" -m "0.2,0.5" -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test12.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test12.tsv > test/tmp/test12.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test12.tsv > test/tmp/test12.tsv.sorted
diff -q test/tmp/test12.exp.tsv.sorted test/tmp/test12.tsv.sorted || die "${testname} diff failed"

testname="Test 13: view ont with -o flag"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -m 0.0 -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam -o test/tmp/test13.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2.tsv > test/tmp/test13.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test13.tsv > test/tmp/test13.tsv.sorted
diff -q test/tmp/test13.exp.tsv.sorted test/tmp/test13.tsv.sorted || die "${testname} diff failed"

testname="Test 14: mod-freq ont with -o flag"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod mod-freq -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam -o test/tmp/test14.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/tmp/test14.tsv > test/tmp/test14.tsv.sorted
diff -q test/tmp/test5.exp.tsv.sorted test/tmp/test14.tsv.sorted || die "${testname} diff failed"

# ======= Extensive tests (prerequisits: buttery-eel, minimap2) - tested on gtgpu =======
# blow5=/home/hasindu/scratch/hg2_prom_lsk114_5khz/chr22/PGXXXX230339_reads_chr22.blow5
# guppybin=/data/suneth/tools/ont-dorado-server/bin
# model=dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac.cfg
# ref=/genome/hg38noAlt.fa
# refidx=/genome/hg38noAlt.idx

# buttery-eel -g $guppybin --port 5000 --use_tcp --device cuda:all --call_mods --config $model -i $blow5 -o "reads.unaligned.sam" || die "${testname} Running buttery-eel failed"
# samtools fastq -@ 32 -TMM,ML "reads.unaligned.sam" | minimap2 -t 32 -x map-ont --sam-hit-only -Y -a -y --secondary=no $refidx - | samtools sort -@ 32 - > "reads.bam" || die "${testname} Running samtools, minimap2 failed"
# samtools index -@ 32 "reads.bam" || die "${testname} Running samtools index failed"

# testname="Ext. Test 1: view extensive test using PGXXXX230339_reads_chr22"
# expected=/data/suneth/do_not_delete/ext_test_minimod/ext.test1.tsv
# echo -e "${BLUE}${testname}${NC}"
# ex  ./minimod view -m 0.0 -t 32 $ref reads.bam > test/tmp/ext.test1.tsv  || die "${testname} Running the tool failed"
# diff -q $expected test/tmp/ext.test1.tsv || die "${testname} diff failed"

# testname="Ext. Test 2: mod-freq extensive test using PGXXXX230339_reads_chr22"
# expected=/data/suneth/do_not_delete/ext_test_minimod/ext.test2.tsv
# echo -e "${BLUE}${testname}${NC}"
# # ex  ./minimod mod-freq -t 32 $ref reads.bam > test/tmp/ext.test2.tsv || die "${testname} Running the tool failed"
# sort -k1,1 -k2,2n -k4,4 $expected > test/expected/ext.test2.tsv.sorted
# sort -k1,1 -k2,2n -k4,4 test/tmp/ext.test2.tsv > test/tmp/ext.test2.tsv.sorted
# diff -q test/expected/ext.test2.tsv.sorted test/tmp/ext.test2.tsv.sorted || die "${testname} diff failed"
# ============================================================================================

echo -e "${GREEN}ALL TESTS PASSED !${NC}"