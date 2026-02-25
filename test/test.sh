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

#**** OLD TESTS - will be removed in future ****

testname="Accuracy Test: freq results correlation with modkit and truthset"
exp_modkit_corr=0.97 # update this if the expected correlation changes
exp_truth_corr=0.85 # update this if the expected correlation changes
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -t 8 -b test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/accu.freq.bedmethyl  || die "${testname} Running the tool failed"
corr=`./test/compare.py test/data/accu.mk.pileup.bedmethyl test/tmp/accu.freq.bedmethyl`
if (( $(echo "$corr >= $exp_modkit_corr" | bc -l) )); then
    echo -e "${GREEN}Corr with modkit: $corr\tExpected: $exp_modkit_corr\tPassed${NC}\n"
elif (( $(echo "$exp_modkit_corr > $corr" | bc -l) )); then
    echo -e "${RED}Corr with modkit: $corr\tExpected: $exp_modkit_corr\tDecreased${NC}\n"
    die "${testname} Correlation with modkit decreased"
fi

corr=`./test/compare.py test/tmp/truth.tsv test/tmp/accu.freq.bedmethyl`
if (( $(echo "$corr >= $exp_truth_corr" | bc -l) )); then
    echo -e "${GREEN}Corr with truthset: $corr\tExpected: $exp_truth_corr\tPassed${NC}\n"
elif (( $(echo "$exp_truth_corr > $corr" | bc -l) )); then
    echo -e "${RED}Corr with truthset: $corr\tExpected: $exp_truth_corr\tDecreased${NC}\n"
    die "${testname} Correlation with truthset decreased"
fi

# testname="Accuracy Test: minimod view vs modkit extract full"
# ex ./minimod view -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/accu.view.tsv || die "${testname} Running the tool failed"
# mkdir -p test/tmp/view_compare
# ./test/compare_view_mkbed_mmtsv.sh -y test/data/accu.mk.extract.bedmethyl test/tmp/accu.view.tsv test/tmp/view_compare || die "${testname} Comparison failed"

# missing1_lines=$(wc -l < test/tmp/view_compare/missing_in_file2.tsv)
# if [ "$missing1_lines" -gt 0 ]; then
#     die "${testname} minimod view missing records compared to modkit extract full"
# fi

testname="Test 1: view hifi"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[CG] -t 8 test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test1.tsv  || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test1.tsv > test/tmp/test1.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test1.tsv > test/tmp/test1.tsv.sorted
diff -q test/tmp/test1.exp.tsv.sorted test/tmp/test1.tsv.sorted || die "${testname} diff failed"

testname="Test 2: view ont"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[CG] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2.tsv > test/tmp/test2.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2.tsv > test/tmp/test2.tsv.sorted
diff -q test/tmp/test2.exp.tsv.sorted test/tmp/test2.tsv.sorted || die "${testname} diff failed"

testname="Test 2a: view ont with insertions"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[CG] --insertions test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2a.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2a.tsv > test/tmp/test2a.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2a.tsv > test/tmp/test2a.tsv.sorted
diff -q test/tmp/test2a.exp.tsv.sorted test/tmp/test2a.tsv.sorted || die "${testname} diff failed"

testname="Test 2b: view ont with all contexts"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[*] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2b.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2b.tsv > test/tmp/test2b.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2b.tsv > test/tmp/test2b.tsv.sorted
diff -q test/tmp/test2b.exp.tsv.sorted test/tmp/test2b.tsv.sorted || die "${testname} diff failed"

testname="Test 2c: view ont all mod codes with wildcard"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "*" test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2c_wild.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2c_wild.tsv > test/tmp/test2c_wild.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2c_wild.tsv > test/tmp/test2c_wild.tsv.sorted
diff -q test/tmp/test2c_wild.exp.tsv.sorted test/tmp/test2c_wild.tsv.sorted || die "${testname} diff failed"

testname="Test 2c: view ont with haplotypes"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[CG] --haplotypes test/tmp/genome_chr1.fa test/data/hap.bam > test/tmp/test2c.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2c.tsv > test/tmp/test2c.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2c.tsv > test/tmp/test2c.tsv.sorted
diff -q test/tmp/test2c.exp.tsv.sorted test/tmp/test2c.tsv.sorted || die "${testname} diff failed"

testname="Test 2d: view ont with U in reference and t in required context"
echo -e "${BLUE}${testname}${NC}"
sed 's/T/U/g' test/tmp/genome_chr22.fa > test/tmp/genome_chr22_U.fa # replace all T with U in chr22 reference
ex  ./minimod view -t 8 -c m[Ct] test/tmp/genome_chr22_U.fa test/data/example-ont.bam > test/tmp/test2d_U.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2d_U.tsv > test/tmp/test2d_U.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2d_U.tsv > test/tmp/test2d_U.tsv.sorted
diff -q test/tmp/test2d_U.tsv.sorted test/tmp/test2d_U.tsv.sorted || die "${testname} diff failed"

testname="Test 3: freq hifi"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -t 8 test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test3.tsv  || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test3.tsv > test/tmp/test3.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test3.tsv > test/tmp/test3.tsv.sorted
diff -q test/tmp/test3.exp.tsv.sorted test/tmp/test3.tsv.sorted || die "${testname} diff failed"

testname="Test 4: freq hifi bedmethyl output batch size 1"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -b -t 8 -K 1 test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test4.bedmethyl  || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k6,6 test/expected/test4.bedmethyl > test/tmp/test4.bedmethyl.exp.sorted
sort -k1,1 -k2,2n -k6,6 test/tmp/test4.bedmethyl > test/tmp/test4.bedmethyl.sorted
diff -q test/tmp/test4.bedmethyl.exp.sorted test/tmp/test4.bedmethyl.sorted || die "${testname} diff failed"

testname="Test 5: freq ont"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test5.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test5.tsv > test/tmp/test5.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test5.tsv > test/tmp/test5.tsv.sorted
diff -q test/tmp/test5.exp.tsv.sorted test/tmp/test5.tsv.sorted || die "${testname} diff failed"

testname="Test 5a: freq ont with insertions"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -t 8 --insertions test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test5a.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test5a.tsv > test/tmp/test5a.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test5a.tsv > test/tmp/test5a.tsv.sorted
diff -q test/tmp/test5a.exp.tsv.sorted test/tmp/test5a.tsv.sorted || die "${testname} diff failed"

testname="Test 5b: freq ont with all contexts"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -t 8 -c m[*] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test5b.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test5b.tsv > test/tmp/test5b.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test5b.tsv > test/tmp/test5b.tsv.sorted
diff -q test/tmp/test5b.exp.tsv.sorted test/tmp/test5b.tsv.sorted || die "${testname} diff failed"

testname="Test 5c: freq ont with haplotypes"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -t 8 --haplotypes test/tmp/genome_chr1.fa test/data/hap.bam > test/tmp/test5c.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test5c.tsv > test/tmp/test5c.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test5c.tsv > test/tmp/test5c.tsv.sorted
diff -q test/tmp/test5c.exp.tsv.sorted test/tmp/test5c.tsv.sorted || die "${testname} diff failed"

testname="Test 6: freq ont bedmethyl output"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -b -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test6.bedmethyl || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k6,6 test/expected/test6.bedmethyl > test/tmp/test6.bedmethyl.exp.sorted
sort -k1,1 -k2,2n -k6,6 test/tmp/test6.bedmethyl > test/tmp/test6.bedmethyl.sorted
diff -q test/tmp/test6.bedmethyl.exp.sorted test/tmp/test6.bedmethyl.sorted || die "${testname} diff failed"

testname="Test 7: freq ont with mod threshold"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -m 0.8 -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test7.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test7.tsv > test/tmp/test7.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test7.tsv > test/tmp/test7.tsv.sorted
diff -q test/tmp/test7.exp.tsv.sorted test/tmp/test7.tsv.sorted || die "${testname} diff failed"

testname="Test 8: freq ont with mod codes m and h"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c "m,h" -m 0.8,0.8 -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test8.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test8.tsv > test/tmp/test8.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test8.tsv > test/tmp/test8.tsv.sorted
diff -q test/tmp/test8.exp.tsv.sorted test/tmp/test8.tsv.sorted || die "${testname} diff failed"

testname="Test 9: freq ont with mod codes h only"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c "h" -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test9.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test9.tsv > test/tmp/test9.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test9.tsv > test/tmp/test9.tsv.sorted
diff -q test/tmp/test9.exp.tsv.sorted test/tmp/test9.tsv.sorted || die "${testname} diff failed"

testname="Test 10: view ont with mod threshold"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test10.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test10.tsv > test/tmp/test10.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test10.tsv > test/tmp/test10.tsv.sorted
diff -q test/tmp/test10.exp.tsv.sorted test/tmp/test10.tsv.sorted || die "${testname} diff failed"

testname="Test 11: view ont with mod codes m and h with different thresholds"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c "m,h"  -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test11.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test11.tsv > test/tmp/test11.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test11.tsv > test/tmp/test11.tsv.sorted
diff -q test/tmp/test11.exp.tsv.sorted test/tmp/test11.tsv.sorted || die "${testname} diff failed"

testname="Test 12: freq ont with mod codes m and h with different thresholds"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c "m,h" -m "0.8,0.5" -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test12.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test12.tsv > test/tmp/test12.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test12.tsv > test/tmp/test12.tsv.sorted
diff -q test/tmp/test12.exp.tsv.sorted test/tmp/test12.tsv.sorted || die "${testname} diff failed"

testname="Test 13: view ont with -o flag"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam -o test/tmp/test13.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2.tsv > test/tmp/test13.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test13.tsv > test/tmp/test13.tsv.sorted
diff -q test/tmp/test13.exp.tsv.sorted test/tmp/test13.tsv.sorted || die "${testname} diff failed"

testname="Test 14: freq ont with -o flag"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -t 8 test/tmp/genome_chr22.fa test/data/example-ont.bam -o test/tmp/test14.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/tmp/test14.tsv > test/tmp/test14.tsv.sorted
diff -q test/tmp/test5.exp.tsv.sorted test/tmp/test14.tsv.sorted || die "${testname} diff failed"

testname="Test 15: view ont with mod codes e and b"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c e,b -t 8 test/tmp/genome_chr1.fa test/data/eb.bam > test/tmp/test15.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test15.tsv > test/tmp/test15.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test15.tsv > test/tmp/test15.tsv.sorted
diff -q test/tmp/test15.exp.tsv.sorted test/tmp/test15.tsv.sorted || die "${testname} diff failed"

testname="Test 16: freq ont with mod codes e and b"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c e,b -m 0.5 -t 8 test/tmp/genome_chr1.fa test/data/eb.bam > test/tmp/test16.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test16.tsv > test/tmp/test16.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test16.tsv > test/tmp/test16.tsv.sorted
diff -q test/tmp/test16.exp.tsv.sorted test/tmp/test16.tsv.sorted || die "${testname} diff failed"

testname="Test 17: freq ChEBI mod code test with pseudouridine (ChEBI: 17802)"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c "17802[*]" -b -t 8 test/tmp/genome_chr22.fa test/data/dRNA.bam > test/tmp/dRNA.mm.freq.17802.bedmethyl || die "${testname} Running the tool failed"
/storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit pileup --reference test/tmp/genome_chr22.fa test/data/dRNA.bam test/expected/dRNA.mk.pileup.bedmethyl || die "${testname} Running modkit pileup failed"
awk 'NR==1 || $4=="17802"' test/expected/dRNA.mk.pileup.bedmethyl > test/expected/dRNA.mk.pileup.17802.bedmethyl
# test/compare_freq_bed_bed.sh -y test/expected/dRNA.mk.pileup.17802.bedmethyl test/tmp/dRNA.mm.freq.17802.bedmethyl test/tmp/dRNA_17802_compare || die "${testname} Comparison failed"
corr=`python3 test/compare.py test/expected/dRNA.mk.pileup.17802.bedmethyl test/tmp/dRNA.mm.freq.17802.bedmethyl`
echo "Correlation with modkit pileup for ChEBI: 17802: $corr"
if [ $(echo "$corr < 0.97" | bc -l) -eq 1 ]; then
    die "${testname} Correlation with modkit pileup decreased for ChEBI: 17802"
fi


testname="Test 17a: view ChEBI mod code test with pseudouridine (ChEBI: 17802)"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c "17802[*]" -t 8 test/tmp/genome_chr22.fa test/data/dRNA.bam > test/tmp/test17a.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test17a.tsv > test/tmp/test17a.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test17a.tsv > test/tmp/test17a.tsv.sorted
diff -q test/tmp/test17a.exp.tsv.sorted test/tmp/test17a.tsv.sorted || die "${testname} diff failed"

testname="Test 18: summary dRNA"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod summary -t 8 test/data/dRNA.bam > test/tmp/test18.tsv || die "${testname} Running the tool failed"
diff -q test/expected/test18.tsv test/tmp/test18.tsv || die "${testname} diff failed"

testname="Test 19: view RNA aligned to genome. Check if both positive and negative strands present"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c a[A] -t 8 test/tmp/genome_chr22.fa test/data/rna_algn_to_genome.bam > test/tmp/test19.tsv || die "${testname} Running the tool failed"
neg_count=`awk 'NR > 1 && $3 == "-" { print }' test/tmp/test19.tsv | wc -l`
echo "Negative strand alignments count: $neg_count"
pos_count=`awk 'NR > 1 && $3 == "+" { print }' test/tmp/test19.tsv | wc -l`
echo "Positive strand alignments count: $pos_count"
if [ "$neg_count" -ne 359 ] || [ "$pos_count" -ne 450 ]; then
    die "${testname} strand counts do not match expected values"
fi

#**** END of OLD TESTS - will be removed in future ****


#************** example-hifi C|m|? **************
testname="view m[CG] example-hifi.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[CG] -t 8 test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/example-hifi.mm.view.m.CG.tsv  || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/example-hifi.bam test/expected/example-hifi.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/example-hifi.mk.extract.CG.bed > test/expected/example-hifi.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-hifi.mk.extract.m.CG.bed test/tmp/example-hifi.mm.view.m.CG.tsv test/tmp/example_hifi_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_hifi_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[CG] example-hifi_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c h[CG] -t 8 test/tmp/genome_chr22.fa test/data/example-hifi_no_supps.bam > test/tmp/example-hifi.mm.view.h.CG.tsv  || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/example-hifi_no_supps.bam test/expected/example-hifi.mk.extract.CG.bed
# awk 'NR==1 || $14=="h"' test/expected/example-hifi.mk.extract.CG.bed > test/expected/example-hifi.mk.extract.h.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-hifi.mk.extract.h.CG.bed test/tmp/example-hifi.mm.view.h.CG.tsv test/tmp/example_hifi_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_hifi_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[*] example-hifi_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[*] -t 8 test/tmp/genome_chr22.fa test/data/example-hifi_no_supps.bam > test/tmp/example-hifi.mm.view.m.all.tsv  || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-hifi_no_supps.bam test/expected/example-hifi.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || ($14=="m" && (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17))))' test/expected/example-hifi.mk.extract.bed > test/expected/example-hifi.mk.extract.m.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-hifi.mk.extract.m.bed test/tmp/example-hifi.mm.view.m.all.tsv test/tmp/example_hifi_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_hifi_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[*] example-hifi_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c h[*] -t 8 test/tmp/genome_chr22.fa test/data/example-hifi_no_supps.bam > test/tmp/example-hifi.mm.view.h.all.tsv  || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1  --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-hifi_no_supps.bam test/expected/example-hifi.mk.extract.bed
awk 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || ($14=="h" && (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17))))' test/expected/example-hifi.mk.extract.bed > test/expected/example-hifi.mk.extract.h.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-hifi.mk.extract.h.bed test/tmp/example-hifi.mm.view.h.all.tsv test/tmp/example_hifi_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_hifi_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"




#************** example-ont_no_supps C|h|? C|m|? **************
testname="view m[CG] example-ont_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[CG] test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam > test/tmp/example-ont.mm.view.m.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam test/expected/example-ont_no_supps.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/example-ont_no_supps.mk.extract.CG.bed > test/expected/example-ont_no_supps.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont_no_supps.mk.extract.m.CG.bed test/tmp/example-ont_no_supps.mm.view.m.CG.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] example-ont_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[C] test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam > test/tmp/example-ont_no_supps.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam test/expected/example-ont_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/example-ont_no_supps.mk.extract.C.bed > test/expected/example-ont_no_supps.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont_no_supps.mk.extract.m.C.bed test/tmp/example-ont_no_supps.mm.view.m.C.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[CG] example-ont_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c h[CG] test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam > test/tmp/example-ont_no_supps.mm.view.h.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam test/expected/example-ont_no_supps.mk.extract.CG.bed
# awk 'NR==1 || $14=="h"' test/expected/example-ont_no_supps.mk.extract.CG.bed > test/expected/example-ont_no_supps.mk.extract.h.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont_no_supps.mk.extract.h.CG.bed test/tmp/example-ont_no_supps.mm.view.h.CG.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[C] example-ont_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c h[C] test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam > test/tmp/example-ont_no_supps.mm.view.h.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam test/expected/example-ont_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="h"' test/expected/example-ont_no_supps.mk.extract.C.bed > test/expected/example-ont_no_supps.mk.extract.h.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont_no_supps.mk.extract.h.C.bed test/tmp/example-ont_no_supps.mm.view.h.C.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[*] example-ont_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[*] test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam > test/tmp/example-ont_no_supps.mm.view.m.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam test/expected/example-ont_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || ($14=="m" && (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17))))' test/expected/example-ont_no_supps.mk.extract.bed > test/expected/example-ont_no_supps.mk.extract.m.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont_no_supps.mk.extract.m.bed test/tmp/example-ont_no_supps.mm.view.m.all.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[*] example-ont_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c h[*] test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam > test/tmp/example-ont_no_supps.mm.view.h.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1  --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-ont_no_supps.bam test/expected/example-ont_no_supps.mk.extract.bed
awk 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || ($14=="h" && (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17))))' test/expected/example-ont_no_supps.mk.extract.bed > test/expected/example-ont_no_supps.mk.extract.h.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont_no_supps.mk.extract.h.bed test/tmp/example-ont_no_supps.mm.view.h.all.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"


#************** dna_5mCG_5hmCG_mm_chr22_no_supps C|h|? C|m|? **************
testname="view m[CG] dna_5mCG_5hmCG_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[CG] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.m.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.CG.bed > test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.m.CG.bed test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.m.CG.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[CG] dna_5mCG_5hmCG_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c h[CG] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.h.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.CG.bed
# awk 'NR==1 || $14=="h"' test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.CG.bed > test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.h.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.h.CG.bed test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.h.CG.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] dna_5mCG_5hmCG_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[C] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.C.bed > test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.m.C.bed test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.m.C.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[C] dna_5mCG_5hmCG_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c h[C] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.h.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="h"' test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.C.bed > test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.h.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.h.C.bed test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.h.C.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[*] dna_5mCG_5hmCG_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[*] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.m.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || ($14=="m" && (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17))))'  test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.bed > test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.m.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.m.bed test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.m.all.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * dna_5mCG_5hmCG_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22_no_supps.bam test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))'  test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.bed > test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22_no_supps.mk.extract.matching.bed test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps.mm.view.all.all.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** dna_4mC_5mC_mm_chr22_no_supps C|21839|. C|m|. **************
testname="view m[CG] dna_4mC_5mC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[CG] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22_no_supps.bam > test/tmp/dna_4mC_5mC_mm_chr22_no_supps.mm.view.m.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22_no_supps.bam test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.CG.bed > test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.m.CG.bed test/tmp/dna_4mC_5mC_mm_chr22_no_supps.mm.view.m.CG.tsv test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] dna_4mC_5mC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[C] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22_no_supps.bam > test/tmp/dna_4mC_5mC_mm_chr22_no_supps.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22_no_supps.bam test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.C.bed > test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.m.C.bed test/tmp/dna_4mC_5mC_mm_chr22_no_supps.mm.view.m.C.tsv test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 21839[C] dna_4mC_5mC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c 21839[C] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22_no_supps.bam > test/tmp/dna_4mC_5mC_mm_chr22_no_supps.mm.view.21839.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22_no_supps.bam test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="21839"' test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.C.bed > test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.21839.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.21839.C.bed test/tmp/dna_4mC_5mC_mm_chr22_no_supps.mm.view.21839.C.tsv test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * dna_4mC_5mC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22_no_supps.bam > test/tmp/dna_4mC_5mC_mm_chr22_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22_no_supps.bam test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.bed > test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_4mC_5mC_mm_chr22_no_supps.mk.extract.matching.bed test/tmp/dna_4mC_5mC_mm_chr22_no_supps.mm.view.all.all.tsv test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** dna_5mCG_5hmCG_mm_chr22_no_supps C|h|. C|m|. **************
testname="view m[CG] dna_5mC_5hmC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[CG] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam > test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.m.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.CG.bed > test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.m.CG.bed test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.m.CG.tsv test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[CG] dna_5mC_5hmC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c h[CG] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam > test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.h.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.CG.bed
# awk 'NR==1 || $14=="h"' test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.CG.bed > test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.h.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.h.CG.bed test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.h.CG.tsv test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] dna_5mC_5hmC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[C] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam > test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.C.bed > test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.m.C.bed test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.m.C.tsv test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[C] dna_5mC_5hmC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c h[C] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam > test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.h.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="h"' test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.C.bed > test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.h.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.h.C.bed test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.h.C.tsv test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[*] dna_5mC_5hmC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c m[*] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam > test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.m.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || ($14=="m" && (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17))))' test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.bed > test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.m.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.m.bed test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.m.all.tsv test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * dna_5mC_5hmC_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam > test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22_no_supps.bam test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.bed > test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22_no_supps.mk.extract.matching.bed test/tmp/dna_5mC_5hmC_mm_chr22_no_supps.mm.view.all.all.tsv test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** dna_6mA_mm_chr22_no_supps A|a|. **************
testname="view a[A] dna_6mA_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c a[A] test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22_no_supps.bam > test/tmp/dna_6mA_mm_chr22_no_supps.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22_no_supps.bam test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.A.bed > test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.a.A.bed test/tmp/dna_6mA_mm_chr22_no_supps.mm.view.a.A.tsv test/tmp/dna_6mA_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view a[*] dna_6mA_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "a[*]" test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22_no_supps.bam > test/tmp/dna_6mA_mm_chr22_no_supps.mm.view.a.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22_no_supps.bam test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || ($14=="a" && (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17))))' test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.bed > test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.a.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.a.bed test/tmp/dna_6mA_mm_chr22_no_supps.mm.view.a.all.tsv test/tmp/dna_6mA_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * dna_6mA_mm_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22_no_supps.bam > test/tmp/dna_6mA_mm_chr22_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22_no_supps.bam test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.bed > test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_6mA_mm_chr22_no_supps.mk.extract.matching.bed test/tmp/dna_6mA_mm_chr22_no_supps.mm.view.all.all.tsv test/tmp/dna_6mA_mm_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"




# #************** rna_2OmeG_mm_hg38_chr22_no_supps G|19229|. **************
testname="view 19229[G] rna_2OmeG_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "19229[G]" test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22_no_supps.bam > test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps.mm.view.19229.G.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif G 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22_no_supps.bam test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.G.bed
# awk 'NR==1 || $14=="19229"' test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.G.bed > test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.19229.G.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.19229.G.bed test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps.mm.view.19229.G.tsv test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"


testname="view 19229[*] rna_2OmeG_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "19229[*]" test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22_no_supps.bam > test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps.mm.view.19229.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22_no_supps.bam test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || ($14=="19229" && (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17))))' test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.bed > test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.19229.all.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.19229.all.bed test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps.mm.view.19229.all.tsv test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_2OmeG_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22_no_supps.bam > test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22_no_supps.bam test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.bed > test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_hg38_chr22_no_supps.mk.extract.matching.bed test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps.mm.view.all.all.tsv test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps A|69426|. A|a|. A|17596|. **************
testname="view 69426[A] rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "69426[A]" test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mm.view.69426.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.A.bed
# awk 'NR==1 || $14=="69426"' test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.69426.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.69426.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mm.view.69426.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view a[A] rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "a[A]" test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.a.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mm.view.a.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 17596[A] rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "17596[A]" test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mm.view.17596.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.A.bed
# awk 'NR==1 || $14=="17596"' test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.17596.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.17596.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mm.view.17596.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.bam test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.bed > test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mk.extract.matching.bed test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps.mm.view.all.all.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"






#************** rna_m5C_2OmeC_mm_hg38_chr22_no_supps C|19228|. C|m|. **************
testname="view 19228[C] rna_m5C_2OmeC_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "19228[C]" test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.bam > test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mm.view.19228.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.bam test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="19228"' test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.C.bed > test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.19228.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.19228.C.bed test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mm.view.19228.C.tsv test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] rna_m5C_2OmeC_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "m[C]" test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.bam > test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.bam test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.m.C.bed
# awk 'NR==1 || $14=="m"' test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.m.C.bed > test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.19228.mC.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.19228.mC.bed test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mm.view.m.C.tsv test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_m5C_2OmeC_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.bam > test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.bam test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.bed > test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mk.extract.matching.bed test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps.mm.view.all.all.tsv test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"






#************** rna_m6A_DRACH_mm_hg38_chr22_no_supps A|a|? ***************
testname="view a[A] rna_m6A_DRACH_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "a[A]" test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22_no_supps.bam > test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22_no_supps.bam test/expected/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mk.extract.A.bed > test/expected/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mk.extract.a.A.bed test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mm.view.a.A.tsv test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_m6A_DRACH_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22_no_supps.bam > test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22_no_supps.bam test/expected/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mk.extract.bed > test/expected/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mk.extract.matching.bed test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps.mm.view.all.all.tsv test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_pseU_2OmeU_mm_hg38_chr22_no_supps T|19227|. T|17802|. ***************
testname="view 19227[T] rna_pseU_2OmeU_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "19227[T]" test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.bam > test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mm.view.19227.T.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif T 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.bam test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.T.bed
# awk 'NR==1 || $14=="19227"' test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.T.bed > test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.19227.T.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.19227.T.bed test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mm.view.19227.T.tsv test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 17802[T] rna_pseU_2OmeU_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "17802[T]" test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.bam > test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mm.view.17802.T.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif T 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.bam test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.T.bed
# awk 'NR==1 || $14=="17802"' test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.T.bed > test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.17802.T.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.17802.T.bed test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mm.view.17802.T.tsv test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_pseU_2OmeU_mm_hg38_chr22_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.bam > test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --kmer-size 1 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.bam test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.bed > test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mk.extract.matching.bed test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps.mm.view.all.all.tsv test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_2OmeG_mm_trans_ENST00000249299.7_no_supps G|19229|.**************
testname="view 19229[G] rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "19229[G]" test/data/transcript_ENST00000249299.7.fa test/data/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mm.view.19229.G.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --motif G 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mk.extract.G.bed
# awk 'NR==1 || $14=="19229"' test/expected/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mk.extract.G.bed > test/expected/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mk.extract.19229.G.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mk.extract.19229.G.bed test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mm.view.19229.G.tsv test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps_view_compare_mG || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps_view_compare_mG/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps_view_compare_mG/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps_view_compare_mG/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 19229[*] rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "19229[*]" test/data/transcript_ENST00000249299.7.fa test/data/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mm.view.19229.all.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --kmer-size 1 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mk.extract.bed > test/expected/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mk.extract.19229.all.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mk.extract.19229.all.bed test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps.mm.view.19229.all.tsv test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps A|69426|. A|a|. A|17596|. **************
testname="view 69426[A] rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "69426[A]" test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mm.view.69426.A.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --motif A 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.A.bed
# awk 'NR==1 || $14=="69426"' test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.69426.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.69426.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mm.view.69426.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view a[A] rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "a[A]" test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --motif A 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.a.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mm.view.a.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 17596[A] rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "17596[A]" test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mm.view.17596.A.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --motif A 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.A.bed
# awk 'NR==1 || $14=="17596"' test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.17596.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mk.extract.17596.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps.mm.view.17596.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps C|19228|. C|m|. **************
testname="view 19228[C] rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "19228[C]" test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mm.view.19228.C.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --motif C 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="19228"' test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.C.bed > test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.19228.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.19228.C.bed test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mm.view.19228.C.tsv test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "m[C]" test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --motif C 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.C.bed > test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.m.C.bed test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mm.view.m.C.tsv test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"


testname="view * rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --kmer-size 1 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.bed > test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mk.extract.matching.bed test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps.mm.view.all.all.tsv test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"


#************** rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps A|a|? ***************
testname="view a[A] rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "a[A]" test/data/transcript_ENST00000249299.7.fa test/data/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --motif A 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mk.extract.A.bed > test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mk.extract.a.A.bed test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mm.view.a.A.tsv test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/data/transcript_ENST00000249299.7.fa test/data/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --kmer-size 1 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mk.extract.bed > test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mk.extract.matching.bed test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps.mm.view.all.all.tsv test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"




#************** rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps T|19227|. T|17802|. ***************
testname="view 19227[T] rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "19227[T]" test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mm.view.19227.T.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --motif T 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.T.bed
# awk 'NR==1 || $14=="19227"' test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.T.bed > test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.19227.T.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.19227.T.bed test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mm.view.19227.T.tsv test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 17802[T] rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c "17802[T]" test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mm.view.17802.T.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --motif T 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.T.bed
# awk 'NR==1 || $14=="17802"' test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.T.bed > test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.17802.T.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.17802.T.bed test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mm.view.17802.T.tsv test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -t 8 -c '*' test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.bam > test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /storage/suneth/install/dist_modkit_v0.5.1_8fa79e3/modkit extract full -t 32 --kmer-size 1 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.bam test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.bed
# awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || (($21==16 && toupper($16)==c($17)) || ($21!=16 && toupper($16)==toupper($17)))' test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.bed > test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.matching.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mk.extract.matching.bed test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps.mm.view.all.all.tsv test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_no_supps_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"