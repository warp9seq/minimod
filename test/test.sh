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

#**** OLD TESTS ****

testname="Accuracy Test: freq results correlation with modkit and truthset"
exp_modkit_corr=0.97 # update this if the expected correlation changes
exp_truth_corr=0.85 # update this if the expected correlation changes
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -b test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/accu.freq.bedmethyl  || die "${testname} Running the tool failed"
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

testname="Test 1: view hifi"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[CG] test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test1.tsv  || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test1.tsv > test/tmp/test1.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test1.tsv > test/tmp/test1.tsv.sorted
diff -q test/tmp/test1.exp.tsv.sorted test/tmp/test1.tsv.sorted || die "${testname} diff failed"

testname="Test 2: view ont"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[CG] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2.tsv > test/tmp/test2.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2.tsv > test/tmp/test2.tsv.sorted
diff -q test/tmp/test2.exp.tsv.sorted test/tmp/test2.tsv.sorted || die "${testname} diff failed"

testname="Test 2a: view ont with insertions"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[CG] --insertions test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2a.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2a.tsv > test/tmp/test2a.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2a.tsv > test/tmp/test2a.tsv.sorted
diff -q test/tmp/test2a.exp.tsv.sorted test/tmp/test2a.tsv.sorted || die "${testname} diff failed"

testname="Test 2b: view ont with all contexts"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[*] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2b.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2b.tsv > test/tmp/test2b.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2b.tsv > test/tmp/test2b.tsv.sorted
diff -q test/tmp/test2b.exp.tsv.sorted test/tmp/test2b.tsv.sorted || die "${testname} diff failed"

testname="Test 2c: view ont all mod codes with wildcard"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c "*" test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test2c_wild.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2c_wild.tsv > test/tmp/test2c_wild.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2c_wild.tsv > test/tmp/test2c_wild.tsv.sorted
diff -q test/tmp/test2c_wild.exp.tsv.sorted test/tmp/test2c_wild.tsv.sorted || die "${testname} diff failed"

testname="Test 2c: view ont with haplotypes"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[CG] --haplotypes test/tmp/genome_chr1.fa test/data/hap.bam > test/tmp/test2c.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2c.tsv > test/tmp/test2c.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2c.tsv > test/tmp/test2c.tsv.sorted
diff -q test/tmp/test2c.exp.tsv.sorted test/tmp/test2c.tsv.sorted || die "${testname} diff failed"

testname="Test 2d: view ont with U in reference and t in required context"
echo -e "${BLUE}${testname}${NC}"
sed 's/T/U/g' test/tmp/genome_chr22.fa > test/tmp/genome_chr22_U.fa # replace all T with U in chr22 reference
ex  ./minimod view -c m[Ct] test/tmp/genome_chr22_U.fa test/data/example-ont.bam > test/tmp/test2d_U.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2d_U.tsv > test/tmp/test2d_U.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test2d_U.tsv > test/tmp/test2d_U.tsv.sorted
diff -q test/tmp/test2d_U.tsv.sorted test/tmp/test2d_U.tsv.sorted || die "${testname} diff failed"

testname="Test 3: freq hifi"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test3.tsv  || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test3.tsv > test/tmp/test3.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test3.tsv > test/tmp/test3.tsv.sorted
diff -q test/tmp/test3.exp.tsv.sorted test/tmp/test3.tsv.sorted || die "${testname} diff failed"

testname="Test 4: freq hifi bedmethyl output batch size 1"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -b -K 1 test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/test4.bedmethyl  || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k6,6 test/expected/test4.bedmethyl > test/tmp/test4.bedmethyl.exp.sorted
sort -k1,1 -k2,2n -k6,6 test/tmp/test4.bedmethyl > test/tmp/test4.bedmethyl.sorted
diff -q test/tmp/test4.bedmethyl.exp.sorted test/tmp/test4.bedmethyl.sorted || die "${testname} diff failed"

testname="Test 5: freq ont"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test5.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test5.tsv > test/tmp/test5.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test5.tsv > test/tmp/test5.tsv.sorted
diff -q test/tmp/test5.exp.tsv.sorted test/tmp/test5.tsv.sorted || die "${testname} diff failed"

testname="Test 5a: freq ont with insertions"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq --insertions test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test5a.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test5a.tsv > test/tmp/test5a.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test5a.tsv > test/tmp/test5a.tsv.sorted
diff -q test/tmp/test5a.exp.tsv.sorted test/tmp/test5a.tsv.sorted || die "${testname} diff failed"

testname="Test 5b: freq ont with all contexts"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c m[*] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test5b.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test5b.tsv > test/tmp/test5b.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test5b.tsv > test/tmp/test5b.tsv.sorted
diff -q test/tmp/test5b.exp.tsv.sorted test/tmp/test5b.tsv.sorted || die "${testname} diff failed"

testname="Test 5c: freq ont with haplotypes"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq --haplotypes test/tmp/genome_chr1.fa test/data/hap.bam > test/tmp/test5c.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test5c.tsv > test/tmp/test5c.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test5c.tsv > test/tmp/test5c.tsv.sorted
diff -q test/tmp/test5c.exp.tsv.sorted test/tmp/test5c.tsv.sorted || die "${testname} diff failed"

testname="Test 6: freq ont bedmethyl output"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -b test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test6.bedmethyl || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k6,6 test/expected/test6.bedmethyl > test/tmp/test6.bedmethyl.exp.sorted
sort -k1,1 -k2,2n -k6,6 test/tmp/test6.bedmethyl > test/tmp/test6.bedmethyl.sorted
diff -q test/tmp/test6.bedmethyl.exp.sorted test/tmp/test6.bedmethyl.sorted || die "${testname} diff failed"

testname="Test 7: freq ont with mod threshold"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -m 0.8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test7.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test7.tsv > test/tmp/test7.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test7.tsv > test/tmp/test7.tsv.sorted
diff -q test/tmp/test7.exp.tsv.sorted test/tmp/test7.tsv.sorted || die "${testname} diff failed"

testname="Test 8: freq ont with mod codes m and h"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c "m,h" -m 0.8,0.8 test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test8.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test8.tsv > test/tmp/test8.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test8.tsv > test/tmp/test8.tsv.sorted
diff -q test/tmp/test8.exp.tsv.sorted test/tmp/test8.tsv.sorted || die "${testname} diff failed"

testname="Test 9: freq ont with mod codes h only"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c "h" test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test9.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test9.tsv > test/tmp/test9.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test9.tsv > test/tmp/test9.tsv.sorted
diff -q test/tmp/test9.exp.tsv.sorted test/tmp/test9.tsv.sorted || die "${testname} diff failed"

testname="Test 10: view ont with mod threshold"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test10.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test10.tsv > test/tmp/test10.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test10.tsv > test/tmp/test10.tsv.sorted
diff -q test/tmp/test10.exp.tsv.sorted test/tmp/test10.tsv.sorted || die "${testname} diff failed"

testname="Test 11: view ont with mod codes m and h with different thresholds"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c "m,h"  test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test11.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test11.tsv > test/tmp/test11.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test11.tsv > test/tmp/test11.tsv.sorted
diff -q test/tmp/test11.exp.tsv.sorted test/tmp/test11.tsv.sorted || die "${testname} diff failed"

testname="Test 12: freq ont with mod codes m and h with different thresholds"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c "m,h" -m "0.8,0.5" test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/test12.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test12.tsv > test/tmp/test12.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test12.tsv > test/tmp/test12.tsv.sorted
diff -q test/tmp/test12.exp.tsv.sorted test/tmp/test12.tsv.sorted || die "${testname} diff failed"

testname="Test 13: view ont with -o flag"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view test/tmp/genome_chr22.fa test/data/example-ont.bam -o test/tmp/test13.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test2.tsv > test/tmp/test13.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test13.tsv > test/tmp/test13.tsv.sorted
diff -q test/tmp/test13.exp.tsv.sorted test/tmp/test13.tsv.sorted || die "${testname} diff failed"

testname="Test 14: freq ont with -o flag"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq test/tmp/genome_chr22.fa test/data/example-ont.bam -o test/tmp/test14.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/tmp/test14.tsv > test/tmp/test14.tsv.sorted
diff -q test/tmp/test5.exp.tsv.sorted test/tmp/test14.tsv.sorted || die "${testname} diff failed"

testname="Test 15: view ont with mod codes e and b"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c e,b test/tmp/genome_chr1.fa test/data/eb.bam > test/tmp/test15.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test15.tsv > test/tmp/test15.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test15.tsv > test/tmp/test15.tsv.sorted
diff -q test/tmp/test15.exp.tsv.sorted test/tmp/test15.tsv.sorted || die "${testname} diff failed"

testname="Test 16: freq ont with mod codes e and b"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c e,b -m 0.5 test/tmp/genome_chr1.fa test/data/eb.bam > test/tmp/test16.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k4,4 test/expected/test16.tsv > test/tmp/test16.exp.tsv.sorted
sort -k1,1 -k2,2n -k4,4 test/tmp/test16.tsv > test/tmp/test16.tsv.sorted
diff -q test/tmp/test16.exp.tsv.sorted test/tmp/test16.tsv.sorted || die "${testname} diff failed"

testname="Test 17: freq ChEBI mod code test with pseudouridine (ChEBI: 17802)"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -c "17802[*]" -b test/tmp/genome_chr22.fa test/data/dRNA.bam > test/tmp/dRNA.mm.freq.17802.bedmethyl || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --reference test/tmp/genome_chr22.fa test/data/dRNA.bam test/expected/dRNA.mk.pileup.bedmethyl || die "${testname} Running modkit pileup failed"
# awk 'NR==1 || $4=="17802"' test/expected/dRNA.mk.pileup.bedmethyl > test/expected/dRNA.mk.pileup.17802.bedmethyl
corr=`python3 test/compare.py test/expected/dRNA.mk.pileup.17802.bedmethyl test/tmp/dRNA.mm.freq.17802.bedmethyl`
echo "Correlation with modkit pileup for ChEBI: 17802: $corr"
if [ $(echo "$corr < 0.97" | bc -l) -eq 1 ]; then
    die "${testname} Correlation with modkit pileup decreased for ChEBI: 17802"
fi

testname="Test 17a: view ChEBI mod code test with pseudouridine (ChEBI: 17802)"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c "17802[*]" test/tmp/genome_chr22.fa test/data/dRNA.bam > test/tmp/test17a.tsv || die "${testname} Running the tool failed"
sort -k1,1 -k2,2n -k3,3 -k6,6 test/expected/test17a.tsv > test/tmp/test17a.exp.tsv.sorted
sort -k1,1 -k2,2n -k3,3 -k6,6 test/tmp/test17a.tsv > test/tmp/test17a.tsv.sorted
diff -q test/tmp/test17a.exp.tsv.sorted test/tmp/test17a.tsv.sorted || die "${testname} diff failed"

testname="Test 18: summary dRNA"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod summary test/data/dRNA.bam > test/tmp/test18.tsv || die "${testname} Running the tool failed"
diff -q test/expected/test18.tsv test/tmp/test18.tsv || die "${testname} diff failed"

testname="Test 19: view RNA aligned to genome. Check if both positive and negative strands present"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c a[A] test/tmp/genome_chr22.fa test/data/rna_algn_to_genome.bam > test/tmp/test19.tsv || die "${testname} Running the tool failed"
neg_count=`awk 'NR > 1 && $3 == "-" { print }' test/tmp/test19.tsv | wc -l`
echo "Negative strand alignments count: $neg_count"
pos_count=`awk 'NR > 1 && $3 == "+" { print }' test/tmp/test19.tsv | wc -l`
echo "Positive strand alignments count: $pos_count"
if [ "$neg_count" -ne 359 ] || [ "$pos_count" -ne 450 ]; then
    die "${testname} strand counts do not match expected values"
fi

#**** END of OLD TESTS ****





#************** example-hifi C|m|? **************
testname="view m[CG] example-hifi.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[CG] --skip-supplementary test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/example-hifi.mm.view.m.CG.tsv  || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/example-hifi.bam test/expected/example-hifi.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/example-hifi.mk.extract.CG.bed > test/expected/example-hifi.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-hifi.mk.extract.m.CG.bed test/tmp/example-hifi.mm.view.m.CG.tsv test/tmp/example_hifi_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_hifi_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[CG] example-hifi.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c h[CG] --skip-supplementary test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/example-hifi.mm.view.h.CG.tsv  || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/example-hifi.bam test/expected/example-hifi.mk.extract.CG.bed
# awk 'NR==1 || $14=="h"' test/expected/example-hifi.mk.extract.CG.bed > test/expected/example-hifi.mk.extract.h.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-hifi.mk.extract.h.CG.bed test/tmp/example-hifi.mm.view.h.CG.tsv test/tmp/example_hifi_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_hifi_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[*] example-hifi.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c m[*] --skip-supplementary test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/example-hifi.mm.view.m.all.tsv  || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-hifi.bam test/expected/example-hifi.mk.extract.bed
# awk 'NR==1 || $14=="m"' test/expected/example-hifi.mk.extract.bed > test/expected/example-hifi.mk.extract.m.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-hifi.mk.extract.m.bed test/tmp/example-hifi.mm.view.m.all.tsv test/tmp/example_hifi_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_hifi_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[*] example-hifi.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view -c h[*] --skip-supplementary test/tmp/genome_chr22.fa test/data/example-hifi.bam > test/tmp/example-hifi.mm.view.h.all.tsv  || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-hifi.bam test/expected/example-hifi.mk.extract.bed
# awk 'NR==1 || $14=="h"' test/expected/example-hifi.mk.extract.bed > test/expected/example-hifi.mk.extract.h.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-hifi.mk.extract.h.bed test/tmp/example-hifi.mm.view.h.all.tsv test/tmp/example_hifi_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_hifi_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_hifi_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"






#************** example-ont C|h|? C|m|? **************
testname="view m[CG] example-ont.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[CG] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/example-ont.mm.view.m.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/example-ont.bam test/expected/example-ont.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/example-ont.mk.extract.CG.bed > test/expected/example-ont.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont.mk.extract.m.CG.bed test/tmp/example-ont.mm.view.m.CG.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] example-ont.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[C] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/example-ont.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-ont.bam test/expected/example-ont.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/example-ont.mk.extract.C.bed > test/expected/example-ont.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont.mk.extract.m.C.bed test/tmp/example-ont.mm.view.m.C.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[CG] example-ont.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c h[CG] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/example-ont.mm.view.h.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/example-ont.bam test/expected/example-ont.mk.extract.CG.bed
# awk 'NR==1 || $14=="h"' test/expected/example-ont.mk.extract.CG.bed > test/expected/example-ont.mk.extract.h.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont.mk.extract.h.CG.bed test/tmp/example-ont.mm.view.h.CG.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[C] example-ont.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c h[C] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/example-ont.mm.view.h.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-ont.bam test/expected/example-ont.mk.extract.C.bed
# awk 'NR==1 || $14=="h"' test/expected/example-ont.mk.extract.C.bed > test/expected/example-ont.mk.extract.h.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont.mk.extract.h.C.bed test/tmp/example-ont.mm.view.h.C.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[*] example-ont.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[*] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/example-ont.mm.view.m.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-ont.bam test/expected/example-ont.mk.extract.bed
# awk 'NR==1 || $14=="m"' test/expected/example-ont.mk.extract.bed > test/expected/example-ont.mk.extract.m.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont.mk.extract.m.bed test/tmp/example-ont.mm.view.m.all.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[*] example-ont.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c h[*] test/tmp/genome_chr22.fa test/data/example-ont.bam > test/tmp/example-ont.mm.view.h.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full  --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/example-ont.bam test/expected/example-ont.mk.extract.bed
# awk 'NR==1 || $14=="h"' test/expected/example-ont.mk.extract.bed > test/expected/example-ont.mk.extract.h.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/example-ont.mk.extract.h.bed test/tmp/example-ont.mm.view.h.all.tsv test/tmp/example_ont_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/example_ont_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/example_ont_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"






#************** dna_5mCG_5hmCG_mm_chr22 C|h|? C|m|? **************
testname="view m[CG] dna_5mCG_5hmCG_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[CG] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.m.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.CG.bed > test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.m.CG.bed test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.m.CG.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[CG] dna_5mCG_5hmCG_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c h[CG] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.h.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.CG.bed
# awk 'NR==1 || $14=="h"' test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.CG.bed > test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.h.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.h.CG.bed test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.h.CG.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] dna_5mCG_5hmCG_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[C] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.C.bed > test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.m.C.bed test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.m.C.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[C] dna_5mCG_5hmCG_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c h[C] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.h.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.C.bed
# awk 'NR==1 || $14=="h"' test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.C.bed > test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.h.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.h.C.bed test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.h.C.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[*] dna_5mCG_5hmCG_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[*] test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.m.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.bed > test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.m.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.m.bed test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.m.all.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * dna_5mCG_5hmCG_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_chr22.mk.extract.bed test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.view.all.all.tsv test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[CG] --allow-secondary dna_5mCG_5hmCG_mm_with_secondary_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
# /data/suneth/install/samtools-1.23/samtools sort -n -o test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22_namesort.bam test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22.bam || die "${testname} Sorting the BAM file failed"
# /data/suneth/install/samtools-1.23/samtools fixmate -M test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22_namesort.bam test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN.bam || die "${testname} Running samtools fixmate failed"
ex  ./minimod view --allow-secondary test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN.bam > test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22.mm.view.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --allow-non-primary --cpg --force --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN.bam test/expected/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN.mk.extract.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN.mk.extract.bed > test/expected/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN.mk.extract.m.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN.mk.extract.m.bed test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN.mm.view.tsv test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_MN_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq m[CG] dna_5mCG_5hmCG_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -b --skip-supplementary test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.freq.m.CG.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --cpg --region chr22 --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/expected/dna_5mCG_5hmCG_mm_chr22.mk.pileup.CG.bed || die "${testname} Running modkit pileup failed"
# awk 'NR==1 || $4=="m"' test/expected/dna_5mCG_5hmCG_mm_chr22.mk.pileup.CG.bed > test/expected/dna_5mCG_5hmCG_mm_chr22.mk.pileup.m.CG.bed
corr=`test/compare.py test/expected/dna_5mCG_5hmCG_mm_chr22.mk.pileup.m.CG.bed test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.freq.m.CG.bed` || die "${testname} Comparison failed"
echo "Correlation of freq m[CG] with modkit pileup: $corr"
[ "$(echo "$corr < 0.999" | bc -l)" -eq 1 ] && die "${testname} Correlation of freq m[CG] with modkit pileup is less than 0.999"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq * dna_5mCG_5hmCG_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq -b --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.freq.all.alt.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --region chr22 --reference test/tmp/genome_chr22.fa test/data/dna_5mCG_5hmCG_mm_chr22.bam test/expected/dna_5mCG_5hmCG_mm_chr22.mk.pileup.bed || die "${testname} Running modkit pileup failed"
corr=`test/compare.py test/expected/dna_5mCG_5hmCG_mm_chr22.mk.pileup.bed test/tmp/dna_5mCG_5hmCG_mm_chr22.mm.freq.all.alt.bed` || die "${testname} Comparison failed"
echo "Correlation of freq * with modkit pileup: $corr"
[ "$(echo "$corr < 0.9708" | bc -l)" -eq 1 ] && die "${testname} Correlation of freq * with modkit pileup is less than 0.9708"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="summary dna_5mCG_5hmCG_mm_with_secondary_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod summary test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary.tsv || die "${testname} Running the tool failed"
ex  ./minimod summary --allow-secondary test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary_sec.tsv || die "${testname} Running the tool failed"
ex  ./minimod summary --skip-supplementary test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary_nosup.tsv || die "${testname} Running the tool failed"
ex  ./minimod summary --allow-secondary --skip-supplementary test/data/dna_5mCG_5hmCG_mm_with_secondary_chr22.bam > test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_summar_sec_nosup.tsv || die "${testname} Running the tool failed"
diff -q test/expected/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary.tsv test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary.tsv || die "${testname} Summary output does not match expected"
diff -q test/expected/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary_sec.tsv test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary_sec.tsv || die "${testname} Summary output with --allow-secondary does not match expected"
diff -q test/expected/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary_nosup.tsv test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary_nosup.tsv || die "${testname} Summary output with --skip-supplementary does not match expected"
diff -q test/expected/dna_5mCG_5hmCG_mm_with_secondary_chr22_summary_sec_nosup.tsv test/tmp/dna_5mCG_5hmCG_mm_with_secondary_chr22_summar_sec_nosup.tsv || die "${testname} Summary output with --allow-secondary and --skip-supplementary does not match expected"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** dna_4mC_5mC_mm_chr22 C|21839|. C|m|. **************
testname="view m[CG] dna_4mC_5mC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[CG] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam > test/tmp/dna_4mC_5mC_mm_chr22.mm.view.m.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/expected/dna_4mC_5mC_mm_chr22.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_4mC_5mC_mm_chr22.mk.extract.CG.bed > test/expected/dna_4mC_5mC_mm_chr22.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_4mC_5mC_mm_chr22.mk.extract.m.CG.bed test/tmp/dna_4mC_5mC_mm_chr22.mm.view.m.CG.tsv test/tmp/dna_4mC_5mC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] dna_4mC_5mC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[C] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam > test/tmp/dna_4mC_5mC_mm_chr22.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/expected/dna_4mC_5mC_mm_chr22.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_4mC_5mC_mm_chr22.mk.extract.C.bed > test/expected/dna_4mC_5mC_mm_chr22.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_4mC_5mC_mm_chr22.mk.extract.m.C.bed test/tmp/dna_4mC_5mC_mm_chr22.mm.view.m.C.tsv test/tmp/dna_4mC_5mC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 21839[C] dna_4mC_5mC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c 21839[C] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam > test/tmp/dna_4mC_5mC_mm_chr22.mm.view.21839.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/expected/dna_4mC_5mC_mm_chr22.mk.extract.C.bed
# awk 'NR==1 || $14=="21839"' test/expected/dna_4mC_5mC_mm_chr22.mk.extract.C.bed > test/expected/dna_4mC_5mC_mm_chr22.mk.extract.21839.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_4mC_5mC_mm_chr22.mk.extract.21839.C.bed test/tmp/dna_4mC_5mC_mm_chr22.mm.view.21839.C.tsv test/tmp/dna_4mC_5mC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * dna_4mC_5mC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam > test/tmp/dna_4mC_5mC_mm_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/expected/dna_4mC_5mC_mm_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_4mC_5mC_mm_chr22.mk.extract.bed test/tmp/dna_4mC_5mC_mm_chr22.mm.view.all.all.tsv test/tmp/dna_4mC_5mC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq m[C] dna_4mC_5mC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq --skip-supplementary -b -c m[C] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam > test/tmp/dna_4mC_5mC_mm_chr22.mm.freq.m.C.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --motif C 0 --region chr22 --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.C.bed || die "${testname} Running modkit pileup failed"
awk 'NR==1 || $4=="m"' test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.C.bed > test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.m.C.bed
corr=`test/compare.py test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.m.C.bed test/tmp/dna_4mC_5mC_mm_chr22.mm.freq.m.C.bed` || die "${testname} Comparison failed"
echo "Correlation of freq m[C] with modkit pileup: $corr"
[ "$(echo "$corr < 0.985" | bc -l)" -eq 1 ] && die "${testname} Correlation of freq m[C] with modkit pileup is less than 0.985"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq * dna_4mC_5mC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq --skip-supplementary -b -c '*' test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam > test/tmp/dna_4mC_5mC_mm_chr22.mm.freq.all.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --region chr22 --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.bed || die "${testname} Running modkit pileup failed"
corr=`test/compare.py test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.bed test/tmp/dna_4mC_5mC_mm_chr22.mm.freq.all.bed` || die "${testname} Comparison failed"
echo "Correlation of freq * with modkit pileup: $corr"
[ "$(echo "$corr < 0.998" | bc -l)" -eq 1 ] && die "${testname} Correlation of freq * with modkit pileup is less than 0.998"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq m[CG] dna_4mC_5mC_mm_chr22.bam compare with freq.sh output"
./minimod view --skip-supplementary -c m[CG] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam > test/tmp/dna_4mC_5mC_mm_chr22.mm.view.m.CG.tsv || die "Running minimod view failed"
test/freq.sh m 0.8 test/tmp/dna_4mC_5mC_mm_chr22.mm.view.m.CG.tsv > test/tmp/dna_4mC_5mC_mm_chr22.freqscript.view.m.CG.tsv || die "Running freq.sh on minimod view output failed"
ex  ./minimod freq --skip-supplementary -b -c m[CG] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam > test/tmp/dna_4mC_5mC_mm_chr22.mm.freq.m.CG.bed || die "${testname} Running the tool failed"
test/compare_freq_mmbed_scripttsv.sh -y test/tmp/dna_4mC_5mC_mm_chr22.mm.freq.m.CG.bed test/tmp/dna_4mC_5mC_mm_chr22.freqscript.view.m.CG.tsv test/tmp/dna_4mC_5mC_mm_chr22_freq_compare_script || die "${testname} Comparison of minimod freq output with freq.sh output failed"
echo "debug logs"
wc -l test/tmp/dna_4mC_5mC_mm_chr22_freq_compare_script/missing_in_file1.tsv
wc -l test/tmp/dna_4mC_5mC_mm_chr22_freq_compare_script/missing_in_file2.tsv
wc -l test/tmp/dna_4mC_5mC_mm_chr22_freq_compare_script/large_freq_diff.tsv
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_freq_compare_script/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod freq missing records compared to freq.sh output"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_freq_compare_script/missing_in_file2.tsv)" -gt 1 ] && die "${testname} freq.sh output missing records compared to minimod freq"
[ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_freq_compare_script/large_freq_diff.tsv)" -gt 1 ] && die "${testname} Records with large freq diff between minimod freq and freq.sh output"
echo -e "${GREEN}${testname} passed!${NC}\n"


# THIS IS TEST IS COMMENTED OUT because minimod can't match modkit's 3 way classification oputput
# testname="freq m[CG] dna_4mC_5mC_mm_chr22.bam using compare_freq_bed_bed.sh"
# echo -e "${BLUE}${testname}${NC}"
# ex  ./minimod freq --skip-supplementary -m 0.8828125 -b -c m[CG] test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam > test/tmp/dna_4mC_5mC_mm_chr22.mm.freq.m.CG.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --cpg --region chr22 --reference test/tmp/genome_chr22.fa test/data/dna_4mC_5mC_mm_chr22.bam test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.CG.bed || die "${testname} Running modkit pileup failed"
# awk 'NR==1 || $4=="m"' test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.CG.bed > test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.m.CG.bed
# test/compare_freq_bed_bed.sh -y test/expected/dna_4mC_5mC_mm_chr22.mk.pileup.m.CG.bed test/tmp/dna_4mC_5mC_mm_chr22.mm.freq.m.CG.bed test/tmp/dna_4mC_5mC_mm_chr22_freq_compare || die "${testname} Comparison failed"
# [ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_freq_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod freq missing records compared to modkit pileup"
# [ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_freq_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit pileup missing records compared to minimod freq"
# [ "$(wc -l < test/tmp/dna_4mC_5mC_mm_chr22_freq_compare/large_freq_diff.tsv)" -gt 1 ] && die "${testname} Records with large freq diff between minimod freq and modkit pileup"
# echo -e "${GREEN}${testname} passed!${NC}\n"






#************** dna_5mC_5hmC_mm_chr22 C|h|. C|m|. **************
testname="view m[CG] dna_5mC_5hmC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[CG] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam > test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.m.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.CG.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.CG.bed > test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.m.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.m.CG.bed test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.m.CG.tsv test/tmp/dna_5mC_5hmC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[CG] dna_5mC_5hmC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c h[CG] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam > test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.h.CG.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --cpg --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.CG.bed
# awk 'NR==1 || $14=="h"' test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.CG.bed > test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.h.CG.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.h.CG.bed test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.h.CG.tsv test/tmp/dna_5mC_5hmC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] dna_5mC_5hmC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[C] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam > test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.C.bed > test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.m.C.bed test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.m.C.tsv test/tmp/dna_5mC_5hmC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view h[C] dna_5mC_5hmC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c h[C] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam > test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.h.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.C.bed
# awk 'NR==1 || $14=="h"' test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.C.bed > test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.h.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.h.C.bed test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.h.C.tsv test/tmp/dna_5mC_5hmC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[*] dna_5mC_5hmC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c m[*] test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam > test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.m.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.bed
# awk 'NR==1 || $14=="m"' test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.bed > test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.m.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.m.bed test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.m.all.tsv test/tmp/dna_5mC_5hmC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * dna_5mC_5hmC_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam > test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_5mC_5hmC_mm_chr22.bam test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_5mC_5hmC_mm_chr22.mk.extract.bed test/tmp/dna_5mC_5hmC_mm_chr22.mm.view.all.all.tsv test/tmp/dna_5mC_5hmC_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_5mC_5hmC_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** dna_6mA_mm_chr22 A|a|. **************
testname="view a[A] dna_6mA_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c a[A] test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam > test/tmp/dna_6mA_mm_chr22.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam test/expected/dna_6mA_mm_chr22.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/dna_6mA_mm_chr22.mk.extract.A.bed > test/expected/dna_6mA_mm_chr22.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_6mA_mm_chr22.mk.extract.a.A.bed test/tmp/dna_6mA_mm_chr22.mm.view.a.A.tsv test/tmp/dna_6mA_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view a[*] dna_6mA_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "a[*]" test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam > test/tmp/dna_6mA_mm_chr22.mm.view.a.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam test/expected/dna_6mA_mm_chr22.mk.extract.bed
# awk 'NR==1 || $14=="a"' test/expected/dna_6mA_mm_chr22.mk.extract.bed > test/expected/dna_6mA_mm_chr22.mk.extract.a.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_6mA_mm_chr22.mk.extract.a.bed test/tmp/dna_6mA_mm_chr22.mm.view.a.all.tsv test/tmp/dna_6mA_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * dna_6mA_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam > test/tmp/dna_6mA_mm_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam test/expected/dna_6mA_mm_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/dna_6mA_mm_chr22.mk.extract.bed test/tmp/dna_6mA_mm_chr22.mm.view.all.all.tsv test/tmp/dna_6mA_mm_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/dna_6mA_mm_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq a[A] dna_6mA_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq --skip-supplementary -b -c a[A] test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam > test/tmp/dna_6mA_mm_chr22.mm.freq.a.A.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --motif A 0 --region chr22 --reference test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam test/expected/dna_6mA_mm_chr22.mk.pileup.A.bed || die "${testname} Running modkit pileup failed"
# awk 'NR==1 || $4=="a"' test/expected/dna_6mA_mm_chr22.mk.pileup.A.bed > test/expected/dna_6mA_mm_chr22.mk.pileup.a.A.bed
corr=`test/compare.py test/expected/dna_6mA_mm_chr22.mk.pileup.a.A.bed test/tmp/dna_6mA_mm_chr22.mm.freq.a.A.bed` || die "${testname} Comparison failed"
echo "Correlation of freq a[A] with modkit pileup: $corr"
[ "$(echo "$corr < 0.988" | bc -l)" -eq 1 ] && die "${testname} Correlation of freq a[A] with modkit pileup is less than 0.988"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq * dna_6mA_mm_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq --skip-supplementary -b -c '*' test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam > test/tmp/dna_6mA_mm_chr22.mm.freq.all.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --region chr22 --reference test/tmp/genome_chr22.fa test/data/dna_6mA_mm_chr22.bam test/expected/dna_6mA_mm_chr22.mk.pileup.bed || die "${testname} Running modkit pileup failed"
corr=`test/compare.py test/expected/dna_6mA_mm_chr22.mk.pileup.bed test/tmp/dna_6mA_mm_chr22.mm.freq.all.bed` || die "${testname} Comparison failed" 
echo "Correlation of freq * with modkit pileup: $corr"
[ "$(echo "$corr < 0.989" | bc -l)" -eq 1 ] && die "${testname} Correlation of freq * with modkit pileup is less than 0.989"
echo -e "${GREEN}${testname} passed!${NC}\n"





# #************** rna_2OmeG_mm_hg38_chr22 G|19229|. **************
testname="view 19229[G] rna_2OmeG_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "19229[G]" test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22.bam > test/tmp/rna_2OmeG_mm_hg38_chr22.mm.view.19229.G.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif G 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22.bam test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.G.bed
# awk 'NR==1 || $14=="19229"' test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.G.bed > test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.19229.G.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.19229.G.bed test/tmp/rna_2OmeG_mm_hg38_chr22.mm.view.19229.G.tsv test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 19229[*] rna_2OmeG_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "19229[*]" test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22.bam > test/tmp/rna_2OmeG_mm_hg38_chr22.mm.view.19229.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22.bam test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.bed
# awk 'NR==1 || $14=="19229"' test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.bed > test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.19229.all.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.19229.all.bed test/tmp/rna_2OmeG_mm_hg38_chr22.mm.view.19229.all.tsv test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_2OmeG_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22.bam > test/tmp/rna_2OmeG_mm_hg38_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_2OmeG_mm_hg38_chr22.bam test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_hg38_chr22.mk.extract.bed test/tmp/rna_2OmeG_mm_hg38_chr22.mm.view.all.all.tsv test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_inosine_m6A_2OmeA_mm_hg38_chr22 A|69426|. A|a|. A|17596|. **************
testname="view 69426[A] rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "69426[A]" test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mm.view.69426.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.A.bed
# awk 'NR==1 || $14=="69426"' test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.69426.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.69426.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mm.view.69426.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view a[A] rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "a[A]" test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.a.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mm.view.a.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 17596[A] rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "17596[A]" test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mm.view.17596.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.A.bed
# awk 'NR==1 || $14=="17596"' test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.17596.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.17596.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mm.view.17596.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_inosine_m6A_2OmeA_mm_hg38_chr22.bam test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mk.extract.bed test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22.mm.view.all.all.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"






#************** rna_m5C_2OmeC_mm_hg38_chr22 C|19228|. C|m|. **************
testname="view 19228[C] rna_m5C_2OmeC_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "19228[C]" test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22.bam > test/tmp/rna_m5C_2OmeC_mm_hg38_chr22.mm.view.19228.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22.bam test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.C.bed
# awk 'NR==1 || $14=="19228"' test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.C.bed > test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.19228.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.19228.C.bed test/tmp/rna_m5C_2OmeC_mm_hg38_chr22.mm.view.19228.C.tsv test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] rna_m5C_2OmeC_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "m[C]" test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22.bam > test/tmp/rna_m5C_2OmeC_mm_hg38_chr22.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif C 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22.bam test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.m.C.bed
# awk 'NR==1 || $14=="m"' test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.m.C.bed > test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.19228.mC.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.19228.mC.bed test/tmp/rna_m5C_2OmeC_mm_hg38_chr22.mm.view.m.C.tsv test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_m5C_2OmeC_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22.bam > test/tmp/rna_m5C_2OmeC_mm_hg38_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_mm_hg38_chr22.bam test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_hg38_chr22.mk.extract.bed test/tmp/rna_m5C_2OmeC_mm_hg38_chr22.mm.view.all.all.tsv test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"






#************** rna_m6A_DRACH_mm_hg38_chr22 A|a|? ***************
testname="view a[A] rna_m6A_DRACH_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "a[A]" test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam > test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.extract.A.bed > test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.extract.a.A.bed test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.view.a.A.tsv test/tmp/rna_m6A_DRACH_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_m6A_DRACH_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam > test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.extract.bed test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.view.all.all.tsv test/tmp/rna_m6A_DRACH_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq a[A] rna_m6A_DRACH_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq --skip-supplementary -b -c "a[A]" test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam > test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.a.A.freq.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --motif A 0 --region chr22 --reference test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup.A.bed || die "${testname} Running modkit pileup failed"
# awk 'NR==1 || $4=="a"' test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup.A.bed > test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup.a.A.bed
corr=`test/compare.py test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup.a.A.bed test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.a.A.freq.bed` || die "${testname} Correlation comparison failed"
[ "$(echo "$corr < 0.995" | bc -l)" -eq 1 ] && die "${testname} Correlation between minimod freq and modkit pileup is less than 0.995: $corr"
echo "Correlation of freq a[A] with modkit pileup: $corr" 
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq * rna_m6A_DRACH_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq --skip-supplementary -b -c '*' test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam > test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.all.all.freq.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --region chr22 --reference test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup.bed || die "${testname} Running modkit pileup failed"
corr=`test/compare.py test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup.bed test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.all.all.freq.bed` || die "${testname} Correlation comparison failed"
[ "$(echo "$corr < 0.995" | bc -l)" -eq 1 ] && die "${testname} Correlation between minimod freq and modkit pileup is less than 0.995: $corr"
echo "Correlation of freq * with modkit pileup: $corr"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq a[A] rna_m6A_DRACH_mm_hg38_chr22.bam using compare_freq_bed_bed.sh"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod freq --skip-supplementary -m 0.9 -b -c "a[A]" test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam > test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.a.A.freq_0.96.bed || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit pileup --filter-threshold A:0.9 --motif A 0 --region chr22 --reference test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup_0.96.A.bed || die "${testname} Running modkit pileup failed"
# awk 'NR==1 || $4=="a"' test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup_0.96.A.bed > test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup_0.96.a.A.bed
test/compare_freq_bed_bed.sh -y test/expected/rna_m6A_DRACH_mm_hg38_chr22.mk.pileup_0.96.a.A.bed test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.a.A.freq_0.96.bed test/tmp/rna_m6A_DRACH_mm_hg38_chr22_freq_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_freq_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod freq missing records compared to modkit pileup"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_freq_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit pileup missing records compared to minimod freq"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_freq_compare/large_freq_diff.tsv)" -gt 1 ] && die "${testname} Records with large freq diff between minimod freq and modkit pileup"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="freq a[A] rna_m6A_DRACH_mm_hg38_chr22.bam compare with freq.sh output"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "a[A]" test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam > test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.view.a.A.tsv || die "${testname} Running minimod view failed"
test/freq.sh a 0.8 test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.view.a.A.tsv > test/tmp/rna_m6A_DRACH_mm_hg38_chr22.freqscript.a.A.tsv || die "${testname} Running freq.sh failed"
ex  ./minimod freq --skip-supplementary -b -c "a[A]" test/tmp/genome_chr22.fa test/data/rna_m6A_DRACH_mm_hg38_chr22.bam > test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.a.A.freq.bed || die "${testname} Running the tool failed"
test/compare_freq_mmbed_scripttsv.sh -y test/tmp/rna_m6A_DRACH_mm_hg38_chr22.mm.a.A.freq.bed test/tmp/rna_m6A_DRACH_mm_hg38_chr22.freqscript.a.A.tsv test/tmp/rna_m6A_DRACH_mm_hg38_chr22_freq_compare_script || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_freq_compare_script/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod freq missing records compared to freq.sh output"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_freq_compare_script/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} freq.sh output missing records compared to minimod freq"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_hg38_chr22_freq_compare_script/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large freq diff between minimod freq and freq.sh output"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_pseU_2OmeU_mm_hg38_chr22 T|19227|. T|17802|. ***************
testname="view 19227[T] rna_pseU_2OmeU_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "19227[T]" test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22.bam > test/tmp/rna_pseU_2OmeU_mm_hg38_chr22.mm.view.19227.T.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif T 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22.bam test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.T.bed
# awk 'NR==1 || $14=="19227"' test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.T.bed > test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.19227.T.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.19227.T.bed test/tmp/rna_pseU_2OmeU_mm_hg38_chr22.mm.view.19227.T.tsv test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 17802[T] rna_pseU_2OmeU_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "17802[T]" test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22.bam > test/tmp/rna_pseU_2OmeU_mm_hg38_chr22.mm.view.17802.T.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif T 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22.bam test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.T.bed
# awk 'NR==1 || $14=="17802"' test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.T.bed > test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.17802.T.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.17802.T.bed test/tmp/rna_pseU_2OmeU_mm_hg38_chr22.mm.view.17802.T.tsv test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_pseU_2OmeU_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22.bam > test/tmp/rna_pseU_2OmeU_mm_hg38_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_pseU_2OmeU_mm_hg38_chr22.bam test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_hg38_chr22.mk.extract.bed test/tmp/rna_pseU_2OmeU_mm_hg38_chr22.mm.view.all.all.tsv test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22 A|17596|. C|m|. A|a|. C|19228|. T|19227|. A|69426|. G|19229|. T|17802|. ***************
testname="view 17596[A] rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "17596[A]" test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.bam > test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mm.view.17596.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --motif A 0 --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.bam test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mk.extract.A.bed
# awk 'NR==1 || $14=="17596"' test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mk.extract.A.bed > test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mk.extract.17596.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mk.extract.17596.A.bed test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mm.view.17596.A.tsv test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view  --skip-supplementary -c '*' test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.bam > test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full --mapped-only --force --reference test/tmp/genome_chr22.fa test/data/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.bam test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mk.extract.bed test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22.mm.view.all.all.tsv test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22_view_compare/missing_in_file2.tsv)" -gt 1 ]&& die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_hg38_chr22_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_2OmeG_mm_trans_ENST00000249299.7 G|19229|.**************
testname="view 19229[G] rna_2OmeG_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "19229[G]" test/data/transcript_ENST00000249299.7.fa test/data/rna_2OmeG_mm_trans_ENST00000249299.7.bam > test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7.mm.view.19229.G.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif G 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_2OmeG_mm_trans_ENST00000249299.7.bam test/expected/rna_2OmeG_mm_trans_ENST00000249299.7.mk.extract.G.bed
# awk 'NR==1 || $14=="19229"' test/expected/rna_2OmeG_mm_trans_ENST00000249299.7.mk.extract.G.bed > test/expected/rna_2OmeG_mm_trans_ENST00000249299.7.mk.extract.19229.G.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_trans_ENST00000249299.7.mk.extract.19229.G.bed test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7.mm.view.19229.G.tsv test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_view_compare_mG || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_view_compare_mG/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_view_compare_mG/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_view_compare_mG/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 19229[*] rna_2OmeG_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "19229[*]" test/data/transcript_ENST00000249299.7.fa test/data/rna_2OmeG_mm_trans_ENST00000249299.7.bam > test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7.mm.view.19229.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_2OmeG_mm_trans_ENST00000249299.7.bam test/expected/rna_2OmeG_mm_trans_ENST00000249299.7.mk.extract.bed
# awk 'NR==1 || $14=="19229"' test/expected/rna_2OmeG_mm_trans_ENST00000249299.7.mk.extract.bed > test/expected/rna_2OmeG_mm_trans_ENST00000249299.7.mk.extract.19229.all.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_2OmeG_mm_trans_ENST00000249299.7.mk.extract.19229.all.bed test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7.mm.view.19229.all.tsv test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_2OmeG_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7 A|69426|. A|a|. A|17596|. **************
testname="view 69426[A] rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "69426[A]" test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mm.view.69426.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif A 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.bam test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.A.bed
# awk 'NR==1 || $14=="69426"' test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.69426.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.69426.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mm.view.69426.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view a[A] rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "a[A]" test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif A 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.bam test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.a.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mm.view.a.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 17596[A] rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "17596[A]" test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.bam > test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mm.view.17596.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif A 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.bam test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.A.bed
# awk 'NR==1 || $14=="17596"' test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.A.bed > test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.17596.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mk.extract.17596.A.bed test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7.mm.view.17596.A.tsv test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_inosine_m6A_2OmeA_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_m5C_2OmeC_mm_trans_ENST00000249299.7 C|19228|. C|m|. **************
testname="view 19228[C] rna_m5C_2OmeC_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "19228[C]" test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.bam > test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mm.view.19228.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif C 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.bam test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.C.bed
# awk 'NR==1 || $14=="19228"' test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.C.bed > test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.19228.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.19228.C.bed test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mm.view.19228.C.tsv test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view m[C] rna_m5C_2OmeC_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "m[C]" test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.bam > test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mm.view.m.C.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif C 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.bam test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.C.bed
# awk 'NR==1 || $14=="m"' test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.C.bed > test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.m.C.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.m.C.bed test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mm.view.m.C.tsv test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_m5C_2OmeC_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.bam > test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.bam test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mk.extract.bed test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7.mm.view.all.all.tsv test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_m6A_DRACH_mm_trans_ENST00000249299.7 A|a|? ***************
testname="view a[A] rna_m6A_DRACH_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "a[A]" test/data/transcript_ENST00000249299.7.fa test/data/rna_m6A_DRACH_mm_trans_ENST00000249299.7.bam > test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mm.view.a.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif A 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m6A_DRACH_mm_trans_ENST00000249299.7.bam test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mk.extract.A.bed
# awk 'NR==1 || $14=="a"' test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mk.extract.A.bed > test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mk.extract.a.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mk.extract.a.A.bed test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mm.view.a.A.tsv test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_m6A_DRACH_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/data/transcript_ENST00000249299.7.fa test/data/rna_m6A_DRACH_mm_trans_ENST00000249299.7.bam > test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m6A_DRACH_mm_trans_ENST00000249299.7.bam test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mk.extract.bed
# awk 'NR==1 || $14=="a"' test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mk.extract.bed > test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mk.extract.a.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mk.extract.a.bed test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7.mm.view.all.all.tsv test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m6A_DRACH_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_pseU_2OmeU_mm_trans_ENST00000249299.7 T|19227|. T|17802|. ***************
testname="view 19227[T] rna_pseU_2OmeU_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "19227[T]" test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.bam > test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mm.view.19227.T.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif T 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.bam test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.T.bed
# awk 'NR==1 || $14=="19227"' test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.T.bed > test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.19227.T.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.19227.T.bed test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mm.view.19227.T.tsv test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view 17802[T] rna_pseU_2OmeU_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c "17802[T]" test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.bam > test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mm.view.17802.T.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif T 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.bam test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.T.bed
# awk 'NR==1 || $14=="17802"' test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.T.bed > test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.17802.T.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.17802.T.bed test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mm.view.17802.T.tsv test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_pseU_2OmeU_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view --skip-supplementary -c '*' test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.bam > test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.bam test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mk.extract.bed test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7.mm.view.all.all.tsv test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_pseU_2OmeU_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"





#************** rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7 A|17596|. C|m|. A|a|. C|19228|. T|19227|. A|69426|. G|19229|. T|17802|. ***************
testname="view 69426[A] rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view  --skip-supplementary -c "69426[A]" test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.bam > test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mm.view.69426.A.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --motif A 0 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.bam test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mk.extract.A.bed
# awk 'NR==1 || $14=="69426"' test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mk.extract.A.bed > test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mk.extract.69426.A.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mk.extract.69426.A.bed test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mm.view.69426.A.tsv test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"

testname="view * rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.bam"
echo -e "${BLUE}${testname}${NC}"
ex  ./minimod view  --skip-supplementary -c '*' test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.bam > test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mm.view.all.all.tsv || die "${testname} Running the tool failed"
# /install/modkit-v0.5.1/modkit extract full -t 32 --mapped-only --force --reference test/data/transcript_ENST00000249299.7.fa test/data/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.bam test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mk.extract.bed
test/compare_view_mkbed_mmtsv.sh -y test/expected/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mk.extract.bed test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7.mm.view.all.all.tsv test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7_view_compare || die "${testname} Comparison failed"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7_view_compare/missing_in_file1.tsv)" -gt 1 ] && die "${testname} minimod view missing records compared to modkit extract full"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7_view_compare/missing_in_file2.tsv)" -gt 1 ] && die "${testname} modkit extract full missing records compared to minimod view"
[ "$(wc -l < test/tmp/rna_m5C_2OmeC_inosine_m6A_2OmeA_pseU_2OmeU_2OmeG_mm_trans_ENST00000249299.7_view_compare/large_prob_diff.tsv)" -gt 1 ] && die "${testname} Records with large prob diff between minimod view and modkit extract full"
echo -e "${GREEN}${testname} passed!${NC}\n"
