#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
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

echo "Test 1"
ex  ./minimod view test/r84088_20230609_025659_1_A01_fmr1.bam > test/tmp.txt  || die "Running the tool failed"
# diff -q test/example.exp test/tmp.txt || die "diff failed"

echo "Test 2"
ex  ./minimod view test/alignment.bam | awk 'NR>1{print $2"\t"$3}' | sort -n -k 1 > test/alignment_actual.tsv || die "Running the tool failed"
diff -q test/alignment_actual.tsv test/alignment_expected.tsv || die "diff failed: Invalid alignment positions"

echo "Test 3"
ex  ./minimod view test/example-ont.bam | awk 'NR>1{print $2"\t"$11}' | sort -n -k 1 > test/example-ont_prob_actual.tsv || die "Running the tool failed"
diff -q test/example-ont_prob_actual.tsv test/example-ont_prob_expected.tsv || die "diff failed: Invalid probabilities"

echo "Tests passed"