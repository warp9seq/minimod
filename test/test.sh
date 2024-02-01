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
ex  ./minimod view test/alignment.bam | awk '{print $2"\t"$3}' | sort -k2,2n > test/alignment_actual.tsv
diff -q test/alignment_actual.tsv test/alignment_expected.tsv || die "diff failed: Invalid alignment positions"

echo "Tests passed"