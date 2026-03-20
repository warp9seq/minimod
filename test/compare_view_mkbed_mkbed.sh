#!/usr/bin/env bash
#
# Usage: ./compare_view_mkbed_mkbed.sh file1.bedmethyl file2.tsv output_dir
#
# This script compares two files based on a constructed key and outputs
# the matched records, records with large differences, and records
# missing in either file.

usage() {
    echo "Usage: $0 file1.bedmethyl file2.bedmethyl output_dir"
    echo "Options:"
    echo "  -y : Overwrite existing output files without prompt."
}

# --- Parse Options ---
while getopts ":y" opt; do
    case $opt in
        y)
            OVERWRITE=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

# --- Configuration ---
TOLERANCE=0.002 # Define the probability difference tolerance
# --- Argument Check ---
if [[ $# -ne 3 ]]; then
    usage
    exit 1
fi
f1="$1"
f2="$2"
out_dir=$3
# Check if extension of files are correct. file1 should be .bed or .bedmethyl and file2 should be .bed or .bedmethyl
if [[ "$f1" != *.bed  && "$f1" != *.bedmethyl ]] || [[ "$f2" != *.bed  && "$f2" != *.bedmethyl ]]; then
    echo "Error: file1 should have .bed or .bedmethyl extension and file2 should have .bed or .bedmethyl extension."
    exit 1
fi
# Extract base names without extension
base1=$(basename "$f1" .bedmethyl)
base2=$(basename "$f2" .bedmethyl)

if [[ "$f1" == *.bed ]]; then
    base1=$(basename "$f1" .bed)
fi

mkdir -p "$out_dir"
script_name=$(basename "$0" .sh)
# --- Output Files ---
OUT_MISSING_1="$out_dir/missing_in_file1.tsv" # Keys in file2 but not in file1
OUT_MISSING_2="$out_dir/missing_in_file2.tsv" # Keys in file1 but not in file2
OUT_MATCH="$out_dir/in_both.tsv"              # Keys in both with small prob difference
OUT_LARGE_DIFF="$out_dir/large_prob_diff.tsv" # Keys in both with large prob difference
# Remove output files if they already exist to avoid appending to old data. Ask for confirmation.
if [[ -z "$OVERWRITE" && ( -f "$OUT_MISSING_1" || -f "$OUT_MISSING_2" || -f "$OUT_MATCH" || -f "$OUT_LARGE_DIFF" ) ]]; then
    read -p "Output files already exist. Overwrite? (y/n) " choice
    if [[ "$choice" != "y" ]]; then
        echo "Exiting without overwriting files."
        exit 1
    fi
fi
rm -f "$OUT_MISSING_1" "$OUT_MISSING_2" "$OUT_MATCH" "$OUT_LARGE_DIFF"
# --- AWK Processing ---
awk -F'\t' -v tol="$TOLERANCE" \
    -v out_match="$OUT_MATCH" \
    -v out_large_diff="$OUT_LARGE_DIFF" \
    -v out_missing1="$OUT_MISSING_1" \
    -v out_missing2="$OUT_MISSING_2" '
    # Helper function to trim whitespace and carriage returns
    function clean(x) {
        gsub(/\r/, "", x)
        gsub(/^ +| +$/, "", x)
        return x
    }
    # Set headers for output files at the beginning
    BEGIN {
        print "key\tprob\tdiff" > out_match
        print "key\tprob_file1\tprob_file2\tdiff" > out_large_diff
        print "key\tprob_file2" > out_missing1
        print "key\tprob_file1" > out_missing2
    }
    # First pass: read file1 into an associative array `prob`
    # NR==FNR is true only while reading the first file
    NR==FNR {
        # skip header lines
        if (FNR == 1) next
        for (i=1; i<=NF; i++) $i = clean($i)
        # Construct the key from the first 6 columns
        key = $4 "\t" $3 "\t" $6 "\t" $1 "\t" $2 "\t" $14
        prob[key] = $13
        next
    }

    # Second pass: process file2 and compare against the `prob` array
    {
        # skip header lines
        if (FNR == 1) next
        for (i=1; i<=NF; i++) $i = clean($i)
        # Construct the key from file2 columns in the order that matches file1
        key1 = $4 "\t" $3 "\t" $6 "\t" $1 "\t" $2 "\t" $14
        if (key1 in prob) {
            # Key found in both files: calculate the difference
            diff = prob[key1] - $13
            if (diff < 0) diff = -diff # Absolute difference
            if (diff <= tol) {
                # Difference is within tolerance: write to matched.tsv
                print key1 "\t" prob[key1] "\t" diff > out_match
            } else {
                # Difference is too large: write to large_prob_diff.tsv
                print key1 "\t" prob[key1] "\t" $13 "\t" diff > out_large_diff
            }
            # Delete the key from the array to mark it as found.
            # At the end, only keys from file1 not in file2 will remain.
            delete prob[key1]
        } else {
            # Key is in file2 but was not in file1
            print key1 "\t" $13 > out_missing1
        }
    }
    # END block: runs after all lines from all files have been processed
    END {
        # Any keys remaining in the `prob` array are from file1
        # and were not found in file2.
        for (key in prob) {
            print key "\t" prob[key] > out_missing2
        }
    }
' "$f1" "$f2"
# --- Completion Message ---
echo "Input files:"
echo "  - file1: $f1"
echo "  - file2: $f2"
echo "Results written to:"
echo "  - $OUT_MATCH"
echo "  - $OUT_LARGE_DIFF"
echo "  - $OUT_MISSING_1"
echo "  - $OUT_MISSING_2"
