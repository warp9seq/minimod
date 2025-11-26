#!/usr/bin/env bash
#
# Usage: ./compare_files.sh file1.tsv file2.bedmethyl out_dir
#
# This script compares two files based on a constructed key and outputs
# the matched records, records with large differences, and records
# missing in either file.

# --- Configuration ---
f1="$1"
f2="$2"
out_dir="$3"
TOLERANCE=0.05 # Define the probability difference tolerance

# --- Argument Check ---
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 file1.tsv file2.bedmethyl out_dir"
    exit 1
fi

# Check if extension of files are correct. file1 should be .tsv and file2 should be .bedmethyl or .bed 
if [[ "$f1" != *.tsv ]] || [[ "$f2" != *.bedmethyl && "$f2" != *.bed ]]; then
    echo "Error: file1 should have .tsv extension and file2 should have .bedmethyl or .bed extension."
    exit 1
fi

# Extract base names without extension
base1=$(basename "$f1" .tsv)
base2=$(basename "$f2" .bedmethyl)
if [[ "$f2" == *.bed ]]; then
    base2=$(basename "$f2" .bed)
fi

mkdir -p "$out_dir"
script_name=$(basename "$0" .sh)

# --- Output Files ---
OUT_MISSING_1="$out_dir/missing_in_${base1}.tsv" # Keys in file2 but not in file1
OUT_MISSING_2="$out_dir/missing_in_${base2}.tsv" # Keys in file1 but not in file2
OUT_MATCH="$out_dir/in_both.tsv"              # Keys in both with small prob difference
OUT_LARGE_DIFF="$out_dir/large_prob_diff.tsv" # Keys in both with large prob difference

# # Remove output files if they already exist to avoid appending to old data. Ask for confirmation.
# if [[ -f "$OUT_MISSING_1" || -f "$OUT_MISSING_2" || -f "$OUT_MATCH" || -f "$OUT_LARGE_DIFF" ]]; then
#     read -p "Output files already exist. Overwrite? (y/n) " choice
#     if [[ "$choice" != "y" ]]; then
#         echo "Exiting without overwriting files."
#         exit 1
#     fi
# fi
rm -f "$OUT_MISSING_1" "$OUT_MISSING_2" "$OUT_MATCH" "$OUT_LARGE_DIFF"
touch "$OUT_MISSING_1" "$OUT_MISSING_2" "$OUT_MATCH" "$OUT_LARGE_DIFF"

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
        print "key\tprob1\tprob2\tdiff" > out_large_diff
    }

    # First pass: read file1 into an associative array `prob`
    # NR==FNR is true only while reading the first file
    NR==FNR {
        if (FNR == 1) next # Skip header
        for (i=1; i<=NF; i++) $i = clean($i)
        # Construct the key from the first 6 columns
        key = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
        prob[key] = $7
        next
    }

    # Skip the header of the second file
    FNR==1 { next }

    # Second pass: process file2 and compare against the `prob` array
    {
        for (i=1; i<=NF; i++) $i = clean($i)
        # Construct the key from file2 columns in the order that matches file1
        key1 = $4 "\t" $3 "\t" $7 "\t" $1 "\t" $2 "\t" $12

        if (key1 in prob) {
            # Key found in both files: calculate the difference
            diff = prob[key1] - $11
            if (diff < 0) diff = -diff # Absolute difference

            if (diff <= tol) {
                # Difference is within tolerance: write to matched.tsv
                print key1 "\t" prob[key1] "\t" diff > out_match
            } else {
                # Difference is too large: write to large_prob_diff.tsv
                print key1 "\t" prob[key1] "\t" $11 "\t" diff > out_large_diff
            }
            # Delete the key from the array to mark it as found.
            # At the end, only keys from file1 not in file2 will remain.
            delete prob[key1]
        } else {
            # Key is in file2 but was not in file1
            print key1 > out_missing1
        }
    }

    # END block: runs after all lines from all files have been processed
    END {
        # Any keys remaining in the `prob` array are from file1
        # and were not found in file2.
        for (key in prob) {
            print key > out_missing2
        }
    }
' "$f1" "$f2"

# --- Completion Message ---
echo "Results written to:"
echo "  - $OUT_MATCH"
echo "  - $OUT_LARGE_DIFF"
echo "  - $OUT_MISSING_1"
echo "  - $OUT_MISSING_2"