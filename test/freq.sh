#!/bin/bash

# Author: Suneth Samarsinghe (2026)

# This script calculates the frequency of modified bases from a minimod view.tsv file.
# How:
# Aggregate records in view.tsv to a map with key->value pairs : (ref_contig, ref_pos, strand, mod_code) -> (n_called, n_mod)
# if mod_prob >= threshold, then n_mod += 1 and n_called += 1
# if mod_prob < 1 - threshold, then n_called += 1
# else, ignore the record (ambiguous)
# Finally, calculate the frequency as n_mod / (n_mod+n_called) for each key and print the results

# Sample minimod view.tsv file:
# ref_contig	ref_pos	strand	read_id	read_pos	mod_code	mod_prob
# chr22	19979864	+	m84088_230609_030819_s1/55512555/ccs	14	m	0.708984
# chr22	19979882	+	m84088_230609_030819_s1/55512555/ccs	32	m	0.947266

die() {
    echo "$1" >&2
    exit 1
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    die "Usage: $0 <mod_code> <thresh> minimod_view.tsv"
fi

MOD_CODE="$1"
THRESHOLD="$2"
FILE="$3"

if [ ! -f "$FILE" ]; then
    die "File not found: $FILE"
fi

LOWER_THRESHOLD=$(echo "(255.5/256.0) - $THRESHOLD" | bc -l)

declare -A freq_map

{
    # Skip the header line and process the rest.
    read -r _
    while IFS=$'\t' read -r ref_contig ref_pos strand read_id read_pos mod_code mod_prob; do
        if [ -z "$ref_contig" ] || [ "$mod_code" != "$MOD_CODE" ]; then
            continue
        fi

        key="$ref_contig"$'\t'"$ref_pos"$'\t'"$strand"$'\t'"$mod_code"

        if [ -z "${freq_map[$key]}" ]; then
            freq_map["$key"]="0:0"
        fi

        IFS=':' read -r n_mod n_called <<< "${freq_map[$key]}"

        if (( $(echo "$mod_prob >= $THRESHOLD" | bc -l) )); then
            n_mod=$((n_mod + 1))
            n_called=$((n_called + 1))
        elif (( $(echo "$mod_prob <= $LOWER_THRESHOLD" | bc -l) )); then
            n_called=$((n_called + 1))
        else
            continue
        fi

        freq_map["$key"]="$n_mod:$n_called"
    done
} < "$FILE"


echo -e "contig\tref_pos\tstrand\tmod_code\tfrequency"
for key in "${!freq_map[@]}"; do
    IFS=':' read -r n_mod n_called <<< "${freq_map[$key]}"
    if [ "$n_called" -gt 0 ]; then
        frequency=$(echo "scale=2; $n_mod / $n_called" | bc -l)
        echo "$key	$frequency"
    fi
done

