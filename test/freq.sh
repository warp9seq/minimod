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

awk -F'\t' -v mod_code="$MOD_CODE" -v threshold="$THRESHOLD" '
BEGIN {
    OFS = "\t"
    lower_threshold = 1.0 - threshold
    print "contig", "ref_pos", "strand", "mod_code", "frequency"
}

# Skip the header row.
NR == 1 {
    next
}

{
    ref_contig = $1
    ref_pos = $2
    strand = $3
    rec_mod_code = $6
    mod_prob = $7 + 0

    if (ref_contig == "" || rec_mod_code != mod_code) {
        next
    }

    key = ref_contig OFS ref_pos OFS strand OFS rec_mod_code

    if (mod_prob >= threshold) {
        n_mod[key]++
        n_called[key]++
    } else if (mod_prob <= lower_threshold) {
        n_called[key]++
    }
}

END {
    for (key in n_called) {
        if (n_called[key] > 0) {
            split(key, fields, OFS)
            printf "%s\t%s\t%s\t%s\t%.2f\n", fields[1], fields[2], fields[3], fields[4], n_mod[key] / n_called[key]
        }
    }
}
' "$FILE"

