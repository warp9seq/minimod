#!/bin/bash

# usage: bash validate.sh <minimod_varfreq_1.bed> <minimod_varfreq_2.bed>
# output: 
# <minimod_varfreq_1.bed> :<number of modified sites>
# <minimod_varfreq_2.bed> :<number of modified sites>
# union :<number of unique sites in the union of the two files>
# intersection :<number of unique sites in the intersection of the two files> (percentage intersection/union)

# bed format (contig, start, end, mod, n_called, strand, start, end, color, m_called, freq, var_pos, ref_base, alt_base, cg_offset)
# chr22	11211064	11211065	m	7	+	11211064	11211065	255,0,0	7	0.000000	11211064	T	C	0
# chr22	11211065	11211066	m	4	+	11211065	11211066	255,0,0	4	0.000000	11211064	T	C	0
# chr22	11211065	11211066	m	3	-	11211065	11211066	255,0,0	3	100.000000	11211064	T	C	0

# intersetct should match contig, start, end, mod, strand, ref_base, cg_offset

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: bash validate.sh <minimod_varfreq_1.bed> <minimod_varfreq_2.bed>" >&2
    exit 1
fi

f1="$1"
f2="$2"

for f in "$f1" "$f2"; do
    if [ ! -f "$f" ]; then
        echo "error: file not found: $f" >&2
        exit 1
    fi
done

# Key columns matching intersection criteria:
# contig(1), start(2), end(3), mod(4), strand(6), ref_base(13), cg_offset(15)
key1=$(mktemp)
key2=$(mktemp)
trap 'rm -f "$key1" "$key2"' EXIT

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$6,$13,$15}' "$f1" | sort -u > "$key1"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$6,$13,$15}' "$f2" | sort -u > "$key2"

n1=$(wc -l < "$key1")
n2=$(wc -l < "$key2")
union=$(sort -u "$key1" "$key2" | wc -l)
intersection=$(comm -12 "$key1" "$key2" | wc -l)

if [ "$union" -gt 0 ]; then
    pct=$(awk -v i="$intersection" -v u="$union" 'BEGIN{printf "%.2f", (i/u)*100}')
else
    pct="0.00"
fi

echo "${f1} :${n1}"
echo "${f2} :${n2}"
echo "union :${union}"
echo "intersection :${intersection} (${pct}%)"

