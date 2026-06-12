#!/bin/bash
# Build region fixtures by slicing a source BAM + VCF to the reads/variants
# overlapping a user region. Original chromosome names and absolute coordinates
# are preserved so the fixture can be loaded directly in IGV against the same
# reference. The reference used at test time is test/tmp/genome_chr<N>.fa,
# downloaded by test_var.sh (no per-region FASTA is emitted).
#
# Usage: SRC_BAM=... SRC_VCF=... SRC_REF=... ./test/regions/extract.sh [label ...]
# Without args, all regions in regions.tsv are processed; otherwise only the
# named labels are. SRC_REF is used only to generate expected outputs locally.

set -euo pipefail

SRC_BAM=${SRC_BAM:?must set SRC_BAM}
SRC_VCF=${SRC_VCF:?must set SRC_VCF}
SRC_REF=${SRC_REF:?must set SRC_REF}

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
MINIMOD=$ROOT/minimod
TSV=$ROOT/test/regions/regions.tsv
DATA_DIR=$ROOT/test/data/varmod
EXP_DIR=$ROOT/test/expected/varmod
PAD=${PAD:-500}
mkdir -p "$DATA_DIR" "$EXP_DIR"

want_label() {
    # called as: want_label "$label" "${WANTED[@]}"
    local L=$1; shift
    [ $# -eq 0 ] && return 0
    for arg in "$@"; do [ "$L" = "$arg" ] && return 0; done
    return 1
}

extract_one() {
    local label=$1 chrom=$2 start=$3 end=$4 notes=$5

    # query window — small pad so the slice picks up variants on the boundary
    local query_start=$((start - PAD))
    local query_end=$((end + PAD))
    [ $query_start -lt 1 ] && query_start=1

    local bam=$DATA_DIR/region_${label}.bam
    local vcf=$DATA_DIR/region_${label}.vcf

    # confirm at least one read overlaps the region; otherwise skip
    local nr_in_window
    nr_in_window=$(samtools view -c "$SRC_BAM" "${chrom}:${query_start}-${query_end}" 2>/dev/null || echo 0)
    if [ "$nr_in_window" = "0" ]; then
        echo "[extract] $label: no reads in $chrom:$query_start-$query_end — skipping"
        return
    fi

    echo "[extract] $label ${chrom}:${query_start}-${query_end}"

    # 1. BAM: slice to the region; preserve original chrom/POS and the full @SQ header.
    samtools view -bh "$SRC_BAM" "${chrom}:${query_start}-${query_end}" 2>/dev/null > "$bam"
    samtools index "$bam"

    # 2. VCF: keep header verbatim, filter records to the region.
    awk -v chrom="$chrom" -v ss=$query_start -v se=$query_end '
        /^#/ { print; next }
        $1 == chrom && $2 >= ss && $2 <= se { print }
    ' "$SRC_VCF" > "$vcf"

    # 3. Generate expected outputs against the slice + original reference.
    "$MINIMOD" varview --haplotypes -b -c "m,h" "$SRC_REF" "$bam" "$vcf" 2>/dev/null \
        > "$EXP_DIR/region_${label}.mh.varview.bedmethyl"
    "$MINIMOD" varfreq --haplotypes -b -c "m,h" "$SRC_REF" "$bam" "$vcf" 2>/dev/null \
        > "$EXP_DIR/region_${label}.mh.varfreq.bedmethyl"

    local nv=$(($(wc -l < "$vcf") - $(grep -c '^#' "$vcf")))
    local nr=$(samtools view -c "$bam" 2>/dev/null)
    local nx=$(wc -l < "$EXP_DIR/region_${label}.mh.varview.bedmethyl")
    local nf=$(wc -l < "$EXP_DIR/region_${label}.mh.varfreq.bedmethyl")
    echo "         vcf=${nv} reads=${nr} varview_rows=$nx varfreq_rows=$nf"
}

WANTED=("$@")

while IFS=$'\t' read -r label chrom start end notes; do
    [[ "$label" =~ ^#.* ]] && continue
    [ -z "$label" ] && continue
    if want_label "$label" "${WANTED[@]}"; then
        extract_one "$label" "$chrom" "$start" "$end" "$notes"
    fi
done < "$TSV"
