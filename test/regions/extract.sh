#!/bin/bash
# Build self-contained region fixtures by slicing a source BAM/VCF/ref.
# Each region from regions.tsv gets a renamed chromosome (=label) so the
# fixture has no dependency on the full reference.
#
# Usage: SRC_BAM=... SRC_VCF=... SRC_REF=... ./test/regions/extract.sh [label ...]
# Without args, all regions in regions.tsv are processed; otherwise only the
# named labels are.

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

    # query window (user start/end + small pad) to find reads that touch the region.
    local query_start=$((start - PAD))
    local query_end=$((end + PAD))
    [ $query_start -lt 1 ] && query_start=1

    # find the actual ALIGNMENT EXTENT of reads overlapping the query window.
    # Long nanopore reads can span tens of kb; the slice must cover the whole extent
    # of any read we keep, otherwise the BAM/CIGAR won't line up with the renamed chrom.
    local extent
    extent=$(samtools view "$SRC_BAM" "${chrom}:${query_start}-${query_end}" 2>/dev/null \
        | awk -v ws=$query_start -v we=$query_end '
            BEGIN{min=0; max=0}
            {
                cigar=$6; ref=0
                while (match(cigar, /[0-9]+[MDNSHI=X]/)) {
                    op=substr(cigar, RSTART, RLENGTH)
                    n=substr(op, 1, length(op)-1)+0
                    c=substr(op, length(op))
                    if (c=="M"||c=="D"||c=="N"||c=="="||c=="X") ref+=n
                    cigar=substr(cigar, RSTART+RLENGTH)
                }
                e=$4+ref-1
                if (min==0 || $4<min) min=$4
                if (e>max) max=e
            }
            END { if (min==0) {print "0 0"} else {print min, max} }')
    local slice_start=${extent% *}
    local slice_end=${extent#* }
    if [ "$slice_start" = "0" ]; then
        echo "[extract] $label: no reads in $chrom:$query_start-$query_end — skipping"
        return
    fi
    local slice_len=$((slice_end - slice_start + 1))

    local new_chrom="$label"
    local bam=$DATA_DIR/region_${label}.bam
    local vcf=$DATA_DIR/region_${label}.vcf
    local fa=$DATA_DIR/region_${label}.fa

    echo "[extract] $label ${chrom}:${slice_start}-${slice_end} (${slice_len} bp, renamed to $new_chrom)"

    # 1. FASTA: extract slice, rename chromosome
    samtools faidx "$SRC_REF" "${chrom}:${slice_start}-${slice_end}" \
        | awk -v new="$new_chrom" 'NR==1{print ">"new; next} {print}' > "$fa"
    samtools faidx "$fa"

    # 2. VCF: keep header but rewrite contig + records inside window
    awk -v chrom="$chrom" -v new="$new_chrom" -v ss=$slice_start -v se=$slice_end -v sl=$slice_len '
        BEGIN{OFS="\t"; printed_contig=0}
        /^##contig=<ID=/ {
            # only emit our one contig line, replace any matching original
            if (!printed_contig) { printf("##contig=<ID=%s,length=%d>\n", new, sl); printed_contig=1 }
            next
        }
        /^#/ { print; next }
        $1 != chrom { next }
        {
            pos = $2
            if (pos < ss || pos > se) next
            $1 = new
            $2 = pos - ss + 1
            print
        }' "$SRC_VCF" > "$vcf"

    # 3. BAM: rewrite @SQ and shift RNAME/POS/RNEXT/PNEXT.
    # Query with the SAME window used for extent calc so all reads' POS >= slice_start.
    samtools view -h "$SRC_BAM" "${chrom}:${query_start}-${query_end}" 2>/dev/null \
        | awk -v chrom="$chrom" -v new="$new_chrom" -v ss=$slice_start -v sl=$slice_len '
            BEGIN{OFS="\t"; emitted_sq=0}
            /^@HD/ { print; next }
            /^@SQ/ {
                if (!emitted_sq) { printf("@SQ\tSN:%s\tLN:%d\n", new, sl); emitted_sq=1 }
                next
            }
            /^@/ { print; next }
            $3 == chrom {
                $3 = new
                $4 = $4 - ss + 1
                if ($7 == chrom) $7 = new
                if ($7 == new || $7 == "=") $8 = $8 - ss + 1
                print
            }' > "${bam}.sam"
    samtools view -b "${bam}.sam" 2>/dev/null > "${bam}.unsorted"
    samtools sort -o "$bam" "${bam}.unsorted" 2>/dev/null
    rm -f "${bam}.sam" "${bam}.unsorted"
    samtools index "$bam"

    # 4. Generate expected outputs against the slice
    "$MINIMOD" varview --haplotypes -c "m,h" "$fa" "$bam" "$vcf" 2>/dev/null \
        > "$EXP_DIR/region_${label}.mh.varview.tsv"
    "$MINIMOD" varfreq --haplotypes -b -c "m,h" "$fa" "$bam" "$vcf" 2>/dev/null \
        > "$EXP_DIR/region_${label}.mh.varfreq.bedmethyl"

    local nv=$(($(wc -l < "$vcf") - $(grep -c '^#' "$vcf")))
    local nr=$(samtools view -c "$bam" 2>/dev/null)
    local nx=$(wc -l < "$EXP_DIR/region_${label}.mh.varview.tsv")
    local nf=$(wc -l < "$EXP_DIR/region_${label}.mh.varfreq.bedmethyl")
    echo "         vcf=${nv} reads=${nr} varview_rows=$((nx - 1)) varfreq_rows=$nf"
}

WANTED=("$@")

# read TSV, skip comments + blanks
while IFS=$'\t' read -r label chrom start end notes; do
    [[ "$label" =~ ^#.* ]] && continue
    [ -z "$label" ] && continue
    if want_label "$label" "${WANTED[@]}"; then
        extract_one "$label" "$chrom" "$start" "$end" "$notes"
    fi
done < "$TSV"
