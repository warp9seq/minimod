#!/usr/bin/env python3

# Usage: test/varfreq_context.py freq.bedmethyl varfreq.bedmethyl min_window [--bed out.bed]

import sys
import gzip
from bisect import bisect_left
from statistics import mean, stdev

argv = sys.argv[1:]
bed_file = None
if "--bed" in argv:
    i = argv.index("--bed")
    if i + 1 >= len(argv):
        print("--bed requires an output file path")
        sys.exit(1)
    bed_file = argv[i + 1]
    del argv[i:i + 2]

if len(argv) != 3:
    print("Usage: {} freq.bedmethyl varfreq.bedmethyl min_window [--bed out.bed]".format(sys.argv[0]))
    sys.exit(1)

ALL_VARTYPES = ("SNP", "INS", "DEL", "MNP")

freq_file = argv[0]
varfreq_file = argv[1]
min_window = int(argv[2])

def var_type_of(ref_allele, alt_allele):
    if len(ref_allele) == 1 and len(alt_allele) == 1:
        return "SNP"
    if len(ref_allele) > len(alt_allele):
        return "DEL"
    if len(ref_allele) < len(alt_allele):
        return "INS"
    return "MNP"

def site_pct(strand_pcts):
    return mean(strand_pcts)

def site_pcts(sites):
    return [site_pct(strand_pcts) for strand_pcts in sites]

def stats(pcts):
    if not pcts:
        return (0, None, None)
    return (len(pcts), mean(pcts), stdev(pcts) if len(pcts) > 1 else None)

def fmt(x):
    return "NA" if x is None else "{:.1f}".format(x)

def freq_load(fn):
    # (contig, hap, mod_code) -> (sites, [(site, n_called, n_mod), ...])
    per_site = {}
    opener = gzip.open if fn.endswith(".gz") else open
    with opener(fn, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig, mod_code, strand, hap = parts[0], parts[3], parts[5], parts[11]
            pos, n_called, freq = int(parts[1]), int(parts[4]), float(parts[10])

            if hap == "*":
                continue

            site_pos = pos if strand == "+" else pos - 1
            d = per_site.setdefault((contig, hap, mod_code), {})
            d.setdefault(site_pos, []).append(freq)

    out = {}
    for key, d in per_site.items():
        sites = sorted(d)
        out[key] = (sites, [site_pct(d[site]) for site in sites])
    return out

def cpg_key_offset(pos, offset, strand, var_pos, ref_len, alt_len):
    if pos == var_pos and offset >= 1:
        o = offset + 1
    elif offset == ref_len:
        o = alt_len + 1
    else:
        o = offset + 1
    return o if strand == "+" else o - 1

def varfreq_load(fn):
    groups = {}
    opener = gzip.open if fn.endswith(".gz") else open
    with opener(fn, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig, mod_code, strand, hap = parts[0], parts[3], parts[5], parts[17]
            pos, n_called, freq = int(parts[1]), int(parts[4]), float(parts[10])
            var_pos, ref_allele, alt_allele = int(parts[12]), parts[14], parts[15]
            offset_str = parts[16]
            if offset_str == "*":
                continue
            offset = int(parts[16])

            if hap == "*":
                continue

            cpg_key = cpg_key_offset(pos, offset, strand, var_pos, len(ref_allele), len(alt_allele))
            # keep only CpGs whose C lies inside the ALT allele
            if not 1 <= cpg_key <= len(alt_allele):
                continue

            # check if CG is in the ALT allele
            if alt_allele[cpg_key - 1:cpg_key + 1].upper() != "CG":
                continue

            key = (contig, var_pos, ref_allele, alt_allele, hap, mod_code)
            g = groups.get(key)
            if g is None:
                g = [var_type_of(ref_allele, alt_allele), {}]
                groups[key] = g
            g[1].setdefault(cpg_key, []).append(freq)
    return groups

def background(freqs, contig, hap, mod_code, var_pos, var_end, radius):
    entry = freqs.get((contig, hap, mod_code))
    if entry is None:
        return ([], [])
    sites, pcts = entry
    lo = bisect_left(sites, var_pos - radius)
    mid_l = bisect_left(sites, var_pos)
    mid_r = bisect_left(sites, var_end)
    hi = bisect_left(sites, var_end + radius)
    return (pcts[lo:mid_l], pcts[mid_r:hi])

CATEGORIES = ("meth_in_unmeth", "unmeth_in_meth", "equal")
COLOUR = {"meth_in_unmeth": "255,0,0", "unmeth_in_meth": "0,0,255", "equal": "128,128,128"}

freqs = freq_load(freq_file)
varfreqs = varfreq_load(varfreq_file)

print("# freq    : {}".format(freq_file), file=sys.stderr)
print("# varfreq : {}".format(varfreq_file), file=sys.stderr)
print("# min_window={} bp  per-CpG freq = mean of its + and - strand freqs  allele level = unweighted mean over the CpGs in ALT".format(
    min_window), file=sys.stderr)
print("# meth_in_unmeth if var_meth_mean > both flank means  unmeth_in_meth if a flank mean > var_meth_mean  equal if it ties the larger flank mean", file=sys.stderr)
print("# meth_mean_diff = var_meth_mean - the larger flank mean, so its sign follows the category", file=sys.stderr)

def new_bucket():
    return {c: [] for c in CATEGORIES}
buckets = {t: new_bucket() for t in ALL_VARTYPES}
skipped_no_bg = {t: 0 for t in ALL_VARTYPES}

for key in sorted(varfreqs):
    contig, var_pos, ref_allele, alt_allele, hap, mod_code = key
    var_type, var_sites = varfreqs[key]
    var_end = var_pos + len(ref_allele)

    var_ncpg, var_mean, var_sd = stats(site_pcts(var_sites.values()))

    radius = max(min_window, len(alt_allele))/2

    bg_l_pcts, bg_r_pcts = background(freqs, contig, hap, mod_code, var_pos, var_end, radius)
    bg_l_ncpg, bg_l_mean, _ = stats(bg_l_pcts)
    bg_r_ncpg, bg_r_mean, _ = stats(bg_r_pcts)
    bg_ncpg, bg_mean, bg_sd = stats(bg_l_pcts + bg_r_pcts)

    flanks = [m for m in (bg_l_mean, bg_r_mean) if m is not None]
    if not flanks:
        skipped_no_bg[var_type] += 1
        continue

    meth_diff = var_mean - bg_mean

    rec = (contig, var_pos, ref_allele, alt_allele, hap, mod_code,
           var_ncpg, var_mean, var_sd,
           bg_ncpg, bg_mean, bg_sd, bg_l_ncpg, bg_l_mean, bg_r_ncpg, bg_r_mean, meth_diff, radius)

    if meth_diff > 0:
        buckets[var_type]["meth_in_unmeth"].append(rec)
    elif meth_diff < 0:
        buckets[var_type]["unmeth_in_meth"].append(rec)
    else:
        buckets[var_type]["equal"].append(rec)

COLUMNS = ["chrom", "start", "end", "var_type", "ref", "alt", "hap", "mod", "category",
           "alt_len", "var_ncpg", "var_meth_mean", "var_meth_sd",
           "bg_span", "bg_ncpg", "bg_meth_mean", "bg_meth_sd",
           "bg_l_ncpg", "bg_l_meth_mean", "bg_r_ncpg", "bg_r_meth_mean", "meth_mean_diff"]

tsv_rows = []
for var_type in ALL_VARTYPES:
    b = buckets[var_type]
    for category in CATEGORIES:
        for (contig, var_pos, ref_allele, alt_allele, hap, mod_code,
             var_ncpg, var_mean, var_sd,
             bg_ncpg, bg_mean, bg_sd, bg_l_ncpg, bg_l_mean, bg_r_ncpg, bg_r_mean,
             meth_diff, radius) in b[category]:
            tsv_rows.append((contig, var_pos, var_pos + len(ref_allele), var_type,
                             ref_allele, alt_allele, hap, mod_code, category,
                             len(alt_allele), var_ncpg, fmt(var_mean), fmt(var_sd),
                             2 * radius, bg_ncpg, fmt(bg_mean), fmt(bg_sd),
                             bg_l_ncpg, fmt(bg_l_mean), bg_r_ncpg, fmt(bg_r_mean),
                             fmt(meth_diff)))

tsv_rows.sort(key=lambda r: (r[0], r[1], r[2], r[8]))

print("\t".join(COLUMNS))
for r in tsv_rows:
    print("\t".join(str(f) for f in r))

for var_type in ALL_VARTYPES:
    b = buckets[var_type]
    print("# {}: meth_in_unmeth={} unmeth_in_meth={} equal={} no_bg_cpg={}".format(
        var_type, len(b["meth_in_unmeth"]), len(b["unmeth_in_meth"]), len(b["equal"]),
        skipped_no_bg[var_type]), file=sys.stderr)

def write_bed(fn, buckets):
    rows = []
    for var_type in ALL_VARTYPES:
        for category in CATEGORIES:
            for (contig, var_pos, ref_allele, alt_allele, hap, mod_code,
                 var_ncpg, var_mean, var_sd,
                 bg_ncpg, bg_mean, bg_sd, bg_l_ncpg, bg_l_mean, bg_r_ncpg, bg_r_mean,
                 meth_diff, radius) in buckets[var_type][category]:
                start = var_pos
                end = var_pos + len(ref_allele)
                name = "{}_{}_{}>{}_hap{}_{}_v{:.0f}/b{}".format(
                    var_type, category, ref_allele, alt_allele, hap, mod_code, var_mean,
                    "NA" if bg_mean is None else "{:.0f}".format(bg_mean))
                score = min(1000, int(round(var_mean * 10)))
                rows.append((contig, start, end, name, score, COLOUR[category]))

    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(fn, "w") as out:
        out.write('track name="varfreq_context" description="variant vs background methylation" itemRgb="On"\n')
        for contig, start, end, name, score, colour in rows:
            out.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\n".format(
                contig, start, end, name, score, start, end, colour))
    return len(rows)

if bed_file is not None:
    n = write_bed(bed_file, buckets)
    print("# wrote {} feature(s) to {}".format(n, bed_file), file=sys.stderr)
