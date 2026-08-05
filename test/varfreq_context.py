#!/usr/bin/env python3

# Usage: test/varfreq_context.py freq.bedmethyl varfreq.bedmethyl radius meth_freq_thresh var_meth_pct var_unmeth_pct bg_meth_pct bg_unmeth_pct [--bed out.bed]

import sys
import gzip
from bisect import bisect_left, bisect_right

argv = sys.argv[1:]
bed_file = None
if "--bed" in argv:
    i = argv.index("--bed")
    if i + 1 >= len(argv):
        print("--bed requires an output file path")
        sys.exit(1)
    bed_file = argv[i + 1]
    del argv[i:i + 2]

if len(argv) != 8:
    print("Usage: {} freq.bedmethyl varfreq.bedmethyl radius meth_freq_thresh var_meth_pct var_unmeth_pct bg_meth_pct bg_unmeth_pct [--bed out.bed]".format(sys.argv[0]))
    sys.exit(1)

ALL_VARTYPES = ("SNP", "INS", "DEL", "MNP")

freq_file = argv[0]
varfreq_file = argv[1]
radius = int(argv[2])
meth_freq_thresh = float(argv[3])
var_meth_pct = float(argv[4])
var_unmeth_pct = float(argv[5])
bg_meth_pct = float(argv[6])
bg_unmeth_pct = float(argv[7])

def var_type_of(ref_allele, alt_allele):
    if len(ref_allele) == 1 and len(alt_allele) == 1:
        return "SNP"
    if len(ref_allele) > len(alt_allele):
        return "DEL"
    if len(ref_allele) < len(alt_allele):
        return "INS"
    return "MNP"

def freq_load(fn):
    # (contig, hap, mod_code) -> (sites, [(site, n_called, n_mod), ...])
    per_site = {}
    opener = gzip.open if fn.endswith(".gz") else open
    with opener(fn, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig, mod_code, strand, hap = parts[0], parts[3], parts[5], parts[11]
            pos, n_called, freq = int(parts[1]), int(parts[4]), float(parts[10])

            n_mod = n_called * freq / 100.0

            if hap == "*":
                continue

            site_pos = pos if strand == "+" else pos - 1
            d = per_site.setdefault((contig, hap, mod_code), {})
            s = d.get(site_pos)
            if s:
                s[0] += n_called
                s[1] += n_mod
            else:
                d[site_pos] = [n_called, n_mod]

    out = {}
    for key, d in per_site.items():
        sites = sorted(d)
        out[key] = (sites, [(site, d[site][0], d[site][1]) for site in sites])
    return out

def varfreq_load(fn):
    groups = {}
    opener = gzip.open if fn.endswith(".gz") else open
    with opener(fn, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig, mod_code, strand, hap = parts[0], parts[3], parts[5], parts[16]
            pos, n_called, freq = int(parts[1]), int(parts[4]), float(parts[10])
            var_pos, ref_allele, alt_allele = int(parts[11]), parts[13], parts[14]

            n_mod = n_called * freq / 100.0

            if hap == "*":
                continue

            site_pos = pos if strand == "+" else pos - 1
            key = (contig, var_pos, ref_allele, alt_allele, hap, mod_code)
            g = groups.get(key)
            if g is None:
                g = [var_type_of(ref_allele, alt_allele), {}]
                groups[key] = g
            s = g[1].get(site_pos)
            if s:
                s[0] += n_called
                s[1] += n_mod
            else:
                g[1][site_pos] = [n_called, n_mod]
    return groups

def meth_site_count(sites_data):
    return sum(1 for n_called, n_mod in sites_data
               if n_called > 0 and n_mod / n_called > meth_freq_thresh)

def background(freqs, contig, hap, mod_code, var_pos):
    entry = freqs.get((contig, hap, mod_code))
    if entry is None:
        return ((0, 0), (0, 0))
    sites, records = entry
    lo = bisect_left(sites, var_pos - radius)
    mid_l = bisect_left(sites, var_pos)
    mid_r = bisect_right(sites, var_pos)
    hi = bisect_right(sites, var_pos + radius)
    left_meth = meth_site_count((records[i][1], records[i][2]) for i in range(lo, mid_l))
    right_meth = meth_site_count((records[i][1], records[i][2]) for i in range(mid_r, hi))
    return ((left_meth, mid_l - lo), (right_meth, hi - mid_r))

CATEGORIES = ("meth_in_unmeth", "unmeth_in_meth")
COLOUR = {"meth_in_unmeth": "255,0,0", "unmeth_in_meth": "0,0,255"}

freqs = freq_load(freq_file)
varfreqs = varfreq_load(varfreq_file)

print("# freq    : {}".format(freq_file), file=sys.stderr)
print("# varfreq : {}".format(varfreq_file), file=sys.stderr)
print("# radius={} bp  site methylated if freq>{}  variant methylated if meth/len >{}%  variant unmethylated if meth/len <{}%  background methylated if meth/len >{}%  background unmethylated if meth/len <{}%".format(
    radius, meth_freq_thresh, var_meth_pct, var_unmeth_pct, bg_meth_pct, bg_unmeth_pct), file=sys.stderr)

def new_bucket():
    return {c: [] for c in CATEGORIES}
buckets = {t: new_bucket() for t in ALL_VARTYPES}
ignored = {t: {"meth_in_meth": 0, "unmeth_in_unmeth": 0, "ambiguous": 0} for t in ALL_VARTYPES}

for key in sorted(varfreqs):
    contig, var_pos, ref_allele, alt_allele, hap, mod_code = key
    var_type, var_sites = varfreqs[key]
    b = buckets[var_type]

    var_ncpg = len(var_sites)
    (bg_l_nmeth, bg_l_ncpg), (bg_r_nmeth, bg_r_ncpg) = background(freqs, contig, hap, mod_code, var_pos)
    bg_nmeth = bg_l_nmeth + bg_r_nmeth
    bg_ncpg = bg_l_ncpg + bg_r_ncpg

    var_nmeth = meth_site_count(var_sites.values())
    var_pct = 100.0 * var_nmeth / len(alt_allele)
    bg_l_pct = 100.0 * bg_l_nmeth / radius
    bg_r_pct = 100.0 * bg_r_nmeth / radius
    bg_pct = 100.0 * bg_nmeth / (2 * radius)
    var_methylated = var_pct > var_meth_pct
    var_unmethylated = var_pct < var_unmeth_pct
    bg_methylated = bg_l_pct > bg_meth_pct or bg_r_pct > bg_meth_pct
    bg_unmethylated = bg_l_pct < bg_unmeth_pct and bg_r_pct < bg_unmeth_pct

    rec = (contig, var_pos, ref_allele, alt_allele, hap, mod_code,
           var_pct, var_nmeth, var_ncpg, bg_pct, bg_l_pct, bg_r_pct, bg_nmeth, bg_ncpg)

    if var_methylated and bg_unmethylated:
        b["meth_in_unmeth"].append(rec)
    elif var_unmethylated and bg_methylated:
        b["unmeth_in_meth"].append(rec)
    elif var_methylated and bg_methylated:
        ignored[var_type]["meth_in_meth"] += 1
    elif var_unmethylated and bg_unmethylated:
        ignored[var_type]["unmeth_in_unmeth"] += 1
    else:
        ignored[var_type]["ambiguous"] += 1

COLUMNS = ["chrom", "start", "end", "var_type", "ref", "alt", "hap", "mod",
           "category", "var_meth", "var_len", "var_meth_pct", "var_cpg",
           "bg_meth", "bg_len", "bg_meth_pct", "bg_l_meth_pct", "bg_r_meth_pct", "bg_cpg"]

tsv_rows = []
for var_type in ALL_VARTYPES:
    b = buckets[var_type]
    for category in CATEGORIES:
        for (contig, var_pos, ref_allele, alt_allele, hap, mod_code,
             var_pct, var_nmeth, var_ncpg, bg_pct, bg_l_pct, bg_r_pct, bg_nmeth, bg_ncpg) in b[category]:
            tsv_rows.append((contig, var_pos, var_pos + len(ref_allele), var_type,
                             ref_allele, alt_allele, hap, mod_code, category,
                             var_nmeth, len(alt_allele), var_pct, var_ncpg,
                             bg_nmeth, 2 * radius, bg_pct, bg_l_pct, bg_r_pct, bg_ncpg))

tsv_rows.sort(key=lambda r: (r[0], r[1], r[2], r[8]))

print("\t".join(COLUMNS))
for r in tsv_rows:
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{:.1f}\t{}".format(*r))

for var_type in ALL_VARTYPES:
    b = buckets[var_type]
    ig = ignored[var_type]
    print("# {}: meth_in_unmeth={} unmeth_in_meth={} meth_in_meth={} unmeth_in_unmeth={} ambiguous={}".format(
        var_type, len(b["meth_in_unmeth"]), len(b["unmeth_in_meth"]),
        ig["meth_in_meth"], ig["unmeth_in_unmeth"], ig["ambiguous"]), file=sys.stderr)

def write_bed(fn, buckets):
    rows = []
    for var_type in ALL_VARTYPES:
        for category in CATEGORIES:
            for (contig, var_pos, ref_allele, alt_allele, hap, mod_code,
                 var_pct, var_nmeth, var_ncpg, bg_pct, bg_l_pct, bg_r_pct, bg_nmeth, bg_ncpg) in buckets[var_type][category]:
                start = var_pos
                end = var_pos + len(ref_allele)
                name = "{}_{}_{}>{}_hap{}_{}_v{:.0f}/b{:.0f}".format(
                    var_type, category, ref_allele, alt_allele, hap, mod_code, var_pct, bg_pct)
                score = min(1000, int(round(var_pct * 10)))
                rows.append((contig, start, end, name, score, COLOUR[category]))

    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(fn, "w") as out:
        out.write('track name="varfreq_context" description="variant vs background methylation contrasts" itemRgb="On"\n')
        for contig, start, end, name, score, colour in rows:
            out.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\n".format(
                contig, start, end, name, score, start, end, colour))
    return len(rows)

if bed_file is not None:
    n = write_bed(bed_file, buckets)
    print("# wrote {} contrast feature(s) to {}".format(n, bed_file), file=sys.stderr)
