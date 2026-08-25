#!/usr/bin/env python3

# Usage: test/varfreq_average.py varfreq.bedmethyl

import sys
import gzip
from statistics import mean

if len(sys.argv) != 2:
    print("Usage: {} varfreq.bedmethyl".format(sys.argv[0]))
    sys.exit(1)

varfreq_file = sys.argv[1]

def c_index(pos, offset, strand, var_pos, ref_len, alt_len):

    if pos == var_pos and offset >= 1:
        o_val = offset + 1
    elif offset == ref_len:
        o_val = alt_len + 1
    else:
        o_val = offset + 1
    if strand != "+":
        o_val -= 1
    return o_val - 1

variants = {}

opener = gzip.open if varfreq_file.endswith(".gz") else open
with opener(varfreq_file, "rt") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        contig, mod_code, strand = parts[0], parts[3], parts[5]
        pos, freq = int(parts[1]), float(parts[10])
        var_pos, gt, ref_allele, alt_allele = int(parts[12]), parts[13], parts[14], parts[15]
        offset = parts[16]
        hap = parts[17] if len(parts) > 17 else "*"

        if hap == "*" or offset == "*":
            continue

        c_idx = c_index(pos, int(offset), strand, var_pos, len(ref_allele), len(alt_allele))

        key = (contig, var_pos, ref_allele, alt_allele, gt, mod_code)
        variants.setdefault(key, {}).setdefault(hap, {}).setdefault(c_idx, []).append(freq)

print("\t".join(["chrom", "start", "end", "ref", "alt", "gt", "mod", "hap", "n_cpg", "meth_mean"]))

for key in variants:
    contig, var_pos, ref_allele, alt_allele, gt, mod_code = key
    haps = variants[key]
    prefix = [contig, var_pos, var_pos + len(ref_allele), ref_allele, alt_allele, gt, mod_code]

    hap_means = []
    for hap in haps:
        sites = haps[hap]
        meth_mean = mean([mean(freqs) for freqs in sites.values()])
        hap_means.append(meth_mean)
        print("\t".join(str(f) for f in prefix + [hap, len(sites), "{:.1f}".format(meth_mean)]))

    if len(haps) > 1:
        n_cpg = len(set(c_idx for sites in haps.values() for c_idx in sites))
        print("\t".join(str(f) for f in prefix + ["*", n_cpg, "{:.1f}".format(mean(hap_means))]))
