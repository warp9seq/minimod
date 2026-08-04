#!/usr/bin/env python3

# Usage: test/varfreq_summary.py variant.vcf minimod.freq.haplotagged.bedmethyl minimod.varfreq.haplotagged.bedmethyl [reference.fa]
#
# The optional reference.fa enables the reference-CpG leak check (see
# reference_cpg_in_varfreq_check). Without it that check is skipped.

import sys
import gzip

if len(sys.argv) not in (4, 5):
    print("Usage: {} variant.vcf minimod.freq.haplotagged.bedmethyl minimod.varfreq.haplotagged.bedmethyl [reference.fa]".format(sys.argv[0]))
    sys.exit(1)

def vcf_load(fn):
    d = {}
    opener = gzip.open if fn.endswith(".gz") else open
    with opener(fn, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            contig = parts[0]
            pos = int(parts[1])
            ref_allele = parts[3]
            alt_alleles = parts[4].split(",")
            n_alt = len(alt_alleles)
            gt = parts[9].strip().split(":")[0]

            for alt_allele in alt_alleles:
                var_type = ""
                if alt_allele == ".":
                    var_type = "REF"
                elif(len(ref_allele) == 1 and len(alt_allele) == 1):
                    var_type = "SNP"
                elif(len(ref_allele) == len(alt_allele)):
                    var_type = "MNP"
                elif(len(ref_allele) > len(alt_allele)):
                    var_type = "DEL"
                elif(len(ref_allele) < len(alt_allele)):
                    var_type = "INS"
                else:
                    print("WARNING: skipping {}:{} because of unknown variant type {} vs {}".format(contig, pos, ref_allele, alt_allele), file=sys.stderr)
                    continue
                d[(contig, pos)] = (var_type, ref_allele, alt_allele, gt)
    return d

def freq_load(fn):
    d = {}
    opener = gzip.open if fn.endswith(".gz") else open
    with opener(fn, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")

            contig = parts[0]
            pos = int(parts[1])
            mod_code = parts[3]
            n_called = int(parts[4])
            strand = parts[5]
            freq = float(parts[10])
            n_mod = n_called * freq / 100.0
            hap = parts[11]

            if hap == "*":
                continue;

            if strand == "+":
                key = (contig, pos, mod_code, hap)
            else:
                key = (contig, pos - 1, mod_code, hap)

            if key in d:
                v = d[key]
                d[key] = (v[0] + n_called, v[1] + n_mod)
            else:
                d[key] = (n_called, n_mod)
    return d

def varfreq_load(fn):
    d = {}
    opener = gzip.open if fn.endswith(".gz") else open
    with opener(fn, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")

            contig = parts[0]
            pos = int(parts[1])
            mod_code = parts[3]
            n_called = int(parts[4])
            strand = parts[5]
            freq = float(parts[10])
            n_mod = n_called * freq / 100.0
            var_pos = int(parts[11])
            gt = parts[12]
            ref_allele = parts[13]
            alt_allele = parts[14]
            offset = int(parts[15])
            hap = parts[16]

            if hap == "*":
                continue;
            
            var_type = ""
            if(len(ref_allele) == 1 and len(alt_allele) == 1):
                var_type = "SNP"
            elif(len(ref_allele) > len(alt_allele)):
                var_type = "DEL"
            elif(len(ref_allele) < len(alt_allele)):
                var_type = "INS"
            elif(len(ref_allele) == len(alt_allele)):
                var_type = "MNP"
            else:
                print("WARNING: skipping {}:{} because of unknown variant type {} vs {}".format(contig, pos, ref_allele, alt_allele), file=sys.stderr)
                continue
                
            if strand == "+":
                key = (contig, pos, mod_code, hap)
            else:
                key = (contig, pos - 1, mod_code, hap)

            if key in d:
                vars = d[key]
                if var_pos in vars:
                    v = vars[var_pos]
                    d[key][var_pos] = (var_type, ref_allele, alt_allele, offset, gt, v[5] + n_called, v[6] + n_mod)
                else:
                    d[key][var_pos] = (var_type, ref_allele, alt_allele, offset, gt, n_called, n_mod)
            else:
                v = {}
                v[var_pos] = (var_type, ref_allele, alt_allele, offset, gt, n_called, n_mod)
                d[key] = v
    return d

def vcf_variant_summary(variants):
    summary = {}
    for key in variants:
        var_type = variants[key][0]
        if var_type not in summary:
            summary[var_type] = 0
        summary[var_type] += 1
    
    print("## Variant summary:")
    total=0
    for var_type in summary:
        print("{}: {}".format(var_type, summary[var_type]))
        total += summary[var_type]
    print("Total: {}".format(total))

def vcf_gt_summary(variants):
    summary = {}
    for key in variants:
        gt = variants[key][3]
        if gt not in summary:
            summary[gt] = 0
        summary[gt] += 1

    print("## Genotype summary:")
    total=0
    for gt in summary:
        print("{}: {}".format(gt, summary[gt]))
        total += summary[gt]
    print("Total: {}".format(total))

def freq_mod_summary(freqs):
    mod_summary = {}
    for key in freqs:
        mod_code = key[2]
        if mod_code not in mod_summary:
            mod_summary[mod_code] = 0
        mod_summary[mod_code] += 1
    
    print("## Modification summary:")
    total=0
    for mod_code in mod_summary:
        print("{}: {}".format(mod_code, mod_summary[mod_code]))
        total += mod_summary[mod_code]
    print("Total: {}".format(total))

def freq_hap_summary(freqs):
    hap_summary = {}
    for key in freqs:
        hap = key[3]
        if hap not in hap_summary:
            hap_summary[hap] = 0
        hap_summary[hap] += 1
    
    print("## Haplotype summary:")
    total=0
    for hap in hap_summary:
        print("{}: {}".format(hap, hap_summary[hap]))
        total += hap_summary[hap]
    print("Total: {}".format(total))

def varfreq_mod_summary(varfreqs):
    mod_summary = {}
    for key in varfreqs:
        mod_code = key[2]
        if mod_code not in mod_summary:
            mod_summary[mod_code] = 0
        mod_summary[mod_code] += 1
    
    print("## Modification summary:")
    total=0
    for mod_code in mod_summary:
        print("{}: {}".format(mod_code, mod_summary[mod_code]))
        total += mod_summary[mod_code]
    print("Total: {}".format(total))

def varfreq_hap_summary(varfreqs):
    hap_summary = {}
    for key in varfreqs:
        hap = key[3]
        if hap not in hap_summary:
            hap_summary[hap] = 0
        hap_summary[hap] += 1
    
    print("## Haplotype summary:")
    total=0
    for hap in hap_summary:
        print("{}: {}".format(hap, hap_summary[hap]))
        total += hap_summary[hap]
    print("Total: {}".format(total))

def varfreq_vartype_summary(varfreqs):
    vartype_summary = {}
    for key in varfreqs:
        for var_pos in varfreqs[key]:
            var_type = varfreqs[key][var_pos][0]
            if var_type not in vartype_summary:
                vartype_summary[var_type] = 0
            vartype_summary[var_type] += 1
    
    print("## Variant type summary:")
    total=0
    for var_type in vartype_summary:
        print("{}: {}".format(var_type, vartype_summary[var_type]))
        total += vartype_summary[var_type]
    print("Total: {}".format(total))

def varfreq_offset_summary(varfreqs):
    offset_summary = {}
    for key in varfreqs:
        for var_pos in varfreqs[key]:
            offset = varfreqs[key][var_pos][3]
            if offset not in offset_summary:
                offset_summary[offset] = 0
            offset_summary[offset] += 1
    
    print("## Offset summary:")
    total=0
    for offset in sorted(offset_summary.keys()):
        print("{}: {}".format(offset, offset_summary[offset]))
        total += offset_summary[offset]
    print("Total: {}".format(total))

def varfreq_gt_summary(varfreqs):
    gt_summary = {}
    for key in varfreqs:
        for var_pos in varfreqs[key]:
            gt = varfreqs[key][var_pos][4]
            if gt not in gt_summary:
                gt_summary[gt] = 0
            gt_summary[gt] += 1
    
    print("## Genotype summary:")
    total=0
    for gt in sorted(gt_summary.keys()):
        print("{}: {}".format(gt, gt_summary[gt]))
        total += gt_summary[gt]
    print("Total: {}".format(total))

def cpg_gain_summary(varfreq_file):
    cats = ["SNP_start", "SNP_end",
            "DEL_start", "DEL_end",
            "INS_start", "INS_within", "INS_end",
            "MNP_start", "MNP_within", "MNP_end",
            "uncategorized"]
    summary = {c: 0 for c in cats}

    opener = gzip.open if varfreq_file.endswith(".gz") else open
    with opener(varfreq_file, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")

            contig = parts[0]
            pos = int(parts[1])
            strand = parts[5]
            ref_allele = parts[13]
            alt_allele = parts[14]
            c_off = int(parts[15])
            hap = parts[16]

            if hap == "*":
                continue
            if strand != "+":
                continue

            alt_len = len(alt_allele)
            if len(ref_allele) == 1 and alt_len == 1:
                var_type = "SNP"
            elif len(ref_allele) > alt_len:
                var_type = "DEL"
            elif len(ref_allele) < alt_len:
                var_type = "INS"
            else:
                var_type = "MNP"

            if c_off == -1:
                loc = "start"
            elif c_off == alt_len - 1:
                loc = "end"
            elif 0 <= c_off < alt_len - 1:
                loc = "within"
            else:
                loc = None

            if loc is None:
                summary["uncategorized"] += 1
                print("WARNING: uncategorized CpG gain: {}:{} {}->{} c_off={} hap={}".format(
                    contig, pos, ref_allele, alt_allele, c_off, hap), file=sys.stderr)
                continue

            bucket = "{}_{}".format(var_type, loc)
            summary.setdefault(bucket, 0) 
            summary[bucket] += 1

    total = 0
    for key in cats:
        print("{}: {}".format(key, summary[key]))
        total += summary[key]
    for key in summary:
        if key not in cats:
            print("{}: {}".format(key, summary[key]))
            total += summary[key]
    print("Total: {}".format(total))

def load_reference(fn):
    seqs = {}
    name = None
    chunks = []
    opener = gzip.open if fn.endswith(".gz") else open
    with opener(fn, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs

def reference_cpg_in_varfreq_check(varfreq_file, reference):
    leaks = []
    checked = 0
    opener = gzip.open if varfreq_file.endswith(".gz") else open
    with opener(varfreq_file, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[5] != "+" or parts[16] == "*":
                continue
            contig = parts[0]
            pos = int(parts[1])
            mod_code = parts[3]
            var_pos = int(parts[11])
            ref_allele = parts[13]
            alt_allele = parts[14]
            offset = int(parts[15])
            hap = parts[16]

            if offset != pos - var_pos:
                continue
            seq = reference.get(contig)
            if seq is None or pos + 1 >= len(seq):
                continue

            ref_len = len(ref_allele)
            alt_len = len(alt_allele)
            g_idx = offset + 2
            g_in = g_idx > ref_len and g_idx <= alt_len
            g_rp = var_pos + ref_len + (g_idx - alt_len - 1) if g_idx > alt_len else var_pos - 1 + g_idx

            checked += 1
            
            if (not g_in and g_rp == pos + 1
                    and seq[pos] in "Cc" and seq[pos + 1] in "Gg"):
                leaks.append((contig, pos, mod_code, hap, ref_allele, alt_allele))

    print("reference-mapped '+' varfreq CpGs checked: {}".format(checked))
    print("  genuine reference-CpG leaks: {}".format(len(leaks)))
    for contig, pos, mod_code, hap, ref_allele, alt_allele in leaks:
        print("  LEAK {}:{} mod={} hap={} {}->{}".format(contig, pos, mod_code, hap, ref_allele, alt_allele))

vcf_file = sys.argv[1]
freq_file = sys.argv[2]
varfreq_file = sys.argv[3]
reference_file = sys.argv[4] if len(sys.argv) == 5 else None

print("# vcf: {}".format(vcf_file))
variants = vcf_load(vcf_file)
vcf_variant_summary(variants)
vcf_gt_summary(variants)

print()

print("# freq : {}".format(freq_file))
freqs = freq_load(freq_file)
freq_mod_summary(freqs)
freq_hap_summary(freqs)

print()

print("# varfreq : {}".format(varfreq_file))
varfreqs = varfreq_load(varfreq_file)
varfreq_mod_summary(varfreqs)
varfreq_hap_summary(varfreqs)
varfreq_vartype_summary(varfreqs)
# varfreq_offset_summary(varfreqs)
varfreq_gt_summary(varfreqs)

print()

print("## CpG gain summary:")
cpg_gain_summary(varfreq_file)

print()

print("## Reference CpG leak check:")
if reference_file is None:
    print("(skipped: pass reference.fa as the 4th argument to run this check)")
else:
    reference = load_reference(reference_file)
    reference_cpg_in_varfreq_check(varfreq_file, reference)

