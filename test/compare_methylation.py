#! /usr/bin/env python

#---------------------------------------------------------
# Copyright 2015 Ontario Institute for Cancer Research
# Written by Jared Simpson (jared.simpson@oicr.on.ca)
#---------------------------------------------------------
# Obtained from https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html

import sys
import csv
import argparse

def make_key(c, s, e):
    return c + ":" + str(s) + "-" + str(e)

def parse_key(key):
    return key.replace("-", ":").split(":")

class MethylationStats:
    def __init__(self, num_reads, num_methylated, atype):
        self.num_reads = num_reads
        self.num_methylated_reads = num_methylated
        self.analysis_type = atype

    def accumulate(self, num_reads, num_methylated):
        self.num_reads += num_reads
        self.num_methylated_reads += num_methylated

    def methylation_frequency(self):
        return float(self.num_methylated_reads) / self.num_reads

def update_stats(collection, key, num_reads, num_methylated_reads, atype):
    if key not in collection:
        collection[key] = MethylationStats(num_reads, num_methylated_reads, atype)
    else:
        collection[key].accumulate(num_reads, num_methylated_reads)

def load_mmtsv(filename):
    out = dict()
    csv_reader = csv.DictReader(open(filename), delimiter='\t')

    for record in csv_reader:
        chrom = record["contig"]
        start = int(record["start"])
        strand = record["strand"]
        n_mod = int(record["n_mod"])
        n_called = int(record["n_called"])
        # n_skipped = int(record["n_skipped"])

        # accumulate on forward strand
        if strand == "+":
            key = make_key(chrom, str(start), str(start))
        else:
            key = make_key(chrom, str(start - 1), str(start - 1))

        update_stats(out, key, n_called, n_mod, "mm.tsv")

    return out

def load_tsv(filename):
    out = dict()
    csv_reader = csv.DictReader(open(filename), delimiter='\t')

    for record in csv_reader:
        key = make_key(record["chromosome"], record["start"], record["end"])

        # skip non-singleton, for now
        if int(record["num_motifs_in_group"]) > 1:
            continue

        num_reads = int(record["called_sites"])
        methylated_reads = int(record["called_sites_methylated"])
        out[key] = MethylationStats(num_reads, methylated_reads, "tsv")

    return out

def load_bedmethyl(filename):
    out = dict()
    fh = open(filename)
    for line in fh:
        fields = line.rstrip().split()
        chromosome = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[5]
        num_reads = int(fields[9])
        percent_methylated = float(fields[10])
        # The bedMethyl percentage is a rounded print of methylated/num_reads, so
        # the product lands just below the true integer (e.g. 33.333333% of 3 is
        # 0.99999999). Truncating with int() would silently drop that read.
        methylated_reads = round((percent_methylated / 100) * num_reads)
        key = ""

        # accumulate on forward strand
        if strand == "+":
            key = make_key(chromosome, str(start), str(start))
        else:
            key = make_key(chromosome, str(start - 1), str(start - 1))

        update_stats(out, key, num_reads, methylated_reads, "bedmethyl")

    return out

def load_bed(filename):
    out = dict()
    fh = open(filename)
    for line in fh:
        fields = line.rstrip().split()
        chromosome = fields[0]
        start = int(fields[1])
        strand = fields[5]
        num_reads = 1
        methylated_reads = 1
        key = ""

        # accumulate on forward strand
        if strand == "+":
            key = make_key(chromosome, str(start), str(start))
        else:
            key = make_key(chromosome, str(start - 1), str(start - 1))

        update_stats(out, key, num_reads, methylated_reads, "bed")

    return out

def load_pb_cpg_tool_bed(filename):
    out = dict()
    fh = open(filename)
    for line in fh:
        fields = line.rstrip().split()
        chromosome = fields[0]
        start = int(fields[1])
        modified_count = float(fields[6])
        unmodified_count = float(fields[7])
        num_reads = modified_count + unmodified_count

        key = make_key(chromosome, str(start), str(start))

        out[key] = MethylationStats(num_reads, modified_count, "pb-cpg-tool-bed")

    return out

def load_rmbase_bed(filename):
    out = dict()
    fh = open(filename)
    for line in fh:
        fields = line.rstrip().split()
        chromosome = fields[0]
        start = int(fields[1])
        
        strand = fields[5]
        modified_count = int(fields[7])
        num_reads = 4
        
        key = ""
        # accumulate on forward strand
        if strand == "+":
            key = make_key(chromosome, str(start), str(start))
        else:
            key = make_key(chromosome, str(start - 1), str(start - 1))

        out[key] = MethylationStats(num_reads, modified_count, "rmbase-bed")

    return out

def load_cov(filename):
    out = dict()
    fh = open(filename)
    for line in fh:
        fields = line.rstrip().split()
        chromosome = fields[0]
        start = int(fields[1])
        # end = int(fields[2])
        # percent_methylated = float(fields[3])
        methylated_reads = int(fields[4])
        unmethylated_reads = int(fields[5])
        key = make_key(chromosome, str(start), str(start))
        num_reads = methylated_reads + unmethylated_reads
        out[key] = MethylationStats(num_reads, methylated_reads, "cov")

    return out

# Load the file of methylation frequency based on the filename
def load_methylation(filename):
    if filename.find("bedmethyl") != -1: #nanopolish bedmethyl output
        return load_bedmethyl(filename)
    if filename.find("pb.bed") != -1: #pg-cpg-tool output bed file
        return load_pb_cpg_tool_bed(filename)
    elif filename.find("rmbase.bed") != -1: # data from https://rna.sysu.edu.cn/rmbase3
        return load_rmbase_bed(filename)
    if filename.find("bed") != -1:
        return load_bed(filename)
    elif filename.find("mm.tsv") != -1: #minimod meth_freq output
        return load_mmtsv(filename)
    elif filename.find("tsv") != -1:
        return load_tsv(filename)
    elif filename.find("cov") != -1:
        return load_cov(filename)
    else:
        sys.stderr.write("ERROR: unknown methylation file format. Suppprted ones are .tsv and .bedmethyl\n" % filename)
        sys.exit(1)

set1 = load_methylation(sys.argv[1])
set2 = load_methylation(sys.argv[2])

output = 0
print("key\tdepth_1\tfrequency_1\tdepth_2\tfrequency_2")
for key in set1:
    if key in set2:

        d1 = set1[key]
        d2 = set2[key]

        if d1.num_reads == 0 or d2.num_reads == 0:
            continue
        print("%s\t%d\t%.4f\t%d\t%.4f" % (key, d1.num_reads, d1.methylation_frequency(), d2.num_reads, d2.methylation_frequency()))
        output += 1

sys.stderr.write("set1 sites: %d set2 sites: %d output: %d\n" % (len(set1), len(set2), output))