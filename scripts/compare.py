#! /usr/bin/env python
#! pip install numpy
#---------------------------------------------------------
# Copyright 2015 Ontario Institute for Cancer Research
# Written by Jared Simpson (jared.simpson@oicr.on.ca)
#---------------------------------------------------------
# Obtained from https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html
#---------------------------------------------------------
# Updated by Suneth Samarasinghe (suneth@unsw.edu.au) in 2024
#---------------------------------------------------------

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
        num_reads = float(fields[9])
        percent_methylated = float(fields[10])
        methylated_reads = int( (percent_methylated / 100) * num_reads)
        key = ""

        # accumulate on forward strand
        if strand == "+":
            key = make_key(chromosome, str(start), str(start))
        else:
            key = make_key(chromosome, str(start - 1), str(start - 1))

        update_stats(out, key, num_reads, methylated_reads, "bedmethyl")

    return out

# Load the file of methylation frequency based on the filename
def load_methylation(filename):
    if filename.find("bedmethyl") != -1:
        return load_bedmethyl(filename)
    elif filename.find("tsv") != -1:
        return load_tsv(filename)
    else:
        sys.stderr.write("ERROR: unknown methylation file format. Suppprted ones are .tsv and .bedmethyl\n" % filename)
        sys.exit(1)

# Pearson correlation coefficient
def pearson_correlation(x, y):
    n = len(x)
    if n == 0:
        return 0
    sum_x = float(sum(x))
    sum_y = float(sum(y))
    sum_x_sq = sum(xi*xi for xi in x)
    sum_y_sq = sum(yi*yi for yi in y)
    psum = sum(xi*yi for xi, yi in zip(x, y))
    num = psum - (sum_x * sum_y/n)
    den = ((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n)) ** 0.5
    if den == 0: return 0
    return num / den

set1 = load_methylation(sys.argv[1])
set2 = load_methylation(sys.argv[2])

list1 = []
list2 = []

for key in set1:
    if key in set2:
        list1.append(set1[key].methylation_frequency())
        list2.append(set2[key].methylation_frequency())

print(pearson_correlation(list1, list2))