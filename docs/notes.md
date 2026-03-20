# Notes

## Information about the implementations

- From >=v0.5.0, we print a warning if a certain modification and a context is untested.  For versions before this, only m[CG] for DNA (genome) and a[A] for transcriptome were tested.

## Major changes from v0.4.0 to v0.5.0

- For <=v0.4.0, when non CG contexts were requested, both minimod view and freq were reporting modifications from erroneous contexts. We have fixed and tested this issue from >= v0.5.0.

- For <=v0.4.0, we allowed primary, secondary, supplementary alignments when viewing and calculating frequencies. From >=v0.5.0, we only consider primary and supplementary alignments by default. We introduced --allow-secondary option to enable considering secondary alignments. Minimod still errors out if hard-clipping is found.

- For <=v0.4.0, when --insertions option is used, erroneous reference positions were included in both freq and view outputs. We have fixed it in >=v0.5.0.

- For <=v0.4.0, status of skipped bases (low probability C+m. or unknown C+m? in MM tag) were not handled. We have fixed it in >=v0.5.0. When skipped as low probability, minimod consider them to have 0 modification probability. Otherwise, when skipped as unknown, they are ignored.

- For >=v0.5.0, mod_prob in the output of view subtool will slightly change due to a change on how probability conversion from unsigned 8-bit N to float probability p. Previously, it was p = N/255.0 and now we use p = (N+0.5)/256.0 as the conversion function.

- For >=v0.5.0, when a context is provided (ex: -c h[CG]) option, we output modifications where the modified read base matches aligned reference base. However, when the context is * (ex: -c a[*]), this comparison is ignored and modifications at both matched and mismatched positions are output. Note that only the modified base is compared with reference base, not the whole context.


## Modkit consistency
Tool versions we used for comparisons are modkit 0.5.1 and minimod 0.5.0

- Minimod by default outputs m modified bases in CG context in reference that are mapped and matched with reference base. 

  Modkit by default outputs all modification without matching the read base with reference base when the context is not specified using --cpg or --motif options.

- Minimod by default ignores secondary alignments and uses primary and supplementary alignments. Using --allow-secondary options can allow secondary alignments. minimod does not require the MN tag to allow non-primary alignments. Further, minimod errors out when hard-clipping is detected. 

  Modkit by default extract ignores non-primary alignments and --allow-non-primary option can allow secondary and supplementary alignments if a valid MN tag is found (https://github.com/nanoporetech/modkit/blob/481e3c9e7930f3f499eadf1ef441606f33e6881c/book/src/intro_extract.md#note-on-non-primary-alignments).


### Minimod view vs Modkit extract full

Following pair of commands using minimod v0.5.0 and modkit 0.5.1 should give the same output. To compare modkit's extract bed with minimod's view tsv, test/compare_view_mkbed_mmtsv.sh script can be used. Let's assume reads.bam contains mapped and unmapped, primary, secondary and supplementary alignments.

When the context and modification type is unknown
```bash
minimod view -c '*' --skip-supplementary  ref.fa reads.bam > mm_view.tsv

modkit extract full --mapped-only reads.bam mk_extract.bed

test/compare_view_mkbed_mmtsv.sh mk_extract.bed mm_view.tsv out_dir
```

When the context is known and modification type is unknown
```bash
minimod view -c '*[A]' --skip-supplementary ref.fa reads.bam > mm_view_A.tsv

modkit extract full --motif A 0 --mapped-only --reference ref.fa reads.bam mk_extract_A.bed

test/compare_view_mkbed_mmtsv.sh mk_extract_A.bed mm_view_A.tsv out_dir
```

```bash
# alternative way to filter context from modkit output using awk instead of using --motif A 0 
modkit extract full --mapped-only --kmer-size 1 --reference ref.fa reads.bam mk_extract.bed
awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || $16=="A" && ($21==16 ? toupper($16)==c($17) : toupper($16)==toupper($17))' mk_extract.bed > mk_extract_A.bed
```

When the context and modification type is known
```bash
minimod view -c 'a[A]' --skip-supplementary ref.fa reads.bam > mm_view_aA.tsv

modkit extract full --motif A 0 --mapped-only --reference ref.fa reads.bam mk_extract_A.bed
awk 'NR==1 || $14=="a"' mk_extract_A.bed > mk_extract_aA.bed

test/compare_view_mkbed_mmtsv.sh mk_extract_aA.bed mm_view_aA.tsv out_dir
```

```bash
# alternative way to filter context from modkit output using awk instead of using --motif A 0 
modkit extract full --mapped-only --kmer-size 1 --reference ref.fa reads.bam mk_extract.bed
awk 'NR==1 || $14=="a"' mk_extract.bed > mk_extract_a.bed
awk -F '\t' 'function c(b){b=toupper(b);return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b} NR==1 || $16=="A" && ($21==16 ? toupper($16)==c($17) : toupper($16)==toupper($17))' mk_extract_a.bed > mk_extract_aA.bed
```

### Minimod freq vs Modkit pileup

While the default behaviours mentioned [above](#modkit-consistency) are applied to minimod freq and modkit pileup tools, following default behaviours can cause decrepencies between the outputs from two tools.

- Minimod freq by default uses 0.8 as the modification threshold (which used to determine if a base if modified).

  Modkit pileup by default computes the modification threshold programatically observing the data unless the threshold is provided by the user.

- Minimod freq uses a modification filtering method as explained in [here ](../README.md#modification-threshold) in README.md.

  Modkit pileup uses two-way, three-way base modification calls as explained [here](https://github.com/nanoporetech/modkit/blob/v0.5.1-rc1/book/src/filtering_details.md)

Following pair of commands using minimod v0.5.0 and modkit 0.5.1 should give similar outputs. To compute Pearson correlation coefficient between outputs from modkit's pileup and minimod's freq, test/compare.py Python script can be used. Let's assume reads.bam contains mapped and unmapped, primary, secondary and supplementary alignments.

When the context and modification type is unknown
```bash
minimod freq -b -c '*' --skip-supplementary ref.fa reads.bam > mm_freq.bed

modkit pileup --reference ref.fa reads.bam mk_pileup.bed

test/compare.py mk_pileup.bed mm_freq.bed
```

When the context is known and modification type is unknown
```bash
minimod freq -b -c '*[A]' --skip-supplementary ref.fa reads.bam > mm_freq_A.bed

modkit pileup --motif A 0 --reference ref.fa reads.bam mk_pileup_A.bed

test/compare.py mk_pileup_A.bed mm_freq_A.bed
```

When the context and modification type is known
```bash
minimod freq -b -c 'a[A]' --skip-supplementary ref.fa reads.bam > mm_freq_aA.bed

modkit pileup --motif A 0 --reference ref.fa reads.bam mk_pileup_A.bed
awk 'NR==1 || $4=="a"' mk_pileup_A.bed > mk_pileup_aA.bed # not needed if mod+basecall model is for only type a modifications (ex: m6A_DRACH)

test/compare.py mk_pileup_A.bed mm_freq_aA.bed
```

When mod+basecalled using a 2-way classification model such as m6A_DRACH (classifies into modified as m6A or canonical A), minimod's and modkit's filtering methods results in same outputs if the modification threshold is matched using -m in minimod or --filter-threshold in modkit.
```bash
minimod freq --skip-supplementary -m 0.9 -b -c "a[A]" ref.fa rna_m6A_DRACH.bam > mm_freq_aA.bed

modkit pileup --filter-threshold A:0.9 --motif A 0 --reference ref.fa rna_m6A_DRACH.bam mk_pileup_aA.bed

test/compare_freq_bed_bed.sh mm_freq_aA.bed mk_pileup_aA.bed out_dir
```

## Philosophy

Minimod is intended to be kept simple. Its main use case is for reference-based analysis.

## Rationale for decisions

- Why reference is still taken for a[*] ?
  - Sanity check if the BAM header contigs match to the reference
- Why unsigned 8-bit probability, N in MM tag converted to float probability, p using p = (N+0.5)/256 formula ?
  - According to SAMtags specification, N represents a probability range N/256 to (N+1)/256. Therefore, It is reasonable to choose the mid of that range when converting back to float probability.
