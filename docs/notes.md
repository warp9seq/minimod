# Notes


## Information about the implementations

- From >=v0.5.0, we print a warning if a certain modfication and a context is untested.  For versions before this, only m[CG] for DNA (genome) and a[A] for transcriptome were tested.

## Breaking changes from v0.4.0 to v0.5.0

- Before v0.5.0, when non CG contexts were requested, both minimod view and freq were reporting modifications from errornous contexts. We have fixed and tested this issue from >= v0.5.0.

- From >=v0.5.0, when a context is provided (ex: -c h[CG]) option, we output modifications where the modified read base matches aligned reference base. However, when the context is * (ex: -c a[*]), this comparison is ignored and modifications at both matched and mismatched positions are output. Note that only the modified base is compared with reference base, not the whole context.

- Before v0.4.0, we allowed primary, secondary, suppelmentary alignments when viewing and calculating frequencies. From >=v0.5.0, we only consider primary, supplementary alignments by default. We introduced --secondary option to enable considering secondary alignments. Minimod still errors out if hard-clipping is found.

- Before v0.4.0, when --insertions option is used, errornous reference positions were included in both freq and view outputs. We have fixed it in >=v0.5.0.

## Modkit consistency
Tool versions we used for comparisons are modkit 0.5.1 and minimod 0.5.0

- modkit by default outputs all modification without matching the read base with reference base when the context is not specified using --cpg or --motif options.
- minimod by default outputs m modified bases in CG context in reference that are mapped and matched with reference base.
- modkit by default extract ignores non-primary alignments and --allow-non-primary option can allow secondary and supplementary alignments if a valid MN tag is found (https://github.com/nanoporetech/modkit/blob/481e3c9e7930f3f499eadf1ef441606f33e6881c/book/src/intro_extract.md#note-on-non-primary-alignments).
- minimod by default ignores secondary alignments and uses primary and supplementary alignments. Using --secondary options can allow them secondary alignments. minimod does not require the MN tag to allow non-primary alignments. Further, minimod errors out when hard-clipping is detected.


Following pair of commands using minimod v0.5.0 and modkit 0.5.1 should give the same output. (To compare modkit's extract bed with minimod's view tsv, test/compare_view_mkbed_mmtsv.sh script can be used)

Let's assume reads.bam containes mapped and unmapped, primary, secondary and supplementary alignments.

When the context and modification type is unknown
```bash
modkit extract full --mapped-only reads.bam mk_extract.bed
minimod view -c '*' --skip-supplementary  ref.fa reads.bam > mm_view.tsv
test/compare_view_mkbed_mmtsv.sh mk_extract.bed mm_view.tsv out_dir
```

When the context is known and modification type is unknown
```bash
modkit extract full --motif A 0 --mapped-only --reference ref.fa reads.bam mk_extract_A.bed
minimod view -c '*[A]' --skip-supplementary ref.fa reads.bam > mm_view_A.tsv
test/compare_view_mkbed_mmtsv.sh mk_extract_A.bed mm_view.tsv out_dir
```

When the context and modification type is known
```bash
modkit extract full --motif A 0 --mapped-only --reference ref.fa reads.bam mk_extract_A.bed
awk 'NR==1 || $14=="a"' mk_extract_A.bed > mk_extract_aA.bed
minimod view -c 'a[A]' --skip-supplementary ref.fa reads.bam > mm_view_aA.tsv
test/compare_view_mkbed_mmtsv.sh mk_extract_a.bed mm_view_aA.tsv out_dir
```

## Philosophy

Minimod is intended to be kept simple. Its main use case is for reference-based analysis.

## Ratinale for decisions

Why reference is still take for a[*]?
Sanity check if the BAM header contigs match to the reference

## modkit v0.5.0 equivalence commands