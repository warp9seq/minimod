# Notes


## Information about the implementations

- From >=v0.5.0, we print a warning if a certain modfication and a context is untested.  For versions before this, only m[CG] for DNA (genome) and a[A] for transcriptome were tested.

## Breaking changes from v0.4.0 to v0.5.0

- Before v0.5.0, when non CG contexts were requested, both minimod view and freq were reporting modifications from errornous contexts. We have fixed and tested this issue from >= v0.5.0.

- From >=v0.5.0, we compare if modified read base matches aligned reference base and output only the matching bases by default. However, this comparison can be avoided using --include-alt-alleles option. (Note that only the modified base is compared, not the whole context.)

- Before v0.4.0, we allowed primary, secondary, suppelmentary alignments when viewing and calculating frequencies. From >=v0.5.0, we only consider primary, supplementary alignments by default. We introduced --secondary option to enable considering secondary alignments. Minimod still errors out if hard-clipping is found.

- Before v0.4.0, when --insertions option is used, errornous reference positions were included in both freq and view outputs. We have fixed it in >=v0.5.0.

## Modkit consistency
Tool versions we used for comparisons are modkit 0.5.1 and minimod 0.5.0

- Comparing modified base with reference base
    - minimod by default outputs modified bases that match with reference base.
    - modkit output modified bases without matching with reference base when the context is not specified using --cpg or --motif options.
    - to match minimod's behaviour with modkit's, use --include-alt-alleles option with minimod. 
- Which reads are used for computations
    - By default modkit extract ignores non-primary alignments and --allow-non-primary option can allow secondary and supplementary alignments if a valid MN tag is found (https://github.com/nanoporetech/modkit/blob/481e3c9e7930f3f499eadf1ef441606f33e6881c/book/src/intro_extract.md#note-on-non-primary-alignments). 
    - minimod ignores secondary alignments by default and uses primary and supplementary alignments. Using --secondary options can allow them secondary alignments. minimod does not require the MN tag to allow non-primary alignments. Further, minimod errors out when hard-clipping is detected.
