# minimod

A simple base modification kit.

## Installation
### Building from source
```bash
git clone https://github.com/warp9seq/minimod
cd minimod
./scripts/install-hts.sh
make
```

## Usage
Usage information can be printed using ```minimod -h``` command.
```bash
Usage: minimod <command> [options]

command:
         view          view base modifications
         meth_freq     call methylation and output methylation frequency
```

## Subtools
### view
Print base modification details to the standard output in tsv format. Following command writes the output to a file.
```bash
basic options:
   -r FILE                    reference genome fasta file
   -h                         help
   -o FILE                    output to file [stdout]
   --version                  print version
```

Example
```bash
minimod view -r ref.fa reads.bam > meth.tsv
```
<!-- | Field    | Definition    |
|----------|-------------|
| 1. read_id | name of the read |
| 2. read_pos | 0 based position of the base on read sequence |
| 3. ref_pos   | 0 based position of the base on reference sequence |
| 4. chrom | chromosome name |
| 5. mod_strand |  |
| 6. ref_strand |  |
| 7. ref_mod_strand |  |
| 8. fw_soft_clipped_start |  |
| 9. fw_soft_clipped_end |  |
| 10. read_length |  |
| 11. mod_qual |  |
| 12. mod_code |  |
| 13. base_qual |  |
| 14. ref_kmer |  |
| 15. query_kmer |  |
| 16. canonical_base |  |
| 17. modified_primary_base |  |
| 18. inferred |  |
| 19. flag |  | -->

### meth_freq
Print methylation frequencies to the standard output in tsv format. Following command writes the output to a file.
```bash
basic options:
   -r FILE                    reference genome fasta file
   -h                         help
   -b                         output in bedMethyl format
   -o FILE                    output to file [stdout]
   --version                  print version
```

Example

```bash
minimod meth_freq -r ref.fa reads.bam > methfreq.tsv
```

Output fields
| Field    | Definition    |
|----------|-------------|
| 1. chrom | choromosome |
| 2. start | position (0-based) of the base |
| 3. end   | position (0-based) of the base |
| 4. depth | number of reads covering the base |
| 5. n_mod | number of reads with probability (>0.2) for base modification |
| 6. n_called | number of reads called for base modification |
| 7. n_skipped | number of reads skipped the base as having less likelihood for modification |
| 8. freq | n_mod/n_called ratio |
| 9. mod_code | modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf) |
| 10. strand | strand (+/-) where the base modification is observed |
| 11. ref_base | base charater as in reference |

#### Sample methfreq.tsv
```
chrom	start	end	depth	n_mod	n_called	n_skipped	freq	mod_code	strand	ref_base
chr22	20016337	20016337	5	0	5	0	0.000000	m	-	G
chr22	20016594	20016594	2	0	2	0	0.000000	m	+	C
chr22	20017045	20017045	1	0	1	0	0.000000	m	+	C
chr22	19970705	19970705	1	0	1	0	0.000000	m	+	C
chr22	19981716	19981716	1	1	1	0	1.000000	m	+	C
chr22	20020909	20020909	3	0	3	0	0.000000	m	-	G
chr22	19995719	19995719	4	2	4	0	0.500000	m	+	C
chr22	20017060	20017060	1	0	1	0	0.000000	m	+	C
chr22	19971259	19971259	1	1	1	0	1.000000	m	+	C
```

#### Sample methfreq.bedmethyl
```
chr22	20016337	20016337	m	5	-	20016337	20016337	255,0,0	5	0.000000
chr22	20016594	20016594	m	2	+	20016594	20016594	255,0,0	2	0.000000
chr22	20017045	20017045	m	1	+	20017045	20017045	255,0,0	1	0.000000
chr22	19970705	19970705	m	1	+	19970705	19970705	255,0,0	1	0.000000
chr22	19981716	19981716	m	1	+	19981716	19981716	255,0,0	1	1.000000
chr22	20020909	20020909	m	3	-	20020909	20020909	255,0,0	3	0.000000
chr22	19995719	19995719	m	4	+	19995719	19995719	255,0,0	4	0.500000
chr22	20017060	20017060	m	1	+	20017060	20017060	255,0,0	1	0.000000
chr22	19971259	19971259	m	1	+	19971259	19971259	255,0,0	1	1.000000
```