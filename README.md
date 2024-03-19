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
minimod meth_freq reads.bam > meth.tsv

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
minimod meth_freq reads.bam > methfreq.tsv
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

#### Sample methfreq.tsv
```
chrom	start	end	depth	n_mod	n_called	n_skipped	freq	mod_code	strand
chr22	20026776	20026776	1	1	1	0	1.000000	m	-
chr22	20016594	20016594	2	0	2	0	0.000000	m	+
chr22	20019069	20019069	1	1	1	0	1.000000	m	+
chr22	19970705	19970705	1	0	1	0	0.000000	m	+
chr22	19981716	19981716	1	1	1	0	1.000000	m	+
chr22	20020909	20020909	3	0	3	0	0.000000	m	-
chr22	19988672	19988672	2	2	2	0	1.000000	m	-
chr22	20017060	20017060	1	0	1	0	0.000000	m	+
chr22	20016854	20016854	5	0	2	0	0.000000	m	-
```