# minimod

A simple base modification kit. It takes a alignment BAM file and the reference FASTA as input, and outputs a base modifications (TSV) and base modification frequency (TSV or bedmethyl).

Minimod reads base modification information encoded under MM:Z and ML:B:C SAM tags and relies on the [SAMtags](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf) specification.

## Important !
Make sure that following requirements are met for each step.

### Base-calling
- Use a basecalling model trained to identify modified bases.

   Example: Base-calling a slow5 file using [buttery-eel](https://github.com/Psy-Fer/buttery-eel)
   ```
   buttery-eel --call_mods --config dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac.cfg -i reads.slow5 -o reads.sam -g path/to/guppy/bin
   samtools fastq -TMM,ML reads.sam > reads.fastq
   ```

### Aligning
- Avoid unmapped reads
- Avoid secondary alignments
- Use soft clipping for supplementary alignments

   Corresponding minimap2 flags are as follows.
   | Minimap2 Flag | Description |
   |-|-|
   |--sam-hit-only| Avoid unmapped reads.|
   |-Y | Use soft clipping for supplementary alignments.|
   |--secondary=no| Avoid secondary alignments |

   Example: aligning ONT reads using [minimap2](https://github.com/lh3/minimap2)
   ```
   minimap2 -ax map-ont --sam-hit-only -Y --secondary=no ref.idx reads.fastq
   ```

# Installation
## Building from source
```bash
git clone https://github.com/warp9seq/minimod
cd minimod
./scripts/install-hts.sh
make
```

# Usage
Usage information can be printed using ```minimod -h``` command.
```bash
Usage: minimod <command> [options]

command:
         view          view base modifications
         meth_freq     call methylation and output methylation frequency
```

## view
Print base modification details to the standard output in tsv format. Following command writes the output to a file.
```bash
basic options:
   -r FILE                    reference genome fasta file
   -m FLOAT                   mmin modification threshold (inclusive, range 0.0 to 1.0) [0.0]
   -h                         help
   -o FILE                    output to file [stdout]
   --version                  print version
```

**Example**
```bash
minimod view -r ref.fa reads.bam > meth.tsv
```

**Output fields**
| Field    | Type | Definition    |
|----------|-------------|-------------|
| 1. read_id | str | name of the read |
| 2. read_pos | int | position (0-based) of the base in read |
| 3. ref_pos   | int | position (0-based) of the base in reference |
| 4. chrom | str | chromosome |
| 5. mod_strand | char | strand (+/-) where the base modification was observed (reported by sequencer) |
| 6. read_strand | char | strand (+/-) of the read |
| 7. read_length | int | length of read |
| 8. mod_prob | float | probabiliy (0.0-1.0) of base modification |
| 9. mod_code | char | base modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf)  |
| 10. base_qual | int | basecall quality score (phred) |
| 11. ref_base | char | base as in reference  |
| 12. read_base | char | base as in read  |
| 13. mod_base | char | base as in modification tag |
| 14. flag | int | FLAG from alignment record |

**Sample meth.tsv output**
```bash
read_id	read_pos	ref_pos	chrom	mod_strand	read_strand	read_length	mod_qual	mod_code	base_qual	ref_base	read_base	mod_base	flag
m84088_230609_030819_s1/81921122/ccs	16597	20005265	chr22	+	+	16932	0.909804	m	40	c	C	C	0
m84088_230609_030819_s1/81921122/ccs	16610	20005278	chr22	+	+	16932	0.317647	m	40	c	C	C	0
m84088_230609_030819_s1/81921122/ccs	16876	20005544	chr22	+	+	16932	0.901961	m	40	c	C	C	0
m84088_230609_030819_s1/81921122/ccs	16888	20005556	chr22	+	+	16932	0.854902	m	35	c	C	C	0
m84088_230609_030819_s1/81921122/ccs	16916	20005584	chr22	+	+	16932	0.078431	m	40	c	C	C	0
m84088_230609_030819_s1/46597904/ccs	91	20001933	chr22	+	-	11990	0.494118	m	40	G	A	C	16
m84088_230609_030819_s1/46597904/ccs	123	20001901	chr22	+	-	11990	0.227451	m	40	G	A	C	16
m84088_230609_030819_s1/46597904/ccs	201	20001823	chr22	+	-	11990	0.054902	m	40	G	A	C	16
m84088_230609_030819_s1/46597904/ccs	423	20001601	chr22	+	-	11990	0.988235	m	40	g	C	C	16
m84088_230609_030819_s1/46597904/ccs	493	20001531	chr22	+	-	11990	0.003922	m	40	g	C	C	16
```

## meth_freq
Print methylation frequencies to the standard output in tsv format. Following command writes the output to a file.
```bash
basic options:
   -r FILE                    reference genome fasta file
   -b                         output in bedMethyl format
   -m FLOAT                   min modification threshold (inclusive, range 0.0 to 1.0) [0.2]
   -h                         help
   -o FILE                    output to file [stdout]
   --version                  print version
```

**Example**

```bash
minimod meth_freq -r ref.fa reads.bam > methfreq.tsv
```

**Output fields**
| Field    | Type | Definition    |
|----------|-------------|-------------|
| 1. chrom | str | choromosome |
| 2. start | int | position (0-based) of the base |
| 3. end   | int | position (0-based) of the base |
| 4. depth | int | number of reads covering the base |
| 5. n_mod | int | number of reads with probability (>0.2) for base modification |
| 6. n_called | int | number of reads called for base modification |
| 7. n_skipped | int | number of reads skipped the base as having less likelihood (modification status is '.' as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf)) |
| 8. freq | float | n_mod/n_called ratio |
| 9. mod_code | char | base modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf) |
| 10. strand | char | strand (+/-) where the base modification was observed (reported by sequencer) |
| 11. ref_base | char | base as in reference |

**Sample methfreq.tsv output**
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

**Sample methfreq.bedmethyl output**
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