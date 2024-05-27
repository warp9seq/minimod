# minimod

A simple base modification kit. It takes a alignment BAM file and the reference FASTA as input, and outputs a base modifications (TSV) and base modification frequency (TSV or bedmethyl).

Minimod reads base modification information encoded under MM:Z and ML:B:C SAM tags and relies on the [SAMtags](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf) specification.

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
         meth-freq     call methylation and output methylation frequency
```

# Examples
```bash
# base modification details tsv format
minimod view -r ref.fa reads.bam > meth.tsv

# methylation frequencies tsv format.
minimod meth-freq -r ref.fa reads.bam > methfreq.tsv

# methylation frequencies bed format.
minimod meth-freq -r ref.fa -b reads.bam > methfreq.tsv
```

# minimod view
```bash
minimod view -r ref.fa reads.bam > meth.tsv
```
Print base modification details to the standard output in tsv format. Following command writes the output to a file.
```bash
basic options:
   -r FILE                    reference genome fasta file
   -m FLOAT                   min modification threshold (inclusive, range 0.0 to 1.0) [0.0]
   -h                         help
   -o FILE                    output to file [stdout]
   --version                  print version
```

**Sample meth.tsv output**
```bash
ref_contig	ref_pos	strand	read_id	read_pos	mod_code	mod_prob
chr22	19979864	+	m84088_230609_030819_s1/55512555/ccs	14	m	0.709804
chr22	19979882	+	m84088_230609_030819_s1/55512555/ccs	32	m	0.949020
chr22	19979885	+	m84088_230609_030819_s1/55512555/ccs	35	m	0.980392
chr22	19979888	+	m84088_230609_030819_s1/55512555/ccs	38	m	0.780392
chr22	19979900	+	m84088_230609_030819_s1/55512555/ccs	50	m	0.623529
chr22	19979902	+	m84088_230609_030819_s1/55512555/ccs	52	m	0.992157
chr22	19979929	+	m84088_230609_030819_s1/55512555/ccs	79	m	0.941176
chr22	19979939	+	m84088_230609_030819_s1/55512555/ccs	89	m	0.141176
chr22	19979948	+	m84088_230609_030819_s1/55512555/ccs	98	m	0.623529
```

| Field    | Type | Definition    |
|----------|-------------|-------------|
| 1. ref_contig | str | chromosome |
| 2. ref_pos   | int | position (0-based) of the base in reference |
| 3. strand | char | strand (+/-) where the base modification was observed (reported by sequencer) |
| 4. read_id | str | name of the read |
| 5. read_pos | int | position (0-based) of the base in read |
| 6. mod_code | char | base modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf)  |
| 7. mod_prob | float | probabiliy (0.0-1.0) of base modification |

# minimod meth-freq
```bash
minimod meth-freq -r ref.fa reads.bam > methfreq.tsv
```
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

**Sample methfreq.tsv output**
```
contig	start	end	strand	n_called	n_mod	freq	mod_code
chr22	20016337	20016337	+	5	0	0.000000	m
chr22	20016594	20016594	+	2	0	0.000000	m
chr22	20017045	20017045	+	1	0	0.000000	m
chr22	19970705	19970705	+	1	0	0.000000	m
chr22	19981716	19981716	+	1	1	1.000000	m
chr22	20020909	20020909	+	3	0	0.000000	m
chr22	19995719	19995719	+	4	2	0.500000	m
chr22	20017060	20017060	+	1	0	0.000000	m
chr22	19971259	19971259	+	1	1	1.000000	m
```
| Field    | Type | Definition    |
|----------|-------------|-------------|
| 1. contig | str | choromosome |
| 2. start | int | position (0-based) of the base |
| 3. end   | int | position (0-based) of the base |
| 4. strand | char | strand (+/-) where the base modification was observed (reported by sequencer) |
| 5. n_called | int | number of reads called for base modification |
| 6. n_mod | int | number of reads with probability (>0.2) for base modification |
| 7. freq | float | n_mod/n_called ratio |
| 8. mod_code | char | base modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf) |

**Sample methfreq.bedmethyl output**
```
chr22	20016337	20016338	m	5	+	20016337	20016337	255,0,0	5	0.000000
chr22	20016594	20016595	m	2	+	20016594	20016594	255,0,0	2	0.000000
chr22	20017045	20017046	m	1	+	20017045	20017045	255,0,0	1	0.000000
chr22	19970705	19970706	m	1	+	19970705	19970705	255,0,0	1	0.000000
chr22	19981716	19981717	m	1	+	19981716	19981716	255,0,0	1	1.000000
chr22	20020909	20020910	m	3	+	20020909	20020909	255,0,0	3	0.000000
chr22	19995719	19995720	m	4	+	19995719	19995719	255,0,0	4	0.500000
chr22	20017060	20017061	m	1	+	20017060	20017060	255,0,0	1	0.000000
chr22	19971259	19971260	m	1	+	19971259	19971259	255,0,0	1	1.000000
chr22	19973437	19973438	m	1	+	19973437	19973437	255,0,0	1	1.000000
```
| Field    | Type | Definition    |
|----------|-------------|-------------|
| 1. contig | str | choromosome |
| 2. start | int | position (0-based) of the base |
| 3. end   | int | position (0-based) of the base |
| 4. mod_code | char | base modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf) |
| 5. n_mod | int | number of reads with probability (>0.2) for base modification |
| 6. strand | char | strand (+/-) where the base modification was observed (reported by sequencer) |
| 7. start | int | = field 2 |
| 8. end   | int | = field 3 |
| 9. n_mod | int | = field 5 |
| 10. freq | float | n_mod/n_called ratio |

# Important !
Make sure that following requirements are met for each step.

## Base-calling
- Use a basecalling model trained to identify modified bases.

   Example: Base-calling a slow5 file using [buttery-eel](https://github.com/Psy-Fer/buttery-eel)
   ```
   buttery-eel --call_mods --config dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac.cfg -i reads.slow5 -o reads.sam -g path/to/guppy/bin
   samtools fastq -TMM,ML reads.sam > reads.fastq
   ```

## Aligning
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