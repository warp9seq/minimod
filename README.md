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
         mod-freq      output base modifications frequencies
```

# Examples
```bash
# modification of type m in tsv format
minimod view ref.fa reads.bam > mods.tsv

# modification frequencies of type m in tsv format (default threshold 0.2)
minimod mod-freq ref.fa reads.bam > modfreqs.tsv

# modification frequencies of type m (5-methylcytosine) in bedmethyl format (default threshold 0.2)
minimod mod-freq -b ref.fa reads.bam > modfreqs.bedmethyl

# modification frequencies  of types m (5-methylcytosine) and h (5-hydroxymethylcytosine) with thresholds 0.2 and 0.3 respectively in tsv format
minimod mod-freq -c mh -m 0.2,0.3 ref.fa reads.bam > mods.tsv
```
Click here to see [how threshold is used in minimod](#modification-codes-and-threshold)

# minimod view
```bash
minimod view ref.fa reads.bam > mods.tsv
```
This writes all base modifications (default modification code "m" and modification threshold 0.2) to a file (mods.tsv) in tsv format. Sample output is given below.
```
Usage: minimod view ref.fa reads.bam

basic options:
   -c STR                     modification code(s) (ex. m , h or mh) [m]
   -t INT                     number of processing threads [8]
   -K INT                     batch size (max number of reads loaded at once) [512]
   -B FLOAT[K/M/G]            max number of bytes loaded at once [20.0M]
   -h                         help
   -p INT                     print progress every INT seconds (0: per batch) [0]
   -o FILE                    output file [stdout]
   --verbose INT              verbosity level [4]
   --version                  print version
```

**Sample mods.tsv output**
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
| 3. strand | char | strand (+/-) of the read |
| 4. read_id | str | name of the read |
| 5. read_pos | int | position (0-based) of the base in read |
| 6. mod_code | char | base modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf)  |
| 7. mod_prob | float | probability (0.0-1.0) of base modification |

# minimod mod-freq
```bash
minimod mod-freq ref.fa reads.bam > modfreqs.tsv
```
This writes base modification frequencies (default modification code "m" and modification threshold 0.2) to a file (modfreqs.tsv) file in tsv format. Sample output is given below.
```
Usage: minimod mod-freq ref.fa reads.bam

basic options:
   -b                         output in bedMethyl format [not set]
   -c STR                     modification codes (ex. m , h or mh) [m]
   -m FLOAT                   min modification threshold(s). Comma separated values for each modification code given in -c [0.2]
   -t INT                     number of processing threads [8]
   -K INT                     batch size (max number of reads loaded at once) [512]
   -B FLOAT[K/M/G]            max number of bytes loaded at once [20.0M]
   -h                         help
   -p INT                     print progress every INT seconds (0: per batch) [0]
   -o FILE                    output file [stdout]
   --verbose INT              verbosity level [4]
   --version                  print version
```

**Sample modfreqs.tsv output**
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
| 4. strand | char | strand (+/-) of the read |
| 5. n_called | int | number of reads called for base modification |
| 6. n_mod | int | number of reads with probability (>0.2) for base modification |
| 7. freq | float | n_mod/n_called ratio |
| 8. mod_code | char | base modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf) |

**Sample modfreqs.bedmethyl output**
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
| 6. strand | char | strand (+/-) of the read |
| 7. start | int | = field 2 |
| 8. end   | int | = field 3 |
| 9. n_mod | int | = field 5 |
| 10. freq | float | n_mod/n_called ratio |

# Modification codes and threshold

Base modification codes and thresholds can be set for mod-freq tools using -c and -m flags respectively.


01. 5mC modification frequencies with threshold 0.2
> ```
>   minimod mod-freq ref.fa reads.bam -c m -m 0.2
> ```
![Fig](docs/figs/m_threshold.png)
> ```
> If p(5mC) >=  0.8 (threshold),       called(5mC) and modified(5mC)
> If p(5mC) <=  0.2 (1-threshold),     called(5mC)
> ```

03. 5mC and 5hmC base modification frequencies with thresholds 0.2, 0.5 respectively
> ```
>   minimod mod-freq ref.fa reads.bam -c mh -m 0.2,0.3
> ```
![Fig](docs/figs/m_h_threshold.png)
> ```
> If p(5mC)  >=  0.2 (threshold),      called(5mC) and modified(5mC)
> If p(5mC)  <=  0.8 (1-threshold),    called(5mC) 
>
> If p(5hmC) >=  0.3 (threshold),      called(5hmC) and modified(5hmC)
> If p(5hmC) <=  0.7 (1-threshold),    called(5hmC)
> ```

# Important !
Make sure that following requirements are met for each step in base modification calling pipeline.

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