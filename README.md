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
# view all modifications in tsv format (default mod code: m, context:CG)
minimod view ref.fa reads.bam > mods.tsv

# modification frequencies in tsv format (default mod code: m and threshold: 0.8 context:CG)
minimod mod-freq ref.fa reads.bam > modfreqs.tsv

# modification frequencies in bedmethyl format (default mod code: m and threshold: 0.8 context:CG)
minimod mod-freq -b ref.fa reads.bam > modfreqs.bedmethyl

# modification frequencies of multiple types ( m (5-methylcytosine) and h (5-hydroxymethylcytosine) in CG context with thresholds 0.8 and 0.7 respectively )
minimod mod-freq -c m[CG],h[CG] -m 0.8,0.7 ref.fa reads.bam > mods.tsv
```
- See [how modification codes can be specified?](#modification-codes)
- See [how threshold is used in minimod?](#modification-threshold)

# minimod view
```bash
minimod view ref.fa reads.bam > mods.tsv
```
This writes all base modifications (default modification code "m") to a file (mods.tsv) in tsv format. Sample output is given below.
```bash
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
   --insertions               enable modifications in insertions
```

- See [how to consider inserted modified bases?](#modified-bases-in-insertions)

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
| 8. ins_offset | int | offset of inserted base from ref_pos (only output when --insertions is specified) |

# minimod mod-freq
```bash
minimod mod-freq ref.fa reads.bam > modfreqs.tsv
```
This writes base modification frequencies (default modification code "m" in CG context with modification threshold 0.8) to a file (modfreqs.tsv) file in tsv format.
```bash
Usage: minimod mod-freq ref.fa reads.bam

basic options:
   -b                         output in bedMethyl format [not set]
   -c STR                     modification codes (ex. m , h or mh) [m]
   -m FLOAT                   min modification threshold(s). Comma separated values for each modification code given in -c [0.8]
   -t INT                     number of processing threads [8]
   -K INT                     batch size (max number of reads loaded at once) [512]
   -B FLOAT[K/M/G]            max number of bytes loaded at once [20.0M]
   -h                         help
   -p INT                     print progress every INT seconds (0: per batch) [0]
   -o FILE                    output file [stdout]
   --verbose INT              verbosity level [4]
   --version                  print version
   --insertions               enable modifications in insertions
```

**Sample modfreqs.tsv output**
```bash
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
| 6. n_mod | int | number of reads with base modification |
| 7. freq | float | n_mod/n_called ratio |
| 8. mod_code | char | base modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf) |
| 9. ins_offset | int | offset of inserted base from ref_pos (only output when --insertions is specified)

**Sample modfreqs.bedmethyl output**
```bash
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
| 5. n_mod | int | number of reads with base modification |
| 6. strand | char | strand (+/-) of the read |
| 7. start | int | = field 2 |
| 8. end   | int | = field 3 |
| 9. n_mod | int | = field 5 |
| 10. freq | float | n_mod/n_called ratio |

# Modification codes
Base modification codes can be set for both view and mod-freq tool using -c option.
 Here is an example command to explain all possible context formats
```bash
minimod view -c a(C),h(CG),m,a(*) ref.fa reads.bam
minimod mod-freq -c a(C),h(CG),m,a(*) ref.fa reads.bam
```
Consider following modification 
- type a modifications of all A bases
- type h modifications in CG context (CpG sites)
- type m modifications in default CG context
- type a modifications in all contexts

# Modification threshold
Base modification threshold can be set for mod-freq tool using -m option.


01. 5mC modification(default context :CG) frequencies with threshold 0.8
> ```
>   minimod mod-freq -c m -m 0.8 ref.fa reads.bam
> ```
![Fig](docs/figs/m_threshold.png)
> ```
> If p(5mC) >=  0.8 (threshold),       called(5mC) and modified(5mC)
> If p(5mC) <=  0.2 (1-threshold),     called(5mC)
> else,                                ignored as ambiguous
>
> mod_freq(5mC) = total_modified(5mC)/total_called(5mC)
> ```

03. 5mC and 5hmC base modification(default context :CG) frequencies with thresholds 0.8, 0.7 respectively
> ```
>   minimod mod-freq -c mh -m 0.8,0.7 ref.fa reads.bam
> ```
![Fig](docs/figs/m_h_threshold.png)
> ```
> If p(5mC)  >=  0.8 (threshold),      called(5mC) and modified(5mC)
> If p(5mC)  <=  0.2 (1-threshold),    called(5mC)
> else,                                ambiguous
>
> mod_freq(5mC) = total_modified(5mC)/total_called(5mC)
> ```

> ```
> If p(5hmC) >=  0.7 (threshold),      called(5hmC) and modified(5hmC)
> If p(5hmC) <=  0.3 (1-threshold),    called(5hmC)
> else,                                ignored as ambiguous
>
> mod_freq(5hmC) = total_modified(5hmC)/total_called(5hmC)
> ```

# Modified bases in insertions
minimod can handle insterted modified bases(where canonical base in not in reference) by specifiying --insertions flag for both mod-freq and view tools.

Specifying --insertions will add an extra ins_offset column to the output. Sample outputs are given below.

**Sample output of view with --insertions**
```bash
ref_contig	ref_pos	strand	read_id	read_pos	mod_code	mod_prob	ins_offset
chr22	19966677	+	89870c83-8790-419f-acf8-8a8e93a0f3c9	186	m	0.972549	0
chr22	19966761	+	89870c83-8790-419f-acf8-8a8e93a0f3c9	266	m	0.505882	0
chr22	19966774	+	89870c83-8790-419f-acf8-8a8e93a0f3c9	279	m	0.949020	0
chr22	19966782	+	89870c83-8790-419f-acf8-8a8e93a0f3c9	287	m	0.972549	0
chr22	19966804	+	89870c83-8790-419f-acf8-8a8e93a0f3c9	309	m	0.286275	0
chr22	19966873	+	89870c83-8790-419f-acf8-8a8e93a0f3c9	377	m	0.003922	0
chr22	19966962	+	89870c83-8790-419f-acf8-8a8e93a0f3c9	463	m	0.988235	0
chr22	19967025	+	89870c83-8790-419f-acf8-8a8e93a0f3c9	526	m	0.952941	0
chr22	19967029	+	89870c83-8790-419f-acf8-8a8e93a0f3c9	530	m	0.886275	0
```

**Sample output of mod-freq with --insertions**
```bash
contig	start	end	strand	n_called	n_mod	freq	mod_code	ins_offset
chr22	20016024	20016024	+	2	0	0.000000	m	0
chr22	20017059	20017059	-	4	0	0.000000	m	0
chr22	19989035	19989035	+	3	2	0.666667	m	0
chr22	20016841	20016841	-	5	0	0.000000	m	0
chr22	20016249	20016249	-	1	0	0.000000	m	0
chr22	20003531	20003531	-	1	0	0.000000	m	0
chr22	19975732	19975732	+	1	0	0.000000	m	0
chr22	19995038	19995038	-	1	0	0.000000	m	0
chr22	19999619	19999619	-	9	6	0.666667	m	0
```

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

# Limitations / Future Improvements
- Does not support alignment BAM files with modification codes as numeric ChEBI codes in MM tag
- Status of skipped bases (encoded as . or ? in MM tag) are ignored