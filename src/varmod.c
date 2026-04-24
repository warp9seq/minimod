/**
 * @file varmod.c
 * @brief variant-aware modification tags
 */

#include "varmod.h"
#include "misc.h"
#include "error.h"
#include "khash.h"
#include "ref.h"
#include "ksort.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>


// Fallback for systems/compilers that don't expose drand48
#ifndef drand48
#define drand48() ((double)rand() / RAND_MAX)
#endif

extern const char *get_mm_tag_ptr(bam1_t *record);
extern uint8_t *get_ml_tag(bam1_t *record, uint32_t *len_ptr);
extern uint8_t get_hp_tag(bam1_t *record);

#define WILDCARD_STR "*"
#define THRESH_UINT8_TO_DBL(x) ((double)( (x + 0.5) / 256.0 ))
#define IS_ALPHA(c) (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')
#define IS_DIGIT(c) (c >= '0' && c <= '9')


// zero-allocation comparator
static int cmp_key_fast(const char *key_a, const char *key_b) {
    // find the first tab (end of the contig string)
    const char *tab_a = strchr(key_a, '\t');
    const char *tab_b = strchr(key_b, '\t');

    // calculate the length of the contig portions
    size_t len_a = tab_a ? (size_t)(tab_a - key_a) : strlen(key_a);
    size_t len_b = tab_b ? (size_t)(tab_b - key_b) : strlen(key_b);

    // compare the contigs up to the length of the shorter one
    size_t min_len = len_a < len_b ? len_a : len_b;
    int cmp = strncmp(key_a, key_b, min_len);

    // if contigs are different, we have our answer
    if (cmp != 0) return cmp;
    
    // if the prefixes match but lengths differ, the shorter one comes first
    if (len_a != len_b) return (len_a < len_b) ? -1 : 1;

    // contigs are exactly identical. compare start positions.
    // if there is no tab, they are identical strings.
    if (!tab_a || !tab_b) return 0;

    // atoi automatically stops reading when it hits the next tab
    int start_a = atoi(tab_a + 1);
    int start_b = atoi(tab_b + 1);

    return start_a - start_b;
}

#define varfreq_kv_lt(a, b) (cmp_key_fast((a).key, (b).key) < 0)
#define varview_kv_lt(a, b) (cmp_key_fast((a).key, (b).key) < 0)

KSORT_INIT(varfreq, varfreq_kv_t, varfreq_kv_lt)
KSORT_INIT(varview, varview_kv_t, varview_kv_lt)


static const int valid_bases[256] = { ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1, ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1, ['u'] = 1, ['n'] = 1 };
static const int valid_strands[256] = { ['+'] = 1, ['-'] = 1 };
static const int base_idx_lookup[256] = { ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['U'] = 3, ['N'] = 4, ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3, ['u'] = 3, ['n'] = 4 };
static const char base_complement_lookup[256] = { ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A', ['N'] = 'N', ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a', ['n'] = 'n' };
static const char* tested_cases[] = {"m[CG]", "h[CG]"};


static char* make_key(const char *chrom, int pos, uint16_t ins_offset, char * mod_code, char strand, int haplotype){
    int start_strlen = snprintf(NULL, 0, "%d", pos);
    int offset_strlen = snprintf(NULL, 0, "%d", ins_offset);
    int mod_code_strlen = strlen(mod_code);
    int haplotype_strlen = snprintf(NULL, 0, "%d", haplotype);
    int key_strlen = strlen(chrom) + start_strlen  + offset_strlen + mod_code_strlen + haplotype_strlen + 7;

    char* key = (char *)malloc(key_strlen * sizeof(char));
    MALLOC_CHK(key);
    snprintf(key, key_strlen, "%s\t%d\t%c\t%s\t%u\t%d", chrom, pos, strand, mod_code, ins_offset, haplotype);
    return key;
}

static void decode_key(char *key, char **chrom, int *pos, uint16_t * ins_offset, char **mod_code, char *strand, int *haplotype){
    char* token = strtok(key, "\t");
    *chrom = calloc(strlen(token)+1, sizeof(char));
    MALLOC_CHK(*chrom);
    strcpy(*chrom, token);

    *pos = atoi(strtok(NULL, "\t"));
    *strand = strtok(NULL, "\t")[0];
    
    token = strtok(NULL, "\t");
    *mod_code = calloc(strlen(token)+1, sizeof(char));
    MALLOC_CHK(*mod_code);
    strcpy(*mod_code, token);

    *ins_offset = strtoul(strtok(NULL, "\t"), NULL, 10);
    *haplotype = atoi(strtok(NULL, "\t"));
}


void find_cg_contexts_in_sequence(const char* seq, int seq_len, int** offsets, int* offsets_len, int* offsets_cap) {
    *offsets_len = 0;
    *offsets_cap = 16;
    *offsets = (int*)malloc(sizeof(int) * (*offsets_cap));
    MALLOC_CHK(*offsets);
    
    for(int i = 0; i < seq_len - 1; i++) {
        if((seq[i] == 'C' || seq[i] == 'c') && (seq[i+1] == 'G' || seq[i+1] == 'g')) {
            if(*offsets_len >= *offsets_cap) {
                *offsets_cap *= 2;
                *offsets = (int*)realloc(*offsets, sizeof(int) * (*offsets_cap));
                MALLOC_CHK(*offsets);
            }
            (*offsets)[(*offsets_len)++] = i;
        }
    }
}

void load_var_map(const char* vcf_file, khash_t(varm)* var_map) {
    
    htsFile *vcf_fp = hts_open(vcf_file, "r");
    if(vcf_fp == NULL) {
        ERROR("Failed to open VCF file: %s\n", vcf_file);
        exit(EXIT_FAILURE);
    }
    
    bcf_hdr_t * vcf_hdr = bcf_hdr_read(vcf_fp);
    if(vcf_hdr == NULL) {
        ERROR("Failed to read VCF header: %s\n", vcf_file);
        exit(EXIT_FAILURE);
    }

    bcf1_t *rec = bcf_init();
    
    while(bcf_read(vcf_fp, vcf_hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);

        // fprintf(stderr, "Processing variant: %s:%d %s -> ", bcf_hdr_id2name(vcf_hdr, rec->rid), rec->pos, rec->d.allele[0]);
        // for(int i = 1; i < rec->n_allele; i++) {
        //     fprintf(stderr, "%s ", rec->d.allele[i]);
        // }
        // fprintf(stderr, "\n");
        
        const char *contig = bcf_hdr_id2name(vcf_hdr, rec->rid);
        int32_t pos = rec->pos;
        const char *ref_allele = rec->d.allele[0];
        int ref_len = strlen(ref_allele);
    
        char * contig_dup = (malloc(sizeof(char) * (strlen(contig) + 1)));
        MALLOC_CHK(contig_dup);
        strcpy(contig_dup, contig);
        
        int ret;
        khint_t k = kh_put(varm, var_map, contig_dup, &ret);
        vars_t * vars = kh_value(var_map, k);
        
        // already exists
        if(ret == 0) {
            free(contig_dup);
        } else {
            vars = (vars_t *)malloc(sizeof(vars_t));
            MALLOC_CHK(vars);
            vars->vars_len = 0;
            vars->vars_cap = 1;
            vars->vars = (var_t*)malloc(sizeof(var_t) * vars->vars_cap);
            MALLOC_CHK(vars->vars);
            kh_value(var_map, k) = vars;
        }
        
        for(int i = 1; i < rec->n_allele; i++) {
            const char *alt_allele = rec->d.allele[i];

            if(vars->vars_len >= vars->vars_cap) {
                vars->vars_cap *= 2;
                vars->vars = (var_t*)realloc(vars->vars, sizeof(var_t) * vars->vars_cap);
                MALLOC_CHK(vars->vars);
            }
            var_t var;
            var.start = pos;
            strcpy(var.ref_allele, ref_allele);
            var.ref_len = ref_allele[0] == '.' ? 0 : strlen(ref_allele);
            strcpy(var.alt_allele, alt_allele);
            var.alt_len = alt_allele[0] == '.' ? 0 : strlen(alt_allele);
            // TODO: Refer vcf spec. and handle indels at terminal positions of the contig.
            vars->vars[vars->vars_len++] = var;

        }
    }
    
    bcf_destroy(rec);
    bcf_hdr_destroy(vcf_hdr);
    hts_close(vcf_fp);
    
}

void destroy_var_map(khash_t(varm)* var_map) {
    khint_t k;
    for (k = kh_begin(var_map); k != kh_end(var_map); ++k) {
        if (kh_exist(var_map, k)) {
            char * key = (char*) kh_key(var_map, k);
            vars_t *vars = kh_value(var_map, k);
            free(key);
            free(vars->vars);
            free(vars);
        }
    }
    kh_destroy(varm, var_map);
}


static void get_aln(core_t * core, db_t *db, bam_hdr_t *hdr, bam1_t *record, int bam_i){
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = (tid >= 0) ? hdr->target_name[tid] : "*";
    int32_t pos = record->core.pos;
    int32_t end = bam_endpos(record);

    const char *qname = bam_get_qname(record);

    int8_t rev = bam_is_rev(record);

    uint32_t *cigar = bam_get_cigar(record);
    uint32_t n_cigar = record->core.n_cigar;

    int seq_len = record->core.l_qseq;

    ref_t *ref = get_ref(tname);
        ASSERT_MSG(ref != NULL, "Contig %s not found in reference provided\n", tname);
  
    int read_pos = 0;
    int ref_pos = pos;

    int * aligned_pairs = db->aln[bam_i];
    //fill the aligned_pairs array with -1
    for(int i=0;i<seq_len;i++){
        aligned_pairs[i] = -1;
    }

    for(int i=0;i<seq_len;i++){
        db->ins[bam_i][i] = -1;
        db->ins_offset[bam_i][i] = 0;
    }

    for (uint32_t ci = 0; ci < n_cigar; ++ci) {
        uint32_t c = cigar[ci];
        if(rev) {
            c = cigar[n_cigar - ci - 1];
        }
        int cigar_len = bam_cigar_oplen(c);
        int cigar_op = bam_cigar_op(c);

        // Set the amount that the ref/read positions should be incremented
        // based on the cigar operation
        int read_inc = 0;
        int ref_inc = 0;

        // Process match between the read and the reference
        int8_t is_aligned = 0, is_inserted = 0;
        if(cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            is_aligned = 1;
            read_inc = 1;
            ref_inc = 1;
        } else if(cigar_op == BAM_CDEL) {
            ref_inc = 1;
        } else if(cigar_op == BAM_CREF_SKIP) {
            // end the current segment and start a new one
            //out.push_back(AlignedSegment());
            ref_inc = 1;
        } else if(cigar_op == BAM_CINS) {
            read_inc = 1;
            is_inserted = 1;
        } else if(cigar_op == BAM_CSOFT_CLIP) {
            read_inc = 1;
        } else if(cigar_op == BAM_CHARD_CLIP) { // TODO: use MN tag (seq len at the time MM value was last written) to check this?
            read_inc = 0;
            ERROR("Hard clipping found in %s and they are not supported.\nTry following workarounds.\n\t01. Filter out non-primary alignments\n\t\tsamtools view -h -F 2308 reads.bam -o primary_reads.bam\n\t02. Use minimap2 with -Y to use soft clipping for suplimentary alignments.\n", qname); 
            exit(EXIT_FAILURE);
        } else {
            ERROR("Unhandled CIGAR OPT Cigar: %d\n", cigar_op);
            exit(EXIT_FAILURE);
        }

        // Iterate over the pairs of aligned bases
        for(int j = 0; j < cigar_len; ++j) {
            if(is_aligned) {
                ASSERT_MSG(read_pos < seq_len, "read_pos:%d seq_len:%d\n", read_pos, seq_len);
                int start = ref_pos;
                if(rev) {
                    start = pos + end - ref_pos - 1;
                }
                aligned_pairs[read_pos] = start;

                ASSERT_MSG(ref_pos >= 0 && ref_pos < ref->ref_seq_length, "ref_pos:%d ref_len:%d\n", ref_pos, ref->ref_seq_length);
                ASSERT_MSG(ref->ref_seq_length == hdr->target_len[tid], "ref_len:%d target_len:%d\n", ref->ref_seq_length, hdr->target_len[tid]);
            }

            if(is_inserted) {
                ASSERT_MSG(read_pos < seq_len, "read_pos:%d seq_len:%d\n", read_pos, seq_len);
                int start = ref_pos-1;
                int offset = j+1;
                if(rev) {
                    start = pos + end - ref_pos - 1;
                    offset = cigar_len - j;
                }
                db->ins[bam_i][read_pos] = start;
                db->ins_offset[bam_i][read_pos] = offset;
            }

            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }
}


void update_varfreq_map(khash_t(varfreqm) *varfreq_map, const char *tname, int ref_pos, int ins_offset, char *mod_code, char strand, int haplotype, int is_called, int is_mod) {
    char * key = make_key(tname, ref_pos, ins_offset, mod_code, strand, haplotype);
    khiter_t k = kh_get(varfreqm, varfreq_map, key);
    if (k == kh_end(varfreq_map)) { // not found, add
        varfreq_t * varfreq = (varfreq_t *)malloc(sizeof(varfreq_t));
        MALLOC_CHK(varfreq);
        varfreq->n_called = is_called;
        varfreq->n_mod = is_mod;
        int ret;
        k = kh_put(varfreqm, varfreq_map, key, &ret);
        kh_value(varfreq_map, k) = varfreq;
    } else { // found, update
        free(key);
        varfreq_t * varfreq = kh_value(varfreq_map, k);
        varfreq->n_called += is_called;
        varfreq->n_mod += is_mod;
        // check if varfreq->n_called overflows
        if(varfreq->n_called == 0){
            ERROR("n_called overflowed for key %s. Please report this issue.", key);
            exit(EXIT_FAILURE);
        }
    }

    if(haplotype != -1) {
        char * key = make_key(tname, ref_pos, ins_offset, mod_code, strand, -1); // aggregate all haplotypes
        khiter_t k = kh_get(varfreqm, varfreq_map, key);
        if (k == kh_end(varfreq_map)) { // not found, add
            varfreq_t * varfreq = (varfreq_t *)malloc(sizeof(varfreq_t));
            MALLOC_CHK(varfreq);
            varfreq->n_called = is_called;
            varfreq->n_mod = is_mod;
            int ret;
            k = kh_put(varfreqm, varfreq_map, key, &ret);
            kh_value(varfreq_map, k) = varfreq;
        } else { // found, update
            free(key);
            varfreq_t * varfreq = kh_value(varfreq_map, k);
            varfreq->n_called += is_called;
            varfreq->n_mod += is_mod;
            // check if varfreq->n_called overflows
            if(varfreq->n_called == 0){
                ERROR("n_called overflowed for key %s. Please report this issue.", key);
                exit(EXIT_FAILURE);
            }
        }
    }
}

void add_varview_entry(khash_t(varviewm) *varview_map, const char *tname, int ref_pos, int ins_offset, char *mod_code, char strand, int haplotype, uint8_t mod_prob, int read_pos, var_t var) {

    char *key = make_key(tname, ref_pos, ins_offset, mod_code, strand, haplotype);
    khiter_t k = kh_get(varviewm, varview_map, key);
    if (k == kh_end(varview_map)) { // not found, add
        varview_t *varview = (varview_t *)malloc(sizeof(varview_t));
        MALLOC_CHK(varview);
        varview->mod_prob = mod_prob;
        varview->read_pos = read_pos;
        varview->var = var;
        int ret;
        k = kh_put(varviewm, varview_map, key, &ret);
        kh_value(varview_map, k) = varview;
    } else { // found, update
        free(key);
    }
}

void varviewfreq_single(core_t * core, db_t *db, int32_t bam_i) {
   bam1_t *record = db->bam_recs[bam_i];
    // const char *qname = bam_get_qname(record);
    int8_t rev = bam_is_rev(record);
    bam_hdr_t *hdr = core->bam_hdr;
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = (tid >= 0) ? hdr->target_name[tid] : "*";
    uint8_t *seq = bam_get_seq(record);
    uint32_t seq_len = record->core.l_qseq;
    char strand = rev ? '-' : '+';
    const char *mm_string = db->mm[bam_i];
    uint32_t ml_len = db->ml_lens[bam_i];
    uint8_t *ml = db->ml[bam_i];
    int haplotype = core->opt.haplotypes ? get_hp_tag(record) : -1;
    int * aln_pairs = db->aln[bam_i];

    // get the aligned positions and insertions
    get_aln(core, db, hdr, record, bam_i);
    
    // 5 int arrays to keep base pos of A, C, G, T, N bases.
    // A: 0, C: 1, G: 2, T: 3, U:4, N: 5
    // so that, nth base of A is at base_pos[0][n] and so on.
    int **bases_pos = db->bases_pos[bam_i];
    int bases_pos_lens[N_BASES] = {0};
    memset(db->mod_codes[bam_i], 0, core->opt.n_mods);

    int i;
    for(i=0;i<seq_len;i++){
        int base_char = seq_nt16_str[bam_seqi(seq, i)];
        int idx = base_idx_lookup[(int)base_char];
        bases_pos[idx][bases_pos_lens[idx]++] = i;
    }

    int mm_str_len = strlen(mm_string);
    i = 0;
    int ml_start_idx = 0;

    char modbase;
    // char mod_strand;
    char * mod_codes = db->mod_codes[bam_i];
    int mod_codes_len;
    int * skip_counts = db->skip_counts[bam_i];
    int skip_counts_len;
    char status_flag;

    while (i < mm_str_len) {
        // reset skip counts and mod codes
        skip_counts_len = 0;
        mod_codes_len = 0;

        // set default status flag to '.' (when not present or '.' in the MM string)
        status_flag = '.';

        // get base
        if(i < mm_str_len) {
            ASSERT_MSG(valid_bases[(int)mm_string[i]], "Invalid base:%c\n", mm_string[i]);
            modbase = mm_string[i] == 'U' ? 'T' : mm_string[i]; // convert U to T
            i++;
        }

        // get strand
        if(i < mm_str_len) {
            ASSERT_MSG(valid_strands[(int)mm_string[i]], "Invalid strand:%c\n", mm_string[i]);
            // mod_strand = mm_string[i];
            i++;
        }

        // get base modification codes. can handle multiple codes giver as chars. TO-DO: handle when given as a ChEBI id
        int j = 0;
        int has_nums = 0;
        int has_alpha = 0;
        while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';' && mm_string[i] != '?' && mm_string[i] != '.') {

            // ASSERT_MSG(valid_mod_codes[(int)mm_string[i]], "Invalid base modification code:%c\n", mm_string[i]);

            if(IS_DIGIT(mm_string[i])) {
                has_nums = 1;
            } else if(IS_ALPHA(mm_string[i])) {
                has_alpha = 1;
            } else {
                ERROR("Invalid base modification code:%c. Modification codes should be either numeric or alphabetic.\n", mm_string[i]);
                exit(EXIT_FAILURE);
            }

            if(j >= db->mod_codes_cap[bam_i]) {
                db->mod_codes_cap[bam_i] *= 2;
                db->mod_codes[bam_i] = (char *)realloc(db->mod_codes[bam_i], sizeof(char) * (db->mod_codes_cap[bam_i] + 1)); // +1 for null terminator
                MALLOC_CHK(db->mod_codes[bam_i]);
            }
            mod_codes = db->mod_codes[bam_i];

            mod_codes[j] = mm_string[i];
            j++;
            i++;
        }
        mod_codes[j] = '\0';
        mod_codes_len = j;

        if(has_nums) {
            mod_codes_len = 1; // if chebi id is given, then only one code is present
        }

        // validate mod codes
        ASSERT_MSG(mod_codes_len>0, "Invalid modification codes:%s. Modification codes cannot be empty.\n", mod_codes);
        ASSERT_MSG((has_nums && has_alpha) == 0, "Invalid modification codes:%s. Modification codes should be either numeric or alphabetic, not both.\n", mod_codes);
        
        // get modification status flag
        if(i < mm_str_len && ( mm_string[i] == '?' || mm_string[i] == '.' )) {
            status_flag = mm_string[i];
            i++;
        } else { // if not present, set to '.'
            status_flag = '.';
        }

        // get skip counts
        int k = 0;
        while (i < mm_str_len && mm_string[i] != ';') {

            // skip if a comma
            if(i < mm_str_len && mm_string[i] == ',') {
                i++;
                continue;
            }

            char skip_count_str[10];
            int l = 0;
            while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';') {
                skip_count_str[l] = mm_string[i];
                i++;
                l++;
                assert(l < 10); // if this fails, use dynamic allocation for skip_count_str
            }
            skip_count_str[l] = '\0';
            ASSERT_MSG(l > 0, "Invalid skip count:%d.\n", l);
            sscanf(skip_count_str, "%d", &skip_counts[k]);
            ASSERT_MSG(skip_counts[k] >= 0, "Skip count cannot be negative: %d.\n", skip_counts[k]);
            
            k++;
        }
        skip_counts_len = k;
        i++;

        char mb = rev? base_complement_lookup[(int)modbase] : modbase;
        int idx = base_idx_lookup[(int)mb];
        
        int ml_idx = ml_start_idx;
        int base_rank = -1; // 0-based rank
        for(int c=0; c<skip_counts_len; c++) {
            base_rank += skip_counts[c] + 1;

            int read_pos;

            if (modbase == 'N') {
                if(rev) {
                    read_pos = seq_len - base_rank - 1;
                } else {
                    read_pos = base_rank;
                }
            } else {
                if(rev) {
                    read_pos = bases_pos[idx][bases_pos_lens[idx] - base_rank - 1];
                } else {
                    read_pos = bases_pos[idx][base_rank];
                }
            }

            ASSERT_MSG(read_pos>=0 && read_pos < seq_len, "Read pos cannot exceed seq len. read_pos: %d seq_len: %d\n", read_pos, seq_len);

            int fastq_read_pos = rev ? (seq_len - read_pos -1) : read_pos;

            int ref_pos = aln_pairs[fastq_read_pos];

            if(ref_pos == -1 && db->ins[bam_i][fastq_read_pos] == -1) { // not aligned nor insertion
                if(mod_codes_len > 0) {
                    ml_idx = ml_start_idx + c*mod_codes_len + mod_codes_len - 1;
                }
                continue;
            }

            int ins_start = db->ins[bam_i][fastq_read_pos];
            int ins_offset = db->ins_offset[bam_i][fastq_read_pos];

            // char out_strand = strand;
            // if(mod_strand == '-') {
            //     out_strand = strand == '+' ? '-' : '+';
            // }
            
            // mod prob per each mod code.
            for(int m=0; m<mod_codes_len; m++) {                
                ml_idx = ml_start_idx + c*mod_codes_len + m;

                // check required mod codes
                khint_t mk;

                mk = kh_get(modcodesm, core->opt.modcodes_map, WILDCARD_STR); // check for wildcard first
                char * mod_code = NULL;
                if (has_nums) { // chebi id
                    mod_code = mod_codes;
                } else { // not chebi id, need to check for each mod code
                    mod_code = &(mod_codes[m]);
                }
                if (mk != kh_end(core->opt.modcodes_map)) { // wildcard present, all mod codes are required
                    // do nothing, just proceed
                } else {
                    mk = kh_get(modcodesm, core->opt.modcodes_map, mod_code);
                    if(mk == kh_end(core->opt.modcodes_map)) continue; // mod code not required
                }

                modcodem_t *req_mod = kh_value(core->opt.modcodes_map, mk);

                int is_in_context = 0;
                var_t var;
                khint_t ck = kh_get(varm, core->var_map, tname);
                if(ck != kh_end(core->var_map)) {
                    vars_t *vars = kh_value(core->var_map, ck);
                    
                    for(int v=0; v<vars->vars_len; v++) {
                        var = vars->vars[v];

                        int ins_pos = ins_start + ins_offset;
                        int var_end = var.start + var.alt_len;
                        if(ref_pos == -1 && ins_pos != -1 && ins_pos >= var.start && ins_pos < var_end) {
                            is_in_context = 1;
                            break;
                        } else if(ref_pos != -1 && ref_pos >= var.start && ref_pos < var_end) {
                            is_in_context = 1;
                            break;
                        }
                        
                        // if(var.alt_len > var.ref_len) { // insertion
                        //     fprintf(stderr, "Checking var at %s:%d ref_len:%d alt_len:%d for read %s at ref_pos:%d ins_pos:%d ins_offset:%d\n", tname, var.start, var.ref_len, var.alt_len, bam_get_qname(record), ref_pos, ins_start + ins_offset, ins_offset);
                        //     int ins_pos = ins_start + ins_offset;
                        //     int var_end = var.start + var.alt_len;
                        //     if(ref_pos == -1 && ins_pos >= var.start && ins_pos < var_end) {
                        //         is_in_context = 1;
                        //         break;
                        //     }
                        // } else if(var.alt_len < var.ref_len) { // deletion
                        //     break; // currently not handling deletions
                        // } else { // SNP
                        //     break; // currently not handling SNPs
                        // }
                    }
                }

                if (is_in_context) { // in context and mod_base matches reference
                } else {
                    continue;
                }

                ASSERT_MSG(ml_idx<ml_len, "read_id:%s mod prob index mismatch. ml_idx:%d ml_len:%d \n", bam_get_qname(record), ml_idx, ml_len);
                uint8_t mod_prob = ml[ml_idx];
                ASSERT_MSG(mod_prob <= 255 && mod_prob>=0, "Invalid mod_prob:%d\n", mod_prob);

                
                if(core->opt.subtool == FREQ) {
                    uint8_t is_mod = 0, is_called = 0;
                    double thresh = req_mod->thresh;
                    double mod_prob_dbl = THRESH_UINT8_TO_DBL(mod_prob);
                    
                    if(mod_prob_dbl >= thresh){ // modified with mod_code
                        is_called = 1;
                        is_mod = 1;
                    } else if(mod_prob_dbl <= 1 - thresh){ // not modified with mod_code
                        is_called = 1;
                    } else { // ambiguous
                        continue;
                    }
                    
                    update_varfreq_map(db->varfreq_maps[bam_i], tname, ref_pos==-1?ins_start:ref_pos, ins_offset, mod_code, strand, haplotype, is_called, is_mod);
                } else if (core->opt.subtool == VARVIEW) {
                    add_varview_entry(db->varview_maps[bam_i], tname, ref_pos==-1?ins_start:ref_pos, ins_offset, mod_code, strand, haplotype, mod_prob, fastq_read_pos, var);
                }
            }

        }
        if(skip_counts_len > 0) ml_start_idx = ml_idx + 1;

        
        // skipped bases
        if (status_flag == '.') {
            int skip_base_rank = -1; // 0-based rank
            int prev_skip_base_rank = -1; 
            for(int c=0; c<skip_counts_len; c++) {
                skip_base_rank += skip_counts[c] + 1;
            
                // modification probability is 0 for skipped bases
                for(int s=prev_skip_base_rank+1; s<skip_base_rank; s++) {
                    int skip_read_pos;
                    if (modbase == 'N') {
                        if(rev) {
                            skip_read_pos = seq_len - s - 1;
                        } else {
                            skip_read_pos = s;
                        }
                    } else {
                        if(rev) {
                            skip_read_pos = bases_pos[idx][bases_pos_lens[idx] - s - 1];
                        } else {
                            skip_read_pos = bases_pos[idx][s];
                        }
                    }

                    ASSERT_MSG(skip_read_pos>=0 && skip_read_pos < seq_len, "Read pos cannot exceed seq len. read_pos: %d seq_len: %d\n", skip_read_pos, seq_len);

                    int skip_fastq_read_pos = rev ? (seq_len - skip_read_pos -1) : skip_read_pos;

                    int skip_ref_pos = aln_pairs[skip_fastq_read_pos];

                    if(skip_ref_pos == -1 && db->ins[bam_i][skip_fastq_read_pos] == -1) { // not aligned nor insertion
                        if(mod_codes_len > 0) {
                            ml_idx = ml_start_idx + c*mod_codes_len + mod_codes_len - 1;
                        }
                        continue;
                    }

                    int ins_start = db->ins[bam_i][skip_read_pos];
                    int ins_offset = db->ins_offset[bam_i][skip_read_pos];

                    // mod prob per each mod code.
                    for(int m=0; m<mod_codes_len; m++) {

                        // check required mod codes
                        khint_t mk;

                        mk = kh_get(modcodesm, core->opt.modcodes_map, WILDCARD_STR); // check for wildcard first
                        char * mod_code = NULL;
                        if (has_nums) { // chebi id
                            mod_code = mod_codes;
                        } else { // not chebi id, need to check for each mod code
                            mod_code = &(mod_codes[m]);
                        }
                        if (mk != kh_end(core->opt.modcodes_map)) { // wildcard present, all mod codes are required
                            // do nothing, just proceed
                        } else {
                            mk = kh_get(modcodesm, core->opt.modcodes_map, mod_code);
                            if(mk == kh_end(core->opt.modcodes_map)) {
                                continue; // mod code not required
                            }
                        }

                        int skip_is_in_context = 0;
                        var_t var;
                        khint_t ck = kh_get(varm, core->var_map, tname);
                        if(ck != kh_end(core->var_map)) {
                            vars_t *vars = kh_value(core->var_map, ck);
                            for(int v=0; v<vars->vars_len; v++) {
                                var = vars->vars[v];
                                int ins_pos = ins_start + ins_offset;
                                int var_end = var.start + var.alt_len;
                                if(skip_ref_pos == -1 && ins_pos != -1 && ins_pos >= var.start && ins_pos < var_end) {
                                    skip_is_in_context = 1;
                                    break;
                                } else if(skip_ref_pos != -1 && skip_ref_pos >= var.start && skip_ref_pos < var_end) {
                                    skip_is_in_context = 1;
                                    break;
                                }
                                // if(var.alt_len > var.ref_len) { // insertion
                                //     int ins_pos = ins_start + ins_offset;
                                //     int var_end = var.start + var.alt_len;
                                //     if(skip_ref_pos == -1 && ins_pos >= var.start && ins_pos < var_end) {
                                //         skip_is_in_context = 1;
                                //         break;
                                //     }
                                // } else if(var.alt_len < var.ref_len) { // deletion
                                //     break; // currently not handling deletions
                                // } else { // SNP
                                //     break; // currently not handling SNPs
                                // }
                            }
                        }

                        if (skip_is_in_context) { // in context and mod_base matches reference
                        } else {
                            continue;
                        }

                        if(core->opt.subtool == FREQ) {
                            uint8_t is_mod = 0, is_called = 1; // skipped bases are called as unmodified
                            update_varfreq_map(db->varfreq_maps[bam_i], tname, skip_ref_pos==-1?ins_start:skip_ref_pos, ins_offset, mod_code, strand, haplotype, is_called, is_mod);
                        } else if (core->opt.subtool == VARVIEW) {
                            add_varview_entry(db->varview_maps[bam_i], tname, skip_ref_pos==-1?ins_start:skip_ref_pos, ins_offset, mod_code, strand, haplotype, 0, skip_fastq_read_pos, var);
                        }
                    }
                }
                prev_skip_base_rank = skip_base_rank;
            }

            // handle skipped bases after the last skip count
            for(int s=prev_skip_base_rank+1; s<bases_pos_lens[idx]; s++) {
                int skip_read_pos;
                if (modbase == 'N') {
                    if(rev) {
                        skip_read_pos = seq_len - s - 1;
                    } else {
                        skip_read_pos = s;
                    }
                } else {
                    if(rev) {
                        skip_read_pos = bases_pos[idx][bases_pos_lens[idx] - s - 1];
                    } else {
                        skip_read_pos = bases_pos[idx][s];
                    }
                }

                ASSERT_MSG(skip_read_pos>=0 && skip_read_pos < seq_len, "Read pos cannot exceed seq len. read_pos: %d seq_len: %d\n", skip_read_pos, seq_len);

                int skip_fastq_read_pos = rev ? (seq_len - skip_read_pos -1) : skip_read_pos;

                int skip_ref_pos = aln_pairs[skip_fastq_read_pos];
                skip_ref_pos = aln_pairs[skip_fastq_read_pos];


                int ins_start = db->ins[bam_i][skip_fastq_read_pos];
                int ins_offset = db->ins_offset[bam_i][skip_fastq_read_pos];

                // mod prob per each mod code.
                for(int m=0; m<mod_codes_len; m++) {

                    // check required mod codes
                    khint_t mk;

                    mk = kh_get(modcodesm, core->opt.modcodes_map, WILDCARD_STR); // check for wildcard first
                    char * mod_code = NULL;
                    if (has_nums) { // chebi id
                        mod_code = mod_codes;
                    } else { // not chebi id, need to check for each mod code
                        mod_code = &(mod_codes[m]);
                    }
                    if (mk != kh_end(core->opt.modcodes_map)) { // wildcard present, all mod codes are required
                        // do nothing, just proceed
                    } else {
                        mk = kh_get(modcodesm, core->opt.modcodes_map, mod_code);
                        if(mk == kh_end(core->opt.modcodes_map)) {
                            continue; // mod code not required
                        }
                    }

                    int skip_is_in_context = 0;
                    var_t var;
                    khint_t ck = kh_get(varm, core->var_map, tname);
                    if(ck != kh_end(core->var_map)) {
                        vars_t *vars = kh_value(core->var_map, ck);
                        
                        for(int v=0; v<vars->vars_len; v++) {
                            for(int v=0; v<vars->vars_len; v++) {
                                var = vars->vars[v];
                                int ins_pos = ins_start + ins_offset;
                                int var_end = var.start + var.alt_len;
                                if(skip_ref_pos == -1 && ins_pos != -1 && ins_pos >= var.start && ins_pos < var_end) {
                                    skip_is_in_context = 1;
                                    break;
                                } else if(skip_ref_pos != -1 && skip_ref_pos >= var.start && skip_ref_pos < var_end) {
                                    skip_is_in_context = 1;
                                    break;
                                }
                                // if(var.alt_len > var.ref_len) { // insertion
                                //     int ins_pos = ins_start + ins_offset;
                                //     int var_end = var.start + var.alt_len;
                                //     if(skip_ref_pos == -1 && ins_pos >= var.start && ins_pos < var_end) {
                                //         skip_is_in_context = 1;
                                //         break;
                                //     }
                                // } else if(var.alt_len < var.ref_len) { // deletion
                                //     break; // currently not handling deletions
                                // } else { // SNP
                                //     break; // currently not handling SNPs
                                // }
                            }
                        }
                    }

                    if (skip_is_in_context) { // in context and mod_base matches reference
                    } else {
                        continue;
                    }

                    if(core->opt.subtool == FREQ) {
                        uint8_t is_mod = 0, is_called = 1; // skipped bases are called as unmodified
                        update_varfreq_map(db->varfreq_maps[bam_i], tname, skip_ref_pos==-1?ins_start:skip_ref_pos, ins_offset, mod_code, strand, haplotype, is_called, is_mod);
                    } else if (core->opt.subtool == VARVIEW) {
                        add_varview_entry(db->varview_maps[bam_i], tname, skip_ref_pos==-1?ins_start:skip_ref_pos, ins_offset, mod_code, strand, haplotype, 0, skip_fastq_read_pos, var);
                    }
                }
            }
        
        }
        
    }
}


void warn_untested_cases_var(opt_t * opt) {
    int n_tested_cases = sizeof(tested_cases) / sizeof(tested_cases[0]);
    for(int i=0; i < opt->n_mods; i++) {
        khint_t i;
        for(i=kh_begin(opt->modcodes_map); i < kh_end(opt->modcodes_map); ++i) {
            if (!kh_exist(opt->modcodes_map, i)) continue;
            char * mod_code = (char *) kh_key(opt->modcodes_map, i);
            char * context = kh_value(opt->modcodes_map, i)->context;

            char * mod_code_with_context = (char *)malloc(strlen(mod_code) + strlen(context) + 3); // for null terminator and brackets
            MALLOC_CHK(mod_code_with_context);
            snprintf(mod_code_with_context, strlen(mod_code) + strlen(context) + 3, "%s[%s]", mod_code, context);

            int is_tested = false;
            for(int j=0; j < n_tested_cases; j++){
                if(strcmp(mod_code_with_context, tested_cases[j]) == 0){
                    is_tested = true;
                    break;
                }
            }
            if(!is_tested){
                WARNING("Modification code with context %s has not been tested.", mod_code_with_context);
            }
            free(mod_code_with_context);
        }
    }
}

void print_varview_header(core_t* core) {
    char * common = "ref_contig\tref_pos\tstrand\tread_id\tread_pos\tmod_code\tmod_prob\tins_offset\tref_allele\talt_allele";
    char * ins_offset = "";
    char * haplotype = "";
    if(core->opt.haplotypes){
        haplotype = "\thaplotype";
    }

    fprintf(core->opt.output_fp, "%s%s%s\n", common, ins_offset, haplotype);
}

void print_varview_output(core_t* core, db_t* db) {
    FILE *out_fp = core->opt.output_fp;
    int do_insertions = core->opt.subtool == VARVIEW;
    int do_haplotypes = core->opt.haplotypes == 1;

    int is_bed = true;

    // Reusable buffer
    int max_arr_capacity = 0;
    varview_kv_t *sorted_arr = NULL;

    for(int i = 0; i < db->n_bam_recs; i++) {
        bam1_t *record = db->bam_recs[i];
        const char *qname = bam_get_qname(record);
        khash_t(varviewm) *varview_map = db->varview_maps[i];
        khint_t map_size = kh_size(varview_map);

        if (map_size == 0) continue;

        if (map_size > max_arr_capacity) {
            max_arr_capacity = map_size;
            sorted_arr = (varview_kv_t *)realloc(sorted_arr, sizeof(varview_kv_t) * max_arr_capacity);
            MALLOC_CHK(sorted_arr);
        }

        int size = 0;
        for (khint_t k = kh_begin(varview_map); k != kh_end(varview_map); k++) {
            if (kh_exist(varview_map, k)) {
                sorted_arr[size].key = (char *)kh_key(varview_map, k);
                sorted_arr[size].view = kh_value(varview_map, k);
                size++;
            }
        }

        // qsort(sorted_arr, size, sizeof(varview_kv_t), cmp_varview_kv);
        ks_introsort_varview(size, sorted_arr);

        for (int j = 0; j < size; j++) {
            varview_t* varview = sorted_arr[j].view;
            char *tname = NULL;
            int ref_pos;
            uint16_t ins_offset;
            char *mod_code;
            char strand;
            int haplotype;
            char * key = sorted_arr[j].key;
            decode_key(key, &tname, &ref_pos, &ins_offset, &mod_code, &strand, &haplotype);

            fprintf(out_fp, "%s\t%d\t%c\t%s\t%d\t%s\t%f", tname, ref_pos, strand, qname, varview->read_pos, mod_code, THRESH_UINT8_TO_DBL(varview->mod_prob));
            if(do_insertions){
                fprintf(out_fp, "\t%d", db->ins_offset[i][varview->read_pos]);
            }
            if(do_haplotypes){
                fprintf(out_fp, "\t%d", haplotype);
            }
            fprintf(out_fp, "%s\t%s", varview->var.ref_allele, varview->var.alt_allele);
            fputc('\n', out_fp);
            free(tname);
            free(mod_code);
        }
    }

    if (sorted_arr) {
        free(sorted_arr);
    }

    if(out_fp != stdout){
        fclose(out_fp);
    }
}
