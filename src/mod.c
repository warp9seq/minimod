/**
 * @file pulse.c
 * @brief modification tags

MIT License

Copyright (c) 2023 Hasindu Gamaarachchi
Copyright (c) 2024 Suneth Samarasinghe

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


******************************************************************************/

#include "mod.h"
#include "misc.h"
#include "error.h"
#include "khash.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>

// #include <sys/wait.h>
// #include <unistd.h>

#define INIT_MODS 100
#define INIT_SKIP_COUNTS 100
#define INIT_MOD_CODES 2
#define INIT_BASE_POS 100
#define INIT_MOD_BASES 2
#define MOD_THRESHOLD 0.2

typedef struct {
    char base;
    char strand;
    char * mod_codes;
    int mod_codes_cap;
    int mod_codes_len;
    int * skip_counts;
    int skip_counts_cap;
    int skip_counts_len;
    char status_flag;
} mod_tag_t; // free in free_mod_tags()

typedef struct {
    char mod_code;
    int mod_strand;
    double mod_prob;
    char mod_base;
    void * base;
} mod_t;

typedef struct {
    char base;
    uint8_t qual;
    const char * chrom;
    int ref_pos;
    mod_t * mods;
    int mods_cap;
    int mods_len;
    char strand;
    int * is_skipped; // for 5 mod codes;
    int depth;
    int * is_called; // for 5 mod codes;
} base_t; // free in free_mods_per_base()

enum MOD_CODES {
    MOD_5mC = 'm',
    MOD_5hmC = 'h',
    MOD_5fC = 'f',
    MOD_5caC = 'c',
    MOD_xC = 'C'
};

typedef struct {
    char * chrom;
    int start;
    int end;
    int depth;
    int n_mod;
    int n_called;
    int n_skipped;
    double freq;
    char base;
    char mod_base;
    char mod_code;
    char mod_strand;
    char strand;
} stat_t;

KHASH_MAP_INIT_STR(str, stat_t *);

khash_t(str)* stats_map;
char* make_key(const char *chrom, int start, int end, char mod_code, char strand){
    int start_strlen = snprintf(NULL, 0, "%d", start);
    int end_strlen = snprintf(NULL, 0, "%d", end);
    int key_strlen = strlen(chrom) + start_strlen + end_strlen  + 7;
    
    char* key = (char *)malloc(key_strlen * sizeof(char));
    MALLOC_CHK(key);
    snprintf(key, key_strlen, "%s\t%d\t%d\t%c\t%c", chrom, start, end, mod_code, strand);
    return key;
}

static inline int isValidBase(char ch) {
    ch = toupper(ch);
    return (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'U' || ch == 'N');
}

static inline int isValidStrand(char ch) {
    return (ch == '+' || ch == '-');
}

static inline int isValidModificationCode(char ch) {
    return (ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'z');
}

static inline int die(const char *format, ...) {
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(EXIT_FAILURE);
}

/*
* Extracts the modifications from the MM string
* @param mm_string MM string
* @param len_mods pointer to the variable to store the number of modifications
* @return pointer to the array of mod_tag_t
*/
static mod_tag_t *extract_mods(const char *mm_string, uint32_t *len_mods) {
    if (mm_string == NULL || strlen(mm_string) == 0) {
        WARNING("%s","Empty MM string. Continuing\n");
        return NULL;
    }

    // allocate initial memory for modifications
    int mods_cap = INIT_MODS;
    mod_tag_t * mod_tags = (mod_tag_t *) malloc(mods_cap * sizeof(mod_tag_t));
    MALLOC_CHK(mod_tags);
    int num_mods = 0;

    int mm_str_len = strlen(mm_string);
    int i = 0;
    while (i < mm_str_len) {

        if (num_mods >= INIT_MODS) {
            mods_cap *= 2;
            mod_tags = (mod_tag_t *) realloc(mod_tags, mods_cap * sizeof(mod_tag_t));
            MALLOC_CHK(mod_tags);
        }

        mod_tag_t current_mod;
        memset(&current_mod, 0, sizeof(mod_tag_t));

        // allocate initial memory for skip counts
        current_mod.skip_counts_cap = INIT_SKIP_COUNTS;
        current_mod.skip_counts = (int *) malloc(current_mod.skip_counts_cap * sizeof(int));
        MALLOC_CHK(current_mod.skip_counts);

        // allocate initial memory for modification codes
        current_mod.mod_codes_cap = INIT_MOD_CODES;
        current_mod.mod_codes = (char *) malloc(current_mod.mod_codes_cap * sizeof(char));
        MALLOC_CHK(current_mod.mod_codes);

        // set default status flag to '.' (when not present or '.' in the MM string)
        current_mod.status_flag = '.';

        // get base
        if(i < mm_str_len) {
            current_mod.base = mm_string[i];
            ASSERT_MSG(isValidBase(current_mod.base), "Invalid base:%c\n", current_mod.base);
            i++;
        }

        // get strand
        if(i < mm_str_len) {
            current_mod.strand = mm_string[i];
            ASSERT_MSG(current_mod.strand == '+' || current_mod.strand == '-', "Invalid strand:%c\n", current_mod.strand);
            i++;
        }

        // get base modification codes. can handle multiple codes giver as chars. TO-DO: handle when given as a ChEBI id
        int j = 0;
        while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';' && mm_string[i] != '?' && mm_string[i] != '.') {

            if (j >= current_mod.mod_codes_cap) {
                current_mod.mod_codes_cap *= 2;
                current_mod.mod_codes = (char *) realloc(current_mod.mod_codes, current_mod.mod_codes_cap * sizeof(char));
                MALLOC_CHK(current_mod.mod_codes);
            }

            ASSERT_MSG(isValidModificationCode(mm_string[i]), "Invalid base modification code:%c\n", mm_string[i]);

            current_mod.mod_codes[j] = mm_string[i];
            j++;

            i++;
        }
        current_mod.mod_codes[j] = '\0';
        current_mod.mod_codes_len = j;

        // get modification status flag
        if(i < mm_str_len && ( mm_string[i] == '?' || mm_string[i] == '.' )) {
            current_mod.status_flag = mm_string[i];
            i++;
        } else { // if not present, set to '.'
            current_mod.status_flag = '.';
        }

        // get skip counts
        int k = 0;
        while (i < mm_str_len && mm_string[i] != ';') {

            // skip if a comma
            if(i < mm_str_len && mm_string[i] == ',') {
                i++;
                continue;
            }

            if (k >= current_mod.skip_counts_cap) {
                current_mod.skip_counts_cap *= 2;
                current_mod.skip_counts = (int *) realloc(current_mod.skip_counts, current_mod.skip_counts_cap * sizeof(int));
                MALLOC_CHK(current_mod.skip_counts);
            }


            char skip_count_str[10];
            int l = 0;
            while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';') {
                // printf("mm_string[%d]: %c\n", i, mm_string[i]);
                skip_count_str[l] = mm_string[i];
                i++;
                l++;
                assert(l < 10); // if this fails, use dynamic allocation for skip_count_str
            }
            skip_count_str[l] = '\0';
            // printf("skip_count_str: %s\n", skip_count_str);
            ASSERT_MSG(l > 0, "invalid skip count:%d.\n", l);
            sscanf(skip_count_str, "%d", &current_mod.skip_counts[k]);
            ASSERT_MSG(current_mod.skip_counts[k] >= 0, "skip count cannot be negative: %d.\n", current_mod.skip_counts[k]);
            
            k++;
        }
        current_mod.skip_counts_len = k;
        i++;

        if(current_mod.skip_counts_len == 0) { // no skip counts, no modification
            free(current_mod.skip_counts);
            free(current_mod.mod_codes);
            continue;
        }

        mod_tags[num_mods] = current_mod;
        num_mods++;
    }

    *len_mods = num_mods;

    return mod_tags;

}

// function to print an array of given type
void print_array(void *array, int len, char type){
    if(type == 'i'){
        int *arr = (int *)array;
        for(int i=0;i<len;i++){
            fprintf(stderr, "%d,", arr[i]);
        }
    }else if(type == 'c'){
        char *arr = (char *)array;
        for(int i=0;i<len;i++){
            fprintf(stderr, "%c,", arr[i]);
        }
    }
    fprintf(stderr, "\n");
}

int base_to_idx(char base){
    if(base == 'A'){
        return 0;
    }else if(base == 'C'){
        return 1;
    }else if(base == 'G'){
        return 2;
    }else if(base == 'T'){
        return 3;
    }
    return 4;
}

int mod_code_to_idx(char mod_code){
    switch (mod_code) {
        case MOD_5mC:
            return 0;
        case MOD_5hmC:
            return 1;
        case MOD_5fC:
            return 2;
        case MOD_5caC:
            return 3;
        case MOD_xC:
            return 4;
        default:
            return -1;
    }
}

char base_complement(char base){
    if(base == 'A'){
        return 'T';
    }else if(base == 'C'){
        return 'G';
    }else if(base == 'G'){
        return 'C';
    }else if(base == 'T'){
        return 'A';
    }
    return base;
}

static void update_stats(base_t *bases, uint32_t seq_len, khash_t(str)* stats){
    for(int i=0;i<seq_len;i++){
        base_t base = bases[i];
        if(base.ref_pos == -1){ //not aligned
            continue;
        }
        for(int j=0;j<base.mods_len;j++){
            mod_t mod = base.mods[j];

            char *key = make_key(base.chrom, base.ref_pos, base.ref_pos, mod.mod_code, base.strand);
            khiter_t k = kh_get(str, stats, key);
            if (k == kh_end(stats)) {
                stat_t * stat = (stat_t *)malloc(sizeof(stat_t));
                MALLOC_CHK(stat);
                stat->chrom = (char *)malloc(strlen(base.chrom)+1);
                MALLOC_CHK(stat->chrom);
                strcpy(stat->chrom, base.chrom);
                stat->start = base.ref_pos;
                stat->end = base.ref_pos;
                stat->mod_code = mod.mod_code;
                stat->strand = base.strand;

                stat->mod_base = mod.mod_base;
                stat->base = base.base; //TO-do; take from ref
                stat->n_called = base.is_called[mod_code_to_idx(mod.mod_code)];
                stat->n_skipped = base.is_skipped[mod_code_to_idx(mod.mod_code)];
                stat->n_mod = mod.mod_prob >= MOD_THRESHOLD ? 1 : 0;
                stat->mod_strand = mod.mod_strand;
                stat->depth = base.depth;

                int ret;
                k = kh_put(str, stats, key, &ret);
                kh_value(stats, k) = stat;
            } else {
                free(key);
                stat_t * stat = kh_value(stats, k);
                stat->n_called += base.is_called[mod_code_to_idx(mod.mod_code)];
                stat->n_skipped += base.is_skipped[mod_code_to_idx(mod.mod_code)];
                stat->n_mod += mod.mod_prob >= MOD_THRESHOLD ? 1 : 0;
                stat->depth += base.depth;
            }
        }
    }
}

static int * get_aln(bam_hdr_t *hdr, bam1_t *record){
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    int32_t pos = record->core.pos;
    int32_t end = bam_endpos(record);
    const char *qname = bam_get_qname(record);

    int8_t rev = bam_is_rev(record);
    ASSERT_MSG(!(record->core.flag & BAM_FUNMAP), "Unmapped read %s\n", qname);

    uint32_t *cigar = bam_get_cigar(record);
    uint32_t n_cigar = record->core.n_cigar;

    int seq_len = record->core.l_qseq;
    
    int read_pos = 0;
    int ref_pos = pos;

    int * aligned_pairs = (int *)malloc(sizeof(int)*seq_len);
    MALLOC_CHK(aligned_pairs);

    //fill the aligned_pairs array with -1
    for(int i=0;i<seq_len;i++){
        aligned_pairs[i] = -1;
    }

    //fprintf(stderr,"n cigar: %d\n", n_cigar);
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
        int8_t is_aligned = 0;
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
        } else if(cigar_op == BAM_CSOFT_CLIP) {
            read_inc = 1;
        } else if(cigar_op == BAM_CHARD_CLIP) {
            read_inc = 0;
            ERROR("Hard clipping(%d) not supported. Use minimap2 with -Y to use soft clipping for suplimentary alignment.\n", cigar_op);
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
            }

            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }

    return aligned_pairs;
}

static base_t * get_bases(mod_tag_t *mod_tags, uint32_t mods_len, uint8_t * ml, uint32_t ml_len, int * aln_pairs, bam_hdr_t *hdr, bam1_t *record){

    const char *qname = bam_get_qname(record);
    int8_t rev = bam_is_rev(record);

    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";

    uint8_t *seq = bam_get_seq(record);
    uint32_t seq_len = record->core.l_qseq;

    char strand = rev ? '-' : '+';

    if(seq_len == 0){
        ERROR("Sequence length is 0 for read %s. Found %d mod_tags", qname, mods_len);
        exit(EXIT_FAILURE);
    }

    base_t *bases = (base_t *)malloc(sizeof(base_t)*seq_len);
    MALLOC_CHK(bases);
    for(int i=0;i<seq_len;i++){
        bases[i].mods_cap = INIT_MOD_BASES;
        bases[i].mods = (mod_t *)malloc(sizeof(mod_t)*bases[i].mods_cap);
        MALLOC_CHK(bases[i].mods);
        bases[i].mods_len = 0;
        bases[i].ref_pos = aln_pairs[i];
        bases[i].depth = aln_pairs[i] == -1 ? 0 : 1; // count towards depth only if aligned
        bases[i].chrom = tname;
        bases[i].qual = bam_get_qual(record)[i];
        bases[i].strand = strand;
        bases[i].base = seq_nt16_str[bam_seqi(seq, i)];
        bases[i].is_skipped = (int *)malloc(sizeof(int)*5);
        MALLOC_CHK(bases[i].is_skipped);
        bases[i].is_called = (int *)malloc(sizeof(int)*5);
        MALLOC_CHK(bases[i].is_called);
        for(int j=0;j<5;j++){
            bases[i].is_skipped[j] = 0;
            bases[i].is_called[j] = 0;
        }
    }

    // 5 int arrays to keep base pos of A, C, G, T, N bases.
    // A: 0, C: 1, G: 2, T: 3, N: 4
    // so that, nth base of A is at base_pos[0][n] and so on.
    int * bases_pos[5];
    int bases_pos_lens[5];
    int bases_pos_caps[5];
    for(int i=0;i<5;i++){
        bases_pos_caps[i] = INIT_BASE_POS;
        bases_pos[i] = (int *)malloc(sizeof(int)*bases_pos_caps[i]);
        MALLOC_CHK(bases_pos[i]);
        bases_pos_lens[i] = 0;
    }

    // keep record of base pos of A, C, G, T, N bases
    for(int i=0; i<seq_len; i++){
        int base_char = seq_nt16_str[bam_seqi(seq, i)];
        int idx = base_to_idx(base_char);
        if(bases_pos_lens[idx] >= bases_pos_caps[idx]){
            bases_pos_caps[idx] *= 2;
            bases_pos[idx] = (int *)realloc(bases_pos[idx], sizeof(int)*bases_pos_caps[idx]);
            MALLOC_CHK(bases_pos[idx]);
        }
        bases_pos[idx][bases_pos_lens[idx]] = i;
        bases_pos_lens[idx]++;
        
    }


    // ASSERT_MSG(mods_len >= ml_len, "Probaility array len mismatch. mods_len:%d ml_len:%d\n", mods_len, ml_len);
    
    int ml_start_idx = 0;
    // go through mod_tags    
    for(int i=0; i<mods_len; i++) {
        mod_tag_t mod = mod_tags[i];
        int base_rank = -1;

        int prev_read_pos = -1;

        int ml_idx = ml_start_idx;
        for(int j=0; j<mod.skip_counts_len; j++) {
            base_rank += mod.skip_counts[j] + 1;
            // fprintf(stderr, "qname:%s seq_len:%d mod.base:%c seq.strand:%c mod.strand:%c rank:%d\n", qname, seq_len, mod.base, strand, mod.strand, base_rank);
            
            char mod_base;
            int idx;
            int read_pos;

            if(rev){
                mod_base = base_complement(mod.base);
                idx = base_to_idx(mod_base);
            } else {
                mod_base = mod.base;
                idx = base_to_idx(mod_base);
            }
            
            // print_array(bases_pos_lens, 5, 'i');
            if(base_rank >= bases_pos_lens[idx]) {
                WARNING("%d th base of %c not found in SEQ. %c base count is %d read_id:%s seq_len:%d mod.base:%c mod_codes:%s\n", base_rank, mod_base, mod_base, bases_pos_lens[idx], qname, seq_len, mod.base, mod.mod_codes);
                continue;
            }
            ASSERT_MSG(base_rank < bases_pos_lens[idx], "%d th base of %c not found in SEQ. %c base count is %d read_id:%s seq_len:%d mod.base:%c mod_codes:%s\n", base_rank, mod_base, mod_base, bases_pos_lens[idx], qname, seq_len, mod.base, mod.mod_codes);
            
            if(rev) {
                read_pos = bases_pos[idx][bases_pos_lens[idx] - base_rank - 1];
                read_pos = seq_len - read_pos - 1;
            } else {
                read_pos = bases_pos[idx][base_rank];
            }

            ASSERT_MSG(read_pos < seq_len, "Base pos cannot exceed seq len. read_pos: %d seq_len: %d\n", read_pos, seq_len);
                        
            // mod prob per each mod code. TO-DO: need to change when code is ChEBI id
            for(int k=0; k<mod.mod_codes_len; k++) {
                // get the mod prob
                ml_idx = ml_start_idx + j*mod.mod_codes_len + k;
                ASSERT_MSG(ml_idx<ml_len, "ml_idx:%d ml_len:%d\n", ml_idx, ml_len);
                uint8_t mod_prob_scaled = ml[ml_idx];
                ASSERT_MSG(mod_prob_scaled <= 255 && mod_prob_scaled>=0, "mod_prob_scaled:%d\n", mod_prob_scaled);
                double mod_prob = (mod_prob_scaled+1)/256.0;

                // add to mods
                if(bases[read_pos].mods_len >= bases[read_pos].mods_cap){
                    bases[read_pos].mods_cap *= 2;
                    bases[read_pos].mods = (mod_t *)realloc(bases[read_pos].mods, sizeof(mod_t)*bases[read_pos].mods_cap);
                    MALLOC_CHK(bases[read_pos].mods);
                }

                bases[read_pos].mods[bases[read_pos].mods_len].base = &bases[read_pos];
                bases[read_pos].mods[bases[read_pos].mods_len].mod_code = mod.mod_codes[k];
                bases[read_pos].mods[bases[read_pos].mods_len].mod_strand = mod.strand;
                bases[read_pos].mods[bases[read_pos].mods_len].mod_prob = mod_prob;
                bases[read_pos].mods[bases[read_pos].mods_len].mod_base = mod.base;
                bases[read_pos].is_called[mod_code_to_idx(mod.mod_codes[k])] = 1;
                bases[read_pos].mods_len++;

                if(mod.status_flag=='.' && prev_read_pos != -1 && prev_read_pos+1 < read_pos){
                    for(int l=prev_read_pos+1;l<read_pos;l++){
                        bases[read_pos].is_skipped[mod_code_to_idx(mod.mod_codes[k])] = 1;
                    }
                }

            }

            prev_read_pos = read_pos;

        }
        ml_start_idx = ml_idx + 1;
    }

    // free base_pos
    for(int i=0;i<5;i++){
        free(bases_pos[i]);
    }

    return bases;

}

static stat_t ** get_stats(khash_t(str)* stats_map, uint32_t *meth_freqs_len){
    uint32_t len = 0;
    stat_t ** stats = (stat_t **)malloc(sizeof(stat_t *)*kh_size(stats_map));
    MALLOC_CHK(stats);
    for (khiter_t k = kh_begin(stats_map); k != kh_end(stats_map); ++k) {
        if (kh_exist(stats_map, k)) {
            stat_t * stat = kh_value(stats_map, k);
            stat->freq = (double)stat->n_mod/stat->n_called;
            stats[len] = stat;
            len++;
        }
    }
    *meth_freqs_len = len;
    return stats;
}


const char *get_mm_tag_ptr(bam1_t *record){

    const char* tag = "MM";
    // get the mm
    uint8_t *data = bam_aux_get(record, tag);
    if(data == NULL){
        WARNING("%s tag not found in read %s",tag, bam_get_qname(record));
        return NULL;
        //exit(EXIT_FAILURE);
    }

    const char *mm_str = bam_aux2Z(data);
    if(mm_str == NULL){
        WARNING("%s tag could not be decoded for %s. Is it type Z?",tag, bam_get_qname(record));
        exit(EXIT_FAILURE);
    }

    return mm_str;

}


uint8_t *get_ml_tag(bam1_t *record, uint32_t *len_ptr){

    const char* tag = "ML";
    // get the mm
    uint8_t *data = bam_aux_get(record, tag);
    if(data == NULL){
        WARNING("%s tag not found in read %s",tag, bam_get_qname(record));
        return NULL;
        //exit(EXIT_FAILURE);
    }

    // check if the type of the tag is array of integers
    const char aux_type = data[0];
    if (aux_type != 'B') {
        ERROR("%s tag is not of type B in read %s",tag, bam_get_qname(record));
        exit(EXIT_FAILURE);
    }

    // get the array len
    uint32_t len = bam_auxB_len(data);

    // check if the array type is uint8_t
    const char aux_array_type = data[1];
    if (aux_array_type != 'C') {
        ERROR("%s array tag type '%c' is not of type C in read %s",tag, aux_array_type, bam_get_qname(record));
        exit(EXIT_FAILURE);
    }

    // get the actual stuff
    uint8_t *array = (uint8_t *)malloc(sizeof(uint8_t)*len); //can be optimised for premalloced array as arg
    MALLOC_CHK(array);

    //fprintf(stderr,"%s: ",tag);
    for(int i=0;i<len;i++){
        array[i] = bam_auxB2i(data,i);
        //fprintf(stderr, "%d,", array[i]);
    }
    //fprintf(stderr, "\n");

    //set the length
    *len_ptr = len;

    return array;

}

static void print_meth_call_hdr(){
    printf("read_id\tread_pos\tref_pos\tchrom\tmod_strand\tref_strand\tref_mod_strand\tfw_soft_clipped_start\tfw_soft_clipped_end\tread_length\tmod_qual\tmod_code\tbase_qual\tref_kmer\tquery_kmer\tcanonical_base\tmodified_primary_base\tinferred\tflag\n");
}

static void print_meth_freq_hdr(FILE * output_file){
    fprintf(output_file, "chrom\tstart\tend\tdepth\tn_mod\tn_called\tn_skipped\tfreq\tmod_code\tstrand\n");
}

static void print_mods(base_t *bases, uint32_t seq_len, bam_hdr_t *hdr, bam1_t *record, enum MOD_CODES print_mod_code){

    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *qname = bam_get_qname(record);

    uint16_t flag = record->core.flag;


    for(int i=0;i<seq_len;i++){
        for(int j=0;j<bases[i].mods_len;j++){
            mod_t mod = bases[i].mods[j];
            base_t base = bases[i];
            if(print_mod_code !='*' && mod.mod_code != print_mod_code){
                continue;
            }
            fprintf(stdout, "%s\t%d\t%d\t%s\t%c\t%c\t%c\t%d\t%d\t%d\t%f\t%c\t%d\t%s\t%s\t%c\t%c\t%c\t%d\n", qname, i, base.ref_pos, base.chrom, mod.mod_strand, base.strand, '*', -1, -1, seq_len, mod.mod_prob, mod.mod_code, base.qual, "*", "*", '*', mod.mod_base, '*', flag);
        }
    }
}

static void print_meth_freq(FILE * output_file, stat_t ** stats, uint32_t seq_len, enum MOD_CODES print_mod_code){
    print_meth_freq_hdr(output_file);
    for(int i=0;i<seq_len;i++){
        stat_t * stat = stats[i];
        if(print_mod_code !='*' && stat->mod_code != print_mod_code){
            continue;
        }
        fprintf(output_file, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%c\t%c\n", stat->chrom, stat->start, stat->end, stat->depth, stat->n_mod, stat->n_called, stat->n_skipped, stat->freq, stat->mod_code, stat->strand);
    }

}

static void free_mod_tags(mod_tag_t *mod_tags, uint32_t len){
    for(int i=0;i<len;i++){
        free(mod_tags[i].mod_codes);
        free(mod_tags[i].skip_counts);
    }
    free(mod_tags);
}

static void free_bases(base_t *bases, uint32_t len){
    for(int i=0;i<len;i++){
        free(bases[i].mods);
    }
    free(bases);
}

static void free_stats_map(khash_t(str)* stats_map){
    //free keys
    for (khiter_t k = kh_begin(stats_map); k != kh_end(stats_map); ++k) {
        if (kh_exist(stats_map, k)) {
            stat_t * stat = kh_value(stats_map, k);
            free((char *) stat->chrom);
            free((char *)kh_key(stats_map, k));
            free(kh_value(stats_map, k));
        }
    }
    kh_destroy(str, stats_map);
}

void simple_meth_view(core_t* core){

    print_meth_call_hdr();
    
    bam1_t *record = bam_init1();
    while(sam_itr_next(core->bam_fp, core->itr, record) >= 0){

        uint32_t seq_len = record->core.l_qseq;

        if(seq_len==0){
            continue;
        }

        const char *mm = get_mm_tag_ptr(record);
        uint32_t ml_len;
        uint8_t *ml = get_ml_tag(record, &ml_len);

        if(ml == NULL || ml_len <= 0){
            free(ml);
            continue;
        }

        uint32_t mods_len = 0;
        mod_tag_t *mod_tags = extract_mods(mm, &mods_len);

        bam_hdr_t *hdr = core->bam_hdr;

        int * aln_pairs = get_aln(hdr, record);
        base_t * bases = get_bases(mod_tags, mods_len, ml, ml_len, aln_pairs, hdr, record);
        
        print_mods(bases, seq_len, hdr, record, '*');

        free_bases(bases, seq_len);
        free(aln_pairs);
        free_mod_tags(mod_tags, mods_len);
        free(ml);

    }
    
    bam_destroy1(record);
    return;
}

void meth_freq(core_t* core){

    bam1_t *record = bam_init1();
    while(sam_itr_next(core->bam_fp, core->itr, record) >= 0){

        uint32_t seq_len = record->core.l_qseq;

        if(seq_len==0){
            continue;
        }

        const char *mm = get_mm_tag_ptr(record);
        uint32_t ml_len;
        uint8_t *ml = get_ml_tag(record, &ml_len);

        if(ml == NULL || ml_len <= 0){
            free(ml);
            continue;
        }

        uint32_t mods_len = 0;
        mod_tag_t *mod_tags = extract_mods(mm, &mods_len);

        bam_hdr_t *hdr = core->bam_hdr;

        int * aln_pairs = get_aln(hdr, record);

        base_t * bases = get_bases(mod_tags, mods_len, ml, ml_len, aln_pairs, hdr, record);
        update_stats(bases, seq_len, stats_map);

        free_bases(bases, seq_len);
        free(aln_pairs);
        free_mod_tags(mod_tags, mods_len);
        free(ml);
    }

    bam_destroy1(record);
    return;
}

void init_mod(){
    stats_map = kh_init(str);
}

void destroy_mod(){
    free_stats_map(stats_map);
}

void print_stats(FILE * output_file){
    uint32_t meth_freqs_len = 0;
    stat_t ** stats = get_stats(stats_map, &meth_freqs_len);
    print_meth_freq(output_file, stats, meth_freqs_len, MOD_5mC);
    free(stats);
}