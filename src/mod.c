/**
 * @file pulse.c
 * @brief modification tags

MIT License

Copyright (c) 2023 Hasindu Gamaarachchi
Copyright (c) 2023 Suneth Samarasinghe

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
// #include "misc.h"
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
} modbase_t;

typedef struct {
    modbase_t * mod_bases;
    int mod_bases_cap;
    int mod_bases_len;
} mod_t; // free in free_mods_per_base()

typedef struct {
    int ref_pos;
    int read_pos;
    char base;
    char mod_strand;
    char mod_code;
    double mod_prob;
} meth_t;

enum MOD_CODES {
    MOD_ALL = 'A',
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
    int n_reads;
    int n_methylated;
    double freq;
} meth_freq_t;

KHASH_MAP_INIT_STR(str, meth_freq_t);

char* make_key(char *chrom, int start, int end){
    int start_strlen = snprintf(NULL, 0, "%d", start);
    int end_strlen = snprintf(NULL, 0, "%d", end);
    int key_strlen = strlen(chrom) + start_strlen + end_strlen + 3;
    char* key = (char *)malloc(key_strlen * sizeof(char));
    MALLOC_CHK(key);
    snprintf(key, key_strlen, "%s\t%d\t%d", chrom, start, end);
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
* @param mod_code code of the modification to extract. If MOD_ALL, all modifications are extracted
* @return pointer to the array of mod_tag_t
*/
static mod_tag_t *extract_mods(const char *mm_string, uint32_t *len_mods, enum MOD_CODES mod_code) {
    if (mm_string == NULL || strlen(mm_string) == 0) {
        ERROR("%s","Error: Empty MM string.\n");
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

            if(mod_code == MOD_ALL || mm_string[i] == mod_code) {
                current_mod.mod_codes[j] = mm_string[i];
                j++;
            }

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

static void print_meth_call_hdr(){
    printf("chrom\tref_pos\tread_name\tread_pos\tstrand\tbase\tmod_strand\tmod_code\tmod_prob\tflag\n");
}

static void print_meth_freq_hdr(){
    printf("chrom\tstart\tend\tn_reads\tn_methylated\tfreq\n");
}

static meth_t * call_meth(mod_t * mods, uint32_t prob_len, bam_hdr_t *hdr, bam1_t *record, uint32_t * meths_len){
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";
    int32_t pos = record->core.pos;
    int32_t end = bam_endpos(record);
    const char *qname = bam_get_qname(record);

    int8_t rev = bam_is_rev(record);
    const char strand = rev ? '-' : '+';
    ASSERT_MSG(!(record->core.flag & BAM_FUNMAP), "Unmapped read %s\n", qname);

    uint32_t len = 0;
    meth_t * meths = (meth_t *)malloc(sizeof(meth_t)*prob_len);
    MALLOC_CHK(meths);

    if(!rev){

        uint32_t *cigar = bam_get_cigar(record);
        uint32_t n_cigar = record->core.n_cigar;

        int read_pos = 0;
        int ref_pos = pos;

        //fprintf(stderr,"n cigar: %d\n", n_cigar);

        for (uint32_t ci = 0; ci < n_cigar; ++ci) {

            const uint32_t c = cigar[ci];
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
            } else {
                ERROR("Unhandled CIGAR OPT Cigar: %d\n", cigar_op);
                exit(EXIT_FAILURE);
            }

            // Iterate over the pairs of aligned bases
            for(int j = 0; j < cigar_len; ++j) {
                if(is_aligned) {
                    assert(read_pos < prob_len);
                    char base = seq_nt16_str[bam_seqi(bam_get_seq(record), read_pos)];
                    mod_t mod = mods[read_pos];
                    for(int k=0;k<mod.mod_bases_len;k++){
                        modbase_t mod_base = mod.mod_bases[k];
                        meths[len].ref_pos = ref_pos;
                        meths[len].read_pos = read_pos;
                        meths[len].base = base;
                        meths[len].mod_strand = mod_base.mod_strand;
                        meths[len].mod_code = mod_base.mod_code;
                        meths[len].mod_prob = mod_base.mod_prob;
                        len++;
                        // fprintf(stdout, "%s\t%d\t%s\t%d\t%c\t%c\t%c\t%c\t%f\n", tname, ref_pos, qname, read_pos, strand, base, mod_base.mod_strand, mod_base.mod_code, mod_base.mod_prob);
                    }
                }

                // increment
                read_pos += read_inc;
                ref_pos += ref_inc;
            }
        }
        
    } 

    *meths_len = len;

    ASSERT_MSG(len <= prob_len, "len:%d prob_len:%d\n", len, prob_len);

    return meths;
}

static mod_t * get_mods_per_base(mod_tag_t *mod_tags, uint32_t mods_len, uint8_t * ml, uint32_t ml_len, bam_hdr_t *hdr, bam1_t *record){

    const char *qname = bam_get_qname(record);
    int8_t rev = bam_is_rev(record);

    uint8_t *seq = bam_get_seq(record);
    uint32_t seq_len = record->core.l_qseq;

    if(seq_len == 0){
        ERROR("Sequence length is 0 for read %s. Found %d mod_tags", qname, mods_len);
        exit(EXIT_FAILURE);
    }

    mod_t *mods = (mod_t *)malloc(sizeof(mod_t)*seq_len);
    MALLOC_CHK(mods);
    for(int i=0;i<seq_len;i++){
        mods[i].mod_bases_cap = INIT_MOD_BASES;
        mods[i].mod_bases = (modbase_t *)malloc(sizeof(modbase_t)*mods[i].mod_bases_cap);
        MALLOC_CHK(mods[i].mod_bases);
        mods[i].mod_bases_len = 0;
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
    
    int ml_idx = 0;
    // go through mod_tags    
    for(int i=0; i<mods_len; i++) {
        mod_tag_t mod = mod_tags[i];
        int base_rank = -1;


        for(int j=0; j<mod.skip_counts_len; j++) {
            base_rank += mod.skip_counts[j] + 1;
            // fprintf(stderr, "qname:%s seq_len:%d mod.base:%c seq.strand:%c mod.strand:%c rank:%d\n", qname, seq_len, mod.base, strand, mod.strand, base_rank);
            
            char mod_base;
            int idx;
            int base_pos;

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
                base_pos = bases_pos[idx][bases_pos_lens[idx] - base_rank - 1];
                base_pos = seq_len - base_pos - 1;
            } else {
                base_pos = bases_pos[idx][base_rank];
            }
            
            ASSERT_MSG(base_pos < seq_len, "Base pos cannot exceed seq len. base_pos: %d seq_len: %d\n", base_pos, seq_len); 
            
            ml_idx += j*mod.mod_codes_len;
            // mod prob per each mod code. TO-DO: need to change when code is ChEBI id
            for(int k=0; k<mod.mod_codes_len; k++) {
                // get the mod prob
                ml_idx += k;
                uint8_t mod_prob_scaled = ml[j*mod.mod_codes_len + k];
                double mod_prob = (double)(mod_prob_scaled+1)/256.0;
                
                // add to mods
                if(mods[base_pos].mod_bases_len >= mods[base_pos].mod_bases_cap){
                    mods[base_pos].mod_bases_cap *= 2;
                    mods[base_pos].mod_bases = (modbase_t *)realloc(mods[base_pos].mod_bases, sizeof(modbase_t)*mods[base_pos].mod_bases_cap);
                    MALLOC_CHK(mods[base_pos].mod_bases);
                }

                mods[base_pos].mod_bases[mods[base_pos].mod_bases_len].mod_code = mod.mod_codes[k];
                mods[base_pos].mod_bases[mods[base_pos].mod_bases_len].mod_strand = mod.strand;
                mods[base_pos].mod_bases[mods[base_pos].mod_bases_len].mod_prob = mod_prob;
                mods[base_pos].mod_bases_len++;

            }

        }
        ml_idx++;
    }

    // free base_pos
    for(int i=0;i<5;i++){
        free(bases_pos[i]);
    }

    return mods;

}

static void update_meth_freq(khash_t(str) *meth_freqs, meth_freq_t meth_freq){
    char *key = make_key(meth_freq.chrom, meth_freq.start, meth_freq.end);
    khiter_t k = kh_get(str, meth_freqs, key);
    if (k == kh_end(meth_freqs)) {
        int ret;
        k = kh_put(str, meth_freqs, key, &ret);
        kh_value(meth_freqs, k) = meth_freq;
    } else {
        meth_freq_t *mf = &kh_value(meth_freqs, k);
        mf->n_reads += meth_freq.n_reads;
        mf->n_methylated += meth_freq.n_methylated;
        mf->freq = (double)mf->n_methylated/mf->n_reads;
    }
}

static void call_freqs(meth_t * meths, uint32_t meths_len, khash_t(str)* freq_map,  bam_hdr_t *hdr, bam1_t *record){
    uint32_t len = 0;

    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";
    const char *qname = bam_get_qname(record);

    for(int i=0;i<meths_len;i++) {
        meth_t meth = meths[i];
        meth_freq_t meth_freq;
        meth_freq.chrom = tname;
        meth_freq.start = meth.ref_pos;
        meth_freq.end = meth.ref_pos;
        meth_freq.n_reads = 1; // revist this
        meth_freq.n_methylated = meth.mod_prob > 0.0 ? 1 : 0; // revist this
        meth_freq.freq = 0.0;
        update_meth_freq(freq_map, meth_freq);
    }

}

static meth_freq_t * get_meth_freqs(khash_t(str)* freq_map, uint32_t *meth_freqs_len){
    uint32_t len = 0;
    meth_freq_t * meth_freqs = (meth_freq_t *)malloc(sizeof(meth_freq_t)*kh_size(freq_map));
    MALLOC_CHK(meth_freqs);
    for (khiter_t k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
        if (kh_exist(freq_map, k)) {
            meth_freqs[len] = kh_value(freq_map, k);
            len++;
        }
    }
    *meth_freqs_len = len;
    return meth_freqs;
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


static void print_ml_array(const uint8_t *array, uint32_t len, bam1_t *record){

    fprintf(stdout, "%s\t%d\t%s\t", bam_get_qname(record),record->core.pos, "ML");
    for(int i=0;i<len;i++){
        fprintf(stdout, "%d,", array[i]);
    }
    fprintf(stdout, "\n");

}

static void print_mm_array(const char *mm, uint32_t len, bam1_t *record){

    fprintf(stdout, "%s\t%d\t%s\t", bam_get_qname(record),record->core.pos, "MM");
    for(int i=0;i<len;i++){
        fprintf(stdout, "%c", mm[i]);
    }
    fprintf(stdout, "\n");

}

static void print_meths(meth_t *meths, uint32_t len, bam_hdr_t *hdr, bam1_t *record){

    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";
    const char *qname = bam_get_qname(record);

    uint16_t flag = record->core.flag;

    int8_t rev = bam_is_rev(record);
    const char strand = rev ? '-' : '+';

    for(int i=0;i<len;i++){
        fprintf(stdout, "%s\t%d\t%s\t%d\t%c\t%c\t%c\t%c\t%f\t%d\n", tname, meths[i].ref_pos, qname, meths[i].read_pos, strand, meths[i].base, meths[i].mod_strand, meths[i].mod_code, meths[i].mod_prob, flag);
    }
}

static void print_meth_freq(meth_freq_t * meth_freqs, uint32_t len){

    for(int i=0;i<len;i++){
        fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%f\n", meth_freqs[i].chrom, meth_freqs[i].start, meth_freqs[i].end, meth_freqs[i].n_reads, meth_freqs[i].n_methylated, meth_freqs[i].freq);
    }

}

static void free_mod_tags(mod_tag_t *mod_tags, uint32_t len){
    for(int i=0;i<len;i++){
        free(mod_tags[i].mod_codes);
        free(mod_tags[i].skip_counts);
    }
    free(mod_tags);
}

static void free_mods_per_base(mod_t *mods, uint32_t len){
    for(int i=0;i<len;i++){
        free(mods[i].mod_bases);
    }
    free(mods);
}

void simple_meth_view(core_t* core){

    print_meth_call_hdr();

    bam1_t *record = bam_init1();
    while(sam_itr_next(core->bam_fp, core->itr, record) >= 0){


        const char *mm = get_mm_tag_ptr(record);
        uint32_t ml_len;
        uint8_t *ml = get_ml_tag(record, &ml_len);

        uint32_t mods_len = 0;
        mod_tag_t *mod_tags = extract_mods(mm, &mods_len, MOD_ALL);

        bam_hdr_t *hdr = core->bam_hdr;

        if(ml==NULL){
            continue;
        }

        if(ml_len==0){ // no mod probs, therefore no mods
            ASSERT_MSG(mods_len==0, "Mods len should be 0 when ml_len is 0. mods_len:%d ml_len:%d\n", mods_len, ml_len);
        }

        mod_t * mods_per_base = get_mods_per_base(mod_tags, mods_len, ml, ml_len, hdr, record);
        uint32_t mods_per_base_len = record->core.l_qseq;

        uint32_t meths_len = 0;
        meth_t * meths = call_meth(mods_per_base, mods_per_base_len, hdr, record, &meths_len);
        
        // print_ml_array(ml, len, record);
        // print_mm_array(mm, len, record);
        // print_mods(mod_tags, mods_len, hdr, record);

        print_meths(meths, meths_len, hdr, record);
        
        free(ml);
        free_mod_tags(mod_tags, mods_len);
        free_mods_per_base(mods_per_base, mods_per_base_len);
        free(meths);


    }
    bam_destroy1(record);
    return;
}

void meth_freq(core_t* core){

    print_meth_freq_hdr();
    khash_t(str)* freq_map = kh_init(str);

    bam1_t *record = bam_init1();
    while(sam_itr_next(core->bam_fp, core->itr, record) >= 0){


        const char *mm = get_mm_tag_ptr(record);
        uint32_t ml_len;
        uint8_t *ml = get_ml_tag(record, &ml_len);

        uint32_t mods_len = 0;
        mod_tag_t *mod_tags = extract_mods(mm, &mods_len, MOD_ALL);

        bam_hdr_t *hdr = core->bam_hdr;

        if(ml==NULL){
            continue;
        }

        if(ml_len==0){ // no mod probs, therefore no mods
            ASSERT_MSG(mods_len==0, "Mods len should be 0 when ml_len is 0. mods_len:%d ml_len:%d\n", mods_len, ml_len);
        }

        mod_t * mods_per_base = get_mods_per_base(mod_tags, mods_len, ml, ml_len, hdr, record);
        uint32_t mods_per_base_len = record->core.l_qseq;

        uint32_t meths_len = 0;
        meth_t * meths = call_meth(mods_per_base, mods_per_base_len, hdr, record, &meths_len);

        call_freqs(meths, meths_len, freq_map, hdr, record);
        
        free(ml);
        free_mod_tags(mod_tags, mods_len);
        free_mods_per_base(mods_per_base, mods_per_base_len);
        free(meths);

    }

    uint32_t meth_freqs_len = 0;
    meth_freq_t * meth_freqs = get_meth_freqs(freq_map, &meth_freqs_len);

    print_meth_freq(meth_freqs, meth_freqs_len);

    kh_destroy(str, freq_map);

    bam_destroy1(record);
    return;
}


