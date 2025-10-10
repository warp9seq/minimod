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
#include "ref.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <stdbool.h>
#include <pthread.h>

#define IS_ALPHA(c) (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')
#define IS_DIGIT(c) (c >= '0' && c <= '9')

#define KEY_SIZE 4

static const int valid_bases[256] = { ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1, ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1, ['u'] = 1, ['n'] = 1 };
static const int valid_strands[256] = { ['+'] = 1, ['-'] = 1 };
static const int base_idx_lookup[256] = { ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['U'] = 4, ['N'] = 5, ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3, ['u'] = 4, ['n'] = 5 };
static const char base_complement_lookup[256] = { ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A', ['N'] = 'N', ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a', ['n'] = 'n' };
static const char* default_context[256] = { ['*'] = "*", ['m'] = "CG", ['h'] = "CG", ['f'] = "C", ['c'] = "C", ['C'] = "C", ['g'] = "T", ['e'] = "T", ['b'] = "T", ['T'] = "T", ['U'] = "U", ['a'] = "A", ['A'] = "A", ['o'] = "G", ['G'] = "G", ['n'] = "N", ['N'] = "N" };

static inline char* get_default_context(char* mod_code) {
    if(strlen(mod_code) == 1) {
        char c = mod_code[0];
        if (default_context[(int)c] != NULL) {
            return (char*)default_context[(int)c];
        }
    }

    return "CG"; // default context for unknown modification codes
}

static inline int die(const char *format, ...) {
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(EXIT_FAILURE);
}

const char *get_mm_tag_ptr(bam1_t *record){

    const char* tag = "MM";
    // get the mm
    uint8_t *data = bam_aux_get(record, tag);
    if(data == NULL){
        LOG_TRACE("%s tag not found in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    const char *mm_str = bam_aux2Z(data);
    if(mm_str == NULL){
        LOG_TRACE("%s tag could not be decoded for %s. Is it type Z?",tag, bam_get_qname(record));
        return NULL;
    }

    return mm_str;
}

uint8_t *get_ml_tag(bam1_t *record, uint32_t *len_ptr){

    const char* tag = "ML";
    // get the mm
    uint8_t *data = bam_aux_get(record, tag);
    if(data == NULL){
        LOG_TRACE("%s tag not found in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    // check if the type of the tag is array of integers
    const char aux_type = data[0];
    if (aux_type != 'B') {
        LOG_TRACE("%s tag is not of type B in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    // get the array len
    uint32_t len = bam_auxB_len(data);
    if(len == 0){
        LOG_TRACE("%s tag array length is 0 in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    // check if the array type is uint8_t
    const char aux_array_type = data[1];
    if (aux_array_type != 'C') {
        LOG_TRACE("%s array tag type '%c' is not of type C in read %s",tag, aux_array_type, bam_get_qname(record));
        return NULL;
    }

    // get the actual stuff
    uint8_t *array = (uint8_t *)malloc(sizeof(uint8_t)*len); //can be optimised for premalloced array as arg
    MALLOC_CHK(array);

    for(int i=0;i<len;i++){
        array[i] = bam_auxB2i(data,i);
    }

    //set the length
    *len_ptr = len;

    return array;
}

// get the haplotype integer from the HP tag
uint8_t get_hp_tag(bam1_t *record){
    
        const char* tag = "HP";
        // get the HP tag
        uint8_t *data = bam_aux_get(record, tag);
        if(data == NULL){
            LOG_TRACE("%s tag not found in read %s",tag, bam_get_qname(record));
            return 0;
        }
    
        // get the integer value
        uint8_t hp = (uint8_t) bam_aux2i(data);
    
        return hp;
}

void parse_mod_codes(opt_t *opt) {

    uint8_t n_codes = 0;
    int i=0;

    i=0;
    while(1){
        if(opt->mod_codes_str[i] == '\0'){ // end of string
            break;
        }
        int has_nums = 0;
        int has_alpha = 0;

        int mod_code_cap = 2;
        char * mod_code = (char *)malloc(sizeof(char) * mod_code_cap);
        MALLOC_CHK(mod_code);
        
        int j = 0;
        while(opt->mod_codes_str[i] != ',' && opt->mod_codes_str[i] != '[' && opt->mod_codes_str[i] != '\0'){
            if(IS_ALPHA(opt->mod_codes_str[i]) || opt->mod_codes_str[i] == '*'){
                has_alpha = 1;
            } else if (IS_DIGIT(opt->mod_codes_str[i])){
                has_nums = 1;
            } else {
                ERROR("Invalid character %c in modification code in -c argument", opt->mod_codes_str[i]);
                exit(EXIT_FAILURE);
            }

            if(j >= mod_code_cap){
                mod_code_cap *= 2;
                mod_code = (char *)realloc(mod_code, sizeof(char) * mod_code_cap);
                MALLOC_CHK(mod_code);
            }

            mod_code[j] = opt->mod_codes_str[i];
            j++;
            i++;
        }
        mod_code[j] = '\0';

        if(has_alpha && has_nums){
            ERROR("Modification code %s cannot contain both letters and numbers in -c argument", mod_code);
            exit(EXIT_FAILURE);
        }
        
        int context_cap = 3;
        char * context = (char *)malloc(context_cap * sizeof(char) + 1);
        MALLOC_CHK(context);

        if(opt->mod_codes_str[i] == '['){ // context given
            i++;
            int j = 0;
            int is_star = 0;
            while(opt->mod_codes_str[i] != ']' || opt->mod_codes_str[i] == '\0'){
                if(opt->mod_codes_str[i] == '*'){
                    is_star = 1;
                } else if(valid_bases[(int)opt->mod_codes_str[i]]==0){
                    ERROR("Invalid character %c in context for modification code %s in -c argument", opt->mod_codes_str[i], mod_code);
                    exit(EXIT_FAILURE);
                }
                if(j >= context_cap){
                    context_cap *= 2;
                    context = (char *)realloc(context, context_cap * sizeof(char) + 1);
                    MALLOC_CHK(context);
                }
                context[j] = toupper(opt->mod_codes_str[i]);
                i++;
                j++;
            }
            if(opt->mod_codes_str[i] == '\0'){
                ERROR("Context not closed with a ] for modification code %s in -c argument", mod_code);
                exit(EXIT_FAILURE);
            }
            context[j] = '\0';
            if(is_star && j > 1){
                ERROR("Invalid context for modification code %s. * should be the only character within [ and ] in -c argument", mod_code);
                exit(EXIT_FAILURE);
            }

            i++;
            if(opt->mod_codes_str[i] == ','){ // more modification codes
                i++;
            }
        } else if(opt->mod_codes_str[i] == ','){ // context is *
            char * default_ctx = get_default_context(mod_code);
            strcpy(context, default_ctx);
            INFO("Context not provided for modification code %s. Using %s in -c argument", mod_code, context);
            
            i++;
        } else if(opt->mod_codes_str[i] == '\0'){
            char * default_ctx = get_default_context(mod_code);
            strcpy(context, default_ctx);
            INFO("Context not provided for modification code %s in -c argument. Using %s", mod_code, context);
        } else {
            ERROR("Invalid character %c after modification code %s in -c argument", opt->mod_codes_str[i], mod_code);
            exit(EXIT_FAILURE);
        }

        modcodem_t * map_entry = (modcodem_t *)malloc(sizeof(modcodem_t));
        MALLOC_CHK(map_entry);
        map_entry->context = context;
        map_entry->index = n_codes;

        int ret;
        khint_t k = kh_put(modcodesm, opt->modcodes_map, mod_code, &ret);
        if (ret == -1) {
            ERROR("Failed to insert modification code %s into map", mod_code);
            exit(EXIT_FAILURE);
        }
        if (ret == 0) {
            ERROR("Duplicate modification code %s found in -c argument", mod_code);
            exit(EXIT_FAILURE);
        }
        kh_value(opt->modcodes_map, k) = map_entry;

        n_codes++;
    }

    opt->n_mods = n_codes;
}

void parse_mod_threshes(opt_t * opt) {
    int i=0;
    int n_thresh = 0;
    int thresh_str_cap = 1;
    double d = 0.0;

    while(opt->mod_threshes_str[i] != '\0'){
        char * thresh_str = (char *)malloc(thresh_str_cap * sizeof(char) + 1);
        MALLOC_CHK(thresh_str);
        int j = 0;
        while(opt->mod_threshes_str[i] != ',' && opt->mod_threshes_str[i] != '\0'){
            if(j >= thresh_str_cap){
                thresh_str_cap *= 2;
                thresh_str = (char *)realloc(thresh_str, thresh_str_cap * sizeof(char) + 1);
                MALLOC_CHK(thresh_str);
            }
            thresh_str[j] = opt->mod_threshes_str[i];
            i++;
            j++;
        }
        thresh_str[j] = '\0';

        errno = 0;
        d = atof(thresh_str);

        if(errno != 0){
            ERROR("Invalid threshold. You entered %s",thresh_str);
            exit(EXIT_FAILURE);
        }
        
        if(d<0 || d>1){
            ERROR("Modification threshold should be in the range 0.0 to 1.0. You entered %f",d);
            exit(EXIT_FAILURE);
        }

        khint_t m;
        
        for(m = kh_begin(opt->modcodes_map); m < kh_end(opt->modcodes_map); ++m) {
            if (kh_exist(opt->modcodes_map, m)) {
                char * key = (char *) kh_key(opt->modcodes_map, m);
                modcodem_t *mod_code_entry = kh_value(opt->modcodes_map, m);
                if(mod_code_entry->index != n_thresh){
                    continue; // skip if the index does not match
                }
                INFO("Modification code: %s, Context: %s, Threshold: %f", key, mod_code_entry->context, d);
                mod_code_entry->thresh = (uint8_t)(d * 255); // convert to 0-255 range
            }
        }

        free(thresh_str);
        n_thresh++;
        if(opt->mod_threshes_str[i] == '\0'){
            break;
        }
        i++;
    }

    if(n_thresh==1){
        khint_t i;
        for(i=kh_begin(opt->modcodes_map); i < kh_end(opt->modcodes_map); ++i) {
            if (!kh_exist(opt->modcodes_map, i)) continue;
            modcodem_t *mod_code_map = kh_value(opt->modcodes_map, i);
            char * mod_code = (char *) kh_key(opt->modcodes_map, i);
            mod_code_map->thresh = (uint8_t)(d * 255); // set the same threshold for all codes
            INFO("Modification code: %s, Context: %s, Threshold: %f", mod_code, mod_code_map->context, d);
        }
    } else if(n_thresh != opt->n_mods){
        ERROR("Number of modification codes and thresholds do not match. Codes:%d, Thresholds:%d",opt->n_mods,n_thresh);
        exit(EXIT_FAILURE);
    }
}

char* make_key(const char *chrom, int pos, uint16_t ins_offset, char * mod_code, char strand, int haplotype){
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

void decode_key(char *key, char **chrom, int *pos, uint16_t * ins_offset, char **mod_code, char *strand, int *haplotype){
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

void print_view_options(opt_t *opt) {
    khint_t i;
    for(i=kh_begin(opt->modcodes_map); i < kh_end(opt->modcodes_map); ++i) {
        if (!kh_exist(opt->modcodes_map, i)) continue;
        modcodem_t *mod_code_map = kh_value(opt->modcodes_map, i);
        char * mod_code = (char *) kh_key(opt->modcodes_map, i);
        INFO("Modification code: %s, Context: %s", mod_code, mod_code_map->context);
    }
}

void print_view_header(core_t* core) {
    char * common = "ref_contig\tref_pos\tstrand\tread_id\tread_pos\tmod_code\tmod_prob";
    char * ins_offset = "";
    char * haplotype = "";
    if(core->opt.insertions){
        ins_offset = "\tins_offset";
    }
    if(core->opt.haplotypes){
        haplotype = "\thaplotype";
    }

    fprintf(core->opt.output_fp, "%s%s%s\n", common, ins_offset, haplotype);
}

void print_view_output(core_t* core, db_t* db) {
    
    for(int i =0; i < db->n_bam_recs; i++) {
        bam1_t *record = db->bam_recs[i];
        const char *qname = bam_get_qname(record);
        khash_t(viewm) *view_map = db->view_maps[i];

        khint_t k;
        for (k = kh_begin(view_map); k != kh_end(view_map); ++k) {
            if (kh_exist(view_map, k)) {
                view_t* view = kh_value(view_map, k);
                char *tname = NULL;
                int ref_pos;
                uint16_t ins_offset;
                char *mod_code;
                char strand;
                int haplotype;
                char * key = (char *) kh_key(view_map, k);
                decode_key(key, &tname, &ref_pos, &ins_offset, &mod_code, &strand, &haplotype);

                fprintf(core->opt.output_fp, "%s\t%d\t%c\t%s\t%d\t%s\t%f", tname, ref_pos, strand, qname, view->read_pos, mod_code, view->mod_prob/255.0);
                if(core->opt.insertions==1){
                    fprintf(core->opt.output_fp, "\t%d", db->ins_offset[i][view->read_pos]);
                }
                if(core->opt.haplotypes==1){
                    fprintf(core->opt.output_fp, "\t%d", haplotype);
                }
                fprintf(core->opt.output_fp, "\n");

                free(tname);
                free(mod_code);
            }
        }
    }

    

    if(core->opt.output_fp != stdout){
        fclose(core->opt.output_fp);
    }
}

void print_freq_header(core_t * core) {
    if(!core->opt.bedmethyl_out) { // tsv output header, no header for bedmethyl
        char * common = "contig\tstart\tend\tstrand\tn_called\tn_mod\tfreq\tmod_code";
        char * ins_offset = "";
        char * haplotype = "";
        if(core->opt.insertions){
            ins_offset = "\tins_offset";
        }
        if(core->opt.haplotypes){
            haplotype = "\thaplotype";
        }

        fprintf(core->opt.output_fp, "%s%s%s\n", common, ins_offset, haplotype);
    }
}

/* Split tab-delimited keys for sorting*/
char** split_key(char* key, int size) {
    char** tok = (char**)malloc(sizeof(char*) * size);
    MALLOC_CHK(tok);
    char* cpy = (char*)malloc((strlen(key) + 1) * sizeof(char));
    MALLOC_CHK(cpy);
    strcpy(cpy, key);
    char* cpy_start = cpy;

    char* t = strtok(cpy, "\t");
    if (t) {
        char * tmp = (char *)malloc((strlen(t) + 1) * sizeof(char));
        MALLOC_CHK(tmp);
        strcpy(tmp, t);
        tok[0] = tmp;
    }

    for (int i = 1; i < size; i++) {
        char* t = strtok(NULL, "\t");
        if (t) {
            char * tmp = (char *)malloc((strlen(t) + 1) * sizeof(char));
            MALLOC_CHK(tmp);
            strcpy(tmp, t);
            tok[i] = tmp;
        }
    }

    free(cpy_start);

    return tok;
}

/* Compare two keys for sorting*/
int cmp_key(const void* a, const void* b) {
    char* key_a = *(char**)a;
    char* key_b = *(char**)b;

    char** toks_a = split_key(key_a, KEY_SIZE);
    char** toks_b = split_key(key_b, KEY_SIZE);

    int chromosome_neq = strcmp(toks_a[0], toks_b[0]);

    if (chromosome_neq) {
        for (int i = 0; i < KEY_SIZE; i++) {
            free(toks_a[i]);
            free(toks_b[i]);
        }
        free(toks_a);
        free(toks_b);
        return chromosome_neq;
    }

    int start_a = atoi(toks_a[1]);
    int start_b = atoi(toks_b[1]);
    int end_a = atoi(toks_a[2]);
    int end_b = atoi(toks_b[2]);

    for (int i = 0; i < KEY_SIZE; i++) {
        free(toks_a[i]);
        free(toks_b[i]);
    }
    free(toks_a);
    free(toks_b);

    if (start_a == start_b) {
        return end_a - end_b;
    }

    return start_a - start_b;
}

void print_freq_output(core_t * core) {
    khash_t(freqm) *freq_map = core->freq_map;

    // sort the keys
    char ** sorted_keys = (char **)malloc(sizeof(char*) * kh_size(freq_map));
    MALLOC_CHK(sorted_keys);
    int size = 0;
    for (khint_t k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
        if (kh_exist(freq_map, k)) {
            sorted_keys[size++] = (char *)kh_key(freq_map, k);
        }
    }
    qsort(sorted_keys, size, sizeof(char*), cmp_key);

    if(core->opt.bedmethyl_out) {
        // chrom, start, end, mod_code, n_called, strand, start, end, "255,0,0",  n_called, freq
        khint_t k;
        for (int i = 0; i < size; i++) {
            k = kh_get(freqm, freq_map, sorted_keys[i]);
            if (kh_exist(freq_map, k)) {
                freq_t* freq = kh_value(freq_map, k);
                double freq_value = (double)freq->n_mod*100/freq->n_called;
                char *contig = NULL;
                int ref_pos;
                uint16_t ins_offset;
                char *mod_code;
                char strand;
                int haplotype;
                char * key = (char *) kh_key(freq_map, k);
                decode_key(key, &contig, &ref_pos, &ins_offset, &mod_code, &strand, &haplotype);
                int end = ref_pos+1;
                fprintf(core->opt.output_fp, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t255,0,0\t%d\t%f\n", contig, ref_pos, end, mod_code, freq->n_called, strand, ref_pos, end, freq->n_called, freq_value);
                free(contig);
                free(mod_code);
            }
        }
        
    } else {
        // contig, start, end, strand, n_called, n_mod, freq, mod_code
        khint_t k;
        for (int i = 0; i < size; i++) {
            k = kh_get(freqm, freq_map, sorted_keys[i]);
            if (kh_exist(freq_map, k)) {
                freq_t* freq = kh_value(freq_map, k);
                double freq_value = (double)freq->n_mod/freq->n_called;
                char * contig = NULL;
                int ref_pos;
                uint16_t ins_offset;
                char * mod_code;
                char strand;
                int haplotype;
                char * key = (char *) kh_key(freq_map, k);
                decode_key(key, &contig, &ref_pos, &ins_offset, &mod_code, &strand, &haplotype);

                fprintf(core->opt.output_fp, "%s\t%d\t%d\t%c\t%d\t%d\t%f\t%s", contig, ref_pos, ref_pos, strand, freq->n_called, freq->n_mod, freq_value, mod_code);

                if(core->opt.insertions){
                    fprintf(core->opt.output_fp, "\t%d", ins_offset);
                } 
                if(core->opt.haplotypes) {
                    if(haplotype == -1){
                        fprintf(core->opt.output_fp, "\t*");
                    } else {
                        fprintf(core->opt.output_fp, "\t%d", haplotype);
                    }
                }
                fprintf(core->opt.output_fp, "\n");
                free(contig);
                free(mod_code);
            }
        }
    }

    if(core->opt.output_fp != stdout){
        fclose(core->opt.output_fp);
    }
    free(sorted_keys);
}

void destroy_freq_map(khash_t(freqm)* freq_map){
    khint_t k;
    for (k = kh_begin(freq_map); k != kh_end(freq_map); k++) {
        if (kh_exist(freq_map, k)) {
            char *key = (char *) kh_key(freq_map, k);
            freq_t* freq = kh_value(freq_map, k);
            free(freq);
            free(key);
        }
    }
    kh_destroy(freqm, freq_map);
}

void merge_freq_maps(core_t* core, db_t* db) {
    khash_t(freqm) *core_map = core->freq_map;
    for(int i=0;i<db->n_bam_recs;i++){
        khash_t(freqm) *rec_map = db->freq_maps[i];
        khint_t k;
        for (k = kh_begin(rec_map); k != kh_end(rec_map); ++k) {
            if (kh_exist(rec_map, k)) {
                freq_t *db_freq = kh_value(rec_map, k);
                char *key = (char *) kh_key(rec_map, k);
                freq_t *core_freq;
                khint_t core_k = kh_get(freqm, core_map, key);
                if (core_k == kh_end(core_map)) {
                    // key not found in core map, insert it
                    int ret;
                    core_k = kh_put(freqm, core_map, key, &ret);
                    core_freq = (freq_t *)malloc(sizeof(freq_t));
                    MALLOC_CHK(core_freq);
                    kh_value(core_map, core_k) = core_freq;
                    core_freq->n_called = db_freq->n_called;
                    core_freq->n_mod = db_freq->n_mod;
                } else {
                    // key found in core map
                    core_freq = kh_value(core_map, core_k);
                    core_freq->n_called += db_freq->n_called;
                    core_freq->n_mod += db_freq->n_mod;
                    free(key); // free the key as it is already copied in core_freq
                }
            }
        }
    }
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

    if(core->opt.insertions){
        for(int i=0;i<seq_len;i++){
            db->ins[bam_i][i] = -1;
            db->ins_offset[bam_i][i] = 0;
        }
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
        } else if(cigar_op == BAM_CHARD_CLIP) {
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

            if(core->opt.insertions && is_inserted) {
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

static void update_freq_map(khash_t(freqm) *freq_map, const char *tname, int ref_pos, int ins_offset, char *mod_code, char strand, int haplotype, int is_called, int is_mod) {
    char * key = make_key(tname, ref_pos, ins_offset, mod_code, strand, haplotype);
    khiter_t k = kh_get(freqm, freq_map, key);
    if (k == kh_end(freq_map)) { // not found, add
        freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
        MALLOC_CHK(freq);
        freq->n_called = is_called;
        freq->n_mod = is_mod;
        int ret;
        k = kh_put(freqm, freq_map, key, &ret);
        kh_value(freq_map, k) = freq;
    } else { // found, update
        free(key);
        freq_t * freq = kh_value(freq_map, k);
        freq->n_called += is_called;
        freq->n_mod += is_mod;
        // check if freq->n_called overflows
        if(freq->n_called == 0){
            ERROR("n_called overflowed for key %s. Please report this issue.", key);
            exit(EXIT_FAILURE);
        }
    }

    if(haplotype != -1) {
        char * key = make_key(tname, ref_pos, ins_offset, mod_code, strand, -1); // aggregate all haplotypes
        khiter_t k = kh_get(freqm, freq_map, key);
        if (k == kh_end(freq_map)) { // not found, add
            freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
            MALLOC_CHK(freq);
            freq->n_called = is_called;
            freq->n_mod = is_mod;
            int ret;
            k = kh_put(freqm, freq_map, key, &ret);
            kh_value(freq_map, k) = freq;
        } else { // found, update
            free(key);
            freq_t * freq = kh_value(freq_map, k);
            freq->n_called += is_called;
            freq->n_mod += is_mod;
            // check if freq->n_called overflows
            if(freq->n_called == 0){
                ERROR("n_called overflowed for key %s. Please report this issue.", key);
                exit(EXIT_FAILURE);
            }
        }
    }
}

static void add_view_entry(khash_t(viewm) *view_map, const char *tname, int ref_pos, int ins_offset, char *mod_code, char strand, int haplotype, uint8_t mod_prob, int read_pos) {

    char *key = make_key(tname, ref_pos, ins_offset, mod_code, strand, haplotype);
    khiter_t k = kh_get(viewm, view_map, key);
    if (k == kh_end(view_map)) { // not found, add
        view_t *view = (view_t *)malloc(sizeof(view_t));
        MALLOC_CHK(view);
        view->mod_prob = mod_prob;
        view->read_pos = read_pos;
        int ret;
        k = kh_put(viewm, view_map, key, &ret);
        kh_value(view_map, k) = view;
    } else { // found, update
        free(key);
        ERROR("Duplicate entry found for key %s. This should not happen. Please report this issue.", kh_key(view_map, k));
        exit(EXIT_FAILURE);
    }
}

void freq_view_single(core_t * core, db_t *db, int32_t bam_i) {
    bam1_t *record = db->bam_recs[bam_i];
    const char *qname = bam_get_qname(record);
    int8_t rev = bam_is_rev(record);
    bam_hdr_t *hdr = core->bam_hdr;
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = (tid >= 0) ? hdr->target_name[tid] : "*";
    uint8_t *seq = bam_get_seq(record);
    uint32_t seq_len = record->core.l_qseq;
    char strand = rev ? '-' : '+';
    ref_t *ref = get_ref(tname);
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
    // char mod_strand;  // commented for now. might need to revisit
    char * mod_codes = db->mod_codes[bam_i];
    int mod_codes_len;
    int * skip_counts = db->skip_counts[bam_i];
    int skip_counts_len;
    // char status_flag;

    while (i < mm_str_len) {
        // reset skip counts and mod codes
        skip_counts_len = 0;
        mod_codes_len = 0;

        // set default status flag to '.' (when not present or '.' in the MM string)
        // status_flag = '.';

        // get base
        if(i < mm_str_len) {
            ASSERT_MSG(valid_bases[(int)mm_string[i]], "Invalid base:%c\n", mm_string[i]);
            modbase = mm_string[i];
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
            // status_flag = mm_string[i];
            i++;
        } else { // if not present, set to '.'
            // status_flag = '.';
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

        if(skip_counts_len == 0) { // no skip counts, no modification
            continue;
        }

        int base_rank = -1;

        int ml_idx = ml_start_idx;
        for(int c=0; c<skip_counts_len; c++) {
            base_rank += skip_counts[c] + 1;
            char mb;
            int idx;
            int read_pos;

            if (modbase == 'N') {
                if(rev) {
                    read_pos = seq_len - base_rank - 1;
                } else {
                    read_pos = base_rank;
                }
                
            } else {
                if(rev) {
                    mb = base_complement_lookup[(int)modbase];
                } else {
                    mb = modbase;
                }

                idx = base_idx_lookup[(int)mb];

                // print_array(bases_pos_lens, 5, 'i');
                if(base_rank >= bases_pos_lens[idx]) {
                    WARNING("%d th base of %c not found in SEQ. %c base count is %d read_id:%s seq_len:%d mod.base:%c mod_codes:%s\n", base_rank, mb, mb, bases_pos_lens[idx], qname, seq_len, modbase, mod_codes);
                    continue;
                }
                ASSERT_MSG(base_rank < bases_pos_lens[idx], "%d th base of %c not found in SEQ. %c base count is %d read_id:%s seq_len:%d mod.base:%c mod_codes:%s\n", base_rank, mb, mb, bases_pos_lens[idx], qname, seq_len, modbase, mod_codes);
                
                if(rev) {
                    read_pos = seq_len - bases_pos[idx][bases_pos_lens[idx] - base_rank - 1] - 1;
                } else {
                    read_pos = bases_pos[idx][base_rank];
                }
            }

            ASSERT_MSG(read_pos>=0 && read_pos < seq_len, "Read pos cannot exceed seq len. read_pos: %d seq_len: %d\n", read_pos, seq_len);

            int ref_pos = aln_pairs[read_pos];
            if(core->opt.insertions) {
                ref_pos = ref_pos == -1 ? db->ins[bam_i][read_pos] : ref_pos;
            }

            if(ref_pos == -1) { // not aligned nor insertion
                if(mod_codes_len > 0) {
                    ml_idx = ml_start_idx + c*mod_codes_len + mod_codes_len - 1;
                }
                continue;
            }

            // mod prob per each mod code.
            for(int m=0; m<mod_codes_len; m++) {                
                ml_idx = ml_start_idx + c*mod_codes_len + m;

                // check required mod codes
                khint_t mk;

                mk = kh_get(modcodesm, core->opt.modcodes_map, "*"); // check for wildcard first
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

                if(core->opt.insertions) { // no need to check context for insertions

                } else if ((!rev && ref->is_context[req_mod->index][ref_pos] == 1) || (rev && ref->is_context[req_mod->index][ref_pos - 1] == 1)) { // in context
                } else {
                    continue;
                }

                ASSERT_MSG(ml_idx<ml_len, "Mod prob index mismatch. ml_idx:%d ml_len:%d\n", ml_idx, ml_len);
                uint8_t mod_prob = ml[ml_idx];
                ASSERT_MSG(mod_prob <= 255 && mod_prob>=0, "Invalid mod_prob:%d\n", mod_prob);

                int ins_offset = core->opt.insertions ? db->ins_offset[bam_i][read_pos] : 0;
                if(core->opt.subtool == FREQ) {
                    uint8_t is_mod = 0, is_called = 0;
                    uint8_t thresh = req_mod->thresh;
                    
                    if(mod_prob >= thresh){ // modified with mod_code
                        is_called = 1;
                        is_mod = 1;
                    } else if(mod_prob <= 255-thresh){ // not modified with mod_code
                        is_called = 1;
                    } else { // ambiguous
                        continue;
                    }
                    
                    update_freq_map(db->freq_maps[bam_i], tname, ref_pos, ins_offset, mod_code, strand, haplotype, is_called, is_mod);
                } else if (core->opt.subtool == VIEW) {
                    add_view_entry(db->view_maps[bam_i], tname, ref_pos, ins_offset, mod_code, strand, haplotype, mod_prob, read_pos);
                }
            }

        }
        ml_start_idx = ml_idx + 1;

    }
}
