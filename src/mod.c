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

static const int valid_bases[256] = { ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1, ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1, ['u'] = 1, ['n'] = 1 };
static const int valid_strands[256] = { ['+'] = 1, ['-'] = 1 };
static const int base_idx_lookup[256] = { ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['U'] = 4, ['N'] = 5, ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3, ['u'] = 4, ['n'] = 5 };
static const char base_complement_lookup[256] = { ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A', ['N'] = 'N', ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a', ['n'] = 'n' };
static const uint8_t valid_mod_codes[256] = { ['a'] = 1, ['b'] = 1, ['c'] = 1, ['e'] = 1, ['f'] = 1, ['g'] = 1, ['h'] = 1, ['m'] = 1, ['n'] = 1, ['o'] = 1,  ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1 };
static int mod_code_idx[256]; // address by mod code to get the index
static const char* default_context[256] = { ['m'] = "CG", ['h'] = "CG", ['f'] = "C", ['c'] = "C", ['C'] = "C", ['g'] = "T", ['e'] = "T", ['b'] = "T", ['T'] = "T", ['U'] = "U", ['a'] = "A", ['A'] = "A", ['o'] = "G", ['G'] = "G", ['n'] = "N", ['N'] = "N" };

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

    for(i=0;i<256;i++){ // reset the mod_code_idx array to -1 (not requested)
        mod_code_idx[i] = -1;
    }

    i=0;
    while(1){
        if(opt->mod_codes_str[i] == '\0'){ // end of string
            break;
        }
        
        char c = opt->mod_codes_str[i];
        if(valid_mod_codes[(int)c] == 0){
            ERROR("Invalid modification code %c", c);
            exit(EXIT_FAILURE);
        }
        if(mod_code_idx[(int)c] != -1){ // already set
            ERROR("Duplicate modification code %c", c);
            exit(EXIT_FAILURE);
        }
        
        mod_code_idx[(int)c] = n_codes;
        opt->req_mod_codes[n_codes] = c;
        
        i++;

        int context_cap = 2;
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
                    ERROR("Invalid character %c in context for modification code %c", opt->mod_codes_str[i], c);
                    exit(EXIT_FAILURE);
                }
                if(j >= context_cap){
                    context_cap *= 2;
                    context = (char *)realloc(context, context_cap * sizeof(char) + 1);
                    MALLOC_CHK(context);
                }
                context[j] = opt->mod_codes_str[i];
                i++;
                j++;
            }
            if(opt->mod_codes_str[i] == '\0'){
                ERROR("Context not closed with a ] for modification code %c", c);
                exit(EXIT_FAILURE);
            }
            context[j] = '\0';
            if(is_star && j > 1){
                ERROR("Invalid context for modification code %c. * should be the only character within [ and ]", c);
                exit(EXIT_FAILURE);
            }
            opt->req_mod_contexts[n_codes] = context;
            i++;
            if(opt->mod_codes_str[i] == ','){ // more modification codes
                i++;
            }
        } else if(opt->mod_codes_str[i] == ','){ // context is *
            INFO("Context not provided for modification code %c. Using %s", c, default_context[(int)c]);
            strcpy(context, default_context[(int)c]);
            opt->req_mod_contexts[n_codes] = context;
            i++;
        } else if(opt->mod_codes_str[i] == '\0'){
            INFO("Context not provided for modification code %c. Using %s", c, default_context[(int)c]);
            strcpy(context, default_context[(int)c]);
            opt->req_mod_contexts[n_codes] = context;
        } else {
            ERROR("Invalid character %c after modification code %c", opt->mod_codes_str[i], c);
            exit(EXIT_FAILURE);
        }
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

        INFO("Modification code: %c, Context: %s, Threshold: %f", opt->req_mod_codes[n_thresh], opt->req_mod_contexts[n_thresh], d);

        opt->req_threshes[n_thresh] = d*255;
        
        free(thresh_str);
        n_thresh++;
        if(opt->mod_threshes_str[i] == '\0'){
            break;
        }
        i++;
    }

    if(n_thresh==1){
        for(int i=1;i<opt->n_mods;i++){
            opt->req_threshes[i] = opt->req_threshes[0];
            INFO("Modification code: %c, Context: %s, Threshold: %f", opt->req_mod_codes[i], opt->req_mod_contexts[i], d);
        }
    } else if(n_thresh != opt->n_mods){
        ERROR("Number of modification codes and thresholds do not match. Codes:%d, Thresholds:%d",opt->n_mods,n_thresh);
        exit(EXIT_FAILURE);
    }
}


void print_view_options(opt_t *opt) {
    int i;
    for(i=0; i < opt->n_mods; i++){
        INFO("Modification code: %c, Context: %s", opt->req_mod_codes[i], opt->req_mod_contexts[i]);
    }
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

char* make_key(const char *chrom, int pos, uint16_t ins_offset, char mod_code, char strand, int haplotype){
    int start_strlen = snprintf(NULL, 0, "%d", pos);
    int offset_strlen = snprintf(NULL, 0, "%d", ins_offset);
    int haplotype_strlen = snprintf(NULL, 0, "%d", haplotype);
    int key_strlen = strlen(chrom) + start_strlen  + offset_strlen + haplotype_strlen + 8;
    
    char* key = (char *)malloc(key_strlen * sizeof(char));
    MALLOC_CHK(key);
    snprintf(key, key_strlen, "%s\t%d\t%u\t%c\t%c\t%d", chrom, pos, ins_offset, mod_code, strand, haplotype);
    return key;
}

void decode_key(char *key, char **chrom, int *pos, uint16_t * ins_offset, char *mod_code, char *strand, int *haplotype){
    char* token = strtok(key, "\t");
    *chrom = calloc(strlen(token)+1, sizeof(char));
    MALLOC_CHK(*chrom);
    strcpy(*chrom, token);
    *pos = atoi(strtok(NULL, "\t"));
    *ins_offset = strtoul(strtok(NULL, "\t"), NULL, 10);
    *mod_code = strtok(NULL, "\t")[0];
    *strand = strtok(NULL, "\t")[0];
    *haplotype = atoi(strtok(NULL, "\t"));
}

void update_freq_map(core_t * core, db_t * db) {
    bam_hdr_t * hdr = core->bam_hdr;
    khash_t(freqm) *freq_map = core->freq_map;

    for(int i=0;i<db->n_bam_recs;i++){
        bam1_t *record = db->bam_recs[i];
        int8_t rev = bam_is_rev(record);
        int32_t tid = record->core.tid;
        assert(tid < hdr->n_targets);
        const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";
        uint32_t seq_len = record->core.l_qseq;
        char strand = rev ? '-' : '+';

        for(int j=0;j<core->opt.n_mods;j++){
            char mod_code = core->opt.req_mod_codes[j];
            for(int seq_i=0;seq_i<seq_len;seq_i++){
                
                if(db->mod_prob[i][j][seq_i] == -1){ // not called
                    continue;
                }

                if(core->opt.insertions && db->aln[i][seq_i] == -1 && db->ins[i][seq_i] == -1){ // not aligned, not inserted
                    continue;
                } else if(!core->opt.insertions && db->aln[i][seq_i] == -1){ // not aligned
                    continue;
                }

                int ref_pos = db->aln[i][seq_i] == -1 ? db->ins[i][seq_i] : db->aln[i][seq_i];
                
                uint8_t is_mod = 0, is_called = 0;
                uint8_t thresh = core->opt.req_threshes[j];
                int mod_prob = db->mod_prob[i][j][seq_i];
                
                if(mod_prob >= thresh){ // modified with mod_code
                    is_called = 1;
                    is_mod = 1;
                } else if(mod_prob <= 255-thresh){ // not modified with mod_code
                    is_called = 1;
                } else { // ambiguous
                    continue;
                }

                int ins_offset = 0;
                int haplotype = 0;
                char *key = NULL;
                if(core->opt.insertions){
                    ins_offset = db->ins_offset[i][seq_i];
                }
                if(core->opt.haplotypes){
                    haplotype = db->haplotypes[i];
                }

                key = make_key(tname, ref_pos, ins_offset, mod_code, strand, haplotype);

                khiter_t k = kh_get(freqm, freq_map, key);
                if (k == kh_end(freq_map)) { // not found, add to map
                    freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
                    MALLOC_CHK(freq);
                    freq->n_called = is_called;
                    freq->n_mod = is_mod;
                    int ret;
                    k = kh_put(freqm, freq_map, key, &ret);
                    kh_value(freq_map, k) = freq;
                } else { // found, update the values
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

                if(core->opt.haplotypes) {
                    // total counts
                    key = make_key(tname, ref_pos, ins_offset, mod_code, strand, -1);

                    k = kh_get(freqm, freq_map, key);
                    if (k == kh_end(freq_map)) { // not found, add to map
                        freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
                        MALLOC_CHK(freq);
                        freq->n_called = is_called;
                        freq->n_mod = is_mod;
                        int ret;
                        k = kh_put(freqm, freq_map, key, &ret);
                        kh_value(freq_map, k) = freq;
                    } else { // found, update the values
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
        }
    }
}

void print_freq_tsv_header(core_t * core) {
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

void print_freq_output(core_t * core) {
    khash_t(freqm) *freq_map = core->freq_map;

    if(core->opt.bedmethyl_out) {
        // chrom, start, end, mod_code, n_called, strand, start, end, "255,0,0",  n_called, freq
        khint_t k;
        for (k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
            if (kh_exist(freq_map, k)) {
                freq_t* freq = kh_value(freq_map, k);
                double freq_value = (double)freq->n_mod*100/freq->n_called;
                char *contig = NULL;
                int ref_pos;
                uint16_t ins_offset;
                char mod_code;
                char strand;
                int haplotype;
                char * key = (char *) kh_key(freq_map, k);
                decode_key(key, &contig, &ref_pos, &ins_offset, &mod_code, &strand, &haplotype);
                int end = ref_pos+1;
                fprintf(core->opt.output_fp, "%s\t%d\t%d\t%c\t%d\t%c\t%d\t%d\t255,0,0\t%d\t%f\n", contig, ref_pos, end, mod_code, freq->n_called, strand, ref_pos, end, freq->n_called, freq_value);
                free(contig);
            }
        }
        
    } else {
        // contig, start, end, strand, n_called, n_mod, freq, mod_code
        khint_t k;
        for (k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
            if (kh_exist(freq_map, k)) {
                freq_t* freq = kh_value(freq_map, k);
                double freq_value = (double)freq->n_mod/freq->n_called;
                char * contig = NULL;
                int ref_pos;
                uint16_t ins_offset;
                char mod_code;
                char strand;
                int haplotype;
                char * key = (char *) kh_key(freq_map, k);
                decode_key(key, &contig, &ref_pos, &ins_offset, &mod_code, &strand, &haplotype);

                fprintf(core->opt.output_fp, "%s\t%d\t%d\t%c\t%d\t%d\t%f\t%c", contig, ref_pos, ref_pos, strand, freq->n_called, freq->n_mod, freq_value, mod_code);

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
            }
        }
    }

    if(core->opt.output_fp != stdout){
        fclose(core->opt.output_fp);
    }
}

static void get_aln(core_t * core, int ** aln, int ** ins, int ** ins_offset, bam_hdr_t *hdr, bam1_t *record){
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

    int * aligned_pairs = *aln;
    int * ins_pos = NULL;
    int * ins_off = NULL;
    //fill the aligned_pairs array with -1
    for(int i=0;i<seq_len;i++){
        aligned_pairs[i] = -1;
    }

    if(core->opt.insertions){
        ins_pos = *ins;
        ins_off = *ins_offset;
        for(int i=0;i<seq_len;i++){
            ins_pos[i] = -1;
            ins_off[i] = 0;
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
                ins_pos[read_pos] = start;
                ins_off[read_pos] = offset;
            }

            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }
}

static void get_bases(core_t * core, db_t *db, int32_t bam_i, const char *mm_string, uint8_t *ml, uint32_t ml_len, int *aln_pairs, int * ins, int * ins_offsets, bam_hdr_t *hdr, bam1_t *record) {
    const char *qname = bam_get_qname(record);
    int8_t rev = bam_is_rev(record);
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = (tid >= 0) ? hdr->target_name[tid] : "*";
    uint8_t *seq = bam_get_seq(record);
    uint32_t seq_len = record->core.l_qseq;

    ref_t *ref = get_ref(tname);

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
        
        for(int j=0;j<core->opt.n_mods;j++){
            db->mod_prob[bam_i][j][i] = -1;
        }
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
        while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';' && mm_string[i] != '?' && mm_string[i] != '.') {

            ASSERT_MSG(valid_mod_codes[(int)mm_string[i]], "Invalid base modification code:%c\n", mm_string[i]);

            if(j >= db->mod_codes_cap[bam_i]) {
                db->mod_codes_cap[bam_i] *= 2;
                db->mod_codes[bam_i] = (char *)realloc(db->mod_codes[bam_i], sizeof(char) * (db->mod_codes_cap[bam_i] + 1)); // +1 for null terminator
                MALLOC_CHK(db->mod_codes[bam_i]);
            }

            mod_codes[j] = mm_string[i];
            j++;

            i++;
        }
        mod_codes[j] = '\0';
        mod_codes_len = j;

        ASSERT_MSG(mod_codes_len>0 && mod_codes_len <= 16, "mod_codes_len:%d\n", mod_codes_len);

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
                ref_pos = ref_pos == -1 ? ins[read_pos] : ref_pos;
            }

            if(ref_pos == -1) { // not aligned nor insertion
                if(mod_codes_len > 0) {
                    ml_idx = ml_start_idx + c*mod_codes_len + mod_codes_len - 1;
                }
                continue;
            }

            // mod prob per each mod code. TO-DO: need to change when code is ChEBI id
            for(int m=0; m<mod_codes_len; m++) {
                char mod_code = mod_codes[m];
                
                ml_idx = ml_start_idx + c*mod_codes_len + m;

                if(mod_code_idx[(int)mod_code]==-1){ // mod code not required
                    continue;
                }

                if(core->opt.insertions) { // no need to check context for insertions

                } else if ((!rev && ref->is_context[mod_code_idx[(int)mod_code]][ref_pos] == 1) || (rev && ref->is_context[mod_code_idx[(int)mod_code]][ref_pos - 1] == 1)) { // in context
                } else {
                    continue;
                }

                ASSERT_MSG(ml_idx<ml_len, "Mod prob index mismatch. ml_idx:%d ml_len:%d\n", ml_idx, ml_len);
                uint8_t mod_prob = ml[ml_idx];
                ASSERT_MSG(mod_prob <= 255 && mod_prob>=0, "Invalid mod_prob:%d\n", mod_prob);

                if(db->mod_prob[bam_i][mod_code_idx[(int)mod_code]][read_pos] != -1) {
                    WARNING("Multiple modification calls for the same base in read:%s ref_pos:%d mod_code:%c read_pos:%d\nKeeping the call with higher modification probability. Ignoring the other calls.\n", qname, ref_pos, mod_code, read_pos);
                    if(mod_prob > db->mod_prob[bam_i][mod_code_idx[(int)mod_code]][read_pos]) {
                        db->mod_prob[bam_i][mod_code_idx[(int)mod_code]][read_pos] = mod_prob;
                    }
                } else {
                    db->mod_prob[bam_i][mod_code_idx[(int)mod_code]][read_pos] = mod_prob;
                }      
            }

        }
        ml_start_idx = ml_idx + 1;

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
    for(int i=0;i<db->n_bam_recs;i++){
        bam1_t *record = db->bam_recs[i];
        bam_hdr_t * hdr = core->bam_hdr;
        const char *qname = bam_get_qname(record);
        int32_t tid = record->core.tid;
        assert(tid < hdr->n_targets);
        const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";
        int32_t seq_len = record->core.l_qseq;

        int8_t rev = bam_is_rev(record);
        char strand = rev ? '-' : '+';

        for(int j=0;j<core->opt.n_mods;j++){
            char mod_code = core->opt.req_mod_codes[j];
            for(int seq_i=0;seq_i<seq_len;seq_i++){
                if(db->mod_prob[i][j][seq_i] == -1){ // no modification
                    continue;
                }

                if(core->opt.insertions && db->aln[i][seq_i] == -1 && db->ins[i][seq_i] == -1){ // not aligned, not inserted
                    continue;
                } else if(!core->opt.insertions && db->aln[i][seq_i] == -1){ // not aligned
                    continue;
                }

                int ref_pos = db->aln[i][seq_i] == -1 ? db->ins[i][seq_i] : db->aln[i][seq_i];

                fprintf(core->opt.output_fp, "%s\t%d\t%c\t%s\t%d\t%c\t%f", tname, ref_pos, strand, qname, seq_i, mod_code, db->mod_prob[i][j][seq_i]/255.0);
                if(core->opt.insertions==1){
                    fprintf(core->opt.output_fp, "\t%d", db->ins_offset[i][seq_i]);
                }
                if(core->opt.haplotypes==1){
                    fprintf(core->opt.output_fp, "\t%d", db->haplotypes[i]);
                }
                fprintf(core->opt.output_fp, "\n");
            }
        }
    }

}


void modbases_single(core_t* core, db_t* db, int32_t i) {
    bam1_t *record = db->bam_recs[i];

    const char *mm = db->mm[i];
    uint32_t ml_len = db->ml_lens[i];
    uint8_t *ml = db->ml[i];

    bam_hdr_t *hdr = core->bam_hdr;

    if(core->opt.haplotypes) {
        db->haplotypes[i] = get_hp_tag(record);
    }

    if(core->opt.insertions) {
        get_aln(core, &(db->aln[i]), &(db->ins[i]), &(db->ins_offset[i]), hdr, record);
        get_bases(core, db, i, mm, ml, ml_len, db->aln[i], db->ins[i], db->ins_offset[i], hdr, record);
    } else {
        get_aln(core, &(db->aln[i]), NULL, NULL, hdr, record);
        get_bases(core, db, i, mm, ml, ml_len, db->aln[i], NULL, NULL, hdr, record);
    }

}
