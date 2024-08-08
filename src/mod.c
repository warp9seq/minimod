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

#define INIT_MODS 100
#define INIT_SKIP_COUNTS 100
#define INIT_MOD_CODES 1
#define INIT_BASE_POS 100
#define INIT_MOD_BASES 1
#define N_BASES 6 // A, C, G, T, N, U

static const int valid_bases[256] = {
    ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1,
    ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1, ['u'] = 1, ['n'] = 1
};

static const int valid_strands[256] = {
    ['+'] = 1,
    ['-'] = 1
};

static const int valid_mod_codes[256] = {
    // ['0'] = 1, ['1'] = 1, ['2'] = 1, ['3'] = 1, ['4'] = 1, ['5'] = 1, ['6'] = 1, ['7'] = 1, ['8'] = 1, ['9'] = 1, // for ChEBI ids
    ['a'] = 1, ['b'] = 1, ['c'] = 1, ['e'] = 1, ['f'] = 1, ['g'] = 1, ['h'] = 1, ['m'] = 1, ['n'] = 1, ['o'] = 1, 
    ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1
};

static const int base_idx_lookup[256] = {
    ['A'] = 0,
    ['C'] = 1,
    ['G'] = 2,
    ['T'] = 3,
    ['U'] = 4,
    ['N'] = 5,
    ['a'] = 0,
    ['c'] = 1,
    ['g'] = 2,
    ['t'] = 3,
    ['u'] = 4,
    ['n'] = 5,
};

static const char base_complement_lookup[256] = {
    ['A'] = 'T',
    ['C'] = 'G',
    ['G'] = 'C',
    ['T'] = 'A',
    ['U'] = 'A',
    ['N'] = 'N',
    ['a'] = 't',
    ['c'] = 'g',
    ['g'] = 'c',
    ['t'] = 'a',
    ['u'] = 'a',
    ['n'] = 'n',
};

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
        WARNING("%s tag not found in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    const char *mm_str = bam_aux2Z(data);
    if(mm_str == NULL){
        WARNING("%s tag could not be decoded for %s. Is it type Z?",tag, bam_get_qname(record));
        return NULL;
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
    }

    // check if the type of the tag is array of integers
    const char aux_type = data[0];
    if (aux_type != 'B') {
        ERROR("%s tag is not of type B in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    // get the array len
    uint32_t len = bam_auxB_len(data);
    if(len == 0){
        WARNING("%s tag array length is 0 in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    // check if the array type is uint8_t
    const char aux_array_type = data[1];
    if (aux_array_type != 'C') {
        ERROR("%s array tag type '%c' is not of type C in read %s",tag, aux_array_type, bam_get_qname(record));
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

static int is_required_mod_code(char mod_code, char * mod_codes){

    if(mod_codes[0] == '*') {
        return 1;
    }

    int i=0;
    while(mod_codes[i] != '\0'){
        if(mod_codes[i] == mod_code){
            return 1;
        }
        i++;
    }
    return 0;
}

static int is_required_mod_code_and_thresh(char mod_code, char * mod_codes, double mod_prob, double * mod_threshes){
    
        if(mod_codes[0] == '*' && mod_prob >= mod_threshes[0]) {
            return 1;
        }
    
        int i=0;
        while(mod_codes[i] != '\0'){
            if(mod_codes[i] == mod_code && mod_prob >= mod_threshes[i]){
                return 1;
            }
            i++;
        }
        return 0;
}

void destroy_freq_map(khash_t(freqm)* freq_map){
    khint_t k;
    for (k = kh_begin(freq_map); k != kh_end(freq_map); k++) {
        if (kh_exist(freq_map, k)) {
            char *key = (char *) kh_key(freq_map, k);
            freq_t* freq = kh_value(freq_map, k);
            free((char*)freq->contig);
            free(freq);
            free(key);
        }
    }
    kh_destroy(freqm, freq_map);
}

char* make_key(const char *chrom, int pos, char mod_code, char strand){
    int start_strlen = snprintf(NULL, 0, "%d", pos);
    int key_strlen = strlen(chrom) + start_strlen  + 7;
    
    char* key = (char *)malloc(key_strlen * sizeof(char));
    MALLOC_CHK(key);
    snprintf(key, key_strlen, "%s\t%d\t%c\t%c", chrom, pos, mod_code, strand);
    return key;
}

static double mod_thresh(char mod_code, char * mod_codes, double * mod_threshes){
    if(mod_codes[0] == '*') {
        return mod_threshes[0];
    }

    int i=0;
    while(mod_codes[i] != '\0'){
        if(mod_codes[i] == mod_code){
            return mod_threshes[i];
        }
        i++;
    }
    return -1;
}
// freq_map is used by multiple threads, so need to lock it
void update_freq_map(core_t * core, db_t * db) {

    for(int i=0;i<db->n_bam_recs;i++){
        bam1_t *record = db->bam_recs[i];
        int8_t rev = bam_is_rev(record);

        bam_hdr_t * hdr = core->bam_hdr;

        int32_t tid = record->core.tid;
        assert(tid < hdr->n_targets);
        const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";

        uint32_t seq_len = record->core.l_qseq;

        char strand = rev ? '-' : '+';

        modbase_t * modbases = db->modbases[i];

        for(int seq_i=0;seq_i<seq_len;seq_i++){
            modbase_t base = modbases[seq_i];
            for(int j=0;j<base.mods_len;j++){
                mod_t mod = base.mods[j];

                if(base.is_aln == 0 || base.is_cpg == 0 || is_required_mod_code(mod.mod_code, core->opt.mod_codes) == 0) {
                    continue;
                }

                char *key = make_key(tname, base.ref_pos, mod.mod_code, strand);
                khash_t(freqm) *freq_map = core->freq_map;
                khiter_t k = kh_get(freqm, freq_map, key);
                if (k == kh_end(freq_map)) { // not found, add to map
                    freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
                    MALLOC_CHK(freq);
                    char* contig = (char *)calloc(strlen(tname)+1, sizeof(char));
                    MALLOC_CHK(contig);
                    strcpy(contig, tname);
                    freq->contig = contig;
                    freq->ref_pos = base.ref_pos;
                    freq->mod_code = mod.mod_code;

                    freq->n_called = 1;
                    freq->n_mod = (mod.mod_prob/255.0) >= mod_thresh(mod.mod_code, core->opt.mod_codes, core->opt.mod_threshes) ? 1 : 0;
                    freq->strand = mod.mod_strand;
                    int ret;
                    k = kh_put(freqm, freq_map, key, &ret);
                    kh_value(freq_map, k) = freq;
                } else { // found, update the values
                    free(key);
                    freq_t * freq = kh_value(freq_map, k);
                    freq->n_called += 1;
                    freq->n_mod += (mod.mod_prob/255.0) >= mod_thresh(mod.mod_code, core->opt.mod_codes, core->opt.mod_threshes) ? 1 : 0;
                }
            }
        }
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
                double freq_value = (double)freq->n_mod/freq->n_called*100;
                int end = freq->ref_pos+1;
                fprintf(core->opt.output_fp, "%s\t%d\t%d\t%c\t%d\t%c\t%d\t%d\t255,0,0\t%d\t%f\n", freq->contig, freq->ref_pos, end, freq->mod_code, freq->n_called, freq->strand, freq->ref_pos, end, freq->n_called, freq_value);
            }
        }
        
    } else {
        // contig, start, end, strand, n_called, n_mod, freq, mod_code
        fprintf(core->opt.output_fp, "contig\tstart\tend\tstrand\tn_called\tn_mod\tfreq\tmod_code\n");
        khint_t k;
        for (k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
            if (kh_exist(freq_map, k)) {
                freq_t* freq = kh_value(freq_map, k);
                double freq_value = (double)freq->n_mod/freq->n_called;
                fprintf(core->opt.output_fp, "%s\t%d\t%d\t%c\t%d\t%d\t%f\t%c\n", freq->contig, freq->ref_pos, freq->ref_pos, freq->strand, freq->n_called, freq->n_mod, freq_value, freq->mod_code);
            }
        }
    }

    if(core->opt.output_fp != stdout){
        fclose(core->opt.output_fp);
    }
}

int * get_aln(bam_hdr_t *hdr, bam1_t *record){
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    int32_t pos = record->core.pos;
    int32_t end = bam_endpos(record);

    int8_t rev = bam_is_rev(record);

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

static void get_bases(const char *mm_string, modbase_t * bases, uint8_t * ml, uint32_t ml_len, int * aln_pairs, bam_hdr_t *hdr, bam1_t *record){

    const char *qname = bam_get_qname(record);
    int8_t rev = bam_is_rev(record);

    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";

    uint8_t *seq = bam_get_seq(record);
    uint32_t seq_len = record->core.l_qseq;

    char strand = rev ? '-' : '+';

    // 5 int arrays to keep base pos of A, C, G, T, N bases.
    // A: 0, C: 1, G: 2, T: 3, U:4, N: 5
    // so that, nth base of A is at base_pos[0][n] and so on.
    int * bases_pos[N_BASES];
    int bases_pos_lens[N_BASES];
    int bases_pos_caps[N_BASES];
    int i;
    for(i=0;i<N_BASES;i++){
        bases_pos_caps[i] = INIT_BASE_POS;
        bases_pos[i] = (int *)malloc(sizeof(int)*bases_pos_caps[i]);
        MALLOC_CHK(bases_pos[i]);
        bases_pos_lens[i] = 0;
    }

    // keep record of base pos of A, C, G, T, N bases
    for(i=0; i<seq_len; i++){
        int base_char = seq_nt16_str[bam_seqi(seq, i)];
        int idx = base_idx_lookup[(int)base_char];
        if(bases_pos_lens[idx] >= bases_pos_caps[idx]){
            bases_pos_caps[idx] *= 2;
            bases_pos[idx] = (int *)realloc(bases_pos[idx], sizeof(int)*bases_pos_caps[idx]);
            MALLOC_CHK(bases_pos[idx]);
        }
        bases_pos[idx][bases_pos_lens[idx]] = i;
        bases_pos_lens[idx]++;
        
    }

    for(i=0;i<seq_len;i++){
        bases[i].mods_cap = INIT_MOD_BASES;
        bases[i].mods_len = 0;
        bases[i].ref_pos = aln_pairs[i];
        bases[i].base = seq_nt16_str[bam_seqi(seq, i)];
        bases[i].is_aln = aln_pairs[i] == -1 ? 0 : 1;
        bases[i].is_cpg = 0;
        bases[i].base = seq_nt16_str[bam_seqi(seq, i)];
        // check if the base belongs to a cpg site using the ref
        if(bases[i].is_aln){ // if aligned
            int ref_pos = bases[i].ref_pos;

            ref_t *ref = get_ref(tname);

            ASSERT_MSG(ref != NULL, "Contig %s not found in ref_map\n", tname);
            ASSERT_MSG(ref_pos >= 0 && ref_pos < ref->ref_seq_length, "ref_pos:%d ref_len:%d\n", ref_pos, ref->ref_seq_length);
            ASSERT_MSG(ref->ref_seq_length == hdr->target_len[tid], "ref_len:%d target_len:%d\n", ref->ref_seq_length, hdr->target_len[tid]);
                        
            // check if the base is a CpG site
            char * ref_seq = ref->forward;
            if(ref_pos+1 < ref->ref_seq_length && strand=='+' && ref_seq[ref_pos] == 'C' && ref_seq[ref_pos+1] == 'G'){
                bases[i].is_cpg = 1;
            } else if(ref_pos > 0 && strand=='-' && ref_seq[ref_pos] == 'G' && ref_seq[ref_pos-1] == 'C'){
                bases[i].is_cpg = 1;
            }

        }
    }

    int mm_str_len = strlen(mm_string);
    i = 0;
    int ml_start_idx = 0;

    char modbase;
    char mod_strand;
    char * mod_codes;
    int mod_codes_cap;
    int mod_codes_len;
    int * skip_counts;
    int skip_counts_cap;
    int skip_counts_len;
    // char status_flag;

    // allocate initial memory for skip counts
    skip_counts_cap = INIT_SKIP_COUNTS;
    skip_counts = (int *) malloc(skip_counts_cap * sizeof(int));
    MALLOC_CHK(skip_counts);

    // allocate initial memory for modification codes
    mod_codes_cap = INIT_MOD_CODES;
    mod_codes = (char *) malloc(mod_codes_cap * sizeof(char));
    MALLOC_CHK(mod_codes);
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
            mod_strand = mm_string[i];
            i++;
        }

        // get base modification codes. can handle multiple codes giver as chars. TO-DO: handle when given as a ChEBI id
        int j = 0;
        while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';' && mm_string[i] != '?' && mm_string[i] != '.') {

            if (j >= mod_codes_cap) {
                mod_codes_cap *= 2;
                mod_codes = (char *) realloc(mod_codes, mod_codes_cap * sizeof(char));
                MALLOC_CHK(mod_codes);
            }

            ASSERT_MSG(valid_mod_codes[(int)mm_string[i]], "Invalid base modification code:%c\n", mm_string[i]);
            mod_codes[j] = mm_string[i];
            j++;

            i++;
        }
        // mod_codes[j] = '\0';
        mod_codes_len = j;

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

            if (k >= skip_counts_cap) {
                skip_counts_cap *= 2;
                skip_counts = (int *) realloc(skip_counts, skip_counts_cap * sizeof(int));
                MALLOC_CHK(skip_counts);
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
            ASSERT_MSG(l > 0, "invalid skip count:%d.\n", l);
            sscanf(skip_count_str, "%d", &skip_counts[k]);
            ASSERT_MSG(skip_counts[k] >= 0, "skip count cannot be negative: %d.\n", skip_counts[k]);
            
            k++;
        }
        skip_counts_len = k;
        i++;

        if(skip_counts_len == 0) { // no skip counts, no modification
            free(skip_counts);
            free(mod_codes);
            continue;
        }
 
        int base_rank = -1;

        int ml_idx = ml_start_idx;
        for(int c=0; c<skip_counts_len; c++) {
            base_rank += skip_counts[c] + 1;
            char mb;
            int idx;
            int read_pos;

            if(rev){
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
                read_pos = bases_pos[idx][bases_pos_lens[idx] - base_rank - 1];
                read_pos = seq_len - read_pos - 1;
            } else {
                read_pos = bases_pos[idx][base_rank];
            }

            ASSERT_MSG(read_pos < seq_len, "Base pos cannot exceed seq len. read_pos: %d seq_len: %d\n", read_pos, seq_len);

            modbase_t* base = &bases[read_pos];

            if(base->mods_len == 0 && mod_codes_len > 0) {
                base->mods = (mod_t *)malloc(sizeof(mod_t)*base->mods_cap);
                MALLOC_CHK(base->mods);
            }
            
            int mod_i = base->mods_len;
            // mod prob per each mod code. TO-DO: need to change when code is ChEBI id
            for(int m=0; m<mod_codes_len; m++) {
                // get the mod prob
                ml_idx = ml_start_idx + c*mod_codes_len + m;
                ASSERT_MSG(ml_idx<ml_len, "ml_idx:%d ml_len:%d\n", ml_idx, ml_len);
                uint8_t mod_prob = ml[ml_idx];
                ASSERT_MSG(mod_prob <= 255 && mod_prob>=0, "mod_prob:%d\n", mod_prob);

                // add to mods
                if(mod_i >= base->mods_cap){
                    base->mods_cap *= 2;
                    base->mods = (mod_t *)realloc(base->mods, sizeof(mod_t)*base->mods_cap);
                    MALLOC_CHK(base->mods);
                }

                base->mods[mod_i].mod_code = mod_codes[m];
                base->mods[mod_i].mod_strand = mod_strand;
                base->mods[mod_i].mod_prob = mod_prob;
                base->mods[mod_i].mod_base = modbase;
                mod_i++;
            }
            base->mods_len = mod_i;

        }
        ml_start_idx = ml_idx + 1;

    }

    free(skip_counts);
    free(mod_codes);

    // free base_pos
    for(int b=0;b<N_BASES;b++){
        free(bases_pos[b]);
    }

}

void print_view_output(core_t* core, db_t* db) {
    fprintf(core->opt.output_fp, "ref_contig\tref_pos\tstrand\tread_id\tread_pos\tmod_code\tmod_prob\n");
    for(int i=0;i<db->n_bam_recs;i++){
        bam1_t *record = db->bam_recs[i];
        

        bam_hdr_t * hdr = core->bam_hdr;
        const char *qname = bam_get_qname(record);
        int32_t tid = record->core.tid;
        assert(tid < hdr->n_targets);
        const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";
        int32_t seq_len = record->core.l_qseq;

        modbase_t * modbases = db->modbases[i];

        for(int seq_i=0;seq_i<seq_len;seq_i++){
            for(int j=0;j<modbases[seq_i].mods_len;j++){
                mod_t mod = modbases[seq_i].mods[j];
                modbase_t base = modbases[seq_i];
                
                double mod_prob = mod.mod_prob/255.0;

                if(base.is_aln == 0 ||  is_required_mod_code_and_thresh(mod.mod_code, core->opt.mod_codes, mod_prob, core->opt.mod_threshes) == 0){
                    continue;
                }
                fprintf(core->opt.output_fp, "%s\t%d\t%c\t%s\t%d\t%c\t%f\n", tname, base.ref_pos, mod.mod_strand, qname, seq_i, mod.mod_code, mod_prob);
            
                
            }
        }
    }

}

// static void free_mod_tags(mod_tag_t *mod_tags, uint32_t len){
//     for(int i=0;i<len;i++){
//         free(mod_tags[i].mod_codes);
//         free(mod_tags[i].skip_counts);
//     }
//     free(mod_tags);
// }

void modbases_single(core_t* core, db_t* db, int32_t i) {
    bam1_t *record = db->bam_recs[i];

    const char *mm = db->mm[i];
    uint32_t ml_len = db->ml_lens[i];
    uint8_t *ml = db->ml[i];

    bam_hdr_t *hdr = core->bam_hdr;

    int * aln_pairs = db->aln[i];
    get_bases(mm, db->modbases[i], ml, ml_len, aln_pairs, hdr, record);

}
