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

static const int valid_bases[256] = {
    ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1,
    ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1, ['u'] = 1, ['n'] = 1
};

static const int valid_strands[256] = {
    ['+'] = 1,
    ['-'] = 1
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

//C:mhfcC, T:gebT, U:U, A:aA, G:oG, N:nN
static const int mod_code_to_idx[256] = {
    ['m'] = 0, ['h'] = 1, ['f'] = 2, ['c'] = 3, ['C'] = 4,
    ['g'] = 5, ['e'] = 6, ['b'] = 7, ['T'] = 8,
    ['U'] = 9,
    ['a'] = 10, ['A'] = 11,
    ['o'] = 12, ['G'] = 13,
    ['n'] = 14, ['N'] = 15
};

static const char idx_to_mod_code[16] = {
    [0] = 'm', [1] = 'h', [2] = 'f', [3] = 'c', [4] = 'C',
    [5] = 'g', [6] = 'e', [7] = 'b', [8] = 'T',
    [9] = 'U',
    [10] = 'a', [11] = 'A',
    [12] = 'o', [13] = 'G',
    [14] = 'n', [15] = 'N'
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

void decode_key(char *key, char **chrom, int *pos, char *mod_code, char *strand){
    char* token = strtok(key, "\t");
    *chrom = calloc(strlen(token)+1, sizeof(char));
    MALLOC_CHK(*chrom);
    strcpy(*chrom, token);
    *pos = atoi(strtok(NULL, "\t"));
    *mod_code = strtok(NULL, "\t")[0];
    *strand = strtok(NULL, "\t")[0];
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

        modbase_t ** bases = db->modbases[i];

        for(int j=0;j<N_MODS;j++){
            for(int seq_i=0;seq_i<seq_len;seq_i++){
                modbase_t * base = &bases[j][seq_i];

                uint8_t mod_prob = base->mods_prob;
                if(mod_prob == 0){ // no modification
                    continue;
                }

                if(base->is_aln_cpg != 2){ // not aln_cpg
                    continue;
                }

                char mod_code = idx_to_mod_code[j];
                if(is_required_mod_code(mod_code, core->opt.mod_codes) == 0){
                    continue;
                }

                char *key = make_key(tname, base->ref_pos, mod_code, strand);            
                khiter_t k = kh_get(freqm, freq_map, key);
                if (k == kh_end(freq_map)) { // not found, add to map
                    freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
                    MALLOC_CHK(freq);
                    freq->mod_code = mod_code;
                    freq->n_called = 1;
                    freq->n_mod = (mod_prob/255.0) >= mod_thresh(mod_code, core->opt.mod_codes, core->opt.mod_threshes) ? 1 : 0;
                    int ret;
                    k = kh_put(freqm, freq_map, key, &ret);
                    kh_value(freq_map, k) = freq;
                } else { // found, update the values
                    free(key);
                    freq_t * freq = kh_value(freq_map, k);
                    freq->n_called += 1;
                    freq->n_mod += (mod_prob/255.0) >= mod_thresh(mod_code, core->opt.mod_codes, core->opt.mod_threshes) ? 1 : 0;
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
                double freq_value = (double)freq->n_mod*100/freq->n_called;
                char *contig = NULL;
                int ref_pos;
                char mod_code;
                char strand;
                char * key = (char *) kh_key(freq_map, k);
                decode_key(key, &contig, &ref_pos, &mod_code, &strand);
                int end = ref_pos+1;
                fprintf(core->opt.output_fp, "%s\t%d\t%d\t%c\t%d\t%c\t%d\t%d\t255,0,0\t%d\t%f\n", contig, ref_pos, end, mod_code, freq->n_called, strand, ref_pos, end, freq->n_called, freq_value);
                free(contig);
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
                char * contig = NULL;
                int ref_pos;
                char mod_code;
                char strand;
                char * key = (char *) kh_key(freq_map, k);
                decode_key(key, &contig, &ref_pos, &mod_code, &strand);
                fprintf(core->opt.output_fp, "%s\t%d\t%d\t%c\t%d\t%d\t%f\t%c\n", contig, ref_pos, ref_pos, strand, freq->n_called, freq->n_mod, freq_value, freq->mod_code);
                free(contig);
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

static void get_bases(db_t *db, int32_t bam_i, const char *mm_string, uint8_t *ml, uint32_t ml_len, int *aln_pairs, bam_hdr_t *hdr, bam1_t *record) {
    const char *qname = bam_get_qname(record);
    int8_t rev = bam_is_rev(record);
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = (tid >= 0) ? hdr->target_name[tid] : "*";
    uint8_t *seq = bam_get_seq(record);
    uint32_t seq_len = record->core.l_qseq;
    char strand = rev ? '-' : '+';

    // 5 int arrays to keep base pos of A, C, G, T, N bases.
    // A: 0, C: 1, G: 2, T: 3, U:4, N: 5
    // so that, nth base of A is at base_pos[0][n] and so on.
    int **bases_pos = db->bases_pos[bam_i];
    int bases_pos_lens[N_BASES] = {0};
    memset(db->mod_codes[bam_i], 0, N_MODS);

    int i;
    for(i=0;i<seq_len;i++){
        modbase_t base;;
        int base_char = seq_nt16_str[bam_seqi(seq, i)];
        int idx = base_idx_lookup[(int)base_char];
        bases_pos[idx][bases_pos_lens[idx]++] = i;

        base.ref_pos = aln_pairs[i];
        base.is_aln_cpg = aln_pairs[i] == -1 ? 0 : 1;
        // check if the base belongs to a cpg site using the ref
        if(base.is_aln_cpg){ // if aligned
            int ref_pos = base.ref_pos;
            ref_t *ref = get_ref(tname);

            ASSERT_MSG(ref != NULL, "Contig %s not found in ref_map\n", tname);
            ASSERT_MSG(ref_pos >= 0 && ref_pos < ref->ref_seq_length, "ref_pos:%d ref_len:%d\n", ref_pos, ref->ref_seq_length);
            ASSERT_MSG(ref->ref_seq_length == hdr->target_len[tid], "ref_len:%d target_len:%d\n", ref->ref_seq_length, hdr->target_len[tid]);

            // check if the base is a CpG site
            char * ref_seq = ref->forward;
            if ((ref_pos + 1 < ref->ref_seq_length && strand == '+' && ref_seq[ref_pos] == 'C' && ref_seq[ref_pos + 1] == 'G') ||
                (ref_pos > 0 && strand == '-' && ref_seq[ref_pos] == 'G' && ref_seq[ref_pos - 1] == 'C')) {
                base.is_aln_cpg = 2;
            }
        }

        base.mods_prob = 0;
        for(int j=0;j<N_MODS;j++){
            db->modbases[bam_i][j][i] = base;
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
            mod_codes[j] = mm_string[i];
            j++;

            i++;
        }
        // mod_codes[j] = '\0';
        mod_codes_len = j;

        ASSERT_MSG(mod_codes_len <= 16, "mod_codes_len:%d\n", mod_codes_len);

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
            ASSERT_MSG(l > 0, "invalid skip count:%d.\n", l);
            sscanf(skip_count_str, "%d", &skip_counts[k]);
            ASSERT_MSG(skip_counts[k] >= 0, "skip count cannot be negative: %d.\n", skip_counts[k]);
            
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
                read_pos = bases_pos[idx][bases_pos_lens[idx] - base_rank - 1];
                read_pos = seq_len - read_pos - 1;
            } else {
                read_pos = bases_pos[idx][base_rank];
            }

            ASSERT_MSG(read_pos>=0 && read_pos < seq_len, "Base pos cannot exceed seq len. read_pos: %d seq_len: %d\n", read_pos, seq_len);

            // mod prob per each mod code. TO-DO: need to change when code is ChEBI id
            for(int m=0; m<mod_codes_len; m++) {
                // get the mod prob
                ml_idx = ml_start_idx + c*mod_codes_len + m;
                ASSERT_MSG(ml_idx<ml_len, "ml_idx:%d ml_len:%d\n", ml_idx, ml_len);
                uint8_t mod_prob = ml[ml_idx];
                ASSERT_MSG(mod_prob <= 255 && mod_prob>=0, "mod_prob:%d\n", mod_prob);

                db->modbases[bam_i][mod_code_to_idx[(int)mod_codes[m]]][read_pos].mods_prob = mod_prob;

            }

        }
        ml_start_idx = ml_idx + 1;

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

        int8_t rev = bam_is_rev(record);
        char strand = rev ? '-' : '+';

        modbase_t ** bases = db->modbases[i];
        for(int j=0;j<N_MODS;j++){
            for(int seq_i=0;seq_i<seq_len;seq_i++){
                modbase_t base = bases[j][seq_i];
            
                uint8_t mod_prob = base.mods_prob;
                if(mod_prob == 0){ // no modification
                    continue;
                }
                char mod_code = idx_to_mod_code[j];
                double prob = mod_prob/255.0;

                if(base.is_aln_cpg != 2 ||  is_required_mod_code_and_thresh(mod_code, core->opt.mod_codes, prob, core->opt.mod_threshes) == 0){
                    continue;
                }
                fprintf(core->opt.output_fp, "%s\t%d\t%c\t%s\t%d\t%c\t%f\n", tname, base.ref_pos, strand, qname, seq_i, mod_code, prob);
            
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

    int * aln_pairs = db->aln[i];
    get_bases(db, i, mm, ml, ml_len, aln_pairs, hdr, record);

}
