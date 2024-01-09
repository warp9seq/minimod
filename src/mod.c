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
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>

// #include <sys/wait.h>
// #include <unistd.h>

#define INIT_MODS 100
#define INIT_SKIP_COUNTS 100
#define INIT_MODS_CODES 3

typedef struct {
    char base;
    char strand;
    char * mod_codes;
    int mod_codes_cap;
    int num_mod_codes;
    int * skip_counts;
    int skip_counts_cap;
    int num_skips;
    char status_flag;
    uint8_t * probs;
} mod_t;

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

// Function to extract methylated C base modification regions from MM string
static mod_t *extract_mods(const char *mm_string, const uint8_t *ml, uint32_t *len_mods) {
    if (mm_string == NULL || strlen(mm_string) == 0) {
        ERROR("%s","Error: Empty MM string.\n");
        return NULL;
    }

    // allocate initial memory for modifications
    int mods_cap = INIT_MODS;
    mod_t * mods = (mod_t *) malloc(mods_cap * sizeof(mod_t));
    MALLOC_CHK(mods);
    int num_mods = 0;

    int mm_str_len = strlen(mm_string);
    int i = 0;
    int probs_i = 0;
    while (i < mm_str_len) {

        if (num_mods >= INIT_MODS) {
            mods_cap *= 2;
            mods = (mod_t *) realloc(mods, mods_cap * sizeof(mod_t));
            MALLOC_CHK(mods);
            // die("Error: Too many modifications than INIT_MODS=%d.\n", INIT_MODS);
        }

        mod_t current_mod;
        memset(&current_mod, 0, sizeof(mod_t));

        // allocate initial memory for skip counts
        current_mod.skip_counts_cap = INIT_SKIP_COUNTS;
        current_mod.skip_counts = (int *) malloc(current_mod.skip_counts_cap * sizeof(int));
        MALLOC_CHK(current_mod.skip_counts);

        // allocate initial memory for probabilities
        current_mod.probs = (uint8_t *) malloc(current_mod.skip_counts_cap * sizeof(uint8_t));
        MALLOC_CHK(current_mod.probs);

        // allocate initial memory for modification codes
        current_mod.mod_codes_cap = INIT_MODS_CODES;
        current_mod.mod_codes = (char *) malloc(current_mod.mod_codes_cap * sizeof(char));
        MALLOC_CHK(current_mod.mod_codes);

        // set default status flag to '.' (when not present or '.' in the MM string)
        current_mod.status_flag = '.';

        // get base
        current_mod.base = mm_string[i];
        if (!isValidBase(current_mod.base))
            die("Error: Invalid base:%c\n", current_mod.base);
        i++;

        // get strand
        current_mod.strand = mm_string[i];
        if (!isValidStrand(current_mod.strand))
            die("Error: Invalid strand:%c\n", current_mod.strand);
        i++;

        // get base modification codes
        int j = 0;
        while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';' && mm_string[i] != '?' && mm_string[i] != '.') {

            if (j >= current_mod.mod_codes_cap) {
                current_mod.mod_codes_cap *= 2;
                current_mod.mod_codes = (char *) realloc(current_mod.mod_codes, current_mod.mod_codes_cap * sizeof(char));
                MALLOC_CHK(current_mod.mod_codes);
                // die("Error: Too many mod codes than INIT_MODS_CODES=.\n", INIT_MODS_CODES);
            }

            if (!isValidModificationCode(mm_string[i]))
                die("Error: Invalid base modification code:%c\n", mm_string[i]);
            current_mod.mod_codes[j] = mm_string[i];

            i++;
            j++;
        }
        current_mod.num_mod_codes = j;

        // get modification status flag
        if (mm_string[i] == '?' || mm_string[i] == '.') {
            current_mod.status_flag = mm_string[i];
            i++;
        } else if (mm_string[i] == ',') { // if not present, set to '.'
            current_mod.status_flag = '.';
            i++;
        } else if (mm_string[i] == ';') { // end of current mod, no skip counts
            current_mod.status_flag = '.';
        }

        // get skip counts
        int k = 0;
        while (i < mm_str_len && mm_string[i] != ';') {

            if (k >= current_mod.skip_counts_cap) {
                current_mod.skip_counts_cap *= 2;
                current_mod.skip_counts = (int *) realloc(current_mod.skip_counts, current_mod.skip_counts_cap * sizeof(int));
                MALLOC_CHK(current_mod.skip_counts);

                current_mod.probs = (uint8_t *) realloc(current_mod.probs, current_mod.skip_counts_cap * sizeof(uint8_t));
                MALLOC_CHK(current_mod.probs);
                // die("Error: Too many skip counts than INIT_SKIP_COUNTS=%d.\n", INIT_SKIP_COUNTS);
            }


            char skip_count_str[10];
            int l = 0;
            while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';') {
                // printf("mm_string[%d]: %c\n", i, mm_string[i]);
                skip_count_str[l] = mm_string[i];
                i++;
                l++;
            }
            skip_count_str[l] = '\0';
            // printf("skip_count_str: %s\n", skip_count_str);
            current_mod.skip_counts[k] = atoi(skip_count_str);
            current_mod.probs[k] = ml[probs_i++];
            i++;
            k++;
        }
        current_mod.num_skips = k;

        mods[num_mods] = current_mod;
        num_mods++;

        /*
        printf("Base: %c\n", current_mod.base);
        printf("Strand: %c\n", current_mod.strand);
        printf("Modification codes: %s\n", current_mod.mod_codes);
        printf("Status flag: %c\n", current_mod.status_flag);
        printf("Skip counts: ");
        for (int m = 0; m < current_mod.num_skips; m++) {
            printf("%d ", current_mod.skip_counts[m]);
        }
        printf("\n");
        */


    }

    *len_mods = num_mods;

    return mods;

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
    if(len==0){
        ERROR("%s array length is 0 in read %s",tag, bam_get_qname(record));
        exit(EXIT_FAILURE);
    }

    // get the sequence length and check if they match to array len
    // int seq_len = record->core.l_qseq;
    // if(seq_len != len){
    // int seq_len = record->core.l_qseq;
    //     ERROR("%d array length does not match sequence length %d in read %s. Hard clipped alignments? If aligner is minimap2, make sure you use -Y option.", len, seq_len, bam_get_qname(record));
    // }

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

static void print_mods(mod_t *mods, uint32_t len, bam_hdr_t *hdr, bam1_t *record){

    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = tid >= 0 ? hdr->target_name[tid] : ".";
    int32_t pos = record->core.pos;
    int32_t end = bam_endpos(record);
    const char *qname = bam_get_qname(record);

    int8_t rev = bam_is_rev(record);
    const char strand = rev ? '-' : '+';
    assert(!(record->core.flag & BAM_FUNMAP));

    printf("chromosome\tstrand\tstart\tend\tread_name\n");

    for(int i=0;i<len;i++){
        fprintf(stdout, "%s\t%c\t%d\t%d\t%s\n", tname, strand, pos, end, qname);
    }

}

static void print_hdr(){
    printf("chromosome\tstrand\tstart\tend\tread_name\tprob");
}

void simple_meth_view(core_t* core){

    bam1_t *record = bam_init1();
    while(sam_itr_next(core->bam_fp, core->itr, record) >= 0){


        const char *mm = get_mm_tag_ptr(record);
        uint32_t len;
        const uint8_t *ml = get_ml_tag(record, &len);

        uint32_t mods_len = 0;
        mod_t *mods = extract_mods(mm, ml, &mods_len);

        bam_hdr_t *hdr = core->bam_hdr;

        if(ml==NULL){
            continue;
        }

        print_ml_array(ml, len, record);
        print_mm_array(mm, len, record);
        print_mods(mods, mods_len, hdr, record);
        
        //free(ri);
        //free(rp);


    }
    bam_destroy1(record);
    return;
}



