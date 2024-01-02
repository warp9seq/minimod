/**
 * @file pulse.c
 * @brief modification tags

MIT License

Copyright (c) 2023 Hasindu Gamaarachchi (hasindu@unsw.edu.au)

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

#define MAX_MODIFICATIONS 100
#define MAX_SKIP_COUNTS 1000
#define MAX_MODIFICATIONS_CODES 3

typedef struct {
    char base; 
    char strand;
    char modification_codes[MAX_MODIFICATIONS_CODES];
    int skip_counts[MAX_SKIP_COUNTS];
    int num_skips;
    char status_flag;
} Modification;

Modification mods[MAX_MODIFICATIONS];

int isValidBase(char ch) {
    ch = toupper(ch);
    return (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'U' || ch == 'N');
}

int isValidStrand(char ch) {
    return (ch == '+' || ch == '-');
}

int isValidModificationCode(char ch) {
    return (ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'z');
}

int die(const char *format, ...) {
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(EXIT_FAILURE);
}

// Function to extract methylated C base modification regions from MM string
void extractModifications(const char *mm_string) {
    if (mm_string == NULL || strlen(mm_string) == 0) {
        printf("Error: Empty MM string.\n");
        return;
    }
    
    int num_mods = 0;

    int mm_str_len = strlen(mm_string);
    int i = 0;
    while (i < mm_str_len) {

        (num_mods < MAX_MODIFICATIONS) || die("Error: Too many modifications than MAX_MODIFICATIONS=%d.\n", MAX_MODIFICATIONS);

        Modification current_mod;
        memset(&current_mod, 0, sizeof(Modification));

        // set default status flag to '.' (when not present or '.' in the MM string)
        current_mod.status_flag = '.';

        // get base
        current_mod.base = mm_string[i];
        isValidBase(current_mod.base) || die("Error: Invalid base:%c\n", current_mod.base);
        i++;

        // get strand
        current_mod.strand = mm_string[i];
        isValidStrand(current_mod.strand) || die("Error: Invalid strand:%c\n", current_mod.strand);
        i++;

        // get base modification codes
        int j = 0;
        while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';') {

            (j < MAX_MODIFICATIONS_CODES) || die("Error: Too many modification codes than MAX_MODIFICATIONS_CODES=.\n", MAX_MODIFICATIONS_CODES);

            if (j>0 && (mm_string[i] == '?' || mm_string[i] == '.')) {
                current_mod.status_flag = mm_string[i];
                i+=2; // skip the comma
                j++;
                current_mod.modification_codes[j] = '\0';
                break;
            }
            
            isValidModificationCode(mm_string[i]) || die("Error: Invalid base modification code:%c\n", mm_string[i]);
            current_mod.modification_codes[j] = mm_string[i];
            
            i++;
            j++;
        }

        // get skip counts
        int k = 0;
        while (i < mm_str_len && mm_string[i] != ';') {

            (k < MAX_SKIP_COUNTS) || die("Error: Too many skip counts than MAXMAX_SKIP_COUNTS=%d.\n", MAX_SKIP_COUNTS);

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
            i++;
            k++;
        }
        current_mod.num_skips = k;

        mods[num_mods] = current_mod;
        num_mods++;

        /*
        printf("Base: %c\n", current_mod.base);
        printf("Strand: %c\n", current_mod.strand);
        printf("Modification codes: %s\n", current_mod.modification_codes);
        printf("Status flag: %c\n", current_mod.status_flag);
        printf("Skip counts: ");
        for (int m = 0; m < current_mod.num_skips; m++) {
            printf("%d ", current_mod.skip_counts[m]);
        }
        printf("\n");
        */
        

    }

}

void *get_mm_tag(bam1_t *record, uint32_t *len_ptr){

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

    fprintf(stdout, "%s\t%s\t%s\n", bam_get_qname(record), tag, mm_str);

    extractModifications(mm_str);


    return NULL;

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


static void print_ml_array(uint8_t *array, uint32_t len, bam1_t *record){

    fprintf(stdout, "%s\t%d\t%s\t", bam_get_qname(record),record->core.pos, "ML");
    for(int i=0;i<len;i++){
        fprintf(stdout, "%d,", array[i]);
    }
    fprintf(stdout, "\n");

}


void simple_meth_view(core_t* core){

    bam1_t *record = bam_init1();
    while(sam_itr_next(core->bam_fp, core->itr, record) >= 0){

        uint32_t len;
        uint16_t *mm = get_mm_tag(record, &len);
        uint8_t *ml = get_ml_tag(record, &len);



        if(ml==NULL){
            continue;
        }
        //print_mm_array(mm, len, record);
        print_ml_array(ml, len, record);


        //free(mm);
        free(ml);
        //free(ri);
        //free(rp);


    }
    bam_destroy1(record);
    return;
}



