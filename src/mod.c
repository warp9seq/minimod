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
#include <ctype.h>

// #include <sys/wait.h>
// #include <unistd.h>

#define MAX_MODIFICATIONS 100
#define MAX_SKIP_COUNTS 256
#define MAX_MODIFICATIONS_TYPES 3

typedef struct {
    char base; 
    char strand;
    char modification_type[MAX_MODIFICATIONS_TYPES];
    int skip_counts[MAX_SKIP_COUNTS];
    int num_skips;
    char status_flag;
} Modification;

Modification mods[MAX_MODIFICATIONS];

// Function to check if a character is a valid base (A, C, G, T, U, N)
int isValidBase(char ch) {
    ch = toupper(ch);
    return (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'U' || ch == 'N');
}

// Function to extract methylated C base modification regions from MM string
void extractModifications(const char *mm_string) {
    if (mm_string == NULL || strlen(mm_string) == 0) {
        printf("Error: Empty MM string.\n");
        return;
    }

    const char * delim_semi = ";";
    const char * delim_comma = ",";
    char *token_semi, *token_comma;

    for(token_semi = strtok(mm_string, delim_semi); token_semi != NULL; token_semi = strtok(NULL, delim_semi)) {
        Modification current_mod;
        memset(&current_mod, 0, sizeof(Modification));

        // get first token
        token_comma = strtok(token_semi, delim_comma);
        if (token_comma == NULL) {
            printf("Error: Invalid modification.\n");
            return;
        }

        if (strlen(token_comma) < 3) {
            printf("Error: Invalid modification.\n");
            return;
        }

        if (token_comma[3] == '?' || token_comma[3] == '.') {
            current_mod.status_flag = token_comma[3];
            token_comma[3] = '\0';
        } else {
            current_mod.status_flag = '\0';
        }

        // scan base, strand and modification type
        int parsed = sscanf(token_comma, "%c%c%s", &current_mod.base, &current_mod.strand, current_mod.modification_type);
        if (parsed != 3) {
            printf("Error: Invalid modification.\n");
            return;
        }

        // check base
        if (isValidBase(current_mod.base) == 0) {
            printf("Error: Invalid base.\n");
            return;
        }

        // check strand
        if (current_mod.strand != '+' && current_mod.strand != '-') {
            printf("Error: Invalid strand.\n");
            return;
        }

        // get skip counts
        token_comma = strtok(NULL, delim_comma);
        while(token_comma != NULL) {

            if (isdigit(token_comma[0]) == 0) {
                printf("Error: Invalid skip count.\n");
                return;
            }

            if (current_mod.num_skips < MAX_SKIP_COUNTS) {
                current_mod.skip_counts[current_mod.num_skips] = atoi(token_comma);
            } else {
                printf("Error: Too many skip counts.\n");
                return;
            }

            token_comma = strtok(NULL, delim_comma);
        }
                
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



