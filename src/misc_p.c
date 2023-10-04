/* @file misc.c
**
** miscellaneous definitions and inline functions
** @@
******************************************************************************/
#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "error.h"

//read bed file
char **read_bed_regions(char *bedfile, int64_t *count){

    FILE *bedfp = fopen(bedfile,"r");
    F_CHK(bedfp,bedfile);

    char* buffer = (char*)malloc(sizeof(char) * (100)); //READ+newline+nullcharacter
    MALLOC_CHK(buffer);

    int64_t reg_capcacity = 1024;
    int64_t reg_i = 0;
    char **reg_list = (char **)malloc(reg_capcacity * sizeof(char *));
    MALLOC_CHK(reg_list);


    size_t bufferSize = 100;
    ssize_t readlinebytes = 0;
    int64_t line_no = 0;


    while ((readlinebytes = getline(&buffer, &bufferSize, bedfp)) != -1) {

        char *ref = (char *)malloc(sizeof(char)*readlinebytes);
        MALLOC_CHK(ref);
        int64_t beg=-1;
        int64_t end=-1;

        //TODO can optimised though strtok etc later
        int ret=sscanf(buffer,"%s\t%ld\t%ld",ref,&beg, &end);
        if(ret!=3 || end<beg){
            ERROR("Malformed bed entry at line %ld",line_no);
            exit(EXIT_FAILURE);
        }

        if(reg_i>=reg_capcacity){
            if(reg_capcacity>1000000){
                WARNING("The region bed file has over %ld regions. To reduce memory usage, you may consider merging bed regions.",reg_i);
            }
            reg_capcacity=reg_capcacity*2;
            reg_list = (char **)realloc((void *)reg_list,reg_capcacity * sizeof(char *));
            MALLOC_CHK(reg_list);

        }

        reg_list[reg_i] = (char *)malloc(sizeof(char)*readlinebytes);
        sprintf(reg_list[reg_i],"%s:%ld-%ld",ref, beg, end);
        reg_i++;


        free(ref);



        line_no++;
    }

    fclose(bedfp);
    free(buffer);
    *count = reg_i;

    return reg_list;
}
