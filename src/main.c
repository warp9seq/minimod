/**
 * @file main.c
 * @brief entry point
 * @author Hasindu Gamaarachchi (hasindu@unsw.edu.au)

MIT License

Copyright (c) 2019 Hasindu Gamaarachchi (hasindu@unsw.edu.au)

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "error.h"
#include "misc.h"
#include "invar.h"

int depth_main(int argc, char* argv[]);

int print_usage(FILE *fp_help){

    fprintf(fp_help,"Usage: invar <command> [options]\n\n");
    fprintf(fp_help,"command:\n");
    fprintf(fp_help,"         depth      do something\n");
    fprintf(fp_help,"         subtool2      do something\n");

    if(fp_help==stderr){
        return(EXIT_FAILURE);
    } else if(fp_help==stdout){
        return(EXIT_SUCCESS);
    } else {
        return(EXIT_FAILURE);
    }

}

int main(int argc, char* argv[]){

    double realtime0 = realtime();

    int ret=1;

    if(argc<2){
        return print_usage(stderr);
    } else if (strcmp(argv[1],"depth")==0){
        ret=depth_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"subtool2")==0){
        ret=depth_main(argc-1, argv+1);
    } else if(strcmp(argv[1],"--version")==0 || strcmp(argv[1],"-V")==0){
        fprintf(stdout,"invar %s\n",INVAR_VERSION);
        exit(EXIT_SUCCESS);
    } else if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0){
        return print_usage(stdout);
    } else{
        fprintf(stderr,"[invar] Unrecognised command %s\n",argv[1]);
        return print_usage(stderr);
    }

    fprintf(stderr,"[%s] Version: %s\n", __func__, INVAR_VERSION);
    fprintf(stderr, "[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU time: %.3f sec; Peak RAM: %.3f GB\n\n",
            __func__, realtime() - realtime0, cputime(),peakrss() / 1024.0 / 1024.0 / 1024.0);

    return ret;
}
