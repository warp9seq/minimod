/**
 * @file meth_freq.c
 * @brief entry point to meth_freq
 * @author Hasindu Gamaarachchi (hasindu@unsw.edu.au)
 * @author Suneth Samarasinghe (suneth@unsw.edu.au)

MIT License

Copyright (c) 2019 Hasindu Gamaarachchi (hasindu@unsw.edu.au)
Copyright (c) 2024 Suneth Samarasinghe (suneth@unsw.edu.au)

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

#include "minimod.h"
#include "error.h"
#include "misc.h"
#include "meth.h"
#include <assert.h>
#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static struct option long_options[] = {
    {"reference", required_argument, 0, 'r'},      //0 reference genome fasta file [required]
    {"bedmethyl", no_argument, 0, 'b'},            //1 output in bedMethyl format
    {"mod_thresh", required_argument, 0, 'm'},     //2 min modification threshold 0.0 to 1.0 [0.2]
    {"threads", required_argument, 0, 't'},        //3 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //4 batchsize - number of reads loaded at once [512]
    {"max-bytes", required_argument, 0, 'B'},      //5 batchsize - number of bytes loaded at once
    {"verbose", required_argument, 0, 'v'},        //6 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //7
    {"version", no_argument, 0, 'V'},              //8
    {"output",required_argument, 0, 'o'},          //9 output to a file [stdout]
    {"debug-break",required_argument, 0, 0},       //10 break after processing the first batch (used for debugging)
    {"profile-cpu",required_argument, 0, 0},       //11 perform section by section (used for profiling - for CPU only)
    {"accel",required_argument, 0, 0},             //12 accelerator
    {"expand",no_argument, 0, 0},                  //13 expand view
    {0, 0, 0, 0}};


static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: minimod meth-freq reads.bam\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -r FILE                    reference genome fasta file\n");
    fprintf(fp_help,"   -b                         output in bedMethyl format\n");
    fprintf(fp_help,"   -m FLOAT                   min modification threshold (inclusive, range 0.0 to 1.0) [0.2]\n");
    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bytes loaded at once [%.1fM]\n",opt.batch_size_bytes/(float)(1000*1000));
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   -o FILE                    output to file [stdout]\n");
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");

    fprintf(fp_help,"\nadvanced options:\n");
    fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");
    fprintf(fp_help,"   --profile-cpu=yes|no       process section by section (used for profiling on CPU)\n");
#ifdef HAVE_ACC
    fprintf(fp_help,"   --accel=yes|no             Running on accelerator [%s]\n",(opt.flag&minimod_ACC?"yes":"no"));
#endif

}



int meth_freq_main(int argc, char* argv[]) {

    double realtime0 = realtime();

    const char* optstring = "m:r:t:B:K:v:o:hVb";

    int longindex = 0;
    int32_t c = -1;

    const char *bamfile = NULL;
    const char *reffile = NULL;
    int bedmethyl_out = 0;
    double mod_thresh = 0.2;

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c == 'B') {
            opt.batch_size_bytes = mm_parse_num(optarg);
            if(opt.batch_size_bytes<=0){
                ERROR("%s","Maximum number of bytes should be larger than 0.");
                exit(EXIT_FAILURE);
            }
        } else if (c == 'K') {
            opt.batch_size = atoi(optarg);
            if (opt.batch_size < 1) {
                ERROR("Batch size should larger than 0. You entered %d",opt.batch_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 't') {
            opt.num_thread = atoi(optarg);
            if (opt.num_thread < 1) {
                ERROR("Number of threads should larger than 0. You entered %d", opt.num_thread);
                exit(EXIT_FAILURE);
            }
        } else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum log_level_opt)v);
        } else if (c=='V'){
            fprintf(stdout,"minimod %s\n",MINIMOD_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if (c=='r'){
            reffile = optarg;
        } else if (c=='m'){
            mod_thresh = atof(optarg);
        } else if (c=='b'){
            bedmethyl_out = 1;
        }else if(c == 0 && longindex == 10){ //debug break
            opt.debug_break = atoi(optarg);
        } else if(c == 0 && longindex == 11){ //sectional benchmark todo : warning for gpu mode
            yes_or_no(&opt.flag, MINIMOD_PRF, long_options[longindex].name, optarg, 1);
        } else if(c == 0 && longindex == 12){ //accel
        #ifdef HAVE_ACC
            yes_or_no(&opt.flag, minimod_ACC, long_options[longindex].name, optarg, 1);
        #else
            WARNING("%s", "--accel has no effect when compiled for the CPU");
        #endif
        } else if(c == 0 && longindex == 13){ //expand output
            yes_or_no(&opt.flag, MINIMOD_EXP, long_options[longindex].name, "yes", 1);
        }
    }

    // No arguments given
    if (argc - optind != 1 || fp_help == stdout) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    bamfile = argv[optind];

    if (bamfile == NULL) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    if (reffile == NULL) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    //initialise the core data structure
    core_t* core = init_core(bamfile, opt, realtime0);

    init_meth(reffile, mod_thresh);

    meth_freq(core);
    print_stats(stdout, bedmethyl_out);

    destroy_mod();

    //free the core data structure
    free_core(core,opt);

    return 0;
}
