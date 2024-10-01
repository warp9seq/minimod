/**
 * @file mod_freq.c
 * @brief entry point to mod_freq
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
#include "mod.h"
#include "error.h"
#include "misc.h"
#include "ref.h"
#include <assert.h>
#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static struct option long_options[] = {
    {"bedmethyl", no_argument, 0, 'b'},            //0 output in bedMethyl format
    {"mod_codes", required_argument, 0, 'c'},      //1 modification codes (ex. m , h or mh) [m]
    {"mod_thresh", required_argument, 0, 'm'},     //2 min modification threshold 0.0 to 1.0 [0.8]
    {"threads", required_argument, 0, 't'},        //3 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //4 batchsize - number of reads loaded at once [512]
    {"max-bytes", required_argument, 0, 'B'},      //5 batchsize - number of bytes loaded at once
    {"verbose", required_argument, 0, 'v'},        //6 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //7
    {"version", no_argument, 0, 'V'},              //8
    {"prog-interval",required_argument, 0, 'p'},   //9 progress interval
    {"debug-break",required_argument, 0, 0},       //10 break after processing the first batch (used for debugging)
    {"profile-cpu",required_argument, 0, 0},       //11 perform section by section (used for profiling - for CPU only)
    {"accel",required_argument, 0, 0},             //12 accelerator
    {"expand",no_argument, 0, 0},                  //13 expand view
    {"output",required_argument, 0, 'o'},          //14 output file
    {0, 0, 0, 0}};


static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: minimod mod-freq ref.fa reads.bam\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -b                         output in bedMethyl format [%s]\n", (opt.bedmethyl_out?"yes":"not set"));
    fprintf(fp_help,"   -c STR                     modification codes (ex. m , h or mh) [%s]\n", opt.req_mod_codes);
    fprintf(fp_help,"   -m FLOAT                   min modification threshold(s). Comma separated values for each modification code given in -c [%s]\n", opt.req_threshes);
    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bytes loaded at once [%.1fM]\n",opt.batch_size_bytes/(float)(1000*1000));
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   -p INT                     print progress every INT seconds (0: per batch) [%d]\n", opt.progress_interval);
    fprintf(fp_help,"   -o FILE                    output file [%s]\n", opt.output_file==NULL?"stdout":opt.output_file);
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");

    fprintf(fp_help,"\nadvanced options:\n");
    fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");
    fprintf(fp_help,"   --profile-cpu=yes|no       process section by section (used for profiling on CPU)\n");
#ifdef HAVE_ACC
    fprintf(fp_help,"   --accel=yes|no             Running on accelerator [%s]\n",(opt.flag&minimod_ACC?"yes":"no"));
#endif

}

int mod_freq_main(int argc, char* argv[]) {

    double realtime0 = realtime();
    double realtime_prog = realtime();

    const char* optstring = "m:c:t:B:K:v:p:o:hVb";

    int longindex = 0;
    int32_t c = -1;

    char *ref_file = NULL;
    char *bam_file = NULL;
    char *mod_codes_str = NULL;
    char *mod_threshes_str = NULL;

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults
    opt.subtool = MOD_FREQ;

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
        } else if (c=='p'){
            if (atoi(optarg) < 0) {
                ERROR("Progress interval should be 0 or positive. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.progress_interval = atoi(optarg);
        } else if (c=='o'){
            FILE *fp = fopen(optarg, "w");
            if (fp == NULL) {
                ERROR("Cannot open file %s for writing", optarg);
                exit(EXIT_FAILURE);
            }
            opt.output_file = optarg;
            opt.output_fp = fp;
        } else if (c=='V'){
            fprintf(stdout,"minimod %s\n",MINIMOD_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if (c=='m'){
            mod_threshes_str = optarg;
        } else if (c=='c') {
            mod_codes_str = optarg;
        } else if (c=='b'){
            opt.bedmethyl_out = 1;
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
        } else {
            print_help_msg(fp_help, opt);
            if(fp_help == stdout){
                exit(EXIT_SUCCESS);
            }
            exit(EXIT_FAILURE);
        }
    }

    if(mod_codes_str==NULL || strlen(mod_codes_str)==0){
        INFO("%s", "Modification codes not provided. Using default modification code m");
        mod_codes_str = "m";
    }
    
    if(mod_threshes_str==NULL || strlen(mod_threshes_str)==0){
        INFO("%s", "Modification threshold not provided. Using default threshold 0.8");
        mod_threshes_str = "0.8";
    } 
    
    parse_mod_codes(&opt, mod_codes_str);
    parse_mod_threshes(&opt, mod_threshes_str);

    // No arguments given
    if (argc - optind != 2 || fp_help == stdout) {
        WARNING("%s","Missing arguments");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    ref_file = argv[optind];
    bam_file = argv[optind+1];

    if (ref_file == NULL) {
        WARNING("%s","Reference file not provided");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    if (bam_file == NULL) {
        WARNING("%s","BAM file not provided");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    //load the reference genome, get the contexts, and destroy the reference
    double realtime1 = realtime();
    fprintf(stderr, "[%s] Loading reference genome %s\n", __func__, ref_file);
    load_ref(ref_file);
    fprintf(stderr, "[%s] Reference genome loaded in %.3f sec\n", __func__, realtime()-realtime1);

    double realtime2 = realtime();
    fprintf(stderr, "[%s] Loading contexts in reference\n", __func__);
    load_ref_contexts(opt.req_mod_codes, opt.n_mods, opt.req_mod_contexts);
    fprintf(stderr, "[%s] Reference contexts loaded in %.3f sec\n", __func__, realtime()-realtime2);

    destroy_ref_forward();

    //initialise the core data structure
    core_t* core = init_core(bam_file, opt, realtime0);

    int32_t counter=0;

    //initialise a databatch
    db_t* db = init_db(core);

    print_freq_tsv_header(core);

    ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bytes};
    while (status.num_reads >= core->opt.batch_size || status.num_bytes>=core->opt.batch_size_bytes) {

        //load a databatch
        status = load_db(core, db);

        //process the data batch
        process_db(core, db);

        //write the output
        output_db(core, db);

        free_db_tmp(core, db);

        //print progress
        if(opt.progress_interval<=0 || realtime()-realtime_prog > opt.progress_interval){
            fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bytes) processed\t%d Entries (%.1fM bytes) skipped\n", __func__,
                    realtime() - realtime0, cputime() / (realtime() - realtime0),
                    (db->n_bam_recs), (db->sum_bytes)/(1000.0*1000.0),
                    db->skipped_reads,db->skipped_reads_bytes/(1000.0*1000.0));
            realtime_prog = realtime();
        }

        //check if 90% of total reads are skipped
        if(core->skipped_reads>0.9*core->total_reads){
            WARNING("%s","90% of the reads are skipped. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
        }

        if(opt.debug_break==counter){
            break;
        }
        counter++;
    }


    output_core(core);

    destroy_ref(opt.n_mods);

    // free the databatch
    free_db(core, db);

    fprintf(stderr, "[%s] total entries: %ld", __func__,(long)core->total_reads);
    fprintf(stderr,"\n[%s] total bytes: %.1f M",__func__,core->sum_bytes/(float)(1000*1000));
    fprintf(stderr,"\n[%s] total skipped entries: %ld",__func__,(long)core->skipped_reads);
    fprintf(stderr,"\n[%s] total skipped bytes: %.1f M",__func__,core->skipped_reads_bytes/(float)(1000*1000));
    fprintf(stderr,"\n[%s] total processed entries: %ld",__func__,(long)(core->total_reads-core->skipped_reads));
    fprintf(stderr,"\n[%s] total processed bytes: %.1f M",__func__,(core->sum_bytes-core->skipped_reads_bytes)/(float)(1000*1000));

    fprintf(stderr, "\n[%s] Data loading time: %.3f sec", __func__,core->load_db_time);
    fprintf(stderr, "\n[%s] Data processing time: %.3f sec", __func__,core->process_db_time);
    if((core->opt.flag&MINIMOD_PRF)|| core->opt.flag & MINIMOD_ACC){
            fprintf(stderr, "\n[%s]     - Parse time: %.3f sec",__func__, core->parse_time);
            fprintf(stderr, "\n[%s]     - Calc time: %.3f sec",__func__, core->calc_time);
    }
    fprintf(stderr, "\n[%s] Data output time: %.3f sec", __func__,core->output_time);

    fprintf(stderr,"\n");

    //free the core data structure
    free_core(core,opt);

    free_opt(&opt);

    return 0;
}
