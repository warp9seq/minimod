/**
 * @file freq.c
 * @brief entry point to freq
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
    {"mod_codes", required_argument, 0, 'c'},      //1 modification codes (eg. m, h or mh) [m]
    {"mod_thresh", required_argument, 0, 'm'},     //2 min modification threshold 0.0 to 1.0 [0.8]
    {"threads", required_argument, 0, 't'},        //3 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //4 batchsize - number of reads loaded at once [512]
    {"max-bytes", required_argument, 0, 'B'},      //5 batchsize - number of bytes loaded at once
    {"verbose", required_argument, 0, 'v'},        //6 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //7
    {"version", no_argument, 0, 'V'},              //8
    {"prog-interval",required_argument, 0, 'p'},   //9 progress interval
    {"debug-break",required_argument, 0, 0},       //10 break after processing the first batch (used for debugging)
    {"output",required_argument, 0, 'o'},          //11 output file
    {"insertions",no_argument, 0, 0},              //12 enable modifications in insertions
    {"haplotypes",no_argument, 0, 0},              //13 enable haplotype mode
    {0, 0, 0, 0}};


static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: minimod freq ref.fa reads.bam\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -b                         output in bedMethyl format [%s]\n", (opt.bedmethyl_out?"yes":"not set"));
    fprintf(fp_help,"   -c STR                     modification code(s) (eg. m, h or mh or as ChEBI) [%s]\n", opt.mod_codes_str);
    fprintf(fp_help,"   -m FLOAT                   min modification threshold(s). Comma separated values for each modification code given in -c [%s]\n", opt.mod_threshes_str);
    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bases loaded at once [%.1fM]\n",opt.batch_size_bases/(float)(1000*1000));
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   -p INT                     print progress every INT seconds (0: per batch) [%d]\n", opt.progress_interval);
    fprintf(fp_help,"   -o FILE                    output file [%s]\n", opt.output_file==NULL?"stdout":opt.output_file);
    fprintf(fp_help,"   --insertions               enable modifications in insertions [%s]\n", (opt.insertions?"yes":"no"));
    fprintf(fp_help,"   --haplotypes               enable haplotype mode [%s]\n", (opt.haplotypes?"yes":"no"));
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");

    fprintf(fp_help,"\nadvanced options:\n");
    fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");

}

//function that processes a databatch - for pthreads when I/O and processing are interleaved
void* pthread_processor(void* voidargs) {
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    double realtime0=core->realtime0;

    double realtime_prog = realtime();

    //process
    process_db(core, db);

    //print progress
    int32_t skipped_reads = db->total_reads-db->n_bam_recs;
    int64_t skipped_bytes = db->total_bytes-db->processed_bytes;
    if(core->opt.progress_interval<=0 || realtime()-realtime_prog > core->opt.progress_interval){
        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bytes) processed\t%d Entries (%.1fM bytes) skipped\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (db->n_bam_recs), (db->total_bytes)/(1000.0*1000.0),
                skipped_reads,skipped_bytes/(1000.0*1000.0));
        realtime_prog = realtime();
    }

    //need to inform the output thread that we completed the processing
    pthread_mutex_lock(&args->mutex);
    pthread_cond_signal(&args->cond);
    args->finished=1;
    pthread_mutex_unlock(&args->mutex);

    if(get_log_level() > LOG_VERB){
        fprintf(stderr, "[%s::%.3f*%.2f] Signal sent by processor thread!\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));
    }

    pthread_exit(0);
}

//function that prints the output and free - for pthreads when I/O and processing are interleaved
void* pthread_post_processor(void* voidargs){
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    double realtime0=core->realtime0;

    //wait until the processing thread has informed us
    pthread_mutex_lock(&args->mutex);
    while(args->finished==0){
        pthread_cond_wait(&args->cond, &args->mutex);
    }
    pthread_mutex_unlock(&args->mutex);

    if(get_log_level() > LOG_VERB){
        fprintf(stderr, "[%s::%.3f*%.2f] Signal got by post-processor thread!\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));
    }

    //output and free
    merge_db(core, db);

    //check if 90% of total reads are skipped
    int32_t skipped_reads = core->total_reads-core->processed_reads;
    if(skipped_reads>0.9*core->total_reads){
        WARNING("%s","90% of the reads are skipped. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
    }
    if(skipped_reads == core->total_reads){
        ERROR("%s","All reads are skipped. Quitting. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
    }

    free_db_tmp(core, db);
    free_db(core, db);
    free(args);
    pthread_exit(0);
}

int freq_main(int argc, char* argv[]) {

    double realtime0 = realtime();

    const char* optstring = "m:c:t:B:K:v:p:o:hVb";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults
    opt.subtool = FREQ;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c == 'B') {
            opt.batch_size_bases = mm_parse_num(optarg);
            if(opt.batch_size_bases<=0){
                ERROR("%s","Maximum number of bases should be larger than 0.");
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
            opt.mod_threshes_str = (char *)malloc(strlen(optarg)+1);
            MALLOC_CHK(opt.mod_threshes_str);
            strcpy(opt.mod_threshes_str,optarg);
        } else if (c=='c') {
            opt.mod_codes_str = optarg;
        } else if (c=='b'){
            opt.bedmethyl_out = 1;
        }else if(c == 0 && longindex == 10){ //debug break
            opt.debug_break = atoi(optarg);
        }else if(c == 0 && longindex == 11){ //output file
            FILE *fp = fopen(optarg, "w");
            if (fp == NULL) {
                ERROR("Cannot open file %s for writing", optarg);
                exit(EXIT_FAILURE);
            }
            opt.output_file = optarg;
            opt.output_fp = fp;
        } else if(c == 0 && longindex == 12){ //insertions
            opt.insertions = 1;
        } else if(c == 0 && longindex == 13){ //haplotypes
            opt.haplotypes = 1;
        } else {
            print_help_msg(fp_help, opt);
            if(fp_help == stdout){
                exit(EXIT_SUCCESS);
            }
            exit(EXIT_FAILURE);
        }
    }

    if(opt.mod_codes_str==NULL || strlen(opt.mod_codes_str)==0){
        INFO("%s", "Modification codes not provided. Using default modification code m");
        opt.mod_codes_str = "m";
    }

    parse_mod_codes(&opt);

    if(opt.mod_threshes_str==NULL || strlen(opt.mod_threshes_str)==0){
        INFO("%s", "Modification threshold not provided. Using default threshold 0.8");
        
        opt.mod_threshes_str = (char *)malloc(opt.n_mods * 4 * sizeof(char) + 1);
        MALLOC_CHK(opt.mod_threshes_str);
        memset(opt.mod_threshes_str,0,opt.n_mods * 4 * sizeof(char)+1);

        char * thresh_str = "0.8";

        for(int i=0;i<opt.n_mods;i++){
            strcat(opt.mod_threshes_str,thresh_str);
            if(i<opt.n_mods-1){
                strcat(opt.mod_threshes_str,",");
            }
        }
    }
    
    parse_mod_threshes(&opt);
    
    // No arguments given
    if (argc - optind != 2 || fp_help == stdout) {
        WARNING("%s","Missing arguments");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    opt.ref_file = argv[optind];
    opt.bam_file = argv[optind+1];

    if (opt.ref_file == NULL) {
        WARNING("%s","Reference file not provided");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    if (opt.bam_file == NULL) {
        WARNING("%s","BAM file not provided");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    // check if the bam file exists
    if (access(opt.bam_file, F_OK) == -1) {
        ERROR("BAM file %s does not exist", opt.bam_file);
        exit(EXIT_FAILURE);
    }

    //load the reference genome, get the contexts, and destroy the reference
    double realtime1 = realtime();
    fprintf(stderr, "[%s] Loading reference genome %s\n", __func__, opt.ref_file);
    load_ref(opt.ref_file);
    fprintf(stderr, "[%s] Reference genome loaded in %.3f sec\n", __func__, realtime()-realtime1);

    double realtime2 = realtime();
    fprintf(stderr, "[%s] Loading contexts in reference\n", __func__);
    
    char** mod_contexts = (char**)malloc(opt.n_mods * sizeof(char*));
    MALLOC_CHK(mod_contexts);
    for (khint_t i = kh_begin(opt.modcodes_map); i < kh_end(opt.modcodes_map); ++i) {
        if (!kh_exist(opt.modcodes_map, i)) continue;
        modcodem_t *mod_code_map = kh_value(opt.modcodes_map, i);
        mod_contexts[mod_code_map->index] = mod_code_map->context;
    }
    load_ref_contexts(opt.n_mods, mod_contexts);
    free(mod_contexts);
    fprintf(stderr, "[%s] Reference contexts loaded in %.3f sec\n", __func__, realtime()-realtime2);

    destroy_ref_forward();

    //initialise the core data structure
    core_t* core = init_core(opt, realtime0);

    int32_t counter=0;

    print_freq_header(core);

#ifdef IO_PROC_NO_INTERLEAVE

    double realtime_prog = realtime();

    //initialise a databatch
    db_t* db = init_db(core);

    ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bases};
    while (status.num_reads >= core->opt.batch_size || status.num_bases>=core->opt.batch_size_bases) {

        //load a databatch
        status = load_db(core, db);

        //process the data batch
        process_db(core, db);

        //merge into map
        merge_db(core, db);

        free_db_tmp(core, db);

        //print progress
        int32_t skipped_reads = db->total_reads-db->n_bam_recs;
        int64_t skipped_bytes = db->total_bytes-db->processed_bytes;
        if(opt.progress_interval<=0 || realtime()-realtime_prog > opt.progress_interval){
            fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bytes) processed\t%d Entries (%.1fM bytes) skipped\n", __func__,
                    realtime() - realtime0, cputime() / (realtime() - realtime0),
                    (db->n_bam_recs), (db->total_bytes)/(1000.0*1000.0),
                    skipped_reads,skipped_bytes/(1000.0*1000.0));
            realtime_prog = realtime();
        }

        //check if 90% of total reads are skipped
        skipped_reads = core->total_reads-core->processed_reads;
        if(skipped_reads>0.9*core->total_reads){
            WARNING("%s","90% of the reads are skipped. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
        }
        if(skipped_reads == core->total_reads){
            ERROR("%s","All reads are skipped. Quitting. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
        }
        

        if(opt.debug_break==counter){
            break;
        }
        counter++;
    }

    free_db(core, db);

#else //IO_PROC_INTERLEAVE

    ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bases};
    int8_t first_flag_p=0;
    int8_t first_flag_pp=0;
    pthread_t tid_p; //process thread
    pthread_t tid_pp; //post-process thread

    while (status.num_reads >= core->opt.batch_size || status.num_bases>=core->opt.batch_size_bases) {

        //init and load a databatch
        db_t* db = init_db(core);
        status = load_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        if(first_flag_p){ //if not the first time of the "process" wait for the previous "process"
            int ret = pthread_join(tid_p, NULL);
            NEG_CHK(ret);
            if(get_log_level() > LOG_VERB){
                fprintf(stderr, "[%s::%.3f*%.2f] Joined to processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_p);
            }
        }
        first_flag_p=1;

        //set up args
        pthread_arg2_t *pt_arg = (pthread_arg2_t*)malloc(sizeof(pthread_arg2_t));
        pt_arg->core=core;
        pt_arg->db=db;
        pthread_cond_init(&pt_arg->cond, NULL);
        pthread_mutex_init(&pt_arg->mutex, NULL);
        pt_arg->finished = 0;

        //process thread launch
        int ret = pthread_create(&tid_p, NULL, pthread_processor,
                                (void*)(pt_arg));
        NEG_CHK(ret);
        if(get_log_level() > LOG_VERB){
            fprintf(stderr, "[%s::%.3f*%.2f] Spawned processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_p);
        }

        if(first_flag_pp){ //if not the first time of the post-process wait for the previous post-process
            int ret = pthread_join(tid_pp, NULL);
            NEG_CHK(ret);
            if(get_log_level() > LOG_VERB){
                fprintf(stderr, "[%s::%.3f*%.2f] Joined to post-processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_pp);
            }
        }
        first_flag_pp=1;

        //post-process thread launch (output and freeing thread)
        ret = pthread_create(&tid_pp, NULL, pthread_post_processor,
                                (void*)(pt_arg));
        NEG_CHK(ret);
        if(get_log_level() > LOG_VERB){
            fprintf(stderr, "[%s::%.3f*%.2f] Spawned post-processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_pp);
        }

        if(opt.debug_break==counter){
            break;
        }
        counter++;
    }

    //final round
    int ret = pthread_join(tid_p, NULL);
    NEG_CHK(ret);
    if(get_log_level() > LOG_VERB){
        fprintf(stderr, "[%s::%.3f*%.2f] Joined to last processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_p);
    }
    ret = pthread_join(tid_pp, NULL);
    NEG_CHK(ret);
    if(get_log_level() > LOG_VERB){
    fprintf(stderr, "[%s::%.3f*%.2f] Joined to last post-processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_pp);
    }

#endif

    output_core(core);

    destroy_ref(opt.n_mods);

    fprintf(stderr, "[%s] total entries: %ld", __func__,(long)core->total_reads);
    fprintf(stderr,"\n[%s] total bytes: %.1f M",__func__,core->total_bytes/(float)(1000*1000));
    fprintf(stderr,"\n[%s] total skipped entries: %ld",__func__,(long)(core->total_reads-core->processed_reads));
    fprintf(stderr,"\n[%s] total skipped bytes: %.1f M",__func__,(core->total_bytes-core->processed_bytes)/(float)(1000*1000));
    fprintf(stderr,"\n[%s] total processed entries: %ld",__func__,(long)core->processed_reads);
    fprintf(stderr,"\n[%s] total processed bytes: %.1f M",__func__,(core->processed_bytes)/(float)(1000*1000));

    fprintf(stderr, "\n[%s] Data loading time: %.3f sec", __func__,core->load_db_time);
    fprintf(stderr, "\n[%s] Data processing time: %.3f sec", __func__,core->process_db_time);
    fprintf(stderr, "\n[%s] Data merging time: %.3f sec", __func__,core->merge_db_time);
    fprintf(stderr, "\n[%s] Data output time: %.3f sec", __func__,core->output_time);

    fprintf(stderr,"\n");

    //free the core data structure
    free_core(core,opt);

    free_opt(&opt);

    return 0;
}
