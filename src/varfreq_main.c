/**
 * @file varfreq_main.c
 * @brief entry point to variant-aware modification frequency
 * @author Suneth Samarasinghe (imsuneth@gmail.com)

MIT License

Copyright (c) 2026 Suneth Samarasinghe (imsuneth@gmail.com)

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
#include "varmod.h"
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
    {"mod_codes", required_argument, 0, 'c'},      //0 modification codes (eg. m, h or mh) [m]
    {"mod_thresh", required_argument, 0, 'm'},     //1 min modification threshold 0.0 to 1.0 [0.8]
    {"threads", required_argument, 0, 't'},        //2 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //3 batchsize - number of reads loaded at once [512]
    {"max-bytes", required_argument, 0, 'B'},      //4 batchsize - number of bytes loaded at once
    {"verbose", required_argument, 0, 'v'},        //5 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //6
    {"version", no_argument, 0, 'V'},              //7
    {"prog-interval",required_argument, 0, 'p'},   //8 progress interval
    {"debug-break",required_argument, 0, 0},       //9 break after processing the first batch (used for debugging)
    {"output",required_argument, 0, 'o'},          //10 output file
    {"insertions",no_argument, 0, 0},              //11 enable modifications in insertions
    {"haplotypes",no_argument, 0, 0},              //12 enable haplotype mode
    {"allow-secondary",no_argument, 0, 0},         //13 enable secondary alignments
    {"skip-supplementary",no_argument, 0, 0},      //14 skip supplementary alignments
    {0, 0, 0, 0}};


static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: minimod varfreq ref.fa reads.bam variants.vcf\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -c STR                     modification code(s) (eg. m, h or mh or as ChEBI) [%s]\n", opt.mod_codes_str==NULL?"m":opt.mod_codes_str);
    fprintf(fp_help,"   -m FLOAT                   min modification threshold(s). Comma separated values for each modification code given in -c [%s]\n", opt.mod_threshes_str==NULL?"0.8":opt.mod_threshes_str);
    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bases loaded at once [%.1fM]\n",opt.batch_size_bases/(float)(1000*1000));
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   -p INT                     print progress every INT seconds (0: per batch) [%d]\n", opt.progress_interval);
    fprintf(fp_help,"   -o FILE                    output file [%s]\n", opt.output_file==NULL?"stdout":opt.output_file);
    fprintf(fp_help,"   --insertions               output modifications in insertions [%s]\n", (opt.insertions?"yes":"no"));
    fprintf(fp_help,"   --haplotypes               output haplotypes [%s]\n", (opt.haplotypes?"yes":"no"));
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");
    fprintf(fp_help,"   --allow-secondary          allow secondary alignments [%s]\n", (opt.allow_secondary?"yes":"no"));
    fprintf(fp_help,"   --skip-supplementary       skip supplementary alignments [%s]\n", (opt.skip_supplementary?"yes":"no"));

    fprintf(fp_help,"\nadvanced options:\n");
    fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");
}

void* pthread_processor_varfreq(void* voidargs) {
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    double realtime0=core->realtime0;

    double realtime_prog = realtime();

    process_db(core, db);

    int32_t skipped_reads = db->total_reads-db->n_bam_recs;
    int64_t skipped_bytes = db->total_bytes-db->processed_bytes;
    if(core->opt.progress_interval<=0 || realtime()-realtime_prog > core->opt.progress_interval){
        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bytes) processed\t%d Entries (%.1fM bytes) skipped\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (db->n_bam_recs), (db->total_bytes)/(1000.0*1000.0),
                skipped_reads,skipped_bytes/(1000.0*1000.0));
        realtime_prog = realtime();
    }

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

void* pthread_post_processor_varfreq(void* voidargs){
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    double realtime0=core->realtime0;

    pthread_mutex_lock(&args->mutex);
    while(args->finished==0){
        pthread_cond_wait(&args->cond, &args->mutex);
    }
    pthread_mutex_unlock(&args->mutex);

    if(get_log_level() > LOG_VERB){
        fprintf(stderr, "[%s::%.3f*%.2f] Signal got by post-processor thread!\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));
    }

    merge_db(core, db);

    int32_t skipped_reads = core->total_reads-core->processed_reads;
    if(skipped_reads>0.9*(int32_t)core->total_reads){
        WARNING("%s","90% of the reads are skipped. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
    }
    if(skipped_reads == (int32_t)core->total_reads){
        ERROR("%s","All reads are skipped. Quitting. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
    }

    free_db_tmp(core, db);
    free_db(core, db);
    free(args);
    pthread_exit(0);
}


int varfreq_main(int argc, char* argv[]) {

    double realtime0 = realtime();

    const char* optstring = "m:c:t:B:K:v:p:o:hV";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt);
    opt.subtool = VARFREQ;

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
        } else if(c == 0 && longindex == 9){
            opt.debug_break = atoi(optarg);
        } else if(c == 0 && longindex == 10){
            FILE *fp = fopen(optarg, "w");
            if (fp == NULL) {
                ERROR("Cannot open file %s for writing", optarg);
                exit(EXIT_FAILURE);
            }
            opt.output_file = optarg;
            opt.output_fp = fp;
        } else if(c == 0 && longindex == 11){
            opt.insertions = 1;
        } else if(c == 0 && longindex == 12){
            opt.haplotypes = 1;
        } else if(c == 0 && longindex == 13){
            opt.allow_secondary = 1;
        } else if(c == 0 && longindex == 14){
            opt.skip_supplementary = 1;
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
    warn_untested_cases_var(&opt);

    if(opt.mod_threshes_str==NULL || strlen(opt.mod_threshes_str)==0){
        INFO("%s", "Modification threshold not provided. Using default threshold 0.8");

        opt.mod_threshes_str = (char *)malloc(opt.n_mods * 4 * sizeof(char) + 1);
        MALLOC_CHK(opt.mod_threshes_str);
        memset(opt.mod_threshes_str,0,opt.n_mods * 4 * sizeof(char)+1);

        char * thresh_str = "0.8";
        for(int i=0;i<opt.n_mods;i++){
            strcat(opt.mod_threshes_str,thresh_str);
            if(i<opt.n_mods-1) strcat(opt.mod_threshes_str,",");
        }
    }

    parse_mod_threshes(&opt);
    print_view_options(&opt);

    if (argc - optind != 3 || fp_help == stdout) {
        WARNING("%s","Missing arguments. Expected: ref.fa reads.bam variants.vcf");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    opt.ref_file = argv[optind];
    opt.bam_file = argv[optind+1];
    const char* vcf_file = argv[optind+2];

    if(access(opt.ref_file, F_OK) == -1) {
        ERROR("Reference file %s does not exist", opt.ref_file);
        exit(EXIT_FAILURE);
    }

    if (access(opt.bam_file, F_OK) == -1) {
        ERROR("BAM file %s does not exist", opt.bam_file);
        exit(EXIT_FAILURE);
    }

    if (access(vcf_file, F_OK) == -1) {
        ERROR("VCF file %s does not exist", vcf_file);
        exit(EXIT_FAILURE);
    }

    double realtime1 = realtime();
    fprintf(stderr, "[%s] Loading reference genome %s\n", __func__, opt.ref_file);
    load_ref(opt.ref_file);
    fprintf(stderr, "[%s] Reference genome loaded in %.3f sec\n", __func__, realtime()-realtime1);

    core_t* core = init_core(opt, realtime0);

    double realtime3 = realtime();
    fprintf(stderr, "[%s] Loading VCF file %s\n", __func__, vcf_file);
    load_var_map(vcf_file, core->var_map);
    fprintf(stderr, "[%s] VCF file loaded in %.3f sec\n", __func__, realtime()-realtime3);

    int32_t counter=0;

    print_varfreq_header(core);

#ifdef IO_PROC_NO_INTERLEAVE

    double realtime_prog = realtime();

    db_t* db = init_db(core);

    ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bases};
    while (status.num_reads >= core->opt.batch_size || status.num_bases>=core->opt.batch_size_bases) {

        status = load_db(core, db);
        process_db(core, db);
        merge_db(core, db);
        free_db_tmp(core, db);

        int32_t skipped_reads = db->total_reads-db->n_bam_recs;
        int64_t skipped_bytes = db->total_bytes-db->processed_bytes;
        if(opt.progress_interval<=0 || realtime()-realtime_prog > opt.progress_interval){
            fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bytes) processed\t%d Entries (%.1fM bytes) skipped\n", __func__,
                    realtime() - realtime0, cputime() / (realtime() - realtime0),
                    (db->n_bam_recs), (db->total_bytes)/(1000.0*1000.0),
                    skipped_reads,skipped_bytes/(1000.0*1000.0));
            realtime_prog = realtime();
        }

        skipped_reads = core->total_reads-core->processed_reads;
        if(skipped_reads>0.9*(int32_t)core->total_reads){
            WARNING("%s","90% of the reads are skipped. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
        }
        if(skipped_reads == (int32_t)core->total_reads){
            ERROR("%s","All reads are skipped. Quitting. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
        }

        if(opt.debug_break==counter) break;
        counter++;
    }

    free_db(core, db);

#else //IO_PROC_INTERLEAVE

    ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bases};
    int8_t first_flag_p=0;
    int8_t first_flag_pp=0;
    pthread_t tid_p;
    pthread_t tid_pp;

    while (status.num_reads >= core->opt.batch_size || status.num_bases>=core->opt.batch_size_bases) {

        db_t* db = init_db(core);
        status = load_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        if(first_flag_p){
            int ret = pthread_join(tid_p, NULL);
            NEG_CHK(ret);
            if(get_log_level() > LOG_VERB){
                fprintf(stderr, "[%s::%.3f*%.2f] Joined to processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_p);
            }
        }
        first_flag_p=1;

        pthread_arg2_t *pt_arg = (pthread_arg2_t*)malloc(sizeof(pthread_arg2_t));
        pt_arg->core=core;
        pt_arg->db=db;
        pthread_cond_init(&pt_arg->cond, NULL);
        pthread_mutex_init(&pt_arg->mutex, NULL);
        pt_arg->finished = 0;

        int ret = pthread_create(&tid_p, NULL, pthread_processor_varfreq, (void*)(pt_arg));
        NEG_CHK(ret);
        if(get_log_level() > LOG_VERB){
            fprintf(stderr, "[%s::%.3f*%.2f] Spawned processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_p);
        }

        if(first_flag_pp){
            int ret = pthread_join(tid_pp, NULL);
            NEG_CHK(ret);
            if(get_log_level() > LOG_VERB){
                fprintf(stderr, "[%s::%.3f*%.2f] Joined to post-processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_pp);
            }
        }
        first_flag_pp=1;

        ret = pthread_create(&tid_pp, NULL, pthread_post_processor_varfreq, (void*)(pt_arg));
        NEG_CHK(ret);
        if(get_log_level() > LOG_VERB){
            fprintf(stderr, "[%s::%.3f*%.2f] Spawned post-processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_pp);
        }

        if(opt.debug_break==counter) break;
        counter++;
    }

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

    destroy_var_map(core->var_map);
    destroy_ref_wo_context(opt.n_mods);

    fprintf(stderr, "[%s] total entries: %ld", __func__,(long)core->total_reads);
    fprintf(stderr,"\n[%s] total bytes: %.1f M",__func__,core->total_bytes/(float)(1000*1000));
    fprintf(stderr,"\n[%s] total skipped entries: %ld",__func__,(long)(core->total_reads-core->processed_reads));
    fprintf(stderr,"\n[%s] total skipped bytes: %.1f M",__func__,(core->total_bytes-core->processed_bytes)/(float)(1000*1000));
    fprintf(stderr,"\n[%s] total processed entries: %ld",__func__,(long)core->processed_reads);
    fprintf(stderr,"\n[%s] total processed bytes: %.1f M",__func__,(core->processed_bytes)/(float)(1000*1000));

    fprintf(stderr, "\n[%s] Data loading time: %.3f sec", __func__,core->load_db_time);
    fprintf(stderr, "\n[%s] Data processing time: %.3f sec", __func__,core->process_db_time);
    fprintf(stderr, "\n[%s] Data merging time: %.3f sec", __func__,core->merge_db_time);
    fprintf(stderr, "\n[%s] Data sorting time: %.3f sec", __func__,core->sort_time);
    fprintf(stderr, "\n[%s] Data output time: %.3f sec", __func__,core->output_time);
    fprintf(stderr,"\n");

    free_core(core,opt);
    free_opt(&opt);

    return 0;
}
