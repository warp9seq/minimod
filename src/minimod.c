/**
 * @file minimod.c
 * @brief common functions for minimod
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

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "minimod.h"
#include "mod.h"
#include "misc.h"
#include "error.h"
#include "khash.h"

#include <sys/wait.h>
#include <unistd.h>

/* initialise the core data structure */
core_t* init_core(const char *bamfilename, opt_t opt,double realtime0) {

    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    core->opt = opt;

    //realtime0
    core->realtime0=realtime0;

    core->load_db_time=0;
    core->process_db_time=0;
    core->output_time=0;

    core->sum_bytes=0;
    core->total_reads=0;
    core->skipped_reads=0;
    core->skipped_reads_bytes=0;

    // load bam file
    core->bam_fp = sam_open(bamfilename, "r");
    NULL_CHK(core->bam_fp);

    // load bam index file
    core->bam_idx = sam_index_load(core->bam_fp, bamfilename);
    if(core->bam_idx==NULL){
        ERROR("could not load the .bai index file for %s", bamfilename);
        fprintf(stderr, "Please run 'samtools index %s'\n", bamfilename);
        exit(EXIT_FAILURE);
    }

    // read the bam header
    core->bam_hdr = sam_hdr_read(core->bam_fp);
    NULL_CHK(core->bam_hdr);

    // If processing a region of the genome, get clipping coordinates
    core->clip_start = -1;
    core->clip_end = -1;

    core->reg_list=NULL; //region list is NULL by default
    core->reg_i=0;
    core->reg_n=0;
    core->itr = NULL;

    if(opt.region_str == NULL){
        core->itr = sam_itr_queryi(core->bam_idx, HTS_IDX_START, 0, 0);
        if(core->itr==NULL){
            ERROR("%s","sam_itr_queryi failed. A problem with the BAM index?");
            exit(EXIT_FAILURE);
        }
    }
    else{
        //determine if .bed
        int region_str_len = strlen(opt.region_str);
        if(region_str_len>=4 && strcmp(&(opt.region_str[region_str_len-4]),".bed")==0 ){
            VERBOSE("Fetching the list of regions from file: %s", opt.region_str);
            WARNING("%s", "Loading region windows from a bed file is an experimental option and not yet throughly tested.");
            WARNING("%s", "When loading windows from a bed file, output is based on reads that are unclipped. Also, there may be repeated entries when regions overlap.");
            int64_t count=0;
            char **reg_l = read_bed_regions(opt.region_str, &count);
            core->reg_list = reg_l;
            core->reg_i = 0;
            core->reg_n = count;
        }
        else{
            VERBOSE("Iterating over region: %s\n", opt.region_str);
            core->itr = sam_itr_querys(core->bam_idx, core->bam_hdr, opt.region_str);
            if(core->itr==NULL){
                ERROR("sam_itr_querys failed. Please check if the region string you entered [%s] is valid",opt.region_str);
                exit(EXIT_FAILURE);
            }
            hts_parse_reg(opt.region_str, &(core->clip_start) , &(core->clip_end));
        }
    }

    if (opt.subtool == MOD_FREQ) {
        core->freq_map = kh_init(freqm);
    }


#ifdef HAVE_ACC
    if (core->opt.flag & MINIMOD_ACC) {
        VERBOSE("%s","Initialising accelator");
    }
#endif


    return core;
}

/* free the core data structure */
void free_core(core_t* core,opt_t opt) {

    if(core->itr){
        sam_itr_destroy(core->itr);
    }
    if(core->reg_list){
        for(int64_t i=0;i<core->reg_n;i++){
            free(core->reg_list[i]);
        }
        free(core->reg_list);
    }

    bam_hdr_destroy(core->bam_hdr);
    hts_idx_destroy(core->bam_idx);
    sam_close(core->bam_fp);

    if (opt.subtool == MOD_FREQ) {
        destroy_freq_map(core->freq_map);
    }

    free(opt.mod_threshes);

#ifdef HAVE_ACC
    if (core->opt.flag & MINIMOD_ACC) {
        VERBOSE("%s","Freeing accelator");
    }
#endif

    free(core);
}

/* initialise a data batch */
db_t* init_db(core_t* core) {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->cap_bam_recs = core->opt.batch_size;
    db->n_bam_recs = 0;

    db->bam_recs = (bam1_t**)(malloc(sizeof(bam1_t*) * db->cap_bam_recs));
    MALLOC_CHK(db->bam_recs);

    db->mm = (const char**)(malloc(sizeof(char*) * db->cap_bam_recs));
    MALLOC_CHK(db->mm);
    db->ml_lens = (uint32_t*)(malloc(sizeof(uint32_t) * db->cap_bam_recs));
    MALLOC_CHK(db->ml_lens);
    db->ml = (uint8_t**)(malloc(sizeof(uint8_t*) * db->cap_bam_recs));
    MALLOC_CHK(db->ml);
    db->aln = (int**)(malloc(sizeof(int*) * db->cap_bam_recs));
    MALLOC_CHK(db->aln);
    db->bases_pos = (int***)(malloc(sizeof(int**) * db->cap_bam_recs));
    MALLOC_CHK(db->bases_pos);
    db->bases_pos_len = (int**)(malloc(sizeof(int*) * db->cap_bam_recs));
    MALLOC_CHK(db->bases_pos_len);
    db->skip_counts = (int**)(malloc(sizeof(int*) * db->cap_bam_recs));
    MALLOC_CHK(db->skip_counts);
    db->skip_counts_len = (int*)(calloc(db->cap_bam_recs,sizeof(int)));
    MALLOC_CHK(db->skip_counts_len);
    db->mod_codes = (char**)(malloc(sizeof(char*) * db->cap_bam_recs));
    MALLOC_CHK(db->mod_codes);


    int32_t i = 0;
    for (i = 0; i < db->cap_bam_recs; ++i) {
        db->bam_recs[i] = bam_init1();
        NULL_CHK(db->bam_recs[i]);

        db->bases_pos[i] = (int**)malloc(sizeof(int*)*N_BASES);
        MALLOC_CHK(db->bases_pos[i]);

        db->bases_pos_len[i] = (int*)malloc(sizeof(int)*N_BASES);
        MALLOC_CHK(db->bases_pos_len[i]);

        for(int b=0;b<N_BASES;b++){
            db->bases_pos_len[i][b] = 0;
        }

        db->mod_codes[i] = (char*)malloc(sizeof(char)*N_MODS);
        MALLOC_CHK(db->mod_codes[i]);
    }

    db->modbases = (modbase_t**)malloc(sizeof(modbase_t*) * db->cap_bam_recs);
    MALLOC_CHK(db->modbases);

    db->modbases_len = (int32_t*)calloc(db->cap_bam_recs,sizeof(int32_t));
    MALLOC_CHK(db->modbases_len);
    
    db->means = (double*)calloc(db->cap_bam_recs,sizeof(double));
    MALLOC_CHK(db->means);

    db->sum_bytes=0;
    db->skipped_reads=0;
    db->skipped_reads_bytes=0;

    return db;
}

/* load a data batch from disk */
ret_status_t load_db(core_t* core, db_t* db) {

    double load_start = realtime();

    db->n_bam_recs = 0;
    db->sum_bytes = 0;
    db->skipped_reads = 0;
    db->skipped_reads_bytes = 0;

    ret_status_t status={0,0};
    while (db->n_bam_recs < db->cap_bam_recs && db->sum_bytes<core->opt.batch_size_bytes) {
        if(sam_itr_next(core->bam_fp, core->itr, db->bam_recs[db->n_bam_recs])>=0){
            bam1_t* rec = db->bam_recs[db->n_bam_recs];

            db->sum_bytes += rec->l_data;
            
            if(rec->core.flag & BAM_FUNMAP){
                db->skipped_reads++;
                db->skipped_reads_bytes += rec->l_data;
                WARNING("Skipping unmapped read %s",bam_get_qname(rec));
                continue;
            }

            if(rec->core.l_qseq == 0){
                db->skipped_reads++;
                db->skipped_reads_bytes += rec->l_data;
                WARNING("Skipping read with 0 length %s",bam_get_qname(rec));
                continue;
            }

            const char *mm = get_mm_tag_ptr(rec);

            if(mm == NULL){
                db->skipped_reads++;
                db->skipped_reads_bytes += rec->l_data;
                WARNING("Skipping read %s with empty MM tag", bam_get_qname(rec));
                continue;
            }

            uint32_t ml_len;
            uint8_t *ml = get_ml_tag(rec, &ml_len);

            if(ml == NULL){
                db->skipped_reads++;
                db->skipped_reads_bytes += rec->l_data;
                continue;
            }

            if(db->modbases_len[db->n_bam_recs] == 0) {
                db->modbases[db->n_bam_recs] = (modbase_t*)malloc(sizeof(modbase_t)*rec->core.l_qseq);
                MALLOC_CHK(db->modbases[db->n_bam_recs]);
                db->modbases_len[db->n_bam_recs] = rec->core.l_qseq;

                for(int seq_i=0;seq_i<rec->core.l_qseq;seq_i++){
                    db->modbases[db->n_bam_recs][seq_i].mods = (mod_t*)malloc(sizeof(mod_t)*N_MODS);
                    MALLOC_CHK(db->modbases[db->n_bam_recs][seq_i].mods);
                }
            } else if(db->modbases_len[db->n_bam_recs] < rec->core.l_qseq){
                db->modbases[db->n_bam_recs] = (modbase_t*)realloc(db->modbases[db->n_bam_recs],sizeof(modbase_t)*rec->core.l_qseq);
                MALLOC_CHK(db->modbases[db->n_bam_recs]);

                for(int seq_i=db->modbases_len[db->n_bam_recs];seq_i<rec->core.l_qseq;seq_i++){
                    db->modbases[db->n_bam_recs][seq_i].mods = (mod_t*)malloc(sizeof(mod_t)*N_MODS);
                    MALLOC_CHK(db->modbases[db->n_bam_recs][seq_i].mods);
                }
                db->modbases_len[db->n_bam_recs] = rec->core.l_qseq;
            }

            for(int i=0;i<N_BASES;i++){
                if(db->bases_pos_len[db->n_bam_recs][i] == 0){
                    db->bases_pos[db->n_bam_recs][i] = (int*)malloc(sizeof(int)*rec->core.l_qseq);
                    MALLOC_CHK(db->bases_pos[db->n_bam_recs][i]);
                    db->bases_pos_len[db->n_bam_recs][i] = rec->core.l_qseq;
                }else if(db->bases_pos_len[db->n_bam_recs][i] < rec->core.l_qseq){
                    db->bases_pos[db->n_bam_recs][i] = (int*)realloc(db->bases_pos[db->n_bam_recs][i],sizeof(int)*rec->core.l_qseq);
                    MALLOC_CHK(db->bases_pos[db->n_bam_recs][i]);
                    db->bases_pos_len[db->n_bam_recs][i] = rec->core.l_qseq;
                }
            }
            
            if(db->skip_counts_len[db->n_bam_recs] == 0){
                db->skip_counts[db->n_bam_recs] = (int*)malloc(sizeof(int)*rec->core.l_qseq);
                MALLOC_CHK(db->skip_counts[db->n_bam_recs]);
                db->skip_counts_len[db->n_bam_recs] = rec->core.l_qseq;
            }else if(db->skip_counts_len[db->n_bam_recs] < rec->core.l_qseq){
                db->skip_counts[db->n_bam_recs] = (int*)realloc(db->skip_counts[db->n_bam_recs],sizeof(int)*rec->core.l_qseq);
                MALLOC_CHK(db->skip_counts[db->n_bam_recs]);
                db->skip_counts_len[db->n_bam_recs] = rec->core.l_qseq;
            }

            db->mm[db->n_bam_recs] = mm;
            db->ml_lens[db->n_bam_recs] = ml_len;
            db->ml[db->n_bam_recs] = ml;
            db->aln[db->n_bam_recs] = get_aln(core->bam_hdr, rec);
            db->n_bam_recs++;
        }else{
            break;
        }
    }

    status.num_reads=db->n_bam_recs;
    status.num_bytes=db->sum_bytes;

    double load_end = realtime();
    core->load_db_time += (load_end-load_start);

    return status;
}


void parse_single(core_t* core,db_t* db, int32_t i){



}

#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))


void work_per_single_read(core_t* core,db_t* db, int32_t i){

    if(core->opt.subtool==VIEW || core->opt.subtool==MOD_FREQ){
        modbases_single(core,db,i);
    }else{
        ERROR("Unknown subtool %d", core->opt.subtool);
    }

}

void mean_db(core_t* core, db_t* db) {
#ifdef HAVE_ACC
    if (core->opt.flag & MINIMOD_ACC) {
        VERBOSE("%s","Aligning reads with accel");
        work_db(core,db,mean_single);
    }
#endif

    if (!(core->opt.flag & MINIMOD_ACC)) {
        //fprintf(stderr, "cpu\n");
        // work_db(core,db,mean_single);
    }
}


void process_db(core_t* core,db_t* db){
    double proc_start = realtime();

    if(core->opt.flag & MINIMOD_PRF || core->opt.flag & MINIMOD_ACC){
        double a = realtime();
        work_db(core,db,parse_single);
        double b = realtime();
        core->parse_time += (b-a);

        a = realtime();
        mean_db(core,db);
        b = realtime();
        core->calc_time += (b-a);

    } else {
        work_db(core, db, work_per_single_read);
    }

    double proc_end = realtime();
    core->process_db_time += (proc_end-proc_start);
}


/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db) {

    double output_start = realtime();
    if(core->opt.subtool == VIEW){
        print_view_output(core, db);
    } else if(core->opt.subtool == MOD_FREQ){
        update_freq_map(core, db);
    }

    core->sum_bytes += db->sum_bytes;
    core->total_reads += db->n_bam_recs;
    core->skipped_reads += db->skipped_reads;
    core->skipped_reads_bytes += db->skipped_reads_bytes;

    //core->read_index = core->read_index + db->n_bam_recs;
    double output_end = realtime();
    core->output_time += (output_end-output_start);

}

void output_core(core_t* core) {
    double output_start = realtime();

    if(core->opt.subtool == MOD_FREQ){
        print_freq_output(core);
    }

    double output_end = realtime();
    core->output_time += (output_end-output_start);
}

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_bam_recs; i++) {        
        free(db->ml[i]);
        free(db->aln[i]);
    }
}

/* completely free a data batch */
void free_db(core_t* core, db_t* db) {

    int32_t i = 0;
    
    // free the rest of the records
    for (i = 0; i < db->cap_bam_recs; i++) {
        if(db->skip_counts_len[i]>0){
            free(db->skip_counts[i]);
        }
        if(db->modbases_len[i]>0){
            for(int seq_i=0;seq_i<db->modbases_len[i];seq_i++){
                free(db->modbases[i][seq_i].mods);
            }
            free(db->modbases[i]);
        }
        
        for(int b=0;b<N_BASES;b++){
            if(db->bases_pos_len[i][b]>0){
                free(db->bases_pos[i][b]);
            }
        }

        free(db->mod_codes[i]);

        free(db->bases_pos[i]);
        free(db->bases_pos_len[i]);
        bam_destroy1(db->bam_recs[i]);
    }

    free(db->skip_counts);
    free(db->skip_counts_len);
    free(db->mod_codes);
    free(db->modbases);
    free(db->modbases_len);
    free(db->bases_pos);
    free(db->bases_pos_len);
    free(db->ml_lens);
    free(db->mm);
    free(db->ml);
    free(db->aln);
    free(db->bam_recs);
    free(db->means);
    free(db);
}

/* initialise user specified options */
void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->batch_size = 512;
    opt->batch_size_bytes = 20*1000*1000;
    opt->num_thread = 8;

    opt->debug_break=-1;

    opt->mod_codes = "m";
    opt->output_fp = stdout;

#ifdef HAVE_ACC
    opt->flag |= MINIMOD_ACC;
#endif

}
