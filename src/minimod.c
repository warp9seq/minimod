/**
 * @file minimod.c
 * @brief common functions for minimod
 * @author Hasindu Gamaarachchi (hasindu@unsw.edu.au)
 *         Suneth Samarasinghe (imsuneth@gmail.com)

MIT License

Copyright (c) 2019 Hasindu Gamaarachchi (hasindu@unsw.edu.au)
Copyright (c) 2024 Suneth Samarasinghe (imsuneth@gmail.com)

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
#include "ref.h"

#include <sys/wait.h>
#include <unistd.h>

/* initialise the core data structure */
core_t* init_core(opt_t opt,double realtime0) {

    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    core->opt = opt;

    //realtime0
    core->realtime0=realtime0;

    core->load_db_time=0;
    core->process_db_time=0;
    core->output_time=0;
    core->merge_db_time=0;

    core->total_bytes=0;
    core->total_reads=0;
    core->processed_reads=0;
    core->processed_bytes=0;

    // load bam file
    core->bam_fp = sam_open(opt.bam_file, "r");
    NULL_CHK(core->bam_fp);

    if(opt.num_thread > 1){
        hts_set_threads(core->bam_fp, opt.num_thread);
    }

    // load bam index file
    core->bam_idx = sam_index_load(core->bam_fp, opt.bam_file);
    if(core->bam_idx==NULL){
        ERROR("could not load the .bai index file for %s", opt.bam_file);
        fprintf(stderr, "Please run 'samtools index %s'\n", opt.bam_file);
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

    if (opt.subtool == FREQ) {
        core->freq_map = kh_init(freqm);
    }
    
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

    if (opt.subtool == FREQ) {
        destroy_freq_map(core->freq_map);
    }

    free(core);
}

/* initialise a data batch */
db_t* init_db(core_t* core) {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->cap_bam_recs = core->opt.batch_size;
    db->n_bam_recs = 0;
    db->processed_bytes=0;
    db->total_reads=0;
    db->total_bytes=0;

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

    if(core->opt.insertions) {
        db->ins = (int**)(malloc(sizeof(int*) * db->cap_bam_recs));
        MALLOC_CHK(db->ins);
        db->ins_offset = (int**)(malloc(sizeof(int*) * db->cap_bam_recs));
        MALLOC_CHK(db->ins_offset);
    }

    if(core->opt.haplotypes){
        db->haplotypes = (int*)(malloc(sizeof(int) * db->cap_bam_recs));
        MALLOC_CHK(db->haplotypes);
    }
    
    db->bases_pos = (int***)(malloc(sizeof(int**) * db->cap_bam_recs));
    MALLOC_CHK(db->bases_pos);
    db->skip_counts = (int**)(malloc(sizeof(int*) * db->cap_bam_recs));
    MALLOC_CHK(db->skip_counts);
    db->mod_codes = (char**)(malloc(sizeof(char*) * db->cap_bam_recs));
    MALLOC_CHK(db->mod_codes);
    db->mod_codes_cap = (uint8_t*)(malloc(sizeof(uint8_t) * db->cap_bam_recs));
    MALLOC_CHK(db->mod_codes_cap);
    db->mod_prob = (int***)(malloc(sizeof(int**) * db->cap_bam_recs));
    MALLOC_CHK(db->mod_prob);



    int32_t i = 0;
    for (i = 0; i < db->cap_bam_recs; ++i) {
        db->bam_recs[i] = bam_init1();
        NULL_CHK(db->bam_recs[i]);

        db->bases_pos[i] = (int**)malloc(sizeof(int*)*N_BASES);
        MALLOC_CHK(db->bases_pos[i]);

        db->mod_codes[i] = (char*)malloc(sizeof(char)*(core->opt.n_mods+1)); //+1 for null terminator
        MALLOC_CHK(db->mod_codes[i]);

        db->mod_codes_cap[i] = core->opt.n_mods;

        db->mod_prob[i] = (int**)malloc(sizeof(int*)*core->opt.n_mods);
        MALLOC_CHK(db->mod_prob[i]);
    }

    db->means = (double*)calloc(db->cap_bam_recs,sizeof(double));
    MALLOC_CHK(db->means);

    return db;
}

/* load a data batch from disk */
ret_status_t load_db(core_t* core, db_t* db) {

    double load_start = realtime();

    // unset previous counts
    db->n_bam_recs = 0;
    db->processed_bytes = 0;
    db->total_reads = 0;
    db->total_bytes = 0;

    ret_status_t status = {0, 0};
    int32_t i;
    bam1_t* rec;

    while (db->n_bam_recs < db->cap_bam_recs && db->processed_bytes < core->opt.batch_size_bases) {
        if (sam_itr_next(core->bam_fp, core->itr, db->bam_recs[db->n_bam_recs]) < 0) {
            break;
        }
        
        i = db->n_bam_recs;
        rec = db->bam_recs[i];

        db->total_reads++;
        db->total_bytes += rec->l_data;
        
        if(rec->core.flag & BAM_FUNMAP){
                LOG_TRACE("Skipping unmapped read %s",bam_get_qname(rec));
                continue;
            }

            if(rec->core.l_qseq == 0){
                LOG_TRACE("Skipping read with 0 length %s",bam_get_qname(rec));
                continue;
            }

        const char *mm = get_mm_tag_ptr(rec);
        if (!mm) {
            LOG_TRACE("Skipping read %s with empty MM tag", bam_get_qname(rec));
            continue;
        }

        uint32_t ml_len;
        uint8_t *ml = get_ml_tag(rec, &ml_len);
        if (!ml) {
            continue;
        }

        db->aln[i] = (int*)malloc(sizeof(int)*rec->core.l_qseq);
        MALLOC_CHK(db->aln[i]);

        if(core->opt.insertions) {
            db->ins[i] = (int*)malloc(sizeof(int)*rec->core.l_qseq);
            MALLOC_CHK(db->ins[i]);
            db->ins_offset[i] = (int*)malloc(sizeof(int)*rec->core.l_qseq);
            MALLOC_CHK(db->ins_offset[i]);
        }

        for(int j=0;j<core->opt.n_mods;j++){
            db->mod_prob[i][j] = (int*)malloc(sizeof(int)*rec->core.l_qseq);
            MALLOC_CHK(db->mod_prob[i][j]);
        }

        for(int j=0;j<N_BASES;j++){
            db->bases_pos[i][j] = (int*)malloc(sizeof(int)*rec->core.l_qseq);
            MALLOC_CHK(db->bases_pos[i][j]);
        }

        db->skip_counts[i] = (int*)malloc(sizeof(int)*rec->core.l_qseq);
        MALLOC_CHK(db->skip_counts[i]);

        db->mm[i] = mm;
        db->ml_lens[i] = ml_len;
        db->ml[i] = ml;

        db->n_bam_recs++;
        db->processed_bytes += rec->l_data;
    }

    status.num_reads = db->n_bam_recs;
    status.num_bases = db->processed_bytes;

    core->load_db_time += realtime() - load_start;

    return status;
}

void work_per_single_read(core_t* core,db_t* db, int32_t i){
    modbases_single(core, db, i);
}

void process_db(core_t* core,db_t* db){
    double proc_start = realtime();

    work_db(core, db, work_per_single_read);

    core->process_db_time += (realtime()-proc_start);
}


/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db) {

    double output_start = realtime();

    print_view_output(core, db);

    core->total_reads += db->total_reads;
    core->total_bytes += db->total_bytes;
    core->processed_reads += db->n_bam_recs;
    core->processed_bytes += db->processed_bytes;

    core->output_time += (realtime()-output_start);

}

void merge_db(core_t* core, db_t* db) {

    double merge_start = realtime();

    update_freq_map(core, db);

    core->total_reads += db->total_reads;
    core->total_bytes += db->total_bytes;
    core->processed_reads += db->n_bam_recs;
    core->processed_bytes += db->processed_bytes;

    core->merge_db_time += (realtime()-merge_start);

}

void output_core(core_t* core) {
    double output_start = realtime();

    if(core->opt.subtool == FREQ){
        print_freq_output(core);
    }

    core->output_time += (realtime()-output_start);
}

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(core_t* core, db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_bam_recs; i++) {        
        free(db->ml[i]);
        free(db->aln[i]);
        if(core->opt.insertions) {
            free(db->ins[i]);
            free(db->ins_offset[i]);
        }
        for(int j=0;j<core->opt.n_mods;j++){
            free(db->mod_prob[i][j]);
        }

        for(int b=0;b<N_BASES;b++){
            free(db->bases_pos[i][b]);
        }

        free(db->skip_counts[i]);

    }
}

/* completely free a data batch */
void free_db(core_t* core, db_t* db) {

    int32_t i = 0;
    
    // free the rest of the records
    for (i = 0; i < db->cap_bam_recs; i++) {
        free(db->mod_prob[i]);
        free(db->mod_codes[i]);
        free(db->bases_pos[i]);
        bam_destroy1(db->bam_recs[i]);
    }

    free(db->skip_counts);
    free(db->mod_codes);
    free(db->mod_codes_cap);
    free(db->mod_prob);
    free(db->bases_pos);
    free(db->ml_lens);
    free(db->mm);
    free(db->ml);
    free(db->aln);
    if(core->opt.insertions) {
        free(db->ins);
        free(db->ins_offset);
    }
    if(core->opt.haplotypes){
        free(db->haplotypes);
    }
    free(db->bam_recs);
    free(db->means);
    free(db);
}

/* initialise user specified options */
void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->batch_size = 512;
    opt->batch_size_bases = 20*1000*1000;
    opt->num_thread = 8;

    opt->debug_break=-1;

    opt->output_fp = stdout;
    opt->progress_interval = 0;
    opt->output_file = NULL;
    opt->bam_file = NULL;
    opt->ref_file = NULL;
    opt->mod_codes_str = NULL;
    opt->mod_threshes_str = NULL;

    opt->modcodes_map = kh_init(modcodesm);

#ifdef HAVE_ACC
    opt->flag |= MINIMOD_ACC;
#endif

}

/* free user specified options */
void free_opt(opt_t* opt) {
    free(opt->mod_threshes_str);
    khint_t i;
    for (i = kh_begin(opt->modcodes_map); i < kh_end(opt->modcodes_map); ++i) {
        if (kh_exist(opt->modcodes_map, i)) {
            modcodem_t *modcode = kh_value(opt->modcodes_map, i);
            char * key = (char*) kh_key(opt->modcodes_map, i);
            free(key);
            free(modcode->context);
            free(modcode);
        }
    }
    kh_destroy(modcodesm, opt->modcodes_map);
}
