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

#ifndef MINIMOD_H
#define MINIMOD_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include "khash.h"
#include <pthread.h>
#ifndef MINIMOD_VERSION
#define MINIMOD_VERSION "dev"
#endif
/*******************************************************
 * flags related to the user specified options (opt_t) *
 *******************************************************/

#define WORK_STEAL 1 //simple work stealing enabled or not (no work stealing mean no load balancing)
#define STEAL_THRESH 1 //stealing threshold

//set if input, processing and output are not to be interleaved (serial mode) - useful for debugging
// #define IO_PROC_NO_INTERLEAVE 1

#define N_BASES 6 // A, C, G, T, N, U

/* input modification code structure */
typedef struct {
    int index;
    char * context; //context string
    uint8_t thresh; //threshold value (0-255)
} modcodem_t;

/* map of required modification codes to their contexts and thresholds */
KHASH_MAP_INIT_STR(modcodesm, modcodem_t *);

/* frequency map entry */
typedef struct {
    uint32_t n_called;
    uint32_t n_mod;
} freq_t;

typedef struct {
    uint8_t mod_prob; //modification probability (0-255)
    int read_pos; //read position of the base
} view_t;

/* frequency map */
KHASH_MAP_INIT_STR(freqm, freq_t *);

/* view map */
KHASH_MAP_INIT_STR(viewm, view_t *);

enum subtool {VIEW=0, FREQ=1};

/* user specified options */
typedef struct {

    int32_t batch_size;         //max reads loaded at once: K
    int64_t batch_size_bases;   //max bytes loaded at once: B

    int32_t num_thread; //t
    int32_t debug_break;

    char *region_str; //the region string in format chr:start-end

    uint8_t bedmethyl_out; //output in bedMethyl format, only for freq
    char *mod_codes_str;
    char *mod_threshes_str;
    khash_t(modcodesm) *modcodes_map;
    char * bam_file;
    char * ref_file;
    char* output_file;
    FILE* output_fp;
    int progress_interval;

    uint8_t subtool; //0:view, 1:freq

    uint8_t n_mods;
    uint8_t insertions; //is insertions enabled, add ins column to the output
    uint8_t haplotypes; //is haplotypes enabled, add haplotype column to the output

} opt_t;


/* a batch of read data (dynamic data based on the reads) */
typedef struct {
    //bam records
    bam1_t** bam_recs;
    int32_t cap_bam_recs;
    int32_t n_bam_recs; // number of reads processed

    //mod tags
    const char ** mm;
    uint32_t * ml_lens;
    uint8_t ** ml;

    // alignment
    int ** aln; // aln[rec_i][read_pos] = ref_pos
    int ** ins; // ins[rec_i][read_pos] = ins_pos
    int ** ins_offset; // ins_offset[rec_i][read_pos] = ins_offset
    int *** bases_pos; // bases_pos[rec_i][base_i] = read_pos
    int ** skip_counts; // skip_counts[rec_i][read_pos] = skip_count
    char ** mod_codes; // mod_codes[rec_i][mod_i] = mod_code
    uint8_t * mod_codes_cap; // mod_codes_cap[rec_i] = mod_codes_cap

    double *means;

    //stats
    int32_t total_reads; //number of reads in the bam file
    int64_t total_bytes; //number of bytes in the bam file
    int64_t processed_bytes; //number of bytes processed

    khash_t(freqm)** freq_maps; // frequency map per record, only for FREQ subtool
    khash_t(viewm)** view_maps; // view map per record, only for VIEW subtool

} db_t;


/* core data structure (mostly static data throughout the program lifetime) */
typedef struct {

    // options
    opt_t opt;

    // bam file related
    htsFile* bam_fp;
    hts_idx_t* bam_idx;
    bam_hdr_t* bam_hdr;
    hts_itr_t* itr;

    //multi region related
    char **reg_list; //the list of regions
    int64_t reg_n;   //number of regions in list
    int64_t reg_i;   //current region being processed

    //clipping coordinates
    int32_t clip_start;
    int32_t clip_end;

    //realtime0
    double realtime0;

    double load_db_time;
    double process_db_time;
    double merge_db_time;
    double output_time;

    //stats //set by output_db
    uint32_t total_reads; //total number entries in the bam file
    uint64_t total_bytes; //total number of bytes in the bam file
    uint32_t processed_reads; //total number of reads processed
    uint64_t processed_bytes; //total number of bytes processed

    khash_t(freqm)* freq_map;

} core_t;


/* argument wrapper for the multithreaded framework used for data processing */
typedef struct {
    core_t* core;
    db_t* db;
    int32_t starti;
    int32_t endi;
    void (*func)(core_t*,db_t*,int);
    int32_t thread_index;
#ifdef WORK_STEAL
    void *all_pthread_args;
#endif
} pthread_arg_t;

/* argument wrapper for multithreaded framework used for input/processing/output interleaving */
typedef struct {
    core_t* core;
    db_t* db;
    //conditional variable for notifying the processing to the output threads
    pthread_cond_t cond;
    pthread_mutex_t mutex;
    int8_t finished;
} pthread_arg2_t;

/* return status by the load_db - used for termination when all the data is processed */
typedef struct {
    int32_t num_reads;
    int64_t num_bases;
} ret_status_t;

/******************************************
 * function prototype for major functions *
 ******************************************/

/* initialise user specified options */
void init_opt(opt_t* opt);

/* initialise the core data structure */
core_t* init_core(opt_t opt, double realtime0);

/* initialise a data batch */
db_t* init_db(core_t* core);

/* load a data batch from disk */
ret_status_t load_db(core_t* dg, db_t* db);

/* process a single read in the given batch db */
void work_per_single_read(core_t* core,db_t* db, int32_t i);

/* process all reads in the given batch db */
void work_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int));

/* process a data batch */
void process_db(core_t* core, db_t* db);

/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db);

/* merge db data into the core map */
void merge_db(core_t* core, db_t* db);

/* write the output for a all processed data batches */
void output_core(core_t* core);

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(core_t* core, db_t* db);

/* completely free a data batch */
void free_db(core_t* core, db_t* db);

/* free the core data structure */
void free_core(core_t* core,opt_t opt);

/* free user specified options */
void free_opt(opt_t* opt);

#endif
