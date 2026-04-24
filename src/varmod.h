/**
 * @file mod.h
 * @brief modification tags

MIT License

Copyright (c) 2023 Hasindu Gamaarachchi (hasindu@unsw.edu.au)
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

#ifndef VARMOD_H
#define VARMOD_H

#include "minimod.h"

typedef struct {
    char *key;
    varfreq_t *freq;
} varfreq_kv_t;

typedef struct {
    char *key;
    varview_t *view;
} varview_kv_t;


void varviewfreq_single(core_t * core, db_t *db, int32_t bam_i);
void destroy_var_map(khash_t(varm)* var_map);
void load_var_map(const char* vcf_file, khash_t(varm)* var_map);
void warn_untested_cases_var(opt_t * opt);
void print_varview_header(core_t* core);
void print_varview_output(core_t* core, db_t* db);

#endif
