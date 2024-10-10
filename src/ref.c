/**
 * @file ref.c
 * @brief reference genome loading

MIT License

Copyright (c) 2024 Suneth Samarasinghe

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

#include <zlib.h>
#include <string.h>

#include "ref.h"
#include "error.h"
#include "kseq.h"
#include "khash.h"

KSEQ_INIT(gzFile, gzread);
KHASH_MAP_INIT_STR(refm, ref_t *);
khash_t(refm)* ref_map;

void load_ref(const char * genome) {
    gzFile fp;
    kseq_t *seq;
    int l;
    
    fp = gzopen(genome, "r");
    F_CHK(fp, genome);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    ref_map = kh_init(refm);
    
    while ((l = kseq_read(seq)) >= 0) {
        ASSERT_MSG(l == (int) seq->seq.l, "Sequence length mismatch: %d vs %d", l, (int) seq->seq.l);

        // initialize ref
        ref_t * ref = (ref_t *) malloc(sizeof(ref_t));
        MALLOC_CHK(ref);
        char * ref_name = (char *) malloc(strlen(seq->name.s) + 1);
        MALLOC_CHK(ref_name);
        strcpy(ref_name, seq->name.s);
        ref->ref_seq_length = seq->seq.l;
        ref->forward = (char *) malloc(seq->seq.l + 1);
        MALLOC_CHK(ref->forward);
        strcpy(ref->forward, seq->seq.s);

        int ret;
        khiter_t k = kh_put(refm, ref_map, ref_name, &ret);
        kh_value(ref_map, k) = ref;
    }

    
    kseq_destroy(seq);
    gzclose(fp);

}

static void search_context_kmp(const char* pat, const char* txt, uint8_t* result) {
    int M = strlen(pat);
    int N = strlen(txt);
    int* lps = (int*)malloc(M * sizeof(int));
    int len = 0;
    lps[0] = 0;
    int i = 1;
    while (i < M) {
        if (pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        }
        else {
            if (len != 0) {
                len = lps[len - 1];
            }
            else {
                lps[i] = 0;
                i++;
            }
        }
    }

    int count = 0;
    i = 0;
    int j = 0; 
  
    while ((N - i) >= (M - j)) {
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }
        if (j == M) {
            for (int k = 0; k < 1; k++) {
                result[i - j + k] = 1;
            }
            (count)++;
            j = lps[j - 1];
        }
        else if (i < N && pat[j] != txt[i]) {
            if (j != 0) j = lps[j - 1];
            else i = i + 1;
        }
    }
    free(lps);
}

int has_chr(const char * chr) {
    khiter_t k = kh_get(refm, ref_map, chr);
    return k != kh_end(ref_map);
}

ref_t * get_ref(const char * chr) {
    khiter_t k = kh_get(refm, ref_map, chr);
    if (k == kh_end(ref_map)) {
        return NULL;
    }
    return kh_value(ref_map, k);
}

void load_ref_contexts(int n_mod_codes, char ** mod_contexts) {
    for (khiter_t k = kh_begin(ref_map); k != kh_end(ref_map); ++k) {
        if (kh_exist(ref_map, k)) {
            ref_t * ref = kh_value(ref_map, k);
            ref->is_context = (uint8_t **) malloc(n_mod_codes * sizeof(uint8_t *));
            MALLOC_CHK(ref->is_context);
            for (int i = 0; i < n_mod_codes; i++) {
                ref->is_context[i] = (uint8_t *) calloc(ref->ref_seq_length, sizeof(uint8_t));
                MALLOC_CHK(ref->is_context[i]);

                if(strcmp(mod_contexts[i], "*") == 0){ // if context is *, fill is_context with 1
                    for(int j=0;j<ref->ref_seq_length;j++){
                        ref->is_context[i][j] = 1;
                    }
                } else {
                    search_context_kmp(mod_contexts[i], ref->forward, ref->is_context[i]);
                }
            }
        }
    }
}

void destroy_ref_forward() {
    khiter_t k;
    for (k = kh_begin(ref_map); k != kh_end(ref_map); ++k) {
        if (kh_exist(ref_map, k)) {
            ref_t * ref = kh_value(ref_map, k);
            free(ref->forward);
        }
    }
}

void destroy_ref(int n_mod_codes) {
    khiter_t k;
    for (k = kh_begin(ref_map); k != kh_end(ref_map); ++k) {
        if (kh_exist(ref_map, k)) {
            ref_t * ref = kh_value(ref_map, k);
            for (int i = 0; i < n_mod_codes; i++) {
                free(ref->is_context[i]);
            }
            char * ref_name = (char *) kh_key(ref_map, k);
            free(ref_name);
            free(ref->is_context);
            free(ref);
        }
    }
    kh_destroy(refm, ref_map);
}
