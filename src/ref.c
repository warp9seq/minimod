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

int has_chr(const char * chr) {
    khiter_t k = kh_get(refm, ref_map, chr);
    return k != kh_end(ref_map);
}

ref_t * get_ref(const char * chr) {
    khiter_t k = kh_get(refm, ref_map, chr);
    return kh_value(ref_map, k);
}

void destroy_ref() {
    khiter_t k;
    for (k = kh_begin(ref_map); k != kh_end(ref_map); ++k) {
        if (kh_exist(ref_map, k)) {
            ref_t * ref = kh_value(ref_map, k);
            free(ref->forward);
            free(ref);
        }
    }
    kh_destroy(refm, ref_map);
}
