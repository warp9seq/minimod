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
KSEQ_INIT(gzFile, gzread)

ref_t * load_ref(const char * genome) {
    gzFile fp;
    kseq_t *seq;
    int l;
    
    fp = gzopen(genome, "r");
    F_CHK(fp, genome);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    // initialize ref
    ref_t *ref = (ref_t *) malloc(sizeof(ref_t));
    int cap = 1;
    ref->ref_lengths = (int32_t *) malloc(cap * sizeof(int32_t));
    ref->ref_names = (char **) malloc(cap * sizeof(char *));
    ref->ref_seq_lengths = (int32_t *) malloc(cap * sizeof(int32_t));
    ref->forward = (char **) malloc(cap * sizeof(char *));
    ref->num_ref = 0;
    
    int i=0;
    while ((l = kseq_read(seq)) >= 0) {
        ASSERT_MSG(l == (int) seq->seq.l, "Sequence length mismatch: %d vs %d", l, (int) seq->seq.l);
        if (i == cap) {
            cap <<= 1;
            ref->ref_names = (char **) realloc(ref->ref_names, cap * sizeof(char *));
            ref->ref_lengths = (int32_t *) realloc(ref->ref_lengths, cap * sizeof(int32_t));
            ref->ref_seq_lengths = (int32_t *) realloc(ref->ref_seq_lengths, cap * sizeof(int32_t));
            ref->forward = (char **) realloc(ref->forward, cap * sizeof(char *));
        }
        ref->ref_lengths[i] = seq->seq.l;
        ref->ref_names[i] = (char *) malloc(strlen(seq->name.s) + 1);
        MALLOC_CHK(ref->ref_names[i]);
        strcpy(ref->ref_names[i], seq->name.s);
        ref->ref_seq_lengths[i] = seq->seq.l;
        ref->forward[i] = (char *) malloc(seq->seq.l + 1);
        MALLOC_CHK(ref->forward[i]);
        strcpy(ref->forward[i], seq->seq.s);

        i++;
    }
    ref->num_ref = i;    
    
    kseq_destroy(seq);
    gzclose(fp);

    return ref;
}

void destroy_ref(ref_t * ref) {
    for(int i = 0; i < ref->num_ref; i++) {
        free(ref->ref_names[i]);
        free(ref->forward[i]);
    }
    free(ref->ref_names);
    free(ref->ref_seq_lengths);
    free(ref->forward);
    free(ref);
}
