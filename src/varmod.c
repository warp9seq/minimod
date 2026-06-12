/**
 * @file varmod.c
 * @brief variant-aware modification tags
 */

#include "varmod.h"
#include "misc.h"
#include "error.h"
#include "khash.h"
#include "ref.h"
#include "ksort.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>


// Fallback for systems/compilers that don't expose drand48
#ifndef drand48
#define drand48() ((double)rand() / RAND_MAX)
#endif

extern const char *get_mm_tag_ptr(bam1_t *record);
extern uint8_t *get_ml_tag(bam1_t *record, uint32_t *len_ptr);
extern uint8_t get_hp_tag(bam1_t *record);

#define WILDCARD_STR "*"
#define THRESH_UINT8_TO_DBL(x) ((double)( (x + 0.5) / 256.0 ))
#define IS_ALPHA(c) (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')
#define IS_DIGIT(c) (c >= '0' && c <= '9')
// sentinel for offset=-1 (CpG C is one base upstream of the variant); stored as UINT16_MAX in uint16_t
#define OFFSET_NEG1 UINT16_MAX
// compute relative offset from variant position; returns OFFSET_NEG1 for the -1 case
#define REF_OFFSET(ref_pos, var_pos) (((ref_pos) < (var_pos)) ? OFFSET_NEG1 : (uint16_t)((ref_pos) - (var_pos)))
// convert stored offset to signed int for output
#define OFFSET_TO_INT(x) ((x) == OFFSET_NEG1 ? -1 : (int)(x))


// zero-allocation comparator
// key format: chrom\tpos\tstrand\tmod_code\toffset\thaplotype
static int cmp_key_fast(const char *key_a, const char *key_b) {
    // compare contig
    const char *tab_a = strchr(key_a, '\t');
    const char *tab_b = strchr(key_b, '\t');
    size_t len_a = tab_a ? (size_t)(tab_a - key_a) : strlen(key_a);
    size_t len_b = tab_b ? (size_t)(tab_b - key_b) : strlen(key_b);
    size_t min_len = len_a < len_b ? len_a : len_b;
    int cmp = strncmp(key_a, key_b, min_len);
    if (cmp != 0) return cmp;
    if (len_a != len_b) return (len_a < len_b) ? -1 : 1;
    if (!tab_a || !tab_b) return 0;

    // compare pos
    int pos_a = atoi(tab_a + 1);
    int pos_b = atoi(tab_b + 1);
    if (pos_a != pos_b) return (pos_a > pos_b) - (pos_a < pos_b);

    // skip to strand
    const char *p_a = strchr(tab_a + 1, '\t');
    const char *p_b = strchr(tab_b + 1, '\t');
    if (!p_a || !p_b) return 0;

    // compare strand (single char)
    if (p_a[1] != p_b[1]) return (p_a[1] > p_b[1]) - (p_a[1] < p_b[1]);

    // skip to mod_code
    const char *mc_a = strchr(p_a + 1, '\t');
    const char *mc_b = strchr(p_b + 1, '\t');
    if (!mc_a || !mc_b) return 0;

    // compare mod_code (string up to next tab)
    const char *mc_end_a = strchr(mc_a + 1, '\t');
    const char *mc_end_b = strchr(mc_b + 1, '\t');
    size_t mc_len_a = mc_end_a ? (size_t)(mc_end_a - mc_a - 1) : strlen(mc_a + 1);
    size_t mc_len_b = mc_end_b ? (size_t)(mc_end_b - mc_b - 1) : strlen(mc_b + 1);
    size_t mc_min = mc_len_a < mc_len_b ? mc_len_a : mc_len_b;
    cmp = strncmp(mc_a + 1, mc_b + 1, mc_min);
    if (cmp != 0) return cmp;
    if (mc_len_a != mc_len_b) return (mc_len_a < mc_len_b) ? -1 : 1;

    // compare offset (numeric, may be negative)
    if (!mc_end_a || !mc_end_b) return 0;
    int off_a = atoi(mc_end_a + 1);
    int off_b = atoi(mc_end_b + 1);
    return (off_a > off_b) - (off_a < off_b);
}

#define varfreq_kv_lt(a, b) (cmp_key_fast((a).key, (b).key) < 0)
#define varview_kv_lt(a, b) (cmp_key_fast((a).key, (b).key) < 0)

KSORT_INIT(varfreq, varfreq_kv_t, varfreq_kv_lt)
KSORT_INIT(varview, varview_kv_t, varview_kv_lt)

static int cmp_cg_entry(const void *a, const void *b) {
    const cg_entry_t *ea = a, *eb = b;
    return (ea->ref_cg_pos > eb->ref_cg_pos) - (ea->ref_cg_pos < eb->ref_cg_pos);
}

static inline int lower_bound_cg(cg_entry_t *entries, int len, int pos) {
    int lo = 0, hi = len;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (entries[mid].ref_cg_pos < pos) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static const int valid_bases[256] = { ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1, ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1, ['u'] = 1, ['n'] = 1 };
static const int valid_strands[256] = { ['+'] = 1, ['-'] = 1 };
static const int base_idx_lookup[256] = { ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['U'] = 3, ['N'] = 4, ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3, ['u'] = 3, ['n'] = 4 };
static const char base_complement_lookup[256] = { ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A', ['N'] = 'N', ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a', ['n'] = 'n' };
static const char* tested_cases[] = {"m[CG]", "h[CG]"};


static char* make_key(const char *chrom, int pos, uint16_t ins_offset, char * mod_code, char strand, int haplotype){
    int start_strlen = snprintf(NULL, 0, "%d", pos);
    int offset_strlen = snprintf(NULL, 0, "%d", ins_offset);
    int mod_code_strlen = strlen(mod_code);
    int haplotype_strlen = snprintf(NULL, 0, "%d", haplotype);
    int key_strlen = strlen(chrom) + start_strlen  + offset_strlen + mod_code_strlen + haplotype_strlen + 7;

    char* key = (char *)malloc(key_strlen * sizeof(char));
    MALLOC_CHK(key);
    snprintf(key, key_strlen, "%s\t%d\t%c\t%s\t%d\t%d", chrom, pos, strand, mod_code, OFFSET_TO_INT(ins_offset), haplotype);
    return key;
}

static void decode_key(char *key, char **chrom, int *pos, uint16_t * ins_offset, char **mod_code, char *strand, int *haplotype){
    char* token = strtok(key, "\t");
    *chrom = calloc(strlen(token)+1, sizeof(char));
    MALLOC_CHK(*chrom);
    strcpy(*chrom, token);

    *pos = atoi(strtok(NULL, "\t"));
    *strand = strtok(NULL, "\t")[0];
    
    token = strtok(NULL, "\t");
    *mod_code = calloc(strlen(token)+1, sizeof(char));
    MALLOC_CHK(*mod_code);
    strcpy(*mod_code, token);

    long v = strtol(strtok(NULL, "\t"), NULL, 10);
    *ins_offset = (v < 0) ? OFFSET_NEG1 : (uint16_t)v;
    *haplotype = atoi(strtok(NULL, "\t"));
}


int * find_cg_contexts_in_sequence(const char* seq, int seq_len, int* offsets_len) {
    int offsets_cap = 2;
    int * offsets = (int*)malloc(sizeof(int) * offsets_cap);
    MALLOC_CHK(offsets);
    int offsets_len_local = 0;
    for(int i = 0; i < seq_len - 1; i++) {
        if((seq[i] == 'C' || seq[i] == 'c') && (seq[i+1] == 'G' || seq[i+1] == 'g')) {
            if(offsets_len_local >= offsets_cap) {
                offsets_cap *= 2;
                offsets = (int*)realloc(offsets, sizeof(int) * offsets_cap);
                MALLOC_CHK(offsets);
            }
            offsets[offsets_len_local++] = i; // C
            offsets[offsets_len_local++] = i + 1; // G

        }
    }
    *offsets_len = offsets_len_local;
    return offsets;
}

int * find_cg_contexts_in_reverse_sequence(const char* seq, int seq_len, int* offsets_len) {
    int offsets_cap = 2;
    int * offsets = (int*)malloc(sizeof(int) * offsets_cap);
    MALLOC_CHK(offsets);
    int offsets_len_local = 0;
    for(int i = seq_len - 1; i > 0; i--) {
        if((seq[i] == 'G' || seq[i] == 'g') && (seq[i-1] == 'C' || seq[i-1] == 'c')) {
            if(offsets_len_local >= offsets_cap) {
                offsets_cap *= 2;
                offsets = (int*)realloc(offsets, sizeof(int) * offsets_cap);
                MALLOC_CHK(offsets);
            }
            offsets[offsets_len_local++] = seq_len - i - 1; // convert to reverse position
        }
    }
    *offsets_len = offsets_len_local;
    return offsets;
}

// check if the given CG positions already exists in the reference. we only omit CG formed by variants.
static inline int is_reference_cpg(const char *ref_seq, int ref_seq_length,
                                   int c_ref_pos, int g_ref_pos,
                                   int c_is_ins, int g_is_ins) {
    if (c_is_ins || g_is_ins) return 0;
    if (g_ref_pos != c_ref_pos + 1) return 0;
    if (c_ref_pos < 0 || c_ref_pos + 1 >= ref_seq_length) return 0;
    char c = ref_seq[c_ref_pos], g = ref_seq[c_ref_pos + 1];
    return ((c == 'C' || c == 'c') && (g == 'G' || g == 'g'));
}

void load_var_map(const char* vcf_file, const char* sample_name, khash_t(varm)* var_map, int haplotypes) {

    htsFile *vcf_fp = hts_open(vcf_file, "r");
    if(vcf_fp == NULL) {
        ERROR("Failed to open VCF file: %s\n", vcf_file);
        exit(EXIT_FAILURE);
    }

    bcf_hdr_t * vcf_hdr = bcf_hdr_read(vcf_fp);
    if(vcf_hdr == NULL) {
        ERROR("Failed to read VCF header: %s\n", vcf_file);
        exit(EXIT_FAILURE);
    }

    // resolve sample index: NULL -> sample 0 if any sample exists, else -1 (no GT info)
    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    int sample_idx = -1;
    if (n_samples > 0) {
        if (sample_name != NULL) {
            sample_idx = bcf_hdr_id2int(vcf_hdr, BCF_DT_SAMPLE, sample_name);
            if (sample_idx < 0) {
                ERROR("Sample '%s' not found in VCF %s", sample_name, vcf_file);
                exit(EXIT_FAILURE);
            }
        } else {
            sample_idx = 0;
        }
    }

    bcf1_t *rec = bcf_init();

    int32_t *gt_arr = NULL;
    int n_gt_arr = 0;

    while(bcf_read(vcf_fp, vcf_hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR | BCF_UN_FMT);

        const char *contig = bcf_hdr_id2name(vcf_hdr, rec->rid);
        int32_t pos = rec->pos;
        char * ref_allele = rec->d.allele[0];
        if(ref_allele[0] == '.') {
            ERROR("Invalid VCF record at %s:%d: REF allele is '.'\n", contig, pos + 1);
            exit(EXIT_FAILURE);
        }
        int ref_len = strlen(ref_allele);

        // assign haplotype of ALT alleles from the chosen sample's GT
        // -1 : ALT not present in this sample, skip
        // 0 : ALT is present, any haplotype(unphased, hom-alt, or missing GT)
        // 1+ : phased and ALT is present only on that haplotype slot
        int8_t alt_hap[256];
        for (int i = 0; i < rec->n_allele && i < 256; i++) alt_hap[i] = -1;

        // GT string of this sample
        char gt_str[64];
        gt_str[0] = '.'; gt_str[1] = '\0';

        int n_gt = (sample_idx >= 0) ? bcf_get_format_int32(vcf_hdr, rec, "GT", &gt_arr, &n_gt_arr) : -1;

        if (n_gt > 0 && sample_idx >= 0) {
            int ploidy = n_gt / n_samples;
            int32_t *gt = gt_arr + sample_idx * ploidy;

            // count how many slots carry each ALT, whether all separators are '|', and build the GT string
            int slot_carrying[256] = {0};
            int total_carriers = 0;
            int all_phased = 1; // first slot has no separator; check slots 1..ploidy-1
            int any_missing = 0;
            int gt_p = 0;
            for (int j = 0; j < ploidy; j++) {
                if (gt[j] == bcf_int32_vector_end) { ploidy = j; break; }
                if (j > 0 && gt_p < (int)sizeof(gt_str) - 4) {
                    gt_str[gt_p++] = bcf_gt_is_phased(gt[j]) ? '|' : '/';
                }
                if (bcf_gt_is_missing(gt[j])) {
                    if (gt_p < (int)sizeof(gt_str) - 1) gt_str[gt_p++] = '.';
                    any_missing = 1;
                    continue;
                }
                if (j >= 1 && !bcf_gt_is_phased(gt[j])) all_phased = 0;
                int a = bcf_gt_allele(gt[j]);
                gt_p += snprintf(gt_str + gt_p, sizeof(gt_str) - gt_p, "%d", a);
                if (a > 0 && a < 256) {
                    slot_carrying[a]++;
                    total_carriers++;
                }
            }
            gt_str[gt_p] = '\0';

            // any-missing GT means we can't trust phase even if separators are '|'
            if (any_missing) all_phased = 0;

            for (int i = 1; i < rec->n_allele && i < 256; i++) {
                if (slot_carrying[i] == 0) {
                    alt_hap[i] = -1; // ALT not present in this sample
                } else if (slot_carrying[i] > 1 || !all_phased || ploidy < 2) {
                    alt_hap[i] = 0;  // hom-alt, unphased, or haploid -> any read
                } else {
                    // exactly one slot carries this ALT in a phased call
                    for (int j = 0; j < ploidy; j++) {
                        if (gt[j] == bcf_int32_vector_end || bcf_gt_is_missing(gt[j])) continue;
                        if (bcf_gt_allele(gt[j]) == i) { alt_hap[i] = (int8_t)(j + 1); break; }
                    }
                }
            }

            // if total carriers is 0 the record is hom-ref for this sample; skip all ALTs
            if (total_carriers == 0) continue;
        } else {
            // no GT info available, all ALTs present, any haplotype
            for (int i = 1; i < rec->n_allele && i < 256; i++) alt_hap[i] = 0;
        }

        ref_t *ref = get_ref(contig);
        ASSERT_MSG(ref != NULL, "Contig %s not found in reference provided\n", contig);
        const char *ref_seq = ref->forward;
    
        char * contig_dup = (malloc(sizeof(char) * (strlen(contig) + 1)));
        MALLOC_CHK(contig_dup);
        strcpy(contig_dup, contig);
        
        int ret;
        khint_t k = kh_put(varm, var_map, contig_dup, &ret);
        vars_t * vars;

        // already exists
        if(ret == 0) {
            free(contig_dup);
            vars = kh_value(var_map, k);
        } else {
            vars = (vars_t *)malloc(sizeof(vars_t));
            MALLOC_CHK(vars);
            vars->vars_len = 0;
            vars->vars_cap = 1;
            vars->vars = (var_t*)malloc(sizeof(var_t) * vars->vars_cap);
            MALLOC_CHK(vars->vars);
            vars->cg_entries_len = 0;
            vars->cg_entries_cap = 4;
            vars->cg_entries = (cg_entry_t*)malloc(sizeof(cg_entry_t) * vars->cg_entries_cap);
            MALLOC_CHK(vars->cg_entries);
            kh_value(var_map, k) = vars;
        }
        
        for(int i = 1; i < rec->n_allele; i++) {
            char * alt_allele = rec->d.allele[i];
            if(alt_allele[0] == '.') continue; // reference call, no variant
            if(i >= 256 || alt_hap[i] < 0) continue; // ALT not present in chosen sample

            int alt_len = strlen(alt_allele);

            char * before_site = (char*)malloc(sizeof(char) * (ref_len + 3));
            MALLOC_CHK(before_site);

            char prev_base = (pos > 0) ? ref_seq[pos-1] : 'N';
            snprintf(before_site, ref_len + 3, "%c%s%c", prev_base, ref_allele, ref_seq[pos + ref_len]);

            // create site string = base before REF + ALT + base after REF
            char *after_site = (char*)malloc(sizeof(char) * (alt_len + 3));
            MALLOC_CHK(after_site);

            snprintf(after_site, alt_len + 3, "%c%s%c", prev_base, alt_allele, ref->forward[pos + ref_len]);

            int n_cg_offsets = 0;
            int * cg_offsets = find_cg_contexts_in_sequence(after_site, strlen(after_site), &n_cg_offsets);

            if(vars->vars_len >= vars->vars_cap) {
                vars->vars_cap *= 2;
                vars->vars = (var_t*)realloc(vars->vars, sizeof(var_t) * vars->vars_cap);
                MALLOC_CHK(vars->vars);
            }
            var_t var;
            var.pos = pos;
            var.ref_len = ref_len;
            var.cg_offsets_len = n_cg_offsets;
            var.cg_offsets = cg_offsets;
            var.ref_allele = (char*)malloc(sizeof(char) * (ref_len + 1));
            MALLOC_CHK(var.ref_allele);
            strcpy(var.ref_allele, ref_allele);
            var.alt_allele = (char*)malloc(sizeof(char) * (alt_len + 1));
            MALLOC_CHK(var.alt_allele);
            strcpy(var.alt_allele, alt_allele);
            var.before_site = before_site;
            var.after_site = after_site;
            var.hap = alt_hap[i];
            var.gt = (char*)malloc(strlen(gt_str) + 1);
            MALLOC_CHK(var.gt);
            strcpy(var.gt, gt_str);
            vars->vars[vars->vars_len++] = var;

            // build CG-position index for O(log N) lookup during processing
            int var_idx = vars->vars_len - 1;
            for (int o = 0; o < n_cg_offsets; o++) {
                int o_val = cg_offsets[o];
                int8_t is_ins = (o_val > ref_len && o_val <= alt_len) ? 1 : 0;

                // skip if CG is a ref CG
                {
                    int c_idx = (after_site[o_val] == 'C' || after_site[o_val] == 'c') ? o_val : o_val - 1;
                    int g_idx = c_idx + 1;
                    int c_rp = (c_idx > alt_len) ? pos + ref_len + (c_idx - alt_len - 1) : pos - 1 + c_idx;
                    int g_rp = (g_idx > alt_len) ? pos + ref_len + (g_idx - alt_len - 1) : pos - 1 + g_idx;
                    int c_in = (c_idx > ref_len && c_idx <= alt_len);
                    int g_in = (g_idx > ref_len && g_idx <= alt_len);
                    if (is_reference_cpg(ref_seq, ref->ref_seq_length, c_rp, g_rp, c_in, g_in)) continue;
                }

                if (vars->cg_entries_len >= vars->cg_entries_cap) {
                    vars->cg_entries_cap *= 2;
                    vars->cg_entries = (cg_entry_t*)realloc(vars->cg_entries, sizeof(cg_entry_t) * vars->cg_entries_cap);
                    MALLOC_CHK(vars->cg_entries);
                }

                int ref_cg_pos = (o_val > alt_len)
                    ? pos + ref_len + (o_val - alt_len - 1)
                    : pos - 1 + o_val;
                vars->cg_entries[vars->cg_entries_len].ref_cg_pos = ref_cg_pos;
                vars->cg_entries[vars->cg_entries_len].var_idx = var_idx;
                vars->cg_entries[vars->cg_entries_len].is_insertion_only = is_ins;
                vars->cg_entries[vars->cg_entries_len].is_compound = 0;
                vars->cg_entries_len++;
            }

            // TODO: Refer vcf spec. and handle indels at terminal positions of the contig.

        }

    }
    
    free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(vcf_hdr);
    hts_close(vcf_fp);


    if (haplotypes) {
    for (khint_t k = kh_begin(var_map); k != kh_end(var_map); ++k) {
        if (!kh_exist(var_map, k)) continue;
        const char *contig = (const char *) kh_key(var_map, k);
        ref_t *ref = get_ref(contig);
        if (ref == NULL) continue;
        const char *ref_seq = ref->forward;
        vars_t *vars = kh_value(var_map, k);

        int8_t max_hap = 0;
        for (int i = 0; i < vars->vars_len; i++) {
            if (vars->vars[i].hap > max_hap) max_hap = vars->vars[i].hap;
        }
        if (max_hap < 1) continue;

        int idx_cap = 16;
        int *hap_idx = (int*)malloc(sizeof(int) * idx_cap);
        MALLOC_CHK(hap_idx);

        for (int hap = 1; hap <= max_hap; hap++) {
            int n_hap = 0;
            for (int i = 0; i < vars->vars_len; i++) {
                if (vars->vars[i].hap != hap) continue;
                if (n_hap >= idx_cap) { idx_cap *= 2; hap_idx = (int*)realloc(hap_idx, sizeof(int)*idx_cap); MALLOC_CHK(hap_idx); }
                hap_idx[n_hap++] = i;
            }
            if (n_hap < 2) continue;

            // insertion sort by pos
            for (int i = 1; i < n_hap; i++) {
                int key = hap_idx[i];
                int j = i - 1;
                while (j >= 0 && vars->vars[hap_idx[j]].pos > vars->vars[key].pos) {
                    hap_idx[j+1] = hap_idx[j]; j--;
                }
                hap_idx[j+1] = key;
            }

            // drop variants that overlap their predecessor (rare in normal VCFs); compact in place
            int w = 0;
            for (int i = 0; i < n_hap; i++) {
                if (w > 0 && vars->vars[hap_idx[i]].pos < vars->vars[hap_idx[w-1]].pos + vars->vars[hap_idx[w-1]].ref_len) {
                    continue; // overlapping; skip
                }
                hap_idx[w++] = hap_idx[i];
            }
            n_hap = w;
            if (n_hap < 2) continue;

            // span [span_start, span_end) covers 1 ref base before the first variant and 1 after the last
            int span_start = vars->vars[hap_idx[0]].pos - 1;
            int last_end = vars->vars[hap_idx[n_hap-1]].pos + vars->vars[hap_idx[n_hap-1]].ref_len;
            int span_end = last_end + 1;
            if (span_start < 0) span_start = 0;
            if (span_end > ref->ref_seq_length) span_end = ref->ref_seq_length;

            // worst-case patched length = ref span + total inserted bases
            long total_ins = 0;
            for (int v = 0; v < n_hap; v++) {
                int al = (int)strlen(vars->vars[hap_idx[v]].alt_allele);
                int rl = vars->vars[hap_idx[v]].ref_len;
                if (al > rl) total_ins += (al - rl);
            }
            long cap = (long)(span_end - span_start) + total_ins + 8;

            char    *hap_seq     = (char*)   malloc(sizeof(char)    * cap);
            int     *hap_ref_pos = (int*)    malloc(sizeof(int)     * cap);
            int8_t  *hap_is_ins  = (int8_t*) malloc(sizeof(int8_t)  * cap);
            int     *hap_owner   = (int*)    malloc(sizeof(int)     * cap); // var_idx, or -1 for ref base
            MALLOC_CHK(hap_seq); MALLOC_CHK(hap_ref_pos); MALLOC_CHK(hap_is_ins); MALLOC_CHK(hap_owner);

            long len = 0;
            int rp = span_start;
            int vptr = 0;
            while (rp < span_end) {
                if (vptr < n_hap && rp == vars->vars[hap_idx[vptr]].pos) {
                    int vi = hap_idx[vptr];
                    var_t *vt = &vars->vars[vi];
                    int alt_len = (int)strlen(vt->alt_allele);
                    for (int a = 0; a < alt_len; a++) {
                        hap_seq[len] = vt->alt_allele[a];
                        // ref_cg_pos convention matches the per-variant after_site path:
                        // alt index a uses vt->pos + a, is_ins=1 when a >= ref_len.
                        hap_ref_pos[len] = vt->pos + a;
                        hap_is_ins[len] = (a >= vt->ref_len) ? 1 : 0;
                        hap_owner[len] = vi;
                        len++;
                    }
                    rp = vt->pos + vt->ref_len; // skip the replaced ref bases
                    vptr++;
                } else {
                    hap_seq[len] = ref_seq[rp];
                    hap_ref_pos[len] = rp;
                    hap_is_ins[len] = 0;
                    hap_owner[len] = -1;
                    len++;
                    rp++;
                }
            }

            for (long p = 0; p + 1 < len; p++) {
                char cc = hap_seq[p], gc = hap_seq[p+1];
                if ((cc != 'C' && cc != 'c') || (gc != 'G' && gc != 'g')) continue;
                // skip CpGs already present in the reference; only emit variant-introduced ones
                if (is_reference_cpg(ref_seq, ref->ref_seq_length,
                                     hap_ref_pos[p], hap_ref_pos[p+1],
                                     hap_is_ins[p], hap_is_ins[p+1])) continue;
                int ref_cg_pos = hap_ref_pos[p];
                int8_t is_ins = hap_is_ins[p];
                // attribute to the owning variant of the C; if C is a ref base, fall back to the
                // owner of the G, else the nearest preceding variant base.
                int attr = hap_owner[p];
                if (attr < 0) attr = hap_owner[p+1];
                if (attr < 0) { for (long q = p; q >= 0; q--) { if (hap_owner[q] >= 0) { attr = hap_owner[q]; break; } } }
                if (attr < 0) attr = hap_idx[0];

                // dedup against existing entries at same (ref_cg_pos, is_ins) on this hap;
                // for insertion entries also require matching var.pos (matches scan filter)
                int dup = 0;
                for (int e = 0; e < vars->cg_entries_len; e++) {
                    if (vars->cg_entries[e].ref_cg_pos != ref_cg_pos) continue;
                    if (vars->cg_entries[e].is_insertion_only != is_ins) continue;
                    if (is_ins && vars->vars[vars->cg_entries[e].var_idx].pos != vars->vars[attr].pos) continue;
                    if (vars->vars[vars->cg_entries[e].var_idx].hap == hap) { dup = 1; break; }
                }
                if (dup) continue;

                if (vars->cg_entries_len >= vars->cg_entries_cap) {
                    vars->cg_entries_cap *= 2;
                    vars->cg_entries = (cg_entry_t*)realloc(vars->cg_entries, sizeof(cg_entry_t) * vars->cg_entries_cap);
                    MALLOC_CHK(vars->cg_entries);
                }
                vars->cg_entries[vars->cg_entries_len].ref_cg_pos = ref_cg_pos;
                vars->cg_entries[vars->cg_entries_len].var_idx = attr;
                vars->cg_entries[vars->cg_entries_len].is_insertion_only = is_ins;
                vars->cg_entries[vars->cg_entries_len].is_compound = 1;
                vars->cg_entries_len++;
            }

            free(hap_seq); free(hap_ref_pos); free(hap_is_ins); free(hap_owner);
        }
        free(hap_idx);
    }
    } // if (haplotypes)

    // sort each contig's CG-position index for binary search during processing
    for (khint_t k = kh_begin(var_map); k != kh_end(var_map); ++k) {
        if (kh_exist(var_map, k)) {
            vars_t *vars = kh_value(var_map, k);
            qsort(vars->cg_entries, vars->cg_entries_len, sizeof(cg_entry_t), cmp_cg_entry);
        }
    }
}

void destroy_var_map(khash_t(varm)* var_map) {
    khint_t k;
    for (k = kh_begin(var_map); k != kh_end(var_map); ++k) {
        if (kh_exist(var_map, k)) {
            char * key = (char*) kh_key(var_map, k);
            vars_t *vars = kh_value(var_map, k);
            free(key);
            for (int i = 0; i < vars->vars_len; i++) {
                free(vars->vars[i].before_site);
                free(vars->vars[i].after_site);
                free(vars->vars[i].ref_allele);
                free(vars->vars[i].alt_allele);
                free(vars->vars[i].gt);
                free(vars->vars[i].cg_offsets);
            }
            free(vars->vars);
            free(vars->cg_entries);
            free(vars);
        }
    }
    kh_destroy(varm, var_map);
}


static void get_aln(core_t * core, db_t *db, bam_hdr_t *hdr, bam1_t *record, int bam_i){
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = (tid >= 0) ? hdr->target_name[tid] : "*";
    int32_t pos = record->core.pos;
    int32_t end = bam_endpos(record);

    const char *qname = bam_get_qname(record);

    int8_t rev = bam_is_rev(record);

    uint32_t *cigar = bam_get_cigar(record);
    uint32_t n_cigar = record->core.n_cigar;

    int seq_len = record->core.l_qseq;

    ref_t *ref = get_ref(tname);
        ASSERT_MSG(ref != NULL, "Contig %s not found in reference provided\n", tname);
  
    int read_pos = 0;
    int ref_pos = pos;

    int * aligned_pairs = db->aln[bam_i];
    //fill the aligned_pairs array with -1
    for(int i=0;i<seq_len;i++){
        aligned_pairs[i] = -1;
    }

    for(int i=0;i<seq_len;i++){
        db->ins[bam_i][i] = -1;
        db->ins_offset[bam_i][i] = 0;
    }

    for (uint32_t ci = 0; ci < n_cigar; ++ci) {
        uint32_t c = cigar[ci];
        if(rev) {
            c = cigar[n_cigar - ci - 1];
        }
        int cigar_len = bam_cigar_oplen(c);
        int cigar_op = bam_cigar_op(c);

        // Set the amount that the ref/read positions should be incremented
        // based on the cigar operation
        int read_inc = 0;
        int ref_inc = 0;

        // Process match between the read and the reference
        int8_t is_aligned = 0, is_inserted = 0;
        if(cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            is_aligned = 1;
            read_inc = 1;
            ref_inc = 1;
        } else if(cigar_op == BAM_CDEL) {
            ref_inc = 1;
        } else if(cigar_op == BAM_CREF_SKIP) {
            ref_inc = 1;
        } else if(cigar_op == BAM_CINS) {
            read_inc = 1;
            is_inserted = 1;
        } else if(cigar_op == BAM_CSOFT_CLIP) {
            read_inc = 1;
        } else if(cigar_op == BAM_CHARD_CLIP) { // TODO: use MN tag (seq len at the time MM value was last written) to check this?
            read_inc = 0;
            ERROR("Hard clipping found in %s and they are not supported.\nTry following workarounds.\n\t01. Filter out non-primary alignments\n\t\tsamtools view -h -F 2308 reads.bam -o primary_reads.bam\n\t02. Use minimap2 with -Y to use soft clipping for suplimentary alignments.\n", qname); 
            exit(EXIT_FAILURE);
        } else {
            ERROR("Unhandled CIGAR OPT Cigar: %d\n", cigar_op);
            exit(EXIT_FAILURE);
        }

        // Iterate over the pairs of aligned bases
        for(int j = 0; j < cigar_len; ++j) {
            if(is_aligned) {
                ASSERT_MSG(read_pos < seq_len, "read_pos:%d seq_len:%d\n", read_pos, seq_len);
                int start = ref_pos;
                if(rev) {
                    start = pos + end - ref_pos - 1;
                }
                aligned_pairs[read_pos] = start;

                ASSERT_MSG(ref_pos >= 0 && ref_pos < ref->ref_seq_length, "ref_pos:%d ref_len:%d\n", ref_pos, ref->ref_seq_length);
                ASSERT_MSG(ref->ref_seq_length == hdr->target_len[tid], "ref_len:%d target_len:%d\n", ref->ref_seq_length, hdr->target_len[tid]);
            }

            if(is_inserted) {
                ASSERT_MSG(read_pos < seq_len, "read_pos:%d seq_len:%d\n", read_pos, seq_len);
                int start = ref_pos-1;
                int offset = j+1;
                if(rev) {
                    start = pos + end - ref_pos - 1;
                    offset = cigar_len - j;
                }
                db->ins[bam_i][read_pos] = start;
                db->ins_offset[bam_i][read_pos] = offset;
            }

            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }
}


void update_varfreq_map(khash_t(varfreqm) *varfreq_map, const char *tname, int ref_pos, int ins_offset, char *mod_code, char strand, int haplotype, int is_called, int is_mod, const char *ref_allele, const char *alt_allele, const char *gt, int var_pos) {
    char * key = make_key(tname, ref_pos, ins_offset, mod_code, strand, haplotype);
    khiter_t k = kh_get(varfreqm, varfreq_map, key);
    if (k == kh_end(varfreq_map)) { // not found, add
        varfreq_t * varfreq = (varfreq_t *)malloc(sizeof(varfreq_t));
        MALLOC_CHK(varfreq);
        varfreq->n_called = is_called;
        varfreq->n_mod = is_mod;
        varfreq->ref_allele = ref_allele;
        varfreq->alt_allele = alt_allele;
        varfreq->gt = gt;
        varfreq->var_pos = var_pos;
        int ret;
        k = kh_put(varfreqm, varfreq_map, key, &ret);
        kh_value(varfreq_map, k) = varfreq;
    } else { // found, update
        varfreq_t * varfreq = kh_value(varfreq_map, k);
        varfreq->n_called += is_called;
        varfreq->n_mod += is_mod;
        // check if varfreq->n_called overflows
        if(varfreq->n_called == 0){
            ERROR("n_called overflowed for key %s. Please report this issue.", key);
            exit(EXIT_FAILURE);
        }
        free(key);
    }
}

void add_varview_entry(khash_t(varviewm) *varview_map, const char *tname, int ref_pos, int ins_offset, char *mod_code, char strand, int haplotype, uint8_t mod_prob, int read_pos, var_t var) {

    char *key = make_key(tname, ref_pos, ins_offset, mod_code, strand, haplotype);
    khiter_t k = kh_get(varviewm, varview_map, key);
    if (k == kh_end(varview_map)) { // not found, add
        varview_t *varview = (varview_t *)malloc(sizeof(varview_t));
        MALLOC_CHK(varview);
        varview->mod_prob = mod_prob;
        varview->read_pos = read_pos;
        varview->var = var;
        int ret;
        k = kh_put(varviewm, varview_map, key, &ret);
        kh_value(varview_map, k) = varview;
    } else { // found, update
        free(key);
    }
}

void varviewfreq_single(core_t * core, db_t *db, int32_t bam_i) {
    bam1_t *record = db->bam_recs[bam_i];
    int8_t rev = bam_is_rev(record);
    bam_hdr_t *hdr = core->bam_hdr;
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = (tid >= 0) ? hdr->target_name[tid] : "*";
    uint8_t *seq = bam_get_seq(record);
    uint32_t seq_len = record->core.l_qseq;
    char strand = rev ? '-' : '+';
    const char *mm_string = db->mm[bam_i];
    uint8_t *ml = db->ml[bam_i];
    // phase-aware processing only when --haplotypes is specified with phased VCF and BAM
    uint8_t read_hp = core->opt.haplotypes ? get_hp_tag(record) : 0;
    int haplotype = core->opt.haplotypes ? (int)read_hp : -1;
    int *aln_pairs = db->aln[bam_i];

    ref_t *ref = get_ref(tname);
    ASSERT_MSG(ref != NULL, "Contig %s not found in reference provided\n", tname);

    get_aln(core, db, hdr, record, bam_i);

    int **bases_pos = db->bases_pos[bam_i];
    int bases_pos_lens[N_BASES] = {0};
    memset(db->mod_codes[bam_i], 0, core->opt.n_mods);

    int i;
    for (i = 0; i < (int)seq_len; i++) {
        int base_char = seq_nt16_str[bam_seqi(seq, i)];
        int idx = base_idx_lookup[(int)base_char];
        bases_pos[idx][bases_pos_lens[idx]++] = i;
    }

    khint_t contig_k = kh_get(varm, core->var_map, tname);
    vars_t *vars = (contig_k != kh_end(core->var_map)) ? kh_value(core->var_map, contig_k) : NULL;
    khint_t wc_k = kh_get(modcodesm, core->opt.modcodes_map, WILDCARD_STR);
    int has_wildcard = (wc_k != kh_end(core->opt.modcodes_map));

    int mm_str_len = strlen(mm_string);
    i = 0;
    int ml_start_idx = 0;

    char modbase;
    // char mod_strand;
    char * mod_codes = db->mod_codes[bam_i];
    int mod_codes_len;
    int * skip_counts = db->skip_counts[bam_i];
    int skip_counts_len;
    char status_flag;

    while (i < mm_str_len) {
        // reset skip counts and mod codes
        skip_counts_len = 0;
        mod_codes_len = 0;

        // set default status flag to '.' (when not present or '.' in the MM string)
        status_flag = '.';

        // get base
        if(i < mm_str_len) {
            ASSERT_MSG(valid_bases[(int)mm_string[i]], "Invalid base:%c\n", mm_string[i]);
            modbase = mm_string[i] == 'U' ? 'T' : mm_string[i]; // convert U to T
            i++;
        }

        // get strand
        if(i < mm_str_len) {
            ASSERT_MSG(valid_strands[(int)mm_string[i]], "Invalid strand:%c\n", mm_string[i]);
            // mod_strand = mm_string[i];
            i++;
        }

        // get base modification codes. can handle multiple codes giver as chars. TO-DO: handle when given as a ChEBI id
        int j = 0;
        int has_nums = 0;
        int has_alpha = 0;
        while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';' && mm_string[i] != '?' && mm_string[i] != '.') {

            // ASSERT_MSG(valid_mod_codes[(int)mm_string[i]], "Invalid base modification code:%c\n", mm_string[i]);

            if(IS_DIGIT(mm_string[i])) {
                has_nums = 1;
            } else if(IS_ALPHA(mm_string[i])) {
                has_alpha = 1;
            } else {
                ERROR("Invalid base modification code:%c. Modification codes should be either numeric or alphabetic.\n", mm_string[i]);
                exit(EXIT_FAILURE);
            }

            if(j >= db->mod_codes_cap[bam_i]) {
                db->mod_codes_cap[bam_i] *= 2;
                db->mod_codes[bam_i] = (char *)realloc(db->mod_codes[bam_i], sizeof(char) * (db->mod_codes_cap[bam_i] + 1)); // +1 for null terminator
                MALLOC_CHK(db->mod_codes[bam_i]);
            }
            mod_codes = db->mod_codes[bam_i];

            mod_codes[j] = mm_string[i];
            j++;
            i++;
        }
        mod_codes[j] = '\0';
        mod_codes_len = j;

        if(has_nums) {
            mod_codes_len = 1; // if chebi id is given, then only one code is present
        }

        // validate mod codes
        ASSERT_MSG(mod_codes_len>0, "Invalid modification codes:%s. Modification codes cannot be empty.\n", mod_codes);
        ASSERT_MSG((has_nums && has_alpha) == 0, "Invalid modification codes:%s. Modification codes should be either numeric or alphabetic, not both.\n", mod_codes);
        
        // get modification status flag
        if(i < mm_str_len && ( mm_string[i] == '?' || mm_string[i] == '.' )) {
            status_flag = mm_string[i];
            i++;
        } else { // if not present, set to '.'
            status_flag = '.';
        }

        // get skip counts
        int k = 0;
        while (i < mm_str_len && mm_string[i] != ';') {

            // skip if a comma
            if(i < mm_str_len && mm_string[i] == ',') {
                i++;
                continue;
            }

            char skip_count_str[10];
            int l = 0;
            while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';') {
                skip_count_str[l] = mm_string[i];
                i++;
                l++;
                assert(l < 10); // if this fails, use dynamic allocation for skip_count_str
            }
            skip_count_str[l] = '\0';
            ASSERT_MSG(l > 0, "Invalid skip count:%d.\n", l);
            skip_counts[k] = atoi(skip_count_str);
            ASSERT_MSG(skip_counts[k] >= 0, "Skip count cannot be negative: %d.\n", skip_counts[k]);
            
            k++;
        }
        skip_counts_len = k;
        i++;

        char mb = rev? base_complement_lookup[(int)modbase] : modbase;
        int idx = base_idx_lookup[(int)mb];

        modcodem_t *req_mods[MOD_CODE_LEN];
        int req_mod_valid[MOD_CODE_LEN];
        for (int m = 0; m < mod_codes_len; m++) {
            char *mc = has_nums ? mod_codes : &mod_codes[m];
            khint_t mk = has_wildcard ? wc_k : kh_get(modcodesm, core->opt.modcodes_map, mc);
            if (mk == kh_end(core->opt.modcodes_map)) {
                req_mod_valid[m] = 0;
                req_mods[m] = NULL;
            } else {
                req_mod_valid[m] = 1;
                req_mods[m] = kh_value(core->opt.modcodes_map, mk);
            }
        }

        int ml_idx = ml_start_idx;
        int base_rank = -1; // 0-based rank
        for(int c=0; c<skip_counts_len; c++) {
            base_rank += skip_counts[c] + 1;

            int read_pos;

            if (modbase == 'N') {
                read_pos = rev ? (int)(seq_len - base_rank - 1) : base_rank;
            } else {
                read_pos = rev ? bases_pos[idx][bases_pos_lens[idx] - base_rank - 1] : bases_pos[idx][base_rank];
            }

            ASSERT_MSG(read_pos >= 0 && read_pos < (int)seq_len, "Read pos cannot exceed seq len. read_pos: %d seq_len: %d\n", read_pos, seq_len);

            int fastq_read_pos = rev ? (int)(seq_len - read_pos - 1) : read_pos;
            int ref_pos = aln_pairs[fastq_read_pos];
            int ins_start = db->ins[bam_i][fastq_read_pos];
            int ins_offset = db->ins_offset[bam_i][fastq_read_pos];

            if (ref_pos == -1 && ins_start == -1) {
                if (mod_codes_len > 0) {
                    ml_idx = ml_start_idx + c * mod_codes_len + mod_codes_len - 1;
                }
                continue;
            }

            for (int m = 0; m < mod_codes_len; m++) {
                ml_idx = ml_start_idx + c * mod_codes_len + m;

                if (!req_mod_valid[m]) continue;

                char *mod_code = has_nums ? mod_codes : &mod_codes[m];
                modcodem_t *req_mod = req_mods[m];

                if (vars) {
                    int want_ins = (ref_pos == -1);
                    int lookup_pos = want_ins ? (ins_start + ins_offset) : ref_pos;
                    int out_pos = want_ins ? ins_start : ref_pos;

                    int ei = lower_bound_cg(vars->cg_entries, vars->cg_entries_len, lookup_pos);
                    for (; ei < vars->cg_entries_len && vars->cg_entries[ei].ref_cg_pos == lookup_pos; ei++) {
                        if (vars->cg_entries[ei].is_insertion_only != want_ins) continue;
                        if (vars->cg_entries[ei].is_compound && !core->opt.haplotypes) continue;
                        var_t var = vars->vars[vars->cg_entries[ei].var_idx];
                        if (want_ins && var.pos != ins_start) continue;
                        // phase-aware filter: only when --haplotypes is on; phased ALT (var.hap > 0) only counts reads with matching HP tag
                        if (core->opt.haplotypes && var.hap > 0 && (int)read_hp != var.hap) continue;
                        uint8_t mod_prob = ml[ml_idx];
                        uint16_t offset = want_ins ? (uint16_t)ins_offset : REF_OFFSET(out_pos, var.pos);
                        if (core->opt.subtool == VARVIEW) {
                            add_varview_entry(db->varview_maps[bam_i], tname, out_pos, offset, mod_code, strand, haplotype, mod_prob, fastq_read_pos, var);
                        } else {
                            double mod_prob_dbl = THRESH_UINT8_TO_DBL(mod_prob);
                            double thresh = req_mod->thresh;
                            int is_called = 0, is_mod = 0;
                            if (mod_prob_dbl >= thresh) { is_called = 1; is_mod = 1; }
                            else if (mod_prob_dbl <= 1 - thresh) { is_called = 1; }
                            else continue;
                            update_varfreq_map(db->varfreq_maps[bam_i], tname, out_pos, offset, mod_code, strand, haplotype, is_called, is_mod, var.ref_allele, var.alt_allele, var.gt, var.pos);
                        }
                    }
                }
            }
        }
        if (skip_counts_len > 0) ml_start_idx = ml_idx + 1;

        // Skipped bases (mod_prob = 0, status_flag == '.')
        if (status_flag == '.') {
            int skip_base_rank = -1;
            int prev_skip_base_rank = -1;
            for (int c = 0; c < skip_counts_len; c++) {
                skip_base_rank += skip_counts[c] + 1;

                for (int s = prev_skip_base_rank + 1; s < skip_base_rank; s++) {
                    int skip_read_pos;
                    if (modbase == 'N') {
                        skip_read_pos = rev ? (int)(seq_len - s - 1) : s;
                    } else {
                        skip_read_pos = rev ? bases_pos[idx][bases_pos_lens[idx] - s - 1] : bases_pos[idx][s];
                    }

                    ASSERT_MSG(skip_read_pos >= 0 && skip_read_pos < (int)seq_len, "Read pos cannot exceed seq len. read_pos: %d seq_len: %d\n", skip_read_pos, seq_len);

                    int skip_fastq_read_pos = rev ? (int)(seq_len - skip_read_pos - 1) : skip_read_pos;
                    int skip_ref_pos = aln_pairs[skip_fastq_read_pos];
                    int skip_ins_start = db->ins[bam_i][skip_fastq_read_pos];
                    int skip_ins_offset = db->ins_offset[bam_i][skip_fastq_read_pos];

                    if (skip_ref_pos == -1 && skip_ins_start == -1) continue;

                    for (int m = 0; m < mod_codes_len; m++) {
                        if (!req_mod_valid[m]) continue;
                        char *mod_code = has_nums ? mod_codes : &mod_codes[m];

                        if (vars) {
                            int want_ins = (skip_ref_pos == -1);
                            int lookup_pos = want_ins ? (skip_ins_start + skip_ins_offset) : skip_ref_pos;
                            int out_pos = want_ins ? skip_ins_start : skip_ref_pos;

                            int ei = lower_bound_cg(vars->cg_entries, vars->cg_entries_len, lookup_pos);
                            for (; ei < vars->cg_entries_len && vars->cg_entries[ei].ref_cg_pos == lookup_pos; ei++) {
                                if (vars->cg_entries[ei].is_insertion_only != want_ins) continue;
                                if (vars->cg_entries[ei].is_compound && !core->opt.haplotypes) continue;
                                var_t var = vars->vars[vars->cg_entries[ei].var_idx];
                                if (want_ins && var.pos != skip_ins_start) continue;
                                if (core->opt.haplotypes && var.hap > 0 && (int)read_hp != var.hap) continue;
                                uint16_t offset = want_ins ? (uint16_t)skip_ins_offset : REF_OFFSET(out_pos, var.pos);
                                if (core->opt.subtool == VARVIEW) {
                                    add_varview_entry(db->varview_maps[bam_i], tname, out_pos, offset, mod_code, strand, haplotype, 0, skip_fastq_read_pos, var);
                                } else {
                                    update_varfreq_map(db->varfreq_maps[bam_i], tname, out_pos, offset, mod_code, strand, haplotype, 1, 0, var.ref_allele, var.alt_allele, var.gt, var.pos);
                                }
                            }
                        }
                    }
                }
                prev_skip_base_rank = skip_base_rank;
            }

            // handle skipped bases after the last skip count
            for(int s=prev_skip_base_rank+1; s<bases_pos_lens[idx]; s++) {
                int skip_read_pos;
                if (modbase == 'N') {
                    skip_read_pos = rev ? (int)(seq_len - s - 1) : s;
                } else {
                    skip_read_pos = rev ? bases_pos[idx][bases_pos_lens[idx] - s - 1] : bases_pos[idx][s];
                }

                ASSERT_MSG(skip_read_pos >= 0 && skip_read_pos < (int)seq_len, "Read pos cannot exceed seq len. read_pos: %d seq_len: %d\n", skip_read_pos, seq_len);

                int skip_fastq_read_pos = rev ? (int)(seq_len - skip_read_pos - 1) : skip_read_pos;
                int skip_ref_pos = aln_pairs[skip_fastq_read_pos];
                int skip_ins_start = db->ins[bam_i][skip_fastq_read_pos];
                int skip_ins_offset = db->ins_offset[bam_i][skip_fastq_read_pos];

                if (skip_ref_pos == -1 && skip_ins_start == -1) continue;

                for (int m = 0; m < mod_codes_len; m++) {
                    if (!req_mod_valid[m]) continue;
                    char *mod_code = has_nums ? mod_codes : &mod_codes[m];

                    if (vars) {
                        int want_ins = (skip_ref_pos == -1);
                        int lookup_pos = want_ins ? (skip_ins_start + skip_ins_offset) : skip_ref_pos;
                        int out_pos = want_ins ? skip_ins_start : skip_ref_pos;

                        int ei = lower_bound_cg(vars->cg_entries, vars->cg_entries_len, lookup_pos);
                        for (; ei < vars->cg_entries_len && vars->cg_entries[ei].ref_cg_pos == lookup_pos; ei++) {
                            if (vars->cg_entries[ei].is_insertion_only != want_ins) continue;
                            if (vars->cg_entries[ei].is_compound && !core->opt.haplotypes) continue;
                            var_t var = vars->vars[vars->cg_entries[ei].var_idx];
                            if (want_ins && var.pos != skip_ins_start) continue;
                            uint16_t offset = want_ins ? (uint16_t)skip_ins_offset : REF_OFFSET(out_pos, var.pos);
                            if (core->opt.subtool == VARVIEW) {
                                add_varview_entry(db->varview_maps[bam_i], tname, out_pos, offset, mod_code, strand, haplotype, 0, skip_fastq_read_pos, var);
                            } else {
                                update_varfreq_map(db->varfreq_maps[bam_i], tname, out_pos, offset, mod_code, strand, haplotype, 1, 0, var.ref_allele, var.alt_allele, var.gt, var.pos);
                            }
                        }
                    }
                }
            }
        }
    }
}


void warn_untested_cases_var(opt_t * opt) {
    int n_tested_cases = sizeof(tested_cases) / sizeof(tested_cases[0]);
    khint_t i;
    for(i=kh_begin(opt->modcodes_map); i < kh_end(opt->modcodes_map); ++i) {
        if (!kh_exist(opt->modcodes_map, i)) continue;
        char * mod_code = (char *) kh_key(opt->modcodes_map, i);
        char * context = kh_value(opt->modcodes_map, i)->context;

        char * mod_code_with_context = (char *)malloc(strlen(mod_code) + strlen(context) + 3); // for null terminator and brackets
        MALLOC_CHK(mod_code_with_context);
        snprintf(mod_code_with_context, strlen(mod_code) + strlen(context) + 3, "%s[%s]", mod_code, context);

        int is_tested = false;
        for(int j=0; j < n_tested_cases; j++){
            if(strcmp(mod_code_with_context, tested_cases[j]) == 0){
                is_tested = true;
                break;
            }
        }
        if(!is_tested){
            WARNING("Modification code with context %s has not been tested.", mod_code_with_context);
        }
        free(mod_code_with_context);
    }
}

void print_varview_header(core_t* core) {
    if(core->opt.bedmethyl_out) return;
    char * common = "ref_contig\tref_pos\tstrand\tread_id\tread_pos\tmod_code\tmod_prob\tvar_pos\tvar_gt\tref_allele\talt_allele\toffset";
    char * haplotype = "";
    if(core->opt.haplotypes){
        haplotype = "\thaplotype";
    }
    fprintf(core->opt.output_fp, "%s%s\n", common, haplotype);
}

void print_varview_output(core_t* core, db_t* db) {
    FILE *out_fp = core->opt.output_fp;
    // int do_haplotypes = core->opt.haplotypes == 1;

    int is_bed = core->opt.bedmethyl_out;

    // Reusable buffer
    int max_arr_capacity = 0;
    varview_kv_t *sorted_arr = NULL;

    for(int i = 0; i < db->n_bam_recs; i++) {
        bam1_t *record = db->bam_recs[i];
        const char *qname = bam_get_qname(record);
        khash_t(varviewm) *varview_map = db->varview_maps[i];
        khint_t map_size = kh_size(varview_map);

        if (map_size == 0) continue;

        if (map_size > max_arr_capacity) {
            max_arr_capacity = map_size;
            sorted_arr = (varview_kv_t *)realloc(sorted_arr, sizeof(varview_kv_t) * max_arr_capacity);
            MALLOC_CHK(sorted_arr);
        }

        int size = 0;
        for (khint_t k = kh_begin(varview_map); k != kh_end(varview_map); k++) {
            if (kh_exist(varview_map, k)) {
                sorted_arr[size].key = (char *)kh_key(varview_map, k);
                sorted_arr[size].view = kh_value(varview_map, k);
                size++;
            }
        }

        // qsort(sorted_arr, size, sizeof(varview_kv_t), cmp_varview_kv);
        ks_introsort_varview(size, sorted_arr);

        for (int j = 0; j < size; j++) {
            varview_t* varview = sorted_arr[j].view;
            char *tname = NULL;
            int ref_pos;
            uint16_t ins_offset;
            char *mod_code;
            char strand;
            int haplotype;
            char * key = sorted_arr[j].key;
            decode_key(key, &tname, &ref_pos, &ins_offset, &mod_code, &strand, &haplotype);

            if(is_bed) {
                // contig \t start \t end \t mod_code \t mod_prob \t strand \t start \t end \t 255,0,0 \t 1 \t 1 \t 0 \t qname \t read_pos \t var_pos \t var_gt \t ref_allele \t alt_allele \t offset
                fprintf(out_fp, "%s\t%d\t%d\t%s\t%f\t%c\t%d\t%d\t255,0,0\t1\t1\t0\t%s\t%d\t%d\t%s\t%s\t%s\t%d", tname, ref_pos, ref_pos+1, mod_code, THRESH_UINT8_TO_DBL(varview->mod_prob), strand, ref_pos, ref_pos+1, qname, varview->read_pos, varview->var.pos, varview->var.gt ? varview->var.gt : ".", varview->var.ref_allele, varview->var.alt_allele, OFFSET_TO_INT(ins_offset));

                //print var.before_site
                fprintf(out_fp, "\t%s\t", varview->var.before_site);

                // print var.after_site
                fprintf(out_fp, "%s\t", varview->var.after_site);

                // print var.cg_offsets comma separated
                if(varview->var.cg_offsets_len > 0) {
                    fprintf(out_fp, "\t%d\t", varview->var.cg_offsets_len);
                    for(int o=0; o<varview->var.cg_offsets_len; o++) {
                        fprintf(out_fp, "%d,", varview->var.cg_offsets[o]);
                    }
                }
                fputc('\n', out_fp);
            } else {
                fprintf(out_fp, "%s\t%d\t%c\t%s\t%d\t%s\t%f\t%d\t%s\t%s\t%s\t%d\n", tname, ref_pos, strand, qname, varview->read_pos, mod_code, THRESH_UINT8_TO_DBL(varview->mod_prob), varview->var.pos, varview->var.gt ? varview->var.gt : ".", varview->var.ref_allele, varview->var.alt_allele, OFFSET_TO_INT(ins_offset));
            }
            free(tname);
            free(mod_code);
        }
    }

    if (sorted_arr) {
        free(sorted_arr);
    }
}

void merge_varfreq_maps(core_t* core, db_t* db) {
    khash_t(varfreqm) *core_map = core->varfreq_map;

    for (int i = 0; i < db->n_bam_recs; i++) {
        khash_t(varfreqm) *rec_map = db->varfreq_maps[i];

        if (kh_size(rec_map) == 0) continue;

        for (khint_t k = kh_begin(rec_map); k != kh_end(rec_map); ++k) {
            if (kh_exist(rec_map, k)) {
                char *key = (char *) kh_key(rec_map, k);
                varfreq_t *db_varfreq = kh_value(rec_map, k);

                int ret;
                khint_t core_k = kh_put(varfreqm, core_map, key, &ret);

                if (ret == 0) {
                    varfreq_t *core_varfreq = kh_value(core_map, core_k);
                    core_varfreq->n_called += db_varfreq->n_called;
                    core_varfreq->n_mod += db_varfreq->n_mod;
                } else {
                    kh_value(core_map, core_k) = db_varfreq;
                    kh_del(varfreqm, rec_map, k);
                }
            }
        }
    }
}

void destroy_varfreq_map(khash_t(varfreqm)* varfreq_map) {
    for (khint_t k = kh_begin(varfreq_map); k != kh_end(varfreq_map); k++) {
        if (kh_exist(varfreq_map, k)) {
            free((char *) kh_key(varfreq_map, k));
            free(kh_value(varfreq_map, k));
        }
    }
    kh_destroy(varfreqm, varfreq_map);
}

void print_varfreq_header(core_t* core) {
    if(core->opt.bedmethyl_out) return;
    char * common = "contig\tstart\tend\tstrand\tn_called\tn_mod\tfreq\tmod_code\tvar_pos\tvar_gt\tref_allele\talt_allele";
    char * hp_str = "";
    if(core->opt.haplotypes) hp_str = "\thaplotype";
    fprintf(core->opt.output_fp, "%s\toffset%s\n", common, hp_str);
}

void print_varfreq_output(core_t* core) {
    khash_t(varfreqm) *varfreq_map = core->varfreq_map;
    khint_t map_size = kh_size(varfreq_map);

    if (map_size == 0) return;

    double sort_start = realtime();
    varfreq_kv_t *sorted_arr = (varfreq_kv_t *)malloc(sizeof(varfreq_kv_t) * map_size);
    MALLOC_CHK(sorted_arr);
    int size = 0;
    for (khint_t k = kh_begin(varfreq_map); k != kh_end(varfreq_map); k++) {
        if (kh_exist(varfreq_map, k)) {
            sorted_arr[size].key = (char *)kh_key(varfreq_map, k);
            sorted_arr[size].freq = kh_value(varfreq_map, k);
            size++;
        }
    }
    ks_introsort_varfreq(size, sorted_arr);
    core->sort_time = realtime() - sort_start;

    double output_start = realtime();

    FILE *out_fp = core->opt.output_fp;
    int do_haplotypes = core->opt.haplotypes;

    int is_bed = core->opt.bedmethyl_out;
    char *agg_chrom = NULL;
    int agg_pos = -1;
    uint16_t agg_ins_offset = 0;
    char agg_strand = 0;
    char *agg_mod_code = NULL;
    uint64_t agg_n_called = 0, agg_n_mod = 0;
    int agg_count = 0;
    varfreq_t *agg_ref = NULL;

    for (int i = 0; i <= size; i++) {
        char *contig = NULL;
        int ref_pos = 0;
        uint16_t ins_offset = 0;
        char *mod_code = NULL;
        char strand = 0;
        int haplotype = 0;
        varfreq_t *varfreq = NULL;

        if (i < size) {
            varfreq = sorted_arr[i].freq;
            decode_key(sorted_arr[i].key, &contig, &ref_pos, &ins_offset, &mod_code, &strand, &haplotype);
        }

        //  - haplotypes on : group by (contig,pos,strand,mod,offset), sum over haplotype -> haplotype=* row
        //  - haplotypes off: group by (contig,pos,strand,mod), sum over offset            -> offset=* row
        //    (haplotype is always -1 when off, so it is not part of the grouping)
        int is_new_group = (i == size) || !agg_chrom ||
            strcmp(contig, agg_chrom) != 0 || ref_pos != agg_pos ||
            strand != agg_strand || strcmp(mod_code, agg_mod_code) != 0 ||
            (do_haplotypes && ins_offset != agg_ins_offset);

        //  - haplotypes on : >=2 haplotypes at this (pos,offset) -> haplotype=* row
        //  - haplotypes off: >=2 offsets at this position        -> offset=* row
        if (is_new_group && agg_count >= 2) {
            double avg_freq = (double)agg_n_mod / agg_n_called;
            if (is_bed) {
                int end = agg_pos + 1;
                fprintf(out_fp, "%s\t%d\t%d\t%s\t%llu\t%c\t%d\t%d\t255,0,0\t%llu\t%f\t%d\t%s\t%s\t%s\t",
                    agg_chrom, agg_pos, end, agg_mod_code,
                    (unsigned long long)agg_n_called, agg_strand, agg_pos, end,
                    (unsigned long long)agg_n_called, avg_freq * 100, agg_ref->var_pos,
                    agg_ref->gt ? agg_ref->gt : ".",
                    agg_ref->ref_allele ? agg_ref->ref_allele : ".",
                    agg_ref->alt_allele ? agg_ref->alt_allele : ".");
                // haplotypes on: real offset + haplotype=* ; off: offset=*
                if (do_haplotypes) fprintf(out_fp, "%d\t*\n", OFFSET_TO_INT(agg_ins_offset));
                else fputs("*\n", out_fp);
            } else {
                fprintf(out_fp, "%s\t%d\t%d\t%c\t%llu\t%llu\t%f\t%s\t%d\t%s\t%s\t%s\t",
                    agg_chrom, agg_pos, agg_pos + 1, agg_strand,
                    (unsigned long long)agg_n_called, (unsigned long long)agg_n_mod, avg_freq,
                    agg_mod_code, agg_ref->var_pos,
                    agg_ref->gt ? agg_ref->gt : ".",
                    agg_ref->ref_allele ? agg_ref->ref_allele : ".",
                    agg_ref->alt_allele ? agg_ref->alt_allele : ".");
                if (do_haplotypes) fprintf(out_fp, "%d\t*\n", OFFSET_TO_INT(agg_ins_offset));
                else fputs("*\n", out_fp);
            }
        }

        if (i == size) break;

        if (is_new_group) {
            free(agg_chrom); agg_chrom = contig; contig = NULL;
            free(agg_mod_code); agg_mod_code = mod_code; mod_code = NULL;
            agg_pos = ref_pos; agg_strand = strand;
            agg_ins_offset = ins_offset;
            agg_n_called = 0; agg_n_mod = 0; agg_count = 0; agg_ref = varfreq;
        }
        agg_n_called += varfreq->n_called;
        agg_n_mod += varfreq->n_mod;
        agg_count++;

        double freq_value = (double)varfreq->n_mod / varfreq->n_called;
        if (is_bed) {
            int end = ref_pos + 1;
            fprintf(out_fp, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t255,0,0\t%d\t%f\t%d\t%s\t%s\t%s\t%d",
                agg_chrom, ref_pos, end, agg_mod_code, varfreq->n_called, strand, ref_pos, end,
                varfreq->n_called, freq_value * 100,
                varfreq->var_pos,
                varfreq->gt ? varfreq->gt : ".",
                varfreq->ref_allele ? varfreq->ref_allele : ".",
                varfreq->alt_allele ? varfreq->alt_allele : ".",
                OFFSET_TO_INT(ins_offset));
            if (do_haplotypes) {
                if (haplotype == -1) fputs("\t*", out_fp);
                else fprintf(out_fp, "\t%d", haplotype);
            }
            fputc('\n', out_fp);
        } else {
            fprintf(out_fp, "%s\t%d\t%d\t%c\t%d\t%d\t%f\t%s\t%d\t%s\t%s\t%s\t%d",
                agg_chrom, ref_pos, ref_pos + 1, strand,
                varfreq->n_called, varfreq->n_mod, freq_value, agg_mod_code,
                varfreq->var_pos,
                varfreq->gt ? varfreq->gt : ".",
                varfreq->ref_allele ? varfreq->ref_allele : ".",
                varfreq->alt_allele ? varfreq->alt_allele : ".",
                OFFSET_TO_INT(ins_offset));
            if (do_haplotypes) {
                if (haplotype == -1) fputs("\t*", out_fp);
                else fprintf(out_fp, "\t%d", haplotype);
            }
            fputc('\n', out_fp);
        }
        free(contig);
        free(mod_code);
    }

    free(agg_chrom);
    free(agg_mod_code);

    if(out_fp != stdout) fclose(out_fp);

    free(sorted_arr);
    core->output_time += realtime() - output_start;
}
