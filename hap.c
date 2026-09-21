/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2026 Chenxi Zhou <chnx.zhou@gmail.com>                          *
 *                                                                               *
 * Permission is hereby granted, free of charge, to any person obtaining a copy  *
 * of this software and associated documentation files (the "Software"), to deal *
 * in the Software without restriction, including without limitation the rights  *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
 * copies of the Software, and to permit persons to whom the Software is         *
 * furnished to do so, subject to the following conditions:                      *
 *                                                                               *
 * The above copyright notice and this permission notice shall be included in    *
 * all copies or substantial portions of the Software.                           *
 *                                                                               *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE *
 * SOFTWARE.                                                                     *
 *********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>
#include <pthread.h>

#include "kvec.h"
#include "kthread.h"

#include "hic.h"
#include "sdict.h"
#include "range.h"
#include "agp-spec.h"
#include "misc.h"

#include "interfaces/highs_c_api.h"

#include "busco.h"
#include "overlap.h"
#include "ploidy.h"
#include "hap.h"

static AGP_CT_t LG_AGP_SEQ_COMPONENT_TYPE = AGP_CT_W;
static AGP_CT_t LG_AGP_GAP_COMPONENT_TYPE = AGP_CT_N;
static AGP_LE_t LG_AGP_LINKAGE_EVIDENCE = AGP_LE_ALIGN_GENUS;
static int LG_AGP_GAP_SIZE = DEFAULT_AGP_U_GAP_SIZE;

static inline void write_agp_seq(FILE *fo, char *s_name, uint64 s_beg, uint64 s_end, uint32 b,
    char *c_name, uint32 c_beg, uint32 c_end, uint32 c_oris)
{
    fprintf(fo, "%s\t%llu\t%llu\t%u\t%s\t%s\t%u\t%u\t%s\n", s_name, s_beg, s_end, b,
        agp_component_type_val(LG_AGP_SEQ_COMPONENT_TYPE), c_name, c_beg, c_end, 
        agp_orientation_val(c_oris? AGP_OT_MINUS : AGP_OT_PLUS));
}

static inline void write_agp_gap(FILE *fo, char *s_name, uint64 s_beg, uint64 s_end, uint32 b)
{
    fprintf(fo, "%s\t%llu\t%llu\t%u\t%s\t%d\t%s\t%s\t%s\n", s_name, s_beg, s_end, b,
        agp_component_type_val(LG_AGP_GAP_COMPONENT_TYPE), 
        LG_AGP_GAP_SIZE, 
        agp_gap_type_val(AGP_GT_SCAFFOLD),
        agp_linkage_val(AGP_LG_YES), 
        agp_linkage_evidence_val(LG_AGP_LINKAGE_EVIDENCE));
}

typedef struct { int read, beg, end; } seg_t;
typedef kvec_t(seg_t) seg_vec_t;

static int seg_acmpfunc(const void *a, const void *b)
{
    int x = ((seg_t *) a)->read;
    int y = ((seg_t *) b)->read;
    
    if (x == y) {
        x = ((seg_t *) a)->beg;
        y = ((seg_t *) b)->beg;
    }

    return (x > y) - (x < y);
}

void scf_free(scf_t *scf)
{
    if (!scf) return;
    free(scf->ctgs);
    free(scf->segs);
}

static int scaff_size_dcmpfunc(const void *a, const void *b)
{
    int64 xm, ym;

    xm = ((scf_t *) a)->len;
    ym = ((scf_t *) b)->len;

    return ((xm < ym) - (xm > ym));
}

int scaff_natural_cmpfunc(const void *a, const void *b)
{
    const scf_t *x = (const scf_t *) a;
    const scf_t *y = (const scf_t *) b;

    if (x->grp >= 0 && y->grp < 0) return -1;
    if (x->grp < 0 && y->grp >= 0) return  1;
    
    if (x->grp >= 0 && y->grp >= 0 && x->grp != y->grp)
        return x->grp - y->grp;

    if (x->hap >= 0 && y->hap < 0) return -1;
    if (x->hap < 0 && y->hap >= 0) return  1;
    
    if (x->hap >= 0 && y->hap >= 0 && x->hap != y->hap)
        return x->hap - y->hap;
    
    return (x->len < y->len) - (x->len > y->len);
}

static void busco_summary_report_segments(busco_table_t *buscos, sdict_t *dicts, seg_t *segs, int nseg, ord_u64_t *ords, int print_summary, int max_hist)
{
    if (!buscos) return;

    int i, which;
    seq_range_t *range;
    kvec_t(seq_range_t) ranges = {0, 0, 0};

    for (i = 0; i < nseg; i++) {
        which = ords? ords[i].which : i;
        kv_pushp(seq_range_t, ranges, &range);
        range->seq = segs[which].read;
        range->beg = segs[which].beg;
        range->end = segs[which].end;
    }
    busco_summary_report(buscos, ranges.a, ranges.n, print_summary, max_hist);
    kv_destroy(ranges);
}

static void busco_summary_report_scaffolds(busco_table_t *buscos, sdict_t *dicts, scf_t *scfs, int nscf, ord_u64_t *ords, int print_summary, int max_hist)
{
    if (!buscos) return;

    int i, j, which, nctg;
    uint32 v, *ctgs;
    uint64 *segs;
    seq_range_t *range;
    kvec_t(seq_range_t) ranges = {0, 0, 0};

    for (i = 0; i < nscf; i++) {
        which = ords? ords[i].which : i;
        ctgs = scfs[which].ctgs;
        segs = scfs[which].segs;
        nctg = scfs[which].nctg;
        for (j = 0; j < nctg; j++) {
            v = ctgs[j];
            kv_pushp(seq_range_t, ranges, &range);
            range->seq = v>>1;
            if (v&1) {
                range->beg = dicts->s[v>>1].len - (segs[j]>>32) - (uint32)segs[j];
                range->end = dicts->s[v>>1].len - (segs[j]>>32);
            } else {
                range->beg = segs[j]>>32;
                range->end = (segs[j]>>32) + (uint32)segs[j];
            }
        }
    }
    busco_summary_report(buscos, ranges.a, ranges.n, print_summary, max_hist);
    kv_destroy(ranges);
}

static void haplotype_group_sequence_summary(scf_t *scfs, int nscf, sdict_t *dicts, busco_table_t *buscos, int ploidy, int long_format, FILE *out)
{
    if (!nscf) return;

    int64 i, j, k, slen, nctg;
    int hap;
    ord_u64_t *ords;

    if (out == NULL) out = stderr;

    MYMALLOC(ords, nscf);
    for (i = 0; i < nscf; i++) {
        ords[i].which = i;
        ords[i].event = (uint64) scfs[i].hap;
    }
    qsort(ords, nscf, sizeof(ord_u64_t), ord_u64_acmpfunc);
    for (i = j = 0; i <= nscf; i++) {
        if (i == nscf || ords[i].event != ords[j].event) {
            hap = (int64) ords[j].event;
            slen = 0;
            nctg = 0;
            for (k = j; k < i; k++) {
                nctg += scfs[ords[k].which].nctg;
                slen += scfs[ords[k].which].len;
            }
            if (long_format) {
                if (hap < 0)
                    fprintf(out, "[M::%s] %8lld sequences of %12lld bp not assigned\n", __func__, nctg, slen);
                else
                    fprintf(out, "[M::%s] %8lld sequences of %12lld bp assigned to haplotype group %d\n", __func__, nctg, slen, hap+1);
                // bsuco stats
                busco_summary_report_scaffolds(buscos, dicts, scfs, i-j, ords+j, 1, ploidy);
            } else {
                if (hap < 0)
                    fprintf(out, "[M::%s] %8lld sequences of %12lld bp not assigned", __func__, nctg, slen);
                else
                    fprintf(out, "[M::%s] %8lld sequences of %12lld bp assigned to haplotype group %d", __func__, nctg, slen, hap+1);
                if (buscos) {
                    busco_summary_report_scaffolds(buscos, dicts, scfs, i-j, ords+j, 0, 0);
                    fprintf(out, " [%s]\n", BUSCO_SUMMARY_STRING);
                } else
                    fputc('\n', out);
            }
            j = i;
        }
    }
    free(ords);
}

static double range_coverage(ovl_t *ovls, int64 novl, int beg, int end)
{
    int64 i;
    double cov = 0.;

    for (i = 0; i < novl; i++) {
        if (ovls[i].del)
            continue;
        if (ovls[i].abpos >= end || ovls[i].aepos <= beg)
            continue;
        cov += MIN(ovls[i].aepos, end) - MAX(ovls[i].abpos, beg);
    }

    return cov;
}

static ovl_t *find_b_overlap(ovl_t *ovls, int size, uint32 b)
{
    int left, mid, right;
    left = 0;
    right = size - 1;

    while (left <= right) {
        mid = left + (right - left) / 2;
        if (ovls[mid].bread == b)
            return ovls + mid;
        if (ovls[mid].bread < b)
            left = mid + 1;
        else right = mid - 1;
    }
    return NULL;
}

static int overlap_accorss_position(ovl_t *ovls, uint64 *index, uint32 s1, int p1, uint32 s2, int p2, int min_ext)
{
    ovl_t *ovl = find_b_overlap(ovls + (index[s1>>1] >> 32), (uint32) index[s1>>1], s2>>1);
    if (ovl == NULL || 
        ovl->del || 
        (ovl->arev == ovl->brev) != ((s1&1) == (s2&1)) ||
        p1 - ovl->abpos < min_ext || 
        p2 - ovl->bbpos < min_ext ||
        ovl->aepos - p1 < min_ext || 
        ovl->bepos - p2 < min_ext)
        return 0;
    return 1;
}

typedef struct {
    int a, b;
    double l, v, p;
} hlk_t;

static int hlk_abseq_cmpfunc(const void *a, const void *b)
{
    hlk_t *x = (hlk_t *) a;
    hlk_t *y = (hlk_t *) b;

    if (x->a != y->a)
        return (x->a > y->a) - (x->a < y->a);
    return (x->b > y->b) - (x->b < y->b);
}

static int hlk_pval_cmpfunc(const void *a, const void *b)
{
    hlk_t *x = (hlk_t *) a;
    hlk_t *y = (hlk_t *) b;

    if (x->p != y->p)
        return (x->p > y->p) - (x->p < y->p);
    return (x->v < y->v) - (x->v > y->v);
}

typedef struct {
    int a, b;
    int64 s, e;
} hidx_t;

// same ordering as hic_cmpfunc in hic.c
static int hic_link_cmpfunc(const void *a, const void *b)
{
    const hic_t *x = a, *y = b;
    if (x->aseq != y->aseq)
        return (x->aseq > y->aseq) - (x->aseq < y->aseq);
    if (x->bseq != y->bseq)
        return (x->bseq > y->bseq) - (x->bseq < y->bseq);
    if (x->apos != y->apos)
        return (x->apos > y->apos) - (x->apos < y->apos);
    return (x->bpos > y->bpos) - (x->bpos < y->bpos);
}

static hidx_t *build_hidx(hic_t *hics, int64 nhic, int64 *_nidx)
{
    if (nhic == 0) {
        *_nidx = 0;
        return NULL;
    }

    kvec_t(hidx_t) hidx;
    hidx_t *idx;
    int64 i;
    int a, b;

    kv_init(hidx);
    kv_resize(hidx_t, hidx, 1<<20);

    a = hics[0].aseq;
    b = hics[0].bseq;
    kv_pushp(hidx_t, hidx, &idx);
    idx->a = a;
    idx->b = b;
    idx->s = 0;
    for (i = 1; i < nhic; i++) {
        if (hics[i].aseq != a || hics[i].bseq != b) {
            idx->e = i;
            a = hics[i].aseq;
            b = hics[i].bseq;
            kv_pushp(hidx_t, hidx, &idx);
            idx->a = a;
            idx->b = b;
            idx->s = i;
        }
    }
    idx->e = nhic;

    MYREALLOC(hidx.a, hidx.n);
    *_nidx = hidx.n;

    return hidx.a;
}

static int64 lower_bound_apos(hic_t *hics, int64 low, int64 high, int apos)
{
    int64 mid;
    while (low < high) {
        mid = low + (high - low) / 2;
        if (hics[mid].apos < apos) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    return low;
}

static int64 upper_bound_apos(hic_t *hics, int64 low, int64 high, int apos)
{
    int64 mid;
    while (low < high) {
        mid = low + (high - low) / 2;
        if (hics[mid].apos <= apos) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    return low;
}

// binary search for the (a, b, [ab, ae])
static int64 hidx_bsearch(hic_t *hics, hidx_t *hidx, int64 nidx, int a, int b, int ab, int ae, 
    int64 *_s, int64 *_e)
{
    *_s = 0;
    *_e = 0;
    
    if (nidx == 0 || ab > ae) { 
        return 0;
    }

    int64 low, mid, high, block_start, block_end;
    
    low = 0;
    high = nidx - 1;
    block_start = block_end = -1;
    // find a-b block
    while (low <= high) {
        mid = low + (high - low) / 2;
        if (hidx[mid].a == a && hidx[mid].b == b) {
            block_start = hidx[mid].s;
            block_end = hidx[mid].e;
            break;
        } else if (hidx[mid].a < a || 
            (hidx[mid].a == a && hidx[mid].b < b)) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    // block not found
    if (block_start == -1)
        return 0;

    // find first element >= ab
    *_s = lower_bound_apos(hics, block_start, block_end, ab);

    // find first element > ae
    *_e = upper_bound_apos(hics, *_s, block_end, ae);

    return *_e - *_s;
}

static int64 count_cis_links(hic_t *hics, int64 s, int64 e, uint32 ab, uint32 ae, uint32 bb, uint32 be)
{
    int64 i, cnt;
    uint32 apos, bpos;
    hic_t *h;

    cnt = 0;
    h = hics + s;
    for (i = s; i < e; i++, h++) {
        apos = h->apos;
        bpos = h->bpos;
        // find if overlaps with the given positions
        if (ab <= apos && apos <= ae && 
            bb <= bpos && bpos <= be)
            cnt += h->nhic;
    }

    return cnt;
}

static int64 count_trans_links(hic_t *hics, hidx_t *hidx, int64 nidx,
    int a, int a_wb, int a_we, ovl_t *ovls, uint32 novl, int b_wb, int b_we)
{
    int64 i, s, e, cnt;
    int b, bb, be;

    cnt = 0;
    for (i = 0; i < novl; i++, ovls++) {
        if (ovls->del)
            continue;
        b = ovls->bread;
        if (b == a)
            continue;
        // skip overlaps that don't touch b's window range [b_wb, b_we)
        if (ovls->aepos / HIC_NORM_WINDOW < b_wb || ovls->abpos / HIC_NORM_WINDOW >= b_we)
            continue;
        bb = ovls->bbpos / HIC_NORM_WINDOW;
        be = ovls->bepos / HIC_NORM_WINDOW;
        if (hidx_bsearch(hics, hidx, nidx, b, a, bb, be, &s, &e) > 0)
            cnt += count_cis_links(hics, s, e, bb, be, a_wb, a_we);
    }

    return cnt;
}

typedef struct {
    int64 which;
    int b, bb, be;
    int cnt;
} ovl_cnt_t;

static int ovl_cnt_b_cmpfunc(const void *a, const void *b)
{
    ovl_cnt_t *x = (ovl_cnt_t *) a;
    ovl_cnt_t *y = (ovl_cnt_t *) b;
    if (x->b != y->b)
        return (x->b > y->b) - (x->b < y->b);
    if (x->bb != y->bb)
        return (x->bb > y->bb) - (x->bb < y->bb);
    return (x->be > y->be) - (x->be < y->be);
}

static int ovl_cnt_w_cmpfunc(const void *a, const void *b)
{
    ovl_cnt_t *x = (ovl_cnt_t *) a;
    ovl_cnt_t *y = (ovl_cnt_t *) b;
    return (x->which > y->which) - (x->which < y->which);
}

static void count_olaps(range_t *rngs, uint32 nrng, ovl_cnt_t *ovls, uint32 novl)
{
    if (rngs == NULL || nrng == 0 || 
        ovls == NULL || novl == 0)
        return;
    
    uint32 i, j, s;
    int bb, be, c;

    for (i = s = 0; i < novl; i++) {
        if (ovls[i].cnt < 0) {
            ovls[i].cnt = 0;
            continue; // skip deleted overlaps
        }

        bb = ovls[i].bb;
        be = ovls[i].be + 1;

        while (s < nrng && rngs[s].end <= bb) s++;

        for (j = s; j < nrng; j++) {
            if (rngs[j].beg >= be)
                break;
            c = MIN(be, rngs[j].end) - MAX(bb, rngs[j].beg);
            if (c > 0) ovls[i].cnt += c;
        }
    }
}

static int64 count_trans_olaps(ovl_cnt_t *ovls, uint32 novl)
{
    int64 i, cnt;
    
    cnt = 0;
    for (i = 0; i < novl; i++)
        cnt += ovls[i].cnt;
    return cnt;
}

// Computes the upper-tail p-value P(X >= k) for a Binomial(n, p) distribution
static double calculate_binomial_pvalue(int k, int n, double p)
{
    if (k > n) return 0.0;
    if (k <= 0) return 1.0;
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return 1.0;

    // If our observed links (k) are less than or equal to the expected 
    // random links (n * p), we are on the left side of the bell curve.
    // The p-value is guaranteed to be > 0.5. Since we only care about 
    // highly significant pulls (e.g., p < 0.01), we can instantly bail out.
    if (k <= n * p) return 1.0; 

    double log_pmf_k = lgamma(n + 1.0) 
                     - lgamma(k + 1.0) 
                     - lgamma(n - k + 1.0) 
                     + k * log(p) 
                     + (n - k) * log(1.0 - p);
    
    double current_pmf = exp(log_pmf_k);
    double p_value = current_pmf;

    for (int i = k; i < n; ++i) {
        current_pmf = current_pmf * (n - i) / (i + 1.0) * (p / (1.0 - p));
        
        if (current_pmf / p_value < DBL_EPSILON)
            break;

        p_value += current_pmf;
    }
    
    return (p_value > 1.0) ? 1.0 : p_value;
}

// Computes the upper-tail p-value P(X >= k) for a Poisson(lambda) distribution
static double calculate_poisson_pvalue(int k, double lambda)
{
    if (k <= 0)        return 1.0;
    if (lambda <= 0.0) return 0.0;

    const double sqrt_lambda = sqrt(lambda);
    const double dk = (double) k;
    
    if (dk > lambda + 10.0 * sqrt_lambda) {
        double log_term = dk * log(lambda) - lambda - lgamma(dk + 1.0);
        if (log_term < -708.0) return 0.0; // underflows to 0 in double

        double term = exp(log_term);
        double sum  = term;
        int x;
        for (x = k + 1; term > DBL_EPSILON * sum; x++) {
            term *= lambda / (double)x;
            sum  += term;
        }
        return sum > 1.0 ? 1.0 : sum;
    }

    if (dk < lambda - 10.0 * sqrt_lambda) {
        double term = exp(-lambda);
        double cdf  = term;
        int x;
        for (x = 1; x < k; x++) {
            term *= lambda / (double)x;
            cdf  += term;
        }
        double sf = 1.0 - cdf;
        return sf < 0.0 ? 0.0 : sf;
    }

    double log_mode = (double)(k - 1) * log(lambda) - lambda - lgamma(dk);
    double log_px, cdf_sum = 0.0;
    int x;
    for (x = 0; x < k; x++) {
        log_px = (double)x * log(lambda) - lambda - lgamma((double)x + 1.0);
        cdf_sum += exp(log_px - log_mode);
    }

    double log_cdf = log(cdf_sum) + log_mode;
    if (log_cdf >= 0.0) return 0.0;

    double log_sf = log1p(-exp(log_cdf));
    double sf     = exp(log_sf);

    return sf < 0.0 ? 0.0 : (sf > 1.0 ? 1.0 : sf);
}

static const int MIN_HIC_SIZE = 10000;
static const int MIN_HIC_WINDOW = MIN_HIC_SIZE / HIC_NORM_WINDOW;
static const int MIN_CIS_HLINK = 3;
static const double MAX_MODIFIED_Z = 3.5;
static const double PVAL_THRESHOLD = 0.01;

static void test_hic_linkage_map(hlk_t *hlks, int64 nhlk, ovl_t *ovls, int64 novl, sdict_t *dicts)
{
    hlk_t *hlk;
    ovl_t *ovl;
    ord_dbl_t *hords;
    uint64 *hindex, *oindex;
    int a, b, g, n, qs, qe, nseq, *smarks, *sgroup, *squeue;
    int64 i, j, k, slen, olen;

    // number of sequences
    nseq = dicts->n;

    // build index for quick search
    MYCALLOC(hindex, nseq);
    MYCALLOC(oindex, nseq);
    MYCALLOC(smarks, nseq);
    MYCALLOC(sgroup, nseq);
    MYCALLOC(squeue, nseq);
    MYMALLOC(hords, nhlk);
    if (hindex == NULL || oindex == NULL || smarks == NULL || 
        sgroup == NULL || squeue == NULL || hords == NULL)
        mem_alloc_error("hic linkage map index");

    // build index for hic linkages
    a = hlks->a;
    for (i = 1, j = 0; i < nhlk; i++) {
        if (hlks[i].a != a) {
            hindex[a] = (uint64) j << 32 | (i - j);
            a = hlks[i].a;
            j = i;
        }
    }
    hindex[a] = (uint64) j << 32 | (i - j);
    
    // build index for overlaps
    a = ovls->aread;
    for (i = 1, j = 0; i < novl; i++) {
        if (ovls[i].aread != a) {
            oindex[a] = (uint64) j << 32 | (i - j);
            a = ovls[i].aread;
            j = i;
        }
    }
    oindex[a] = (uint64) j << 32 | (i - j);

    // hic linkages
    for (i = 0; i < nhlk; i++) {
        hlk = hlks + i;
        fprintf(stdout, "H %s\t%12d\t%s\t%12d\t%.3f\t%.3f\t%.3e\n", 
            dicts->s[hlk->a].name, dicts->s[hlk->a].len, 
            dicts->s[hlk->b].name, dicts->s[hlk->b].len, 
            hlk->l, hlk->v, hlk->p);
    }

    // group sequences by hic linkages
    MYBZERO(sgroup, nseq);
    MYBZERO(smarks, nseq);
    g = 0;
    for (i = 0; i < nseq; i++) {
        if (sgroup[i])
            continue;
        qs = qe = 0;
        squeue[qe++] = i;
        smarks[i] = 1;
        while (qs < qe) {
            a = squeue[qs++];
            hlk = hlks + (hindex[a] >> 32);
            n = (uint32) hindex[a];
            for (j = 0; j < n; j++, hlk++) {
                b = hlk->b;
                if (!smarks[b]) {
                    squeue[qe++] = b;
                    smarks[b] = 1;
                }
            }
        }
        g++;
        slen = 0;
        for (j = 0; j < qe; j++) {
            a = squeue[j];
            sgroup[a] = g;
            slen += dicts->s[a].len;
        }
        olen = 0;
        for (j = 0; j < qe; j++) {
            a = squeue[j];
            ovl = ovls + (oindex[a] >> 32);
            n = (uint32) oindex[a];
            for (k = 0; k < n; k++, ovl++) {
                b = ovl->bread;
                if (sgroup[b] == g)
                    olen += ovl->bepos - ovl->bbpos;
            }
        }
        olen /= 2; // each overlap is counted twice

        fprintf(stdout, "G %6d %6d %12lld %12lld\n", g, qe, slen, olen);
    }

    free(hindex);
    free(oindex);
    free(smarks);
    free(sgroup);
    free(squeue);
    free(hords);
}

static hlk_t *build_hic_linkage_map(const char *hic_bfile, int min_qual, sdict_t *dicts, asm_dict_t *break_dict, ovl_t *ovls, int64 novl, int64 *_nhlk)
{
    if (_nhlk) *_nhlk = 0;
    if (hic_bfile == NULL)
        return NULL;

    kvec_t(range_t) rngs;
    kvec_t(hlk_t) hlks;
    range_t *rng, *srngs;
    uint64 *index, *ridx;
    hidx_t *hidx, *idx;
    ovl_t *ov;
    ovl_cnt_t *oc, *tcnts;
    ord_dbl_t *sords;
    hic_t *hics;
    uint8 *smark;
    int64 i, j, n, s, e, nseq, mseq, nidx, stot, htot, nhlk, nhic, *wcnts, *hcnts, *hcnts_beg, *hcnts_end;
    int64 acnt, tcnt, hc_a, hc_b, cish_ab, cisl_a, cisl_b, transh_a, transh_b, transl_a, transl_b;
    double h_pb, p1, p2, p3, cis_l, d_val, p_val, p_thresh, nf, median;
    int a, b, p, w, *clens, *clens_beg, *clens_end, *tlens, *tlens_beg, *tlens_end;
    int a_wb, a_we, b_wb, b_we, best_ea, best_eb;
    int a_wl, b_wl, a_st, b_st, ea, eb, ea_wb, ea_we, eb_wb, eb_we;
    int Wend;

#undef TEST_HIC
#ifdef TEST_HIC
    // read HiC links                                                        
    FILE *fp = fopen(hic_bfile, "r");                                        
    fread(&nhic, sizeof(int64), 1, fp);                                      
    MYMALLOC(hics, nhic);
    if (hics == NULL)
        mem_alloc_error("HiC links array");
    fread(hics, sizeof(hic_t), nhic, fp);
    fclose(fp);
#else
    // read HiC links
    // with a break AGP the binary file is on the raw sequences
    // coordinates are converted to piece coordinates while reading
    nhic = 0;
    hics = break_dict?
        read_hic_from_binary_sd_conversion((char *) hic_bfile, break_dict, HIC_NORM_WINDOW, min_qual, &nhic) :
        read_hic_from_binary((char *) hic_bfile, dicts, HIC_NORM_WINDOW, min_qual, &nhic);
    if (break_dict && hics && nhic > 0) {
        // the sd_conversion reader returns upper-triangle links only
        // add symmetric links to match read_hic_from_binary
        int64 m = nhic;
        MYREALLOC(hics, m * 2);
        if (hics == NULL)
            mem_alloc_error("HiC links array");
        for (i = 0; i < m; i++)
            if (hics[i].aseq != hics[i].bseq)
                hics[nhic++] = (hic_t){hics[i].bseq, hics[i].bpos, hics[i].aseq, hics[i].apos, hics[i].nhic};
        MYREALLOC(hics, nhic);
        qsort(hics, nhic, sizeof(hic_t), hic_link_cmpfunc);
    }
#endif

    if (hics == NULL || nhic == 0) {
        fprintf(stderr, "[M::%s] no valid HiC links found in the input file\n", __func__);
        return NULL;
    }
    
    // sequence number
    nseq = dicts->n;
    
    // maximum sequence length in windows
    mseq = 0;
    for (i = 0; i < nseq; i++) {
        s = (dicts->s[i].len - 1) / HIC_NORM_WINDOW + 1;
        if (s > mseq) mseq = s;
    }

    // allocate memory
    MYCALLOC(tcnts, novl);
    MYCALLOC(wcnts, mseq);
    MYCALLOC(index, nseq);
    MYCALLOC(hcnts, nseq);
    MYCALLOC(hcnts_beg, nseq);
    MYCALLOC(hcnts_end, nseq);
    MYCALLOC(clens, nseq);
    MYCALLOC(clens_beg, nseq);
    MYCALLOC(clens_end, nseq);
    MYCALLOC(tlens, nseq);
    MYCALLOC(tlens_beg, nseq);
    MYCALLOC(tlens_end, nseq);
    MYCALLOC(smark, nseq);
    MYCALLOC(sords, nseq);
    MYCALLOC(ridx, nseq);

    if (tcnts == NULL || wcnts == NULL || index == NULL || 
        hcnts == NULL || hcnts_beg == NULL || hcnts_end == NULL ||
        clens == NULL || clens_beg == NULL || clens_end == NULL ||
        tlens == NULL || tlens_beg == NULL || tlens_end == NULL ||
        smark == NULL || sords == NULL || ridx == NULL)
        mem_alloc_error("hic index array");

    // build overlap index for quick search
    // overlaps are sorted by abseqs
    a = ovls->aread;
    for (i = 1, j = 0; i < novl; i++) {
        if (ovls[i].aread != a) {
            index[a] = (uint64) j << 32 | (i - j);
            a = ovls[i].aread;
            j = i;
        }
    }
    index[a] = (uint64) j << 32 | (i - j);

    // build hic index for quick search
    nidx = 0;
    hidx = build_hidx(hics, nhic, &nidx);

    // estimate Wend
    for (i = 0; i < nhic; i++) {
        if (hics[i].aseq != hics[i].bseq)
            continue;
        wcnts[hics[i].bpos-hics[i].apos] += hics[i].nhic;
    }
    acnt = tcnt = 0;
    for (i = 0; i < mseq; i++)
        tcnt += wcnts[i];
    fprintf(stderr, "[M::%s] total intra hic links: %lld\n", __func__, tcnt);
    Wend = 0;
    nf = (double)((int) ((double) (wcnts[0] - 1) / tcnt * 10)) / 10 + .1;
    for (i = 0; i < mseq; i++) {
        acnt += wcnts[i];
        if (acnt >= nf * tcnt) {
            fprintf(stderr, "[M::%s] %10lld %10lld %10.3f%%\n", __func__, i, acnt, (double) acnt / tcnt * 100);
            nf += .1;
            if (nf < 1.) Wend = i; // using 90% intra links to estimate the end window size
        }
    }
    if (Wend < 100000 / HIC_NORM_WINDOW) // 100kb minimum end window size
        Wend = 100000 / HIC_NORM_WINDOW;
    fprintf(stderr, "[M::%s] using hic window size: %d bp\n", __func__, HIC_NORM_WINDOW);
    fprintf(stderr, "[M::%s] using terminal window: %d kb\n", __func__, Wend * HIC_NORM_WINDOW / 1000);

    // collect windows covered by hic links
    kv_init(rngs);
    kv_resize(range_t, rngs, 1<<20);
    a = hidx[0].a;
    b = rngs.n;
    a_wl = dicts->s[a].len / HIC_NORM_WINDOW;
    for (i = 0; i <= nidx; i++) {
        if (i == nidx || hidx[i].a != a) {
            p = merge_ranges(rngs.a + b, rngs.n - b, 0);
            ridx[a] = (uint64) b << 32 | p;
            rngs.n = b + p;

            if (i == nidx) break;
            a = hidx[i].a;
            b = rngs.n;
            a_wl = dicts->s[a].len / HIC_NORM_WINDOW;
        }
        if (a == hidx[i].b)
            continue;
        s = hidx[i].s;
        e = hidx[i].e;
        kv_pushp(range_t, rngs, &rng);
        p = hics[s].apos;
        rng->beg = p;
        hcnts[a] += hics[s].nhic;
        if (hics[s].apos < Wend)
            hcnts_beg[a] += hics[s].nhic;
        if (hics[s].apos >= a_wl - Wend)
            hcnts_end[a] += hics[s].nhic;
        for (j = s + 1; j < e; j++) {
            if (hics[j].apos != ++p) {
                rng->end = p;
                kv_pushp(range_t, rngs, &rng);
                p = hics[j].apos;
                rng->beg = p;
            }
            hcnts[a] += hics[j].nhic;
            if (hics[j].apos < Wend)
                hcnts_beg[a] += hics[j].nhic;
            if (hics[j].apos >= a_wl - Wend)
                hcnts_end[a] += hics[j].nhic;
        }
        rng->end = ++p;
    }
    srngs = rngs.a;

    // collect stats
    stot = htot = 0;
    for (i = 0; i < nseq; i++) {
        a_wl = dicts->s[i].len / HIC_NORM_WINDOW;
        a_wl = MIN(a_wl, Wend);
        b_wl = dicts->s[i].len / HIC_NORM_WINDOW - Wend;
        b_wl = MAX(b_wl, 0);
        rng = srngs + (ridx[i] >> 32);
        p = (uint32) ridx[i];
        clens[i] = clens_beg[i] = clens_end[i] = 0;
        while (p--) {
            a_wb = rng->beg, a_we = rng->end;
            clens[i] += a_we - a_wb;
            b_wb = a_wb; b_we = MIN(a_we, a_wl);
            if (b_wb < b_we) clens_beg[i] += b_we - b_wb;
            b_wb = MAX(a_wb, b_wl); b_we = a_we;
            if (b_wb < b_we) clens_end[i] += b_we - b_wb;
            rng++;
        }
        stot += clens[i];
        htot += hcnts[i];
    }
    htot /= 2; // each hic is counted twice
    h_pb = (double) htot / stot / stot ; // hic link per base density
    fprintf(stderr, "[M::%s] total windows covered: %lld\n", __func__, stot);
    fprintf(stderr, "[M::%s] total inter hic links: %lld\n", __func__, htot);
    fprintf(stderr, "[M::%s] per-window hic density: %.3e\n", __func__, h_pb);

    // filter short/repetitive sequences and with low hic density
    MYBZERO(smark, nseq);
    j = 0;
    for (i = 0; i < nseq; i++) {
        if (clens[i] < MIN_HIC_WINDOW ||
            hcnts[i] == 0)
            smark[i] = 1;
        else
            sords[j++] = (ord_dbl_t) {i, log2((double) hcnts[i] / clens[i])};
    }
    // MAD filtering
    if (j > 0) {
        qsort(sords, j, sizeof(ord_dbl_t), ord_dbl_acmpfunc);
        median = (j % 2) ? sords[j/2].event : 0.5 * (sords[j/2 - 1].event + sords[j/2].event);
        for (i = 0; i < j; i++)
            sords[i].event = fabs(sords[i].event - median);
        qsort(sords, j, sizeof(ord_dbl_t), ord_dbl_acmpfunc);
        median = (j % 2) ? sords[j/2].event : 0.5 * (sords[j/2 - 1].event + sords[j/2].event);
        for (i = 0; i < j; i++)
            if (sords[i].event * 0.6745 > MAX_MODIFIED_Z * median)
                smark[sords[i].which] = 1;
    }
    
    // pre-compute trans overlaps
    for (i = 0; i < novl; i++)
        tcnts[i] = (ovl_cnt_t) {
            i, 
            ovls[i].bread, 
            ovls[i].bbpos / HIC_NORM_WINDOW, 
            ovls[i].bepos / HIC_NORM_WINDOW, 
            ovls[i].del ? -1 : 0 // mark deleted overlaps here
        };
    // sort by b so we can do a linear scan to compute all trans
    qsort(tcnts, novl, sizeof(ovl_cnt_t), ovl_cnt_b_cmpfunc);
    b = tcnts[0].b;
    for (i = j = 0; i <= novl; i++) {
        if (i == novl || tcnts[i].b != b) {
            count_olaps(srngs + (ridx[b] >> 32), (uint32) ridx[b], tcnts + j, i - j);
            if (i == novl) break;
            j = i;
            b = tcnts[j].b;
        }
    }
    // sort it back so will share index with the original overlaps
    qsort(tcnts, novl, sizeof(ovl_cnt_t), ovl_cnt_w_cmpfunc);
    // count for each sequence
    for (i = 0; i < nseq; i++) {
        n = (uint32) index[i];
        w = dicts->s[i].len / HIC_NORM_WINDOW;
        oc = tcnts + (index[i]>>32);
        tlens[i] = count_trans_olaps(oc, n);
        tlens_beg[i] = tlens_end[i] = 0;
        for (j = 0; j < n; j++) {
            ov = ovls + oc[j].which;
            if (ov->del) continue;
            if (ov->abpos / HIC_NORM_WINDOW < Wend)
                tlens_beg[i] += oc[j].cnt;
            if (ov->aepos / HIC_NORM_WINDOW >= w - Wend)
                tlens_end[i] += oc[j].cnt;
        }
    }

    // calculate hic linkages
    kv_init(hlks);
    kv_resize(hlk_t, hlks, 1<<20);
    for (i = 0; i < nidx; i++) {
        idx = hidx + i;
        a = idx->a;
        b = idx->b;
        if (a >= b || smark[a] || smark[b])
            continue;
        if (find_b_overlap(ovls + (index[a] >> 32), (uint32) index[a], b))
            continue; // skip overlapped sequences

        // check four beg/end recombinations to find the one with the most cis links
        a_wl = dicts->s[a].len / HIC_NORM_WINDOW;
        b_wl = dicts->s[b].len / HIC_NORM_WINDOW;
        a_st = (a_wl <= Wend);
        b_st = (b_wl <= Wend);
        best_ea = best_eb = 0;
        cish_ab = -1;
        for (ea = 0; ea <= (a_st ? 0 : 1); ea++) {
            if (ea == 0) { 
                ea_wb = 0;
                ea_we = MIN(a_wl, Wend);
            } else {
                ea_wb = a_wl - Wend;
                ea_we = a_wl;
            }
            for (eb = 0; eb <= (b_st ? 0 : 1); eb++) {
                if (eb == 0) {
                    eb_wb = 0;
                    eb_we = MIN(b_wl, Wend);
                } else {
                    eb_wb = b_wl - Wend;
                    eb_we = b_wl;
                }
                n = count_cis_links(hics, idx->s, idx->e, ea_wb, ea_we - 1, eb_wb, eb_we - 1);
                if (n > cish_ab) {
                    cish_ab = n;
                    best_ea = ea;
                    best_eb = eb;
                }
            }
        }
        // set final end window ranges
        if (best_ea == 0) { 
            a_wb = 0;
            a_we = MIN(a_wl, Wend);
        } else {
            a_wb = a_wl - Wend;
            a_we = a_wl;
        }
        if (best_eb == 0) { 
            b_wb = 0;
            b_we = MIN(b_wl, Wend);
        } else {
            b_wb = b_wl - Wend;
            b_we = b_wl;
        }

        if (cish_ab < MIN_CIS_HLINK)
            continue;

        cisl_a   = best_eb == 0? clens_beg[b] : clens_end[b];
        cisl_b   = best_ea == 0? clens_beg[a] : clens_end[a];
        transl_a = best_eb == 0? tlens_beg[b] : tlens_end[b];
        transl_b = best_ea == 0? tlens_beg[a] : tlens_end[a];
        transh_a = count_trans_links(hics, hidx, nidx, a, a_wb, a_we - 1,
            ovls + (index[b]>>32), (uint32) index[b], b_wb, b_we);
        transh_b = count_trans_links(hics, hidx, nidx, b, b_wb, b_we - 1,
            ovls + (index[a]>>32), (uint32) index[a], a_wb, a_we);

        hc_a = (best_ea == 0) ? hcnts_beg[a] : hcnts_end[a];
        hc_b = (best_eb == 0) ? hcnts_beg[b] : hcnts_end[b];
        p1 = calculate_poisson_pvalue(cish_ab, (double) hc_a * hc_b / htot);

        p2 = (cisl_a + transl_a > 0)
             ? calculate_binomial_pvalue(cish_ab, cish_ab + transh_a, (double) cisl_a / (cisl_a + transl_a))
             : 1.0;
        p3 = (cisl_b + transl_b > 0)
             ? calculate_binomial_pvalue(cish_ab, cish_ab + transh_b, (double) cisl_b / (cisl_b + transl_b))
             : 1.0;
        
        p_val = p1;
        p_val = MAX(p_val, p2);
        p_val = MAX(p_val, p3);

        cis_l = 2.0 * cisl_a * cisl_b / (cisl_a + cisl_b);
        //cis_l = sqrt(cisl_a) * sqrt(cisl_b);
        d_val = cish_ab / sqrt(cisl_a) / sqrt(cisl_b);

        kv_push(hlk_t, hlks, ((hlk_t) {a, b, cis_l, d_val, p_val}));

#ifdef DEBUG_HIC_LINKAGE
        fprintf(stdout, "%s\t%d\t%lld\t%s\t%d\t%lld\t|\t%12.3f\t%12.3f|\t%lld\t%lld\t%lld\t%lld\t|\t%lld\t%lld\t%lld\t%lld\t|\t%.3e\t%.3e\t%.3e\t|\t%.3e\tPG:A:%c\n", 
            dicts->s[a].name, clens[a], hcnts[a], 
            dicts->s[b].name, clens[b], hcnts[b],
            cis_l, d_val,
            cisl_a, cish_ab, transl_a, transh_a, 
            cisl_b, cish_ab, transl_b, transh_b,
            p1, p2, p3, p_val,
            PG_HIC);
#endif
    }

    nhlk = hlks.n;
    if (!nhlk) {
        kv_destroy(hlks);
        hlks.a = NULL;
        goto cleanup;
    }

    // Benjamini-Hochberg adjustment
    qsort(hlks.a, nhlk, sizeof(hlk_t), hlk_pval_cmpfunc);
    p_thresh = PVAL_THRESHOLD / nhlk;
    for (i = nhlk - 1; i >= 0; i--)
        if (hlks.a[i].p <= (double) (i + 1) * p_thresh)
            break;
    p_thresh = hlks.a[i].p;
    nhlk = i + 1;

    fprintf(stderr, "[M::%s] minimum cis-link count: %d\n", __func__, MIN_CIS_HLINK);
    fprintf(stderr, "[M::%s] candidate hic linkages: %zu\n", __func__, hlks.n);
    fprintf(stderr, "[M::%s] maximum p-value: %.3e\n", __func__, PVAL_THRESHOLD);
    fprintf(stderr, "[M::%s] BH adjusted p-value: %.3e\n", __func__, p_thresh);
    fprintf(stderr, "[M::%s] final linkage retained: %lld\n", __func__, nhlk);

    // make symmetric copies
    kv_resize(hlk_t, hlks, nhlk * 2);
    hlk_t *src = hlks.a;
    hlk_t *dst = hlks.a + nhlk;
    for (i = 0; i < nhlk; i++, src++, dst++)
        *dst = (hlk_t) {src->b, src->a, src->l, src->v, src->p};
    nhlk *= 2;

    // sort by sequence ids
    qsort(hlks.a, nhlk, sizeof(hlk_t), hlk_abseq_cmpfunc);
    
    MYREALLOC(hlks.a, nhlk);
    if (_nhlk) *_nhlk = nhlk;

    // run tests
    //test_hic_linkage_map(hlks.a, nhlk, ovls, novl, dicts);

cleanup:
    free(hics);
    free(hidx);
    free(ridx);
    free(smark);
    free(clens);
    free(clens_beg);
    free(clens_end);
    free(tlens);
    free(tlens_beg);
    free(tlens_end);
    free(hcnts);
    free(hcnts_beg);
    free(hcnts_end);
    free(tcnts);
    free(wcnts);
    free(srngs);
    free(sords);
    free(index);

    return hlks.a;
}

static const double BACKBONE_MIN_RECRUIT = 0.8;
static const double BACKBONE_MIN_CONFIDENCE_SCORE = 0.7;

static const double BACKBONE_MIN_OVL_SIZE = 5.0; // log10(100000)
static const double BACKBONE_MAX_OVL_SIZE = 6.0; // log10(1000000)
static const double BACKBONE_LOG_DIFF_OVL_SIZE = BACKBONE_MAX_OVL_SIZE - BACKBONE_MIN_OVL_SIZE;
static const double BACKBONE_MIN_OVL_TO_EXT_RATIO_AT_MIN = 1.0;
static const double BACKBONE_MIN_OVL_TO_EXT_RATIO_AT_MAX = 0.1;
static const double BACKBONE_OVL_TO_EXT_RATIO_DIFF = BACKBONE_MIN_OVL_TO_EXT_RATIO_AT_MAX - BACKBONE_MIN_OVL_TO_EXT_RATIO_AT_MIN;
static inline double backbone_min_ovl_to_ext_ratio(double l)
{
    if (l <= 1.0) return BACKBONE_MIN_OVL_TO_EXT_RATIO_AT_MIN;

    double log_L = log10(l);
    if (log_L <= BACKBONE_MIN_OVL_SIZE) return BACKBONE_MIN_OVL_TO_EXT_RATIO_AT_MIN;
    if (log_L >= BACKBONE_MAX_OVL_SIZE) return BACKBONE_MIN_OVL_TO_EXT_RATIO_AT_MAX;
    return BACKBONE_MIN_OVL_TO_EXT_RATIO_AT_MIN + (log_L - BACKBONE_MIN_OVL_SIZE) / BACKBONE_LOG_DIFF_OVL_SIZE * BACKBONE_OVL_TO_EXT_RATIO_DIFF;
}

// number of conflict types
// 1. sequence overlap
// 2. hic links
#define CFT_OVL 0
#define CFT_HIC 1

typedef struct { int64 vals[2]; } cft_v_t;

typedef struct {
    int a, b; // sequences
    uint32 v:31, t:1; // value and type
} cft_t;

typedef kvec_t(cft_t) cft_vec_t;

static inline int cft_v_zero(cft_v_t *v)
{
    return v->vals[CFT_OVL] == 0 && v->vals[CFT_HIC] == 0;
}

static inline void cft_v_reset(cft_v_t *v)
{
    v->vals[CFT_OVL] = v->vals[CFT_HIC] = 0;
}

static inline void cft_v_add(cft_v_t *a, cft_v_t *b)
{
    a->vals[CFT_OVL] += b->vals[CFT_OVL];
    a->vals[CFT_HIC] += b->vals[CFT_HIC];
}

static inline void cft_v_sub(cft_v_t *a, cft_v_t *b)
{
    a->vals[CFT_OVL] -= b->vals[CFT_OVL];
    a->vals[CFT_HIC] -= b->vals[CFT_HIC];
}

static int cft_v_cmpfunc(cft_v_t *a, cft_v_t *b)
{
    // we prefer lower overlaps but higher hic links
    // the entry with larger value is considered better
    // i.e., fewer conflicts
    if (a->vals[CFT_OVL] != b->vals[CFT_OVL])
        return (a->vals[CFT_OVL] > b->vals[CFT_OVL]) - (a->vals[CFT_OVL] < b->vals[CFT_OVL]);
    return (a->vals[CFT_HIC] < b->vals[CFT_HIC]) - (a->vals[CFT_HIC] > b->vals[CFT_HIC]);
}

static int cft_abseqs_cmpfunc(const void *a, const void *b)
{
    cft_t *x = (cft_t *) a;
    cft_t *y = (cft_t *) b;

    if (x->a != y->a)
        return (x->a > y->a) - (x->a < y->a);
    if (x->b != y->b)
        return (x->b > y->b) - (x->b < y->b);
    return (x->t > y->t) - (x->t < y->t);
}

// calculate conflict after removing sequence from group h
// sequence is the first sequence in cfts
// the haplotype group is not actually updated here
// the conflict values are changed in-place in hcft
static inline void cft_v_unset(cft_t *cfts, int ncft, uint8 h, uint8 *haps, cft_v_t *hcft)
{
    if (!h || cfts == NULL || !ncft) return;
    int i;
    for (i = 0; i < ncft; i++, cfts++)
        if (haps[cfts->b] == h) hcft->vals[cfts->t] -= cfts->v;
}

// calculate conflict after adding sequence to group h
// sequence is the first sequence in cfts
// the haplotype group is not actually updated here
// the conflict values are changed in-place in hcft
static inline void cft_v_set(cft_t *cfts, int ncft, uint8 h, uint8 *haps, cft_v_t *hcft)
{
    if (!h || cfts == NULL || !ncft) return;
    int i;
    for (i = 0; i < ncft; i++, cfts++)
        if (haps[cfts->b] == h) hcft->vals[cfts->t] += cfts->v;
}

// calculate conflict after changing sequence from group hc to ht
// sequence is the first sequence in cfts
// the haplotype group is not actually updated here
// the conflict values are changed in-place in hcft
static inline void cft_v_move(cft_t *cfts, int ncft, uint8 hc, uint8 ht, uint8 *haps, cft_v_t *hcft)
{
    if (hc == ht || cfts == NULL || !ncft) return;
    if (!hc) return cft_v_set(cfts, ncft, ht, haps, hcft);
    if (!ht) return cft_v_unset(cfts, ncft, hc, haps, hcft);
    int i;
    for (i = 0; i < ncft; i++, cfts++) {
        if (haps[cfts->b] == hc)
            hcft->vals[cfts->t] -= cfts->v;
        else if (haps[cfts->b] == ht)
            hcft->vals[cfts->t] += cfts->v;
    }
}

// conflicts to all haplotype groups
// hctfs size should be at least (ploidy+1)
static inline void cft_v_all(cft_t *cfts, int ncft, uint8 *haps, cft_v_t *hcfts)
{
    if (cfts == NULL || !ncft) return;
    int i;
    for (i = 0; i < ncft; i++, cfts++)
        hcfts[haps[cfts->b]].vals[cfts->t] += cfts->v;
}

static inline void cft_v_total(cft_t *cfts, int ncft, uint8 *haps, cft_v_t *hcft)
{
    while (ncft--) {
        if (haps[cfts->a] && (haps[cfts->a] == haps[cfts->b]))
            hcft->vals[cfts->t] += cfts->v;
        cfts++;
    }
    // each overlap is counted twice
    hcft->vals[CFT_OVL] /= 2;
    hcft->vals[CFT_HIC] /= 2;
}

// this is for a min-heap
typedef struct {
    cft_v_t cfs;
    uint32 len;
    uint8 hap;
} cft_prior_t;

static void cft_prior_cpyfunc(void *x, void *y)
{
    *(cft_prior_t *) x = *(cft_prior_t *) y;
}

static int fm_prior_cmpfunc(void *x, void *y)
{
    cft_prior_t *a = (cft_prior_t *) x;
    cft_prior_t *b = (cft_prior_t *) y;

    // for selecting next sequence to place
    // we prioritise larger gain, longer sequence, and smaller haplotype
    int res;
    if ((res = cft_v_cmpfunc(&b->cfs, &a->cfs)) != 0)
        return res;
    if (a->len != b->len)
        return (a->len < b->len) - (a->len > b->len);
    if (a->hap != b->hap)
        return (a->hap > b->hap) - (a->hap < b->hap);
    return 0;
}

static int fs_prior_cmpfunc(void *x, void *y)
{
    cft_prior_t *a = (cft_prior_t *) x;
    cft_prior_t *b = (cft_prior_t *) y;

    // for selecting next sequence to solve conflict
    // we prioritise larger gain, smaller sequence and larger haplotype
    int res;
    if ((res = cft_v_cmpfunc(&b->cfs, &a->cfs)) != 0)
        return res;
    if (a->len != b->len)
        return (a->len > b->len) - (a->len < b->len);
    if (a->hap != b->hap)
        return (a->hap < b->hap) - (a->hap > b->hap);
    return 0;
}

#define PQ_INVALID_INDEX -1 

typedef struct {
    int id;
    void *priority;
} pq_elem_t;

static pq_elem_t PQ_ELEM_NULL = (pq_elem_t) {0, NULL};

typedef struct {
    int PQ_HEAP_SIZE;
    int *PQ_ELEM_MAP; // position of element [id] in heap
    pq_elem_t *PQ_HEAP;
    size_t PQ_PRIOR_WIDTH;
    void *PQ_PRIOR;
    void (*PQ_CPYFUNC)(void *, void *);
    int  (*PQ_CMPFUNC)(void *, void *);
} pq_t;

#define PQ_PRIOR_AT(pq, id) ((pq)->PQ_PRIOR + (id) * (pq)->PQ_PRIOR_WIDTH)

void pq_elem_swap(pq_t *pq, int i, int j)
{
    pq_elem_t temp = pq->PQ_HEAP[i];
    pq->PQ_HEAP[i] = pq->PQ_HEAP[j];
    pq->PQ_HEAP[j] = temp;
    pq->PQ_ELEM_MAP[pq->PQ_HEAP[i].id] = i;
    pq->PQ_ELEM_MAP[pq->PQ_HEAP[j].id] = j;
}

static void pq_heapify_up(pq_t *pq, int index)
{
    int parent;
    while (index > 0) {
        parent = (index - 1) / 2;
        if (pq->PQ_CMPFUNC(pq->PQ_HEAP[parent].priority, pq->PQ_HEAP[index].priority) > 0) {
            pq_elem_swap(pq, parent, index);
            index = parent;
        } else break;
    }
}

void pq_heapify_down(pq_t *pq, int index)
{
    int l, r, i;
    while (1) {
        l = 2 * index + 1;
        r = 2 * index + 2;
        i = index;
        if (l < pq->PQ_HEAP_SIZE && pq->PQ_CMPFUNC(pq->PQ_HEAP[l].priority, pq->PQ_HEAP[i].priority) < 0)
            i = l;
        if (r < pq->PQ_HEAP_SIZE && pq->PQ_CMPFUNC(pq->PQ_HEAP[r].priority, pq->PQ_HEAP[i].priority) < 0)
            i = r;
        if (i == index) break;
        pq_elem_swap(pq, index, i);
        index = i;
    }
}

static pq_elem_t pq_pop(pq_t *pq) 
{
    if (pq->PQ_HEAP_SIZE == 0)
        return PQ_ELEM_NULL;
    pq_elem_t top = pq->PQ_HEAP[0];
    pq->PQ_ELEM_MAP[top.id] = PQ_INVALID_INDEX;
    pq->PQ_HEAP_SIZE--;
    if (pq->PQ_HEAP_SIZE > 0) {
        pq->PQ_HEAP[0] = pq->PQ_HEAP[pq->PQ_HEAP_SIZE];
        pq->PQ_ELEM_MAP[pq->PQ_HEAP[0].id] = 0;
        pq_heapify_down(pq, 0);
    }
    return top;
}

static int pq_pop1(pq_t *pq, pq_elem_t *top) 
{
    if (pq->PQ_HEAP_SIZE == 0)
        return 0;
    *top = pq->PQ_HEAP[0];
    pq->PQ_ELEM_MAP[top->id] = PQ_INVALID_INDEX;
    pq->PQ_HEAP_SIZE--;
    if (pq->PQ_HEAP_SIZE > 0) {
        pq->PQ_HEAP[0] = pq->PQ_HEAP[pq->PQ_HEAP_SIZE];
        pq->PQ_ELEM_MAP[pq->PQ_HEAP[0].id] = 0;
        pq_heapify_down(pq, 0);
    }
    return 1;
}

static pq_elem_t pq_peek(pq_t *pq) 
{
    if (pq->PQ_HEAP_SIZE == 0)
        return PQ_ELEM_NULL;
    return pq->PQ_HEAP[0];
}

static pq_elem_t *pq_at(pq_t *pq, int index) 
{
    return pq->PQ_HEAP + index;
}

static int pq_size(pq_t *pq) 
{
    return pq->PQ_HEAP_SIZE;
}

static int pq_is_empty(pq_t *pq) 
{
    return pq->PQ_HEAP_SIZE == 0;
}

static void *pq_prior(pq_t *pq, int id, int check) 
{
    int index = pq->PQ_ELEM_MAP[id];
    if (check && index == PQ_INVALID_INDEX)
        return NULL;
    else return PQ_PRIOR_AT(pq, id);
}

static int pq_exist(pq_t *pq, int id)
{
    return pq->PQ_ELEM_MAP[id] != PQ_INVALID_INDEX;
}

static pq_t *pq_init(int size, size_t width, void (*cpyfunc)(void *, void *), int (*cmpfunc)(void *, void *))
{
    pq_t *pq;
    int i;
    MYCALLOC(pq, 1);
    if (!pq)
        mem_alloc_error("pq_t struct");
    MYCALLOC(pq->PQ_HEAP, size);
    MYCALLOC(pq->PQ_ELEM_MAP, size);
    pq->PQ_PRIOR = malloc(width * size);
    if (pq->PQ_HEAP == NULL || pq->PQ_PRIOR == NULL || 
        pq->PQ_ELEM_MAP == NULL)
        mem_alloc_error("pq arrays");
    for (i = 0; i < size; i++)
        pq->PQ_ELEM_MAP[i] = PQ_INVALID_INDEX;
    pq->PQ_PRIOR_WIDTH = width;
    pq->PQ_CPYFUNC = cpyfunc;
    pq->PQ_CMPFUNC = cmpfunc;
    pq->PQ_HEAP_SIZE = 0;
    return pq;
}

static void pq_destroy(pq_t *pq)
{
    free(pq->PQ_HEAP);
    free(pq->PQ_PRIOR);
    free(pq->PQ_ELEM_MAP);
    free(pq);
}

static void pq_insert(pq_t *pq, int id, void *priority)
{
    int index = pq->PQ_ELEM_MAP[id];
    void *pptr = PQ_PRIOR_AT(pq, id);
    if (index == PQ_INVALID_INDEX) {
        // insert
        pq->PQ_CPYFUNC(pptr, priority); // copy data
        pq->PQ_HEAP[pq->PQ_HEAP_SIZE].id = id;
        pq->PQ_HEAP[pq->PQ_HEAP_SIZE].priority = pptr;
        pq->PQ_ELEM_MAP[id] = pq->PQ_HEAP_SIZE;
        pq_heapify_up(pq, pq->PQ_HEAP_SIZE++);
    } else {
        // update
        if (pq->PQ_CMPFUNC(priority, pptr) < 0) {
            pq->PQ_CPYFUNC(pptr, priority); // copy data
            pq_heapify_up(pq, index);
        } else {
            pq->PQ_CPYFUNC(pptr, priority); // copy data
            pq_heapify_down(pq, index);
        }
    }
}

static void pq_remove(pq_t *pq, int id)
{
    int index = pq->PQ_ELEM_MAP[id];
    if (index != PQ_INVALID_INDEX) {
        pq->PQ_ELEM_MAP[id] = PQ_INVALID_INDEX;
        pq->PQ_HEAP_SIZE--;
        if (index < pq->PQ_HEAP_SIZE) {
            pq->PQ_HEAP[index] = pq->PQ_HEAP[pq->PQ_HEAP_SIZE];
            pq->PQ_ELEM_MAP[pq->PQ_HEAP[index].id] = index;
            pq_heapify_up(pq, index);
            pq_heapify_down(pq, index);
        }
    }
}

typedef struct {
    int seq;
    uint8 hap;
} hs_opt_t;

typedef kvec_t(hs_opt_t) hs_opts_t;

typedef struct {
    int64 max_iter;
    uint8 *haps;
    uint8 *bhap;
    int *smark;
    cft_v_t *hcfts;
    cft_v_t *dcfts;
    pq_t *fm_q;
    pq_t *fs_q;
    hs_opts_t hs_opts;
} vns_data_t;

static vns_data_t *vns_data_init(int size, int ploidy)
{
    vns_data_t *data;
    MYMALLOC(data, 1);
    if (!data)
        mem_alloc_error("vns_data_t struct");
    MYMALLOC(data->haps, size);
    MYMALLOC(data->bhap, size);
    MYMALLOC(data->smark, size);
    MYMALLOC(data->hcfts, ploidy+3);
    MYMALLOC(data->dcfts, size*2);
    data->fm_q = pq_init(size, sizeof(cft_prior_t), cft_prior_cpyfunc, fm_prior_cmpfunc);
    data->fs_q = pq_init(size, sizeof(cft_prior_t), cft_prior_cpyfunc, fs_prior_cmpfunc);
    kv_init(data->hs_opts);
    if (data->haps == NULL || data->bhap == NULL || 
        data->smark == NULL || data->hcfts == NULL ||
        data->fm_q == NULL || data->fs_q == NULL ||
        data->dcfts == NULL)
        mem_alloc_error("vns_data_t arrays");
    data->max_iter = 1000000;
    return data;
}

static void vns_data_destroy(vns_data_t *data)
{
    if (!data) return;
    free(data->haps);
    free(data->bhap);
    free(data->smark);
    free(data->hcfts);
    free(data->dcfts);
    pq_destroy(data->fm_q);
    pq_destroy(data->fs_q);
    kv_destroy(data->hs_opts);
    free(data);
}

typedef struct {
    int nseq;
    uint32 *seqs;
    uint32 *type;
    int64 slen;
} scf_block_t;

typedef struct {
    uint32 seq:31, rev:1;
    uint32 len;
    int grp;
    uint8 hap;
    int64 bpos, epos; // position on the backbone
} hap_info_t;

static int hap_info_acmpfunc(const void *a, const void *b)
{
    const hap_info_t *x = a;
    const hap_info_t *y = b;
    if (x->grp != y->grp)
        return (x->grp > y->grp) - (x->grp < y->grp);
    if (x->hap != y->hap)
        return (x->hap > y->hap) - (x->hap < y->hap);
    if (x->bpos != y->bpos)
        return (x->bpos > y->bpos) - (x->bpos < y->bpos);
    if (x->epos != y->epos)
        return (x->epos > y->epos) - (x->epos < y->epos);
    if (x->len != y->len)
        return (x->len < y->len) - (x->len > y->len);
    return (x->seq > y->seq) - (x->seq < y->seq);
}

static int hap_info_gcmpfunc(const void *a, const void *b)
{
    const hap_info_t *x = a;
    const hap_info_t *y = b;
    if (x->grp != y->grp)
        return (x->grp > y->grp) - (x->grp < y->grp);
    if (x->bpos != y->bpos)
        return (x->bpos > y->bpos) - (x->bpos < y->bpos);
    if (x->epos != y->epos)
        return (x->epos > y->epos) - (x->epos < y->epos);
    if (x->len != y->len)
        return (x->len < y->len) - (x->len > y->len);
    if (x->hap != y->hap)
        return (x->hap > y->hap) - (x->hap < y->hap);
    return (x->seq > y->seq) - (x->seq < y->seq);
}

typedef struct {
    uint32 s;
    ovl_t *o;
    int p;
    int g;
    int ovl, ext;
} ext_info_t;

typedef struct {
    int g;
    double score;
} ext_clus_t;

static int einfo_ext_dcmpfunc(const void *a, const void *b)
{
    const ext_info_t *x = a;
    const ext_info_t *y = b;
    if (x->ext != y->ext)
        return (x->ext < y->ext) - (x->ext > y->ext);
    if (x->ovl != y->ovl)
        return (x->ovl < y->ovl) - (x->ovl > y->ovl);
    return (x->s > y->s) - (x->s < y->s);
}

static int eclus_score_dcmpfunc(const void *a, const void *b)
{
    const ext_clus_t *x = a;
    const ext_clus_t *y = b;
    return (x->score < y->score) - (x->score > y->score);
}

static inline double scaled_b_ovl_size(ovl_t *ovl)
{
    return (double) ovl->blen * log2(ovl->neff + 1.);
}

typedef struct {
        int a, b;
        ovl_t *o;
} grp_ovl_t;

typedef struct {
    uint32 g:31, r:1;
    double w;   // total (symmetric) weight — used for spectral and degree
    double wp;  // directed forward weight: cost of placing g BEFORE this node
                // i.e., wp[i→j] > 0 means evidence says i should come before j;
                // cost contribution = wp * max(0, rank_j - rank_i)  [correct order: free]
                //                   + (w-wp) * max(0, rank_i - rank_j)  [wrong order: penalised]
} grp_edge_t;

typedef kvec_t(grp_edge_t) grp_edge_vec_t;

typedef struct {
    int a, b;
    double cw, pw;
} grp_conflict_t;

typedef struct {
    int seq;
    int64 beg;
    int64 end;
} bb_range_t;

static int bb_range_pos_acmpfunc(const void *a, const void *b)
{
    const bb_range_t *x = a;
    const bb_range_t *y = b;
    if (x->beg != y->beg)
        return (x->beg > y->beg) - (x->beg < y->beg);
    return (x->end > y->end) - (x->end < y->end);
}

static int govl_abseq_cmpfunc(const void *a, const void *b)
{
    grp_ovl_t *x = (grp_ovl_t *) a;
    grp_ovl_t *y = (grp_ovl_t *) b;
    if (x->a != y->a)
        return (x->a > y->a) - (x->a < y->a);
    if (x->b != y->b)
        return (x->b > y->b) - (x->b < y->b);
    return 0;
}

#define SWAP_WIN_MAX 1000
#define SWAP_WIN_MIN  100
#define FM_PASS_MAX  1000

// priority for group-ordering FM swap PQ: keyed by group ga,
// stores its best swap partner (gb) and the resulting cost delta.
// min-heap on delta (most negative = most improving = highest priority).
typedef struct { double delta; int partner; } swap_prior_t;
static void swap_prior_cpyfunc(void *x, void *y) { *(swap_prior_t*)x = *(swap_prior_t*)y; }
static int swap_prior_cmpfunc(void *x, void *y) 
{
    double da = ((swap_prior_t*)x)->delta, db = ((swap_prior_t*)y)->delta;
    return (da > db) - (da < db);
}

// Directed edge cost: for edge (i→j) with forward weight wp and total weight w,
// cost = wp * max(0, rj - ri) + (w - wp) * max(0, ri - rj)
// = correct order is cheap (weighted by wp); wrong order is expensive (weighted by w-wp).
// When wp = w/2 (no directional signal) this reduces to w/2 * |ri - rj| = symmetric cost.
static inline double directed_cost(double w, double wp, int ri, int rj)
{
    int d = rj - ri;
    return d >= 0 ? wp * d : (w - wp) * (-d);
}

// Compute best window-swap partner for group g given current ranks.
// locked may be NULL (then all partners are candidates).
// Returns partner id, writes delta to *out_delta; returns -1 if none found.
static int best_swap_partner(int g, int *grp_rank, int *grp_order, int ngrp,
    grp_edge_vec_t *sadj, const int *locked, int win, double *out_delta)
{
    int ra, lo, hi, rb, gb, rkm, e, km;
    int best_gb;
    double w, wp, delta, best_delta;

    ra = grp_rank[g];
    lo = ra - win; if (lo < 0) lo = 0;
    hi = ra + win + 1; if (hi > ngrp) hi = ngrp;
    best_gb    = -1;
    best_delta = 1e300;
    for (rb = lo; rb < hi; rb++) {
        if (rb == ra) continue;
        gb = grp_order[rb];
        if (locked && locked[gb]) continue;
        delta = 0.;
        for (e = 0; e < (int)sadj[g].n; e++) {
            km = sadj[g].a[e].g; if (km == gb) continue;
            rkm = grp_rank[km];
            w = sadj[g].a[e].w; wp = sadj[g].a[e].wp;
            delta += directed_cost(w, wp, rb, rkm) - directed_cost(w, wp, ra, rkm);
        }
        for (e = 0; e < (int)sadj[gb].n; e++) {
            km = sadj[gb].a[e].g; if (km == g) continue;
            rkm = grp_rank[km];
            w = sadj[gb].a[e].w; wp = sadj[gb].a[e].wp;
            delta -= directed_cost(w, wp, rb, rkm) - directed_cost(w, wp, ra, rkm);
        }
        // For directed costs, swapping g↔gb also changes the cost of the direct {g,gb} edge:
        // old: directed_cost(w, wp, ra, rb); new: directed_cost(w, wp, rb, ra).
        // (Symmetric costs cancel, but directed costs do not.)
        for (e = 0; e < (int)sadj[g].n; e++) {
            if ((int)sadj[g].a[e].g == gb) {
                w = sadj[g].a[e].w; wp = sadj[g].a[e].wp;
                delta += directed_cost(w, wp, rb, ra) - directed_cost(w, wp, ra, rb);
                break;
            }
        }
        if (delta < best_delta) { best_delta = delta; best_gb = gb; }
    }
    *out_delta = best_delta;
    return best_gb;
}

// priority for local within-haplotype FM: keyed by sequence id,
// stores its best target rank and the crossing-count gain (positive = improvement).
// max-heap on gain.
typedef struct { double gain; int best_rank; int hap; } lmove_prior_t;
static void lmove_prior_cpyfunc(void *x, void *y) { *(lmove_prior_t*)x = *(lmove_prior_t*)y; }
static int lmove_prior_cmpfunc(void *x, void *y) 
{
    double da = ((lmove_prior_t*)x)->gain, db = ((lmove_prior_t*)y)->gain;
    return (da > db) - (da < db);   // max-heap: larger gain = higher priority
}

// Compute the signed crossing weight between sequence s (rank p in ha[])
// and sequence ha[p2] (rank p2) across all partner haplotypes.
// W > 0 means moving s past p2 reduces crossings; W < 0 means it increases crossings.
// Partners of s on hj: loop over all hj sequences, check overlap.
// Partners of ha[p2] on hj: same.
// Weight = min(alen, blen) signed by relative rank order on hj.
static double lmove_crossing_weight(int s_seq, int p2_seq,
    hap_info_t **hseqs, int *hcnt, int ploidy_hi, int hi,
    ovl_t *ovls, uint64 *index)
{
    hap_info_t *hb;
    ovl_t *oa, *ob;
    double w, wa, wb, wt;
    int hj, q, r, nb, sgn;

    w = 0.;
    for (hj = 0; hj < ploidy_hi; hj++) {
        if (hj == hi) continue;
        nb = hcnt[hj];
        hb = hseqs[hj];
        // for each pair (qa, qb) where qa overlaps s and qb overlaps p2_seq
        // crossing between (s,qa) and (p2,qb) is removed when rank_hj(qa) > rank_hj(qb)
        // and created when rank_hj(qa) < rank_hj(qb) — if s moves to right of p2
        for (q = 0; q < nb; q++) {
            oa = find_b_overlap(ovls + (index[s_seq] >> 32),
                                (uint32)index[s_seq], hb[q].seq);
            if (!oa || oa->del) continue;
            wa = (double)(oa->alen < oa->blen ? oa->alen : oa->blen);
            for (r = 0; r < nb; r++) {
                ob = find_b_overlap(ovls + (index[p2_seq] >> 32),
                                    (uint32)index[p2_seq], hb[r].seq);
                if (!ob || ob->del) continue;
                wb  = (double)(ob->alen < ob->blen ? ob->alen : ob->blen);
                wt  = wa < wb ? wa : wb;
                sgn = (q > r) - (q < r);
                w  += wt * sgn;
            }
        }
    }
    return w;
}

// Compute the best target rank within [max(0,p-win), min(n-1,p+win)] for sequence
// s_seq at current rank p in haplotype hi.  Returns gain (positive = improvement)
// and writes target rank to *best_r.
static double lmove_best_rank(int s_seq, int p, int hi,
    hap_info_t **hseqs, int *hcnt, int ploidy,
    ovl_t *ovls, uint64 *index, int win, int *best_r, const int *locked)
{
    hap_info_t *ha;
    int n, lo, hi2, r, p2;
    double best_gain, cum, cw;

    ha  = hseqs[hi];
    n   = hcnt[hi];
    lo  = p - win; if (lo < 0) lo = 0;
    hi2 = p + win; if (hi2 >= n) hi2 = n - 1;
    best_gain = 0.;
    *best_r   = p;
    cum = 0.;
    // scan right: p+1 .. hi2
    for (r = p + 1; r <= hi2; r++) {
        p2 = ha[r].seq;
        if (locked && locked[p2]) break;   // stop at locked sequence
        cw  = lmove_crossing_weight(s_seq, p2, hseqs, hcnt, ploidy, hi, ovls, index);
        cum += cw;   // moving right past p2: gain = +W
        if (cum > best_gain) { best_gain = cum; *best_r = r; }
    }
    cum = 0.;
    // scan left: p-1 .. lo
    for (r = p - 1; r >= lo; r--) {
        p2 = ha[r].seq;
        if (locked && locked[p2]) break;
        cw  = lmove_crossing_weight(s_seq, p2, hseqs, hcnt, ploidy, hi, ovls, index);
        cum -= cw;   // moving left past p2: gain = -W
        if (cum > best_gain) { best_gain = cum; *best_r = r; }
    }
    return best_gain;
}

typedef struct { int p; int q; int64 ab, ae, bb, be; ovl_t *ovl; } hap_ovl_pair_t;
static int hap_ovl_pair_cmpfunc(const void *a, const void *b)
{
    const hap_ovl_pair_t *x = a, *y = b;
    if (x->ae != y->ae) return x->ae < y->ae ? -1 : 1;
    if (x->be != y->be) return x->be < y->be ? -1 : 1;
    return 0;
}

typedef struct olink {
    ovl_t *ovl;
    struct olink *next;
} olink_t;

typedef struct gadj {
    struct gadj *next;
    uint32 a, b;
    int ast, bst;
    olink_t *olinks, *olast;
    int *bcov;
    double mscore;
} gadj_t;

static int gadj_neighbor_cmpfunc(const void *a, const void *b)
{
    gadj_t *x = *(gadj_t **)a, *y = *(gadj_t **)b;
    if (x->b != y->b)
        return (x->b > y->b) - (x->b < y->b);
    return (x->a > y->a) - (x->a < y->a);
}

#define GADJ_CHUNK_SIZE 10000
typedef struct gadj_chunk {
    gadj_t adjs[GADJ_CHUNK_SIZE];
    struct gadj_chunk *next;
} gadj_chunk_t;

typedef struct {
    gadj_chunk_t *chunk;
    int idx;
} gadj_pool_t;

static gadj_t *gadj_new(gadj_pool_t *pool)
{
    if (pool->chunk == NULL || pool->idx >= GADJ_CHUNK_SIZE) {
        gadj_chunk_t *new_chunk;
        MYCALLOC(new_chunk, 1);
        new_chunk->next = pool->chunk;
        pool->chunk = new_chunk;
        pool->idx = 0;
    }
    return pool->chunk->adjs + (pool->idx++);
}

static gadj_t *gadj_copy(gadj_pool_t *pool, gadj_t *src)
{
    gadj_t *dst = gadj_new(pool);
    *dst = *src;
    dst->a = src->b;
    dst->b = src->a;
    dst->ast = src->bst;
    dst->bst = src->ast;
    return dst;
}

static void gadj_pool_destroy(gadj_pool_t *pool)
{
    if (!pool) return;

    gadj_t *adj;
    gadj_chunk_t *chunk = pool->chunk;
    int size = pool->idx;
    int i;
    while (chunk) {
        gadj_chunk_t *next = chunk->next;
        adj = chunk->adjs;
        for (i = 0; i < size; i++, adj++)
            if (adj->a < adj->b)
                free(adj->bcov);
        free(chunk);
        chunk = next;
        size = GADJ_CHUNK_SIZE;
    }
    free(pool);
}

typedef gadj_t* hnode_t;

typedef struct {
    // heap for merging groups
    hnode_t *heap;
    int heap_size;
    int heap_capacity;
} HeapSpace;

static HeapSpace *hs_init(int groups, int heap_capacity)
{
    HeapSpace *hs;
    MYCALLOC(hs, 1);
    hs->heap_capacity = heap_capacity;
    hs->heap_size = 0;
    MYMALLOC(hs->heap, heap_capacity);
    return hs;
}

static void hs_heap_push(HeapSpace *hs, hnode_t node)
{
    // ensure capacity
    if (hs->heap_size >= hs->heap_capacity) {
        hs->heap_capacity <<= 1;
        MYREALLOC(hs->heap, hs->heap_capacity);
    }
    
    int i = hs->heap_size++;
    hs->heap[i] = node;
    
    while (i != 0 && hs->heap[(i - 1) / 2]->mscore < hs->heap[i]->mscore) {
        hnode_t temp = hs->heap[i];
        hs->heap[i] = hs->heap[(i - 1) / 2];
        hs->heap[(i - 1) / 2] = temp;
        i = (i - 1) / 2;
    }
}

static hnode_t hs_heap_pop(HeapSpace *hs)
{
    if (hs->heap_size <= 0) return NULL;
    
    if (hs->heap_size == 1) {
        hs->heap_size--;
        return hs->heap[0];
    }
    
    hnode_t root = hs->heap[0];
    hs->heap[0] = hs->heap[hs->heap_size - 1];
    hs->heap_size--;
    
    int i = 0, left, right, best;
    while (1) {
        left = 2 * i + 1;
        right = 2 * i + 2;
        best = i;
        if (left < hs->heap_size && hs->heap[left]->mscore > hs->heap[best]->mscore)
            best = left;
        if (right < hs->heap_size && hs->heap[right]->mscore > hs->heap[best]->mscore)
            best = right;
        if (best != i) {
            hnode_t temp = hs->heap[i];
            hs->heap[i] = hs->heap[best];
            hs->heap[best] = temp;
            i = best;
        } else break;
    }
    return root;
}

void hs_destroy(HeapSpace *hs)
{
    free(hs->heap);
    free(hs);
}

static int merge_coverage(cov_point_t *acovs, int nacov, cov_point_t *bcovs, int nbcov, cov_point_t *covs)
{
    int ai, bi, apos, bpos, pos, acov, bcov, cov;
    int ncov, has_prev, prev_cov;

    acov = bcov = 0;
    ai = bi = 0;
    ncov = 0;
    has_prev = 0;
    while (ai < nacov || bi < nbcov) {
        apos = ai < nacov? acovs[ai].pos : INT32_MAX;
        bpos = bi < nbcov? bcovs[bi].pos : INT32_MAX;
        pos  = apos <= bpos? apos : bpos;

        if (apos == pos) acov = acovs[ai++].cov;
        if (bpos == pos) bcov = bcovs[bi++].cov;

        cov = acov + bcov;
        if (!has_prev || cov != prev_cov) {
            covs[ncov].pos = pos;
            covs[ncov].cov = cov;
            prev_cov = cov;
            has_prev = 1;
            ncov++;
        }
    }
    memcpy(acovs, covs, ncov * sizeof(cov_point_t));

    return ncov;
}

static cov_point_t *merge_offset_covs(cov_point_t *acovs, int nacov, cov_point_t *bcovs, int nbcov, 
    int boff, int brev, int *_nmcov)
{
    cov_point_t *mcovs;
    MYMALLOC(mcovs, (nacov + nbcov));
    if (!mcovs) {
        mem_alloc_error("merged coverage points");
        *_nmcov = 0;
        return NULL;
    }

    int ai, bi, apos, bpos, pos, acov, bcov, nmcov;
    acov = bcov = 0;
    nmcov = 0;
    if (!brev) {
        ai = bi = 0;
        while (ai < nacov || bi < nbcov) {
            apos = ai < nacov? acovs[ai].pos        : INT32_MAX;
            bpos = bi < nbcov? bcovs[bi].pos + boff : INT32_MAX;
            pos  = apos <= bpos? apos : bpos;

            if (apos == pos) acov = acovs[ai++].cov;
            if (bpos == pos) bcov = bcovs[bi++].cov;

            mcovs[nmcov].pos = pos;
            mcovs[nmcov].cov = acov + bcov;
            nmcov++;
        }
    } else {
        ai = 0, bi = nbcov - 1;
        while (ai < nacov || bi >= 0) {
            apos = ai < nacov? acovs[ai].pos     : INT32_MAX;
            bpos = bi >= 0? boff - bcovs[bi].pos : INT32_MAX;
            pos  = apos <= bpos? apos : bpos;

            if (apos == pos) acov = acovs[ai++].cov;
            if (bpos == pos) {
                bcov = (bi > 0) ? bcovs[bi - 1].cov : 0;
                bi--;
            }

            mcovs[nmcov].pos = pos;
            mcovs[nmcov].cov = acov + bcov;
            nmcov++;
        }
    }

    int i, cov, nout;
    i = nout = 0;
    while (i < nmcov) {
        mcovs[nout++] = mcovs[i];
        cov = mcovs[i].cov;
        while (++i < nmcov && mcovs[i].cov == cov);
    }

    *_nmcov = nout;
    return mcovs;
}

static void profile_coverage(cov_point_t *covs, int ncov, int ploidy, int *bcov, int reset)
{
    int i, c, p;
    if (reset)
        MYBZERO(bcov, ploidy+1);
    if (ncov <= 0)
        return;

    p = covs[0].pos;
    c = covs[0].cov;
    for (i = 1; i < ncov; i++) {
        if (c >= ploidy)
            c = ploidy;
        bcov[c] += covs[i].pos - p;
        p = covs[i].pos;
        c = covs[i].cov;
    }
}

static inline int cov_bucket(int cov, int ploidy)
{
    return cov >= ploidy ? ploidy : cov;
}

static void delta_coverage(cov_point_t *covs, int ncov, cov_point_t *tcovs, int tncov, int ploidy, int *bdiv, int reset)
{
    if (reset)
        MYBZERO(bdiv, ploidy+1);
    if (ncov <= 0 || tncov <= 0)
        return;
    
    int ai, bi, apos, bpos, pos, len, acov, bcov, prev_pos;
    ai = bi = 0;
    acov = bcov = 0;
    prev_pos = 0;
    while (ai < ncov || bi < tncov) {
        apos = (ai < ncov)  ? covs[ai].pos  : INT32_MAX;
        bpos = (bi < tncov) ? tcovs[bi].pos : INT32_MAX;
        pos  = apos < bpos ? apos : bpos;
        if (bcov != 0) {
            len = pos - prev_pos;
            bdiv[cov_bucket(acov,        ploidy)] -= len;
            bdiv[cov_bucket(acov + bcov, ploidy)] += len;
        }
        if (apos == pos) acov = covs[ai++].cov;
        if (bpos == pos) bcov = tcovs[bi++].cov;
        prev_pos = pos;
    }
}

static const double LOG_MIN_DELTA_SIZE = 4.0; // log10(10000)
static const double LOG_MAX_DELTA_SIZE = 6.0; // log10(1000000)
static const double LOG_DIFF_DELTA_SIZE = LOG_MAX_DELTA_SIZE - LOG_MIN_DELTA_SIZE;
static const double MIN_FRAC_AT_MIN = .50;
static const double MIN_FRAC_AT_MAX = .10;
static const double MIN_FRAC_DIFF = MIN_FRAC_AT_MAX - MIN_FRAC_AT_MIN;
static inline double min_delta_frac(double l)
{
    if (l <= 1.0) return MIN_FRAC_AT_MIN;

    double log_L = log10(l);
    if (log_L <= LOG_MIN_DELTA_SIZE) return MIN_FRAC_AT_MIN;
    if (log_L >= LOG_MAX_DELTA_SIZE) return MIN_FRAC_AT_MAX;
    return MIN_FRAC_AT_MIN + (log_L - LOG_MIN_DELTA_SIZE) / LOG_DIFF_DELTA_SIZE * MIN_FRAC_DIFF;
}

static int mergeable_group_pair(gadj_t *gadj, cov_point_t **covs, int *ncov, int *abcov, int *bbcov,
    hap_info_t *haps, int ploidy, int *tbcov, cov_point_t *tcovs, int *tcpts, srange_vec_t *srngv)
{
    olink_t *olinks;
    ovl_t *ovl;
    int i, j, a, tncov;
    double alen, blen, mlen, tsum, dfrac, score;

    // collect and sort overlap ranges
    olinks = gadj->olinks;
    srngv->n = 0;
    while (olinks) {
        ovl = olinks->ovl;
        kv_push(srange_t, *srngv, ((srange_t){ovl->aread, ovl->abpos, ovl->aepos}));
        kv_push(srange_t, *srngv, ((srange_t){ovl->bread, ovl->bbpos, ovl->bepos}));
        olinks = olinks->next;
    }
    qsort(srngv->a, srngv->n, sizeof(srange_t), srange_cmpfunc);

    // compute delta coverage after adding overlaps
    MYBZERO(tbcov, ploidy+1);
    a = srngv->a[0].seq;
    for (i = 1, j = 0; i <= srngv->n; i++) {
        if (i == srngv->n || srngv->a[i].seq != a) {
            tncov = calc_coverage_from_intervals(tcovs, tcpts, haps[a].len, srngv->a+j, i-j, NULL, pts_from_sranges);
            delta_coverage(covs[a], ncov[a], tcovs, tncov, ploidy, tbcov, 0);
            if (i < srngv->n)
                a = srngv->a[i].seq;
            j = i;
        }
    }

    alen = blen = tsum = .0;
    for (i = 0; i < ploidy; i++) {
        alen += abcov[i] / (1. + i);
        blen += bbcov[i] / (1. + i);
        tsum += tbcov[i] / (1. + i);
    }
    alen += abcov[ploidy];
    blen += bbcov[ploidy];
    tsum += tbcov[ploidy];

    mlen = MIN(alen, blen);
    tsum *= -1.0;

    // compute merging score
    dfrac = tsum / mlen;
    score = tsum * dfrac;

#ifdef DEBUG_SCAFFOLD_GROUP_MERGE
    // some debug info
    fprintf(stderr, "Group pair: %6d %10.0f %6d %10.0f @ %10.0f %6.3f %10.0f %c\n", gadj->a+1, alen, gadj->b+1, blen, tsum, dfrac, score, dfrac < min_delta_frac(mlen) ? 'N' : 'Y');
    fprintf(stderr, "A Cov Prof:"); for (i = 0; i <= ploidy; i++) fprintf(stderr, " %10d", abcov[i]); fprintf(stderr, "\n");
    fprintf(stderr, "B Cov Prof:"); for (i = 0; i <= ploidy; i++) fprintf(stderr, " %10d", bbcov[i]); fprintf(stderr, "\n");
    fprintf(stderr, "D Cov Prof:"); for (i = 0; i <= ploidy; i++) fprintf(stderr, " %10d", tbcov[i]); fprintf(stderr, "\n");
#endif

    if (dfrac < min_delta_frac(mlen))
        return 0;

    // copy outputs
    MYMALLOC(gadj->bcov, ploidy+1);
    memcpy(gadj->bcov, tbcov, (ploidy+1) * sizeof(int));
    gadj->mscore = score;

    return 1;
}



typedef struct { int a, b; int w; } mst_edge_t;

static int mstedge_w_dcmpfunc(const void *a, const void *b)
{
    const mst_edge_t *x = (const mst_edge_t *) a;
    const mst_edge_t *y = (const mst_edge_t *) b;
    return (x->w < y->w) - (x->w > y->w);
}

typedef struct { int a, b; int n; int64 l; double w; } grp_pair_t;
static int grppair_w_dcmpfunc(const void *a, const void *b)
{
    const grp_pair_t *x = (const grp_pair_t *) a;
    const grp_pair_t *y = (const grp_pair_t *) b;
    return (x->w < y->w) - (x->w > y->w);
}

static int grppair_g_cmpfunc(const void *a, const void *b)
{
    const grp_pair_t *x = (const grp_pair_t *) a;
    const grp_pair_t *y = (const grp_pair_t *) b;
    if (x->a != y->a)
        return (x->a > y->a) - (x->a < y->a);
    else
        return (x->b > y->b) - (x->b < y->b);
}

// MST-based linear ordering of sequences within each scaffold group
// For each group (blks[i]), builds a max spanning tree from intra-group overlaps
// then BFS from a leaf to assign haps[seq].bpos / haps[seq].epos
// skip[i] != 0 (or skip == NULL) controls whether group i is processed
static void mst_order_groups(scf_block_t *blks, int nblk, ovl_t *ovls, uint64 *index, hap_info_t *haps, 
    sdict_t *dicts, int nseq, const int *skip)
{
    kvec_t(mst_edge_t) me;
    ovl_t *ov, *oc, **mst_par_ovl;
    bb_range_t *mst_brngs;
    scf_block_t *blk;
    uint32 *bseqs, bb_r, bb_s;
    int *smap, *uf_r, *mst_par, *mst_deg, *mst_queue;
    int64 bmax, max_epos, me_n_agg, slen;
    int i, j, k, n, a, b, p;
    int bnseq, mst_start, n_mst, nq, root;
    int bb_alen, bb_blen, bb_bgap, bb_egap, bb_bshift, bb_eshift;
    int ia, ib, ra, rb, nb, cur, par, tmp;

    bmax = 0;
    for (i = 0; i < nblk; i++)
        if (blks[i].nseq > bmax)
            bmax = blks[i].nseq;
    if (bmax == 0) return;

    kv_init(me);
    MYMALLOC(uf_r,        bmax);
    MYMALLOC(mst_par,     bmax);
    MYMALLOC(mst_par_ovl, bmax);
    MYMALLOC(mst_deg,     bmax);
    MYMALLOC(mst_queue,   bmax);
    MYMALLOC(mst_brngs,   bmax);
    MYCALLOC(smap,        nseq);
    for (i = 0; i < nseq; i++)
        smap[i] = -1;

    for (i = 0; i < nblk; i++) {
        if (skip && skip[i]) continue;

        blk   = blks + i;
        bnseq = blk->nseq;
        bseqs = blk->seqs;

        if (bnseq == 0) continue;
        if (bnseq == 1) {
            a = bseqs[0] >> 1;
            haps[a].bpos = 0;
            haps[a].epos = dicts->s[a].len;
            continue;
        }

        for (j = 0; j < bnseq; j++)
            smap[bseqs[j]>>1] = j;

        me.n = 0;
        for (j = 0; j < bnseq; j++) {
            a  = bseqs[j] >> 1;
            ov = ovls + (index[a] >> 32);
            n  = (uint32) index[a];
            for (k = 0; k < n; k++, ov++) {
                if (ov->del) continue;
                b  = ov->bread;
                ib = smap[b];
                if (ib < 0 || ib <= j) continue;
                kv_push(mst_edge_t, me, ((mst_edge_t){(int) j, ib, (ov->alen * ov->qual)}));
            }
        }

        qsort(me.a, me.n, sizeof(mst_edge_t), mstedge_w_dcmpfunc);

        for (j = 0; j < bnseq; j++) { uf_r[j] = -1; mst_deg[j] = 0; }
        me_n_agg  = me.n;
        mst_start = me.n;
        n_mst = 0;
        for (j = 0; j < me_n_agg && n_mst < bnseq-1; j++) {
            ia = me.a[j].a; ib = me.a[j].b;
            ra = ia; rb = ib;
            while (uf_r[ra] >= 0) ra = uf_r[ra];
            while (uf_r[rb] >= 0) rb = uf_r[rb];
            if (ra == rb) continue;
            if (-uf_r[ra] < -uf_r[rb]) { tmp = ra; ra = rb; rb = tmp; }
            if (uf_r[ra] == uf_r[rb]) uf_r[ra]--;
            uf_r[rb] = ra;
            kv_push(mst_edge_t, me, me.a[j]);
            mst_deg[ia]++; mst_deg[ib]++;
            n_mst++;
        }

        root = 0;
        for (j = 0; j < bnseq; j++) {
            if (mst_deg[j] == 1) { root = j; break; }
        }
        for (j = 0; j < bnseq; j++) mst_par[j] = bnseq;
        mst_par[root]     = -1;
        mst_par_ovl[root] = NULL;
        mst_queue[0] = root;
        nq = 1;
        for (p = 0; p < nq; p++) {
            cur = mst_queue[p];
            for (j = 0; j < n_mst; j++) {
                ia = me.a[mst_start + j].a;
                ib = me.a[mst_start + j].b;
                nb = -1;
                if (ia == cur && mst_par[ib] == bnseq) nb = ib;
                else if (ib == cur && mst_par[ia] == bnseq) nb = ia;
                if (nb < 0) continue;
                mst_par[nb] = cur;
                mst_par_ovl[nb] = find_b_overlap(
                    ovls + (index[bseqs[nb]>>1] >> 32),
                    (uint32) index[bseqs[nb]>>1],
                    bseqs[cur]>>1);
                mst_queue[nq++] = nb;
            }
        }

        a = bseqs[root] >> 1;
        mst_brngs[root].beg = 0;
        mst_brngs[root].end = dicts->s[a].len;
        haps[a].bpos = 0;
        haps[a].epos = dicts->s[a].len;

        for (p = 1; p < nq; p++) {
            cur = mst_queue[p];
            par = mst_par[cur];
            a = bseqs[cur] >> 1;
            b = bseqs[par] >> 1;
            bb_alen = dicts->s[a].len;
            bb_blen = dicts->s[b].len;
            bb_r    = bseqs[cur] & 1;
            bb_s    = bseqs[par] & 1;
            oc = mst_par_ovl[cur];
            if (oc) {
                if (bb_r) { bb_bgap = bb_alen - oc->aepos; bb_egap = oc->abpos; }
                else      { bb_bgap = oc->abpos; bb_egap = bb_alen - oc->aepos; }
                if (bb_s) { bb_bshift = bb_blen - oc->bepos; bb_eshift = bb_blen - oc->bbpos; }
                else      { bb_bshift = oc->bbpos; bb_eshift = oc->bepos; }
                mst_brngs[cur].beg = mst_brngs[par].beg + bb_bshift - bb_bgap;
                mst_brngs[cur].end = mst_brngs[par].beg + bb_eshift + bb_egap;
            } else {
                mst_brngs[cur].beg = mst_brngs[par].end;
                mst_brngs[cur].end = mst_brngs[par].end + bb_alen;
            }
            haps[a].bpos = mst_brngs[cur].beg;
            haps[a].epos = mst_brngs[cur].end;
        }

        if (nq < bnseq) {
            max_epos = 0;
            for (j = 0; j < bnseq; j++) {
                a = bseqs[j] >> 1;
                if (haps[a].epos > max_epos) max_epos = haps[a].epos;
            }
            for (j = 0; j < bnseq; j++) {
                if (mst_par[j] != bnseq) continue;
                a = bseqs[j] >> 1;
                haps[a].bpos = max_epos;
                haps[a].epos = max_epos + dicts->s[a].len;
                max_epos = haps[a].epos;
            }
        }

        slen = 0;
        for (j = 0; j < bnseq; j++) {
            a = bseqs[j] >> 1;
            if (haps[a].bpos < slen) slen = haps[a].bpos;
        }
        if (slen < 0) {
            for (j = 0; j < bnseq; j++) {
                a = bseqs[j] >> 1;
                haps[a].bpos -= slen;
                haps[a].epos -= slen;
            }
        }

        for (j = 0; j < bnseq; j++) smap[bseqs[j]>>1] = -1;
    }

    free(uf_r);
    free(smap);
    free(mst_par);
    free(mst_deg);
    free(mst_queue);
    free(mst_brngs);
    free(mst_par_ovl);
    kv_destroy(me);
}

static void build_scaffold_partition(ovl_t *ovls, int64 novl, hap_info_t *haps, sdict_t *dicts, 
    busco_table_t *buscos, int ploidy, int min_ext, int with_hap)
{
    if (ovls == NULL || !novl)
        return;

    kvec_t(scf_block_t) blks;
    scf_block_t *blk;
    ovl_t *ovl;
    range_t *rngs;
    cov_point_t **covs, *cov;
    ord_i64_t *sords, *qords;
    kvec_t(ext_info_t) exts;
    kvec_t(ext_clus_t) clus;
    uint64 *index;
    uint32 s, a, b, *cseqs, *ctype;
    uint8 *rmarks, *qmarks;
    int64 i, j, k, l, p, q, g, n, movl, nseq, nrng, slen, ncseq, lowcpy, fovls[2];
    int ga, gb, beg, end, ext, ngrp, *ncov, *grps, *cpts, *label;
    double maxcpy, l_scaled, min_qual = OVL_MIN_QUAL;

    // sequence number
    nseq = dicts->n;
    
    // allocate memory
    MYCALLOC(index, nseq);
    MYMALLOC(covs, nseq);
    MYMALLOC(ncov, nseq);
    MYMALLOC(covs[0], (novl+nseq)*2);
    MYMALLOC(grps, nseq);
    MYMALLOC(sords, nseq);
    MYMALLOC(qords, nseq);
    MYMALLOC(cseqs, nseq);
    MYMALLOC(ctype, nseq);
    MYMALLOC(rmarks, nseq);
    MYMALLOC(qmarks, nseq);
    
    if (covs == NULL || covs[0] == NULL || ncov == NULL || 
        grps == NULL || sords == NULL || qords == NULL ||
        index == NULL || cseqs == NULL || ctype == NULL || 
        rmarks == NULL || qmarks == NULL)
        mem_alloc_error("scaffold partition");
    
    // build sequence orders by length
    for (i = 0; i < nseq; i++)
        sords[i] = (ord_i64_t) {i, dicts->s[i].len};
    qsort(sords, nseq, sizeof(ord_i64_t), ord_i64_dcmpfunc);

    // build overlap index for quick search
    // overlaps are sorted by abseqs
    a = ovls->aread;
    for (i = 1, j = 0; i < novl; i++) {
        if (ovls[i].aread != a) {
            index[a] = (uint64) j << 32 | (i - j);
            a = ovls[i].aread;
            j = i;
        }
    }
    index[a] = (uint64) j << 32 | (i - j);

    // max overlaps
    movl = 0;
    for (i = 0; i < nseq; i++) {
        n = index[i] >> 32;
        if (movl < n)
            movl = n;
    }
    MYMALLOC(rngs, movl);
    MYMALLOC(cpts, movl*2);
    if (rngs == NULL || cpts == NULL)
        mem_alloc_error("scaffold partition");
    // calculate sequence coverage
    for (i = 1; i < nseq; i++)
        covs[i] = covs[i-1] + (uint32) index[i-1] * 2 + 2;
    for (i = 0; i < nseq; i++) {
        ovl = ovls + (index[i] >> 32);
        n = (uint32) index[i];
        ncov[i] = calc_coverage_from_intervals(covs[i], cpts, dicts->s[i].len, ovl, n, &min_qual, pts_from_overlaps);
    }

    // average coverage
    maxcpy = 1. + ploidy;
    for (i = 0; i < nseq; i++) {
        cov = covs[i];
        lowcpy = 0;
        for (j = 1; j < ncov[i]; j++)
            if (cov[j-1].cov < maxcpy)
                lowcpy += cov[j].pos - cov[j-1].pos;
        // mark repeat
        rmarks[i] = lowcpy < dicts->s[i].len * overlap_min_lowcopy_ratio(dicts->s[i].len);
    }

    // build scaffold partition
    kv_init(blks);
    kv_init(exts); kv_resize(ext_info_t, exts, ploidy * 2);
    kv_init(clus); kv_resize(ext_clus_t, clus, ploidy * 2);
    for (i = 0; i < nseq; i++)
        grps[i] = 0;
    ngrp = 0;
    for (i = 0; i < nseq; i++) {
        a = sords[i].which;
        if (grps[a]) continue;
        // build scaffold backbone use a as seed
        grps[a] = ++ngrp;
        if (rmarks[a]) {
            a = a << 1;
            ncseq = 0;
            cseqs[ncseq] = a, ctype[ncseq] = 0, ncseq++;
            goto add_block; // likely repeat, skip
        }
        // left extension
        a = a << 1 | 1;
        ncseq = 0;
        cseqs[ncseq] = a, ctype[ncseq] = 0, ncseq++;
        while (1) {
            ovl = ovls + (index[a>>1] >> 32);
            n = (uint32) index[a>>1];
            exts.n = 0;
            for (j = 0; j < n; j++, ovl++) {
                b = ovl->bread;
                if (grps[b] ||
                    rmarks[b] ||
                    ovl->del || 
                    !(ovl->type & OVL_EXTD) ||
                    ovl->arev != (a & 1) ||
                    ovl->aepos - ovl->abpos < min_ext ||
                    ovl->bepos - ovl->bbpos < min_ext)
                    continue;
                if (ovl->brev) {
                    beg = 0;
                    end = ovl->bbpos;
                } else {
                    beg = ovl->bepos;
                    end = dicts->s[b].len;
                }
                ext = end - beg;
                if (ext < min_ext) continue;
                if (average_range_coverage(covs[a>>1], ncov[a>>1], ovl->abpos, ovl->aepos) > maxcpy)
                    continue;
                if (average_range_coverage(covs[b], ncov[b], ovl->bbpos, ovl->bepos) > maxcpy)
                    continue;
                kv_push(ext_info_t, exts, ((ext_info_t) {b << 1 | ovl->brev, ovl, ovl->brev? end : beg, 0, ovl->blen, ext}));
            }
            // group candidates
            if (!exts.n) break;
            qsort(exts.a, exts.n, sizeof(ext_info_t), einfo_ext_dcmpfunc);
            g = 0;
            for (j = 0; j < exts.n; j++) {
                for (k = 0; k < j; k++) {
                    if (overlap_accorss_position(ovls, index, exts.a[j].s, exts.a[j].p, exts.a[k].s, exts.a[k].p, min_ext)) {
                        if (!exts.a[j].g)
                            exts.a[j].g = exts.a[k].g;
                        else if (exts.a[j].g != exts.a[k].g) {
                            // merge groups
                            p = MIN(exts.a[j].g, exts.a[k].g);
                            q = MAX(exts.a[j].g, exts.a[k].g);
                            for (l = 0; l <= j; l++) {
                                if (exts.a[l].g == q)
                                    exts.a[l].g = p;
                            }
                        }
                    }
                }
                if (!exts.a[j].g) exts.a[j].g = ++g;
            }
            // collect groups
            clus.n = 0;
            g = 0;
            for (j = 0; j < exts.n; j++) {
                if (g < exts.a[j].g) {
                    g = exts.a[j].g;
                    kv_push(ext_clus_t, clus, ((ext_clus_t) {g, exts.a[j].ovl}));
                } else {
                    for (k = 0; k < clus.n; k++)
                        if (clus.a[k].g == exts.a[j].g)
                            break;
                    // k must exist
                    clus.a[k].score += exts.a[j].ovl;
                }
            }
            // score and filter groups
            qsort(clus.a, clus.n, sizeof(ext_clus_t), eclus_score_dcmpfunc);
            if (clus.n > 1 && clus.a[0].score * (1. - BACKBONE_MIN_CONFIDENCE_SCORE) < clus.a[1].score)
                break; // no clear winner
            // find the best candidate in the winning group
            g = clus.a[0].g;
            s = (uint32) -1;
            for (j = 0; j < exts.n; j++) {
                if (exts.a[j].g == g) {
                    l_scaled = scaled_b_ovl_size(exts.a[j].o);
                    if (l_scaled >= exts.a[j].ext * backbone_min_ovl_to_ext_ratio(l_scaled)) {
                        s = exts.a[j].s;
                        break;
                    }
                }
            }
            if (s == (uint32) -1) break;
            // add to backbone
            a = s;
            grps[a>>1] = (a&1)? -ngrp : ngrp;
            cseqs[ncseq] = a, ctype[ncseq] = 0, ncseq++;
        }
        
        // reverse the array
        for (j = 0, k = ncseq-1; j <= k; j++, k--) {
            s = cseqs[j];
            cseqs[j] = cseqs[k]^1;
            cseqs[k] = s^1;
            // all zeros 
            // s = ctype[j];
            // ctype[j] = ctype[k];
            // ctype[k] = s;
        }
        for (j = 0; j < ncseq; j++) {
            s = cseqs[j];
            grps[s>>1] = (s&1)? -ngrp : ngrp;
        }

        // right extension
        a = sords[i].which << 1;
        while (1) {
            ovl = ovls + (index[a>>1] >> 32);
            n = (uint32) index[a>>1];
            exts.n = 0;
            for (j = 0; j < n; j++, ovl++) {
                b = ovl->bread;
                if (grps[b] ||
                    rmarks[b] ||
                    ovl->del || 
                    !(ovl->type & OVL_EXTD) ||
                    ovl->arev != (a & 1) ||
                    ovl->aepos - ovl->abpos < min_ext ||
                    ovl->bepos - ovl->bbpos < min_ext)
                    continue;
                if (ovl->brev) {
                    beg = 0;
                    end = ovl->bbpos;
                } else {
                    beg = ovl->bepos;
                    end = dicts->s[b].len;
                }
                ext = end - beg;
                if (ext < min_ext) continue;
                if (average_range_coverage(covs[a>>1], ncov[a>>1], ovl->abpos, ovl->aepos) > maxcpy)
                    continue;
                if (average_range_coverage(covs[b], ncov[b], ovl->bbpos, ovl->bepos) > maxcpy)
                    continue;
                kv_push(ext_info_t, exts, ((ext_info_t) {b << 1 | ovl->brev, ovl, ovl->brev? end : beg, 0, ovl->blen, ext}));
            }
            // group candidates
            if (!exts.n) break;
            qsort(exts.a, exts.n, sizeof(ext_info_t), einfo_ext_dcmpfunc);
            g = 0;
            for (j = 0; j < exts.n; j++) {
                for (k = 0; k < j; k++) {
                    if (overlap_accorss_position(ovls, index, exts.a[j].s, exts.a[j].p, exts.a[k].s, exts.a[k].p, min_ext)) {
                        if (!exts.a[j].g)
                            exts.a[j].g = exts.a[k].g;
                        else if (exts.a[j].g != exts.a[k].g) {
                            // merge groups
                            p = MIN(exts.a[j].g, exts.a[k].g);
                            q = MAX(exts.a[j].g, exts.a[k].g);
                            for (l = 0; l <= j; l++) {
                                if (exts.a[l].g == q)
                                    exts.a[l].g = p;
                            }
                        }
                    }
                }
                if (!exts.a[j].g) exts.a[j].g = ++g;
            }
            // collect groups
            clus.n = 0;
            g = 0;
            for (j = 0; j < exts.n; j++) {
                if (g < exts.a[j].g) {
                    g = exts.a[j].g;
                    kv_push(ext_clus_t, clus, ((ext_clus_t) {g, exts.a[j].ovl}));
                } else {
                    for (k = 0; k < clus.n; k++)
                        if (clus.a[k].g == exts.a[j].g)
                            break;
                    // k must exist
                    clus.a[k].score += exts.a[j].ovl;
                }
            }
            // score and filter groups
            qsort(clus.a, clus.n, sizeof(ext_clus_t), eclus_score_dcmpfunc);
            if (clus.n > 1 && clus.a[0].score * (1. - BACKBONE_MIN_CONFIDENCE_SCORE) < clus.a[1].score)
                break; // no clear winner
            // find the best candidate in the winning group
            g = clus.a[0].g;
            s = (uint32) -1;
            for (j = 0; j < exts.n; j++) {
                if (exts.a[j].g == g) {
                    l_scaled = scaled_b_ovl_size(exts.a[j].o);
                    if (l_scaled >= exts.a[j].ext * backbone_min_ovl_to_ext_ratio(l_scaled)) {
                        s = exts.a[j].s;
                        break;
                    }
                }
            }
            if (s == (uint32) -1) break;
            // add to backbone
            a = s;
            grps[a>>1] = (a&1)? -ngrp : ngrp;
            cseqs[ncseq] = a, ctype[ncseq] = 0, ncseq++;
        }

#ifdef DEBUG_SCAFFOLD_PARTITION
        slen = 0;
        for (j = 0; j < ncseq; j++)
            slen += dicts->s[cseqs[j]>>1].len;
        fprintf(stderr, "[M::%s] scaffold group %d: %lld sequences [%lld bp]\n", __func__, ngrp, ncseq, slen);
        for (j = 0; j < ncseq; j++) {
            fprintf(stderr, "[M::%s] [%4lld] %s%c %12d\n", __func__, 
                j, dicts->s[cseqs[j]>>1].name, "+-"[cseqs[j]&1], dicts->s[cseqs[j]>>1].len);
        }
#endif

        // recruit sequences from backbone sequences
        MYBZERO(qmarks, nseq);
        p = q = 0;
        g = 0;
        while (1) {
            // collect candidate sequences
            for (j = p; j < ncseq; j++) {
                a = cseqs[j];
                ovl = ovls + (index[a>>1] >> 32);
                n = (uint32) index[a>>1];
                for (k = 0; k < n; k++, ovl++) {
                    b = ovl->bread;
                    if (grps[b] ||
                        rmarks[b] ||
                        qmarks[b] ||
                        ovl->del)
                        continue;
                    qmarks[b] = 1;
                    qords[q++] = (ord_i64_t) {b, dicts->s[b].len};
                }
            }
            if (!q) break; // no candidates
            p = ncseq; // mark the end
            // process candidates
            qsort(qords, q, sizeof(ord_i64_t), ord_i64_dcmpfunc);
            k = 0;
            g++;
            for (j = 0; j < q; j++) {
                a = qords[j].which;
                if (grps[a]) continue;
                ovl = ovls + (index[a] >> 32);
                n = (uint32) index[a];
                nrng = 0;
                fovls[0] = fovls[1] = 0;
                for (l = 0; l < n; l++, ovl++) {
                    b = ovl->bread;
                    if (abs(grps[b]) != ngrp ||
                        ovl->del)
                        continue;
                    if (grps[b] > 0)
                        fovls[ovl->arev!=ovl->brev] += ovl->alen;
                    else
                        fovls[ovl->arev==ovl->brev] += ovl->alen;
                    rngs[nrng++] = (range_t) {ovl->abpos, ovl->aepos};
                }
                if (nrng && rangelist_size(rngs, nrng, 0) >= dicts->s[a].len * BACKBONE_MIN_RECRUIT) {
                    // grps[a] = fovls[1] > fovls[0]? -ngrp : ngrp;
                    cseqs[ncseq] = (a << 1) | (fovls[1] > fovls[0]), ctype[ncseq] = g, ncseq++;
                } else qords[k++] = qords[j]; // keep for next round
            }
            // set groups outside the loop for a strict recruitment order
            for (j = p; j < ncseq; j++) {
                a = cseqs[j];
                grps[a>>1] = (a&1)? -ngrp : ngrp;
            }

            q = k; // set end of candidates for next round
            if (p == ncseq)
                break; // no new sequences added
        }

#ifdef DEBUG_SCAFFOLD_PARTITION
        {
            for (j = 0; j < q; j++) {
                a = qords[j].which;
                if (grps[a]) continue;
                ovl = ovls + (index[a] >> 32);
                n = (uint32) index[a];
                nrng = 0;
                fovls[0] = fovls[1] = 0;
                for (l = 0; l < n; l++, ovl++) {
                    b = ovl->bread;
                    if (abs(grps[b]) != ngrp ||
                        ovl->del)
                        continue;
                    if (grps[b] > 0)
                        fovls[ovl->arev!=ovl->brev] += ovl->alen;
                    else
                        fovls[ovl->arev==ovl->brev] += ovl->alen;
                    rngs[nrng++] = (range_t) {ovl->abpos, ovl->aepos};
                }
                slen = nrng? rangelist_size(rngs, nrng, 0) : 0;
                fprintf(stderr, "[M::%s] skipped candidate %s: %d tot=%lld [%.2f] fwd=%lld rev=%lld\n",
                    __func__, dicts->s[a].name, dicts->s[a].len, slen, (double) slen / dicts->s[a].len, fovls[0], fovls[1]);
            }
        }
#endif

add_block:
        slen = 0;
        for (j = 0; j < ncseq; j++)
            slen += dicts->s[cseqs[j]>>1].len;
        
        kv_pushp(scf_block_t, blks, &blk);
        blk->nseq = ncseq;
        MYMALLOC(blk->seqs, ncseq), memcpy(blk->seqs, cseqs, ncseq * sizeof(uint32));
        MYMALLOC(blk->type, ncseq), memcpy(blk->type, ctype, ncseq * sizeof(uint32));
        blk->slen = slen;

#ifdef DEBUG_SCAFFOLD_PARTITION
        fprintf(stderr, "[M::%s] scaffold group %d after recruitment: %lld sequences [%lld bp]\n", __func__, ngrp, ncseq, slen);
#endif
    }

    // order groups within-group according to backbone positions
    {
        bb_range_t *bbrngs;
        ovl_t *bb_o;
        int64 bp, bmax, bb_cpos, bpos_min;
        uint32 bb_r, bb_v, bb_s, bb_t, bb_n, bb_ww, *bseqs, *btype;
        double bb_bsum, bb_esum, bb_ssum;
        int bb_x, bb_y, bnseq, bb_alen, bb_blen, bb_bshift, bb_eshift, bb_bgap, bb_egap, *smap_blk;
        
        // find largest block for bbrngs allocation
        bmax = 0;
        for (i = 0; i < blks.n; i++)
            if (blks.a[i].nseq > bmax) bmax = blks.a[i].nseq;

        MYMALLOC(bbrngs,    bmax);
        MYCALLOC(smap_blk,  nseq);
        for (i = 0; i < nseq; i++) smap_blk[i] = -1;

        // compute backbone positions within each group
        for (i = 0; i < blks.n; i++) {
            blk = blks.a + i;
            bseqs = blk->seqs;
            btype = blk->type;
            bnseq = blk->nseq;

            // build smap for this group
            for (j = 0; j < bnseq; j++)
                smap_blk[bseqs[j]>>1] = j;

            for (j = 0; j < bnseq; j++) {
                a = bseqs[j] >> 1;
                if (!btype[j]) {
                    // backbone sequence: position follows previous backbone end
                    if (j > 0) {
                        bb_cpos = bbrngs[j-1].end;
                        bb_o = find_b_overlap(ovls + (index[a] >> 32), (uint32) index[a], bseqs[j-1]>>1);
                        if (bb_o) bb_cpos -= bb_o->alen;
                    } else bb_cpos = 0;
                    bbrngs[j] = (bb_range_t){j, bb_cpos, bb_cpos + dicts->s[a].len};
                } else {
                    // recruited sequence: weighted anchor from already-placed sequences
                    bb_r   = bseqs[j] & 1;
                    bb_v   = btype[j];
                    bb_x   = smap_blk[a];
                    bb_alen = dicts->s[a].len;
                    bb_o   = ovls + (index[a] >> 32);
                    bb_n   = (uint32) index[a];
                    bb_bsum = bb_esum = bb_ssum = 0.;
                    bb_t   = UINT32_MAX;
                    for (k = 0; k < bb_n; k++, bb_o++) {
                        if (bb_o->del) continue;
                        b = bb_o->bread;
                        bb_y = smap_blk[b];
                        if (bb_y < 0 || bb_y > bb_x) continue;
                        bb_ww = btype[bb_y];
                        if (bb_ww >= bb_v || bb_ww > bb_t) continue;
                        bb_s = bseqs[bb_y] & 1;
                        if ((bb_r == bb_s) != (bb_o->arev == bb_o->brev)) continue;
                        if (bb_ww < bb_t) { bb_bsum = bb_esum = bb_ssum = 0.; bb_t = bb_ww; }
                        bb_blen  = dicts->s[b].len;
                        if (bb_r) { bb_bgap = bb_alen - bb_o->aepos; bb_egap = bb_o->abpos; }
                        else      { bb_bgap = bb_o->abpos; bb_egap = bb_alen - bb_o->aepos; }
                        if (bb_s) { bb_bshift = bb_blen - bb_o->bepos; bb_eshift = bb_blen - bb_o->bbpos; }
                        else      { bb_bshift = bb_o->bbpos; bb_eshift = bb_o->bepos; }
                        bb_ssum += bb_o->alen;
                        bb_bsum += (double)(bbrngs[bb_y].beg + bb_bshift - bb_bgap) * bb_o->alen;
                        bb_esum += (double)(bbrngs[bb_y].beg + bb_eshift + bb_egap) * bb_o->alen;
                    }
                    assert(bb_ssum > 0.);
                    bbrngs[j] = (bb_range_t){j, (int64)(bb_bsum/bb_ssum), (int64)(bb_esum/bb_ssum)};
                }
                haps[a].grp  = i+1;
                haps[a].rev  = bseqs[j] & 1;
                haps[a].bpos = bbrngs[j].beg;
                haps[a].epos = bbrngs[j].end;
            }

            // safety: ensure all bpos >= 0 (recruited seqs may have negative estimates)
            bpos_min = 0;
            for (j = 0; j < bnseq; j++) {
                bp = haps[bseqs[j]>>1].bpos;
                if (bp < bpos_min) bpos_min = bp;
            }
            if (bpos_min < 0) {
                for (j = 0; j < bnseq; j++) {
                    haps[bseqs[j]>>1].bpos -= bpos_min;
                    haps[bseqs[j]>>1].epos -= bpos_min;
                }
            }

            // clear smap for this group
            for (j = 0; j < bnseq; j++)
                smap_blk[bseqs[j]>>1] = -1;
        }

        free(bbrngs);
        free(smap_blk);
    }

    // QPBO via Graph Cut to orient sequences globally
    {
        grp_ovl_t *govls;
        grp_edge_t *edge;
        grp_edge_vec_t *gadjs;
        kvec_t(grp_conflict_t) confs;
        int64 blk_span, nb, ne;
        int cur, n_in_tree, new_root, *prim_s, *prim_r, *in_tree;
        double w, bw, *pathw, *primw;
        
        MYCALLOC(govls, novl);
        n = 0;
        for (i = 0; i < novl; i++) {
            ovl = ovls + i;
            if (ovl->del) continue;
            a = ovl->aread;
            b = ovl->bread;
            ga = abs(grps[a]);
            gb = abs(grps[b]);
            if (ga && gb && ga < gb) {
                govls[n].a = ga;
                govls[n].b = gb;
                govls[n].o = ovl;
                n++;
            }
        }
        qsort(govls, n, sizeof(grp_ovl_t), govl_abseq_cmpfunc);

        if (n > 0) {
            // add edges between each a-b group pair
            MYCALLOC(gadjs, ngrp);
            ga = govls->a;
            gb = govls->b;
            for (i = j = 0; i <= n; i++) {
                if (i == n || govls[i].a != ga || govls[i].b != gb) {
                    if (ga != gb) {
                        l_scaled = 0.;
                        for (k = j; k < i; k++) {
                            ovl = govls[k].o;
                            if ((ovl->arev == ovl->brev) == (grps[ovl->aread] > 0) == (grps[ovl->bread] > 0))
                                // support same orientation
                                l_scaled += 0.5 * (ovl->alen + ovl->blen) * ovl->qual;
                            else
                                // support different orientation
                                l_scaled -= 0.5 * (ovl->alen + ovl->blen) * ovl->qual;
                        }
                        if (l_scaled > 0.) {
                            kv_push(grp_edge_t, gadjs[ga-1], ((grp_edge_t) {gb-1, 0, l_scaled, 0.}));
                            kv_push(grp_edge_t, gadjs[gb-1], ((grp_edge_t) {ga-1, 0, l_scaled, 0.}));
                        } else if (l_scaled < 0.) {
                            kv_push(grp_edge_t, gadjs[ga-1], ((grp_edge_t) {gb-1, 1, -l_scaled, 0.}));
                            kv_push(grp_edge_t, gadjs[gb-1], ((grp_edge_t) {ga-1, 1, -l_scaled, 0.}));
                        }
                    }
                    if (i < n) {
                        ga = govls[i].a;
                        gb = govls[i].b;
                        j = i;
                    }
                }
            }

            // orient groups using maximum spanning tree (Prim's O(V^2))
            kv_init(confs);
            MYMALLOC(label, ngrp);
            MYMALLOC(pathw, ngrp);
            MYMALLOC(primw, ngrp);
            MYMALLOC(prim_s, ngrp);
            MYCALLOC(prim_r, ngrp);
            MYCALLOC(in_tree, ngrp);
            for (i = 0; i < ngrp; i++) {
                prim_s[i] = -1;
                primw[i]    = -1.0;
                pathw[i]    = -1.0;
                label[i]    = -1;
            }
            n_in_tree = 0;
            while (n_in_tree < ngrp) {
                // pick max-weight frontier edge, or start a new component
                bw = -1.0;
                cur = -1; new_root = -1;
                for (j = 0; j < ngrp; j++) {
                    if (in_tree[j]) continue;
                    if (prim_s[j] < 0) {
                        if (new_root < 0) new_root = j;
                    } else if (primw[j] > bw) {
                        bw = primw[j];
                        cur = j;
                    }
                }
                if (cur < 0) cur = new_root;
                // assign label and tree-path bottleneck weight
                if (prim_s[cur] < 0) {
                    label[cur] = 0;
                    pathw[cur] = DBL_MAX;
                } else {
                    label[cur] = label[prim_s[cur]] ^ prim_r[cur];
                    pathw[cur] = MIN(primw[cur], pathw[prim_s[cur]]);
                }
                in_tree[cur] = 1;
                n_in_tree++;
                // update frontier with neighbors of cur
                for (j = 0; j < gadjs[cur].n; j++) {
                    edge = gadjs[cur].a + j;
                    gb = edge->g;
                    if (in_tree[gb]) continue;
                    if (edge->w > primw[gb]) {
                        primw[gb]    = edge->w;
                        prim_s[gb] = cur;
                        prim_r[gb]    = edge->r;
                    }
                }
            }

            // detect frustrated non-tree edges (edge_w <= cut_w after MST)
            for (ga = 0; ga < ngrp; ga++) {
                for (j = 0; j < gadjs[ga].n; j++) {
                    edge = gadjs[ga].a + j;
                    gb = edge->g;
                    if (ga >= gb) continue;
                    if ((label[ga] ^ (int)edge->r) != label[gb]) {
                        w = MIN(pathw[ga], pathw[gb]);
                        kv_push(grp_conflict_t, confs, ((grp_conflict_t) {ga, gb, edge->w, w}));
                    }
                }
            }

            // flip sequences orientation if necessary
            // also the sequence order in the block
            for (i = 0; i < blks.n; i++) {
                if (label[i] != 1)
                    continue;
                blk = blks.a + i;
                for (j = 0; j < blk->nseq; j++) {
                    s = blk->seqs[j]>>1;
                    grps[s] = -grps[s];
                    blk->seqs[j] ^= 1;
                }
                for (j = 0, k = blk->nseq-1; j <= k; j++, k--) {
                    s = blk->seqs[j];
                    blk->seqs[j] = blk->seqs[k];
                    blk->seqs[k] = s;
                    s = blk->type[j];
                    blk->type[j] = blk->type[k];
                    blk->type[k] = s;
                }
                // reflect haps[].bpos/epos: new_bpos = max_epos - old_epos
                blk_span = 0;
                for (j = 0; j < blk->nseq; j++) {
                    s = blk->seqs[j]>>1;
                    if (haps[s].epos > blk_span)
                        blk_span = haps[s].epos;
                }
                for (j = 0; j < blk->nseq; j++) {
                    s = blk->seqs[j]>>1;
                    nb = haps[s].bpos;
                    ne = haps[s].epos;
                    haps[s].bpos = blk_span - ne;
                    haps[s].epos = blk_span - nb;
                    haps[s].rev  = blk->seqs[j] & 1;
                }
            }

            kv_destroy(confs);
            for (i = 0; i < ngrp; i++)
                kv_destroy(gadjs[i]);
            free(gadjs);
            free(label);
            free(pathw);
            free(primw);
            free(prim_s); 
            free(prim_r); 
            free(in_tree);
        }
        free(govls);
    }

    // ploidy-aware group merging
    {
        grp_ovl_t *govls;
        range_vec_t rngv;
        srange_vec_t srngv;
        HeapSpace *hs;
        olink_t *olinks, *olink;
        gadj_t **gadjs, *gadj, *aadj, **adjs;
        gadj_pool_t *gadjp;
        cov_point_t *tcov;
        int64 npair;
        int tncov, nadj, *tbcov, **bcov, *states, *unions;

        // collect inter-group overlaps
        MYCALLOC(govls, novl);
        npair = 0;
        for (i = 0; i < novl; i++) {
            ovl = ovls + i;
            if (ovl->del) continue;
            a = ovl->aread;
            b = ovl->bread;
            ga = abs(grps[a]);
            gb = abs(grps[b]);
            if (ga && gb && ga < gb) {
                govls[npair].a = ga-1;
                govls[npair].b = gb-1;
                govls[npair].o = ovl;
                npair++;
            }
        }
        qsort(govls, npair, sizeof(grp_ovl_t), govl_abseq_cmpfunc);

        hs = hs_init(ngrp, 1<<15);
        if (npair > 0) {
            // intra-group coverage
            kv_init(rngv);
            kv_resize(range_t, rngv, movl);
            for (i = 0; i < nseq; i++) {
                ovl = ovls + (index[i] >> 32);
                n = (uint32) index[i];
                g = abs(grps[i]);
                rngv.n = 0;
                for (j = 0; j < n; j++, ovl++)
                    if (!ovl->del && abs(grps[ovl->bread]) == g)
                        kv_push(range_t, rngv, ((range_t) {ovl->abpos, ovl->aepos}));
                ncov[i] = calc_coverage_from_intervals(covs[i], cpts, dicts->s[i].len, rngv.a, rngv.n, NULL, pts_from_ranges);
            }

            // profile coverage for each group
            MYCALLOC(bcov, ngrp);
            MYCALLOC(bcov[0], ngrp*(ploidy+1));
            for (i = 1; i < ngrp; i++)
                bcov[i] = bcov[i-1] + (ploidy+1);
            for (i = 0; i < ngrp; i++) {
                blk = blks.a + i;
                MYBZERO(bcov[i], ploidy+1);
                for (j = 0; j < blk->nseq; j++) {
                    a = blk->seqs[j]>>1;
                    profile_coverage(covs[a], ncov[a], ploidy, bcov[i], 0);
                }  
            }

            // for union-find of group merging
            MYCALLOC(unions, ngrp);
            for (i = 0; i < ngrp; i++)
                unions[i] = -1;

            // add edges between each a-b group pair
            kv_init(srngv);
            MYCALLOC(states, ngrp);
            MYCALLOC(gadjs, ngrp);
            MYMALLOC(adjs, ngrp*2);
            MYMALLOC(tcov, (movl+1)*6);
            MYCALLOC(olinks, npair);
            MYMALLOC(tbcov, ploidy+1);
            MYCALLOC(gadjp, 1);
            gadj = gadj_new(gadjp);
            ga = govls->a;
            gb = govls->b;
            for (i = 1, j = 0; i <= npair; i++) {
                if (i == npair || govls[i].a != ga || govls[i].b != gb) {
                    olinks[j].ovl = govls[j].o;
                    for (k = j+1; k < i; k++) {
                        olinks[k-1].next = olinks + k;
                        olinks[k].ovl = govls[k].o;
                    }
                    
                    gadj->a = ga;
                    gadj->b = gb;
                    gadj->ast = states[ga];
                    gadj->bst = states[gb];
                    gadj->olinks = olinks + j;
                    gadj->olast = olinks + i - 1;

                    if (mergeable_group_pair(gadj, covs, ncov, bcov[ga], bcov[gb], haps, ploidy, tbcov, tcov, cpts, &srngv)) {
                        gadj->next = gadjs[ga];
                        gadjs[ga] = gadj;
                        hs_heap_push(hs, gadj);
                        gadj = gadj_copy(gadjp, gadj);
                        gadj->next = gadjs[gb];
                        gadjs[gb] = gadj;
                        gadj = gadj_new(gadjp);
                    } else MYBZERO(gadj, 1);
                    if (i < npair) {
                        ga = govls[i].a;
                        gb = govls[i].b;
                        j = i;
                    }
                }
            }

            while ((gadj = hs_heap_pop(hs))) {
                ga = gadj->a;
                gb = gadj->b;

                if (gadj->ast != states[ga] || gadj->bst != states[gb])
                    continue; // stale entry

                // merge gb into ga
                unions[gb] = ga;
                
                // update sequence profiles
                olink = gadj->olinks;
                srngv.n = 0;
                while (olink) {
                    ovl = olink->ovl;
                    kv_push(srange_t, srngv, ((srange_t){ovl->aread, ovl->abpos, ovl->aepos}));
                    kv_push(srange_t, srngv, ((srange_t){ovl->bread, ovl->bbpos, ovl->bepos}));
                    olink = olink->next;
                }
                qsort(srngv.a, srngv.n, sizeof(srange_t), srange_cmpfunc);
                a = srngv.a[0].seq;
                for (i = 1, j = 0; i <= srngv.n; i++) {
                    if (i == srngv.n || srngv.a[i].seq != a) {
                        tncov = calc_coverage_from_intervals(tcov, cpts, haps[a].len, srngv.a+j, i-j, NULL, pts_from_sranges);
                        ncov[a] = merge_coverage(covs[a], ncov[a], tcov, tncov, tcov + tncov);
                        if (i < srngv.n)
                            a = srngv.a[i].seq;
                        j = i;
                    }
                }

                // update group coverage profiles
                for (i = 0; i <= ploidy; i++)
                    bcov[ga][i] += bcov[gb][i] + gadj->bcov[i];

                // collect active edges
                nadj = 0;
                aadj = gadjs[ga];
                while (aadj) {
                    if ((aadj->a == ga && aadj->b != gb && aadj->ast == states[ga] && aadj->bst == states[aadj->b]) || 
                        (aadj->b == ga && aadj->a != gb && aadj->ast == states[aadj->a] && aadj->bst == states[ga]))
                        adjs[nadj++] = aadj;
                    aadj = aadj->next;
                }
                aadj = gadjs[gb];
                while (aadj) {
                    if ((aadj->a == gb && aadj->b != ga && aadj->ast == states[gb] && aadj->bst == states[aadj->b]) || 
                        (aadj->b == gb && aadj->a != ga && aadj->ast == states[aadj->a] && aadj->bst == states[gb]))
                        adjs[nadj++] = aadj;
                    aadj = aadj->next;
                }

                // update states
                states[ga]++;
                states[gb]++;
                // all old edges are now stale
                gadjs[ga] = NULL;
                // push new edges
                if (nadj) {
                    // make b the neighbor for sorting
                    for (i = 0; i < nadj; i++) {
                        aadj = adjs[i];
                        if (aadj->b == ga) {
                            aadj->b = aadj->a;
                            aadj->a = ga;
                        } else if (aadj->b == gb) {
                            aadj->b = aadj->a;
                            aadj->a = gb;
                        }
                    }
                    // sort by neighbor ids
                    qsort(adjs, nadj, sizeof(gadj_t*), gadj_neighbor_cmpfunc);
                    // process each neighbor
                    for (i = 0; i < nadj; i++) {
                        aadj = adjs[i];
                        gb = aadj->b;
                        if (i < nadj-1 && adjs[i+1]->b == gb) {
                            // shared neighbor
                            // append olink list
                            aadj->olast->next = adjs[i+1]->olinks;
                            aadj->olast = adjs[i+1]->olast;
                            i++;
                        }

                        gadj->a = ga;
                        gadj->b = gb;
                        gadj->ast = states[ga];
                        gadj->bst = states[gb];
                        gadj->olinks = aadj->olinks;
                        gadj->olast = aadj->olast;

                        if (mergeable_group_pair(gadj, covs, ncov, bcov[ga], bcov[gb], haps, ploidy, tbcov, tcov, cpts, &srngv)) {
                            gadj->next = gadjs[ga];
                            gadjs[ga] = gadj;
                            hs_heap_push(hs, gadj);
                            gadj = gadj_copy(gadjp, gadj);
                            gadj->next = gadjs[gb];
                            gadjs[gb] = gadj;
                            gadj = gadj_new(gadjp);
                        } else MYBZERO(gadj, 1);
                    }
                }
            }

            // union groups and assign final group ids
            for (i = 0; i < ngrp; i++) {
                if (unions[i] < 0)
                    continue;
                k = i;
                while (unions[k] >= 0)
                    k = unions[k];
                j = i;
                while (j != k) {
                    unions[j] = k;
                    j = unions[j];
                }
            }

            // update groups
            MYBZERO(sords, ngrp);
            for (i = 0; i < ngrp; i++) {
                k = unions[i] >= 0? unions[i] : i;
                sords[k].event += blks.a[i].nseq;
            }
            // expand memory for sequences in merged groups
            for (i = 0; i < ngrp; i++) {
                if (unions[i] >= 0)
                    continue;
                // number sequences
                n = sords[i].event;
                blk = blks.a + i;
                MYREALLOC(blk->seqs, n);
                MYREALLOC(blk->type, n);
            }
            for (i = 0; i < ngrp; i++) {
                if (unions[i] < 0)
                    continue;
                blk = blks.a + unions[i];
                memcpy(blk->seqs + blk->nseq, blks.a[i].seqs, blks.a[i].nseq * sizeof(uint32));
                memcpy(blk->type + blk->nseq, blks.a[i].type, blks.a[i].nseq * sizeof(uint32));
                blk->nseq += blks.a[i].nseq;
                blk->slen += blks.a[i].slen;
            }
            // compact blocks and free merged blocks
            for (i = j = 0; i < ngrp; i++) {
                if (unions[i] >= 0) {
                    blk = blks.a + i;
                    free(blk->seqs);
                    free(blk->type);
                    continue;
                }
                if (i != j)
                    blks.a[j] = blks.a[i];
                j++;
            }
            ngrp = blks.n = j;
            // update grps[] and haps[].grp/rev
            for (i = 0; i < ngrp; i++) {
                blk = blks.a + i;
                for (j = 0; j < blk->nseq; j++) {
                    s = blk->seqs[j]>>1;
                    haps[s].grp = i+1;
                    haps[s].rev = blk->seqs[j] & 1;
                    grps[s] = haps[s].rev? -(i+1) : (i+1);
                }
            }

            // reorder sequences within each group via MST + BFS ordering
            mst_order_groups(blks.a, ngrp, ovls, index, haps, dicts, nseq, NULL);

            free(olinks);
            free(bcov[0]);
            free(bcov);
            free(tcov);
            free(adjs);
            free(tbcov);
            free(gadjs);
            free(states);
            free(unions);
            kv_destroy(rngv);
            kv_destroy(srngv);
            gadj_pool_destroy(gadjp);
        }
        free(govls);
    }

    // order groups globally
    // minimise Σ adj_w(ga,gb)*|rank(ga)-rank(gb)|
    // (1) spectral ordering: Fiedler vector of normalised Laplacian via power iteration — O(n²·iter)
    // (2) insertion sweep refinement — O(n²) per sweep vs O(n³) for pairwise swap
    {
        scf_block_t *new_blks;
        ord_dbl_t *rords;
        uint32 *bseqs, *btype;
        grp_edge_vec_t *sadj; // sparse adj: sadj[i] = list of (j, w) with adj_w[i,j]>0
        double w;
        int *grp_order, *grp_rank;
        int bnseq;

        MYCALLOC(sadj, ngrp);
        MYMALLOC(rords, ngrp);
        MYMALLOC(new_blks, ngrp);
        MYMALLOC(grp_order, ngrp);
        MYMALLOC(grp_rank, ngrp);
        
        // reorder sequences within each block according to backbone positions
        for (i = 0; i < blks.n; i++) {
            blk = blks.a + i;
            bseqs = blk->seqs;
            btype = blk->type;
            bnseq = blk->nseq;
            for (j = 0; j < bnseq; j++) {
                qords[j].which = j;
                qords[j].event = haps[bseqs[j]>>1].bpos;
            }
            qsort(qords, bnseq, sizeof(ord_i64_t), ord_i64_acmpfunc);
            for (j = 0; j < bnseq; j++) {
                a = qords[j].which;
                cseqs[j] = bseqs[a];
                ctype[j] = btype[a];
            }
            memcpy(bseqs, cseqs, bnseq * sizeof(uint32));
            memcpy(btype, ctype, bnseq * sizeof(uint32));
        }

        // make distance matrix between groups
        {
            double *raw_w, *raw_fwd, *softmax_p, *conf;
            double *row, *prow;
            double wmax, exp_sum, H, H_norm, p_val;
            double pij, pji, rij, v, fij, fji, fwd_ij, wp_ij, wp_ji;
            int64 *blk_minbpos, *blk_maxbpos, spa, spb;
            double fa, fb, da, db, da2, db2, dir_ab, dir_ba;
            int cnt;

            MYCALLOC(raw_w,       (int64)ngrp * ngrp);
            MYCALLOC(raw_fwd,     (int64)ngrp * ngrp);
            MYCALLOC(softmax_p,   (int64)ngrp * ngrp);
            MYCALLOC(conf,        ngrp);
            MYCALLOC(blk_minbpos, ngrp);
            MYCALLOC(blk_maxbpos, ngrp);
            
            for (i = 0; i < ngrp; i++) { blk_minbpos[i] = INT64_MAX; blk_maxbpos[i] = 0; }
            for (i = 0; i < nseq; i++) {
                g = abs(grps[i]);
                if (!g) continue;
                g--;
                if (haps[i].bpos < blk_minbpos[g]) blk_minbpos[g] = haps[i].bpos;
                if (haps[i].epos > blk_maxbpos[g]) blk_maxbpos[g] = haps[i].epos;
            }

            // step 1: accumulate raw = alen * qual and raw_fwd = raw * dir (hard filter: qual >= OVL_MIN_QUAL)
            for (i = 0; i < novl; i++) {
                ovl = ovls + i;
                if (ovl->del) continue;
                if (ovl->qual < OVL_MIN_QUAL) continue;
                a = ovl->aread; b = ovl->bread;
                ga = abs(grps[a]); gb = abs(grps[b]);
                if (!ga || !gb || ga == gb) continue;
                w = 0.5 * (ovl->alen + ovl->blen) * ovl->qual;
                raw_w[(int64)(ga-1)*ngrp + (gb-1)] += w;
                raw_w[(int64)(gb-1)*ngrp + (ga-1)] += w;
                // directional signal: fa = fractional bpos of a in its block, fb = same for b
                {
                    spa = blk_maxbpos[ga-1] - blk_minbpos[ga-1];
                    spb = blk_maxbpos[gb-1] - blk_minbpos[gb-1];
                    fa = (spa > 0) ? (double)(haps[a].bpos - blk_minbpos[ga-1]) / spa : 0.5;
                    fb = (spb > 0) ? (double)(haps[b].bpos - blk_minbpos[gb-1]) / spb : 0.5;
                    // dir_ab > 0: evidence for ga→gb (a is near tail of ga, b near head of gb)
                    da = fa - 0.5; 
                    db = 0.5 - fb;
                    dir_ab = (da > 0. && db > 0.) ? 4. * da * db : 0.;
                    // dir_ba > 0: evidence for gb→ga
                    da2 = fb - 0.5;
                    db2 = 0.5 - fa;
                    dir_ba = (da2 > 0. && db2 > 0.) ? 4. * da2 * db2 : 0.;
                    raw_fwd[(int64)(ga-1)*ngrp + (gb-1)] += w * dir_ab;
                    raw_fwd[(int64)(gb-1)*ngrp + (ga-1)] += w * dir_ba;
                }
            }

            // step 2: per-row softmax(T=max) → softmax_p[i,j]; conf[i] = 1 - H_norm
            // softmax: normalise to [0,1] by wmax then apply temperature T=0.2
            // T controls sharpness: smaller T → sharper; T=0.2 gives exp(-1/0.2)=exp(-5)≈0.007
            // for a value at 0.8*wmax, separating well from values at 0.2*wmax
            for (i = 0; i < ngrp; i++) {
                row  = raw_w     + (int64)i * ngrp;
                prow = softmax_p + (int64)i * ngrp;
                wmax = 0.; cnt = 0;
                for (j = 0; j < ngrp; j++) if (row[j] > wmax) wmax = row[j];
                if (wmax == 0.) { conf[i] = 0.; continue; }
                for (j = 0; j < ngrp; j++) if (row[j] > 0.) cnt++;
                exp_sum = 0.;
                for (j = 0; j < ngrp; j++)
                    if (row[j] > 0.) exp_sum += exp((row[j]/wmax - 1.) / 0.2);
                H = 0.;
                for (j = 0; j < ngrp; j++) {
                    if (row[j] <= 0.) continue;
                    p_val = exp((row[j]/wmax - 1.) / 0.2) / exp_sum;
                    prow[j] = p_val;
                    H -= p_val * log(p_val);
                }
                H_norm = cnt > 1 ? H / log((double)cnt) : 0.;
                conf[i] = 1. - H_norm;
            }

            // step 3: build sparse sadj from gated mutual-attention formula
            // adj_w[i,j]  = conf[i]*conf[j] * p[i,j]*p[j,i] * raw[i,j]  (symmetric)
            // adj_wp[i,j] = adj_w[i,j] * raw_fwd[i,j] / raw[i,j]         (directed: i before j)
            // cost contribution of edge (i,j) given ranks ri, rj:
            //   wp[i,j] * max(0, rj-ri)  +  (w-wp)[i,j] * max(0, ri-rj)
            //   = correct-order: free when wp[i,j] dominates; wrong-order: penalised
            for (i = 0; i < ngrp; i++) {
                for (j = i+1; j < ngrp; j++) {
                    pij = softmax_p[(int64)i*ngrp+j];
                    pji = softmax_p[(int64)j*ngrp+i];
                    rij = 0.5 * (raw_w[(int64)i*ngrp+j] + raw_w[(int64)j*ngrp+i]);
                    v = conf[i] * conf[j] * pij * pji * rij;
                    if (v > 0.) {
                        // fwd fraction for i→j: raw_fwd[i,j] / raw_w[i,j]
                        fij = raw_w[(int64)i*ngrp+j] > 0.
                            ? raw_fwd[(int64)i*ngrp+j] / raw_w[(int64)i*ngrp+j] : 0.5;
                        fji = raw_w[(int64)j*ngrp+i] > 0.
                            ? raw_fwd[(int64)j*ngrp+i] / raw_w[(int64)j*ngrp+i] : 0.5;
                        // average fwd fraction from both sides: fwd_ij ∈ [0,1], large means i before j
                        fwd_ij = 0.5 * (fij + (1. - fji));
                        // directed_cost(w, wp, ri, rj): ri<rj costs wp*d, ri>rj costs (w-wp)*d.
                        // To prefer i before j (fwd_ij large), make the ri<rj branch cheap:
                        //   wp_ij must be SMALL when fwd_ij is large → wp_ij = v*(1-fwd_ij).
                        wp_ij = v * (1. - fwd_ij); // penalty for i before j (low = prefer i first)
                        wp_ji = v * fwd_ij;         // penalty for j before i (low = prefer j first)
                        kv_push(grp_edge_t, sadj[i], ((grp_edge_t){(uint32)j, 0, v, wp_ij}));
                        kv_push(grp_edge_t, sadj[j], ((grp_edge_t){(uint32)i, 0, v, wp_ji}));
                    }
                }
            }
            free(conf);
            free(raw_w);
            free(raw_fwd);
            free(softmax_p);
            free(blk_minbpos);
            free(blk_maxbpos);
        }

        // step A: spectral ordering via power iteration on normalised adjacency
        // finds Fiedler vector of normalised Laplacian = 2nd eigenvector of D^{-1}A
        // by iterating x = D^{-1}A x, projecting out constant eigenvector each step
        {
            double *deg, *x, *xnew;
            double norm, mean, s, diff1, diff2, d1, d2;
            int e, iter;
            MYCALLOC(deg, ngrp); MYCALLOC(x, ngrp); MYCALLOC(xnew, ngrp);
            for (i = 0; i < ngrp; i++)
                for (e = 0; e < (int)sadj[i].n; e++)
                    deg[i] += sadj[i].a[e].w;
            // initialise x with deterministic spread (rank-based)
            for (i = 0; i < ngrp; i++) x[i] = (double)i / ngrp - 0.5;
            norm = 0.;
            for (i = 0; i < ngrp; i++) norm += x[i]*x[i];
            norm = sqrt(norm);
            for (i = 0; i < ngrp; i++) x[i] /= norm;
            for (iter = 0; iter < 500; iter++) {
                // x_new = D^{-1} A x
                for (i = 0; i < ngrp; i++) {
                    if (deg[i] == 0.) { xnew[i] = x[i]; continue; }
                    s = 0.;
                    for (e = 0; e < (int)sadj[i].n; e++) s += sadj[i].a[e].w * x[sadj[i].a[e].g];
                    xnew[i] = s / deg[i];
                }
                // subtract mean (orthogonalise against constant eigenvector 1/√n)
                mean = 0.;
                for (i = 0; i < ngrp; i++) mean += xnew[i];
                mean /= ngrp;
                for (i = 0; i < ngrp; i++) xnew[i] -= mean;
                // normalise
                norm = 0.;
                for (i = 0; i < ngrp; i++) norm += xnew[i]*xnew[i];
                norm = sqrt(norm);
                if (norm > 0.) for (i = 0; i < ngrp; i++) xnew[i] /= norm;
                // convergence: allow sign flip
                diff1 = 0.; diff2 = 0.;
                for (i = 0; i < ngrp; i++) {
                    d1 = xnew[i] - x[i]; d2 = xnew[i] + x[i];
                    diff1 += d1*d1; diff2 += d2*d2;
                }
                memcpy(x, xnew, ngrp * sizeof(double));
                if (diff1 < 1e-12 || diff2 < 1e-12) break;
            }
            for (i = 0; i < ngrp; i++) { rords[i].which = i; rords[i].event = x[i]; }
            qsort(rords, ngrp, sizeof(ord_dbl_t), ord_dbl_acmpfunc);
            for (i = 0; i < ngrp; i++) { grp_order[i] = rords[i].which; grp_rank[grp_order[i]] = i; }

            free(deg); 
            free(x); 
            free(xnew);
        }

        // FM swap refinement using a priority queue keyed by group id.
        // Mirrors refine_haplotype_partition: for each group ga the PQ stores its best
        // swap partner (gb, delta) within SWAP_WIN rank positions.
        //
        // Structure of one FM pass (analogous to the fm_q pass in refine_haplotype_partition):
        //   1. Populate PQ: for every group ga compute best_swap_partner and insert.
        //   2. Pop ga with most-improving (most negative) delta.
        //      Stale check: if stored partner gb is already locked, recompute ga's best
        //      unlocked partner and re-insert (lazy invalidation — same trick as smark in
        //      refine_haplotype_partition).
        //   3. Apply swap(ga, gb): only two ranks change, so the delta formula is exact.
        //      Lock both; remove gb from PQ.
        //   4. Update PQ for all unlocked sadj-neighbours of ga and gb: their deltas
        //      changed because rank[ga] and rank[gb] changed (same as updating neighbours
        //      of a committed sequence in refine_haplotype_partition).
        //   5. Track cumulative cost delta; record best-prefix step.
        //   6. Roll back moves best_step..nmoves-1; if best_step>0 repeat.
        //
        // Cost per pass: O(ngrp/2) swaps × O(W·d) per swap (W=window, d=avg sadj degree)
        //   for neighbor updates; O(log n) PQ ops. Total O(n·W·d·log n) vs O(n²·W·d) before.
        {
            swap_prior_t _sp, *sp = &_sp;
            pq_elem_t fme;
            pq_t *swap_pq;
            int *locked, *move_ga, *move_gb, *upd;
            double *cum_delta;
            int fm_pass, fm_improved, max_iter, swap_win;
            int nmoves, best_step, upd_n, ne, dup, q, e;
            int ga, gb, ra, rb, km, rk, bg;
            double w, wp, bd, bdelta, total_d, best_total;

            max_iter = ngrp * 10;
            if (max_iter < FM_PASS_MAX)
                max_iter = FM_PASS_MAX;
            swap_win = ngrp / 10;
            if (swap_win < SWAP_WIN_MIN)
                swap_win = SWAP_WIN_MIN;
            if (swap_win > SWAP_WIN_MAX)
                swap_win = SWAP_WIN_MAX;

            MYCALLOC(locked,    ngrp);
            MYMALLOC(move_ga,   ngrp / 2 + 1);
            MYMALLOC(move_gb,   ngrp / 2 + 1);
            MYMALLOC(cum_delta, ngrp / 2 + 2);
            MYMALLOC(upd,       ngrp);  // scratch: unique unlocked neighbours to update

            swap_pq = pq_init(ngrp, sizeof(swap_prior_t),
                              swap_prior_cpyfunc, swap_prior_cmpfunc);

            fm_improved = 1;
            for (fm_pass = 1; fm_improved && fm_pass <= max_iter; fm_pass++) {
                fm_improved = 0;
                MYBZERO(locked, ngrp);

                // reduce window 5% every 10 passes, floor at SWAP_WIN_MIN
                if (fm_pass % 10 == 0) {
                    swap_win = swap_win * 95 / 100;
                    if (swap_win < SWAP_WIN_MIN)
                        swap_win = SWAP_WIN_MIN;
                }
                
                // populate: compute best window-swap partner for every group
                for (i = 0; i < ngrp; i++) {
                    bg = best_swap_partner(i, grp_rank, grp_order, ngrp,
                                          sadj, NULL, swap_win, &bd);
                    if (bg >= 0) {
                        sp->delta = bd;
                        sp->partner = bg;
                        pq_insert(swap_pq, i, sp);
                    }
                }

                // FM pass
                nmoves = 0;
                total_d = 0.; best_total = 0.;
                best_step = 0;
                cum_delta[0] = 0.;

                while (pq_pop1(swap_pq, &fme)) {
                    ga = fme.id;
                    gb = ((swap_prior_t *)fme.priority)->partner;
                    bdelta = ((swap_prior_t *)fme.priority)->delta;

                    // stale check: partner locked since this entry was inserted
                    if (locked[gb]) {
                        // recompute ga's best unlocked partner and re-insert
                        bg = best_swap_partner(ga, grp_rank, grp_order, ngrp,
                                               sadj, locked, swap_win, &bd);
                        if (bg >= 0) {
                            sp->delta = bd;
                            sp->partner = bg;
                            pq_insert(swap_pq, ga, sp); // will be re-popped later
                        }
                        continue;
                    }

                    // apply swap: only rank[ga] and rank[gb] change — delta formula exact
                    // Recompute bdelta from current ranks: the stored PQ value may be stale
                    // if intermediate swaps changed ranks of groups in sadj[gb] (those rank
                    // changes affect ga's delta but ga is not a sadj-neighbour of those groups,
                    // so its PQ entry was never refreshed by the neighbour-update step).
                    ra = grp_rank[ga];
                    rb = grp_rank[gb];
                    bdelta = 0.;
                    for (e = 0; e < (int)sadj[ga].n; e++) {
                        km = (int)sadj[ga].a[e].g;
                        if (km == gb) continue;
                        rk = grp_rank[km];
                        w = sadj[ga].a[e].w;
                        wp = sadj[ga].a[e].wp;
                        bdelta += directed_cost(w, wp, rb, rk) - directed_cost(w, wp, ra, rk);
                    }
                    for (e = 0; e < (int)sadj[gb].n; e++) {
                        km = (int)sadj[gb].a[e].g; if (km == ga) continue;
                        rk = grp_rank[km];
                        w = sadj[gb].a[e].w;
                        wp = sadj[gb].a[e].wp;
                        bdelta -= directed_cost(w, wp, rb, rk) - directed_cost(w, wp, ra, rk);
                    }
                    // direct {ga,gb} edge: cost changes when both endpoints swap
                    for (e = 0; e < (int)sadj[ga].n; e++) {
                        if ((int)sadj[ga].a[e].g == gb) {
                            w = sadj[ga].a[e].w;
                            wp = sadj[ga].a[e].wp;
                            bdelta += directed_cost(w, wp, rb, ra) - directed_cost(w, wp, ra, rb);
                            break;
                        }
                    }
                    
                    grp_rank[ga] = rb;
                    grp_rank[gb] = ra;
                    grp_order[ra] = gb; 
                    grp_order[rb] = ga;
                    locked[ga] = locked[gb] = 1;
                    pq_remove(swap_pq, gb); // gb locked, remove from PQ
                    move_ga[nmoves] = ga; 
                    move_gb[nmoves] = gb;
                    total_d += bdelta;
                    cum_delta[++nmoves] = total_d;
                    if (total_d < best_total) { 
                        best_total = total_d; 
                        best_step = nmoves;
                    }

                    // update PQ for unlocked sadj-neighbours of ga and gb.
                    // Collect unique unlocked neighbours first to avoid double-processing.
                    upd_n = 0;
                    for (e = 0; e < (int)sadj[ga].n; e++) {
                        ne = (int)sadj[ga].a[e].g;
                        if (!locked[ne]) upd[upd_n++] = ne;
                    }
                    for (e = 0; e < (int)sadj[gb].n; e++) {
                        ne = (int)sadj[gb].a[e].g;
                        if (locked[ne]) continue;
                        dup = 0;
                        for (q = 0; q < upd_n; q++) {
                            if (upd[q] == ne) { 
                                dup = 1; 
                                break; 
                            }
                        }
                        if (!dup) upd[upd_n++] = ne;
                    }
                    for (q = 0; q < upd_n; q++) {
                        ne = upd[q];
                        bg = best_swap_partner(ne, grp_rank, grp_order, ngrp,
                                               sadj, locked, swap_win, &bd);
                        if (bg >= 0) {
                            sp->delta = bd;
                            sp->partner = bg;
                            pq_insert(swap_pq, ne, sp); // insert-or-update
                        } else {
                            pq_remove(swap_pq, ne); // no valid partner any more
                        }
                    }
                }

                // roll back moves from nmoves-1 down to best_step (in reverse)
                for (k = nmoves - 1; k >= best_step; k--) {
                    ga = move_ga[k]; 
                    gb = move_gb[k];
                    ra = grp_rank[ga]; 
                    rb = grp_rank[gb];
                    grp_rank[ga] = rb; 
                    grp_rank[gb] = ra;
                    grp_order[ra] = gb; 
                    grp_order[rb] = ga;
                }

                if (best_step > 0)
                    fm_improved = 1;
            }

            // apply final permutation
            for (i = 0; i < ngrp; i++) new_blks[i] = blks.a[grp_order[i]];
            memcpy(blks.a, new_blks, ngrp * sizeof(scf_block_t));

            for (i = 0; i < nseq; i++) {
                if (!grps[i]) continue;
                ga = abs(grps[i]) - 1;
                gb = grp_rank[ga] + 1;
                grps[i] = grps[i] > 0 ? gb : -gb;
            }
            for (i = 0; i < nseq; i++)
                if (haps[i].grp) haps[i].grp = abs(grps[i]);

            free(locked); 
            free(move_ga); 
            free(move_gb); 
            free(cum_delta); 
            free(upd);
            pq_destroy(swap_pq);
        }

        for (i = 0; i < ngrp; i++)
            kv_destroy(sadj[i]);
        free(sadj);
        free(grp_order);
        free(grp_rank);
        free(new_blks);
        free(rords);
    }

    // find mergeable group according to LCS blocks and do merge
    if (with_hap) {
        hap_info_t *hseq, **hseqs, *ha, *hb;
        kvec_t(hap_ovl_pair_t) pairs;
        kvec_t(grp_pair_t) gpairs;
        hap_ovl_pair_t pr;
        ovl_t *o;
        int64 *hap_bpos, *hap_epos, *bbpos, pa, pb, pc, ea, eb, ca, cb, ext_a, ext_b;
        double v, *dp, best;
        int *used_a, *used_b, *chain, chain_n, tmp;
        int a, b, ga, gb, ia, ib, ic, ia2, ib2, cx, prev;
        int na, nb, np, hi, hj, m, x, y, t, *tb, *hcnt, *unions, *skips;
        
        MYCALLOC(hap_bpos, ploidy);
        MYCALLOC(hap_epos, ploidy);
        MYMALLOC(hseqs, ploidy);
        MYMALLOC(hseqs[0], nseq);
        MYMALLOC(hcnt,  ploidy);
        MYBZERO(hcnt, ploidy);
        for (i = 0; i < nseq; i++)
            hcnt[haps[i].hap-1]++;
        for (k = 1; k < ploidy; k++)
            hseqs[k] = hseqs[k-1] + hcnt[k-1];
        hseq = hseqs[0];
        memcpy(hseq, haps, nseq * sizeof(hap_info_t));
        for (i = 0; i < (int)blks.n; i++) {
            blk = blks.a + i;
            MYBZERO(hap_epos, ploidy);
            for (j = 0; j < blk->nseq; j++) {
                s = blk->seqs[j]>>1;
                k = haps[s].hap-1;
                if (haps[s].epos > hap_epos[k])
                    hap_epos[k] = haps[s].epos;
                hseq[s].grp   = 0; // disemble group
                hseq[s].bpos += hap_bpos[k];
                hseq[s].epos += hap_bpos[k];
            }
            for (k = 0; k < ploidy; k++)
                hap_bpos[k] += hap_epos[k];
        }
        // this is sorted by hap and then bpos
        qsort(hseq, nseq, sizeof(hap_info_t), hap_info_acmpfunc);
        
        for (k = 0; k < ploidy; k++) {
            n = hcnt[k];
            slen = 0;
            for (i = 0; i < n; i++) {
                hseqs[k][i].bpos = slen;
                hseqs[k][i].epos = slen + hseqs[k][i].len;
                slen += hseqs[k][i].len;
            }
        }

        // for each pair of haplotypes: find max-weight chain of overlapping sequence
        // pairs (i,j): position-based chain DP on actual overlap backbone coordinates
        MYMALLOC(bbpos, nseq);
        MYMALLOC(dp, nseq);
        MYMALLOC(tb, nseq);
        MYMALLOC(used_a, nseq);
        MYMALLOC(used_b, nseq);
        MYMALLOC(chain, nseq);
        kv_init(pairs);
        kv_init(gpairs);
        kv_resize(hap_ovl_pair_t, pairs, nseq);
        m = nseq;
        hseq = hseqs[0];
        for (i = 0; i < nseq; i++)
            bbpos[hseq[i].seq] = hseq[i].bpos;
        for (hi = 1; hi <= ploidy; hi++) {
            na = hcnt[hi-1];
            ha = hseqs[hi-1];
            for (hj = hi + 1; hj <= ploidy; hj++) {
                nb = hcnt[hj-1];
                hb = hseqs[hj-1];
                pairs.n = 0;
                for (p = 0; p < na; p++) {
                    for (q = 0; q < nb; q++) {
                        o = find_b_overlap(ovls + (index[ha[p].seq] >> 32),
                                        (uint32) index[ha[p].seq], hb[q].seq);
                        if (o && !o->del) {
                            pr.p = p; pr.q = q; pr.ovl = o;
                            pr.ab = haps[ha[p].seq].rev ? bbpos[ha[p].seq] + dicts->s[ha[p].seq].len - o->aepos
                                : bbpos[ha[p].seq] + o->abpos;
                            pr.ae = haps[ha[p].seq].rev ? bbpos[ha[p].seq] + dicts->s[ha[p].seq].len - o->abpos
                                : bbpos[ha[p].seq] + o->aepos;
                            pr.bb = haps[hb[q].seq].rev ? bbpos[hb[q].seq] + dicts->s[hb[q].seq].len - o->bepos
                                : bbpos[hb[q].seq] + o->bbpos;
                            pr.be = haps[hb[q].seq].rev ? bbpos[hb[q].seq] + dicts->s[hb[q].seq].len - o->bbpos
                                : bbpos[hb[q].seq] + o->bepos;
                            kv_push(hap_ovl_pair_t, pairs, pr);
                        }
                    }
                }
                qsort(pairs.a, pairs.n, sizeof(hap_ovl_pair_t), hap_ovl_pair_cmpfunc);
                if (m < (int)pairs.n) {
                    m = pairs.n << 1;
                    MYREALLOC(dp, m);
                    MYREALLOC(tb, m);
                    MYREALLOC(chain, m+1);
                }
                // O(m^2) chain DP
                best = 0.; t = -1;
                for (x = 0; x < (int)pairs.n; x++) {
                    ea = pairs.a[x].ae - pairs.a[x].ab;
                    eb = pairs.a[x].be - pairs.a[x].bb;
                    dp[x] = ea < eb ? ea : eb;
                    tb[x] = -1;
                    for (y = 0; y < x; y++) {
                        ext_a = pairs.a[x].ae
                            - (pairs.a[y].ae > pairs.a[x].ab ? pairs.a[y].ae : pairs.a[x].ab);
                        ext_b = pairs.a[x].be
                            - (pairs.a[y].be > pairs.a[x].bb ? pairs.a[y].be : pairs.a[x].bb);
                        if (ext_a > 0 && ext_b > 0) {
                            v = dp[y] + (ext_a < ext_b ? ext_a : ext_b);
                            if (v > dp[x]) { dp[x] = v; tb[x] = y; }
                        }
                    }
                    if (dp[x] > best) { best = dp[x]; t = x; }
                }

                // collect chain in forward order
                tmp = t;
                chain_n = 0;
                while (tmp >= 0) { 
                    chain[chain_n++] = tmp; 
                    tmp = tb[tmp];
                }
                for (x = 0, y = chain_n-1; x < y; x++, y--) { 
                    tmp = chain[x]; 
                    chain[x] = chain[y]; 
                    chain[y] = tmp;
                }

                // mark which ha/hb indices are used in the chain
                MYBZERO(used_a, na);
                MYBZERO(used_b, nb);
                for (x = 0; x < chain_n; x++) {
                    used_a[pairs.a[chain[x]].p] = 1;
                    used_b[pairs.a[chain[x]].q] = 1;
                }

                // merge: ia scans ha[], ib scans hb[], ic scans chain[]
                // advance the pointer whose next event comes first
                ia = ib = ic = 0;
                while (ia < na || ib < nb || ic < chain_n) {
                    // find next event positions
                    pa = pb = pc = INT64_MAX;
                    // next unused ha
                    ia2 = ia;
                    while (ia2 < na && used_a[ia2]) ia2++;
                    if (ia2 < na) pa = bbpos[ha[ia2].seq];
                    // next unused hb
                    ib2 = ib;
                    while (ib2 < nb && used_b[ib2]) ib2++;
                    if (ib2 < nb) pb = bbpos[hb[ib2].seq];
                    // next chain entry (keyed by min of its two positions)
                    if (ic < chain_n) {
                        ca = pairs.a[chain[ic]].ab;
                        cb = pairs.a[chain[ic]].bb;
                        pc = ca < cb ? ca : cb;
                    }
                    if (pc <= pa && pc <= pb) {
                        // print paired chain entry
                        cx = chain[ic++];
                        prev = tb[cx];
                        a = ha[pairs.a[cx].p].seq;
                        b = hb[pairs.a[cx].q].seq;
                        ga = haps[a].grp;
                        gb = haps[b].grp;
                        if (prev < 0) {
                            ext_a = pairs.a[cx].ae - pairs.a[cx].ab;
                            ext_b = pairs.a[cx].be - pairs.a[cx].bb;
                        } else {
                            ext_a = pairs.a[cx].ae
                                - (pairs.a[prev].ae > pairs.a[cx].ab ? pairs.a[prev].ae : pairs.a[cx].ab);
                            ext_b = pairs.a[cx].be
                                - (pairs.a[prev].be > pairs.a[cx].bb ? pairs.a[prev].be : pairs.a[cx].bb);
                        }
                        // collect inter-group sequence pairs along the chain
                        ext_a = MIN(ext_a, ext_b);
                        o = find_b_overlap(ovls + (index[a] >> 32), (uint32) index[a], b);
                        if (ga < gb) 
                            kv_push(grp_pair_t, gpairs, ((grp_pair_t){ga, gb, 1, ext_a, ext_a * o->qual}));
                        else if (gb < ga)
                            kv_push(grp_pair_t, gpairs, ((grp_pair_t){gb, ga, 1, ext_a, ext_a * o->qual}));
                        // advance ia/ib past the consumed indices
                        while (ia <= pairs.a[cx].p) ia++;
                        while (ib <= pairs.a[cx].q) ib++;
                    } else if (pa <= pb) {
                        ia = ia2;
                        ia++;
                    } else {
                        ib = ib2;
                        ib++;
                    }
                }
            }
        }
        if (gpairs.n > 0) {
            // sort groups and merge pairs
            qsort(gpairs.a, gpairs.n, sizeof(grp_pair_t), grppair_g_cmpfunc);
            
            ga = gpairs.a[0].a; 
            gb = gpairs.a[0].b;
            np = 0;
            n  = gpairs.n;
            for (i = 1; i < n; i++) {
                if (gpairs.a[i].a == ga && gpairs.a[i].b == gb) {
                    gpairs.a[np].n += gpairs.a[i].n;
                    gpairs.a[np].l += gpairs.a[i].l;
                    gpairs.a[np].w += gpairs.a[i].w;
                } else {
                    ga = gpairs.a[i].a;
                    gb = gpairs.a[i].b;
                    np++;
                    if (np < i)
                        gpairs.a[np] = gpairs.a[i];
                }
            }
            np++;

            // do actual merging
            MYMALLOC(unions, ngrp);
            for (i = 0; i < ngrp; i++)
                unions[i] = -1;
            for (i = 0; i < np; i++) {
                ga = gpairs.a[i].a - 1;
                gb = gpairs.a[i].b - 1;
                slen = blks.a[ga].slen;
                if (slen > blks.a[gb].slen)
                    slen = blks.a[gb].slen;
                for (j = ga+1; j < gb; j++)
                    slen += blks.a[j].slen;
                if (slen * 0.05 > gpairs.a[i].w)
                    continue;
                for (j = ga+1; j <= gb; j++)
                    unions[j] = ga;
            }

            // union groups and assign final group ids
            for (i = 0; i < ngrp; i++) {
                if (unions[i] < 0)
                    continue;
                k = i;
                while (unions[k] >= 0)
                    k = unions[k];
                j = i;
                while (j != k) {
                    unions[j] = k;
                    j = unions[j];
                }
            }

            // update groups
            MYBZERO(sords, ngrp);
            for (i = 0; i < ngrp; i++) {
                k = unions[i] >= 0? unions[i] : i;
                sords[k].which += 1;
                sords[k].event += blks.a[i].nseq;
            }
            // expand memory for sequences in merged groups
            for (i = 0; i < ngrp; i++) {
                if (unions[i] >= 0)
                    continue;
                // number sequences
                n = sords[i].event;
                blk = blks.a + i;
                MYREALLOC(blk->seqs, n);
                MYREALLOC(blk->type, n);
            }
            for (i = 0; i < ngrp; i++) {
                if (unions[i] < 0)
                    continue;
                blk = blks.a + unions[i];
                memcpy(blk->seqs + blk->nseq, blks.a[i].seqs, blks.a[i].nseq * sizeof(uint32));
                memcpy(blk->type + blk->nseq, blks.a[i].type, blks.a[i].nseq * sizeof(uint32));
                blk->nseq += blks.a[i].nseq;
                blk->slen += blks.a[i].slen;
            }

            // reorder sequences in merged groups (skip singletons that weren't actually merged)
            MYMALLOC(skips, ngrp);
            for (i = 0; i < ngrp; i++)
                skips[i] = (sords[i].which <= 1);
            mst_order_groups(blks.a, ngrp, ovls, index, haps, dicts, nseq, skips);
            free(skips);

            // compact blocks and free merged blocks
            for (i = j = 0; i < ngrp; i++) {
                if (unions[i] >= 0) {
                    blk = blks.a + i;
                    free(blk->seqs);
                    free(blk->type);
                    continue;
                }
                if (i != j)
                    blks.a[j] = blks.a[i];
                j++;
            }
            ngrp = blks.n = j;
            // update grps[] and haps[].grp/rev
            for (i = 0; i < ngrp; i++) {
                blk = blks.a + i;
                for (j = 0; j < blk->nseq; j++) {
                    s = blk->seqs[j]>>1;
                    haps[s].grp = i+1;
                    haps[s].rev = blk->seqs[j] & 1;
                    grps[s] = haps[s].rev? -(i+1) : (i+1);
                }
            }

            free(unions);
        }

        free(dp); 
        free(tb); 
        free(bbpos); 
        free(hap_bpos);
        free(hap_epos);
        free(hseqs[0]);
        free(hseqs);
        free(hcnt);
        kv_destroy(pairs); 
        kv_destroy(gpairs);
    }

    // print the final group order and orientation
    if (with_hap && VERBOSE) {
        hap_info_t *ha, *hb;
        range_t *arngs, *brngs;
        kvec_t(hap_ovl_pair_t) pairs;
        hap_ovl_pair_t pr;
        ovl_t *o;
        hap_info_t *hseq, **hseqs;
        int64 *hap_bpos, *hap_epos, *bbpos, ea, eb, ext_a, ext_b;
        int na, nb, ca, cb, hi, hj, m, x, y, *hcnt;
        double *dp, v, best;

        MYCALLOC(hap_bpos, ploidy);
        MYCALLOC(hap_epos, ploidy);
        MYMALLOC(hseqs, ploidy);
        MYMALLOC(hseqs[0], nseq);
        MYMALLOC(hcnt,  ploidy);
        MYBZERO(hcnt, ploidy);
        for (i = 0; i < nseq; i++)
            hcnt[haps[i].hap-1]++;
        for (k = 1; k < ploidy; k++)
            hseqs[k] = hseqs[k-1] + hcnt[k-1];
        hseq = hseqs[0];
        memcpy(hseq, haps, nseq * sizeof(hap_info_t));
        for (i = 0; i < (int)blks.n; i++) {
            blk = blks.a + i;
            MYBZERO(hap_epos, ploidy);
            for (j = 0; j < blk->nseq; j++) {
                s = blk->seqs[j]>>1;
                k = haps[s].hap-1;
                if (haps[s].epos > hap_epos[k])
                    hap_epos[k] = haps[s].epos;
                hseq[s].grp   = 0; // disemble group
                hseq[s].bpos += hap_bpos[k];
                hseq[s].epos += hap_bpos[k];
            }
            for (k = 0; k < ploidy; k++)
                hap_bpos[k] += hap_epos[k];
        }
        // this is sorted by hap and then bpos
        qsort(hseq, nseq, sizeof(hap_info_t), hap_info_acmpfunc);
        
        for (k = 0; k < ploidy; k++) {
            n = hcnt[k];
            slen = 0;
            for (i = 0; i < n; i++) {
                hseqs[k][i].bpos = slen;
                hseqs[k][i].epos = slen + hseqs[k][i].len;
                slen += hseqs[k][i].len;
            }
        }
        
        MYMALLOC(bbpos, nseq);
        MYMALLOC(dp, nseq);
        MYMALLOC(arngs, novl);
        MYMALLOC(brngs, novl);
        kv_init(pairs);
        kv_resize(hap_ovl_pair_t, pairs, nseq);
        m = nseq;
        hseq = hseqs[0];
        for (i = 0; i < nseq; i++)
            bbpos[hseq[i].seq] = hseq[i].bpos;

        fprintf(stderr, "[M::%s] haplotype collinear blocks\n", __func__);
        for (hi = 1; hi <= ploidy; hi++) {
            na = hcnt[hi-1];
            ha = hseqs[hi-1];
            for (hj = hi + 1; hj <= ploidy; hj++) {
                nb = hcnt[hj-1];
                hb = hseqs[hj-1];
                // calculate overlapped bases
                ca = cb = 0;
                for (i = 0; i < novl; i++) {
                    ovl = ovls + i;
                    if (ovl->del) continue;
                    x = ovl->aread;
                    y = ovl->bread;
                    if (haps[x].hap != hi || haps[y].hap != hj)
                        continue;
                    if (haps[x].rev)
                        arngs[ca++] = (range_t){bbpos[x]+dicts->s[x].len-ovl->aepos, bbpos[x]+dicts->s[x].len-ovl->abpos};
                    else
                        arngs[ca++] = (range_t){bbpos[x]+ovl->abpos, bbpos[x]+ovl->aepos};
                    if (haps[y].rev)
                        brngs[cb++] = (range_t){bbpos[y]+dicts->s[y].len-ovl->bepos, bbpos[y]+dicts->s[y].len-ovl->bbpos};
                    else
                        brngs[cb++] = (range_t){bbpos[y]+ovl->bbpos, bbpos[y]+ovl->bepos};
                }
                ca = rangelist_size(arngs, ca, 0);
                cb = rangelist_size(brngs, cb, 0);
                // collect all overlapping pairs and store their hap positions
                pairs.n = 0;
                for (p = 0; p < na; p++) {
                    for (q = 0; q < nb; q++) {
                        o = find_b_overlap(ovls + (index[ha[p].seq] >> 32),
                                        (uint32) index[ha[p].seq], hb[q].seq);
                        if (o && !o->del) {
                            pr.p = p; pr.q = q; pr.ovl = o;
                            pr.ab = ha[p].rev ? bbpos[ha[p].seq] + dicts->s[ha[p].seq].len - o->aepos
                                : bbpos[ha[p].seq] + o->abpos;
                            pr.ae = ha[p].rev ? bbpos[ha[p].seq] + dicts->s[ha[p].seq].len - o->abpos
                                : bbpos[ha[p].seq] + o->aepos;
                            pr.bb = hb[q].rev ? bbpos[hb[q].seq] + dicts->s[hb[q].seq].len - o->bepos
                                : bbpos[hb[q].seq] + o->bbpos;
                            pr.be = hb[q].rev ? bbpos[hb[q].seq] + dicts->s[hb[q].seq].len - o->bbpos
                                : bbpos[hb[q].seq] + o->bepos;
                            kv_push(hap_ovl_pair_t, pairs, pr);
                        }
                    }
                }
                // sort by (ae, be) for correct topological DP order
                qsort(pairs.a, pairs.n, sizeof(hap_ovl_pair_t), hap_ovl_pair_cmpfunc);
                if (m < (int)pairs.n) {
                    m = pairs.n << 1;
                    MYREALLOC(dp, m);
                }
                // O(m^2) chain DP
                best = 0.;
                for (x = 0; x < (int)pairs.n; x++) {
                    ea = pairs.a[x].ae - pairs.a[x].ab;
                    eb = pairs.a[x].be - pairs.a[x].bb;
                    dp[x] = ea < eb ? ea : eb;
                    for (y = 0; y < x; y++) {
                        ext_a = pairs.a[x].ae
                            - (pairs.a[y].ae > pairs.a[x].ab ? pairs.a[y].ae : pairs.a[x].ab);
                        ext_b = pairs.a[x].be
                            - (pairs.a[y].be > pairs.a[x].bb ? pairs.a[y].be : pairs.a[x].bb);
                        if (ext_a > 0 && ext_b > 0) {
                            v = dp[y] + (ext_a < ext_b ? ext_a : ext_b);
                            if (v > dp[x]) dp[x] = v;
                        }
                    }
                    if (dp[x] > best) best = dp[x];
                }
                fprintf(stderr, "[M::%s] LCS overlap H%d-H%d: %10.0f @ %10d %10d\n", __func__, hi, hj, best, ca, cb);
            }
        }
        free(bbpos);
        free(dp);
        free(arngs); 
        free(brngs);
        free(hap_bpos);
        free(hap_epos);
        free(hseqs[0]);
        free(hseqs);
        free(hcnt);
        kv_destroy(pairs);
    }

    // fix bepos so [bpos, bepos) has exactly the sequence length
    for (i = 0; i < nseq; i++)
        haps[i].epos = haps[i].bpos + dicts->s[i].len;

    for (i = 0; i < blks.n; i++) {
        free(blks.a[i].seqs);
        free(blks.a[i].type);
    }
    free(blks.a);
    free(covs[0]);
    free(covs);
    free(ncov);
    free(grps);
    free(rngs);
    free(cpts);
    free(index);
    free(sords);
    free(qords);
    free(cseqs);
    free(ctype);
    free(rmarks);
    free(qmarks);
    kv_destroy(clus);
    kv_destroy(exts);
}

static void print_partition_stats(hap_info_t *seqs, int nseq, cft_t *cfts, uint64 *index, int ploidy)
{
    cft_v_t *hcfts;
    int *hseqs;
    uint8 h;
    cft_t *cft;
    int64 i, c, n, comp, size, *hsize;
    
    // allocate memory
    MYCALLOC(hcfts, ploidy+1);
    MYCALLOC(hseqs, ploidy+1);
    MYCALLOC(hsize, ploidy+1);
    if (hcfts == NULL || hseqs == NULL || hsize == NULL)
        mem_alloc_error("haplotype partition arrays");
    size = 0;
    comp = 0;
    for (i = 0; i < nseq; i++) {
        h = seqs[i].hap;
        cft = cfts + (index[i] >> 32);
        n = (uint32) index[i];
        for (c = 0; c < n; c++, cft++)
            if (seqs[cft->b].hap == h)
                hcfts[h].vals[cft->t] += cft->v;
        size += seqs[i].len;
        comp += (int64) seqs[i].len * seqs[i].hap;
        hseqs[h] += 1;
        hsize[h] += seqs[i].len;
    }
    cft_v_reset(hcfts);
    for (c = 1; c <= ploidy; c++)
        cft_v_add(hcfts, hcfts + c);
    fprintf(stderr, "[M::%s] haplotype partition for %d sequences of %lld bases:\n", __func__, nseq, size);
    for (c = 1; c <= ploidy; c++)
        fprintf(stderr, "[M::%s] [H%lld] %6d %12lld [%12lld %12lld]\n", __func__, 
            c, hseqs[c], hsize[c], hcfts[c].vals[CFT_OVL]/2, hcfts[c].vals[CFT_HIC]/2);
    fprintf(stderr, "[M::%s] total conflitcs: %lld %lld [%lld]\n", __func__, hcfts[0].vals[CFT_OVL]/2, hcfts[0].vals[CFT_HIC]/2, comp);
    
    free(hcfts);
    free(hseqs);
    free(hsize);
}


static void par_prior_cpyfunc(void *x, void *y)
{
    *(int64 *) x = *(int64 *) y;
}

static int par_prior_cmpfunc(void *x, void *y)
{
    int64 a = *(int64 *) x;
    int64 b = *(int64 *) y;
    return (a < b) - (a > b);
}

static int u64_cmpfunc(const void *x, const void *y)
{
    uint64 a = *(uint64 *) x;
    uint64 b = *(uint64 *) y;
    return (a > b) - (a < b);
}

static void build_haplotype_partition_ovl_greedy(cft_t *cfts, uint64 *index, hap_info_t *seqs, int nseq, int ploidy)
{
    ord_i64_t *sords;
    cft_v_t *hcfts, *hcft;
    cft_t *cft;
    uint8 *haps, *chap;
    int64 i, n, size, *hsize;
    int a, c, p, h;
    double covs, *scopy;

    // allocate memory
    MYMALLOC(haps, nseq);
    MYMALLOC(chap, nseq);
    MYMALLOC(sords, nseq);
    MYMALLOC(scopy, nseq);
    MYMALLOC(hcfts, ploidy+1);
    MYMALLOC(hsize, ploidy+1);
    if (haps == NULL || chap == NULL || sords == NULL || 
        scopy == NULL || hcfts == NULL || hsize == NULL)
        mem_alloc_error("haplotype partition arrays");

    // build sequence orders by length
    for (i = 0; i < nseq; i++)
        sords[i] = (ord_i64_t) {i, seqs[i].len};
    qsort(sords, nseq, sizeof(ord_i64_t), ord_i64_dcmpfunc);

    // calculate average copy number of each sequence
    for (i = 0; i < nseq; i++) {
        cft = cfts + (index[i] >> 32);
        n = (uint32) index[i];
        covs = .0;
        for (c = 0; c < n; c++, cft++)
            if (cft->t == CFT_OVL) // sequence overlaps
                covs += cft->v;
        scopy[i] = covs / seqs[i].len;
    }

    // build haplotype partition
    // by greedy graph colouring
    MYBZERO(haps, nseq);
    MYBZERO(hsize, ploidy+1);
    for (i = 0; i < nseq; i++) {
        a = sords[i].which;
        if (haps[a])
            continue;

        // calculate conflicts of the sequence to each haplotype group
        MYBZERO(hcfts, ploidy+1);
        cft = cfts + (index[a] >> 32);
        n = (uint32) index[a];
        cft_v_all(cft, n, haps, hcfts);

        // assign the sequence to the group with least conflicts
        h = 1;
        hcft = hcfts + h;
        size = hsize[h];
        for (p = 2; p <= ploidy; p++) {
            c = cft_v_cmpfunc(hcft, hcfts + p);
            if (c > 0 || 
                (c == 0 && size > hsize[p] && scopy[a] > p-1)) {
                h = p;
                hcft = hcfts + h;
                size = hsize[h];
            }
        }

        haps[a] = h;
        hsize[h] += seqs[a].len; // count the total length of sequences
    }

    // update haplotype partition for each sequence
    for (i = 0; i < nseq; i++)
        seqs[i].hap = haps[i]? haps[i] : 1;

#ifdef DEBUG_HAPLOTYPE_PARTITION
    // some statistics
    fprintf(stderr, "[M::%s] initial global haplotype partition\n", __func__);
    print_partition_stats(seqs, nseq, cfts, index, ploidy);
#endif

    free(haps);
    free(chap);
    free(sords);
    free(scopy);
    free(hcfts);
    free(hsize);
}

#define LOUVAIN_RNG_SEED 42U
#define LOUVAIN_DEFAULT_GRANULARITY 0.5
#define LOUVAIN_HIGH_GRANULARITY 1

static int build_hic_linkage_group_louvain(hlk_t *hlks, int64 nhlk, sdict_t *dicts, int *reps, double granularity, int *comms)
{
    int    *adj_nb_A,    *adj_nb_B,    *adj_nb_cur,    *adj_nb_nxt,    *adj_nb_tmp;
    double *adj_w_A,     *adj_w_B,     *adj_w_cur,     *adj_w_nxt,     *adj_w_tmp;
    int64  *adj_start_A, *adj_start_B, *adj_start_cur, *adj_start_nxt, *adj_start_tmp;
    double *ki_A,        *ki_B,        *ki_cur,        *ki_nxt,        *ki_tmp;
    double *Sigma;
    int    *comm, *node_order, *neigh_comms, *remap, *orig_comm;
    int    *node_list, *comm_start, *cnt_arr, *grp_remap;
    double *neigh_weight;
    int    *has_link;
    int64  i, j, e, pos, nedge_cur, nedge_nxt;
    double m2, ki_comm, score_cur, score, best_score;
    int    u, cv, c, q, ngrp, cur_c, best_c, changed, pass, max_pass, ncomm, level_improved;
    int    n_cur, n_nxt, n_neigh;
    int    nseq;
    unsigned rng;

    // initialise
    ngrp = 0;
    if (nhlk == 0)
        return 0;
    
    nseq = dicts->n;

    MYMALLOC(adj_nb_A,     nhlk);
    MYMALLOC(adj_nb_B,     nhlk);
    MYMALLOC(adj_w_A,      nhlk);
    MYMALLOC(adj_w_B,      nhlk);
    MYCALLOC(adj_start_A,  nseq + 1);
    MYCALLOC(adj_start_B,  nseq + 1);
    MYCALLOC(ki_A,         nseq);
    MYCALLOC(ki_B,         nseq);
    MYCALLOC(Sigma,        nseq);
    MYCALLOC(comm,         nseq);
    MYMALLOC(node_order,   nseq);
    MYCALLOC(neigh_weight, nseq);
    MYMALLOC(neigh_comms,  nseq);
    MYMALLOC(remap,        nseq);
    MYMALLOC(orig_comm,    nseq);
    MYMALLOC(node_list,    nseq);
    MYCALLOC(comm_start,   nseq + 1);
    MYCALLOC(cnt_arr,      nseq);

    for (i = 0; i < nhlk; i++) {
        if (reps && 
            (reps[hlks[i].a] || reps[hlks[i].b])) 
            continue;
        cnt_arr[hlks[i].a]++;
    }
    adj_start_A[0] = 0;
    for (i = 0; i < nseq; i++) 
        adj_start_A[i+1] = adj_start_A[i] + cnt_arr[i];
    MYBZERO(cnt_arr, nseq);
    for (i = 0; i < nhlk; i++) {
        if (reps && 
            (reps[hlks[i].a] || reps[hlks[i].b])) 
            continue;
        u   = hlks[i].a;
        pos = adj_start_A[u] + cnt_arr[u]++;
        adj_nb_A[pos] = hlks[i].b;
        adj_w_A[pos]  = hlks[i].v;
    }
    MYBZERO(cnt_arr, nseq);

    nedge_cur = nhlk;
    n_cur     = nseq;
    m2        = 0.;
    for (i = 0; i < n_cur; i++) {
        ki_A[i] = 0.;
        for (e = adj_start_A[i]; e < adj_start_A[i+1]; e++)
            ki_A[i] += adj_w_A[e];
        m2 += ki_A[i];
    }

    adj_nb_cur    = adj_nb_A;    adj_nb_nxt    = adj_nb_B;
    adj_w_cur     = adj_w_A;     adj_w_nxt     = adj_w_B;
    adj_start_cur = adj_start_A; adj_start_nxt = adj_start_B;
    ki_cur        = ki_A;        ki_nxt        = ki_B;

    for (i = 0; i < nseq; i++) 
        orig_comm[i] = (int)i;
    rng = LOUVAIN_RNG_SEED;

    MYCALLOC(has_link, nseq);
    for (i = 0; i < nseq; i++) 
        has_link[i] = (ki_A[i] > 0.) ? 1 : 0;

    if (m2 <= 0.) goto free_all;

    for (;;) {
        for (i = 0; i < n_cur; i++) {
            comm[i]       = (int)i;
            Sigma[i]      = ki_cur[i];
            node_order[i] = (int)i;
        }
        max_pass = n_cur * 20 + 100;
        level_improved = 0;

        for (pass = 0; pass < max_pass; pass++) {
            /* Fisher-Yates shuffle (Knuth LCG, fixed seed) */
            for (i = (int64)n_cur - 1; i > 0; i--) {
                rng = rng * 1664525U + 1013904223U;
                j = (int64)(rng % (unsigned)(i + 1));
                SWAP(int, node_order[i], node_order[j]);
            }

            changed = 0;
            for (q = 0; q < n_cur; q++) {
                u = node_order[q];
                cur_c = comm[u];
                Sigma[cur_c] -= ki_cur[u];

                n_neigh = 0;
                ki_comm = 0.;
                for (e = adj_start_cur[u]; e < adj_start_cur[u+1]; e++) {
                    c = comm[adj_nb_cur[e]];
                    if (c == cur_c) {
                        ki_comm += adj_w_cur[e];
                    } else {
                        if (neigh_weight[c] == 0.) 
                            neigh_comms[n_neigh++] = c;
                        neigh_weight[c] += adj_w_cur[e];
                    }
                }

                score_cur  = ki_comm - granularity * ki_cur[u] * Sigma[cur_c] / m2;
                best_c     = cur_c;
                best_score = score_cur;
                for (i = 0; i < n_neigh; i++) {
                    c     = neigh_comms[i];
                    score = neigh_weight[c] - granularity * ki_cur[u] * Sigma[c] / m2;
                    if (score > best_score) { 
                        best_score = score; 
                        best_c = c;
                    }
                    neigh_weight[c] = 0.;
                }

                comm[u] = best_c;
                Sigma[best_c] += ki_cur[u];
                if (best_c != cur_c) { 
                    changed = 1; 
                    level_improved = 1;
                }
            }
            if (!changed) break;
        }

        if (!level_improved) break;

        for (i = 0; i < n_cur; i++)
            remap[i] = -1;
        ncomm = 0;
        for (i = 0; i < n_cur; i++)
            if (remap[comm[i]] < 0)
                remap[comm[i]] = ncomm++;

        for (i = 0; i < nseq; i++)
            orig_comm[i] = remap[comm[orig_comm[i]]];

        if (ncomm >= n_cur) break;

        n_nxt = ncomm;

        for (i = 0; i < n_cur; i++) 
            cnt_arr[remap[comm[i]]]++;
        comm_start[0] = 0;
        for (i = 0; i < n_nxt; i++) 
            comm_start[i+1] = comm_start[i] + cnt_arr[i];
        MYBZERO(cnt_arr, n_nxt);
        for (i = 0; i < n_cur; i++) {
            c = remap[comm[i]];
            node_list[comm_start[c] + cnt_arr[c]++] = (int)i;
        }
        MYBZERO(cnt_arr, n_nxt);

        MYBZERO(ki_nxt, n_nxt);
        for (i = 0; i < n_cur; i++)
            ki_nxt[remap[comm[i]]] += ki_cur[i];

        nedge_nxt = 0;
        for (c = 0; c < n_nxt; c++) {
            n_neigh = 0;
            for (i = comm_start[c]; i < comm_start[c+1]; i++) {
                u = node_list[i];
                for (e = adj_start_cur[u]; e < adj_start_cur[u+1]; e++) {
                    cv = remap[comm[adj_nb_cur[e]]];
                    if (cv == c) continue;
                    if (neigh_weight[cv] == 0.) 
                        neigh_comms[n_neigh++] = cv;
                    neigh_weight[cv] += adj_w_cur[e];
                }
            }
            cnt_arr[c] = n_neigh;
            nedge_nxt += n_neigh;
            for (i = 0; i < n_neigh; i++) 
                neigh_weight[neigh_comms[i]] = 0.;
        }
        adj_start_nxt[0] = 0;
        for (i = 0; i < n_nxt; i++) 
            adj_start_nxt[i+1] = adj_start_nxt[i] + cnt_arr[i];
        
        MYBZERO(cnt_arr, n_nxt);
        for (c = 0; c < n_nxt; c++) {
            n_neigh = 0;
            for (i = comm_start[c]; i < comm_start[c+1]; i++) {
                u = node_list[i];
                for (e = adj_start_cur[u]; e < adj_start_cur[u+1]; e++) {
                    cv = remap[comm[adj_nb_cur[e]]];
                    if (cv == c) continue;
                    if (neigh_weight[cv] == 0.) 
                        neigh_comms[n_neigh++] = cv;
                    neigh_weight[cv] += adj_w_cur[e];
                }
            }
            for (i = 0; i < n_neigh; i++) {
                cv  = neigh_comms[i];
                pos = adj_start_nxt[c] + cnt_arr[c]++;
                adj_nb_nxt[pos] = cv;
                adj_w_nxt[pos]  = neigh_weight[cv];
                neigh_weight[cv] = 0.;
            }
        }
        MYBZERO(cnt_arr,    n_nxt);
        MYBZERO(comm_start, n_nxt + 1);

        adj_nb_tmp    = adj_nb_cur;    adj_nb_cur    = adj_nb_nxt;    adj_nb_nxt    = adj_nb_tmp;
        adj_w_tmp     = adj_w_cur;     adj_w_cur     = adj_w_nxt;     adj_w_nxt     = adj_w_tmp;
        adj_start_tmp = adj_start_cur; adj_start_cur = adj_start_nxt; adj_start_nxt = adj_start_tmp;
        ki_tmp        = ki_cur;        ki_cur        = ki_nxt;        ki_nxt        = ki_tmp;
        n_cur         = n_nxt;
        nedge_cur     = nedge_nxt;
        (void)nedge_cur;
    }

    MYMALLOC(grp_remap, nseq);
    for (i = 0; i < nseq; i++) grp_remap[i] = -1;
    for (i = 0; i < nseq; i++) {
        if (!has_link[i]) continue;
        if (grp_remap[orig_comm[i]] < 0) grp_remap[orig_comm[i]] = ++ngrp;
        comms[i] = grp_remap[orig_comm[i]];
    }
    free(grp_remap);

#ifdef DEBUG_HIC_LINKAGE_GROUP
    static const char *PALETTE[12] = {
        "#e6194b", "#3cb44b", "#4363d8", "#f58231",
        "#911eb4", "#42d4f4", "#f032e6", "#bfef45",
        "#fabed4", "#469990", "#dcbeff", "#9a6324"
    };
    int a, b, grp;

    fprintf(stdout, "%c <?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", PG_GML);
    fprintf(stdout, "%c <graphml xmlns=\"http://graphml.graphdrawing.org/graphml\"\n", PG_GML);
    fprintf(stdout, "%c          xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n", PG_GML);
    fprintf(stdout, "%c          xsi:schemaLocation=\"http://graphml.graphdrawing.org/graphml "
            "http://graphml.graphdrawing.org/graphml-1.0.xsd\">\n", PG_GML);
    fprintf(stdout, "%c   <key id=\"label\"   for=\"node\" attr.name=\"label\"   attr.type=\"string\"/>\n", PG_GML);
    fprintf(stdout, "%c   <key id=\"color\"   for=\"node\" attr.name=\"color\"   attr.type=\"string\"/>\n", PG_GML);
    fprintf(stdout, "%c   <key id=\"cluster\" for=\"node\" attr.name=\"cluster\" attr.type=\"int\"/>\n", PG_GML);
    fprintf(stdout, "%c   <key id=\"weight\"  for=\"edge\" attr.name=\"weight\"  attr.type=\"double\"/>\n", PG_GML);
    fprintf(stdout, "%c   <graph id=\"G\" edgedefault=\"undirected\">\n", PG_GML);
    for (i = 0; i < nseq; i++) {
        if (!has_link[i]) continue;
        grp = comms[i];
        fprintf(stdout, "%c     <node id=\"n%d\">\n", PG_GML, (int)i);
        fprintf(stdout, "%c       <data key=\"label\">%s</data>\n",   PG_GML, dicts->s[i].name);
        fprintf(stdout, "%c       <data key=\"color\">%s</data>\n",   PG_GML, (grp > 0) ? PALETTE[(grp - 1) % 12] : "#aaaaaa");
        fprintf(stdout, "%c       <data key=\"cluster\">%d</data>\n", PG_GML, grp);
        fprintf(stdout, "%c     </node>\n", PG_GML);
    }
    for (j = 0; j < nhlk; j++) {
        a = hlks[j].a; 
        b = hlks[j].b;
        if (a >= b || (reps && (reps[a] || reps[b]))) 
            continue;
        fprintf(stdout, "%c     <edge source=\"n%d\" target=\"n%d\">\n", PG_GML, a, b);
        fprintf(stdout, "%c       <data key=\"weight\">%g</data>\n", PG_GML, hlks[j].v);
        fprintf(stdout, "%c     </edge>\n", PG_GML);
    }
    fprintf(stdout, "%c   </graph>\n", PG_GML);
    fprintf(stdout, "%c </graphml>\n", PG_GML);
#endif

free_all:
    free(has_link);
    free(adj_nb_A);    free(adj_nb_B);
    free(adj_w_A);     free(adj_w_B);
    free(adj_start_A); free(adj_start_B);
    free(ki_A);        free(ki_B);
    free(Sigma);
    free(comm);        free(node_order);
    free(neigh_weight); free(neigh_comms);
    free(remap);       free(orig_comm);
    free(node_list);   free(comm_start);  free(cnt_arr);

    return ngrp;
}

#define LEIDEN_RNG_SEED 42U
#define LEIDEN_DEFAULT_GRANULARITY 0.5
#define LEIDEN_HIGH_GRANULARITY 1
#define LEIDEN_THETA 0.05   /* Traag et al. 2019 randomness parameter for refinement
                                merges: smaller -> greedier (closer to argmax),
                                larger -> more uniform/random among valid candidates */

static int build_hic_linkage_group_leiden(hlk_t *hlks, int64 nhlk, sdict_t *dicts, int *reps, double granularity, int *comms)
{
    int    *adj_nb_A,    *adj_nb_B,    *adj_nb_cur,    *adj_nb_nxt,    *adj_nb_tmp;
    double *adj_w_A,     *adj_w_B,     *adj_w_cur,     *adj_w_nxt,     *adj_w_tmp;
    int64  *adj_start_A, *adj_start_B, *adj_start_cur, *adj_start_nxt, *adj_start_tmp;
    double *ki_A,        *ki_B,        *ki_cur,        *ki_nxt,        *ki_tmp;
    double *Sigma;
    int    *comm, *node_order, *neigh_comms, *remap, *orig_comm;
    int    *node_list, *comm_start, *cnt_arr, *grp_remap;
    double *neigh_weight;
    int    *has_link;
    int64  i, j, e, pos, nedge_cur, nedge_nxt;
    double m2, ki_comm, score_cur, score, best_score;
    int    u, cv, c, q, ngrp, cur_c, best_c, changed, pass, max_pass, ncomm;
    int    n_cur, n_nxt, n_neigh;
    int    nseq;
    unsigned rng;

    /* ---- refinement-phase state (this is what makes it Leiden, not Louvain) --
       pcomm/pdeg      : the local-moving partition P, compacted, with K_S per P-community
       rcomm/rsize/rdeg/rext : the refined sub-partition being built inside P, with each
                                refined community's size, total degree and external weight
       rremap/r2p      : compaction of the refined partition + which P-community each
                                refined community belongs to
       prep/seed_*     : carries the non-refined partition P forward as the initial
                                (non-singleton) community assignment for the next level  --- */
    int    *pcomm, *rcomm, *rsize, *rremap, *r2p, *prep, *cand_id;
    double *pdeg, *rdeg, *rext, *cand_gain;
    int    *seed_A, *seed_B, *seed_cur, *seed_nxt, *seed_tmp;
    int64   fstart, fend;
    int     n_refined, n_cand, best_cand;
    double  k_s, e_vs, e_vc, gain, max_gain, tot_w, draw, cum;

    ngrp = 0;
    if (nhlk == 0)
        return 0;
    
    nseq = dicts->n;

    MYMALLOC(adj_nb_A,     nhlk);
    MYMALLOC(adj_nb_B,     nhlk);
    MYMALLOC(adj_w_A,      nhlk);
    MYMALLOC(adj_w_B,      nhlk);
    MYCALLOC(adj_start_A,  nseq + 1);
    MYCALLOC(adj_start_B,  nseq + 1);
    MYCALLOC(ki_A,         nseq);
    MYCALLOC(ki_B,         nseq);
    MYCALLOC(Sigma,        nseq);
    MYCALLOC(comm,         nseq);
    MYMALLOC(node_order,   nseq);
    MYCALLOC(neigh_weight, nseq);
    MYMALLOC(neigh_comms,  nseq);
    MYMALLOC(remap,        nseq);
    MYMALLOC(orig_comm,    nseq);
    MYMALLOC(node_list,    nseq);
    MYCALLOC(comm_start,   nseq + 1);
    MYCALLOC(cnt_arr,      nseq);

    MYMALLOC(pcomm,     nseq);
    MYCALLOC(pdeg,      nseq);
    MYMALLOC(rcomm,     nseq);
    MYMALLOC(rsize,     nseq);
    MYMALLOC(rremap,    nseq);
    MYMALLOC(r2p,       nseq);
    MYMALLOC(prep,      nseq);
    MYCALLOC(rdeg,      nseq);
    MYCALLOC(rext,      nseq);
    MYMALLOC(cand_id,   nseq);
    MYMALLOC(cand_gain, nseq);
    MYMALLOC(seed_A,    nseq);
    MYMALLOC(seed_B,    nseq);

    for (i = 0; i < nhlk; i++) {
        if (reps && 
            (reps[hlks[i].a] || reps[hlks[i].b])) 
            continue;
        cnt_arr[hlks[i].a]++;
    }
    adj_start_A[0] = 0;
    for (i = 0; i < nseq; i++) 
        adj_start_A[i+1] = adj_start_A[i] + cnt_arr[i];
    MYBZERO(cnt_arr, nseq);
    for (i = 0; i < nhlk; i++) {
        if (reps && 
            (reps[hlks[i].a] || reps[hlks[i].b])) 
            continue;
        u   = hlks[i].a;
        pos = adj_start_A[u] + cnt_arr[u]++;
        adj_nb_A[pos] = hlks[i].b;
        adj_w_A[pos]  = hlks[i].v;
    }
    MYBZERO(cnt_arr, nseq);

    nedge_cur = nhlk;
    n_cur     = nseq;
    m2        = 0.;
    for (i = 0; i < n_cur; i++) {
        ki_A[i] = 0.;
        for (e = adj_start_A[i]; e < adj_start_A[i+1]; e++)
            ki_A[i] += adj_w_A[e];
        m2 += ki_A[i];
    }

    adj_nb_cur    = adj_nb_A;    adj_nb_nxt    = adj_nb_B;
    adj_w_cur     = adj_w_A;     adj_w_nxt     = adj_w_B;
    adj_start_cur = adj_start_A; adj_start_nxt = adj_start_B;
    ki_cur        = ki_A;        ki_nxt        = ki_B;

    for (i = 0; i < nseq; i++) 
        orig_comm[i] = (int)i;
    rng = LEIDEN_RNG_SEED;

    MYCALLOC(has_link, nseq);
    for (i = 0; i < nseq; i++) 
        has_link[i] = (ki_A[i] > 0.) ? 1 : 0;

    if (m2 <= 0.) goto free_all;

    seed_cur = seed_A; seed_nxt = seed_B;
    for (i = 0; i < n_cur; i++)
        seed_cur[i] = (int)i;      /* level 0: no prior partition, start singleton */

    for (;;) {
        /* ---- seed this level's communities: singleton at level 0, otherwise the
                previous level's (non-refined) P-communities carried forward ------ */
        for (i = 0; i < n_cur; i++) {
            comm[i]       = seed_cur[i];
            node_order[i] = (int)i;
        }
        MYBZERO(Sigma, nseq);
        for (i = 0; i < n_cur; i++)
            Sigma[comm[i]] += ki_cur[i];

        max_pass = n_cur * 20 + 100;

        for (pass = 0; pass < max_pass; pass++) {
            /* Fisher-Yates shuffle (Knuth LCG, fixed seed) */
            for (i = (int64)n_cur - 1; i > 0; i--) {
                rng = rng * 1664525U + 1013904223U;
                j = (int64)(rng % (unsigned)(i + 1));
                SWAP(int, node_order[i], node_order[j]);
            }

            changed = 0;
            for (q = 0; q < n_cur; q++) {
                u = node_order[q];
                cur_c = comm[u];
                Sigma[cur_c] -= ki_cur[u];

                n_neigh = 0;
                ki_comm = 0.;
                for (e = adj_start_cur[u]; e < adj_start_cur[u+1]; e++) {
                    c = comm[adj_nb_cur[e]];
                    if (c == cur_c) {
                        ki_comm += adj_w_cur[e];
                    } else {
                        if (neigh_weight[c] == 0.) 
                            neigh_comms[n_neigh++] = c;
                        neigh_weight[c] += adj_w_cur[e];
                    }
                }

                score_cur  = ki_comm - granularity * ki_cur[u] * Sigma[cur_c] / m2;
                best_c     = cur_c;
                best_score = score_cur;
                for (i = 0; i < n_neigh; i++) {
                    c     = neigh_comms[i];
                    score = neigh_weight[c] - granularity * ki_cur[u] * Sigma[c] / m2;
                    if (score > best_score) { 
                        best_score = score; 
                        best_c = c;
                    }
                    neigh_weight[c] = 0.;
                }

                comm[u] = best_c;
                Sigma[best_c] += ki_cur[u];
                if (best_c != cur_c) 
                    changed = 1;
            }
            if (!changed) break;
        }

        /* ---- compact the local-moving partition P into 0..ncomm-1 ------------ */
        for (i = 0; i < n_cur; i++)
            remap[i] = -1;
        ncomm = 0;
        for (i = 0; i < n_cur; i++)
            if (remap[comm[i]] < 0)
                remap[comm[i]] = ncomm++;
        for (i = 0; i < n_cur; i++)
            pcomm[i] = remap[comm[i]];

        if (ncomm >= n_cur) {
            /* every node is its own P-community: refinement has nothing to find */
            for (i = 0; i < nseq; i++)
                orig_comm[i] = pcomm[orig_comm[i]];
            break;
        }

        /* K_S per P-community, and bucket current-level nodes by P-community */
        MYBZERO(pdeg, nseq);
        for (i = 0; i < n_cur; i++)
            pdeg[pcomm[i]] += ki_cur[i];

        MYBZERO(cnt_arr, nseq);
        for (i = 0; i < n_cur; i++) cnt_arr[pcomm[i]]++;
        comm_start[0] = 0;
        for (i = 0; i < ncomm; i++) comm_start[i+1] = comm_start[i] + cnt_arr[i];
        MYBZERO(cnt_arr, nseq);
        for (i = 0; i < n_cur; i++) {
            c = pcomm[i];
            node_list[comm_start[c] + cnt_arr[c]++] = (int)i;
        }
        MYBZERO(cnt_arr, nseq);

        /* ==== REFINEMENT ================================================
           Sub-partition each P-community with randomised local merging,
           screened by a well-connectedness threshold. Every community this
           produces induces a connected subgraph in the current-level graph,
           which by induction across levels means every FINAL community is
           connected in the original graph -- the guarantee plain Louvain
           local-moving does not provide (Traag, Waltman & van Eck, "From
           Louvain to Leiden: guaranteeing well-connected communities",
           Sci Rep 2019). ================================================= */
        for (i = 0; i < n_cur; i++) {
            rcomm[i] = (int)i;
            rsize[i] = 1;
            rdeg[i]  = ki_cur[i];
        }

        for (c = 0; c < ncomm; c++) {
            fstart = comm_start[c];
            fend   = comm_start[c+1];
            k_s    = pdeg[c];

            /* initial external weight of every (still-singleton) member w.r.t. S\{member} */
            for (i = fstart; i < fend; i++) {
                u = node_list[i];
                e_vs = 0.;
                for (e = adj_start_cur[u]; e < adj_start_cur[u+1]; e++)
                    if (pcomm[adj_nb_cur[e]] == c)
                        e_vs += adj_w_cur[e];
                rext[u] = e_vs;
            }

            /* shuffle this community's members (continues the same LCG stream) */
            for (i = fend - 1; i > fstart; i--) {
                rng = rng * 1664525U + 1013904223U;
                j = fstart + (int64)(rng % (unsigned)(i - fstart + 1));
                SWAP(int, node_list[i], node_list[j]);
            }

            for (q = fstart; q < fend; q++) {
                u = node_list[q];
                /* eligible only if u is still an untouched singleton: nobody has
                   joined it (rsize==1) AND it hasn't itself moved (rcomm==u) --
                   checking rsize alone is not enough, since a node that moved
                   away keeps a stale rsize of 1 */
                if (rcomm[u] != u || rsize[u] != 1) continue;

                n_neigh = 0;
                e_vs = 0.;
                for (e = adj_start_cur[u]; e < adj_start_cur[u+1]; e++) {
                    if (pcomm[adj_nb_cur[e]] != c) continue;
                    e_vs += adj_w_cur[e];
                    cv = rcomm[adj_nb_cur[e]];
                    if (neigh_weight[cv] == 0.)
                        neigh_comms[n_neigh++] = cv;
                    neigh_weight[cv] += adj_w_cur[e];
                }

                /* node-level well-connectedness screen: is u well connected to
                   the rest of its P-community? if not, it stays singleton */
                if (e_vs < granularity * ki_cur[u] * (k_s - ki_cur[u]) / m2) {
                    for (i = 0; i < n_neigh; i++) neigh_weight[neigh_comms[i]] = 0.;
                    continue;
                }

                /* candidates: well-connected sub-communities of S that u has a
                   direct edge to, plus "stay singleton" (gain 0, always valid) */
                n_cand = 0;
                cand_id[n_cand]   = u;
                cand_gain[n_cand] = 0.;
                n_cand++;
                max_gain = 0.;

                for (i = 0; i < n_neigh; i++) {
                    cv = neigh_comms[i];
                    if (rext[cv] >= granularity * rdeg[cv] * (k_s - rdeg[cv]) / m2) {
                        gain = neigh_weight[cv] - granularity * ki_cur[u] * rdeg[cv] / m2;
                        if (gain >= 0.) {
                            cand_id[n_cand]   = cv;
                            cand_gain[n_cand] = gain;
                            n_cand++;
                            if (gain > max_gain) max_gain = gain;
                        }
                    }
                }

                /* softmax sample among candidates (randomised merge, theta-controlled) */
                tot_w = 0.;
                for (i = 0; i < n_cand; i++) {
                    cand_gain[i] = exp((cand_gain[i] - max_gain) / LEIDEN_THETA);
                    tot_w += cand_gain[i];
                }
                rng  = rng * 1664525U + 1013904223U;
                draw = ((double)rng / 4294967296.0) * tot_w;
                cum  = 0.;
                best_cand = cand_id[n_cand - 1];
                for (i = 0; i < n_cand; i++) {
                    cum += cand_gain[i];
                    if (draw < cum) { best_cand = cand_id[i]; break; }
                }

                if (best_cand != u) {
                    e_vc = neigh_weight[best_cand];
                    rext[best_cand] += e_vs - 2. * e_vc;
                    rdeg[best_cand] += ki_cur[u];
                    rsize[best_cand] += 1;
                    rcomm[u] = best_cand;
                }

                for (i = 0; i < n_neigh; i++) neigh_weight[neigh_comms[i]] = 0.;
            }
        }

        /* ---- compact the refined partition into 0..n_refined-1 --------------- */
        for (i = 0; i < n_cur; i++) rremap[i] = -1;
        n_refined = 0;
        for (i = 0; i < n_cur; i++) {
            u = rcomm[i];
            if (rremap[u] < 0) rremap[u] = n_refined++;
        }
        for (i = 0; i < n_cur; i++)
            r2p[rremap[rcomm[i]]] = pcomm[i];

        for (i = 0; i < nseq; i++)
            orig_comm[i] = rremap[rcomm[orig_comm[i]]];

        if (n_refined >= n_cur) break;    /* refinement found nothing to merge: converged */

        /* ---- seed the NEXT level from P (not from refined singletons): each
                aggregate node starts in whichever P-community it came from,
                using one representative refined-id per P-community -------------- */
        for (i = 0; i < ncomm; i++) prep[i] = -1;
        for (i = 0; i < n_refined; i++)
            if (prep[r2p[i]] < 0) prep[r2p[i]] = i;
        for (i = 0; i < n_refined; i++)
            seed_nxt[i] = prep[r2p[i]];

        n_nxt = n_refined;

        /* ---- bucket current-level nodes by refined community, for aggregation - */
        MYBZERO(cnt_arr, nseq);
        for (i = 0; i < n_cur; i++) cnt_arr[rremap[rcomm[i]]]++;
        comm_start[0] = 0;
        for (i = 0; i < n_nxt; i++) comm_start[i+1] = comm_start[i] + cnt_arr[i];
        MYBZERO(cnt_arr, nseq);
        for (i = 0; i < n_cur; i++) {
            c = rremap[rcomm[i]];
            node_list[comm_start[c] + cnt_arr[c]++] = (int)i;
        }
        MYBZERO(cnt_arr, nseq);

        MYBZERO(ki_nxt, n_nxt);
        for (i = 0; i < n_cur; i++)
            ki_nxt[rremap[rcomm[i]]] += ki_cur[i];

        nedge_nxt = 0;
        for (c = 0; c < n_nxt; c++) {
            n_neigh = 0;
            for (i = comm_start[c]; i < comm_start[c+1]; i++) {
                u = node_list[i];
                for (e = adj_start_cur[u]; e < adj_start_cur[u+1]; e++) {
                    cv = rremap[rcomm[adj_nb_cur[e]]];
                    if (cv == c) continue;
                    if (neigh_weight[cv] == 0.) 
                        neigh_comms[n_neigh++] = cv;
                    neigh_weight[cv] += adj_w_cur[e];
                }
            }
            cnt_arr[c] = n_neigh;
            nedge_nxt += n_neigh;
            for (i = 0; i < n_neigh; i++) 
                neigh_weight[neigh_comms[i]] = 0.;
        }
        adj_start_nxt[0] = 0;
        for (i = 0; i < n_nxt; i++) 
            adj_start_nxt[i+1] = adj_start_nxt[i] + cnt_arr[i];

        MYBZERO(cnt_arr, nseq);
        for (c = 0; c < n_nxt; c++) {
            n_neigh = 0;
            for (i = comm_start[c]; i < comm_start[c+1]; i++) {
                u = node_list[i];
                for (e = adj_start_cur[u]; e < adj_start_cur[u+1]; e++) {
                    cv = rremap[rcomm[adj_nb_cur[e]]];
                    if (cv == c) continue;
                    if (neigh_weight[cv] == 0.) 
                        neigh_comms[n_neigh++] = cv;
                    neigh_weight[cv] += adj_w_cur[e];
                }
            }
            for (i = 0; i < n_neigh; i++) {
                cv  = neigh_comms[i];
                pos = adj_start_nxt[c] + cnt_arr[c]++;
                adj_nb_nxt[pos] = cv;
                adj_w_nxt[pos]  = neigh_weight[cv];
                neigh_weight[cv] = 0.;
            }
        }
        MYBZERO(cnt_arr,    nseq);
        MYBZERO(comm_start, nseq + 1);

        adj_nb_tmp    = adj_nb_cur;    adj_nb_cur    = adj_nb_nxt;    adj_nb_nxt    = adj_nb_tmp;
        adj_w_tmp     = adj_w_cur;     adj_w_cur     = adj_w_nxt;     adj_w_nxt     = adj_w_tmp;
        adj_start_tmp = adj_start_cur; adj_start_cur = adj_start_nxt; adj_start_nxt = adj_start_tmp;
        ki_tmp        = ki_cur;        ki_cur        = ki_nxt;        ki_nxt        = ki_tmp;
        seed_tmp      = seed_cur;      seed_cur      = seed_nxt;      seed_nxt      = seed_tmp;
        n_cur         = n_nxt;
        nedge_cur     = nedge_nxt;
        (void)nedge_cur;
    }

    MYMALLOC(grp_remap, nseq);
    for (i = 0; i < nseq; i++) grp_remap[i] = -1;
    for (i = 0; i < nseq; i++) {
        if (!has_link[i]) continue;
        if (grp_remap[orig_comm[i]] < 0) grp_remap[orig_comm[i]] = ++ngrp;
        comms[i] = grp_remap[orig_comm[i]];
    }
    free(grp_remap);

#ifdef DEBUG_HIC_LINKAGE_GROUP
    static const char *PALETTE[12] = {
        "#e6194b", "#3cb44b", "#4363d8", "#f58231",
        "#911eb4", "#42d4f4", "#f032e6", "#bfef45",
        "#fabed4", "#469990", "#dcbeff", "#9a6324"
    };
    int a, b, grp;

    fprintf(stdout, "%c <?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", PG_GML);
    fprintf(stdout, "%c <graphml xmlns=\"http://graphml.graphdrawing.org/graphml\"\n", PG_GML);
    fprintf(stdout, "%c          xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n", PG_GML);
    fprintf(stdout, "%c          xsi:schemaLocation=\"http://graphml.graphdrawing.org/graphml "
            "http://graphml.graphdrawing.org/graphml-1.0.xsd\">\n", PG_GML);
    fprintf(stdout, "%c   <key id=\"label\"   for=\"node\" attr.name=\"label\"   attr.type=\"string\"/>\n", PG_GML);
    fprintf(stdout, "%c   <key id=\"color\"   for=\"node\" attr.name=\"color\"   attr.type=\"string\"/>\n", PG_GML);
    fprintf(stdout, "%c   <key id=\"cluster\" for=\"node\" attr.name=\"cluster\" attr.type=\"int\"/>\n", PG_GML);
    fprintf(stdout, "%c   <key id=\"weight\"  for=\"edge\" attr.name=\"weight\"  attr.type=\"double\"/>\n", PG_GML);
    fprintf(stdout, "%c   <graph id=\"G\" edgedefault=\"undirected\">\n", PG_GML);
    for (i = 0; i < nseq; i++) {
        if (!has_link[i]) continue;
        grp = comms[i];
        fprintf(stdout, "%c     <node id=\"n%d\">\n", PG_GML, (int)i);
        fprintf(stdout, "%c       <data key=\"label\">%s</data>\n",   PG_GML, dicts->s[i].name);
        fprintf(stdout, "%c       <data key=\"color\">%s</data>\n",   PG_GML, (grp > 0) ? PALETTE[(grp - 1) % 12] : "#aaaaaa");
        fprintf(stdout, "%c       <data key=\"cluster\">%d</data>\n", PG_GML, grp);
        fprintf(stdout, "%c     </node>\n", PG_GML);
    }
    for (j = 0; j < nhlk; j++) {
        a = hlks[j].a; 
        b = hlks[j].b;
        if (a >= b || (reps && (reps[a] || reps[b]))) 
            continue;
        fprintf(stdout, "%c     <edge source=\"n%d\" target=\"n%d\">\n", PG_GML, a, b);
        fprintf(stdout, "%c       <data key=\"weight\">%g</data>\n", PG_GML, hlks[j].v);
        fprintf(stdout, "%c     </edge>\n", PG_GML);
    }
    fprintf(stdout, "%c   </graph>\n", PG_GML);
    fprintf(stdout, "%c </graphml>\n", PG_GML);
#endif

free_all:
    free(has_link);
    free(adj_nb_A);    free(adj_nb_B);
    free(adj_w_A);     free(adj_w_B);
    free(adj_start_A); free(adj_start_B);
    free(ki_A);        free(ki_B);
    free(Sigma);
    free(comm);        free(node_order);
    free(neigh_weight); free(neigh_comms);
    free(remap);       free(orig_comm);
    free(node_list);   free(comm_start);  free(cnt_arr);

    free(pcomm);   free(pdeg);
    free(rcomm);   free(rsize); free(rremap); free(r2p); free(prep);
    free(rdeg);    free(rext);
    free(cand_id); free(cand_gain);
    free(seed_A);  free(seed_B);

    return ngrp;
}

static inline uint8 find_candidate_move(cft_t *cfts, uint64 *index, int seq, uint8 *haps, 
    int ploidy, cft_v_t *hap_conflicts, cft_v_t *delta_conflicts)
{
    cft_v_t *best_conflicts, *hcfts = hap_conflicts + 1;
    cft_t *cft = cfts + (index[seq] >> 32);
    int p, h, bh, n = (uint32) index[seq];
    h = haps[seq];
    MYBZERO(hcfts, ploidy+1);
    cft_v_all(cft, n, haps, hcfts);
    bh = 1 + (h==1);
    best_conflicts = hcfts + bh;
    for (p = bh + 1; p < h; p++) {
        if (cft_v_cmpfunc(hcfts + p, best_conflicts) < 0) {
            bh = p;
            best_conflicts = hcfts + bh;
        }
    }
    for (p = h + 1; p <= ploidy; p++) {
        if (cft_v_cmpfunc(hcfts + p, best_conflicts) < 0) {
            bh = p;
            best_conflicts = hcfts + bh;
        }
    }
    *hap_conflicts = hcfts[h];
    *delta_conflicts = *best_conflicts;
    cft_v_sub(delta_conflicts, hap_conflicts);
    return bh;
}

static inline void update_haplotype(uint8 *haps, hs_opts_t *hs_opts)
{
    int i = hs_opts->n;
    hs_opt_t *opts = hs_opts->a;
    while (i--) haps[opts->seq] = opts->hap, opts++;
}

static void refine_haplotype_partition_ovl_vns(cft_t *cfts, int64 ncft, uint64 *index, hap_info_t *seqs, int nseq, 
    int ploidy, vns_data_t *vns_data)
{
    cft_v_t *best_conflicts, *delta_conflicts, *hcfts, *scfts;
    cft_v_t _tcfts, *tcfts = &_tcfts;
    cft_v_t _bcfts, *bcfts = &_bcfts;
    pq_t *fm_q, *fs_q;
    cft_prior_t _fmp, *fmp = &_fmp;
    pq_elem_t _fme, *fme = &_fme;
    hs_opts_t *hs_opts;
    cft_t *cft;
    uint8 h, *haps, *bhap;
    int64 iter, max_iter;
    int64 tcomp, bcomp;
    int i, j, a, b, c, n, bstep, *smark;

#ifdef DEBUG_REFINE_HAP_PARTITION
    int64 comp_score;
#endif
    
    // collect data
    haps = vns_data->haps;
    bhap = vns_data->bhap;
    fm_q = vns_data->fm_q;
    fs_q = vns_data->fs_q;
    smark = vns_data->smark;
    hs_opts = &vns_data->hs_opts;
    best_conflicts = vns_data->hcfts;
    delta_conflicts = vns_data->dcfts;

    // set max iteration
    // this is arbitrary, but should be enough for convergence
    // bounded by the hyperparameter
    max_iter = (int64) nseq * 3; // * ploidy * 100;
    if (max_iter > vns_data->max_iter)
        max_iter = vns_data->max_iter;

    // some precomputes
    hcfts = best_conflicts + 1;
    scfts = delta_conflicts + nseq;

#ifdef DEBUG_REFINE_HAP_PARTITION
    comp_score = 0;
    for (i = 0; i < nseq; i++)
        comp_score += (int64) seqs[i].len * seqs[i].hap;
#endif

    // record the best solution so far
    for (i = 0; i < nseq; i++)
        bhap[i] = seqs[i].hap;
    cft_v_reset(best_conflicts);
    cft_v_total(cfts, ncft, bhap, best_conflicts);
    
    // add candidate moves
    for (i = 0; i < nseq; i++) {
        h = find_candidate_move(cfts, index, i, bhap, ploidy, hcfts, delta_conflicts+i);
        *fmp = (cft_prior_t) {*hcfts, seqs[i].len, h};
        if (!cft_v_zero(&fmp->cfs) || h < bhap[i]) pq_insert(fs_q, i, fmp);
    }

    // do local search
    iter = 0;
    while (iter++ < max_iter && pq_pop1(fs_q, fme)) {
        // do a fm_pass
        pq_insert(fm_q, fme->id, fme->priority);
        *tcfts = *best_conflicts;
        *bcfts = *best_conflicts;
        scfts[fme->id] = delta_conflicts[fme->id];
        tcomp = 0;
        bcomp = 0;
        bstep = 0;
        hs_opts->n = 0;
        MYBZERO(smark, nseq);
        memcpy(haps, bhap, nseq);
        while (pq_pop1(fm_q, fme)) {
            a = fme->id;
            h = ((cft_prior_t *) (fme->priority))->hap;
            kv_push(hs_opt_t, *hs_opts, ((hs_opt_t) {a, h}));
            tcomp += (int64) seqs[a].len * ((int) h - haps[a]);
            cft_v_add(tcfts, scfts + a);
            if ((c = cft_v_cmpfunc(tcfts, bcfts)) < 0 || 
                (c == 0 && tcomp < bcomp)) {
                *bcfts = *tcfts;
                bcomp = tcomp;
                bstep = hs_opts->n;
            }

            // update candidate and lock
            haps[a] = h;
            smark[a] = 1;

            // update moves of neighbors
            cft = cfts + (index[a] >> 32);
            n = (uint32) index[a];
            for (i = 0; i < n; i++, cft++) {
                b = cft->b;
                if (smark[b])
                    continue;
                h = find_candidate_move(cfts, index, b, haps, ploidy, hcfts, scfts+b);
                *fmp = (cft_prior_t) {*hcfts, seqs[b].len, h};
                if (!cft_v_zero(&fmp->cfs) || h < haps[b])
                    pq_insert(fm_q, b, fmp);
                else pq_remove(fm_q, b);
            }
        }

        if ((c = cft_v_cmpfunc(bcfts, best_conflicts)) < 0 ||
            (c == 0 && bcomp < 0)) {
            *best_conflicts = *bcfts;
            hs_opts->n = bstep; // the best solution in this pass
            update_haplotype(bhap, hs_opts);

#ifdef DEBUG_REFINE_HAP_PARTITION
            comp_score += bcomp;
            fprintf(stderr, "[M::%s] [iter %lld] found better solution with conflicts [%lld %lld] [%lld]\n",
                __func__, iter, best_conflicts->vals[CFT_OVL], best_conflicts->vals[CFT_HIC], comp_score);
#endif

            // also update the candidate queue
            // sequences in the hs_opts and their neighbours
            MYBZERO(smark, nseq);
            for (i = 0; i < bstep; i++) {
                a = hs_opts->a[i].seq;
                h = find_candidate_move(cfts, index, a, bhap, ploidy, hcfts, delta_conflicts+a);
                *fmp = (cft_prior_t) {*hcfts, seqs[a].len, h};
                if (!cft_v_zero(&fmp->cfs) || h < bhap[a])
                    pq_insert(fs_q, a, fmp);
                else pq_remove(fs_q, a);
                smark[a] = 1;
                cft = cfts + (index[a] >> 32);
                n = (uint32) index[a];
                for (j = 0; j < n; j++, cft++) {
                    b = cft->b;
                    if (smark[b])
                        continue;
                    h = find_candidate_move(cfts, index, b, bhap, ploidy, hcfts, delta_conflicts+b);
                    *fmp = (cft_prior_t) {*hcfts, seqs[b].len, h};
                    if (!cft_v_zero(&fmp->cfs) || h < bhap[b])
                        pq_insert(fs_q, b, fmp);
                    else pq_remove(fs_q, b);
                    smark[b] = 1;
                }
            }
        }
    }

#ifdef DEBUG_REFINE_HAP_PARTITION
    fprintf(stderr, "[M::%s] [iter %lld] final solution with conflicts [%lld %lld] [%lld]\n",
        __func__, iter, best_conflicts->vals[CFT_OVL], best_conflicts->vals[CFT_HIC], comp_score);
#endif

    // record the best solution so far
    for (i = 0; i < nseq; i++)
        seqs[i].hap = bhap[i];

#ifdef DEBUG_HAPLOTYPE_PARTITION
    // some statistics
    fprintf(stderr, "[M::%s] haplotype partition after refinement\n", __func__);
    print_partition_stats(seqs, nseq, cfts, index, ploidy);
#endif
}

typedef struct { int a, b; int64 w; } hedge_t;

static int hedge_w_dcmpfunc(const void *x, const void *y)
{
    const hedge_t *a = (const hedge_t *) x;
    const hedge_t *b = (const hedge_t *) y;
    return (a->w < b->w) - (a->w > b->w);
}

static void refine_haplotype_partition_hic_2opt(cft_t *cfts, int64 ncft, uint64 *index, hap_info_t *seqs, int nseq, 
    int ploidy, vns_data_t *vns_data)
{
    kvec_t(hedge_t) hedges;
    kvec_t(int) queue;
    kvec_t(hs_opt_t) best_opts;
    cft_v_t best_cfts, tcfts, bcfts, bcfts_best;
    hedge_t *he;
    cft_t *cft;
    uint8 *haps, *bhap, *visited;
    hs_opts_t *hs_opts;
    int64 i, j, n, hic_w;
    int a, b, ha, hb, improved, d, bstep, bstep_best;
    int has_ovl, s, h_from, h_to, nb, qh, k;

#ifdef DEBUG_REFINE_HAP_PARTITION
    int npass, ncascade;
#endif

    haps    = vns_data->haps;
    bhap    = vns_data->bhap;
    hs_opts = &vns_data->hs_opts;

    for (i = 0; i < nseq; i++)
        bhap[i] = seqs[i].hap;

    // OVL floor — never exceed this
    cft_v_reset(&best_cfts);
    cft_v_total(cfts, (int) ncft, bhap, &best_cfts);

    MYMALLOC(visited, nseq);
    if (visited == NULL)
        mem_alloc_error("hic 2opt visited array");
    kv_init(hedges);
    kv_init(queue);
    kv_init(best_opts);

#ifdef DEBUG_REFINE_HAP_PARTITION
    npass = 0;
#endif
    do {
        improved = 0;
#ifdef DEBUG_REFINE_HAP_PARTITION
        ncascade = 0;
#endif

        // collect inter-haplotype pure-HIC edges
        hedges.n = 0;
        for (a = 0; a < nseq; a++) {
            cft = cfts + (index[a] >> 32);
            n   = (uint32) index[a];
            i   = 0;
            while (i < n) {
                b = cft[i].b;
                if (b <= a) { i++; continue; }
                has_ovl = 0;
                hic_w   = 0;
                j = i;
                while (j < n && cft[j].b == b) {
                    if (cft[j].t == CFT_OVL) has_ovl = 1;
                    else                     hic_w  += cft[j].v;
                    j++;
                }
                if (!has_ovl && hic_w > 0 && bhap[a] != bhap[b])
                    kv_push(hedge_t, hedges, ((hedge_t){a, b, hic_w}));
                i = j;
            }
        }

        if (!hedges.n) break;
        qsort(hedges.a, hedges.n, sizeof(hedge_t), hedge_w_dcmpfunc);

        // for each edge
        for (he = hedges.a; he < hedges.a + hedges.n; he++) {
            a  = he->a;
            b  = he->b;
            ha = bhap[a];
            hb = bhap[b];
            if (ha == hb) continue;  // stale: merged by a previous cascade this pass

            // try both seed directions; keep the best prefix across both
            // d=0: seed a -> hb;  d=1: seed b -> ha
            bcfts_best = best_cfts;  // sentinel — no improvement across either direction
            bstep_best = 0;
            best_opts.n = 0;

            for (d = 0; d < 2; d++) {
                s = (d == 0) ? a : b;   // seed sequence

                memcpy(haps, bhap, nseq);
                MYBZERO(visited, nseq);
                hs_opts->n = 0;
                tcfts = best_cfts;
                bcfts = best_cfts;
                bstep = 0;

                visited[s] = 1;
                kv_push(int, queue, s);
                qh = 0;

                while (qh < (int) queue.n) {
                    s      = queue.a[qh++];
                    h_from = haps[s];
                    h_to   = (h_from == ha) ? hb : ha;

                    // incremental conflict update — call before updating haps[s]
                    cft = cfts + (index[s] >> 32);
                    n   = (uint32) index[s];
                    cft_v_move(cft, (int) n, (uint8) h_from, (uint8) h_to, haps, &tcfts);
                    haps[s] = h_to;

                    kv_push(hs_opt_t, *hs_opts, ((hs_opt_t){s, (uint8) h_to}));

                    // best-prefix: record step if OVL floor is respected and strictly better
                    if (tcfts.vals[CFT_OVL] <= best_cfts.vals[CFT_OVL] &&
                        cft_v_cmpfunc(&tcfts, &bcfts) < 0) {
                        bcfts = tcfts;
                        bstep = hs_opts->n;
                    }

                    // enqueue OVL neighbours now in h_to (they now conflict with s)
                    for (i = 0; i < n; i++) {
                        if (cft[i].t != CFT_OVL) continue;
                        nb = cft[i].b;
                        if (visited[nb])      continue;
                        if (haps[nb] != h_to) continue;
                        visited[nb] = 1;
                        kv_push(int, queue, nb);
                    }
                }
                queue.n = 0;

                // update cross-direction best
                if (bstep > 0 && cft_v_cmpfunc(&bcfts, &bcfts_best) < 0) {
                    bcfts_best = bcfts;
                    bstep_best = bstep;
                    kv_resize(hs_opt_t, best_opts, bstep);
                    memcpy(best_opts.a, hs_opts->a, bstep * sizeof(hs_opt_t));
                    best_opts.n = bstep;
                }
            }

            // commit the best prefix across both directions
            if (bstep_best > 0 && cft_v_cmpfunc(&bcfts_best, &best_cfts) < 0) {
                for (k = 0; k < bstep_best; k++)
                    bhap[best_opts.a[k].seq] = best_opts.a[k].hap;
                best_cfts = bcfts_best;
                improved  = 1;
#ifdef DEBUG_REFINE_HAP_PARTITION
                ncascade++;
#endif
            }
        }
    
#ifdef DEBUG_REFINE_HAP_PARTITION
        npass++;
        fprintf(stderr, "[M::%s] pass %d: %d cascades committed, conflicts [%lld %lld]\n",
            __func__, npass, ncascade, best_cfts.vals[CFT_OVL], best_cfts.vals[CFT_HIC]);
#endif

    } while (improved);

#ifdef DEBUG_REFINE_HAP_PARTITION
    fprintf(stderr, "[M::%s] hic-refine done after %d passes, final conflicts [%lld %lld]\n",
        __func__, npass, best_cfts.vals[CFT_OVL], best_cfts.vals[CFT_HIC]);
#endif

    for (i = 0; i < nseq; i++)
        seqs[i].hap = bhap[i];


#ifdef DEBUG_HAPLOTYPE_PARTITION
    fprintf(stderr, "[M::%s] haplotype partition after HIC 2-opt refinement\n", __func__);
    print_partition_stats(seqs, nseq, cfts, index, ploidy);
#endif

    free(visited);
    kv_destroy(hedges);
    kv_destroy(queue);
    kv_destroy(best_opts);
}

#define MAX_PATH_LEN 10

typedef struct {
    int s[MAX_PATH_LEN];
    double v;
} hs_path_t;

static int hs_path_v_dcmpfunc(const void *a, const void *b)
{
    const hs_path_t *pa = (const hs_path_t *) a;
    const hs_path_t *pb = (const hs_path_t *) b;
    if (pb->v > pa->v) return  1;
    if (pb->v < pa->v) return -1;
    return 0;
}

static inline int _uf_find_hic(int *uf, int x)
{
    while (uf[x] != x) {
        uf[x] = uf[uf[x]];
        x = uf[x];
    }
    return x;
}

static void refine_haplotype_partition_hic_kopt(cft_t *cfts, int64 ncft, uint64 *index,
    hap_info_t *seqs, int nseq, int ploidy, vns_data_t *vns_data)
{
    kvec_t(hedge_t) all_edges;
    kvec_t(hedge_t) comp_edges;
    kvec_t(int)     adj_head;
    kvec_t(int)     adj_next;
    kvec_t(int)     adj_to;
    kvec_t(int64)   adj_w;
    kvec_t(int)     path;
    kvec_t(int)     dfs_stack;
    kvec_t(int8)    dfs_depth;
    kvec_t(int)     queue;
    cft_v_t best_cfts, tcfts, bcfts;
    cft_t *cft;
    hs_path_t hp;
    kvec_t(hs_path_t) paths;
    uint8 *haps, *bhap, *visited, *in_path;
    int   *uf, *uf2, *cls_head, *cls_next, *cls_tail, *cls_size;
    int   *hap_cnt;
    hs_opts_t *hs_opts;
    int64 i, j, n, hic_w;
    int   a, b, ra, rb, s, ps, qs, nb, qh, k, m, e, ei, cur, depth, tmp, nxt, pi;
    int   has_ovl, h_from, h_to, h_maj, h_a, n_minority, pushed;
    int   bstep, improved, path_len;
    double pw;

#ifdef DEBUG_REFINE_HAP_PARTITION
    int npass, ncascade;
#endif

    haps    = vns_data->haps;
    bhap    = vns_data->bhap;
    hs_opts = &vns_data->hs_opts;

    for (i = 0; i < nseq; i++)
        bhap[i] = seqs[i].hap;

    cft_v_reset(&best_cfts);
    cft_v_total(cfts, (int) ncft, bhap, &best_cfts);

    MYMALLOC(uf,        nseq);
    MYMALLOC(cls_head,  nseq);
    MYMALLOC(cls_tail,  nseq);
    MYMALLOC(cls_next,  nseq);
    MYMALLOC(cls_size,  nseq);
    MYMALLOC(visited,   nseq);
    MYMALLOC(in_path,   nseq);
    MYMALLOC(hap_cnt,   ploidy + 1);
    if (!uf || !cls_head || !cls_tail || !cls_next || !cls_size ||
        !visited || !in_path || !hap_cnt)
        mem_alloc_error("hic kopt arrays");

    kv_init(all_edges);
    kv_init(comp_edges);
    kv_init(adj_head);
    kv_init(adj_next);
    kv_init(adj_to);
    kv_init(adj_w);
    kv_init(path);
    kv_init(dfs_stack);
    kv_init(dfs_depth);
    kv_init(queue);
    kv_init(paths);

#ifdef DEBUG_REFINE_HAP_PARTITION
    npass = 0;
#endif

    path_len = MAX_PATH_LEN;
    while (path_len >= 2) {

        improved = 0;
        paths.n = 0;
#ifdef DEBUG_REFINE_HAP_PARTITION
        ncascade = 0;
#endif

        // collect pure-HIC edges and build connected components
        for (i = 0; i < nseq; i++) {
            uf[i]       = i;
            cls_head[i] = i;
            cls_tail[i] = i;
            cls_next[i] = -1;
            cls_size[i] = 1;
        }

        all_edges.n = 0;
        for (a = 0; a < nseq; a++) {
            cft = cfts + (index[a] >> 32);
            n   = (uint32) index[a];
            i   = 0;
            while (i < n) {
                b = cft[i].b;
                if (b <= a) { i++; continue; }
                has_ovl = 0; hic_w = 0;
                j = i;
                while (j < n && cft[j].b == b) {
                    if (cft[j].t == CFT_OVL) has_ovl = 1;
                    else                     hic_w  += cft[j].v;
                    j++;
                }
                if (!has_ovl && hic_w > 0)
                    kv_push(hedge_t, all_edges, ((hedge_t){a, b, hic_w}));
                i = j;
            }
        }
        qsort(all_edges.a, all_edges.n, sizeof(hedge_t), hedge_w_dcmpfunc);

        for (i = 0; i < (int64) all_edges.n; i++) {
            a  = all_edges.a[i].a;
            b  = all_edges.a[i].b;
            ra = _uf_find_hic(uf, a);
            rb = _uf_find_hic(uf, b);
            if (ra == rb) continue;
            if (cls_size[ra] > cls_size[rb]) { tmp = ra; ra = rb; rb = tmp; }
            uf[ra] = rb;
            cls_size[rb] += cls_size[ra];
            cls_next[cls_tail[rb]] = cls_head[ra];
            cls_tail[rb] = cls_tail[ra];
        }

        MYBZERO(visited, nseq);

        // per-component: MST + path enumeration + sequential 2opts
        for (i = 0; i < nseq; i++) {
            ra = _uf_find_hic(uf, i);
            if (ra != i) continue;          // only root of each component
            if (cls_size[i] < 2) continue;  // singleton

            // collect edges belonging to this component
            comp_edges.n = 0;
            for (e = 0; e < (int) all_edges.n; e++) {
                a = all_edges.a[e].a;
                if (_uf_find_hic(uf, a) == i)
                    kv_push(hedge_t, comp_edges, all_edges.a[e]);
            }

            // reset adjacency list entries for nodes in this component
            s = cls_head[i];
            while (s >= 0) {
                kv_a(int, adj_head, s) = -1;
                s = cls_next[s];
            }
            adj_next.n = adj_to.n = adj_w.n = 0;

            // build MST via Kruskal on component edges
            MYMALLOC(uf2, nseq);
            s = cls_head[i];
            while (s >= 0) { uf2[s] = s; s = cls_next[s]; }
            for (e = 0; e < (int) comp_edges.n; e++) {
                a  = comp_edges.a[e].a;
                b  = comp_edges.a[e].b;
                ra = _uf_find_hic(uf2, a);
                rb = _uf_find_hic(uf2, b);
                if (ra == rb) continue;
                if (ra > rb) { tmp = ra; ra = rb; rb = tmp; }
                uf2[ra] = rb;
                // add edge a→b
                ei = adj_to.n;
                kv_push(int,   adj_to,   b);
                kv_push(int64, adj_w,    comp_edges.a[e].w);
                kv_push(int,   adj_next, kv_A(adj_head, a));
                kv_A(adj_head, a) = ei;
                // add edge b→a
                ei = adj_to.n;
                kv_push(int,   adj_to,   a);
                kv_push(int64, adj_w,    comp_edges.a[e].w);
                kv_push(int,   adj_next, kv_A(adj_head, b));
                kv_A(adj_head, b) = ei;
            }
            free(uf2);

            // DFS from each node to enumerate simple paths of exactly path_len nodes
            s = cls_head[i];
            while (s >= 0) {
                MYBZERO(in_path, nseq);
                path.n = dfs_stack.n = dfs_depth.n = 0;

                kv_push(int,  dfs_stack, s);
                kv_push(int8, dfs_depth, 0);

                while (dfs_stack.n > 0) {
                    cur   = kv_A(dfs_stack, dfs_stack.n - 1);
                    depth = kv_A(dfs_depth, dfs_depth.n - 1);

                    if (depth == (int) path.n) {
                        // entering cur: add to path
                        kv_push(int, path, cur);
                        in_path[cur] = 1;
                    } else {
                        // returning: pop path tail
                        in_path[path.a[path.n - 1]] = 0;
                        path.n--;
                        dfs_stack.n--;
                        dfs_depth.n--;
                        continue;
                    }

                    if ((int) path.n == path_len) {
                        // store canonical path (first node < last node avoids duplicate reverses)
                        if (path.a[0] < path.a[path_len - 1]) {
                            pw = 0.0;
                            for (k = 0; k + 1 < path_len; k++)
                                for (e = kv_A(adj_head, path.a[k]); e >= 0; e = kv_A(adj_next, e))
                                    if (kv_A(adj_to, e) == path.a[k + 1]) { pw += (double) kv_A(adj_w, e); break; }
                            for (k = 0; k < path_len; k++) hp.s[k] = path.a[k];
                            hp.v = pw;
                            kv_push(hs_path_t, paths, hp);
                        }
                        // backtrack
                        in_path[path.a[path.n - 1]] = 0;
                        path.n--;
                        dfs_stack.n--;
                        dfs_depth.n--;
                        continue;
                    }

                    // expand neighbours not already on path
                    pushed = 0;
                    for (e = kv_A(adj_head, cur); e >= 0; e = kv_A(adj_next, e)) {
                        nxt = kv_A(adj_to, e);
                        if (in_path[nxt]) continue;
                        kv_push(int,  dfs_stack, nxt);
                        kv_push(int8, dfs_depth, (int8) path.n);
                        pushed = 1;
                    }
                    if (!pushed) {
                        // leaf before target depth: backtrack
                        in_path[path.a[path.n - 1]] = 0;
                        path.n--;
                        dfs_stack.n--;
                        dfs_depth.n--;
                    }
                }

                s = cls_next[s];
            }
        }

        // sort collected paths by total HIC weight descending, then process
        qsort(paths.a, paths.n, sizeof(hs_path_t), hs_path_v_dcmpfunc);

        for (pi = 0; pi < (int) paths.n; pi++) {
            MYBZERO(hap_cnt, ploidy + 1);
            for (k = 0; k < path_len; k++)
                hap_cnt[(int) bhap[paths.a[pi].s[k]]]++;
            h_maj = 1;
            for (k = 2; k <= ploidy; k++)
                if (hap_cnt[k] > hap_cnt[h_maj]) h_maj = k;
            n_minority = path_len - hap_cnt[h_maj];
            if (n_minority == 0) continue;

            for (k = 0; k < path_len; k++) {
                ps  = paths.a[pi].s[k];
                h_a = (int) bhap[ps];
                if (h_a == h_maj) continue;

                memcpy(haps, bhap, nseq);
                MYBZERO(visited, nseq);
                hs_opts->n = 0;
                tcfts = best_cfts;
                bcfts = best_cfts;
                bstep = 0;
                queue.n = 0; qh = 0;

                visited[ps] = 1;
                kv_push(int, queue, ps);

                while (qh < (int) queue.n) {
                    qs     = queue.a[qh++];
                    h_from = haps[qs];
                    h_to   = (h_from == h_a) ? h_maj : h_a;

                    cft = cfts + (index[qs] >> 32);
                    n   = (uint32) index[qs];
                    cft_v_move(cft, (int) n, (uint8) h_from, (uint8) h_to, haps, &tcfts);
                    haps[qs] = h_to;
                    kv_push(hs_opt_t, *hs_opts, ((hs_opt_t){qs, (uint8) h_to}));

                    if (tcfts.vals[CFT_OVL] <= best_cfts.vals[CFT_OVL] &&
                        cft_v_cmpfunc(&tcfts, &bcfts) < 0) {
                        bcfts = tcfts;
                        bstep = hs_opts->n;
                    }

                    for (e = 0; e < (int) n; e++) {
                        if (cft[e].t != CFT_OVL) continue;
                        nb = cft[e].b;
                        if (visited[nb])      continue;
                        if (haps[nb] != h_to) continue;
                        visited[nb] = 1;
                        kv_push(int, queue, nb);
                    }
                }

                if (bstep > 0 && cft_v_cmpfunc(&bcfts, &best_cfts) < 0) {
                    for (m = 0; m < bstep; m++)
                        bhap[hs_opts->a[m].seq] = hs_opts->a[m].hap;
                    best_cfts = bcfts;
                    improved  = 1;
#ifdef DEBUG_REFINE_HAP_PARTITION
                    ncascade++;
#endif
                }
            }
        }

#ifdef DEBUG_REFINE_HAP_PARTITION
        npass++;
        fprintf(stderr, "[M::%s] path_len=%d pass %d: %d commits, conflicts [%lld %lld]\n",
            __func__, path_len, npass, ncascade, best_cfts.vals[CFT_OVL], best_cfts.vals[CFT_HIC]);
#endif

        if (!improved)
            path_len--;
    }

    for (i = 0; i < nseq; i++)
        seqs[i].hap = bhap[i];

#ifdef DEBUG_REFINE_HAP_PARTITION
    fprintf(stderr, "[M::%s] done, final conflicts [%lld %lld]\n",
        __func__, best_cfts.vals[CFT_OVL], best_cfts.vals[CFT_HIC]);
#endif

#ifdef DEBUG_HAPLOTYPE_PARTITION
    fprintf(stderr, "[M::%s] haplotype partition after HIC k-opt refinement\n", __func__);
    print_partition_stats(seqs, nseq, cfts, index, ploidy);
#endif

    free(uf);
    free(cls_head); 
    free(cls_tail); 
    free(cls_next); 
    free(cls_size);
    free(visited); 
    free(in_path); 
    free(hap_cnt);
    kv_destroy(paths);
    kv_destroy(all_edges); 
    kv_destroy(comp_edges);
    kv_destroy(adj_head); 
    kv_destroy(adj_next); 
    kv_destroy(adj_to); 
    kv_destroy(adj_w);
    kv_destroy(path); 
    kv_destroy(dfs_stack); 
    kv_destroy(dfs_depth);
    kv_destroy(queue);
}

static int *mark_repeat_overlaps(ovl_t *ovls, int64 novl, sdict_t *dicts, int ploidy, int *_reps)
{
    if (ovls == NULL || !novl)
        return NULL;

    cov_point_t **covs, *cov;
    ovl_t *ovl;
    uint64 *index;
    uint32 a, b;
    int *reps, *dels;
    int64 i, j, n, nseq, movl, lowa, lowb, sdel, ndel;
    int *ncov, *cpts;
    double min_qual = OVL_MIN_QUAL;
    
    // sequence number
    nseq = dicts->n;
    
    // allocate memory
    MYCALLOC(index, nseq);
    MYMALLOC(covs, nseq);
    MYMALLOC(ncov, nseq);
    MYMALLOC(covs[0], (novl+nseq)*2);
    MYMALLOC(dels, novl);
    if (_reps) reps = _reps;
    else MYMALLOC(reps, nseq);

    if (covs == NULL || covs[0] == NULL || ncov == NULL || 
        index == NULL || reps == NULL || dels == NULL)
        mem_alloc_error("scaffold partition arrays");
    
    // sequence number
    nseq = dicts->n;

    // max copy number for low-copy region
    ploidy += 1;

    // build overlap index for quick search
    // overlaps are sorted by abseqs
    a = ovls->aread;
    for (i = 1, j = 0; i < novl; i++) {
        if (ovls[i].aread != a) {
            index[a] = (uint64) j << 32 | (i - j);
            a = ovls[i].aread;
            j = i;
        }
    }
    index[a] = (uint64) j << 32 | (i - j);
    // max overlaps
    movl = 0;
    for (i = 0; i < nseq; i++) {
        n = index[i] >> 32;
        if (movl < n)
            movl = n;
    }
    MYMALLOC(cpts, movl*2);
    if (cpts == NULL)
        mem_alloc_error("coverage points");
    // calculate sequence coverage
    for (i = 0; i < nseq; i++) {
        if (i)
            covs[i] = covs[i-1] + ncov[i-1];
        ovl = ovls + (index[i] >> 32);
        n = (uint32) index[i];
        ncov[i] = calc_coverage_from_intervals(covs[i], cpts, dicts->s[i].len, ovl, n, &min_qual, pts_from_overlaps);
    }

    // find repeat sequences
    MYBZERO(reps, nseq);
    sdel = ndel = 0;
    for (i = 0; i < nseq; i++) {
        cov = covs[i];
        lowa = 0;
        for (j = 1; j < ncov[i]; j++)
            if (cov[j-1].cov < ploidy)
                lowa += cov[j].pos - cov[j-1].pos;
        // mark repeat
        if (lowa < dicts->s[i].len * overlap_min_lowcopy_ratio(dicts->s[i].len)) {
            reps[i] = 1; // mark repeat sequence
            ndel += 1;
            sdel += dicts->s[i].len;
        }
    }
    if (VERBOSE > 0)
        fprintf(stderr, "[M::%s] marked %lld sequences of %lld bp as repeats\n", __func__, ndel, sdel);
    // mark all overlaps on them as deleted
    MYBZERO(dels, novl);
    ndel = 0;
    for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
        if (reps[ovl->aread] || reps[ovl->bread]) {
            if (!ovl->del) {
                ovl->del = 1; // mark overlap as deleted
                dels[i] = 1; // mark deleted overlaps
                ndel += 1;
            }
        }
    }
    if (VERBOSE > 0)
        fprintf(stderr, "[M::%s] marked %lld associated overlaps\n", __func__, ndel);

    // recalculate coverage after marking repeat sequences
    for (i = 0; i < nseq; i++) {
        if (i)
            covs[i] = covs[i-1] + ncov[i-1];
        ovl = ovls + (index[i] >> 32);
        n = (uint32) index[i];
        ncov[i] = calc_coverage_from_intervals(covs[i], cpts, dicts->s[i].len, ovl, n, &min_qual, pts_from_overlaps);
    }
    free(cpts);

    // find overlaps in repeat regions and mark them as deleted
    ndel = 0;
    for (i = 0; i < novl; i++) {
        ovl = ovls + i;
        a = ovl->aread;
        b = ovl->bread;
        cov = covs[a];
        lowa = 0;
        // find low-copy regions on the overlap region
        for (j = 1; j < ncov[a]; j++) {
            if (cov[j].pos <= ovl->abpos) continue;
            if (cov[j-1].pos >= ovl->aepos) break;
            if (cov[j-1].cov < ploidy)
                lowa += MIN(cov[j].pos, ovl->aepos) - MAX(cov[j-1].pos, ovl->abpos);
        }
        lowb = 0;
        cov = covs[b];
        for (j = 1; j < ncov[b]; j++) {
            if (cov[j].pos <= ovl->bbpos) continue;
            if (cov[j-1].pos >= ovl->bepos) break;
            if (cov[j-1].cov < ploidy)
                lowb += MIN(cov[j].pos, ovl->bepos) - MAX(cov[j-1].pos, ovl->bbpos);
        }

        if (lowa < ovl->alen * overlap_min_lowcopy_ratio(ovl->alen) && 
            lowb < ovl->blen * overlap_min_lowcopy_ratio(ovl->blen))
            if (!ovl->del) {
                ovl->del = 1; // mark overlap as deleted
                dels[i] = 1; // mark deleted overlaps
                ndel += 1;
            }
    }
    if (VERBOSE > 0)
        fprintf(stderr, "[M::%s] marked %lld more overlaps in repeat regions\n", __func__, ndel);
    
    free(index);
    free(covs[0]);
    free(covs);
    free(ncov);
    if (!_reps) free(reps);

    return dels;
}

static inline int parse_line(char *s, char **f, int m)
{
    int n = 0;
    while (*s) {
        while (isspace(*s)) s++;
        if (!*s) break;
        f[n++] = s;
        while (*s && !isspace(*s)) s++;
        if (*s) *s++ = '\0';
        if (n >= m) break;
    }
    return n;
}

typedef struct {
    int a, b;
    int64 v;
} hcft_t;

static int hcft_ab_cmpfunc(const void *x, const void *y)
{   
    const hcft_t *a = (const hcft_t *) x;
    const hcft_t *b = (const hcft_t *) y;
    if (a->a != b->a)
        return (a->a > b->a) - (a->a < b->a);
    return (a->b > b->b) - (a->b < b->b);
}

typedef struct {
    int g, n;
    int *s;
    int64 t;
} grp_info_t;

static int grp_info_cmpfunc(const void *a, const void *b)
{
    const grp_info_t *ga = (const grp_info_t *) a;
    const grp_info_t *gb = (const grp_info_t *) b;
    if (ga->g == 0) return  1; // group 0 is always last
    if (gb->g == 0) return -1; // group 0 is always last
    if (ga->t != gb->t) return (ga->t < gb->t) - (ga->t > gb->t); // sort by total length descending
    if (ga->n != gb->n) return (ga->n < gb->n) - (ga->n > gb->n); // sort by number of sequences descending
    return (ga->g > gb->g) - (ga->g < gb->g); // sort by group id ascending
}

static inline int ccompare(int k, int p, int64 *hcft, int64 *hovl, int64 *hhlk, int64 *hseq) 
{
    int d;
    if ((d = (hcft[k] > hcft[p]) - (hcft[k] < hcft[p]))) return d;
    if ((d = (hovl[k] > hovl[p]) - (hovl[k] < hovl[p]))) return d;
    if ((d = (hhlk[p] > hhlk[k]) - (hhlk[p] < hhlk[k]))) return d;
    return (hseq[k] > hseq[p]) - (hseq[k] < hseq[p]);
}

typedef struct { int u, v; double w; } seed_edge_t;

static int cmp_seed_edge_desc(const void *pa, const void *pb) 
{
    double wa = ((const seed_edge_t *)pa)->w;
    double wb = ((const seed_edge_t *)pb)->w;
    return (wa < wb) - (wa > wb); // descending
}

#define SBIT_WORDS(n)       (((size_t)(n) + 63) / 64)
#define SBIT_ROW(adj,nw,v)  ((adj) + (size_t)(v) * (nw))
#define SBIT_SET(row,i)     ((row)[(i) >> 6] |= (1ULL << ((i) & 63)))

// Bounded backtracking: is there a `need`-clique fully inside bitset `cand`,
// using adjacency `adj` (ns rows of nw words each)? Recursion depth is
// capped at `need` (<= K-2), and each step only walks `cand`'s set bits,
// so cost tracks graph density near the seed edge, not ns itself.
static int seed_clique_in_set(uint64 *adj, size_t nw, uint64 *cand, int need, int *out)
{
    if (need == 0) return 1;

    size_t j, w, xw;
    uint64 bits, *adjx, next[nw];
    int x, xb, bit, popcount;
    
    popcount = 0;
    for (w = 0; w < nw; w++)
        popcount += __builtin_popcountll(cand[w]);
    if (popcount < need) return 0; // not enough candidates left - prune

    for (w = 0; w < nw; w++) {
        bits = cand[w];
        while (bits) {
            bit = __builtin_ctzll(bits);
            x = (int) (w * 64 + (size_t)bit);
            bits &= bits - 1;

            adjx = SBIT_ROW(adj, nw, x);
            for (j = 0; j < nw; j++) 
                next[j] = cand[j] & adjx[j];

            // restrict to indices > x, so each combination is only tried once
            xw = (size_t) x >> 6;
            xb = x & 63;
            if (xb == 63) 
                next[xw] = 0;
            else
                next[xw] &= ~((1ULL << (xb + 1)) - 1);
            for (j = 0; j < xw; j++)
                next[j] = 0;

            out[need - 1] = x;
            if (seed_clique_in_set(adj, nw, next, need - 1, out))
                return 1;
        }
    }
    return 0;
}

// Best bottleneck `target`-clique among local vertices [0, ns).
// `edges` must be pre-sorted descending by weight. On success, fills
// seed_out[0..target-1] with local vertex ids and returns 1.
static int find_seed_clique(int ns, seed_edge_t *edges, int nedges, int target, uint64 *seed_out)
{
    if (target < 2) return 0;

    size_t j, nw;
    uint64 *adj;
    int e, u, v, t, found;

    nw = SBIT_WORDS(ns);
    MYCALLOC(adj, (size_t) ns * nw);
    if (adj == NULL)
        mem_alloc_error("seed clique adjacency bitset");

    found = 0;
    for (e = 0; e < nedges && !found; e++) {
        u = edges[e].u;
        v = edges[e].v;
        SBIT_SET(SBIT_ROW(adj, nw, u), v);
        SBIT_SET(SBIT_ROW(adj, nw, v), u);

        if (target == 2) {
            seed_out[0] = u;
            seed_out[1] = v;
            found = 1;
            break;
        }

        uint64 common[nw];
        uint64 *ru = SBIT_ROW(adj, nw, u), *rv = SBIT_ROW(adj, nw, v);
        for (j = 0; j < nw; j++)
            common[j] = ru[j] & rv[j];

        int rest[target];
        if (seed_clique_in_set(adj, nw, common, target - 2, rest)) {
            seed_out[0] = u; 
            seed_out[1] = v;
            for (t = 0; t < target - 2; t++) 
                seed_out[2 + t] = rest[t];
            found = 1;
        }
    }

    free(adj);
    return found;
}

long   HIGHS_TIME_LIMIT  = 300;   // seconds
double HIGHS_MIP_REL_GAP = 0.001; // 1% relative gap

#define PHASE_UNRESOLVED   0
#define PHASE_ILP_FEASIBLE 1
#define PHASE_ILP_OPTIMUM  2
#define PHASE_HEU_APPROX   4

typedef struct {
    grp_info_t *grps;
    hap_info_t *haps;
    ovl_t *ovls;
    hlk_t *hlks;
    ord_dbl_t *hord;
    int *sidx;
    int K;
    int *ns;
    int *no;
    int *nh;
    int **ols;
    int **hls;
    int *status;
} phase_data_t;

static void phase_data_init(phase_data_t *data, int n)
{
    MYBZERO(data, 1);
    MYCALLOC(data->ns, n);
    MYCALLOC(data->no, n);
    MYCALLOC(data->nh, n);
    MYCALLOC(data->ols, n);
    MYCALLOC(data->hls, n);
    MYCALLOC(data->status, n);

    if (data->ns == NULL || data->no == NULL || data->nh == NULL || 
        data->ols == NULL || data->hls == NULL || data->status == NULL)
        mem_alloc_error("phase_data_t arrays");
}

static void phase_data_destroy(phase_data_t *data)
{
    free(data->ns);
    free(data->no);
    free(data->nh);
    free(data->ols);
    free(data->hls);
    free(data->status);
}

static pthread_mutex_t plock = PTHREAD_MUTEX_INITIALIZER;

static double PHASE_ALPHA = 1.0; // weight for overlap conflict in objective function
static double PHASE_BETA  = 1.0; // weight for hic link support in objective function

static int event_ord_acmpfunc(const void *a, const void *b)
{
    ord_i64_t *x = (ord_i64_t *) a;
    ord_i64_t *y = (ord_i64_t *) b;
    int xm, ym;

    xm = llabs(x->event);
    ym = llabs(y->event);
    
    if (xm == ym) {
        xm = x->event;
        ym = y->event;
    }
    
    return (xm > ym) - (xm < ym);
}

static int u64_dcmp_func(const void *a, const void *b)
{
    uint64 x = *(const uint64 *) a;
    uint64 y = *(const uint64 *) b;
    return (x < y) - (x > y);
}

// record of MIP progress for fallback
// only the last 32 records are kept
typedef struct {
    double running_time;
    double mip_primal_bound;
    double mip_dual_bound;
} mip_track_t;

typedef struct {
	int front, count;
	mip_track_t a[32];
} mip_callback_t;

static inline mip_track_t *mf_push(mip_callback_t *q)
{
	return &q->a[((q->count++) + q->front) & 0x1f];
}

static inline mip_track_t *mf_at(mip_callback_t *q, int i)
{
    if (q->count <= i) return NULL;
    return &q->a[(q->front + i) & 0x1f];
}

static void mip_callback(
    const int callback_type,
    const char *message,
    const HighsCallbackDataOut *out,
    HighsCallbackDataIn *in,
    void *user_data)
{
    mip_callback_t *fallback = (mip_callback_t *)user_data;
    mip_track_t *track = mf_push(fallback);

    if (callback_type == kHighsCallbackMipImprovingSolution) {
        track->mip_primal_bound = out->mip_primal_bound;
        track->running_time = out->running_time;
        track->mip_dual_bound = out->mip_dual_bound;
    }
}

static int MAX_ILP_SEQS = 300;

/*
 * Rolling-horizon MILP phasing.
 *
 * Main ideas:
 *   1. Solve ~target_seq new sequences at a time.
 *   2. Choose the block boundary near target_seq where the number of
 *      pseudo-backbone intervals crossing the boundary is small, while
 *      preferring >= K crossing sequences.
 *   3. After solving a block, choose one already-solved anchor for each
 *      haplotype.  A six-clique is NOT required.
 *   4. Carry these anchors into the next MILP with their x variables fixed.
 *   5. Edges from older fixed sequences that are not carried explicitly
 *      are retained exactly, conditional on their fixed assignment, as
 *      unary costs on the current x variables.
 *
 * The decomposition is heuristic globally: once a sequence leaves the free
 * window its haplotype assignment is not reconsidered.
 */

#define PHASE_ILP_IS_SOLVED(s) \
    ((s) == PHASE_ILP_OPTIMUM || (s) == PHASE_ILP_FEASIBLE)

typedef enum {
    PHASE_MIP_FAILED   = 0,
    PHASE_MIP_FEASIBLE = 1,
    PHASE_MIP_OPTIMAL  = 2
} phase_mip_result_t;

typedef enum {
    PHASE_ROLL_FAILED = 0,
    PHASE_ROLL_APPROX = 1,
    PHASE_ROLL_EXACT  = 2
} phase_roll_result_t;

#define PHASE_ROLL_MIN_FREE 10
#define PHASE_ROLL_GROW_STEP 10
#define PHASE_ROLL_SHRINK_NUM 7
#define PHASE_ROLL_SHRINK_DEN 10

typedef struct {
    int gid;                    /* index in grps[jid].s[] / group-local index */
    long long bpos;
    long long epos;
} phase_ord_t;

typedef struct {
    phase_data_t *data;
    long jid;

    grp_info_t *grps;
    hap_info_t *haps;
    ovl_t *ovls;
    hlk_t *hlks;
    ord_dbl_t *hord;

    int K;
    int ns;
    int nh;
    int no;

    int *sidx;
    int *ols;
    int *hls;
    int *gseq;
} phase_blk_ctx_t;


static int cmp_phase_ord_bpos(const void *a, const void *b)
{
    const phase_ord_t *x = (const phase_ord_t *)a;
    const phase_ord_t *y = (const phase_ord_t *)b;

    if (x->bpos < y->bpos) return -1;
    if (x->bpos > y->bpos) return  1;
    if (x->epos > y->epos) return -1;  /* longer first for identical bpos */
    if (x->epos < y->epos) return  1;
    return x->gid - y->gid;
}


/*
 * Select the exclusive right end of the next free block.
 *
 * We look around start + target_seq and prefer a cut with:
 *   - at least K active pseudo-backbone intervals;
 *   - a small active frontier;
 *   - a block size close to target_seq.
 *
 * cut_pos is the bpos of the first sequence not included in this block.
 */
static int choose_phase_block_end(const phase_blk_ctx_t *ctx,
                                  const phase_ord_t *ord,
                                  int start,
                                  int target_seq,
                                  long long *cut_pos)
{
    int rem = ctx->ns - start;
    int d, lo, hi, r, j;
    int best = -1;
    double best_score = 1e300;

    if (rem <= target_seq) {
        *cut_pos = 0x7fffffffffffffffLL;
        return ctx->ns;
    }

    d = target_seq / 5;
    if (d < 5) d = 5;

    lo = start + target_seq - d;
    hi = start + target_seq + d;

    if (lo < start + ctx->K) lo = start + ctx->K;
    if (hi >= ctx->ns) hi = ctx->ns - 1;

    for (r = lo; r <= hi; r++) {
        long long p;
        int active = 0;
        int nfree = r - start;
        double score;

        /* avoid splitting a set of identical starting coordinates */
        if (r > 0 && ord[r - 1].bpos == ord[r].bpos)
            continue;

        p = ord[r].bpos;

        /*
         * Count all already-started intervals crossing p.  This includes
         * previously solved long sequences as well as sequences in the
         * current block.
         */
        for (j = 0; j < r; j++)
            if (ord[j].epos > p)
                active++;

        if (active < ctx->K)
            continue;

        /*
         * Frontier size is the primary criterion.  The second term keeps
         * the block close to target_seq.
         */
        score = (double)active + 0.20 * fabs((double)nfree - target_seq);

        if (score < best_score) {
            best_score = score;
            best = r;
        }
    }

    /* No geometrically attractive K-covered cut: fall back to target size. */
    if (best < 0) {
        best = start + target_seq;
        if (best >= ctx->ns)
            best = ctx->ns;

        /* Do not split identical bpos values if we can avoid it. */
        while (best < ctx->ns &&
               best > 0 &&
               ord[best - 1].bpos == ord[best].bpos)
            best++;
    }

    if (best >= ctx->ns)
        *cut_pos = 0x7fffffffffffffffLL;
    else
        *cut_pos = ord[best].bpos;

    return best;
}


/*
 * Select one carry-over anchor for each haplotype.
 *
 * Preference is given to solved sequences having the greatest total edge
 * weight to still-unsolved sequences.  If no such edge exists for a
 * haplotype, the sequence extending furthest to the right is used.
 *
 * rank[gid] gives the pseudo-backbone order of group-local sequence gid.
 *
 * Returns the number of haplotypes for which an anchor was found.
 */
static int select_phase_anchors(const phase_blk_ctx_t *ctx,
                                const phase_ord_t *ord,
                                const int *rank,
                                int right,
                                long long cut_pos,
                                int *anchors,
                                double *anchor_score)
{
    double *future_score;
    double best_score[ctx->K];
    long long best_epos[ctx->K];
    int i, e, k, gu, gv, ru, rv;
    int found = 0;

    MYCALLOC(future_score, ctx->ns);
    if (future_score == NULL)
        mem_alloc_error("phase anchor future scores");

    for (k = 0; k < ctx->K; k++) {
        anchors[k] = -1;
        anchor_score[k] = 0.0;
        best_score[k] = -1.0;
        best_epos[k] = -1;
    }

    /* Hi-C edges crossing solved -> unsolved. */
    for (e = 0; e < ctx->nh; e++) {
        hlk_t *hlk = ctx->hlks + ctx->hls[e];
        double w;

        gu = ctx->sidx[hlk->a];
        gv = ctx->sidx[hlk->b];

        if (gu < 0 || gu >= ctx->ns || gv < 0 || gv >= ctx->ns)
            continue;

        ru = rank[gu];
        rv = rank[gv];

        w = hlk->l * HIC_NORM_WINDOW *
            ctx->hord[ctx->hls[e]].event / 1000000.0 * PHASE_BETA;

        if (ru < right && rv >= right)
            future_score[gu] += w;
        else if (rv < right && ru >= right)
            future_score[gv] += w;
    }

    /* OVL edges crossing solved -> unsolved. */
    for (e = 0; e < ctx->no; e++) {
        ovl_t *ovl = ctx->ovls + ctx->ols[e];
        double w;

        if (ovl->del) continue;

        gu = ctx->sidx[ovl->aread];
        gv = ctx->sidx[ovl->bread];

        if (gu < 0 || gu >= ctx->ns || gv < 0 || gv >= ctx->ns)
            continue;

        ru = rank[gu];
        rv = rank[gv];

        w = (double)ovl->alen + ovl->blen;
        w = w * log10(w) / 1000000.0 * PHASE_ALPHA;

        if (ru < right && rv >= right)
            future_score[gu] += w;
        else if (rv < right && ru >= right)
            future_score[gv] += w;
    }

    /*
     * Pick the best solved representative of each haplotype.
     * A positive future_score is preferred; epos breaks ties and supplies
     * a fallback when there is no direct edge into the future.
     */
    for (i = 0; i < right; i++) {
        int gid = ord[i].gid;
        int hap = ctx->haps[ctx->gseq[gid]].hap;

        if (hap < 1 || hap > ctx->K)
            continue;

        k = hap - 1;

        if (anchors[k] < 0 ||
            future_score[gid] > best_score[k] ||
            (future_score[gid] == best_score[k] &&
             ord[i].epos > best_epos[k])) {

            anchors[k] = gid;
            best_score[k] = future_score[gid];
            best_epos[k] = ord[i].epos;
            anchor_score[k] = future_score[gid];
        }
    }

    for (k = 0; k < ctx->K; k++)
        if (anchors[k] >= 0)
            found++;

    if (VERBOSE > 2) {
        pthread_mutex_lock(&plock);
        for (k = 0; k < ctx->K; k++) {
            if (anchors[k] >= 0) {
                hap_info_t *h = ctx->haps + ctx->gseq[anchors[k]];
                fprintf(stderr,
                        "[D::select_phase_anchors] hap=%d gid=%d "
                        "bpos=%lld epos=%lld future_score=%.6g cross=%d\n",
                        k + 1, anchors[k],
                        (long long)h->bpos, (long long)h->epos,
                        anchor_score[k],
                        cut_pos != 0x7fffffffffffffffLL &&
                        (long long)h->epos > cut_pos);
            } else {
                fprintf(stderr,
                        "[D::select_phase_anchors] hap=%d NO ANCHOR\n",
                        k + 1);
            }
        }
        pthread_mutex_unlock(&plock);
    }

    free(future_score);
    return found;
}


/*
 * Solve one rolling-horizon subproblem.
 *
 * Free sequences are ord[start:right].
 * anchors[] contains previously solved group-local sequence indices whose
 * hap assignments are fixed in this MILP.
 *
 * The function writes hap assignments only for the free sequences.
 *
 * Returns 1 if a feasible integer assignment was recovered, otherwise 0.
 */

/*
 * Reduce a failed block while respecting identical bpos values.
 * Returns an exclusive right index; returns start if no smaller block
 * of at least min_free sequences can be formed.
 */
static int shrink_phase_block_end(const phase_ord_t *ord,
                                  int start,
                                  int right,
                                  int min_free)
{
    int n = right - start;
    int new_n;
    int r;

    if (n <= min_free)
        return start;

    new_n = n * PHASE_ROLL_SHRINK_NUM / PHASE_ROLL_SHRINK_DEN;
    if (new_n >= n)
        new_n = n - 1;
    if (new_n < min_free)
        new_n = min_free;

    r = start + new_n;

    /*
     * Prefer moving the cut left across a group with identical starting
     * coordinates.  If that would violate min_free, move right instead.
     */
    while (r > start + min_free &&
           ord[r - 1].bpos == ord[r].bpos)
        r--;

    while (r < right &&
           ord[r - 1].bpos == ord[r].bpos)
        r++;

    if (r >= right)
        return start;

    return r;
}


static phase_mip_result_t solve_phase_subproblem(const phase_blk_ctx_t *ctx,
                                  const phase_ord_t *ord,
                                  int start,
                                  int right,
                                  const int *anchors,
                                  int nanchor,
                                  int first_block)
{
    int nfree = right - start;
    int nloc_max = nanchor + nfree;
    int nloc = 0;
    int *map = NULL;            /* group-local gid -> local MILP index */
    int *local_gid = NULL;      /* local MILP index -> group-local gid */
    int *free_li = NULL;        /* free-sequence index -> local MILP index */
    int *lhic = NULL;           /* indices into ctx->hls[] */
    int *lovl = NULL;           /* indices into ctx->ols[] */
    int nlh = 0, nlo = 0;
    int a, q, i, e, k, gu, gv, lu, lv;
    int wd;
    HighsInt *idxs = NULL;
    double *vals = NULL;
    double *col_costs = NULL;
    double *col_value = NULL, *col_dual = NULL;
    double *row_value = NULL, *row_dual = NULL;
    void *highs = NULL;
    double inf;
    HighsInt offset_x, offset_z, offset_o, tot_vars;
    HighsInt col_idx, num_cols, num_rows;
    HighsInt run_status, model_status;
    int feasible = 0;
    phase_mip_result_t result = PHASE_MIP_FAILED;

    if (nfree <= 0)
        return PHASE_MIP_OPTIMAL;

    MYMALLOC(map, ctx->ns);
    MYMALLOC(local_gid, nloc_max);
    MYMALLOC(free_li, nfree);

    if (ctx->nh > 0) MYMALLOC(lhic, ctx->nh);
    if (ctx->no > 0) MYMALLOC(lovl, ctx->no);

    if (map == NULL || local_gid == NULL || free_li == NULL ||
        (ctx->nh > 0 && lhic == NULL) ||
        (ctx->no > 0 && lovl == NULL))
        mem_alloc_error("rolling phase subproblem arrays");

    for (i = 0; i < ctx->ns; i++)
        map[i] = -1;

    /* Add fixed carry-over anchors first. */
    for (a = 0; a < nanchor; a++) {
        int gid = anchors[a];

        if (gid < 0 || gid >= ctx->ns)
            continue;
        if (map[gid] >= 0)
            continue;
        if (ctx->haps[ctx->gseq[gid]].hap < 1 ||
            ctx->haps[ctx->gseq[gid]].hap > ctx->K)
            continue;

        map[gid] = nloc;
        local_gid[nloc++] = gid;
    }

    /* Add the new/free sequences. */
    for (q = 0; q < nfree; q++) {
        int gid = ord[start + q].gid;

        if (map[gid] >= 0) {
            free_li[q] = map[gid];
            continue;
        }

        map[gid] = nloc;
        free_li[q] = nloc;
        local_gid[nloc++] = gid;
    }

    /* Internal Hi-C edges. */
    for (e = 0; e < ctx->nh; e++) {
        hlk_t *hlk = ctx->hlks + ctx->hls[e];

        gu = ctx->sidx[hlk->a];
        gv = ctx->sidx[hlk->b];
        if (gu < 0 || gu >= ctx->ns || gv < 0 || gv >= ctx->ns)
            continue;

        lu = map[gu];
        lv = map[gv];

        if (lu >= 0 && lv >= 0)
            lhic[nlh++] = e;
    }

    /* Internal OVL edges. */
    for (e = 0; e < ctx->no; e++) {
        ovl_t *ovl = ctx->ovls + ctx->ols[e];

        gu = ctx->sidx[ovl->aread];
        gv = ctx->sidx[ovl->bread];
        if (gu < 0 || gu >= ctx->ns || gv < 0 || gv >= ctx->ns)
            continue;

        lu = map[gu];
        lv = map[gv];

        if (lu >= 0 && lv >= 0)
            lovl[nlo++] = e;
    }

    highs = Highs_create();
    inf = Highs_getInfinity(highs);

    Highs_setBoolOptionValue(highs, "output_flag", 0);
    Highs_setBoolOptionValue(highs, "presolve", 1);

    wd = ctx->K > 3 ? ctx->K : 3;
    MYMALLOC(idxs, wd);
    MYMALLOC(vals, wd);
    if (idxs == NULL || vals == NULL)
        mem_alloc_error("rolling phase row arrays");

    offset_x = 0;
    offset_z = offset_x + (HighsInt)nloc * ctx->K;
    offset_o = offset_z + (HighsInt)nlh * ctx->K;
    tot_vars = offset_o + nlo;

#define BIDX_X(ii, kk) (offset_x + (HighsInt)(ii) * ctx->K + (kk))
#define BIDX_Z(ee, kk) (offset_z + (HighsInt)(ee) * ctx->K + (kk))
#define BIDX_O(ee)     (offset_o + (ee))

    /* x_ik variables. */
    for (i = 0; i < nloc; i++) {
        int gid = local_gid[i];
        int fixed_hap = ctx->haps[ctx->gseq[gid]].hap;

        for (k = 0; k < ctx->K; k++) {
            double lb = 0.0, ub = 1.0;

            /*
             * A sequence preceding start is a carry-over anchor and is
             * already fixed.  New sequences ord[start:right] have hap=0
             * when entering this function.
             */
            if (fixed_hap >= 1 && fixed_hap <= ctx->K) {
                if (k == fixed_hap - 1)
                    lb = ub = 1.0;
                else
                    lb = ub = 0.0;
            }

            col_idx = BIDX_X(i, k);
            Highs_addCol(highs, 0.0, lb, ub, 0, NULL, NULL);
            Highs_changeColIntegrality(highs, col_idx, kHighsVarTypeInteger);
        }
    }

    /* Internal Hi-C z variables. */
    for (e = 0; e < nlh; e++) {
        for (k = 0; k < ctx->K; k++) {
            col_idx = BIDX_Z(e, k);
            Highs_addCol(highs, 0.0, 0.0, 1.0, 0, NULL, NULL);
            Highs_changeColIntegrality(highs, col_idx,
                                       kHighsVarTypeContinuous);
        }
    }

    /* Internal overlap o variables. */
    for (e = 0; e < nlo; e++) {
        col_idx = BIDX_O(e);
        Highs_addCol(highs, 0.0, 0.0, 1.0, 0, NULL, NULL);
        Highs_changeColIntegrality(highs, col_idx,
                                   kHighsVarTypeContinuous);
    }

    /*
     * First block only: use the largest available OVL clique as a symmetry
     * seed, but DO NOT require a K-clique.  If there is no useful clique,
     * fix one arbitrary sequence to haplotype 1, which is pure label
     * symmetry breaking.
     *
     * Subsequent blocks already have fixed carry-over anchors and therefore
     * do not need a new seed clique.
     */
    if (first_block && nanchor == 0) {
        seed_edge_t *sedges = NULL;
        uint64 clique[ctx->K];
        int nsedges = 0;
        int target, clique_size = 0;

        if (nlo > 0) {
            MYMALLOC(sedges, nlo);
            if (sedges == NULL)
                mem_alloc_error("rolling phase seed edges");

            for (e = 0; e < nlo; e++) {
                ovl_t *ovl = ctx->ovls + ctx->ols[lovl[e]];
                if (ovl->del) continue; /* seed graph follows the original code */
                gu = ctx->sidx[ovl->aread];
                gv = ctx->sidx[ovl->bread];

                sedges[nsedges].u = map[gu];
                sedges[nsedges].v = map[gv];
                sedges[nsedges].w =
                    (double)ovl->alen + (double)ovl->blen;
                nsedges++;
            }

            qsort(sedges, nsedges, sizeof(seed_edge_t),
                  cmp_seed_edge_desc);

            for (target = ctx->K;
                 target >= 2 && clique_size == 0;
                 target--) {
                if (find_seed_clique(nloc, sedges, nsedges,
                                     target, clique))
                    clique_size = target;
            }

            if (clique_size > 0) {
                /*
                 * Use the same deterministic label ordering as the full
                 * MILP: longer seed sequences receive lower haplotype
                 * indices.  This makes the single-block rolling model
                 * identical to the full model with respect to seed fixing.
                 */
                for (k = 0; k < clique_size; k++) {
                    int b = (int)clique[k];
                    int gid = local_gid[b];
                    clique[k] =
                        ((uint64)ctx->haps[ctx->gseq[gid]].len << 32) |
                        (uint32)b;
                }
                qsort(clique, clique_size, sizeof(uint64), u64_dcmp_func);

                for (k = 0; k < clique_size; k++) {
                    int b = (uint32)clique[k];

                    Highs_changeColBounds(highs, BIDX_X(b, k),
                                          1.0, 1.0);
                    for (i = 0; i < ctx->K; i++) {
                        if (i == k) continue;
                        Highs_changeColBounds(highs, BIDX_X(b, i),
                                              0.0, 0.0);
                    }
                }
            }

            free(sedges);
        }

        if (clique_size == 0 && nloc > 0) {
            Highs_changeColBounds(highs, BIDX_X(0, 0), 1.0, 1.0);
            for (k = 1; k < ctx->K; k++)
                Highs_changeColBounds(highs, BIDX_X(0, k), 0.0, 0.0);
        }
    }

    /* Exactly one haplotype per local sequence. */
    for (i = 0; i < nloc; i++) {
        for (k = 0; k < ctx->K; k++) {
            idxs[k] = BIDX_X(i, k);
            vals[k] = 1.0;
        }
        Highs_addRow(highs, 1.0, 1.0, ctx->K, idxs, vals);
    }

    /* Internal Hi-C linearisation. */
    for (e = 0; e < nlh; e++) {
        int ge = lhic[e];
        hlk_t *hlk = ctx->hlks + ctx->hls[ge];

        gu = ctx->sidx[hlk->a];
        gv = ctx->sidx[hlk->b];
        lu = map[gu];
        lv = map[gv];

        for (k = 0; k < ctx->K; k++) {
            idxs[0] = BIDX_Z(e, k); vals[0] =  1.0;

            idxs[1] = BIDX_X(lu, k); vals[1] = -1.0;
            Highs_addRow(highs, -inf, 0.0, 2, idxs, vals);

            idxs[1] = BIDX_X(lv, k); vals[1] = -1.0;
            Highs_addRow(highs, -inf, 0.0, 2, idxs, vals);
        }
    }

    /* Internal OVL linearisation. */
    for (e = 0; e < nlo; e++) {
        int ge = lovl[e];
        ovl_t *ovl = ctx->ovls + ctx->ols[ge];

        gu = ctx->sidx[ovl->aread];
        gv = ctx->sidx[ovl->bread];
        lu = map[gu];
        lv = map[gv];

        idxs[0] = BIDX_O(e); vals[0] = 1.0;

        for (k = 0; k < ctx->K; k++) {
            idxs[1] = BIDX_X(lu, k); vals[1] = -1.0;
            idxs[2] = BIDX_X(lv, k); vals[2] = -1.0;
            Highs_addRow(highs, -1.0, inf, 3, idxs, vals);
        }
    }

    free(idxs);
    free(vals);
    idxs = NULL;
    vals = NULL;

    /* Objective. */
    num_cols = Highs_getNumCol(highs);
    MYCALLOC(col_costs, num_cols);
    if (col_costs == NULL)
        mem_alloc_error("rolling phase objective");

    Highs_changeObjectiveSense(highs, -1);

    /* Internal Hi-C terms. */
    for (e = 0; e < nlh; e++) {
        int ge = lhic[e];
        double w = ctx->hlks[ctx->hls[ge]].l *
                   HIC_NORM_WINDOW *
                   ctx->hord[ctx->hls[ge]].event /
                   1000000.0 * PHASE_BETA;

        for (k = 0; k < ctx->K; k++)
            col_costs[BIDX_Z(e, k)] += w;
    }

    /* Internal OVL terms. */
    for (e = 0; e < nlo; e++) {
        int ge = lovl[e];
        ovl_t *ovl = ctx->ovls + ctx->ols[ge];
        double w = (double)ovl->alen + ovl->blen;

        w = w * log10(w) / 1000000.0 * PHASE_ALPHA;
        col_costs[BIDX_O(e)] -= w;
    }

    /*
     * Edges from a local sequence to an older, already-fixed external
     * sequence are exact unary terms conditional on that fixed assignment.
     *
     * Hi-C:
     *      +w [h_local == h_fixed]  -> +w x_local,h_fixed
     *
     * OVL:
     *      -w [h_local == h_fixed]  -> -w x_local,h_fixed
     */
    for (e = 0; e < ctx->nh; e++) {
        hlk_t *hlk = ctx->hlks + ctx->hls[e];
        int ext_gid, loc;
        int hap;
        double w;

        gu = ctx->sidx[hlk->a];
        gv = ctx->sidx[hlk->b];
        if (gu < 0 || gu >= ctx->ns || gv < 0 || gv >= ctx->ns)
            continue;

        lu = map[gu];
        lv = map[gv];

        if ((lu >= 0) == (lv >= 0))
            continue;  /* both local or both external */

        if (lu >= 0) {
            loc = lu;
            ext_gid = gv;
        } else {
            loc = lv;
            ext_gid = gu;
        }

        hap = ctx->haps[ctx->gseq[ext_gid]].hap;
        if (hap < 1 || hap > ctx->K)
            continue;  /* future/unresolved external endpoint */

        w = hlk->l * HIC_NORM_WINDOW *
            ctx->hord[ctx->hls[e]].event /
            1000000.0 * PHASE_BETA;

        col_costs[BIDX_X(loc, hap - 1)] += w;
    }

    for (e = 0; e < ctx->no; e++) {
        ovl_t *ovl = ctx->ovls + ctx->ols[e];
        int ext_gid, loc;
        int hap;
        double w;

        gu = ctx->sidx[ovl->aread];
        gv = ctx->sidx[ovl->bread];
        if (gu < 0 || gu >= ctx->ns || gv < 0 || gv >= ctx->ns)
            continue;

        lu = map[gu];
        lv = map[gv];

        if ((lu >= 0) == (lv >= 0))
            continue;

        if (lu >= 0) {
            loc = lu;
            ext_gid = gv;
        } else {
            loc = lv;
            ext_gid = gu;
        }

        hap = ctx->haps[ctx->gseq[ext_gid]].hap;
        if (hap < 1 || hap > ctx->K)
            continue;

        w = (double)ovl->alen + ovl->blen;
        w = w * log10(w) / 1000000.0 * PHASE_ALPHA;

        col_costs[BIDX_X(loc, hap - 1)] -= w;
    }

    Highs_changeColsCostByRange(highs, 0, num_cols - 1, col_costs);
    free(col_costs);
    col_costs = NULL;

    Highs_setDoubleOptionValue(highs, "mip_rel_gap",
                               HIGHS_MIP_REL_GAP);
    if (HIGHS_TIME_LIMIT > 0)
        Highs_setDoubleOptionValue(highs, "time_limit",
                                   HIGHS_TIME_LIMIT);

    if (VERBOSE > 1) {
        pthread_mutex_lock(&plock);
        fprintf(stderr,
                "[M::%s] rolling MILP: free=%d anchors=%d "
                "local=%d hic=%d ovl=%d vars=%lld rows=%lld\n",
                __func__, nfree, nloc - nfree, nloc, nlh, nlo,
                (long long)tot_vars,
                (long long)Highs_getNumRow(highs));
        pthread_mutex_unlock(&plock);
    }

    run_status = Highs_run(highs);
    model_status = Highs_getModelStatus(highs);

    num_cols = Highs_getNumCol(highs);
    num_rows = Highs_getNumRow(highs);

    MYMALLOC(col_value, num_cols);
    MYMALLOC(col_dual, num_cols);
    MYMALLOC(row_value, num_rows);
    MYMALLOC(row_dual, num_rows);

    if (col_value == NULL || col_dual == NULL ||
        row_value == NULL || row_dual == NULL)
        mem_alloc_error("rolling phase solution arrays");

    /*
     * A time-limit status can still contain a perfectly usable incumbent.
     * Rather than requiring model_status == Optimal, retrieve the final
     * solution and explicitly validate the free x assignments.
     */
    if (run_status >= 0) {
        Highs_getSolution(highs, col_value, col_dual,
                          row_value, row_dual);

        feasible = 1;

        for (q = 0; q < nfree; q++) {
            int li = free_li[q];
            int nk = 0;

            for (k = 0; k < ctx->K; k++)
                if (col_value[BIDX_X(li, k)] > 0.5)
                    nk++;

            if (nk != 1) {
                feasible = 0;
                break;
            }
        }

        if (feasible) {
            for (q = 0; q < nfree; q++) {
                int gid = ord[start + q].gid;
                int li = free_li[q];

                ctx->haps[ctx->gseq[gid]].hap = 0;

                for (k = 0; k < ctx->K; k++) {
                    if (col_value[BIDX_X(li, k)] > 0.5) {
                        ctx->haps[ctx->gseq[gid]].hap = k + 1;
                        break;
                    }
                }
            }

            result = (model_status == kHighsModelStatusOptimal) ?
                     PHASE_MIP_OPTIMAL : PHASE_MIP_FEASIBLE;
        }
    }

    if (VERBOSE > 1) {
        double time = Highs_getRunTime(highs);
        double primal = Highs_getObjectiveValue(highs);
        double dual = 0.0;

        Highs_getDoubleInfoValue(highs, "mip_dual_bound", &dual);

        pthread_mutex_lock(&plock);
        fprintf(stderr,
                "[M::%s] rolling MILP finished: "
                "run_status=%lld model_status=%lld feasible=%d "
                "time=%.2f primal=%.6f dual=%.6f\n",
                __func__,
                (long long)run_status, (long long)model_status,
                feasible, time, primal, dual);
        pthread_mutex_unlock(&plock);
    }

    free(col_value);
    free(col_dual);
    free(row_value);
    free(row_dual);

    Highs_destroy(highs);

    free(map);
    free(local_gid);
    free(free_li);
    free(lhic);
    free(lovl);

#undef BIDX_X
#undef BIDX_Z
#undef BIDX_O

    return result;
}


/*
 * Run rolling-horizon phasing.
 *
 * Returns:
 *   PHASE_ROLL_EXACT:
 *       one unsplit MILP containing every sequence was solved to proven
 *       optimality.  The result is globally optimal for this MILP and no
 *       second full-MILP pass is needed.
 *
 *   PHASE_ROLL_APPROX:
 *       a valid assignment was obtained, but either decomposition was used
 *       or at least one local MILP was not proved optimal.
 *
 *   PHASE_ROLL_FAILED:
 *       no complete rolling assignment could be obtained.
 *
 * This function deliberately does not change data->status[jid]; the outer
 * hybrid_phase_ilp_thread() decides whether the result is globally optimal
 * or only a feasible/suboptimal incumbent.
 */
static phase_roll_result_t hybrid_phase_rolling_solve(void *_data,
                                                       long jid,
                                                       int tid)
{
    phase_blk_ctx_t ctx;
    phase_ord_t *ord = NULL;
    int *rank = NULL;
    int *anchors = NULL;
    double *anchor_score = NULL;

    int target_seq;
    int start, right;
    int nanchor = 0;
    int block_id = 0;
    int first_block = 1;
    int decomposed = 0;
    int all_local_optimal = 1;
    int k;

    (void)tid;

    memset(&ctx, 0, sizeof(ctx));

    ctx.data = (phase_data_t *)_data;
    ctx.jid = jid;
    ctx.grps = ctx.data->grps;
    ctx.haps = ctx.data->haps;
    ctx.ovls = ctx.data->ovls;
    ctx.hlks = ctx.data->hlks;
    ctx.hord = ctx.data->hord;
    ctx.K = ctx.data->K;
    ctx.ns = ctx.data->ns[jid];
    ctx.nh = ctx.data->nh[jid];
    ctx.no = ctx.data->no[jid];
    ctx.sidx = ctx.data->sidx;
    ctx.ols = ctx.data->ols[jid];
    ctx.hls = ctx.data->hls[jid];
    ctx.gseq = ctx.grps[jid].s;

    target_seq = MAX_ILP_SEQS / ctx.K / 10 * 10;
    if (target_seq < PHASE_ROLL_MIN_FREE)
        target_seq = PHASE_ROLL_MIN_FREE;
    if (target_seq > ctx.ns)
        target_seq = ctx.ns;

    MYMALLOC(ord, ctx.ns);
    MYMALLOC(rank, ctx.ns);
    MYMALLOC(anchors, ctx.K);
    MYMALLOC(anchor_score, ctx.K);
    if (ord == NULL || rank == NULL ||
        anchors == NULL || anchor_score == NULL)
        mem_alloc_error("rolling phase order arrays");

    /*
     * Start from a clean assignment.  The rolling solution itself becomes
     * the warm start for the subsequent full MILP.
     */
    for (k = 0; k < ctx.ns; k++) {
        ord[k].gid = k;
        ord[k].bpos = (long long)ctx.haps[ctx.gseq[k]].bpos;
        ord[k].epos = (long long)ctx.haps[ctx.gseq[k]].epos;
        ctx.haps[ctx.gseq[k]].hap = 0;
    }

    qsort(ord, ctx.ns, sizeof(phase_ord_t), cmp_phase_ord_bpos);

    for (k = 0; k < ctx.ns; k++)
        rank[ord[k].gid] = k;
    for (k = 0; k < ctx.K; k++)
        anchors[k] = -1;

    start = 0;

    while (start < ctx.ns) {
        long long cut_pos;
        int desired_right;
        int found = 0;
        int accepted = 0;
        phase_mip_result_t local_result = PHASE_MIP_FAILED;

        block_id++;

        desired_right = choose_phase_block_end(&ctx, ord, start,
                                               target_seq, &cut_pos);
        right = desired_right;

        /*
         * First attempt: solve the desired block.  If HiGHS cannot return
         * any integer incumbent, progressively shrink the block.  This is
         * the primary failure handling for rare difficult local MILPs.
         */
        for (;;) {
            int q;
            int min_free;
            int new_right;

            min_free = PHASE_ROLL_MIN_FREE;
            if (ctx.ns - start < min_free)
                min_free = ctx.ns - start;
            if (first_block && min_free < ctx.K &&
                ctx.ns - start >= ctx.K)
                min_free = ctx.K;

            for (q = start; q < right; q++)
                ctx.haps[ctx.gseq[ord[q].gid]].hap = 0;

            if (VERBOSE > 1) {
                pthread_mutex_lock(&plock);
                fprintf(stderr,
                        "[M::%s] block=%d start=%d right=%d "
                        "free=%d anchors=%d cut=%lld\n",
                        __func__, block_id, start, right,
                        right - start, nanchor, (long long)cut_pos);
                pthread_mutex_unlock(&plock);
            }

            local_result =
                solve_phase_subproblem(&ctx, ord, start, right,
                                       anchors, nanchor, first_block);

            if (local_result != PHASE_MIP_FAILED)
                break;

            new_right =
                shrink_phase_block_end(ord, start, right, min_free);

            if (new_right <= start) {
                if (VERBOSE > 0) {
                    pthread_mutex_lock(&plock);
                    fprintf(stderr,
                            "[W::%s] block %d has no feasible MILP "
                            "incumbent even after shrinking\n",
                            __func__, block_id);
                    pthread_mutex_unlock(&plock);
                }
                goto rolling_fail;
            }

            if (VERBOSE > 0) {
                pthread_mutex_lock(&plock);
                fprintf(stderr,
                        "[W::%s] block %d failed at %d free sequences; "
                        "retrying with %d\n",
                        __func__, block_id, right - start,
                        new_right - start);
                pthread_mutex_unlock(&plock);
            }

            right = new_right;
            if (right < ctx.ns)
                cut_pos = ord[right].bpos;
            else
                cut_pos = 0x7fffffffffffffffLL;
        }

        if (local_result != PHASE_MIP_OPTIMAL)
            all_local_optimal = 0;

        /*
         * If this reaches the end there is no downstream anchor requirement.
         */
        if (right == ctx.ns) {
            found = ctx.K;
            accepted = 1;
        } else {
            found = select_phase_anchors(&ctx, ord, rank, right,
                                         cut_pos, anchors,
                                         anchor_score);
            if (found == ctx.K)
                accepted = 1;
        }

        /*
         * The first accepted block must establish all K labels.  If a
         * successfully solved block does not yet contain all K haplotypes,
         * try progressively larger blocks, up to about 2*target_seq.
         *
         * We do not let a failed larger trial destroy a previously solved
         * prefix; its current free assignments are simply cleared before
         * the next trial.
         */
        if (!accepted && first_block) {
            int max_right = start + 2 * target_seq;
            int grow_right = desired_right;
            int q;

            if (max_right > ctx.ns)
                max_right = ctx.ns;
            if (grow_right < right)
                grow_right = right;

            while (!accepted && grow_right < max_right) {
                phase_mip_result_t grow_result;

                grow_right += PHASE_ROLL_GROW_STEP;
                if (grow_right > max_right)
                    grow_right = max_right;

                while (grow_right < ctx.ns &&
                       grow_right > start &&
                       ord[grow_right - 1].bpos == ord[grow_right].bpos)
                    grow_right++;

                if (grow_right > ctx.ns)
                    grow_right = ctx.ns;

                for (q = start; q < grow_right; q++)
                    ctx.haps[ctx.gseq[ord[q].gid]].hap = 0;

                cut_pos = (grow_right < ctx.ns) ?
                          ord[grow_right].bpos :
                          0x7fffffffffffffffLL;

                if (VERBOSE > 1) {
                    pthread_mutex_lock(&plock);
                    fprintf(stderr,
                            "[M::%s] block=%d lacks all %d haplotypes; "
                            "retrying enlarged block with %d free sequences\n",
                            __func__, block_id, ctx.K, grow_right - start);
                    pthread_mutex_unlock(&plock);
                }

                grow_result =
                    solve_phase_subproblem(&ctx, ord, start, grow_right,
                                           anchors, nanchor, first_block);

                if (grow_result == PHASE_MIP_FAILED)
                    continue;

                if (grow_result != PHASE_MIP_OPTIMAL)
                    all_local_optimal = 0;

                right = grow_right;
                local_result = grow_result;

                if (right == ctx.ns) {
                    found = ctx.K;
                    accepted = 1;
                    break;
                }

                found = select_phase_anchors(&ctx, ord, rank, right,
                                             cut_pos, anchors,
                                             anchor_score);
                if (found == ctx.K)
                    accepted = 1;
            }
        }

        if (!accepted) {
            if (VERBOSE > 0) {
                pthread_mutex_lock(&plock);
                fprintf(stderr,
                        "[W::%s] unable to establish all %d carry-over "
                        "haplotypes after block %d\n",
                        __func__, ctx.K, block_id);
                pthread_mutex_unlock(&plock);
            }
            goto rolling_fail;
        }

        if (right < ctx.ns)
            decomposed = 1;

        if (VERBOSE > 1 && right < ctx.ns) {
            pthread_mutex_lock(&plock);
            fprintf(stderr,
                    "[M::%s] block=%d accepted: solved through rank %d; "
                    "carry anchors:",
                    __func__, block_id, right - 1);
            for (k = 0; k < ctx.K; k++)
                fprintf(stderr, " h%d=%d", k + 1, anchors[k]);
            fprintf(stderr, "\n");
            pthread_mutex_unlock(&plock);
        }

        start = right;
        nanchor = (start < ctx.ns) ? ctx.K : 0;
        first_block = 0;
    }

    free(ord);
    free(rank);
    free(anchors);
    free(anchor_score);

    /*
     * Only a single, unsplit, proven-optimal MILP is an exact/global result.
     * Even if every decomposed block is individually optimal, freezing
     * earlier assignments makes the overall rolling solution heuristic.
     */
    if (!decomposed && block_id == 1 && all_local_optimal)
        return PHASE_ROLL_EXACT;

    return PHASE_ROLL_APPROX;

rolling_fail:
    for (k = 0; k < ctx.ns; k++)
        ctx.haps[ctx.gseq[k]].hap = 0;

    free(ord);
    free(rank);
    free(anchors);
    free(anchor_score);

    return PHASE_ROLL_FAILED;
}



/* Objective value of a complete group-local haplotype assignment. */
static double phase_full_assignment_objective(const phase_blk_ctx_t *ctx,
                                              const int *assign)
{
    double obj = 0.0;
    int e;

    for (e = 0; e < ctx->nh; e++) {
        hlk_t *hlk = ctx->hlks + ctx->hls[e];
        int u = ctx->sidx[hlk->a];
        int v = ctx->sidx[hlk->b];

        if (u < 0 || u >= ctx->ns || v < 0 || v >= ctx->ns)
            continue;

        if (assign[u] == assign[v] &&
            assign[u] >= 1 && assign[u] <= ctx->K) {
            double w =
                hlk->l * HIC_NORM_WINDOW *
                ctx->hord[ctx->hls[e]].event /
                1000000.0 * PHASE_BETA;
            obj += w;
        }
    }

    for (e = 0; e < ctx->no; e++) {
        ovl_t *ovl = ctx->ovls + ctx->ols[e];
        int u = ctx->sidx[ovl->aread];
        int v = ctx->sidx[ovl->bread];

        if (u < 0 || u >= ctx->ns || v < 0 || v >= ctx->ns)
            continue;

        if (assign[u] == assign[v] &&
            assign[u] >= 1 && assign[u] <= ctx->K) {
            double w = (double)ovl->alen + (double)ovl->blen;
            w = w * log10(w) / 1000000.0 * PHASE_ALPHA;
            obj -= w;
        }
    }

    return obj;
}


static int phase_copy_current_assignment(const phase_blk_ctx_t *ctx,
                                         int *assign)
{
    int i;

    for (i = 0; i < ctx->ns; i++) {
        assign[i] = ctx->haps[ctx->gseq[i]].hap;
        if (assign[i] < 1 || assign[i] > ctx->K)
            return 0;
    }

    return 1;
}


/*
 * For a warm-started full MILP, avoid making the rolling solution
 * incompatible with the hard seed labels.
 *
 * We first try to use the same largest overlap clique as the original full
 * solver.  If the rolling solution places those clique members on distinct
 * haplotypes, a global label permutation makes it exactly compatible with
 * the original seed fixing.
 *
 * If they are not distinct, the original clique fixing is not merely a
 * label permutation of the incumbent.  In that rare case we use only one
 * fixed sequence as pure label-symmetry breaking, leaving the full MILP free
 * to improve the rolling partition.
 *
 * Returns:
 *   >=2 : original clique seed applied
 *    1  : one-variable symmetry seed applied
 *    0  : no seed could be applied
 */
static int phase_fix_full_seed(void *highs,
                               phase_blk_ctx_t *ctx,
                               int use_warm_start)
{
    seed_edge_t *sedges = NULL;
    uint64 *clique = NULL;
    int nsedges = 0;
    int target, clique_size = 0;
    int i, k;

#define FSEED_X(ii, kk) ((HighsInt)(ii) * ctx->K + (kk))

    MYMALLOC(sedges, ctx->no);
    MYMALLOC(clique, ctx->K);
    if ((ctx->no > 0 && sedges == NULL) || clique == NULL)
        mem_alloc_error("full MILP seed arrays");

    for (i = 0; i < ctx->no; i++) {
        ovl_t *ovl = ctx->ovls + ctx->ols[i];

        if (ovl->del)
            continue;

        sedges[nsedges].u = ctx->sidx[ovl->aread];
        sedges[nsedges].v = ctx->sidx[ovl->bread];
        sedges[nsedges].w = (double)ovl->alen + (double)ovl->blen;
        nsedges++;
    }

    qsort(sedges, nsedges, sizeof(seed_edge_t), cmp_seed_edge_desc);

    for (target = ctx->K; target >= 2 && clique_size == 0; target--)
        if (find_seed_clique(ctx->ns, sedges, nsedges,
                             target, clique))
            clique_size = target;

    if (clique_size > 0) {
        for (k = 0; k < clique_size; k++) {
            int b = (int)clique[k];
            clique[k] =
                ((uint64)ctx->haps[ctx->gseq[b]].len << 32) |
                (uint32)b;
        }
        qsort(clique, clique_size, sizeof(uint64), u64_dcmp_func);
    }

    if (use_warm_start && clique_size >= 2) {
        int *map = NULL;
        int *used_target = NULL;
        int compatible = 1;

        MYCALLOC(map, ctx->K + 1);
        MYCALLOC(used_target, ctx->K + 1);
        if (map == NULL || used_target == NULL)
            mem_alloc_error("full MILP warm label map");

        for (k = 0; k < clique_size; k++) {
            int b = (uint32)clique[k];
            int old_hap = ctx->haps[ctx->gseq[b]].hap;
            int target_hap = k + 1;

            if (old_hap < 1 || old_hap > ctx->K ||
                map[old_hap] != 0) {
                compatible = 0;
                break;
            }

            map[old_hap] = target_hap;
            used_target[target_hap] = 1;
        }

        if (compatible) {
            int old_hap, target_hap = 1;

            for (old_hap = 1; old_hap <= ctx->K; old_hap++) {
                if (map[old_hap] != 0)
                    continue;

                while (target_hap <= ctx->K &&
                       used_target[target_hap])
                    target_hap++;

                if (target_hap > ctx->K) {
                    compatible = 0;
                    break;
                }

                map[old_hap] = target_hap;
                used_target[target_hap] = 1;
            }
        }

        if (compatible) {
            /*
             * Apply the global label permutation to the rolling incumbent.
             * The objective is label-symmetric, so its value is unchanged.
             */
            for (i = 0; i < ctx->ns; i++) {
                int old_hap = ctx->haps[ctx->gseq[i]].hap;
                if (old_hap >= 1 && old_hap <= ctx->K)
                    ctx->haps[ctx->gseq[i]].hap = map[old_hap];
            }

            for (k = 0; k < clique_size; k++) {
                int b = (uint32)clique[k];

                Highs_changeColBounds(highs, FSEED_X(b, k),
                                      1.0, 1.0);
                for (i = 0; i < ctx->K; i++) {
                    if (i == k)
                        continue;
                    Highs_changeColBounds(highs, FSEED_X(b, i),
                                          0.0, 0.0);
                }
            }

            free(map);
            free(used_target);
            free(sedges);
            free(clique);
/* #undef FSEED_X deferred to function end */
            return clique_size;
        }

        free(map);
        free(used_target);
    }

    if (!use_warm_start && clique_size >= 2) {
        for (k = 0; k < clique_size; k++) {
            int b = (uint32)clique[k];

            Highs_changeColBounds(highs, FSEED_X(b, k), 1.0, 1.0);
            for (i = 0; i < ctx->K; i++) {
                if (i == k)
                    continue;
                Highs_changeColBounds(highs, FSEED_X(b, i),
                                      0.0, 0.0);
            }
        }

        free(sedges);
        free(clique);
/* #undef FSEED_X deferred to function end */
        return clique_size;
    }

    /*
     * Warm-start fallback: one fixed sequence is enough to remove the
     * global label permutation symmetry and can always be made compatible
     * with the rolling incumbent by swapping haplotype labels.
     */
    if (use_warm_start && ctx->ns > 0) {
        int old_hap = ctx->haps[ctx->gseq[0]].hap;

        if (old_hap >= 1 && old_hap <= ctx->K) {
            if (old_hap != 1) {
                for (i = 0; i < ctx->ns; i++) {
                    int h = ctx->haps[ctx->gseq[i]].hap;

                    if (h == 1)
                        ctx->haps[ctx->gseq[i]].hap = old_hap;
                    else if (h == old_hap)
                        ctx->haps[ctx->gseq[i]].hap = 1;
                }
            }

            Highs_changeColBounds(highs, FSEED_X(0, 0), 1.0, 1.0);
            for (k = 1; k < ctx->K; k++)
                Highs_changeColBounds(highs, FSEED_X(0, k),
                                      0.0, 0.0);

            if (VERBOSE > 1) {
                pthread_mutex_lock(&plock);
                fprintf(stderr,
                        "[M::%s] warm start incompatible with full seed "
                        "clique; using one-variable symmetry anchor\n",
                        __func__);
                pthread_mutex_unlock(&plock);
            }

            free(sedges);
            free(clique);
/* #undef FSEED_X deferred to function end */
            return 1;
        }
    }

    free(sedges);
    free(clique);
#undef FSEED_X
    return 0;
}


/*
 * Build a complete primal column solution for the full model.
 * Supplying x, z and o makes the MIP start immediately feasible.
 */
static int phase_fill_full_mip_start(const phase_blk_ctx_t *ctx,
                                     HighsInt offset_x,
                                     HighsInt offset_z,
                                     HighsInt offset_o,
                                     double *start_value)
{
    int i, e;

#define FSTART_X(ii, kk) (offset_x + (HighsInt)(ii) * ctx->K + (kk))
#define FSTART_Z(ee, kk) (offset_z + (HighsInt)(ee) * ctx->K + (kk))
#define FSTART_O(ee)     (offset_o + (ee))

    for (i = 0; i < ctx->ns; i++) {
        int h = ctx->haps[ctx->gseq[i]].hap;

        if (h < 1 || h > ctx->K) {
/* #undef FSTART_X deferred to function end */
/* #undef FSTART_Z deferred to function end */
/* #undef FSTART_O deferred to function end */
            return 0;
        }

        start_value[FSTART_X(i, h - 1)] = 1.0;
    }

    for (e = 0; e < ctx->nh; e++) {
        hlk_t *hlk = ctx->hlks + ctx->hls[e];
        int u = ctx->sidx[hlk->a];
        int v = ctx->sidx[hlk->b];
        int hu, hv;

        if (u < 0 || u >= ctx->ns || v < 0 || v >= ctx->ns)
            continue;

        hu = ctx->haps[ctx->gseq[u]].hap;
        hv = ctx->haps[ctx->gseq[v]].hap;

        if (hu == hv && hu >= 1 && hu <= ctx->K)
            start_value[FSTART_Z(e, hu - 1)] = 1.0;
    }

    for (e = 0; e < ctx->no; e++) {
        ovl_t *ovl = ctx->ovls + ctx->ols[e];
        int u = ctx->sidx[ovl->aread];
        int v = ctx->sidx[ovl->bread];
        int hu, hv;

        if (u < 0 || u >= ctx->ns || v < 0 || v >= ctx->ns)
            continue;

        hu = ctx->haps[ctx->gseq[u]].hap;
        hv = ctx->haps[ctx->gseq[v]].hap;

        if (hu == hv && hu >= 1 && hu <= ctx->K)
            start_value[FSTART_O(e)] = 1.0;
    }

#undef FSTART_X
#undef FSTART_Z
#undef FSTART_O
    return 1;
}



/*
 * Full MILP refinement.
 *
 * If use_warm_start != 0, the current haps[].hap assignment is supplied as
 * a complete feasible MIP start.  The function only overwrites it if HiGHS
 * returns a feasible assignment whose full-model objective is not worse.
 */
static phase_mip_result_t hybrid_phase_full_ilp_solve(void *_data,
                                                       long jid,
                                                       int tid,
                                                       int use_warm_start)
{
    phase_data_t *data = (phase_data_t *)_data;
    phase_blk_ctx_t ctx;
    grp_info_t *grps = data->grps;
    hap_info_t *haps = data->haps;
    ovl_t *ovls = data->ovls;
    hlk_t *hlks = data->hlks;
    ord_dbl_t *hord = data->hord;
    int K = data->K;
    int ns = data->ns[jid];
    int nh = data->nh[jid];
    int no = data->no[jid];
    int *sidx = data->sidx;
    int *ols = data->ols[jid];
    int *hls = data->hls[jid];
    int *gseq = grps[jid].s;

    ovl_t *ovl;
    hlk_t *hlk;
    mip_callback_t _fallback, *fallback = &_fallback;
    mip_track_t *track;
    HighsInt offset_x, offset_z, offset_o, tot_vars;
    HighsInt i, e, k, u, v, col_idx, num_cols, *idxs;
    HighsInt run_status, model_status;
    void *highs;
    double w, inf, *col_costs, *vals;
    double *mip_start = NULL;
    int *warm_assign = NULL;
    double warm_obj = -1e300;
    int wd;
    phase_mip_result_t result = PHASE_MIP_FAILED;

    (void)tid;

    memset(&ctx, 0, sizeof(ctx));
    ctx.data = data;
    ctx.jid = jid;
    ctx.grps = grps;
    ctx.haps = haps;
    ctx.ovls = ovls;
    ctx.hlks = hlks;
    ctx.hord = hord;
    ctx.K = K;
    ctx.ns = ns;
    ctx.nh = nh;
    ctx.no = no;
    ctx.sidx = sidx;
    ctx.ols = ols;
    ctx.hls = hls;
    ctx.gseq = gseq;

    if (use_warm_start) {
        MYMALLOC(warm_assign, ns);
        if (warm_assign == NULL)
            mem_alloc_error("full MILP warm assignment");

        if (!phase_copy_current_assignment(&ctx, warm_assign)) {
            free(warm_assign);
            warm_assign = NULL;
            use_warm_start = 0;
        }
    }

    highs = Highs_create();
    inf = Highs_getInfinity(highs);

    Highs_setBoolOptionValue(highs, "output_flag", 0);
    Highs_setBoolOptionValue(highs, "presolve", 1);

    wd = (K > 3) ? K : 3;
    MYMALLOC(idxs, wd);
    MYMALLOC(vals, wd);
    if (idxs == NULL || vals == NULL)
        mem_alloc_error("full phase MILP row arrays");

    offset_x = 0;
    offset_z = offset_x + (HighsInt)ns * K;
    offset_o = offset_z + (HighsInt)nh * K;
    tot_vars = offset_o + no;

#define IDX_X(i, k) (offset_x + (HighsInt)(i) * K + (k))
#define IDX_Z(e, k) (offset_z + (HighsInt)(e) * K + (k))
#define IDX_O(e)    (offset_o + (e))

    /* x_ik */
    for (i = 0; i < ns; i++) {
        for (k = 0; k < K; k++) {
            col_idx = IDX_X(i, k);
            Highs_addCol(highs, 0.0, 0.0, 1.0, 0, NULL, NULL);
            Highs_changeColIntegrality(highs, col_idx,
                                       kHighsVarTypeInteger);
        }
    }

    /* z_e,k */
    for (e = 0; e < nh; e++) {
        for (k = 0; k < K; k++) {
            col_idx = IDX_Z(e, k);
            Highs_addCol(highs, 0.0, 0.0, 1.0, 0, NULL, NULL);
            Highs_changeColIntegrality(highs, col_idx,
                                       kHighsVarTypeContinuous);
        }
    }

    /* o_e */
    for (e = 0; e < no; e++) {
        col_idx = IDX_O(e);
        Highs_addCol(highs, 0.0, 0.0, 1.0, 0, NULL, NULL);
        Highs_changeColIntegrality(highs, col_idx,
                                   kHighsVarTypeContinuous);
    }

    /*
     * Seed fixing.
     *
     * With no warm start this reproduces the original largest-clique
     * behaviour.  With a warm start, labels are permuted to make the same
     * seed compatible whenever possible; otherwise a single pure-symmetry
     * anchor is used.
     */
    {
        int seed_size = phase_fix_full_seed(highs, &ctx, use_warm_start);

        /*
         * Preserve the original special-case behaviour when there is
         * essentially no usable overlap graph and no rolling incumbent.
         */
        if (!use_warm_start && seed_size < 2) {
            for (i = 0; i < ns; i++)
                haps[gseq[i]].hap = 1;

            free(idxs);
            free(vals);
            free(warm_assign);
            Highs_destroy(highs);

/* #undef IDX_X deferred to function end */
/* #undef IDX_Z deferred to function end */
/* #undef IDX_O deferred to function end */
            return PHASE_MIP_OPTIMAL;
        }
    }

    /* exactly one haplotype per sequence */
    for (i = 0; i < ns; i++) {
        for (k = 0; k < K; k++) {
            idxs[k] = IDX_X(i, k);
            vals[k] = 1.0;
        }
        Highs_addRow(highs, 1.0, 1.0, K, idxs, vals);
    }

    /* Hi-C linearisation */
    for (e = 0; e < nh; e++) {
        hlk = hlks + hls[e];
        u = sidx[hlk->a];
        v = sidx[hlk->b];

        for (k = 0; k < K; k++) {
            idxs[0] = IDX_Z(e, k); vals[0] =  1.0;

            idxs[1] = IDX_X(u, k); vals[1] = -1.0;
            Highs_addRow(highs, -inf, 0.0, 2, idxs, vals);

            idxs[1] = IDX_X(v, k); vals[1] = -1.0;
            Highs_addRow(highs, -inf, 0.0, 2, idxs, vals);
        }
    }

    /* OVL linearisation */
    for (e = 0; e < no; e++) {
        ovl = ovls + ols[e];
        u = sidx[ovl->aread];
        v = sidx[ovl->bread];

        idxs[0] = IDX_O(e); vals[0] = 1.0;

        for (k = 0; k < K; k++) {
            idxs[1] = IDX_X(u, k); vals[1] = -1.0;
            idxs[2] = IDX_X(v, k); vals[2] = -1.0;
            Highs_addRow(highs, -1.0, inf, 3, idxs, vals);
        }
    }

    free(idxs);
    free(vals);

    num_cols = Highs_getNumCol(highs);
    MYCALLOC(col_costs, num_cols);
    if (col_costs == NULL)
        mem_alloc_error("full phasing objective");

    Highs_changeObjectiveSense(highs, -1);

    for (e = 0; e < nh; e++) {
        w = hlks[hls[e]].l * HIC_NORM_WINDOW *
            hord[hls[e]].event / 1000000.0 * PHASE_BETA;

        for (k = 0; k < K; k++)
            col_costs[IDX_Z(e, k)] = +w;
    }

    for (e = 0; e < no; e++) {
        w = (double)ovls[ols[e]].alen + ovls[ols[e]].blen;
        w = w * log10(w) / 1000000.0 * PHASE_ALPHA;
        col_costs[IDX_O(e)] = -w;
    }

    Highs_changeColsCostByRange(highs, 0, num_cols - 1, col_costs);
    free(col_costs);

    /*
     * Set the complete rolling assignment as a MIP start.
     * phase_fix_full_seed() may have globally permuted haplotype labels;
     * rebuild warm_assign after that permutation.
     */
    if (use_warm_start) {
        free(warm_assign);
        warm_assign = NULL;

        MYMALLOC(warm_assign, ns);
        if (warm_assign == NULL)
            mem_alloc_error("full MILP relabelled warm assignment");

        if (phase_copy_current_assignment(&ctx, warm_assign)) {
            HighsInt start_status;

            warm_obj =
                phase_full_assignment_objective(&ctx, warm_assign);

            MYCALLOC(mip_start, num_cols);
            if (mip_start == NULL)
                mem_alloc_error("full MILP start vector");

            if (phase_fill_full_mip_start(&ctx,
                                          offset_x, offset_z, offset_o,
                                          mip_start)) {
                start_status =
                    Highs_setSolution(highs, mip_start,
                                      NULL, NULL, NULL);

                if (VERBOSE > 1) {
                    pthread_mutex_lock(&plock);
                    fprintf(stderr,
                            "[M::%s] full-MILP warm start: "
                            "status=%lld objective=%.6f\n",
                            __func__, (long long)start_status, warm_obj);
                    pthread_mutex_unlock(&plock);
                }
            }

            free(mip_start);
            mip_start = NULL;
        }
    }

    Highs_setDoubleOptionValue(highs, "mip_rel_gap",
                               HIGHS_MIP_REL_GAP);
    if (HIGHS_TIME_LIMIT > 0)
        Highs_setDoubleOptionValue(highs, "time_limit",
                                   HIGHS_TIME_LIMIT);

    if (VERBOSE > 1) {
        pthread_mutex_lock(&plock);
        fprintf(stderr,
                "[M::%s] solving full MILP with %lld variables "
                "and %lld constraints%s\n",
                __func__,
                (long long)tot_vars,
                (long long)Highs_getNumRow(highs),
                use_warm_start ? " (rolling warm start)" : "");
        fprintf(stderr, "[M::%s]    #haps: %d\n", __func__, K);
        fprintf(stderr, "[M::%s]    #seqs: %d\n", __func__, ns);
        fprintf(stderr, "[M::%s]    #ovls: %d\n", __func__, no);
        fprintf(stderr, "[M::%s]    #hlks: %d\n", __func__, nh);
        pthread_mutex_unlock(&plock);
    }

    fallback->front = fallback->count = 0;
    Highs_setCallback(highs, mip_callback, fallback);
    Highs_startCallback(highs, kHighsCallbackMipImprovingSolution);

    run_status = Highs_run(highs);
    model_status = Highs_getModelStatus(highs);

    if (VERBOSE > 1) {
        double time, primal, dual = 0.0;

        pthread_mutex_lock(&plock);
        fprintf(stderr,
                "[M::%s] full HiGHS finished: "
                "run_status=%lld model_status=%lld\n",
                __func__,
                (long long)run_status, (long long)model_status);

        for (i = 0; (track = mf_at(fallback, i)) != NULL; i++)
            fprintf(stderr,
                    "[M::%s]    callback %3lld: "
                    "time=%8.2f primal=%12.6f dual=%12.6f\n",
                    __func__, (long long)i + 1,
                    track->running_time,
                    track->mip_primal_bound,
                    track->mip_dual_bound);

        time = Highs_getRunTime(highs);
        primal = Highs_getObjectiveValue(highs);
        Highs_getDoubleInfoValue(highs, "mip_dual_bound", &dual);

        fprintf(stderr,
                "[M::%s]    final: time=%8.2f "
                "primal=%12.6f dual=%12.6f\n",
                __func__, time, primal, dual);
        pthread_mutex_unlock(&plock);
    }

    if (run_status >= 0) {
        HighsInt num_rows = Highs_getNumRow(highs);
        double *col_value = NULL, *col_dual = NULL;
        double *row_value = NULL, *row_dual = NULL;
        int *candidate = NULL;
        int candidate_feasible = 1;
        double candidate_obj = -1e300;

        MYMALLOC(col_value, num_cols);
        MYMALLOC(col_dual, num_cols);
        MYMALLOC(row_value, num_rows);
        MYMALLOC(row_dual, num_rows);
        MYMALLOC(candidate, ns);

        if (col_value == NULL || col_dual == NULL ||
            row_value == NULL || row_dual == NULL ||
            candidate == NULL)
            mem_alloc_error("full phasing solution arrays");

        Highs_getSolution(highs, col_value, col_dual,
                          row_value, row_dual);

        for (i = 0; i < ns; i++) {
            int nk = 0;
            int h = 0;

            for (k = 0; k < K; k++) {
                if (col_value[IDX_X(i, k)] > 0.5) {
                    nk++;
                    h = (int)k + 1;
                }
            }

            if (nk != 1) {
                candidate_feasible = 0;
                break;
            }

            candidate[i] = h;
        }

        if (candidate_feasible) {
            candidate_obj =
                phase_full_assignment_objective(&ctx, candidate);

            /*
             * With a warm start, never replace the rolling solution by a
             * worse incumbent.  Normally HiGHS itself guarantees this once
             * the start has been accepted, but this guard is inexpensive.
             */
            if (!use_warm_start ||
                warm_obj <= -1e299 ||
                candidate_obj + 1e-9 >= warm_obj) {
                for (i = 0; i < ns; i++)
                    haps[gseq[i]].hap = candidate[i];

                if (model_status == kHighsModelStatusOptimal)
                    result = PHASE_MIP_OPTIMAL;
                else
                    result = PHASE_MIP_FEASIBLE;
            } else {
                result = PHASE_MIP_FEASIBLE;

                if (VERBOSE > 0) {
                    pthread_mutex_lock(&plock);
                    fprintf(stderr,
                            "[W::%s] final full-MILP incumbent %.6f "
                            "is worse than rolling warm start %.6f; "
                            "keeping rolling assignment\n",
                            __func__, candidate_obj, warm_obj);
                    pthread_mutex_unlock(&plock);
                }
            }
        }

        free(col_value);
        free(col_dual);
        free(row_value);
        free(row_dual);
        free(candidate);
    }

    free(warm_assign);
    Highs_destroy(highs);

#undef IDX_X
#undef IDX_Z
#undef IDX_O

    return result;
}



/*
 * Two-stage phasing entry point.
 *
 * 1. Run rolling-horizon MILP.
 * 2. If it was a single unsplit model and was proved optimal, stop.
 * 3. Otherwise use the rolling solution as a warm start for the original
 *    full MILP and allow HiGHS to improve it.
 * 4. If rolling itself fails, try the full MILP without a warm start.
 *
 * Final status:
 *   PHASE_ILP_OPTIMUM  = global/full MILP proved optimal
 *   PHASE_ILP_FEASIBLE = valid incumbent, but no global proof
 *   anything else      = no MILP solution; caller may use heuristic fallback
 */
static void hybrid_phase_ilp_thread(void *_data, long jid, int tid)
{
    phase_data_t *data = (phase_data_t *)_data;
    phase_roll_result_t roll_result;
    phase_mip_result_t full_result;

    if (PHASE_ILP_IS_SOLVED(data->status[jid]))
        return;
    
    roll_result = hybrid_phase_rolling_solve(_data, jid, tid);

    if (roll_result == PHASE_ROLL_EXACT) {
        data->status[jid] = PHASE_ILP_OPTIMUM;

        if (VERBOSE > 1) {
            pthread_mutex_lock(&plock);
            fprintf(stderr,
                    "[M::%s] rolling stage used one unsplit MILP and "
                    "proved it optimal; full refinement skipped\n",
                    __func__);
            pthread_mutex_unlock(&plock);
        }
        return;
    }

    if (roll_result == PHASE_ROLL_APPROX) {
        if (VERBOSE > 1) {
            pthread_mutex_lock(&plock);
            fprintf(stderr,
                    "[M::%s] rolling stage produced a suboptimal/"
                    "decomposed incumbent; starting full-MILP refinement\n",
                    __func__);
            pthread_mutex_unlock(&plock);
        }

        full_result =
            hybrid_phase_full_ilp_solve(_data, jid, tid, 1);

        if (full_result == PHASE_MIP_OPTIMAL) {
            data->status[jid] = PHASE_ILP_OPTIMUM;
            return;
        }

        /*
         * Even if the full refinement fails to return a new incumbent, the
         * rolling assignment remains in haps[].hap and is a valid solution.
         */
        data->status[jid] = PHASE_ILP_FEASIBLE;
        return;
    }

    /*
     * Rolling failed even after adaptive shrinking.  Try the original full
     * MILP once without a warm start.  If this also fails, leave status
     * unchanged so the existing heuristic fallback can take over.
     */
    if (VERBOSE > 0) {
        pthread_mutex_lock(&plock);
        fprintf(stderr,
                "[W::%s] rolling MILP failed; trying full MILP "
                "without warm start\n",
                __func__);
        pthread_mutex_unlock(&plock);
    }

    full_result =
        hybrid_phase_full_ilp_solve(_data, jid, tid, 0);

    if (full_result == PHASE_MIP_OPTIMAL)
        data->status[jid] = PHASE_ILP_OPTIMUM;
    else if (full_result == PHASE_MIP_FEASIBLE)
        data->status[jid] = PHASE_ILP_FEASIBLE;
    else if (VERBOSE > 0) {
        pthread_mutex_lock(&plock);
        fprintf(stderr,
                "[W::%s] both rolling and full MILP failed; "
                "heuristic branch should be tried for this group\n",
                __func__);
        pthread_mutex_unlock(&plock);
    }
}

/*
 * Fallback heuristic for hybrid haplotype phasing.
 *
 * Objective (same signed pairwise objective as the MILP):
 *
 *   sum_HiC w_uv [h_u == h_v] - sum_OVL w_uv [h_u == h_v]
 *
 * Heuristic:
 *   - same overlap-clique anchors as the MILP;
 *   - deterministic + randomized greedy construction;
 *   - exact best-improvement 1-vertex local search;
 *   - repeated perturbation of low-confidence vertices;
 *   - multiple starts;
 *   - any complete assignment left by ILP is used as an extra start.
 *
 */

#ifndef PHASE_HEU_NSTART
#define PHASE_HEU_NSTART 8
#endif

#ifndef PHASE_HEU_NPERTURB
#define PHASE_HEU_NPERTURB 12
#endif

#ifndef PHASE_HEU_PERTURB_FRAC
#define PHASE_HEU_PERTURB_FRAC 0.08
#endif

#ifndef PHASE_HEU_PERTURB_MAX
#define PHASE_HEU_PERTURB_MAX 12
#endif

#ifndef PHASE_HEU_RANDOM_SLACK
#define PHASE_HEU_RANDOM_SLACK 0.05
#endif

#ifndef PHASE_HEU_GAIN_EPS
#define PHASE_HEU_GAIN_EPS 1e-12
#endif

#ifndef PHASE_HEU_MAX_MOVES_FACTOR
#define PHASE_HEU_MAX_MOVES_FACTOR 1000
#endif

#ifndef PHASE_HEU_INCLUDE_DELETED_OVL
#define PHASE_HEU_INCLUDE_DELETED_OVL 1
#endif

#define PHASE_HEU_IS_SOLVED(s) \
    (PHASE_ILP_IS_SOLVED(s) || (s) == PHASE_HEU_APPROX)

typedef struct {
    int u, v;
    double w;
} phase_heu_edge_t;

typedef struct {
    int to, next;
    double w;
} phase_heu_arc_t;

typedef struct {
    int ns, K;
    int nedge, narc;
    phase_heu_edge_t *edge;
    phase_heu_arc_t *arc;
    int *head;
    double *abs_degree;
} phase_heu_graph_t;

typedef struct {
    int v;
    double margin;
} phase_heu_conf_t;


static uint64_t phase_heu_rng_next(uint64_t *s)
{
    uint64_t x = *s;
    if (x == 0) x = UINT64_C(0x9e3779b97f4a7c15);
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *s = x;
    return x * UINT64_C(2685821657736338717);
}

static double phase_heu_rng01(uint64_t *s)
{
    return (phase_heu_rng_next(s) >> 11) *
           (1.0 / 9007199254740992.0);
}

static int phase_heu_rand_int(uint64_t *s, int n)
{
    return n <= 1 ? 0 :
        (int)(phase_heu_rng_next(s) % (uint64_t)n);
}

static void phase_heu_graph_destroy(phase_heu_graph_t *g)
{
    free(g->edge);
    free(g->arc);
    free(g->head);
    free(g->abs_degree);
    memset(g, 0, sizeof(*g));
}

static void phase_heu_graph_add(phase_heu_graph_t *g,
                                int u, int v, double w)
{
    int e;

    if (u < 0 || u >= g->ns || v < 0 || v >= g->ns)
        return;

    e = g->nedge++;
    g->edge[e].u = u;
    g->edge[e].v = v;
    g->edge[e].w = w;

    /* self-edge is an objective constant; omit from move-gain adjacency */
    if (u == v) return;

    g->arc[g->narc].to = v;
    g->arc[g->narc].w = w;
    g->arc[g->narc].next = g->head[u];
    g->head[u] = g->narc++;

    g->arc[g->narc].to = u;
    g->arc[g->narc].w = w;
    g->arc[g->narc].next = g->head[v];
    g->head[v] = g->narc++;

    g->abs_degree[u] += fabs(w);
    g->abs_degree[v] += fabs(w);
}

static int phase_heu_graph_build(phase_heu_graph_t *g,
                                 phase_data_t *data, long jid)
{
    ovl_t *ovls = data->ovls;
    hlk_t *hlks = data->hlks;
    ord_dbl_t *hord = data->hord;
    int ns = data->ns[jid], nh = data->nh[jid], no = data->no[jid];
    int *sidx = data->sidx, *ols = data->ols[jid], *hls = data->hls[jid];
    int maxe = nh + no;
    int e, u, v;
    double w;

    memset(g, 0, sizeof(*g));
    g->ns = ns;
    g->K = data->K;

    g->edge = (phase_heu_edge_t *)
        malloc((size_t)(maxe > 0 ? maxe : 1) * sizeof(*g->edge));
    g->arc = (phase_heu_arc_t *)
        malloc((size_t)(2 * maxe > 0 ? 2 * maxe : 1) * sizeof(*g->arc));
    g->head = (int *)malloc((size_t)ns * sizeof(*g->head));
    g->abs_degree = (double *)calloc((size_t)ns, sizeof(*g->abs_degree));

    if (!g->edge || !g->arc || !g->head || !g->abs_degree) {
        phase_heu_graph_destroy(g);
        return 0;
    }

    for (u = 0; u < ns; ++u) g->head[u] = -1;

    for (e = 0; e < nh; ++e) {
        hlk_t *h = hlks + hls[e];
        u = sidx[h->a];
        v = sidx[h->b];
        w = h->l * HIC_NORM_WINDOW *
            hord[hls[e]].event / 1000000.0 * PHASE_BETA;
        phase_heu_graph_add(g, u, v, +w);
    }

    for (e = 0; e < no; ++e) {
        ovl_t *o = ovls + ols[e];

#if !PHASE_HEU_INCLUDE_DELETED_OVL
        if (o->del) continue;
#endif

        u = sidx[o->aread];
        v = sidx[o->bread];
        w = (double)o->alen + (double)o->blen;
        w = w * log10(w) / 1000000.0 * PHASE_ALPHA;
        phase_heu_graph_add(g, u, v, -w);
    }

    return 1;
}

static double phase_heu_objective(const phase_heu_graph_t *g,
                                  const int *hap)
{
    double z = 0.0;
    int e;

    for (e = 0; e < g->nedge; ++e)
        if (hap[g->edge[e].u] == hap[g->edge[e].v])
            z += g->edge[e].w;

    return z;
}

static void phase_heu_build_scores(const phase_heu_graph_t *g,
                                   const int *hap, double *score)
{
    int i, a;

    memset(score, 0, (size_t)g->ns * g->K * sizeof(*score));

    for (i = 0; i < g->ns; ++i) {
        for (a = g->head[i]; a >= 0; a = g->arc[a].next) {
            int j = g->arc[a].to;
            int k = hap[j];
            if (k >= 0 && k < g->K)
                score[(size_t)i * g->K + k] += g->arc[a].w;
        }
    }
}

static int phase_heu_make_seeds(phase_data_t *data, long jid,
                                const phase_heu_graph_t *g,
                                unsigned char *fixed, int *seed_hap)
{
    grp_info_t *grps = data->grps;
    hap_info_t *haps = data->haps;
    ovl_t *ovls = data->ovls;
    int ns = data->ns[jid], no = data->no[jid], K = data->K;
    int *sidx = data->sidx, *ols = data->ols[jid], *gseq = grps[jid].s;
    seed_edge_t *sedges = NULL;
    uint64 *clique = NULL;
    int nsedges = 0, clique_size = 0;
    int i, k, target;

    memset(fixed, 0, (size_t)ns * sizeof(*fixed));
    for (i = 0; i < ns; ++i) seed_hap[i] = -1;

    if (no > 0) {
        sedges = (seed_edge_t *)malloc((size_t)no * sizeof(*sedges));
        clique = (uint64 *)malloc((size_t)K * sizeof(*clique));
        if (!sedges || !clique) {
            free(sedges);
            free(clique);
            return 0;
        }

        for (i = 0; i < no; ++i) {
            ovl_t *o = ovls + ols[i];
            if (o->del) continue;  /* exactly as full-ILP seed construction */
            sedges[nsedges].u = sidx[o->aread];
            sedges[nsedges].v = sidx[o->bread];
            sedges[nsedges].w = (double)o->alen + (double)o->blen;
            ++nsedges;
        }

        qsort(sedges, (size_t)nsedges, sizeof(*sedges),
              cmp_seed_edge_desc);

        for (target = K; target >= 2 && clique_size == 0; --target)
            if (find_seed_clique(ns, sedges, nsedges, target, clique))
                clique_size = target;

        if (clique_size > 0) {
            for (k = 0; k < clique_size; ++k) {
                int b = (int)clique[k];
                clique[k] =
                    ((uint64)haps[gseq[b]].len << 32) |
                    (uint64)(uint32_t)b;
            }

            qsort(clique, (size_t)clique_size, sizeof(*clique),
                  u64_dcmp_func);

            for (k = 0; k < clique_size; ++k) {
                int b = (int)(uint32_t)clique[k];
                fixed[b] = 1;
                seed_hap[b] = k;
            }
        }
    }

    free(sedges);
    free(clique);

    /*
     * If there is no clique, fix one informative vertex to hap 0.
     * This is only label-symmetry breaking.
     */
    if (clique_size == 0 && ns > 0) {
        int best = 0;
        for (i = 1; i < ns; ++i)
            if (g->abs_degree[i] > g->abs_degree[best])
                best = i;

        fixed[best] = 1;
        seed_hap[best] = 0;
        clique_size = 1;
    }

    return clique_size;
}

/*
 * Relabel a complete incoming assignment to satisfy the seed labels without
 * changing its objective.  Return 0 if the incoming partition is genuinely
 * incompatible with the fixed seed partition.
 */
static int phase_heu_align_to_seeds(int ns, int K, int *hap,
                                    const unsigned char *fixed,
                                    const int *seed_hap)
{
    int *map = (int *)malloc((size_t)K * sizeof(*map));
    unsigned char *used =
        (unsigned char *)calloc((size_t)K, sizeof(*used));
    int i, a, b;

    if (!map || !used) {
        free(map);
        free(used);
        return 0;
    }

    for (a = 0; a < K; ++a) map[a] = -1;

    for (i = 0; i < ns; ++i) {
        if (!fixed[i]) continue;

        a = hap[i];
        b = seed_hap[i];

        if (a < 0 || a >= K || b < 0 || b >= K ||
            (map[a] >= 0 && map[a] != b) ||
            (map[a] < 0 && used[b])) {
            free(map);
            free(used);
            return 0;
        }

        map[a] = b;
        used[b] = 1;
    }

    for (a = 0; a < K; ++a) {
        if (map[a] >= 0) continue;
        for (b = 0; b < K && used[b]; ++b) {}
        if (b == K) {
            free(map);
            free(used);
            return 0;
        }
        map[a] = b;
        used[b] = 1;
    }

    for (i = 0; i < ns; ++i) {
        if (hap[i] < 0 || hap[i] >= K) {
            free(map);
            free(used);
            return 0;
        }
        hap[i] = map[hap[i]];
    }

    free(map);
    free(used);
    return 1;
}

static int phase_heu_choose_hap(const double *s, int K, double conn,
                                int randomized, uint64_t *rng)
{
    int k, bestk = 0;
    double best = s[0];

    for (k = 1; k < K; ++k) {
        if (s[k] > best) {
            best = s[k];
            bestk = k;
        }
    }

    if (!randomized) return bestk;

    {
        double slack = PHASE_HEU_RANDOM_SLACK * (1.0 + conn);
        int ncand = 0, chosen = bestk;

        for (k = 0; k < K; ++k) {
            if (s[k] >= best - slack) {
                ++ncand;
                if (phase_heu_rand_int(rng, ncand) == 0)
                    chosen = k;
            }
        }
        return chosen;
    }
}

static int phase_heu_greedy(const phase_heu_graph_t *g,
                            const unsigned char *fixed,
                            const int *seed_hap,
                            int randomized, uint64_t *rng,
                            int *hap)
{
    double *score =
        (double *)calloc((size_t)g->ns * g->K, sizeof(*score));
    double *conn =
        (double *)calloc((size_t)g->ns, sizeof(*conn));
    unsigned char *assigned =
        (unsigned char *)calloc((size_t)g->ns, sizeof(*assigned));
    int i, k, a, nassigned = 0;

    if (!score || !conn || !assigned) {
        free(score);
        free(conn);
        free(assigned);
        return 0;
    }

    for (i = 0; i < g->ns; ++i) hap[i] = -1;

    for (i = 0; i < g->ns; ++i) {
        if (!fixed[i]) continue;
        hap[i] = seed_hap[i];
        assigned[i] = 1;
        ++nassigned;
    }

    /* contributions from seed vertices */
    for (i = 0; i < g->ns; ++i) {
        if (!assigned[i]) continue;
        k = hap[i];

        for (a = g->head[i]; a >= 0; a = g->arc[a].next) {
            int j = g->arc[a].to;
            if (assigned[j]) continue;
            score[(size_t)j * g->K + k] += g->arc[a].w;
            conn[j] += fabs(g->arc[a].w);
        }
    }

    while (nassigned < g->ns) {
        int best_i = -1;
        double best_priority = -1.0;

        for (i = 0; i < g->ns; ++i) {
            double p;
            if (assigned[i]) continue;

            p = conn[i] > 0.0 ?
                conn[i] : 1e-12 * (1.0 + g->abs_degree[i]);

            if (randomized)
                p *= 0.85 + 0.30 * phase_heu_rng01(rng);

            if (best_i < 0 || p > best_priority) {
                best_i = i;
                best_priority = p;
            }
        }

        if (best_i < 0) break;

        k = phase_heu_choose_hap(
            score + (size_t)best_i * g->K,
            g->K, conn[best_i], randomized, rng);

        hap[best_i] = k;
        assigned[best_i] = 1;
        ++nassigned;

        for (a = g->head[best_i]; a >= 0; a = g->arc[a].next) {
            int j = g->arc[a].to;
            if (assigned[j]) continue;
            score[(size_t)j * g->K + k] += g->arc[a].w;
            conn[j] += fabs(g->arc[a].w);
        }
    }

    free(score);
    free(conn);
    free(assigned);

    return nassigned == g->ns;
}

/*
 * score[i,k] = sum of signed edge weights from i to current hap k.
 * Exact gain of i : old -> k is score[i,k] - score[i,old].
 */
static double phase_heu_local_search(const phase_heu_graph_t *g,
                                     const unsigned char *fixed,
                                     int *hap, double *score,
                                     int *nmove_out)
{
    double obj;
    int nmove = 0;
    int max_moves = PHASE_HEU_MAX_MOVES_FACTOR * g->ns;

    if (max_moves < 1000) max_moves = 1000;

    phase_heu_build_scores(g, hap, score);
    obj = phase_heu_objective(g, hap);

    for (;;) {
        int best_i = -1, best_k = -1;
        double best_gain = PHASE_HEU_GAIN_EPS;
        int i, k;

        for (i = 0; i < g->ns; ++i) {
            int oldk;
            double oldscore;

            if (fixed[i]) continue;

            oldk = hap[i];
            oldscore = score[(size_t)i * g->K + oldk];

            for (k = 0; k < g->K; ++k) {
                double gain;
                if (k == oldk) continue;
                gain = score[(size_t)i * g->K + k] - oldscore;

                if (gain > best_gain) {
                    best_gain = gain;
                    best_i = i;
                    best_k = k;
                }
            }
        }

        if (best_i < 0) break;

        {
            int oldk = hap[best_i];
            int a;
            hap[best_i] = best_k;
            obj += best_gain;
            ++nmove;

            for (a = g->head[best_i]; a >= 0; a = g->arc[a].next) {
                int j = g->arc[a].to;
                double w = g->arc[a].w;
                score[(size_t)j * g->K + oldk] -= w;
                score[(size_t)j * g->K + best_k] += w;
            }
        }

        if (nmove >= max_moves) {
            if (VERBOSE > 1) {
                pthread_mutex_lock(&plock);
                fprintf(stderr,
                        "[W::phase_heu_local_search] move limit %d reached\n",
                        nmove);
                pthread_mutex_unlock(&plock);
            }
            break;
        }
    }

    (void) obj;

    if (nmove_out) *nmove_out = nmove;

    return phase_heu_objective(g, hap);
}

static int phase_heu_cmp_conf(const void *a, const void *b)
{
    const phase_heu_conf_t *x = (const phase_heu_conf_t *)a;
    const phase_heu_conf_t *y = (const phase_heu_conf_t *)b;

    if (x->margin < y->margin) return -1;
    if (x->margin > y->margin) return 1;
    return x->v - y->v;
}

static int phase_heu_perturb(const phase_heu_graph_t *g,
                             const unsigned char *fixed,
                             int *hap, double *score,
                             phase_heu_conf_t *conf,
                             uint64_t *rng)
{
    int i, k, q, nfree = 0;
    int npert, pool;

    phase_heu_build_scores(g, hap, score);

    for (i = 0; i < g->ns; ++i) {
        int oldk;
        double alt = -DBL_MAX;

        if (fixed[i]) continue;

        oldk = hap[i];
        for (k = 0; k < g->K; ++k)
            if (k != oldk &&
                score[(size_t)i * g->K + k] > alt)
                alt = score[(size_t)i * g->K + k];

        conf[nfree].v = i;
        conf[nfree].margin =
            score[(size_t)i * g->K + oldk] - alt;
        ++nfree;
    }

    if (nfree == 0 || g->K <= 1) return 0;

    qsort(conf, (size_t)nfree, sizeof(*conf), phase_heu_cmp_conf);

    npert = (int)ceil(PHASE_HEU_PERTURB_FRAC * nfree);
    if (npert < 2 && nfree >= 2) npert = 2;
    if (npert < 1) npert = 1;
    if (npert > PHASE_HEU_PERTURB_MAX) npert = PHASE_HEU_PERTURB_MAX;
    if (npert > nfree) npert = nfree;

    pool = 3 * npert;
    if (pool > nfree) pool = nfree;

    /* partial Fisher-Yates within lowest-confidence pool */
    for (q = 0; q < npert; ++q) {
        int r = q + phase_heu_rand_int(rng, pool - q);
        phase_heu_conf_t t = conf[q];
        conf[q] = conf[r];
        conf[r] = t;
    }

    for (q = 0; q < npert; ++q) {
        int v = conf[q].v;
        int oldk = hap[v], newk = oldk;

        if (phase_heu_rng01(rng) < 0.70) {
            double best = -DBL_MAX;
            for (k = 0; k < g->K; ++k) {
                double s;
                if (k == oldk) continue;
                s = score[(size_t)v * g->K + k];
                if (newk == oldk || s > best) {
                    best = s;
                    newk = k;
                }
            }
        } else {
            int r = phase_heu_rand_int(rng, g->K - 1);
            newk = r >= oldk ? r + 1 : r;
        }

        hap[v] = newk;
    }

    return npert;
}

static int phase_heu_read_complete_assignment(phase_data_t *data,
                                              long jid, int *hap)
{
    grp_info_t *grps = data->grps;
    hap_info_t *haps = data->haps;
    int *gseq = grps[jid].s;
    int ns = data->ns[jid], K = data->K;
    int i;

    for (i = 0; i < ns; ++i) {
        int h = haps[gseq[i]].hap;
        if (h < 1 || h > K) return 0;
        hap[i] = h - 1;
    }

    return 1;
}

static void phase_heu_store_assignment(phase_data_t *data,
                                       long jid, const int *hap)
{
    grp_info_t *grps = data->grps;
    hap_info_t *haps = data->haps;
    int *gseq = grps[jid].s;
    int ns = data->ns[jid];
    int i;

    for (i = 0; i < ns; ++i)
        haps[gseq[i]].hap = hap[i] + 1;
}


static int phase_heu_solve(phase_data_t *data, long jid,
                           double *best_obj_out)
{
    phase_heu_graph_t g;
    unsigned char *fixed = NULL;
    int *seed_hap = NULL;
    int *work = NULL, *cand = NULL, *rbest = NULL;
    int *gbest = NULL, *incoming = NULL;
    double *score = NULL;
    phase_heu_conf_t *conf = NULL;
    double global_obj = -DBL_MAX;
    int ns = data->ns[jid], K = data->K;
    int nseed, restart, have_incoming;
    uint64_t base_rng;

    if (ns <= 0 || K <= 0) return 0;
    if (!phase_heu_graph_build(&g, data, jid)) return 0;

    fixed = (unsigned char *)calloc((size_t)ns, sizeof(*fixed));
    seed_hap = (int *)malloc((size_t)ns * sizeof(*seed_hap));
    work = (int *)malloc((size_t)ns * sizeof(*work));
    cand = (int *)malloc((size_t)ns * sizeof(*cand));
    rbest = (int *)malloc((size_t)ns * sizeof(*rbest));
    gbest = (int *)malloc((size_t)ns * sizeof(*gbest));
    incoming = (int *)malloc((size_t)ns * sizeof(*incoming));
    score = (double *)malloc((size_t)ns * K * sizeof(*score));
    conf = (phase_heu_conf_t *)malloc((size_t)ns * sizeof(*conf));

    if (!fixed || !seed_hap || !work || !cand || !rbest ||
        !gbest || !incoming || !score || !conf)
        goto fail;

    have_incoming =
        phase_heu_read_complete_assignment(data, jid, incoming);

    nseed = phase_heu_make_seeds(data, jid, &g, fixed, seed_hap);
    if (nseed <= 0 && ns > 0) goto fail;

    if (VERBOSE > 1) {
        pthread_mutex_lock(&plock);
        fprintf(stderr,
                "[M::phase_heu_solve] jid=%ld ns=%d K=%d edges=%d "
                "seeds=%d starts=%d perturb=%d\n",
                jid, ns, K, g.nedge, nseed,
                PHASE_HEU_NSTART, PHASE_HEU_NPERTURB);
        pthread_mutex_unlock(&plock);
    }

    /* ILP/rolling incumbent, if any, is an extra start. */
    if (have_incoming) {
        int nm = 0;

        memcpy(work, incoming, (size_t)ns * sizeof(*work));

        if (phase_heu_align_to_seeds(ns, K, work, fixed, seed_hap)) {
            double obj =
                phase_heu_local_search(&g, fixed, work, score, &nm);

            if (obj > global_obj) {
                global_obj = obj;
                memcpy(gbest, work, (size_t)ns * sizeof(*gbest));
            }

            if (VERBOSE > 2) {
                pthread_mutex_lock(&plock);
                fprintf(stderr,
                        "[D::phase_heu_solve] incoming obj=%.9f moves=%d\n",
                        obj, nm);
                pthread_mutex_unlock(&plock);
            }
        }
    }

    base_rng = UINT64_C(0x6a09e667f3bcc909) ^
        ((uint64_t)(jid + 1) * UINT64_C(0x9e3779b97f4a7c15));

    for (restart = 0; restart < PHASE_HEU_NSTART; ++restart) {
        uint64_t rng = base_rng ^
            ((uint64_t)(restart + 1) * UINT64_C(0xbf58476d1ce4e5b9));
        double robj;
        int nm = 0, round;

        if (!phase_heu_greedy(&g, fixed, seed_hap,
                              restart != 0, &rng, work))
            goto fail;

        robj = phase_heu_local_search(&g, fixed, work, score, &nm);
        memcpy(rbest, work, (size_t)ns * sizeof(*rbest));

        if (robj > global_obj) {
            global_obj = robj;
            memcpy(gbest, rbest, (size_t)ns * sizeof(*gbest));
        }

        if (VERBOSE > 2) {
            pthread_mutex_lock(&plock);
            fprintf(stderr,
                    "[D::phase_heu_solve] start=%d obj=%.9f moves=%d\n",
                    restart + 1, robj, nm);
            pthread_mutex_unlock(&plock);
        }

        for (round = 0; round < PHASE_HEU_NPERTURB; ++round) {
            double obj;
            int np, lm = 0;

            memcpy(cand, rbest, (size_t)ns * sizeof(*cand));

            np = phase_heu_perturb(&g, fixed, cand, score, conf, &rng);
            if (np == 0) break;

            obj = phase_heu_local_search(&g, fixed, cand, score, &lm);

            if (obj > robj + PHASE_HEU_GAIN_EPS) {
                robj = obj;
                memcpy(rbest, cand, (size_t)ns * sizeof(*rbest));
            }

            if (obj > global_obj + PHASE_HEU_GAIN_EPS) {
                global_obj = obj;
                memcpy(gbest, cand, (size_t)ns * sizeof(*gbest));
            }

            if (VERBOSE > 3) {
                pthread_mutex_lock(&plock);
                fprintf(stderr,
                        "[D::phase_heu_solve] start=%d perturb=%d "
                        "npert=%d obj=%.9f moves=%d best=%.9f\n",
                        restart + 1, round + 1, np, obj, lm, global_obj);
                pthread_mutex_unlock(&plock);
            }
        }
    }

    if (global_obj == -DBL_MAX) goto fail;

    phase_heu_store_assignment(data, jid, gbest);
    if (best_obj_out) *best_obj_out = global_obj;

    if (VERBOSE > 1) {
        pthread_mutex_lock(&plock);
        fprintf(stderr,
                "[M::phase_heu_solve] jid=%ld final objective=%.9f\n",
                jid, global_obj);
        pthread_mutex_unlock(&plock);
    }

    free(fixed); free(seed_hap); free(work); free(cand);
    free(rbest); free(gbest); free(incoming); free(score); free(conf);
    phase_heu_graph_destroy(&g);
    return 1;

fail:
    free(fixed); free(seed_hap); free(work); free(cand);
    free(rbest); free(gbest); free(incoming); free(score); free(conf);
    phase_heu_graph_destroy(&g);
    return 0;
}


static void hybrid_phase_heu_thread(void *_data, long jid, int tid)
{
    phase_data_t *data = (phase_data_t *)_data;
    double obj = 0.0;
    int ok;

    (void)tid;

    if (PHASE_HEU_IS_SOLVED(data->status[jid]))
        return;

    ok = phase_heu_solve(data, jid, &obj);

    if (!ok) {
        if (VERBOSE > 0) {
            pthread_mutex_lock(&plock);
            fprintf(stderr,
                    "[E::%s] heuristic phasing failed for jid=%ld\n",
                    __func__, jid);
            pthread_mutex_unlock(&plock);
        }
        return;
    }
    data->status[jid] = PHASE_HEU_APPROX;

    if (VERBOSE > 0) {
        pthread_mutex_lock(&plock);
        fprintf(stderr,
                "[M::%s] heuristic completed for jid=%ld objective=%.9f\n",
                __func__, jid, obj);
        pthread_mutex_unlock(&plock);
    }
}

static void hungarian(int64 **rewards, int K, int8 *perms)
{
    int64 u[K + 1], v[K + 1], minv[K + 1];
    int p[K + 1], way[K + 1], used[K + 1];
    int64 cur, delta;
    int i, j, i0, j0, j1;

    for (i = 0; i <= K; i++) {
        u[i] = 0;
        v[i] = 0;
        p[i] = 0;
        way[i] = 0;
    }

    for (i = 1; i <= K; i++) {
        p[0] = i;
        j0 = 0;

        for (j = 0; j <= K; j++) {
            minv[j] = INT64_MAX;
            used[j] = 0;
        }

        do {
            used[j0] = 1;
            i0 = p[j0];
            delta = INT64_MAX;
            j1 = 0;

            for (j = 1; j <= K; j++) {
                if (!used[j]) {
                    cur = -rewards[i0 - 1][j - 1] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }

            for (j = 0; j <= K; j++) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }

            j0 = j1;
        } while (p[j0] != 0);

        do {
            j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }

    for (j = 1; j <= K; j++) {
        if (p[j] != 0) {
            perms[p[j] - 1] = (int8)(j - 1);
        }
    }
}

typedef struct {
    int64 r, l, s;
} ph_prior_t;

static void ph_prior_cpyfunc(void *x, void *y)
{
    *(ph_prior_t *) x = *(ph_prior_t *) y;
}

static int ph_prior_cmpfunc(void *x, void *y)
{
    ph_prior_t *a = (ph_prior_t *) x;
    ph_prior_t *b = (ph_prior_t *) y;

    // for selecting next to place
    // prioritise larger r, l and s
    if (a->r != b->r)
        return (a->r < b->r) - (a->r > b->r);
    if (a->l != b->l)
        return (a->l < b->l) - (a->l > b->l);
    if (a->s != b->s)
        return (a->s < b->s) - (a->s > b->s);
    return 0;
}


static void build_haplotype_partition_hybrid_core(ovl_t *ovls, int64 novl, hlk_t *hlks, int64 nhlk,
    int n_threads, hap_info_t *haps, sdict_t *dicts, int *reps, int ploidy)
{
    grp_info_t *grp, *grps;
    ord_i64_t *sord;
    ord_dbl_t *hord;
    phase_data_t _data, *data = &_data;
    hlk_t *hlk;
    ovl_t *ovl;
    hcft_t *cft, *cfts;
    uint64 *coff;
    uint32 a, b, nc;
    int64 ncft, hcft[ploidy+1], hseq[ploidy+1], hovl[ploidy+1], hhlk[ploidy+1];
    double p10, p90, idr;
    int *gseqs, *gseq, *mols, *mhls, *sidx;
    int i, j, k, c, n, p, g, ngrp, nseq;

    // sequence number
    nseq = (int) dicts->n;

    // reset haplotype assignments
    for (i = 0; i < nseq; i++)
        haps[i].hap = 0;

    // indexes of ovls and hlks
    MYCALLOC(hord, nhlk);
    if (hord == NULL)
        mem_alloc_error("haplotype partition offsets");

    // hic link normalisation
    for (i = 0; i < nhlk; i++) {
        hord[i].which = i;
        hord[i].event = hlks[i].v;
    }
    qsort(hord, nhlk, sizeof(ord_dbl_t), ord_dbl_acmpfunc);
    p10 = hord[(int) (nhlk * 0.1)].event;
    p90 = hord[(int) (nhlk * 0.9)].event;
    idr = p90 - p10;
    if (idr < 1e-6) idr = 1e-6;
    for (i = 0; i < nhlk; i++) {
        hord[i].event = 0.1 + 0.8 * (hord[i].event - p10) / idr;
        if (hord[i].event < 0) hord[i].event = 0.;
    }
    fprintf(stderr, "[M::%s] hic link normalisation: min=%.3f p10=%.3f p90=%.3f max=%.3f\n", __func__, 
        hord[0].event, p10, p90, hord[nhlk - 1].event);
    if (VERBOSE > 0)
        for (i = 0; i < 100; i++)
            fprintf(stderr, "[M::%s]   p%d = %.3f\n", __func__, i, hord[(int) (nhlk * i / 100.)].event);
    // sort back by index
    qsort(hord, nhlk, sizeof(ord_dbl_t), ord_dbl_wcmpfunc);

    // collect scaffold partition groups
    ngrp = 0;
    for (i = 0; i < nseq; i++)
        if (haps[i].grp > ngrp)
            ngrp = haps[i].grp;
    ngrp += 1;

    MYMALLOC(gseqs, nseq);
    MYCALLOC(grps, ngrp);
    MYMALLOC(sord, nseq);
    if (gseqs == NULL || grps == NULL || sord == NULL)
        mem_alloc_error("haplotype partition groups");

    for (i = 0; i < nseq; i++)
        grps[haps[i].grp].n++;

    gseq = gseqs;
    for (i = 0; i < ngrp; i++) {
        grp = grps + i;
        grp->g = i;
        grp->s = gseq;
        gseq += grp->n;
        grp->n = 0;
    }
    
    for (i = 0; i < nseq; i++) {
        g = haps[i].grp;
        grp = grps + g;
        grp->s[grp->n++] = i;
        grp->t += dicts->s[i].len;
    }

    // within each group sort sequences by start position ascending
    for (i = 0; i < ngrp; i++) {
        grp = grps + i;
        for (j = 0; j < grp->n; j++)
            sord[j] = (ord_i64_t) {grp->s[j], haps[grp->s[j]].bpos};
        qsort(sord, grp->n, sizeof(ord_i64_t), ord_i64_acmpfunc);
        for (j = 0; j < grp->n; j++)
            grp->s[j] = sord[j].which;
    }

    // sort groups by size
    // group 0 is always last (unassigned sequences)
    qsort(grps, ngrp, sizeof(grp_info_t), grp_info_cmpfunc);

    MYMALLOC(sidx, nseq);
    MYMALLOC(mols, novl);
    MYMALLOC(mhls, nhlk);
    if (sidx == NULL || mols == NULL || mhls == NULL)
        mem_alloc_error("haplotype partition arrays");

    // initialise phase_data_t structure for multi-threading
    phase_data_init(data, ngrp);
    data->grps = grps;
    data->haps = haps;
    data->ovls = ovls;
    data->hlks = hlks;
    data->hord = hord;
    data->K = ploidy;
    data->sidx = sidx;

    int ns, no, nh, u, v;
    int *ols = mols, *hls = mhls;
    for (c = 0; c < ngrp; c++) {
        // basic group information
        grp = grps + c;
        gseq = grp->s; // sequences in the group
        g = grp->g; // real group id
        ns = grp->n; // number of sequences in the group
        for (i = 0; i < ns; i++)
            sidx[gseq[i]] = i;
        no = 0;
        for (i = 0, no = 0; i < novl; i++) {
            if (ovls[i].del) continue;
            u = ovls[i].aread;
            v = ovls[i].bread;
            if (u < v && haps[u].grp == g && haps[v].grp == g)
                ols[no++] = i;
        }
        nh = 0;
        for (i = 0; i < nhlk; i++) {
            u = hlks[i].a;
            v = hlks[i].b;
            if (u < v && haps[u].grp == g && haps[v].grp == g)
                hls[nh++] = i;
        }
        
        data->ns[c] = ns;
        data->no[c] = no;
        data->nh[c] = nh;
        data->ols[c] = ols;
        data->hls[c] = hls;
        data->status[c] = PHASE_UNRESOLVED;
        ols += no;
        hls += nh;
    }

    // the actual haplotype phasing
#ifdef DEBUG_HAPLOTYPE_PHASE
    for (c = 0; c < ngrp; c++) {
        hybrid_phase_ilp_thread(data, c, 0);
        hybrid_phase_heu_thread(data, c, 0);
    }
#else
    kt_for(n_threads, hybrid_phase_ilp_thread, data, ngrp);
    kt_for(n_threads, hybrid_phase_heu_thread, data, ngrp);
#endif

    {
        int cnts[ploidy+1];
        int64 lens[ploidy+1];
        for (i = 0; i <= ploidy; i++) {
            cnts[i] = 0;
            lens[i] = 0;
        }
        for (i = 0; i < nseq; i++) {
            k = haps[i].hap;
            cnts[k]++;
            lens[k] += haps[i].len;
        }
        fprintf(stderr, "[M::%s] Haplotype partition summary RAW:\n", __func__);
        for (k = 0; k <= ploidy; k++)
            fprintf(stderr, "[M::%s]    H%d:  %6d seqs %12lld bp\n", __func__, k, cnts[k], lens[k]);
    }
    {
        int64 intra_ovls = 0, inter_ovls = 0;
        for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
            int a, b;
            if (ovl->del)
                continue;
            a = ovl->aread;
            b = ovl->bread;
            // only count overlaps between sequences in the same scaffold group
            if (haps[a].grp != haps[b].grp || haps[a].grp == 0)
                continue;
            if (haps[a].hap > 0 && haps[a].hap == haps[b].hap)
                intra_ovls += (ovl->alen + ovl->blen) / 2;
            else
                inter_ovls += (ovl->alen + ovl->blen) / 2;
        }
        fprintf(stderr, "[M::%s] within-group intra_ovls: %16lld inter_ovls: %16lld\n", __func__, intra_ovls, inter_ovls);
    }

    // process unassigned sequences using greedy colouring
    // number of conflicts
    ncft = nhlk + novl;
    // collect conflict table entries
    MYMALLOC(coff, nseq);
    MYCALLOC(cfts, ncft);
    if (cfts == NULL || coff == NULL)
        mem_alloc_error("haplotype partition offsets");
    cft = cfts;
    for (i = 0, hlk = hlks; i < nhlk; i++, hlk++)
        *(cft++) = (hcft_t) { hlk->a, hlk->b, (int64) (hlk->l * HIC_NORM_WINDOW * hord[i].event) << 8 | CFT_HIC };
    for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
        if (ovl->del) continue;
        *(cft++) = (hcft_t) { (int)ovl->aread, (int)ovl->bread, (int64) (ovl->alen + ovl->blen) << 8 |  CFT_OVL };
    }
    // sort by (a,b) and merge duplicate directed pairs (sum net weight)
    qsort(cfts, ncft, sizeof(hcft_t), hcft_ab_cmpfunc);
    // build offsets
    if (ncft > 0) {
        MYBZERO(coff, nseq);
        a = cfts[0].a;
        for (i = 1, j = 0; i < ncft; i++)
            if (cfts[i].a != a) { coff[a] = (uint64) j << 32 | (i - j); a = cfts[i].a; j = i; }
        coff[a] = (uint64) j << 32 | (ncft - j);
    }
    // now process unassigned sequences and repeats
    n = 0;
    for (i = 0; i < nseq; i++) {
        if (haps[i].hap > 0) continue;
        sord[n++] = (ord_i64_t) {i, haps[i].len};
    }
    qsort(sord, n, sizeof(ord_i64_t), ord_i64_dcmpfunc);
    for (i = 0; i < n; i++) {
        a = sord[i].which;
        if (haps[a].hap > 0) continue;
        cft = cfts + (coff[a] >> 32);
        nc  = (uint32) coff[a];
        MYBZERO(hcft, ploidy+1);
        MYBZERO(hhlk, ploidy+1);
        MYBZERO(hovl, ploidy+1);
        for (k = 0; k < nc; k++) {
            b = cft[k].b;
            if ((p = haps[b].hap) == 0)
                continue;
            if ((cft[k].v & 0xff) == CFT_OVL)
                hovl[p] += cft[k].v;
            else if ((cft[k].v & 0xff) == CFT_HIC && 
                hhlk[p] < cft[k].v)
                hhlk[p] = cft[k].v;
        }
        for (k = 1; k <= ploidy; k++)
            hcft[k] = hovl[k] - hhlk[k];

        p = 1;
        for (k = 2; k <= ploidy; k++)
            if (ccompare(k, p, hcft, hovl, hhlk, hseq) < 0)
                p = k;
        haps[a].hap = p;
        hseq[p] += haps[a].len;
    }

    {
        int cnts[ploidy+1];
        int64 lens[ploidy+1];
        for (i = 0; i <= ploidy; i++) {
            cnts[i] = 0;
            lens[i] = 0;
        }
        for (i = 0; i < nseq; i++) {
            k = haps[i].hap;
            cnts[k]++;
            lens[k] += haps[i].len;
        }
        fprintf(stderr, "[M::%s] Haplotype partition summary:\n", __func__);
        for (k = 0; k <= ploidy; k++)
            fprintf(stderr, "[M::%s]    H%d:  %6d seqs %12lld bp\n", __func__, k, cnts[k], lens[k]);
    }

    /*
     * Group-level phasing.
     *
     * Each scaffold group has already been internally partitioned into K
     * local haplotypes.  The only remaining decision is therefore a
     * permutation of the K local labels for each group.
     *
     * Build a signed group-level linkage graph:
     *
     *      Hi-C  : positive reward for the same final haplotype
     *      OVL   : negative reward for the same final haplotype
     *
     * The weights follow the sequence-level MILP objective (without the
     * common 1e6 scaling, which does not affect the optimum):
     *
     *      w_H = l * HIC_NORM_WINDOW * event * PHASE_BETA
     *      w_O = (alen + blen) * log10(alen + blen) * PHASE_ALPHA
     *
     * Each physical cross-group edge is inserted in both directions so
     * every group can evaluate its own permutation directly.
     */
    {
        hcft_t *gcfts = NULL, *gcft;
        uint64 *gcoff = NULL;
        int64 gncft = 0, gcap;
        int64 same_ovl_before = 0, same_ovl_after = 0, cross_ovl = 0;
        int K = ploidy, K2 = K * K;

        /*
         * Diagnostic before group-level label synchronization.
         */
        for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
            int64 x;
            a = ovl->aread;
            b = ovl->bread;

            if (ovl->del ||
                haps[a].grp <= 0 || haps[b].grp <= 0 ||
                haps[a].hap <= 0 || haps[b].hap <= 0 ||
                haps[a].grp == haps[b].grp)
                continue;

            x = ((int64)ovl->alen + (int64)ovl->blen) / 2;
            cross_ovl += x;
            if (haps[a].hap == haps[b].hap)
                same_ovl_before += x;
        }

        /*
         * At most two directed entries are generated from each original
         * Hi-C/OVL edge.
         */
        gcap = 2 * (nhlk + novl);
        if (gcap < 1)
            gcap = 1;

        MYMALLOC(gcfts, gcap);
        MYCALLOC(gcoff, ngrp);
        if (gcfts == NULL || gcoff == NULL)
            mem_alloc_error("group-level haplotype linkage arrays");

        /*
         * Positive Hi-C group-level links.
         */
        for (i = 0, hlk = hlks; i < nhlk; i++, hlk++) {
            uint32 ga, gb;
            int64 iw;
            double dw;

            a = hlk->a;
            b = hlk->b;

            if (haps[a].grp <= 0 || haps[b].grp <= 0 ||
                haps[a].hap <= 0 || haps[b].hap <= 0 ||
                haps[a].grp == haps[b].grp)
                continue;

            ga = ((haps[a].grp - 1) << 8) | (haps[a].hap - 1);
            gb = ((haps[b].grp - 1) << 8) | (haps[b].hap - 1);

            dw = hlk->l * HIC_NORM_WINDOW * hord[i].event * PHASE_BETA;
            iw = (int64) llround(dw);
            if (iw == 0 && dw > 0.0)
                iw = 1;

            gcfts[gncft++] = (hcft_t) {ga, gb, +iw};
            gcfts[gncft++] = (hcft_t) {gb, ga, +iw};
        }

        /*
         * Negative overlap group-level links.
         */
        for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
            uint32 ga, gb;
            int64 iw;
            double dw, x;

            if (ovl->del)
                continue;

            a = ovl->aread;
            b = ovl->bread;

            if (haps[a].grp <= 0 || haps[b].grp <= 0 ||
                haps[a].hap <= 0 || haps[b].hap <= 0 ||
                haps[a].grp == haps[b].grp)
                continue;

            ga = ((haps[a].grp - 1) << 8) | (haps[a].hap - 1);
            gb = ((haps[b].grp - 1) << 8) | (haps[b].hap - 1);

            x = (double)ovl->alen + (double)ovl->blen;
            dw = x * log10(x) * PHASE_ALPHA;
            iw = (int64) llround(dw);
            if (iw == 0 && dw > 0.0)
                iw = 1;

            gcfts[gncft++] = (hcft_t) {ga, gb, -iw};
            gcfts[gncft++] = (hcft_t) {gb, ga, -iw};
        }

        if (VERBOSE > 0) {
            fprintf(stderr,
                    "[M::%s] cross-group overlaps before synchronization: "
                    "same-hap=%lld / total=%lld (%.2f%%)\n",
                    __func__,
                    (long long)same_ovl_before,
                    (long long)cross_ovl,
                    cross_ovl > 0 ?
                    100.0 * (double)same_ovl_before / (double)cross_ovl :
                    0.0);
        }

        if (gncft > 0) {
            int8 h, ha, hb, hc, *done, *perms, *match;
            int64 **links;
            int64 m;

            /*
             * Merge duplicate directed (group,hap) pairs.  Positive Hi-C and
             * negative OVL support therefore collapse into one net weight.
             */
            qsort(gcfts, gncft, sizeof(hcft_t), hcft_ab_cmpfunc);

            m = 0;
            for (i = 0; i < gncft; i++) {
                if (m > 0 &&
                    gcfts[m-1].a == gcfts[i].a &&
                    gcfts[m-1].b == gcfts[i].b) {
                    gcfts[m-1].v += gcfts[i].v;
                } else {
                    gcfts[m++] = gcfts[i];
                }
            }

            /*
             * Remove pairs whose positive and negative evidence cancels
             * exactly.
             */
            gncft = 0;
            for (i = 0; i < m; i++)
                if (gcfts[i].v != 0)
                    gcfts[gncft++] = gcfts[i];

            if (gncft > 0) {
                /*
                 * Offset table by 0-based real group id (haps[].grp - 1).
                 */
                MYBZERO(gcoff, ngrp);
                a = gcfts[0].a >> 8;
                for (i = 1, j = 0; i < gncft; i++) {
                    uint32 ga = gcfts[i].a >> 8;
                    if (ga != a) {
                        gcoff[a] = (uint64)j << 32 | (uint32)(i - j);
                        a = ga;
                        j = i;
                    }
                }
                gcoff[a] = (uint64)j << 32 | (uint32)(gncft - j);
            }

            MYCALLOC(done, (size_t)ngrp);
            MYCALLOC(perms, (size_t)ngrp * K);
            MYMALLOC(match, (size_t)K);
            MYMALLOC(links, (size_t)K);
            if (done == NULL || perms == NULL ||
                links == NULL || match == NULL)
                mem_alloc_error("group-level haplotype permutation arrays");

            MYMALLOC(links[0], (size_t)K2);
            if (links[0] == NULL)
                mem_alloc_error("group-level haplotype linkage matrix");
            for (i = 1; i < K; i++)
                links[i] = links[i-1] + K;

            /*
             * Initialise all group permutations to identity.  The first
             * group selected in each connected component therefore anchors
             * the arbitrary global haplotype labels.
             */
            for (c = 0; c < ngrp; c++) {
                g = grps[c].g - 1;
                if (g < 0)
                    continue;
                for (k = 0; k < K; k++)
                    perms[g * K + k] = k;
            }

            /*
             * Initial placement.
             *
             * Groups are selected by evidence strength.  Importantly, use
             * sum(abs(weight)), not the signed sum, because strong Hi-C and
             * strong OVL evidence should not cancel when deciding which
             * group is most informative to place next.
             */
            {
                ph_prior_t _pp, *pp = &_pp;
                pq_t *pq;
                pq_elem_t pe;
                int64 strength;
                int8 *flag;
                int *list;

                pq = pq_init(ngrp, sizeof(ph_prior_t),
                             ph_prior_cpyfunc, ph_prior_cmpfunc);
                MYMALLOC(list, ngrp);
                MYMALLOC(flag, ngrp);
                if (list == NULL || flag == NULL || pq == NULL)
                    mem_alloc_error("group-level phasing priority arrays");

                for (c = 0; c < ngrp; c++) {
                    g = grps[c].g - 1;
                    if (g < 0)
                        continue;

                    gcft = gcfts + (gcoff[g] >> 32);
                    nc = (uint32)gcoff[g];
                    strength = 0;
                    for (k = 0; k < nc; k++)
                        strength += llabs(gcft[k].v);

                    *pp = (ph_prior_t) {0, strength, grps[c].t};
                    pq_insert(pq, g, pp);
                }

                while (pq_pop1(pq, &pe)) {
                    g = pe.id;
                    done[g] = 1;

                    /*
                     * Recompute every not-yet-placed neighbouring group
                     * against all groups already placed.
                     */
                    MYBZERO(flag, ngrp);
                    gcft = gcfts + (gcoff[g] >> 32);
                    nc = (uint32)gcoff[g];
                    n = 0;

                    for (i = 0; i < nc; i++) {
                        b = gcft[i].b >> 8;
                        if (done[b] || flag[b])
                            continue;
                        flag[b] = 1;
                        list[n++] = b;
                    }

                    for (i = 0; i < n; i++) {
                        int64 related = 0;

                        a = list[i];
                        gcft = gcfts + (gcoff[a] >> 32);
                        nc = (uint32)gcoff[a];

                        MYBZERO(links[0], K2);

                        for (k = 0; k < nc; k++) {
                            b = gcft[k].b >> 8;
                            if (!done[b])
                                continue;

                            ha = gcft[k].a & 0xff;
                            hb = gcft[k].b & 0xff;
                            hc = perms[b * K + hb];

                            links[ha][hc] += gcft[k].v;
                            related += llabs(gcft[k].v);
                        }

                        hungarian(links, K, match);
                        for (k = 0; k < K; k++)
                            perms[a * K + k] = match[k];

                        /*
                         * 'r' measures how strongly this group is connected
                         * to already placed groups, irrespective of sign.
                         */
                        ph_prior_cpyfunc(pp, pq_prior(pq, a, 1));
                        pp->r = related;
                        pq_insert(pq, a, pp);
                    }
                }

                free(list);
                free(flag);
                pq_destroy(pq);
            }

            /*
             * Coordinate-ascent refinement.
             *
             * With all other group permutations fixed, Hungarian gives the
             * exact best permutation for one group.  Accept only a strictly
             * positive gain, so the signed group-level objective increases
             * monotonically until no single-group permutation can improve it.
             */
            {
                ph_prior_t _pp, *pp = &_pp;
                pq_t *pq;
                pq_elem_t pe;
                int64 gain;
                int8 *permx, *flag;
                int *list;

                pq = pq_init(ngrp, sizeof(ph_prior_t),
                             ph_prior_cpyfunc, ph_prior_cmpfunc);
                MYMALLOC(permx, ngrp * K);
                MYMALLOC(list, ngrp);
                MYMALLOC(flag, ngrp);
                if (pq == NULL || permx == NULL ||
                    list == NULL || flag == NULL)
                    mem_alloc_error("group-level refinement arrays");

                /*
                 * Initial gain for every group.
                 */
                for (i = 0; i < ngrp; i++) {
                    g = grps[i].g - 1;
                    if (g < 0)
                        continue;

                    gcft = gcfts + (gcoff[g] >> 32);
                    nc = (uint32)gcoff[g];

                    MYBZERO(links[0], K2);

                    for (k = 0; k < nc; k++) {
                        b = gcft[k].b >> 8;
                        ha = gcft[k].a & 0xff;
                        hb = gcft[k].b & 0xff;
                        hc = perms[b * K + hb];
                        links[ha][hc] += gcft[k].v;
                    }

                    hungarian(links, K, match);

                    /*
                     * Compare against the CURRENT permutation, not identity.
                     */
                    gain = 0;
                    for (k = 0; k < K; k++)
                        gain += links[k][match[k]]
                              - links[k][perms[g * K + k]];

                    if (gain > 0) {
                        for (k = 0; k < K; k++)
                            permx[g * K + k] = match[k];
                        *pp = (ph_prior_t) {gain, 0, 0};
                        pq_insert(pq, g, pp);
                    }
                }

                while (pq_pop1(pq, &pe)) {
                    g = pe.id;

                    /*
                     * Apply the best currently known improving permutation.
                     */
                    for (k = 0; k < K; k++)
                        perms[g * K + k] = permx[g * K + k];

                    /*
                     * Recompute gains for all neighbouring groups.  Do NOT
                     * use done[] here: all groups are already placed.
                     */
                    MYBZERO(flag, ngrp);
                    gcft = gcfts + (gcoff[g] >> 32);
                    nc = (uint32)gcoff[g];
                    n = 0;

                    for (i = 0; i < nc; i++) {
                        b = gcft[i].b >> 8;
                        if (b == (uint32)g || flag[b])
                            continue;
                        flag[b] = 1;
                        list[n++] = b;
                    }

                    for (i = 0; i < n; i++) {
                        a = list[i];
                        gcft = gcfts + (gcoff[a] >> 32);
                        nc = (uint32)gcoff[a];

                        MYBZERO(links[0], K2);

                        for (k = 0; k < nc; k++) {
                            b = gcft[k].b >> 8;
                            ha = gcft[k].a & 0xff;
                            hb = gcft[k].b & 0xff;
                            hc = perms[b * K + hb];
                            links[ha][hc] += gcft[k].v;
                        }

                        hungarian(links, K, match);

                        gain = 0;
                        for (k = 0; k < K; k++)
                            gain += links[k][match[k]]
                                  - links[k][perms[a * K + k]];

                        if (gain > 0) {
                            for (k = 0; k < K; k++)
                                permx[a * K + k] = match[k];
                            *pp = (ph_prior_t) {gain, 0, 0};
                            pq_insert(pq, a, pp);
                        } else {
                            pq_remove(pq, a);
                        }
                    }
                }

                free(permx);
                free(list);
                free(flag);
                pq_destroy(pq);
            }

            /*
             * Apply the final group permutations to sequence haplotypes.
             */
            for (i = 0; i < nseq; i++) {
                g = haps[i].grp - 1;
                h = haps[i].hap - 1;
                if (g < 0 || h < 0)
                    continue;
                haps[i].hap = perms[g * K + h] + 1;
            }

            free(links[0]);
            free(links);
            free(done);
            free(match);
            free(perms);
        }

        /*
         * Cross-group overlap diagnostic after synchronization.
         */
        for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
            int64 x;

            if (ovl->del)
                continue;

            a = ovl->aread;
            b = ovl->bread;

            if (haps[a].grp <= 0 || haps[b].grp <= 0 ||
                haps[a].hap <= 0 || haps[b].hap <= 0 ||
                haps[a].grp == haps[b].grp)
                continue;

            x = ((int64)ovl->alen + (int64)ovl->blen) / 2;
            if (haps[a].hap == haps[b].hap)
                same_ovl_after += x;
        }

        if (VERBOSE > 0) {
            fprintf(stderr,
                    "[M::%s] cross-group overlaps after synchronization:  "
                    "same-hap=%lld / total=%lld (%.2f%%); delta=%lld\n",
                    __func__,
                    (long long)same_ovl_after,
                    (long long)cross_ovl,
                    cross_ovl > 0 ?
                    100.0 * (double)same_ovl_after / (double)cross_ovl :
                    0.0,
                    (long long)(same_ovl_after - same_ovl_before));
        }

        free(gcfts);
        free(gcoff);
    }

    {
        int cnts[ploidy+1];
        int64 lens[ploidy+1];
        for (i = 0; i <= ploidy; i++) {
            cnts[i] = 0;
            lens[i] = 0;
        }
        for (i = 0; i < nseq; i++) {
            k = haps[i].hap;
            cnts[k]++;
            lens[k] += haps[i].len;
        }
        fprintf(stderr, "[M::%s] Haplotype partition summary:\n", __func__);
        for (k = 0; k <= ploidy; k++)
            fprintf(stderr, "[M::%s]    H%d:  %6d seqs %12lld bp\n", __func__, k, cnts[k], lens[k]);
    }

    free(gseqs);
    free(grps);
    free(sord);
    free(hord);
    free(sidx);
    free(mols);
    free(mhls);
    free(cfts);
    free(coff);
    phase_data_destroy(data);
}

static void build_haplotype_partition_hybrid_core1(ovl_t *ovls, int64 novl, hlk_t *hlks, int64 nhlk,
    int n_threads, hap_info_t *haps, sdict_t *dicts, int *reps, int ploidy)
{
    grp_info_t *grp, *grps;
    ord_i64_t *sord;
    ord_dbl_t *hord;
    phase_data_t _data, *data = &_data;
    hlk_t *hlk;
    ovl_t *ovl;
    hcft_t *cft, *cfts;
    uint64 *coff;
    uint32 a, b, nc;
    int64 ncft, hcft[ploidy+1], hseq[ploidy+1], hovl[ploidy+1], hhlk[ploidy+1];
    double p10, p90, idr;
    int *gseqs, *gseq, *mols, *mhls, *sidx;
    int i, j, k, c, n, p, g, ngrp, nseq;

    // sequence number
    nseq = (int) dicts->n;

    // reset haplotype assignments
    for (i = 0; i < nseq; i++)
        haps[i].hap = 0;

    // indexes of ovls and hlks
    MYCALLOC(hord, nhlk);
    if (hord == NULL)
        mem_alloc_error("haplotype partition offsets");

    // hic link normalisation
    for (i = 0; i < nhlk; i++) {
        hord[i].which = i;
        hord[i].event = hlks[i].v;
    }
    qsort(hord, nhlk, sizeof(ord_dbl_t), ord_dbl_acmpfunc);
    p10 = hord[(int) (nhlk * 0.1)].event;
    p90 = hord[(int) (nhlk * 0.9)].event;
    idr = p90 - p10;
    if (idr < 1e-6) idr = 1e-6;
    for (i = 0; i < nhlk; i++) {
        hord[i].event = 0.1 + 0.8 * (hord[i].event - p10) / idr;
        if (hord[i].event < 0) hord[i].event = 0.;
    }
    fprintf(stderr, "[M::%s] hic link normalisation: min=%.3f p10=%.3f p90=%.3f max=%.3f\n", __func__, 
        hord[0].event, p10, p90, hord[nhlk - 1].event);
    if (VERBOSE > 0)
        for (i = 0; i < 100; i++)
            fprintf(stderr, "[M::%s]   p%d = %.3f\n", __func__, i, hord[(int) (nhlk * i / 100.)].event);
    // sort back by index
    qsort(hord, nhlk, sizeof(ord_dbl_t), ord_dbl_wcmpfunc);

    // collect scaffold partition groups
    ngrp = 0;
    for (i = 0; i < nseq; i++)
        if (haps[i].grp > ngrp)
            ngrp = haps[i].grp;
    ngrp += 1;

    MYMALLOC(gseqs, nseq);
    MYCALLOC(grps, ngrp);
    MYMALLOC(sord, nseq);
    if (gseqs == NULL || grps == NULL || sord == NULL)
        mem_alloc_error("haplotype partition groups");

    for (i = 0; i < nseq; i++)
        grps[haps[i].grp].n++;

    gseq = gseqs;
    for (i = 0; i < ngrp; i++) {
        grp = grps + i;
        grp->g = i;
        grp->s = gseq;
        gseq += grp->n;
        grp->n = 0;
    }
    
    for (i = 0; i < nseq; i++) {
        g = haps[i].grp;
        grp = grps + g;
        grp->s[grp->n++] = i;
        grp->t += dicts->s[i].len;
    }

    // within each group sort sequences by start position ascending
    for (i = 0; i < ngrp; i++) {
        grp = grps + i;
        for (j = 0; j < grp->n; j++)
            sord[j] = (ord_i64_t) {grp->s[j], haps[grp->s[j]].bpos};
        qsort(sord, grp->n, sizeof(ord_i64_t), ord_i64_acmpfunc);
        for (j = 0; j < grp->n; j++)
            grp->s[j] = sord[j].which;
    }

    // sort groups by size
    // group 0 is always last (unassigned sequences)
    qsort(grps, ngrp, sizeof(grp_info_t), grp_info_cmpfunc);

    MYMALLOC(sidx, nseq);
    MYMALLOC(mols, novl);
    MYMALLOC(mhls, nhlk);
    if (sidx == NULL || mols == NULL || mhls == NULL)
        mem_alloc_error("haplotype partition arrays");

    // initialise phase_data_t structure for multi-threading
    phase_data_init(data, ngrp);
    data->grps = grps;
    data->haps = haps;
    data->ovls = ovls;
    data->hlks = hlks;
    data->hord = hord;
    data->K = ploidy;
    data->sidx = sidx;

    int ns, no, nh, u, v;
    int *ols = mols, *hls = mhls;
    for (c = 0; c < ngrp; c++) {
        // basic group information
        grp = grps + c;
        gseq = grp->s; // sequences in the group
        g = grp->g; // real group id
        ns = grp->n; // number of sequences in the group
        for (i = 0; i < ns; i++)
            sidx[gseq[i]] = i;
        no = 0;
        for (i = 0, no = 0; i < novl; i++) {
            if (ovls[i].del) continue;
            u = ovls[i].aread;
            v = ovls[i].bread;
            if (u < v && haps[u].grp == g && haps[v].grp == g)
                ols[no++] = i;
        }
        nh = 0;
        for (i = 0; i < nhlk; i++) {
            u = hlks[i].a;
            v = hlks[i].b;
            if (u < v && haps[u].grp == g && haps[v].grp == g)
                hls[nh++] = i;
        }
        
        data->ns[c] = ns;
        data->no[c] = no;
        data->nh[c] = nh;
        data->ols[c] = ols;
        data->hls[c] = hls;
        data->status[c] = PHASE_UNRESOLVED;
        ols += no;
        hls += nh;
    }

    // the actual haplotype phasing
#ifdef DEBUG_HAPLOTYPE_PHASE
    for (c = 0; c < ngrp; c++) {
        hybrid_phase_ilp_thread(data, c, 0);
        hybrid_phase_heu_thread(data, c, 0);
    }
#else
    kt_for(n_threads, hybrid_phase_ilp_thread, data, ngrp);
    kt_for(n_threads, hybrid_phase_heu_thread, data, ngrp);
#endif

    {
        int cnts[ploidy+1];
        int64 lens[ploidy+1];
        for (i = 0; i <= ploidy; i++) {
            cnts[i] = 0;
            lens[i] = 0;
        }
        for (i = 0; i < nseq; i++) {
            k = haps[i].hap;
            cnts[k]++;
            lens[k] += haps[i].len;
        }
        fprintf(stderr, "[M::%s] Haplotype partition summary RAW:\n", __func__);
        for (k = 0; k <= ploidy; k++)
            fprintf(stderr, "[M::%s]    H%d:  %6d seqs %12lld bp\n", __func__, k, cnts[k], lens[k]);
    }
    {
        int64 intra_ovls = 0, inter_ovls = 0;
        for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
            int a, b;
            if (ovl->del)
                continue;
            a = ovl->aread;
            b = ovl->bread;
            // only count overlaps between sequences in the same scaffold group
            if (haps[a].grp != haps[b].grp || haps[a].grp == 0)
                continue;
            if (haps[a].hap > 0 && haps[a].hap == haps[b].hap)
                intra_ovls += (ovl->alen + ovl->blen) / 2;
            else
                inter_ovls += (ovl->alen + ovl->blen) / 2;
        }
        fprintf(stderr, "[M::%s] within-group intra_ovls: %16lld inter_ovls: %16lld\n", __func__, intra_ovls, inter_ovls);
    }

    // process unassigned sequences using greedy colouring
    // number of conflicts
    ncft = nhlk + novl;
    // collect conflict table entries
    MYMALLOC(coff, nseq);
    MYCALLOC(cfts, ncft);
    if (cfts == NULL || coff == NULL)
        mem_alloc_error("haplotype partition offsets");
    cft = cfts;
    for (i = 0, hlk = hlks; i < nhlk; i++, hlk++)
        *(cft++) = (hcft_t) { hlk->a, hlk->b, (int64) (hlk->l * HIC_NORM_WINDOW * hord[i].event) << 8 | CFT_HIC };
    for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
        if (ovl->del) continue;
        *(cft++) = (hcft_t) { (int)ovl->aread, (int)ovl->bread, (int64) (ovl->alen + ovl->blen) << 8 |  CFT_OVL };
    }
    // sort by (a,b) and merge duplicate directed pairs (sum net weight)
    qsort(cfts, ncft, sizeof(hcft_t), hcft_ab_cmpfunc);
    // build offsets
    if (ncft > 0) {
        MYBZERO(coff, nseq);
        a = cfts[0].a;
        for (i = 1, j = 0; i < ncft; i++)
            if (cfts[i].a != a) { coff[a] = (uint64) j << 32 | (i - j); a = cfts[i].a; j = i; }
        coff[a] = (uint64) j << 32 | (ncft - j);
    }
    // now process unassigned sequences and repeats
    n = 0;
    for (i = 0; i < nseq; i++) {
        if (haps[i].hap > 0) continue;
        sord[n++] = (ord_i64_t) {i, haps[i].len};
    }
    qsort(sord, n, sizeof(ord_i64_t), ord_i64_dcmpfunc);
    for (i = 0; i < n; i++) {
        a = sord[i].which;
        if (haps[a].hap > 0) continue;
        cft = cfts + (coff[a] >> 32);
        nc  = (uint32) coff[a];
        MYBZERO(hcft, ploidy+1);
        MYBZERO(hhlk, ploidy+1);
        MYBZERO(hovl, ploidy+1);
        for (k = 0; k < nc; k++) {
            b = cft[k].b;
            if ((p = haps[b].hap) == 0)
                continue;
            if ((cft[k].v & 0xff) == CFT_OVL)
                hovl[p] += cft[k].v;
            else if ((cft[k].v & 0xff) == CFT_HIC && 
                hhlk[p] < cft[k].v)
                hhlk[p] = cft[k].v;
        }
        for (k = 1; k <= ploidy; k++)
            hcft[k] = hovl[k] - hhlk[k];

        p = 1;
        for (k = 2; k <= ploidy; k++)
            if (ccompare(k, p, hcft, hovl, hhlk, hseq) < 0)
                p = k;
        haps[a].hap = p;
        hseq[p] += haps[a].len;
    }

    {
        int cnts[ploidy+1];
        int64 lens[ploidy+1];
        for (i = 0; i <= ploidy; i++) {
            cnts[i] = 0;
            lens[i] = 0;
        }
        for (i = 0; i < nseq; i++) {
            k = haps[i].hap;
            cnts[k]++;
            lens[k] += haps[i].len;
        }
        fprintf(stderr, "[M::%s] Haplotype partition summary:\n", __func__);
        for (k = 0; k <= ploidy; k++)
            fprintf(stderr, "[M::%s]    H%d:  %6d seqs %12lld bp\n", __func__, k, cnts[k], lens[k]);
    }

    // group-level phasing
    // convert cfts to group-level conflict table
    // only consider CFT_HIC links for group-level phasing
    // (group<<8 | hap)
    n = 0;
    for (i = 0; i < ncft; i++) {
#define USE_HIC_ONLY
#ifdef USE_HIC_ONLY
        if ((cfts[i].v & 0xff) != CFT_HIC)
            continue;
#endif
        a = cfts[i].a;
        b = cfts[i].b;
        if (haps[a].grp <= 0 || 
            haps[b].grp <= 0 ||
            haps[a].hap <= 0 ||
            haps[b].hap <= 0 ||
            haps[a].grp == haps[b].grp)
            continue;
        cfts[n].a = ((haps[a].grp-1) << 8) | (haps[a].hap-1);
        cfts[n].b = ((haps[b].grp-1) << 8) | (haps[b].hap-1);
        cfts[n].v = cfts[i].v >> 8;
#ifndef USE_HIC_ONLY
        if ((cfts[i].v & 0xff) == CFT_OVL)
            cfts[n].v = -cfts[i].v;
#endif
        n++;
    }
    if (n > 0) {
        int8 h, ha, hb, hc, *done, *perms, *match;
        int64 **links;
        int K = ploidy, K2 = K * K;

        MYCALLOC(done, (size_t) ngrp);
        MYCALLOC(perms, (size_t) ngrp * K);
        MYMALLOC(match, (size_t) K);
        MYMALLOC(links, (size_t) K);
        if (done == NULL || perms == NULL || links == NULL || match == NULL)
            mem_alloc_error("haplotype partition arrays");
        MYMALLOC(links[0], (size_t) K2);
        if (links[0] == NULL)
            mem_alloc_error("haplotype partition arrays");
        for (i = 1; i < K; i++)
            links[i] = links[i-1] + K;

        // re-sort by (a,b)
        qsort(cfts, n, sizeof(hcft_t), hcft_ab_cmpfunc);
        // merge duplicates (sum net weight)
        ncft = n;
        a = cfts[0].a;
        b = cfts[0].b;
        n = 0;
        for (i = 1; i < ncft; i++) {
            if (cfts[i].a == a && 
                cfts[i].b == b) {
                cfts[n].v += cfts[i].v;
                continue;
            }
            cfts[++n] = cfts[i];
            a = cfts[i].a;
            b = cfts[i].b;
        }
        ncft = ++n;

        // re-build offsets
        MYBZERO(coff, nseq);
        a = cfts[0].a >> 8;
        for (i = 1, j = 0; i < ncft; i++)
            if ((cfts[i].a >> 8) != a) { coff[a] = (uint64) j << 32 | (i - j); a = cfts[i].a >> 8; j = i; }
        coff[a] = (uint64) j << 32 | (ncft - j);

        // now process each group
        // we have two options here:
        /**
        {   // 1. simply process each group by size in descending order
            g = grps[0].g - 1;
            for (k = 0; k < K; k++)
                perms[g * K + k] = k;
            done[g] = 1;
            for (c = 1; c < ngrp; c++) {
                g = grps[c].g - 1;
                if (g < 0)
                    continue;
                MYBZERO(links[0], K2);
                cft = cfts + (coff[g] >> 32);
                nc  = (uint32) coff[g];
                for (k = 0; k < nc; k++) {
                    if (!done[cft[k].b >> 8])
                        continue;
                    ha = cft[k].a & 0xff;
                    hb = cft[k].b & 0xff;
                    hc = perms[(cft[k].b >> 8) * K + hb];
                    links[ha][hc] += cft[k].v;
                }

                hungarian(links, K, match);

                for (k = 0; k < K; k++)
                    perms[g * K + k] = match[k];
                done[g] = 1;
            }
        }
        **/
        {   // 2. a more sophisticated approach is to process related groups in a greedy manner
            // find the group with most links to other groups as the starting point
            ph_prior_t _pp, *pp = &_pp;
            pq_t *pq;
            pq_elem_t pe;
            int64 l;
            int8 *perm, *flag;
            int *list;

            pq = pq_init(ngrp, sizeof(ph_prior_t), ph_prior_cpyfunc, ph_prior_cmpfunc);
            MYMALLOC(list, ngrp);
            MYMALLOC(flag, ngrp);
            if (list == NULL || flag == NULL || pq == NULL)
                mem_alloc_error("haplotype partition arrays");

            for (c = 0; c < ngrp; c++) {
                g = grps[c].g - 1;
                if (g < 0)
                    continue;
                perm = perms + g * K;
                for (k = 0; k < K; k++)
                    perm[k] = k;
                cft = cfts + (coff[g] >> 32);
                nc  = (uint32) coff[g];
                l   = 0;
                for (k = 0; k < nc; k++)
                    l += cft[k].v;
                *pp = (ph_prior_t) {0, l, grps[c].t};
                pq_insert(pq, g, pp);
            }

            while (pq_pop1(pq, &pe)) {
                g = pe.id;
                // mark as done
                done[g] = 1;
                // update groups with links to this group
                MYBZERO(flag, ngrp);
                cft = cfts + (coff[g] >> 32);
                nc  = (uint32) coff[g];
                n   = 0;
                for (i = 0; i < nc; i++) {
                    b = cft[i].b >> 8;
                    if (done[b] || flag[b])
                        continue;
                    flag[b] = 1;
                    list[n++] = b;
                }
                
                for (i = 0; i < n; i++) {
                    a = list[i];
                    cft = cfts + (coff[a] >> 32);
                    nc  = (uint32) coff[a];
                    MYBZERO(links[0], K2);
                    for (k = 0; k < nc; k++) {
                        if (!done[cft[k].b >> 8])
                            continue;
                        ha = cft[k].a & 0xff;
                        hb = cft[k].b & 0xff;
                        hc = perms[(cft[k].b >> 8) * K + hb];
                        links[ha][hc] += cft[k].v;
                    }
                    hungarian(links, K, match);
                    for (k = 0; k < K; k++)
                        perms[a * K + k] = match[k];
                    l = 0;
                    for (k = 0; k < K; k++)
                        l += links[k][match[k]];
                    // 'a' must be in the pq
                    ph_prior_cpyfunc(pp, pq_prior(pq, a, 1));
                    pp->r = l;
                    pq_insert(pq, a, pp);
                }
            }

            free(list);
            free(flag);
            pq_destroy(pq);
        }

        {   // refine haplotype assignments
            ph_prior_t _pp, *pp = &_pp;
            pq_t *pq;
            pq_elem_t pe;
            int64 l;
            int8 *permx, *flag;
            int *list;
            
            pq = pq_init(ngrp, sizeof(ph_prior_t), ph_prior_cpyfunc, ph_prior_cmpfunc);
            MYMALLOC(permx, ngrp * K);
            MYMALLOC(list, ngrp);
            MYMALLOC(flag, ngrp);
            if (pq == NULL || permx == NULL || list == NULL || flag == NULL)
                mem_alloc_error("haplotype partition arrays");
            
            for (i = 0; i < ngrp; i++) {
                g = grps[i].g - 1;
                if (g < 0)
                    continue;
                cft = cfts + (coff[g] >> 32);
                nc  = (uint32) coff[g];
                MYBZERO(links[0], K2);
                for (k = 0; k < nc; k++) {
                    ha = cft[k].a & 0xff;
                    hb = cft[k].b & 0xff;
                    hc = perms[(cft[k].b >> 8) * K + hb];
                    links[ha][hc] += cft[k].v;
                }

                hungarian(links, K, match);

                l = 0;
                for (k = 0; k < K; k++)
                    l += links[k][match[k]] - links[k][k];
                if (l > 0) {
                    for (k = 0; k < K; k++)
                        permx[g * K + k] = match[k];
                    *pp = (ph_prior_t) {l, -i, 0};
                    pq_insert(pq, g, pp);
                }
            }

            while (pq_pop1(pq, &pe)) {
                g = pe.id;
                // update haplotype assignment
                for (k = 0; k < K; k++)
                    perms[g * K + k] = permx[g * K + k];
                // update groups with links to this group
                MYBZERO(flag, ngrp);
                cft = cfts + (coff[g] >> 32);
                nc  = (uint32) coff[g];
                n   = 0;
                for (i = 0; i < nc; i++) {
                    b = cft[i].b >> 8;
                    if (done[b] || flag[b])
                        continue;
                    flag[b] = 1;
                    list[n++] = b;
                }
                
                for (i = 0; i < n; i++) {
                    a = list[i];
                    cft = cfts + (coff[a] >> 32);
                    nc  = (uint32) coff[a];
                    MYBZERO(links[0], K2);
                    for (k = 0; k < nc; k++) {
                        ha = cft[k].a & 0xff;
                        hb = cft[k].b & 0xff;
                        hc = perms[(cft[k].b >> 8) * K + hb];
                        links[ha][hc] += cft[k].v;
                    }
                    hungarian(links, K, match);
                    l = 0;
                    for (k = 0; k < K; k++)
                        l += links[k][match[k]] - links[k][k];
                    if (l > 0) {
                        for (k = 0; k < K; k++)
                            permx[a * K + k] = match[k];
                        ph_prior_cpyfunc(pp, pq_prior(pq, a, 0));
                        pp->r = l;
                        pq_insert(pq, a, pp);
                    } else pq_remove(pq, a);
                }
            }

            free(permx);
            free(list);
            free(flag);
            pq_destroy(pq);
        }

        for (i = 0; i < nseq; i++) {
            g = haps[i].grp - 1;
            h = haps[i].hap - 1;
            if (g < 0 || h < 0)
                continue;
            haps[i].hap = perms[g * K + h] + 1;
        }

        // fix conflict table to reflect new haplotype assignments
        for (i = 0; i < ncft; i++) {
            a  = cfts[i].a >> 8;
            b  = cfts[i].b >> 8;
            ha = cfts[i].a & 0xff;
            hb = cfts[i].b & 0xff;
            cfts[i].a = (a<<8) | perms[a*K+ha];
            cfts[i].b = (b<<8) | perms[b*K+hb];
        }
        
        free(links[0]);
        free(links);
        free(done);
        free(match);
        free(perms);
    }

    {
        int cnts[ploidy+1];
        int64 lens[ploidy+1];
        for (i = 0; i <= ploidy; i++) {
            cnts[i] = 0;
            lens[i] = 0;
        }
        for (i = 0; i < nseq; i++) {
            k = haps[i].hap;
            cnts[k]++;
            lens[k] += haps[i].len;
        }
        fprintf(stderr, "[M::%s] Haplotype partition summary:\n", __func__);
        for (k = 0; k <= ploidy; k++)
            fprintf(stderr, "[M::%s]    H%d:  %6d seqs %12lld bp\n", __func__, k, cnts[k], lens[k]);
    }

    free(gseqs);
    free(grps);
    free(sord);
    free(hord);
    free(sidx);
    free(mols);
    free(mhls);
    free(cfts);
    free(coff);
    phase_data_destroy(data);
}

static void build_haplotype_partition_hybrid(ovl_t *ovls, int64 novl, hlk_t *hlks, int64 nhlk, int n_threads, 
    sdict_t *dicts, hap_info_t *haps, int ploidy)
{
    cft_t *cfts;
    ovl_t *ovl;
    uint64 *index;
    uint32 a;
    int64 i, j, nseq, ncft, intra_ovls, inter_ovls;
    double minhv, intra_hlks, inter_hlks;
    int *dels, *reps;
    
    // sequence number
    nseq = dicts->n;

    // allocate memory
    MYCALLOC(index, nseq);
    MYCALLOC(cfts, novl + nhlk);
    if (index == NULL || cfts == NULL)
        mem_alloc_error("global haplotype partition arrays");

    // mark repeat sequences and repeat-repeat overlaps up front:
    // both the signed partition and the conflict table exclude them.
    MYCALLOC(reps, nseq);
    dels = NULL;
    //dels = mark_repeat_overlaps(ovls, novl, dicts, ploidy, reps);
    
    // initial haplotype partition based on HiC attraction + overlap repulsion
    build_haplotype_partition_hybrid_core(ovls, novl, hlks, nhlk, n_threads, haps, dicts, reps, ploidy);

    intra_ovls = inter_ovls = 0;
    for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
        if (ovl->del)
            continue;
        if (haps[ovl->aread].hap > 0 && haps[ovl->aread].hap == haps[ovl->bread].hap)
            intra_ovls += (ovl->alen + ovl->blen) / 2;
        else inter_ovls += (ovl->alen + ovl->blen) / 2;
    }
    intra_hlks = inter_hlks = .0;
    for (i = 0; i < nhlk; i++) {
        if (haps[hlks[i].a].hap == haps[hlks[i].b].hap)
            intra_hlks += hlks[i].v;
        else inter_hlks += hlks[i].v;
    }
    fprintf(stderr, "[M::%s] initial greedy haplotype partition\n", __func__);
    fprintf(stderr, "[M::%s] intra_hlks: %16.2f inter_hlks: %16.2f\n", __func__, intra_hlks, inter_hlks);
    fprintf(stderr, "[M::%s] intra_ovls: %16lld inter_ovls: %16lld\n", __func__, intra_ovls, inter_ovls);
/***
    // build the conflict table for local-search refinement:
    // non-repeat overlaps (hard separation floor) plus HiC links.
    ncft = 0;
    for (i = 0, ovl = ovls; i < novl; i++, ovl++)
        if (!ovl->del && !(reps[ovl->aread] || reps[ovl->bread]))
            cfts[ncft++] = (cft_t) {
                ovl->aread,
                ovl->bread,
                (ovl->alen + ovl->blen) / 2,
                CFT_OVL
            };

    // add hic links to the conflict table
    minhv = .001;
    for (i = 0; i < nhlk; i++) {
        if (hlks[i].v < minhv)
            minhv = hlks[i].v;
    }
    for (i = 0; i < nhlk; i++)
        cfts[ncft++] = (cft_t) {
            hlks[i].a, 
            hlks[i].b, 
            (uint32) (hlks[i].v / minhv + 1),
            CFT_HIC
        };
    // sort and build index for quick search
    qsort(cfts, ncft, sizeof(cft_t), cft_abseqs_cmpfunc);
    a = cfts->a;
    for (i = 1, j = 0; i < ncft; i++) {
        if (cfts[i].a != a) {
            index[a] = (uint64) j << 32 | (i - j);
            a = cfts[i].a;
            j = i;
        }
    }
    index[a] = (uint64) j << 32 | (i - j);

    // TODO add an extra step to group repeat sequences

    // refine haplotype partition using hic links
    refine_haplotype_partition_hic_2opt(cfts, ncft, index, haps, nseq, ploidy, vns_data);

    refine_haplotype_partition_hic_kopt(cfts, ncft, index, haps, nseq, ploidy, vns_data);

    intra_ovls = inter_ovls = 0;
    for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
        if (!dels[i] && ovl->del)
            continue;
        if (haps[ovl->aread].hap == haps[ovl->bread].hap)
            intra_ovls += (ovl->alen + ovl->blen) / 2;
        else inter_ovls += (ovl->alen + ovl->blen) / 2;
    }
    intra_hlks = inter_hlks = .0;
    for (i = 0; i < nhlk; i++) {
        if (haps[hlks[i].a].hap == haps[hlks[i].b].hap)
            intra_hlks += hlks[i].v;
        else inter_hlks += hlks[i].v;
    }
    fprintf(stderr, "[M::%s] refined haplotype partition\n", __func__);
    fprintf(stderr, "[M::%s] intra_hlks: %16.2f inter_hlks: %16.2f\n", __func__, intra_hlks, inter_hlks);
    fprintf(stderr, "[M::%s] intra_ovls: %16lld inter_ovls: %16lld\n", __func__, intra_ovls, inter_ovls);

if(0) {
    //char *fs[2] = {"Phase_Test/drRosCani1_ec-20260618.chr1.raw.agp",
    //               "Phase_Test/drRosCani1_ec-20260618.chr1.fix.agp"};
    //char *fs[2] = {"Phase_Test/drRosCani1_ec-20260622.bbseq.scf_scaffolds_final.all_haps.raw.agp",
    //               "Phase_Test/drRosCani1_ec-20260622.bbseq.scf_scaffolds_final.all_haps.fix.agp"};
    //char *fs[1] = {"Phase_Test/drRosCani1_mc.curated.renamed.agp"};
    char *fs[1] = {"Phase_Test/drRosAgre1_mc.curated.renamed.agp"};
    for (int f = 0; f < sizeof(fs) / sizeof(fs[0]); f++) {
        iostream_t *fp = iostream_open(fs[f]);
        char *line;
        char *fields[9];
        while ((line = iostream_getline(fp)) != NULL) {
            if (is_empty_line(line) || parse_line(line, fields, 9) < 9)
                continue;
            if (strcmp(fields[4], "W") == 0)
                haps[sd_get(dicts, fields[5])].hap = atoi(fields[0]+strlen(fields[0])-1);
        }
        intra_ovls = inter_ovls = 0;
        for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
            if (!dels[i] && ovl->del)
                continue;
            if (haps[ovl->aread].hap && haps[ovl->aread].hap == haps[ovl->bread].hap)
                intra_ovls += (ovl->alen + ovl->blen) / 2;
            else inter_ovls += (ovl->alen + ovl->blen) / 2;
        }
        intra_hlks = inter_hlks = .0;
        for (i = 0; i < nhlk; i++) {
            if (haps[hlks[i].a].hap && haps[hlks[i].a].hap == haps[hlks[i].b].hap)
                intra_hlks += hlks[i].v;
            else inter_hlks += hlks[i].v;
        }
        fprintf(stderr, "[M::%s] XXX %d\n", __func__, f+1);
        fprintf(stderr, "[M::%s] intra_hlks: %16.2f inter_hlks: %16.2f\n", __func__, intra_hlks, inter_hlks);
        fprintf(stderr, "[M::%s] intra_ovls: %16lld inter_ovls: %16lld\n", __func__, intra_ovls, inter_ovls);
        iostream_close(fp);
    }

    exit(0);
}
**/
    // restore the deleted overlaps
    if (dels) {
        for (i = 0, ovl = ovls; i < novl; i++, ovl++)
            if (dels[i])
                ovl->del = 0;
    }

    free(index);
    free(cfts);
    free(dels);
    free(reps);

    return;
}

static void build_haplotype_partition_overlap(ovl_t *ovls, int64 novl, sdict_t *dicts, hap_info_t *haps, int ploidy)
{
    cft_t *cfts;
    ovl_t *ovl;
    uint64 *index;
    vns_data_t *vns_data;
    uint32 a;
    int64 i, j, nseq, ncft;
    int64 intra_ovls, inter_ovls;
    
    // sequence number
    nseq = dicts->n;
    
    // allocate memory  
    MYCALLOC(index, nseq);
    MYCALLOC(cfts, novl);
    if (index == NULL || cfts == NULL)
        mem_alloc_error("global haplotype partition arrays");

    vns_data = vns_data_init(nseq, ploidy);

    // add overlaps to the conflict table
    // overlaps are sorted by abseqs
    ncft = 0;
    for (i = 0, ovl = ovls; i < novl; i++, ovl++)
        if (!ovl->del)
            cfts[ncft++] = (cft_t) {
                ovl->aread,
                ovl->bread,
                (ovl->alen + ovl->blen) / 2,
                CFT_OVL
            };
    // sort and build index for quick search
    qsort(cfts, ncft, sizeof(cft_t), cft_abseqs_cmpfunc);
    a = cfts->a;
    for (i = 1, j = 0; i < ncft; i++) {
        if (cfts[i].a != a) {
            index[a] = (uint64) j << 32 | (i - j);
            a = cfts[i].a;
            j = i;
        }
    }
    index[a] = (uint64) j << 32 | (i - j);

    // build global haplotype partition
    build_haplotype_partition_ovl_greedy(cfts, index, haps, nseq, ploidy);

    intra_ovls = inter_ovls = 0;
    for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
        if (ovl->del)
            continue;
        if (haps[ovl->aread].hap == haps[ovl->bread].hap)
            intra_ovls += (ovl->alen + ovl->blen) / 2;
        else inter_ovls += (ovl->alen + ovl->blen) / 2;
    }
    fprintf(stderr, "[M::%s] initial greedy haplotype partition\n", __func__);
    fprintf(stderr, "[M::%s] intra_ovls: %16lld inter_ovls: %16lld\n", __func__, intra_ovls, inter_ovls);

    // refine haplotype partition using overlaps
    refine_haplotype_partition_ovl_vns(cfts, ncft, index, haps, nseq, ploidy, vns_data);
    
    intra_ovls = inter_ovls = 0;
    for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
        if (ovl->del)
            continue;
        if (haps[ovl->aread].hap == haps[ovl->bread].hap)
            intra_ovls += (ovl->alen + ovl->blen) / 2;
        else inter_ovls += (ovl->alen + ovl->blen) / 2;
    }
    fprintf(stderr, "[M::%s] refined haplotype partition\n", __func__);
    fprintf(stderr, "[M::%s] intra_ovls: %16lld inter_ovls: %16lld\n", __func__, intra_ovls, inter_ovls);

    free(index);
    free(cfts);
    vns_data_destroy(vns_data);

    return;
}

#define MIN_VAR_SIZE     100000
#define MIN_VAR_FRAC     0.1
#define MIN_LOW_COV_FRAC 0.7

typedef struct {
    uint32 seq:31, dir:1, pos; // dir stores the diverged side, 0 for left and 1 for right
    ovl_t *ovl;
} sv_t;

static int sv_cmpfunc(const void *a, const void *b)
{
    const sv_t *sa = a, *sb = b;
    if (sa->seq != sb->seq)
        return (sa->seq > sb->seq) - (sa->seq < sb->seq);
    if (sa->dir != sb->dir)
        return (sa->dir > sb->dir) - (sa->dir < sb->dir);
    return (sa->pos > sb->pos) - (sb->pos < sa->pos);
}

static void detect_structural_variants(ovl_t *ovls, int64 novl, sdict_t *dicts, int ploidy)
{
    kvec_t(sv_t) vars;
    sv_t *var;
    uint64 *index;
    ovl_t *ovl;
    cov_point_t **covs;
    uint32 dir;
    int64 i, j, k, nseq;
    int abpos, aepos, bbpos, bepos, alen, blen, acov, bcov;
    int a, b, n, del, sup, nvar, movl, cab, cae, cbb, cbe, max_cov, *ncov, *cpts;
    double min_qual = OVL_MIN_QUAL;
    
    // sequence number
    nseq = dicts->n;

    // allocate memory
    MYCALLOC(index, nseq);
    MYMALLOC(covs, nseq);
    MYMALLOC(ncov, nseq);
    MYMALLOC(covs[0], (novl+nseq)*2);
    if (index == NULL || covs == NULL || ncov == NULL || covs[0] == NULL)
        mem_alloc_error("SV detection arrays");

    // build overlap index for quick search
    // overlaps are sorted by abseqs
    a = ovls->aread;
    for (i = 1, j = 0; i < novl; i++) {
        if (ovls[i].aread != a) {
            index[a] = (uint64) j << 32 | (i - j);
            a = ovls[i].aread;
            j = i;
        }
    }
    index[a] = (uint64) j << 32 | (i - j);

    // max overlaps
    movl = 0;
    for (i = 0; i < nseq; i++) {
        n = index[i] >> 32;
        if (movl < n)
            movl = n;
    }
    MYMALLOC(cpts, movl*2);
    if (cpts == NULL)
        mem_alloc_error("coverage points");
    // calculate sequence coverage
    for (i = 0; i < nseq; i++) {
        if (i)
            covs[i] = covs[i-1] + ncov[i-1];
        ovl = ovls + (index[i] >> 32);
        n = (uint32) index[i];
        ncov[i] = calc_coverage_from_intervals(covs[i], cpts, dicts->s[i].len, ovl, n, &min_qual, pts_from_overlaps);
    }
    free(cpts);


    // maximum coverage in overlapping regions
    max_cov = ploidy - 1;

    // potential structural variants are identified INTL overlaps
    kv_init(vars);
    for (i = 0, ovl = ovls; i < novl; i++, ovl++) {
        if (ovl->del || !(ovl->type & OVL_INTL))
            continue;
        a = ovl->aread;
        b = ovl->bread;
        alen = dicts->s[a].len;
        blen = dicts->s[b].len;
        abpos = ovl->abpos;
        aepos = ovl->aepos;
        bbpos = ovl->bbpos;
        bepos = ovl->bepos;

        acov = total_range_coverage(covs[a], ncov[a], abpos, aepos, max_cov);
        bcov = total_range_coverage(covs[b], ncov[b], bbpos, bepos, max_cov);

        if ((double) ovl->alen * ovl->qual < MIN_VAR_SIZE || 
            (double) ovl->alen * ovl->qual / alen < MIN_VAR_FRAC ||
            (double) ovl->blen * ovl->qual < MIN_VAR_SIZE || 
            (double) ovl->blen * ovl->qual / blen < MIN_VAR_FRAC ||
            (double) acov / ovl->alen < MIN_LOW_COV_FRAC ||
            (double) bcov / ovl->blen < MIN_LOW_COV_FRAC)
            continue;

        cab = ((double) abpos / alen > MIN_VAR_FRAC && abpos > MIN_VAR_SIZE);
        cae = ((double) (alen - aepos) / alen > MIN_VAR_FRAC && (alen - aepos) > MIN_VAR_SIZE);
        cbb = ((double) bbpos / blen > MIN_VAR_FRAC && bbpos > MIN_VAR_SIZE);
        cbe = ((double) (blen - bepos) / blen > MIN_VAR_FRAC && (blen - bepos) > MIN_VAR_SIZE);

        if (ovl->arev == ovl->brev) {
            if (cab && cbb)
                kv_push(sv_t, vars, ((sv_t) {a, 0, abpos, ovl}));
            if (cae && cbe)
                kv_push(sv_t, vars, ((sv_t) {a, 1, aepos, ovl}));
        } else {
            if (cab && cbe)
                kv_push(sv_t, vars, ((sv_t) {a, 0, abpos, ovl}));
            if (cae && cbb)
                kv_push(sv_t, vars, ((sv_t) {a, 1, aepos, ovl}));
        }
    }

    if (!vars.n) goto cleanup;

    // sort by sequence and position
    qsort(vars.a, vars.n, sizeof(sv_t), sv_cmpfunc);

    // check each potential SV
    // we require no other overlaps spanning the breakpoints
    // also there should be overlaps close to the breakpoints in the diverged side to support the SV
    // we do a naive search here
    nvar = 0;
    var = vars.a;
    for (i = 1, j = 0; i <= vars.n; i++) {
        if (i == vars.n || 
            var[i].seq != var[i-1].seq || 
            var[i].dir != var[i-1].dir || 
            var[i].pos  > var[i-1].pos + MIN_VAR_SIZE) {
            // use the middle position
            a     = var[j].seq;
            dir   = var[j].dir;
            abpos = var[(j+i-1)/2].pos;
            ovl = ovls + (index[a] >> 32);
            n = (uint32) index[a];
            del = sup = 0;
            for (k = 0; k < n; k++, ovl++) {
                b = ovl->bread;
                alen = dicts->s[a].len;
                blen = dicts->s[b].len;
                if (ovl->del ||
                    (double) ovl->alen * ovl->qual < MIN_VAR_SIZE || 
                    (double) ovl->alen * ovl->qual / alen < MIN_VAR_FRAC ||
                    (double) ovl->blen * ovl->qual < MIN_VAR_SIZE || 
                    (double) ovl->blen * ovl->qual / blen < MIN_VAR_FRAC)
                    continue;
                // check if the overlap spans the breakpoint
                if (abpos > ovl->abpos + MIN_VAR_SIZE / 2 && 
                    abpos < ovl->aepos - MIN_VAR_SIZE / 2) {
                    del = 1;
                    break;
                }
                // check if there are overlaps close to the breakpoint in the diverged side
                if (!sup) {
                    sup = dir ? 
                        (ovl->abpos > abpos - MIN_VAR_SIZE / 2 && ovl->abpos < abpos + MIN_VAR_SIZE / 2) : 
                        (ovl->aepos > abpos - MIN_VAR_SIZE / 2 && ovl->aepos < abpos + MIN_VAR_SIZE / 2) ;
                }
            }
            if (!del && sup) var[nvar++] = (sv_t) {a, dir, abpos, var[(j+i-1)/2].ovl};
            j = i;
        }
    }
    vars.n = nvar;

    // debug info
#ifdef DEBUG_SV_DETECTION
    int64 slen = 0;
    n = 0;
    for (i = 0; i < vars.n; i++) {
        if (i == 0 || vars.a[i].seq != vars.a[i-1].seq) {
            slen += dicts->s[vars.a[i].seq].len;
            n += 1;
        }
    }
    fprintf(stderr, "[M::%s] detected %d potential SVs in %d sequences of %lld bp\n", __func__, (int) vars.n, n, slen);
    for (i = 0; i < vars.n; i++) {
        ovl = vars.a[i].ovl;
        fprintf(stderr, "[M::%s] potential SV: %s [%c] %10d %s %10d %10d %10d %10d\n", __func__, 
            dicts->s[vars.a[i].seq].name, "LR"[vars.a[i].dir], vars.a[i].pos, 
            dicts->s[ovl->bread].name, ovl->abpos, ovl->aepos, ovl->bbpos, ovl->bepos);
    }
#endif

cleanup:
    free(covs[0]);
    free(covs);
    free(ncov);
    free(index);
    kv_destroy(vars);
}

static void write_yahs_outputs(char *hic_bfile, hap_info_t *haps, int64 nseq, int ploidy, sdict_t *dicts, asm_dict_t *break_dict, char *out_pref)
{
    if (!hic_bfile) return;
    
    kvec_t(int32) bvecs;
    hap_info_t *hseq, *hseqs;
    sdict_t *bdicts;
    sd_seg_t *cseg;
    int64 nscf, *pos, **ends;
    int32 *seg, *seq, *nend, **segs, **seqs;
    FILE *fp1, *fp2, *fp3;
    char *sbuff;
    int64 a, b, e, s, n, m, l, t, r, g, nbsq, slen;
    int i, j, k;
    
    // find backbone sequence to build a pseudo-scaffold for each block
    // then build a map sequences to the backone position
    MYMALLOC(hseqs, nseq);
    if (hseqs == NULL)
        mem_alloc_error("yahs output arrays");
    
    memcpy(hseqs, haps, nseq * sizeof(hap_info_t));
    qsort(hseqs, nseq, sizeof(hap_info_t), hap_info_gcmpfunc);
    
    kv_init(bvecs);
    g = hseqs[0].grp;
    for (i = 1, j = 0; i <= nseq; i++) {
        if (i == nseq || hseqs[i].grp != g) {
            b = hseqs[j].bpos;
            e = hseqs[j].epos;
            for (k = j; k < i; k++)
                if (hseqs[k].epos > e)
                    e = hseqs[k].epos;
            s = j;
            while (b < e) {
                k = -1;
                l =  b;
                while (s < i && hseqs[s].bpos <= b) {
                    if (hseqs[s].epos > l) {
                        k = s;
                        l = hseqs[s].epos;
                    }
                    s++;
                }
                if (k < 0) break;
                kv_push(int32, bvecs, (int32) hseqs[k].seq);
                b = l;
            }

            if (i < nseq)
                g = hseqs[i].grp;
            j = i;
        }
    }

    nscf = hseqs[nseq-1].grp+1;
    nbsq = bvecs.n;

    MYMALLOC(nend, nscf);
    MYMALLOC(ends, nscf);
    MYMALLOC(ends[0], nbsq);
    MYMALLOC(seqs, nscf);
    MYMALLOC(seqs[0], nbsq);
    MYMALLOC(sbuff, strlen(out_pref)+32);
    if (nend == NULL || ends == NULL || ends[0] == NULL || 
        seqs == NULL || seqs[0] == NULL || sbuff == NULL)
        mem_alloc_error("yahs output arrays");
    sprintf(sbuff, "%s.bbseq.agp", out_pref);
    fp1 = fopen(sbuff, "w");
    sprintf(sbuff, "%s.bbscf.agp", out_pref);
    fp2 = fopen(sbuff, "w");
    sprintf(sbuff, "%s.bbpos.txt", out_pref);
    fp3 = fopen(sbuff, "w");
    if (fp1 == NULL || fp2 == NULL || fp3 == NULL) {
        fprintf(stderr, "[E::%s] failed to open output file %s\n", __func__, sbuff);
        exit(EXIT_FAILURE);
    }
    bdicts = sd_init();
    pos = ends[0];
    seq = seqs[0];
    m = 0;
    g = haps[bvecs.a[0]].grp;
    for (i = 1, j = 0; i <= nbsq; i++) {
        if (i == nbsq || haps[bvecs.a[i]].grp != g) {
            ends[g] = pos;
            seqs[g] = seq;
            t = 0;
            l = 0;
            n = 0;
            for (k = j; k < i; k++) {
                if (k > j) { // add gap between sequencesmake
                    write_agp_gap(fp2, sbuff, slen + 1, slen + LG_AGP_GAP_SIZE, ++t);
                    l += LG_AGP_GAP_SIZE;
                }
                a = bvecs.a[k];
                r = haps[a].rev;
                pos[n] = haps[a].epos;
                seq[n] = m++;
                n++;
                // add DICT and AGP entries
                sprintf(sbuff, "s%lld-%lld", g, n);
                slen = n>1? (pos[n-1]-pos[n-2]) : pos[0];
                sd_put(bdicts, sbuff, slen);
                if (r) {
                    e = slen;
                    b = 0;
                } else {
                    e = dicts->s[a].len;
                    b = e - slen;
                }
                if (break_dict) {
                    // compose back onto the raw sequence coordinates
                    cseg = &break_dict->seg[break_dict->s[a].s];
                    write_agp_seq(fp1, sbuff, 1, slen, 1, break_dict->sdict->s[cseg->c>>1].name, cseg->x + b + 1, cseg->x + e, r);
                } else
                    write_agp_seq(fp1, sbuff, 1, slen, 1, dicts->s[a].name, b + 1, e, r);
                sprintf(sbuff, "scf%lld", g);
                write_agp_seq(fp2, sbuff, l+1, l+slen, ++t, bdicts->s[m-1].name, 1, slen, 0);
                l += slen;
            }
            nend[g] = n;
            pos += n;
            seq += n;

            if (i < nbsq)
                g = haps[bvecs.a[i]].grp;
            j = i;
        }
    }
    if (break_dict)
        for (i = 0; i < nseq; i++) {
            hseq = hseqs + i;
            a = hseq->seq;
            cseg = &break_dict->seg[break_dict->s[a].s];
            fprintf(fp3, "%s\t%d\t%c\t%d\t%d\t%lld\n",
                break_dict->sdict->s[cseg->c >> 1].name,
                hseq->len,
                "+-"[hseq->rev],
                hseq->grp,
                hseq->hap,
                cseg->x + hseq->bpos);
        }
    else
        for (i = 0; i < nseq; i++) {
            hseq = hseqs + i;
            fprintf(fp3, "%s\t%d\t%c\t%d\t%d\t%lld\n", 
                dicts->s[hseq->seq].name,
                hseq->len, 
                "+-"[hseq->rev], 
                hseq->grp, 
                hseq->hap, 
                hseq->bpos);
        }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);

    // map sequences to backbone position
    // each sequence will have a segment block, where
    // seg[0]       = sequence length
    // seg[1]       = sequence orientation
    // seg[2]       = number of segment pairs
    // seg[2*i+3]   = end position on sequence
    // seg[2*i+4]   = backbone sequence index
    // to map a position {p} on sequence {i} to backbone position
    // 1. find the first segment pair {s} where seg[2*s+3] > p
    // 2. seg[2*s+4] is the backbone sequence index
    // 3. seg[2*s+3] - p is the distance to the end of the backbone sequence
    n = 0;
    for (i = 0; i < nseq; i++)
        n += nend[haps[i].grp] * 2 + 3;
    MYMALLOC(segs, nseq);
    MYMALLOC(segs[0], n);
    if (segs == NULL || segs[0] == NULL)
        mem_alloc_error("yahs output arrays");
    seg = segs[0];
    for (i = 0; i < nseq; i++) {
        segs[i] = seg;
        seg[0] = haps[i].len;
        seg[1] = haps[i].rev;
        seg += 3;
        // find the first backbone sequence in the block
        g = haps[i].grp;
        pos = ends[g];
        seq = seqs[g];
        n = nend[g];
        b = haps[i].bpos;
        for (j = 0; j < n; j++)
            if (pos[j] > b)
                break;
        e = b + haps[i].len;
        m = 0;
        s = 0;
        while (b < e && j < n) {
            slen = pos[j] - b;
            s += slen;
            seg[m++] = s;
            seg[m++] = seq[j];
            b += slen;
            j++;
        }
        seg[-1] = m / 2;
        seg += m;
    }
    
    // now write binary file for hic links
    sprintf(sbuff, "%s.bbscf-hic.bin", out_pref);
    write_binary_hic_data_pseudo_yahs(hic_bfile, bdicts, break_dict, segs, sbuff);

    free(ends[0]);
    free(ends);
    free(seqs[0]);
    free(seqs);
    free(segs[0]);
    free(segs);
    free(nend);
    free(sbuff);
    free(hseqs);
    kv_destroy(bvecs);
    sd_destroy(bdicts);

    return;
}

scf_t *build_pseudo_scaffolds(ovl_t *ovls, int64 novl, sdict_t *dicts, asm_dict_t *break_dict, busco_table_t *buscos,
    int ploidy, int min_ext, int min_qual, int n_threads, int conf_yahs, char *hic_bfile, 
    char *out_pref, int64 *_nscf)
{
    if (_nscf) *_nscf = 0; 
    if (ovls == NULL || !novl)
        return 0;

    kvec_t(scf_t) scfs;
    scf_t *scf;
    hlk_t *hlks;
    hap_info_t *haps;
    kvec_t(uint32) ctgs;
    kvec_t(uint64) segs;
    int64 i, nseq, nhlk, slen;
    int h, g, nctg;

    // sequence number
    nseq = dicts->n;

    // build linkage map from HiC
    nhlk = 0;
    hlks = build_hic_linkage_map(hic_bfile, min_qual, dicts, break_dict, ovls, novl, &nhlk);

    // initial haplotype group partition
    MYCALLOC(haps, nseq);
    if (haps == NULL)
        mem_alloc_error("haplotype partition array");
    for (i = 0; i < nseq; i++) {
        haps[i].seq = i;
        haps[i].len = dicts->s[i].len;
    }

    // build haplotype partitioned scaffolds
    if (hlks && nhlk) {
        // build scaffold blocks (without haplotype partition)
        build_scaffold_partition(ovls, novl, haps, dicts, buscos, ploidy, min_ext, 0);
        // build haplotype partition using overlaps and HiC links
        build_haplotype_partition_hybrid(ovls, novl, hlks, nhlk, n_threads, dicts, haps, ploidy);
        // rebuild scaffold blocks (with haplotype partition)
        build_scaffold_partition(ovls, novl, haps, dicts, buscos, ploidy, min_ext, 1);
    } else {
        // build haplotype partition using overlaps only
        build_haplotype_partition_overlap(ovls, novl, dicts, haps, ploidy);
        // build scaffold blocks (with haplotype partition)
        build_scaffold_partition(ovls, novl, haps, dicts, buscos, ploidy, min_ext, 1);
    }

    // write outputs for yahs
    if (conf_yahs)
        write_yahs_outputs(hic_bfile, haps, nseq, ploidy, dicts, break_dict, out_pref);

    // build data structure for return
    kv_init(scfs);
    kv_init(ctgs);
    kv_init(segs);
    
    qsort(haps, nseq, sizeof(hap_info_t), hap_info_acmpfunc);
    
    h = haps[0].hap;
    g = haps[0].grp;
    slen = 0;
    ctgs.n = segs.n = 0;
    for (i = 0; i <= nseq; i++) {
        if (i == nseq || haps[i].hap != h || haps[i].grp != g) {
            kv_pushp(scf_t, scfs, &scf);
            nctg = ctgs.n;
            MYMALLOC(scf->ctgs, nctg);
            MYMALLOC(scf->segs, nctg);
            scf->grp = g-1;
            scf->hap = h-1;
            scf->len = slen;
            scf->nctg = nctg;
            memcpy(scf->ctgs, ctgs.a, sizeof(uint32) * nctg);
            memcpy(scf->segs, segs.a, sizeof(uint64) * nctg);
            if (i < nseq) {
                h = haps[i].hap;
                g = haps[i].grp;
                slen = 0;
                ctgs.n = segs.n = 0;
            } else break;
        }
        kv_push(uint32, ctgs, haps[i].seq << 1 | haps[i].rev);
        kv_push(uint64, segs, haps[i].len);
        slen += haps[i].len;
    }
    
    // haplotype group summary stats
    haplotype_group_sequence_summary(scfs.a, scfs.n, dicts, buscos, ploidy, 1, stderr);
    
    kv_destroy(ctgs);
    kv_destroy(segs);
    free(hlks);
    free(haps);
    *_nscf = scfs.n;

    return scfs.a;
}

void write_scf_outputs(scf_t *scfs, int nscf, sdict_t *dicts, asm_dict_t *break_dict, const uint8 opts_out, const char *pref_out)
{
    FILE *fo;
    char *file, *name, *sname;
    int i, j, t, grp, hap, nctg;
    uint32 *ctgs, rid, rev, beg, end, len;
    uint64 *segs, slen;
    scf_t *scf;
    sd_seg_t *cseg;

    MYMALLOC(file, strlen(pref_out) + 35);
    MYMALLOC(name, 128);

    if (opts_out & AGP_OUT) {
        sprintf(file, "%s.grp.agp", pref_out);
        fo = fopen(file, "w");
        for (i = 0; i < nscf; i++) {
            scf = &scfs[i]; 
            grp = scf->grp+1;
            hap = scf->hap+1;
            nctg = scf->nctg;
            ctgs = scf->ctgs;
            segs = scf->segs;
            sprintf(name, "scaffold_%d_h%d", grp, hap);

            slen = 0;
            t = 0;
            for (j = 0; j < nctg; j++) {
                rid = ctgs[j] >> 1;
                rev = ctgs[j] & 1;
                beg = segs[j] >> 32;
                len = (uint32) segs[j];
                end = beg + len;
                if (break_dict) {
                    // compose back onto the raw sequence coordinates
                    cseg = &break_dict->seg[break_dict->s[rid].s];
                    write_agp_seq(fo, name, slen + 1, slen + len, ++t, break_dict->sdict->s[cseg->c>>1].name, cseg->x + beg + 1, cseg->x + end, rev);
                } else
                    write_agp_seq(fo, name, slen + 1, slen + len, ++t, dicts->s[rid].name, beg + 1, end, rev);
                slen += len;
                if (j != nctg - 1) {
                    write_agp_gap(fo, name, slen + 1, slen + LG_AGP_GAP_SIZE, ++t);
                    slen += LG_AGP_GAP_SIZE;
                }
            }
        }
        fclose(fo);
    }

    if (opts_out & GRP_OUT) {
        sprintf(file, "%s.grp.txt", pref_out);
        fo = fopen(file, "w");
        for (i = 0; i < nscf; i++) {
            scf = &scfs[i]; 
            grp = scf->grp+1;
            hap = scf->hap+1;
            nctg = scf->nctg;
            ctgs = scf->ctgs;
            segs = scf->segs;
            for (j = 0; j < nctg; j++) {
                sname = dicts->s[ctgs[j]>>1].name;
                if (break_dict) {
                    // compose back onto the raw sequence coordinates
                    cseg = &break_dict->seg[break_dict->s[ctgs[j]>>1].s];
                    sname = break_dict->sdict->s[cseg->c>>1].name;
                }
                fprintf(fo, "%s %8d %8d %8d %10u\n", sname, grp, hap, i+1, (uint32) segs[j]);
            }
        }
        fclose(fo);
    }

    if (opts_out & PLT_OUT) {
        sprintf(file, "%s.plt.txt", pref_out);
        fo = fopen(file, "w");
        for (i = 0; i < nscf; i++) {
            scf = &scfs[i]; 
            grp = scf->grp+1;
            hap = scf->hap+1;
            nctg = scf->nctg;
            ctgs = scf->ctgs;
            segs = scf->segs;
            for (j = 0; j < nctg; j++) {
                sname = dicts->s[ctgs[j]>>1].name;
                if (break_dict) {
                    // compose back onto the raw sequence coordinates
                    cseg = &break_dict->seg[break_dict->s[ctgs[j]>>1].s];
                    sname = break_dict->sdict->s[cseg->c>>1].name;
                }
                fprintf(fo, "@%s%c\tLG%d.H%d %10u\n", sname, "+-"[ctgs[j]&1], grp, hap, (uint32) segs[j]);
            }
        }
        fclose(fo);
    }

    free(file);
    free(name);
}

