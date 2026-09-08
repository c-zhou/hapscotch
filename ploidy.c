/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2024 Chenxi Zhou <chnx.zhou@gmail.com>                          *
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
#include <math.h>

#include "kvec.h"
#include "sdict.h"
#include "misc.h"
#include "range.h"

#include "ploidy.h"

// upper limit of genome ploidy number to consider, settable via -P/--max-ploidy
int MAX_PLOIDY_NUMBER = 16;

int64 *genome_coverage_histogram(aln_t *alns, int64 naln, sdict_t *dicts, int max_copy, int print_hist)
{
    if (alns == NULL || naln <= 0)
        return NULL;

    if (max_copy <= 0)
        max_copy = MAX_PLOIDY_NUMBER;

    int64 i, j, k, n, amax, ccnt, copy, nseq, seqtot, alntot, *hist;
    uint32 read, prev;
    chord_t *chord;
    point_info_t *point;

#ifdef DEBUG_GENOME_COV_HIST
    int cmax_s, cmax_g;
    kvec_t(uint32) covs_s, covs_g;
    kv_init(covs_s);
    kv_resize(uint32, covs_s, 32);
    kv_init(covs_g);
    kv_resize(uint32, covs_g, 32);
    MYBZERO(covs_g.a, covs_g.m);
    cmax_g = 0;
#endif
    
    nseq = dicts->n;
    read = alns->aread;
    amax = 0;
    for (i = 0, j = 0; i < naln; i++) {
        if (alns[i].aread != read) {
            if (amax < i-j)
                amax = i-j;
            read = alns[i].aread;
            j = i;
        }
    }
    if (amax < i-j) amax = i-j;

    MYMALLOC(chord, amax*2+2);
    MYMALLOC(point, amax*2+1);
    if (chord == NULL || point == NULL)
        mem_alloc_error("chord/space array");

    MYCALLOC(hist, max_copy+1);
    if (hist == NULL)
        mem_alloc_error("hist array");

    alntot = 0;
    prev = 0;
    read = alns->aread;
    for (i = 0, j = 0; i <= naln; i++) {
        if (i == naln || alns[i].aread != read) {
            // sequences with no alignments
            for (k = prev; k < read; k++) {
                hist[0] += dicts->s[k].len;
#ifdef DEBUG_GENOME_COV_HIST
                fprintf(stdout, "%s\t%d\t%d\t%d\t%f\tPG:A:%c\n", dicts->s[k].name, 
                        0, dicts->s[k].len, dicts->s[k].len, 1., PG_HIST);
                covs_g.a[0] += dicts->s[k].len;
#endif
#ifdef DEBUG_GENOME_COV_DEPTH
                fprintf(stdout, "%s\t%d\t%d\t0\tPG:A:%c\n", dicts->s[k].name, 
                        0, dicts->s[k].len, PG_DEPTH);
#endif
            }
            // make chord
            n = (i - j) * 2;
            for (k = 0; k < n; k++)
                point[k].next = &point[k+1];
            point[n].next = NULL;
            ccnt = make_chord(&alns[j], i - j, dicts->s[read].len, chord, point);
            for (k = 1; k < ccnt; k++) {
                copy = chord[k-1].wgt;
                if (copy > max_copy)
                    copy = max_copy;
                hist[copy] += chord[k].pos - chord[k-1].pos;
                alntot += (chord[k].pos - chord[k-1].pos) * copy;
            }
            if (!ccnt)
                hist[0] += dicts->s[read].len;
            
#ifdef DEBUG_GENOME_COV_HIST
            MYBZERO(covs_s.a, covs_s.m);
            cmax_s = 0;
            for (k = 1; k < ccnt; k++) {
                copy = chord[k-1].wgt;
                if (copy >= covs_s.m) {
                    MYREALLOC(covs_s.a, copy * 2);
                    MYBZERO(covs_s.a + covs_s.m, copy * 2 - covs_s.m);
                    covs_s.m = copy * 2;
                }
                covs_s.a[copy] += chord[k].pos - chord[k-1].pos;
                if (copy > cmax_s) cmax_s = copy;
            }
            if (cmax_s > cmax_g) cmax_g = cmax_s;
            if (cmax_s > covs_g.m) {
                MYREALLOC(covs_g.a, cmax_s * 2);
                MYBZERO(covs_g.a + covs_g.m, cmax_s * 2 - covs_g.m);
                covs_g.m = cmax_s * 2;
            }
            for (k = 0; k <= cmax_s; k++) {
                if (covs_s.a[k] > 0) {
                    fprintf(stdout, "%s\t%lld\t%d\t%d\t%f\tPG:A:%c\n", dicts->s[read].name, k, covs_s.a[k], 
                            dicts->s[read].len, (double) covs_s.a[k]/dicts->s[read].len, PG_HIST);
                    covs_g.a[k] += covs_s.a[k];
                }
            }
#endif
#ifdef DEBUG_GENOME_COV_DEPTH
            for (k = 1; k < ccnt; k++)
                fprintf(stdout, "%s\t%d\t%d\t%d\tPG:A:%c\n", dicts->s[read].name, chord[k-1].pos, chord[k].pos, 
                        chord[k-1].wgt, PG_DEPTH);
#endif

            prev = read + 1;
            if (i < naln) {
                read = alns[i].aread;
                j = i;
            }
        }
    }
    for (k = prev; k < nseq; k++) {
        hist[0] += dicts->s[k].len;
#ifdef DEBUG_GENOME_COV_HIST
        fprintf(stdout, "%s\t%d\t%d\t%d\t%f\tPG:A:%c\n", dicts->s[k].name, 
                0, dicts->s[k].len, dicts->s[k].len, 1., PG_HIST);
        covs_g.a[0] += dicts->s[k].len;
#endif
#ifdef DEBUG_GENOME_COV_DEPTH
        fprintf(stdout, "%s\t%d\t%d\t0\tPG:A:%c\n", dicts->s[k].name, 
                0, dicts->s[k].len, PG_DEPTH);
#endif
    }

#ifdef DEBUG_GENOME_COV_HIST
    seqtot = 0;
    for (k = 0, nseq = dicts->n; k < nseq; k++)
        seqtot += dicts->s[k].len;
    for (k = 0; k <= cmax_g; k++)
        if (covs_g.a[k] > 0)
            fprintf(stdout, "genome\t%lld\t%d\t%lld\t%f\tPG:A:%c\n", k, covs_g.a[k], seqtot, 
                    (double) covs_g.a[k]/seqtot, PG_HIST);
    
    kv_destroy(covs_s);
    kv_destroy(covs_g);
#endif

    if (print_hist) {
        seqtot = 0;
        for (k = 0, nseq = dicts->n; k < nseq; k++)
            seqtot += dicts->s[k].len;

        fprintf(stderr, "[M::%s] number of genome sequences: %lld\n", __func__, nseq);
        fprintf(stderr, "[M::%s] number of genome bases: %lld\n", __func__, seqtot);
        fprintf(stderr, "[M::%s] number of genome bases covered: %lld\n", __func__, seqtot - hist[0]);
        fprintf(stderr, "[M::%s] number of bases in alignments: %lld\n", __func__, alntot);
        fprintf(stderr, "[M::%s] coverage histogram:\n", __func__);
        for (k = 0; k < max_copy; k++)
            if (hist[k] > 0)
                fprintf(stderr, "[M::%s] %4lld: %12lld  %5.1f%%\n", __func__, k+1, hist[k], (100.*hist[k])/seqtot);
        if (hist[max_copy] > 0)
            fprintf(stderr, "[M::%s]  >%2d: %12lld  %5.1f%%\n", __func__, max_copy, hist[max_copy], (100.*hist[max_copy])/seqtot);
    }
    
    return hist;
}

void report_genome_coverage_histogram(aln_t *alns, int64 naln, sdict_t *dicts, int max_copy)
{
    int64 *hist;
    hist = genome_coverage_histogram(alns, naln, dicts, max_copy, 1);
    free(hist);
}

int srange_cmpfunc(const void *a, const void *b)
{
    int xm, ym;
    srange_t *x = (srange_t *) a;
    srange_t *y = (srange_t *) b;
    
    xm = x->seq;
    ym = y->seq;

    if (xm == ym) {
        xm = x->beg;
        ym = y->beg;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->end;
        ym = y->end;
    }
    return (xm > ym) - (xm < ym);
}

int64 pts_from_overlaps(void *data, int64 n, void *param, int *points)
{
    ovl_t *ovls = (ovl_t *) data;
    double min_qual = *(double *) param;
    int64 i, npts = 0;
    for (i = 0; i < n; i++) {
        if (ovls[i].del || ovls[i].qual < min_qual)
            continue;
        points[npts++] = -ovls[i].abpos;
        points[npts++] =  ovls[i].aepos;
    }
    return npts;
} 

int64 pts_from_ranges(void *data, int64 n, void *param, int *points)
{
    range_t *rngs = (range_t *) data;
    int64 i, npts = 0;
    for (i = 0; i < n; i++) {
        points[npts++] = -rngs[i].beg;
        points[npts++] =  rngs[i].end;
    }
    return npts;
}

int64 pts_from_sranges(void *data, int64 n, void *param, int *points)
{
    srange_t *rngs = (srange_t *) data;
    int64 i, npts = 0;
    for (i = 0; i < n; i++) {
        points[npts++] = -rngs[i].beg;
        points[npts++] =  rngs[i].end;
    }
    return npts;
}

int calc_coverage_from_intervals(cov_point_t *covs, int *points, int slen, 
    void *data, int64 nd, void *param, int64 (*pts_func) (void *, int64, void *, int *))
{
    cov_point_t *_covs = covs;
    int *_points;
    int64 npts;
    int pos, cov;

    npts = pts_func(data, nd, param, points);

    if (npts <= 0) {
        covs[0].pos = 0;
        covs[0].cov = 0;
        covs[1].pos = slen;
        covs[1].cov = 0;
        return 2;
    }

    qsort(points, npts, sizeof(int), absint_ord_acmpfunc);

    if (points[0] != 0) {
        covs->pos = 0;
        covs->cov = 0;
        covs++;
    }

    cov = 0;
    _points = points + npts;
    while (points < _points) {
        pos = abs(*points);
        cov += *(points++) > 0 ? -1 : 1;
        while (points < _points && abs(*points) == pos)
            cov += *(points++) > 0 ? -1 : 1;
        covs->pos = pos;
        covs->cov = cov;
        covs++;
    }

    if (abs(points[-1]) != slen) {
        covs->pos = slen;
        covs->cov = 0;
        covs++;
    }

    return covs - _covs;
}

double average_range_coverage(cov_point_t *covs, int ncov, int beg, int end)
{
    int i, b, e;
    double cov = 0.;
    
    if (beg + end < covs[ncov-1].pos) {
        // iterate from the left
        for (i = 1; i < ncov; i++) {
            b = covs[i-1].pos;
            e = covs[i].pos;
            if (e >= end) break;
            if (b <= beg) continue;
            cov += (MIN(e, end) - MAX(b, beg)) * covs[i-1].cov;
        }
    } else {
        // iterate from the right
        for (i = ncov - 1; i > 0; i--) {
            b = covs[i-1].pos;
            e = covs[i].pos;
            if (e <= beg) break;
            if (b >= end) continue;
            cov += (MIN(e, end) - MAX(b, beg)) * covs[i-1].cov;
        }
    }
    cov /= end - beg;

    return cov;
}

int total_range_coverage(cov_point_t *covs, int ncov, int beg, int end, int max_cov)
{
    int i, b, e;
    int cov = 0;
    
    if (beg + end < covs[ncov-1].pos) {
        // iterate from the left
        for (i = 1; i < ncov; i++) {
            b = covs[i-1].pos;
            e = covs[i].pos;
            if (e >= end) break;
            if (b <= beg) continue;
            if (covs[i-1].cov > max_cov) continue;
            cov += MIN(e, end) - MAX(b, beg);
        }
    } else {
        // iterate from the right
        for (i = ncov - 1; i > 0; i--) {
            b = covs[i-1].pos;
            e = covs[i].pos;
            if (e <= beg) break;
            if (b >= end) continue;
            if (covs[i-1].cov > max_cov) continue;
            cov += MIN(e, end) - MAX(b, beg);
        }
    }

    return cov;
}

static const double OVERLAP_MIN_LOWCOPY_SIZE = 4.7; // log10(50000)
static const double OVERLAP_MAX_LOWCOPY_SIZE = 6.0; // log10(1000000)
static const double OVERLAP_LOG_DIFF_LOWCOPY_SIZE = OVERLAP_MAX_LOWCOPY_SIZE - OVERLAP_MIN_LOWCOPY_SIZE;
static const double OVERLAP_MIN_LOWCOPY_RATIO_AT_MIN = 0.7;
static const double OVERLAP_MIN_LOWCOPY_RATIO_AT_MAX = 0.3;
static const double OVERLAP_LOWCOPY_RATIO_DIFF = OVERLAP_MIN_LOWCOPY_RATIO_AT_MAX - OVERLAP_MIN_LOWCOPY_RATIO_AT_MIN;
double overlap_min_lowcopy_ratio(double l)
{
    if (l <= 1.0) return OVERLAP_MIN_LOWCOPY_RATIO_AT_MIN;

    double log_L = log10(l);
    if (log_L <= OVERLAP_MIN_LOWCOPY_SIZE) return OVERLAP_MIN_LOWCOPY_RATIO_AT_MIN;
    if (log_L >= OVERLAP_MAX_LOWCOPY_SIZE) return OVERLAP_MIN_LOWCOPY_RATIO_AT_MAX;
    return OVERLAP_MIN_LOWCOPY_RATIO_AT_MIN + (log_L - OVERLAP_MIN_LOWCOPY_SIZE) / OVERLAP_LOG_DIFF_LOWCOPY_SIZE * OVERLAP_LOWCOPY_RATIO_DIFF;
}

int estimate_ploidy_number(ovl_t *ovls, int64 novl, sdict_t *dicts)
{
    if (ovls == NULL || novl <= 0)
        return 1;

    ovl_t *ovl;
    cov_point_t *covs;
    uint64 *index;
    uint32 a;
    int64 i, j, k, l, n, u, d, s, cut, movl, nseq, seqtot, copy, *hist;
    int ploidy, ncov, *cpts;
    double min_qual = OVL_MIN_QUAL;

    // sequence number
    nseq = dicts->n;

    // allocate memory
    MYCALLOC(index, nseq);
    if (index == NULL)
        mem_alloc_error("overlap index array");

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
    MYMALLOC(covs, movl*2 + 2);
    MYCALLOC(hist, MAX_PLOIDY_NUMBER + 1);
    if (cpts == NULL || covs == NULL || hist == NULL)
        mem_alloc_error("coverage points");
    // calculate sequence coverage
    for (i = 0; i < nseq; i++) {
        ovl = ovls + (index[i] >> 32);
        n = (uint32) index[i];
        ncov = calc_coverage_from_intervals(covs, cpts, dicts->s[i].len, ovl, n, &min_qual, pts_from_overlaps);
        for (j = 1; j < ncov; j++) {
            copy = covs[j-1].cov;
            if (copy > MAX_PLOIDY_NUMBER)
                copy = MAX_PLOIDY_NUMBER;
            hist[copy] += covs[j].pos - covs[j-1].pos;
        }
    }
    free(cpts);
    free(covs);
    free(index);
    
    fprintf(stderr, "[M::%s] estimate ploidy number from overlaps\n", __func__);
    
    seqtot = 0;
    for (k = 0, nseq = dicts->n; k < nseq; k++)
        seqtot += dicts->s[k].len;
    fprintf(stderr, "[M::%s] number of genome sequences: %lld\n", __func__, nseq);
    fprintf(stderr, "[M::%s] number of genome bases: %lld\n", __func__, seqtot);
    fprintf(stderr, "[M::%s] number of overlapped bases: %lld\n", __func__, seqtot - hist[0]);
    fprintf(stderr, "[M::%s] coverage histogram:\n", __func__);
    for (k = 0; k < MAX_PLOIDY_NUMBER; k++)
        if (hist[k] > 0)
            fprintf(stderr, "[M::%s] %4lld: %12lld  %5.1f%%\n", __func__, k+1, hist[k], (100.*hist[k])/seqtot);
    if (hist[MAX_PLOIDY_NUMBER] > 0)
        fprintf(stderr, "[M::%s]  >%2d: %12lld  %5.1f%%\n", __func__, MAX_PLOIDY_NUMBER, hist[MAX_PLOIDY_NUMBER], (100.*hist[MAX_PLOIDY_NUMBER])/seqtot);
    
    // find cumulative sum
    for (k = 1; k <= MAX_PLOIDY_NUMBER; k++)
        hist[k] += hist[k-1];
    
    // find minimum ploidy number
    // i.e., coverage at which the cumulative sum reaches 50% of the genome size
    cut = hist[MAX_PLOIDY_NUMBER] * .5;
    for (l = 0; l <= MAX_PLOIDY_NUMBER; l++)
        if (hist[l] >= cut)
            break;
    // find maximum ploidy number
    // i.e., coverage at which the cumulative sum reaches 90% of the genome size
    cut = hist[MAX_PLOIDY_NUMBER] * .9;
    for (u = l; u <= MAX_PLOIDY_NUMBER; u++)
        if (hist[u] >= cut)
            break;
    if (u < MAX_PLOIDY_NUMBER) u += 1; // make it exclusive
    // find the elbow point in (l, u)
    // i.e., maximum distance from line connecting (0, 0) and (u+1, hist[u])
    d = INT64_MIN;
    ploidy = l;
    for (k = l; k < u; k++) {
        s = (u + 1) * hist[k] - (k + 1) * hist[u];
        if (s > d) {
            d = s;
            ploidy = k;
        }
    }
    ploidy += 1; // convert alignment copy number to ploidy

    free(hist);
    
    fprintf(stderr, "[M::%s] estimated ploidy number: %d\n", __func__, ploidy);

    return ploidy;
}

void genome_coverage_summary(rangetree_t *rt, sdict_t *dicts, int max_copy)
{
    int i, n;
    int64 nc, nb, tb, *ns;
    n = dicts->n;
    nc = 0;
    nb = 0;
    MYCALLOC(ns, max_copy+3);
    for (i = 0; i < n; i++) {
        nc += rangetree_coverage(rt+i, 0, INT64_MAX, &ns[max_copy+2], ns, max_copy);
        nb += ns[max_copy+2];
    }
    tb = 0;
    for (i = 0; i < n; i++)
        tb += dicts->s[i].len;
    ns[0] = tb;
    for (i = 1; i <= max_copy+1; i++)
        ns[0] -= ns[i];
    fprintf(stderr, "[M::%s] number of genome sequences: %d\n", __func__, n);
    fprintf(stderr, "[M::%s] number of genome bases: %lld\n", __func__, tb);
    fprintf(stderr, "[M::%s] number of genome bases covered: %lld\n", __func__, nc);
    fprintf(stderr, "[M::%s] number of bases in alignments: %lld\n", __func__, nb);
    fprintf(stderr, "[M::%s] genome coverage summary:\n", __func__);

    for (i = 0; i <= max_copy; i++)
        if (ns[i] > 0)
            fprintf(stderr, "[M::%s] %4d: %12lld  %5.1f%%\n", __func__, i+1, ns[i], (100.*ns[i])/tb);
    if (ns[max_copy+1] > 0)
        fprintf(stderr, "[M::%s]  >%2d: %12lld  %5.1f%%\n", __func__, max_copy+1, ns[max_copy+1], (100.*ns[max_copy+1]/tb));
    free(ns);
}
