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
#include <math.h>
#include <float.h>

#include "kvec.h"

#include "ec.h"
#include "misc.h"

#undef DEBUG_ERROR_CORRECTION

ec_conf_t ec_conf = {
    .bin_size = 1000,
    .min_frag = 10000,
    .med_drop = .3,
    .rec_rate = .8,
    .p_thresh = .01,
};

/* Interface for DynamicPercentile */
typedef struct DynamicPercentile DynamicPercentile_t;
static DynamicPercentile_t *dynamic_percentile_create(size_t max_elements, double percentile);
static void dynamic_percentile_destroy(DynamicPercentile_t *dm);
static void dynamic_percentile_clear(DynamicPercentile_t *dm);
static void dynamic_percentile_insert(DynamicPercentile_t *dm, double val);
static double dynamic_percentile_get(const DynamicPercentile_t *dm);

/* Interface for static array median using quickselect */
static double quick_array_median(double *arr, size_t n);

typedef enum { STATE_BASELINE, STATE_INPEAK } FState;

static void remove_spikes(double *span, int n, double df, double mf, int plat_max)
{
    if (n < 2) return;

    FState state;
    double f, sf, af, base_val;
    int i, j, peak_beg, plat_len;
    
    state = STATE_BASELINE;
    peak_beg = -1;
    plat_len =  0;
    base_val =  0.;
    for (i = 0; i < n; i++) {
        sf = span[i+1] - span[i] - df;
        af = fabs(sf);

        if (state == STATE_BASELINE) {
            if (sf > mf) {
                state = STATE_INPEAK;
                base_val =  span[i];
                peak_beg =  i;
                plat_len =  0;
            }
        } else if (state == STATE_INPEAK) {
            if (af <= mf)
                plat_len++;
            else
                plat_len = 0;

            if (af <= mf && span[i] <= base_val + mf * (i - peak_beg)) {
                // end of peak, return to baseline
                // interpolate the peak region
                f = (span[i] - base_val) / (i - peak_beg);
                for (j = peak_beg + 1; j < i; j++)
                    span[j] = span[j-1] + f;
                state = STATE_BASELINE;
            } else if (plat_len > plat_max) state = STATE_BASELINE;
        }
    }
}

// Computes the lower-tail p-value P(X <= k) for a Poisson(lambda) distribution
static double calculate_poisson_lower_pvalue(int k, double lambda)
{
    if (k < 0)          return 0.0;
    if (lambda <= 0.0)  return 1.0;

    const double sqrt_lambda = sqrt(lambda);
    const double dk = (double) k;

    if (dk > lambda + 10.0 * sqrt_lambda) {
        double log_term = (dk + 1.0) * log(lambda) - lambda - lgamma(dk + 2.0);
        if (log_term < -708.0) return 1.0;

        double term = exp(log_term);
        double sum  = term;
        int x;
        for (x = k + 2; term > DBL_EPSILON * sum; x++) {
            term *= lambda / (double)x;
            sum  += term;
        }
        return sum >= 1.0 ? 0.0 : 1.0 - sum;
    }

    if (dk < lambda - 10.0 * sqrt_lambda) {
        double term = exp(-lambda);
        if (term == 0.0) return 0.0;
        double cdf = term;
        int x;
        for (x = 1; x <= k; x++) {
            term *= lambda / (double)x;
            cdf  += term;
        }
        return cdf > 1.0 ? 1.0 : cdf;
    }

    double log_pk = dk * log(lambda) - lambda - lgamma(dk + 1.0);
    double log_px, sum = 0.0;
    int x;
    for (x = 0; x <= k; x++) {
        log_px = (double)x * log(lambda) - lambda - lgamma((double)x + 1.0);
        sum += exp(log_px - log_pk);
    }
    double log_cdf = log(sum) + log_pk;
    if (log_cdf >= 0.0) return 1.0;
    double cdf = exp(log_cdf);

    return cdf < 0.0 ? 0.0 : (cdf > 1.0 ? 1.0 : cdf);
}

typedef kvec_t(ec_pos_t) ec_pos_v;

static int ec_pval_cmpfunc(const void *a, const void *b)
{
    double pa = ((const ec_pos_t *)a)->pval;
    double pb = ((const ec_pos_t *)b)->pval;
    return (pa > pb) - (pa < pb);
}

static int ec_pos_cmpfunc(const void *a, const void *b)
{
    const ec_pos_t *pa = (const ec_pos_t *)a;
    const ec_pos_t *pb = (const ec_pos_t *)b;
    if (pa->seq != pb->seq)
        return (pa->seq > pb->seq) - (pa->seq < pb->seq);
    return (pa->pos > pb->pos) - (pa->pos < pb->pos);
}

static int call_breaks(double *cnts, double *meds, int n, int min_bins, 
    double med_drop, double rec_rate, int64 *tspace)
{
    int i, j, p, beg, end, ncall, run_length;
    double min_med, min_cnt;

    MYBZERO(tspace, n+1);
    run_length = 0;
    for (i = 0; i <= n ; i++) {
        if (i < n && meds[i] >= rec_rate) {
            run_length++;
        } else {
            if (run_length >= min_bins)
                for (j = i - run_length; j < i; j++)
                    tspace[j] = 1;
            run_length = 0;
        }
    }

    ncall = 0;
    i = 0;
    while (i < n) {
        if (tspace[i] == 0) {
            beg = i;
            while (i < n && tspace[i] == 0) ++i;
            end = i;
            if (beg > 0 && end < n) {
                min_med = meds[beg];
                min_cnt = cnts[beg];
                p = beg;
                for (j = beg + 1; j < end; j++) {
                    if (meds[j] < min_med) {
                        min_med = meds[j];
                        min_cnt = cnts[j];
                        p = j;
                    } else if (meds[j] == min_med && 
                        cnts[j] < min_cnt) {
                        min_cnt = cnts[j];
                        p = j;
                    }
                }
                if (min_med < med_drop)
                    // mark the position p as a call break
                    tspace[ncall++] = p; // this is safe
            }
        } else i++;
    }

    return ncall;
}

ec_pos_t *ec_call_breaks(hic_t *hics, int64 nhic, sdict_t *dicts, int *_ncall)
{
    ec_pos_v _calls = {0, 0, 0}, *calls = &_calls;
    ec_pos_t *call;
    DynamicPercentile_t *dp50;
    int64 i, j, k, p, n, beg, end, nseq, acnt, tcnt;
    int64 mbin, tbin, *tspace, *bins;
    uint32 s;
    double *sspace, *span, **spans, *cnt, *med, *meds[2];
    double f, df, mf, med_drop, rec_rate, p_thresh;
    int a, b, min_bins, max_bins, ncall;
    
    *_ncall = 0;
    if (hics == NULL || nhic <= 0)
        return NULL;
    if (ec_conf.bin_size <= 0 || ec_conf.min_frag < ec_conf.bin_size || 
        ec_conf.med_drop <= 0. || ec_conf.med_drop > 1. || 
        ec_conf.rec_rate <= 0. || ec_conf.rec_rate > 1. ||
        ec_conf.p_thresh <= 0. || ec_conf.p_thresh > 1. )
        return NULL;

    // number of sequences
    nseq = dicts->n;

    // maximum sequence length in windows
    MYCALLOC(bins, nseq);
    if (bins == NULL)
        mem_alloc_error("error correction bins array");
    mbin = tbin = 0;
    for (i = 0; i < nseq; i++) {
        bins[i] = (dicts->s[i].len - 1) / ec_conf.bin_size + 1;
        tbin += bins[i];
        if (bins[i] > mbin) mbin = bins[i];
    }
    // allocate memory
    MYCALLOC(tspace, mbin + 1);
    MYCALLOC(sspace, tbin * 2 + mbin * 2 + nseq + 1);
    MYMALLOC(spans, nseq);
    if (tspace == NULL || sspace == NULL || spans == NULL)
        mem_alloc_error("error correction memory array");
    meds[0] = sspace + tbin + 1;
    meds[1] = meds[0] + mbin;
    spans[0] = meds[1] + mbin;
    for (i = 1; i < nseq; i++)
        spans[i] = spans[i-1] + bins[i-1] + 1;
    // estimate max_distance for intra-hic links
    for (i = 0; i < nhic; i++) {
        if (hics[i].aseq != hics[i].bseq)
            continue;
        tspace[hics[i].bpos-hics[i].apos] += hics[i].nhic;
    }
    acnt = tcnt = 0;
    for (i = 0; i < mbin; i++)
        tcnt += tspace[i];
    fprintf(stderr, "[M::%s] total intra hic links: %lld\n", __func__, tcnt);
    max_bins = 0;
    f = (double)((int) ((double) (tspace[0] - 1) / tcnt * 10)) / 10 + .1;
    for (i = 0; i < mbin; i++) {
        acnt += tspace[i];
        if (acnt >= f * tcnt) {
            f += .1;
            if (f < 1.) max_bins = i; // using 90% intra links to estimate max_bins dist
        }
    }
    max_bins = MAX(max_bins, 1);
    if (max_bins < ec_conf.min_frag / ec_conf.bin_size) // minimum distance
        max_bins = ec_conf.min_frag / ec_conf.bin_size;
    fprintf(stderr, "[M::%s] using hic window size: %d bp\n", __func__, ec_conf.bin_size);
    fprintf(stderr, "[M::%s] using maximum distance: %d kb\n", __func__, max_bins * ec_conf.bin_size / 1000);

    // collect spans
    for (i = 0; i < nhic; i++) {
        if ((s = hics[i].aseq) != hics[i].bseq)
            continue;
        a = hics[i].apos;
        b = hics[i].bpos;
        if (a > b) SWAP(int, a, b);
        if (b - a > 2 && b - a <= max_bins) {
            f = hics[i].nhic;
            spans[s][a+1] += f;
            spans[s][b-1] -= f;
        }
    }

    // finalise spans
    for (i = 0; i < nseq; i++) {
        n = bins[i];
        span = spans[i];
        for (j = 0; j < n; j++)
            span[j+1] += span[j];
    }

    // minimum number of bins for a valid fragment
    min_bins = (ec_conf.min_frag - 1) / ec_conf.bin_size + 1, min_bins = MAX(1, min_bins);

    // remove spikes around abnormal discrete first difference
    tcnt = 0;
    for (i = 0; i < nseq; i++) {
        n = bins[i];
        span = spans[i];
        for (j = 0; j < n; j++)
            sspace[tcnt++] = span[j+1] - span[j];
    }
    df = quick_array_median(sspace, tcnt);
    for (i = 0; i < tcnt; i++)
        sspace[i] = fabs(sspace[i] - df);
    mf = quick_array_median(sspace, tcnt);
    //fprintf(stderr, "[M::%s] median of first difference: %.4f\n", __func__, df);
    //fprintf(stderr, "[M::%s] MAD-based modified z-score: %.4f\n", __func__, mf);

    // apply spike filter
    mf = mf * 3.5 / 0.6745; // outlier threshold 0.6745 * |d - df| > 3.5 * mf
    for (i = 0; i < nseq; i++) {
        n = bins[i];
        span = spans[i];
        remove_spikes(span, n, df, mf, min_bins);
    }

    // call breakpoints
    med_drop = ec_conf.med_drop;
    rec_rate = ec_conf.rec_rate;
    dp50 = dynamic_percentile_create(mbin, 0.5);
    for (i = 0; i < nseq; i++) {
        n = bins[i];
        span = spans[i];

        if (n <= min_bins * 2)
            continue;

        // calculate left- and right-medians
        for (k = 1; k >= -1; k -= 2) {
            if (k > 0) {
                beg = 0;
                end = n;
                med = meds[0];
            } else {
                beg = n - 1;
                end = -1;
                med = meds[1];
            }
            // reset dynamic percentiles for the new sequence
            dynamic_percentile_clear(dp50);
            for (j = beg; j != end; j += k) {
                dynamic_percentile_insert(dp50, span[j]);
                med[j] = dynamic_percentile_get(dp50);
            }
        }

        // collect min-medians
        sspace[0]   = meds[1][1];
        sspace[n-1] = meds[0][n-2];
        for (j = 1; j < n-1; j++)
            sspace[j] = MIN(meds[0][j-1], meds[1][j+1]);
        for (j = 0; j < n; j++) {
            f = sspace[j];
            f = f > 0? span[j] / f : 1.;
            sspace[j] = MIN(f, 1.);
        }

        // do actual detection
        ncall = call_breaks(span, sspace, n, min_bins, med_drop, rec_rate, tspace);
        for (j = 0; j < ncall; j++) {
            kv_pushp(ec_pos_t, *calls, &call);
            p = tspace[j];
            cnt = call->cnts;
            cnt[0] = meds[0][p-1];
            cnt[1] = span[p];
            cnt[2] = meds[1][p+1];
            call->seq = i;
            call->pos = p;
            call->pval = calculate_poisson_lower_pvalue((int)cnt[1], MIN(cnt[0], cnt[2]));
        }

#ifdef DEBUG_ERROR_CORRECTION
        MYBZERO(tspace, n);
        for (j = calls->n - ncall; j < calls->n; j++)
            tspace[calls->a[j].pos] = (int) (-log(calls->a[j].pval + 1e-300));
        for (j = 0; j < n; j++)
            fprintf(stderr, "[M::%s] M/M %s\t%8lld\t%8lld\t%.4f\t%.4f\t%lld\n", __func__, 
                dicts->s[i].name, j, j * ec_conf.bin_size, span[j], sspace[j], tspace[j]);
#endif
    }
    fprintf(stderr, "[M::%s] number of errors called: %zu\n", __func__, calls->n);

    if (calls->n) {
        // Benjamini-Hochberg adjustment
        call  = calls->a;
        ncall = calls->n;
        qsort(call, ncall, sizeof(ec_pos_t), ec_pval_cmpfunc);
        p_thresh = ec_conf.p_thresh / ncall;
        for (i = ncall - 1; i >= 0; i--)
            if (call[i].pval <= (double) (i + 1) * p_thresh)
                break;
        ncall = i + 1;
        if (ncall) {
            p_thresh = call[i].pval;
            // sort back by sequence and position
            qsort(call, ncall, sizeof(ec_pos_t), ec_pos_cmpfunc);
        }
        // some summary statistics
        fprintf(stderr, "[M::%s] maximum p-value: %.3e\n", __func__, ec_conf.p_thresh);
        fprintf(stderr, "[M::%s] BH adjusted p-value: %.3e\n", __func__, p_thresh);
        fprintf(stderr, "[M::%s] final errors called: %d\n", __func__, ncall);
        // set the real number of calls
        calls->n = ncall;
    }

    dynamic_percentile_destroy(dp50);
    free(tspace);
    free(sspace);
    free(spans);
    free(bins);

    if (calls->n == 0) {
        free(calls->a);
        return NULL;
    }

    // convert to the actual sequence position
    for (i = 0; i < calls->n; i++)
        calls->a[i].pos *= ec_conf.bin_size;

    *_ncall = calls->n;
    return calls->a;
}

void ec_write_report(ec_pos_t *calls, int ncall, sdict_t *dicts, FILE *fo)
{
    int64 i;
    fprintf(fo, "#sequence\tposition\tmed-left\tmedian\tmed-right\tpval\n");
    for (i = 0; i < ncall; i++) {
        fprintf(fo, "%s\t%u\t%.0f\t%.0f\t%.0f\t%.4e\n",
            dicts->s[calls[i].seq].name, calls[i].pos,
            calls[i].cnts[0], calls[i].cnts[1], calls[i].cnts[2], calls[i].pval);
    }
}

void ec_write_agp(ec_pos_t *calls, int ncall, sdict_t *dicts, FILE *fo)
{
    int i;
    uint32 s;
    for (s = 0, i = 0; s < dicts->n; s++) {
        uint32 beg = 0, end, part = 0, piece = 0;
        while (i < ncall && calls[i].seq < s)
            i++;
        while (i < ncall && calls[i].seq == s) {
            end = calls[i].pos;
            if (end > beg) {
                fprintf(fo, "%s_ec_%u\t1\t%u\t1\tW\t%s\t%u\t%u\t+\n",
                    dicts->s[s].name, ++piece, end - beg, dicts->s[s].name, beg + 1, end);
                beg = end;
            }
            i++;
        }
        if (beg < dicts->s[s].len) {
            part++;
            fprintf(fo, "%s_ec_%u\t1\t%u\t%u\tW\t%s\t%u\t%u\t+\n",
                dicts->s[s].name, ++piece, dicts->s[s].len - beg, part, dicts->s[s].name,
                beg + 1, dicts->s[s].len);
        }
    }
}

/* *************************************************************************** */
/* *********************** Heap based Dynamic Percentile ********************* */
/* *************************************************************************** */

typedef struct {
    double *data;
    size_t size;
    size_t capacity;
    int is_max_heap;
} Heap;

struct DynamicPercentile {
    Heap max_heap;       // lower portion of numbers
    Heap min_heap;       // upper portion of numbers
    size_t max_elements;
    double percentile;   // 0.0-1.0, e.g. 0.5 = median, 0.8 = 80th percentile
};

static inline void swap_double(double *a, double *b) 
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

static void heap_push(Heap *h, double val) 
{
    if (h->size >= h->capacity) return;

    h->data[h->size] = val;
    size_t parent, idx = h->size++;

    while (idx > 0) {
        parent = (idx - 1) / 2;
        if (h->is_max_heap ? (h->data[idx] > h->data[parent]) : (h->data[idx] < h->data[parent])) {
            swap_double(&h->data[idx], &h->data[parent]);
            idx = parent;
        } else break;
    }
}

static double heap_pop(Heap *h) 
{
    if (h->size == 0) return 0.0;

    double root = h->data[0];
    size_t lft, rht, best, idx = 0;

    h->data[0] = h->data[--h->size];
    while (1) {
        lft = 2 * idx + 1;
        rht = 2 * idx + 2;
        best = idx;

        if (lft < h->size &&
            (h->is_max_heap ? (h->data[lft] > h->data[best]) : (h->data[lft] < h->data[best])))
            best = lft;
        if (rht < h->size &&
            (h->is_max_heap ? (h->data[rht] > h->data[best]) : (h->data[rht] < h->data[best])))
            best = rht;

        if (best != idx) {
            swap_double(&h->data[idx], &h->data[best]);
            idx = best;
        } else break;
    }

    return root;
}

static DynamicPercentile_t *dynamic_percentile_create(size_t max_elements, double percentile) 
{
    if (max_elements == 0) return NULL;
    if (percentile < 0.0 || percentile > 1.0) return NULL;

    DynamicPercentile_t *dp = malloc(sizeof(DynamicPercentile_t));
    if (!dp) return NULL;

    dp->max_elements = max_elements;
    dp->percentile = percentile;

    // max_heap (lower portion) can hold up to ~percentile fraction of elements,
    // min_heap (upper portion) the remaining ~(1-percentile) fraction.
    // +2 gives slack for the transient state right after a push, before rebalance.
    size_t max_cap = (size_t)ceil(percentile * (double)max_elements) + 2;
    size_t min_cap = (size_t)ceil((1.0 - percentile) * (double)max_elements) + 2;
    if (max_cap > max_elements + 2) max_cap = max_elements + 2;
    if (min_cap > max_elements + 2) min_cap = max_elements + 2;

    dp->max_heap.capacity = max_cap;
    dp->max_heap.size = 0;
    dp->max_heap.is_max_heap = 1;
    dp->max_heap.data = malloc(max_cap * sizeof(double));

    dp->min_heap.capacity = min_cap;
    dp->min_heap.size = 0;
    dp->min_heap.is_max_heap = 0;
    dp->min_heap.data = malloc(min_cap * sizeof(double));

    if (!dp->max_heap.data || !dp->min_heap.data) {
        dynamic_percentile_destroy(dp);
        return NULL;
    }

    return dp;
}

static void dynamic_percentile_destroy(DynamicPercentile_t *dp) 
{
    if (!dp) return;
    if (dp->max_heap.data) free(dp->max_heap.data);
    if (dp->min_heap.data) free(dp->min_heap.data);
    free(dp);
}

static void dynamic_percentile_clear(DynamicPercentile_t *dp) 
{
    if (!dp) return;
    dp->max_heap.size = 0;
    dp->min_heap.size = 0;
}

// Number of elements max_heap should hold once there are n total elements,
// so that max_heap.data[0] / min_heap.data[0] straddle the desired rank.
static size_t desired_max_heap_size(double percentile, size_t n) 
{
    if (n == 0) return 0;
    double real_rank = percentile * (double)(n - 1); // 0-indexed rank, "linear" method
    size_t floor_idx = (size_t)floor(real_rank);
    return floor_idx + 1;
}

static void dynamic_percentile_insert(DynamicPercentile_t *dp, double val) 
{
    if (!dp) return;

    // 1. Route insertion based on current max heap root (still correct for any
    //    split ratio - it's just an initial guess; the rebalance below enforces
    //    the actual invariant).
    if (dp->max_heap.size == 0 || val <= dp->max_heap.data[0]) {
        heap_push(&dp->max_heap, val);
    } else {
        heap_push(&dp->min_heap, val);
    }

    // 2. Rebalance to the desired split for the chosen percentile.
    //    At most one element ever needs to move.
    size_t n = dp->max_heap.size + dp->min_heap.size;
    size_t target = desired_max_heap_size(dp->percentile, n);

    if (dp->max_heap.size > target) {
        double v = heap_pop(&dp->max_heap);
        heap_push(&dp->min_heap, v);
    } else if (dp->max_heap.size < target) {
        double v = heap_pop(&dp->min_heap);
        heap_push(&dp->max_heap, v);
    }
}

static double dynamic_percentile_get(const DynamicPercentile_t *dp) 
{
    if (!dp) return NAN;
    size_t n = dp->max_heap.size + dp->min_heap.size;
    if (n == 0) return NAN;

    double real_rank = dp->percentile * (double)(n - 1);
    size_t floor_idx = (size_t)floor(real_rank);
    double frac = real_rank - (double)floor_idx;

    double lower = dp->max_heap.data[0]; // value at rank floor_idx

    if (frac == 0.0 || dp->min_heap.size == 0) {
        return lower;
    }

    double upper = dp->min_heap.data[0]; // value at rank floor_idx + 1
    return lower + frac * (upper - lower);
}

// In-place QuickSelect to find the k-th smallest element (0-indexed)
static double quick_select(double *arr, size_t left, size_t right, size_t k) 
{
    size_t i, j;
    double pivot;

    while (left < right) {
        pivot = arr[right];
        i = left;
        
        for (j = left; j < right; j++) {
            if (arr[j] <= pivot) {
                swap_double(&arr[i], &arr[j]);
                i++;
            }
        }
        swap_double(&arr[i], &arr[right]);

        if (i == k) return arr[i];
        if (i < k) left = i + 1;
        else right = i - 1;
    }

    return arr[left];
}

static double quick_array_median(double *arr, size_t n) 
{
    if (n == 0) return 0.0;
    
    if (n % 2 != 0) {
        return quick_select(arr, 0, n - 1, n / 2);
    } else {
        // For even N, average the two middle values
        double right_mid = quick_select(arr, 0, n - 1, n / 2);
        double left_mid  = quick_select(arr, 0, (n / 2) - 1, (n / 2) - 1);
        return (left_mid + right_mid) / 2.0;
    }
}
