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

 #include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "kvec.h"

#include "ec.h"
#include "misc.h"

#define DEBUG_ERROR_CORRECTION

typedef struct {
    int n;
    double *endpoints;
    double *bias;
    double *diff;
    double *span;
    double *sum;
} ec_track_t;

typedef kvec_t(ec_candidate_t) ec_candidate_vec_t;

static int ec_candidate_p_cmpfunc(const void *a, const void *b)
{
    const ec_candidate_t *x = (const ec_candidate_t *) a;
    const ec_candidate_t *y = (const ec_candidate_t *) b;
    if (x->p != y->p)
        return (x->p > y->p) - (x->p < y->p);
    if (x->seq != y->seq)
        return (x->seq > y->seq) - (x->seq < y->seq);
    return (x->bin > y->bin) - (x->bin < y->bin);
}

static int ec_candidate_pos_cmpfunc(const void *a, const void *b)
{
    const ec_candidate_t *x = (const ec_candidate_t *) a;
    const ec_candidate_t *y = (const ec_candidate_t *) b;
    if (x->seq != y->seq)
        return (x->seq > y->seq) - (x->seq < y->seq);
    return (x->bin > y->bin) - (x->bin < y->bin);
}

void ec_conf_init(ec_conf_t *conf)
{
    conf->bin_size = 1000;
    conf->max_distance = 1000000;
    conf->control_window = 30000;
    conf->min_fragment = 30000;
    conf->min_contacts = 10;
    conf->max_ratio = .35;
    conf->fdr = .01;
}

static void ec_tracks_destroy(ec_track_t *tracks, uint32 n)
{
    uint32 i;
    for (i = 0; i < n; i++) {
        free(tracks[i].endpoints);
        free(tracks[i].bias);
        free(tracks[i].diff);
        free(tracks[i].span);
        free(tracks[i].sum);
    }
    free(tracks);
}

static double ec_sum(double *sum, int beg, int end)
{
    if (end <= beg)
        return 0.;
    return sum[end] - sum[beg];
}

static void ec_normalize_endpoints(ec_track_t *track, int radius)
{
    int i, beg, end;
    double local;

    track->sum[0] = 0.;
    for (i = 0; i < track->n; i++)
        track->sum[i+1] = track->sum[i] + track->endpoints[i];

    for (i = 0; i < track->n; i++) {
        beg = MAX(0, i - radius);
        end = MIN(track->n, i + radius + 1);
        local = ec_sum(track->sum, beg, end) / (end - beg);
        track->bias[i] = MIN_MAX((track->endpoints[i] + 1.) / (local + 1.), .1, 10.);
    }
}

static void ec_adjust_qvalues(ec_candidate_vec_t *cands)
{
    int64 i, n;
    double q;

    n = cands->n;
    if (!n)
        return;
    qsort(cands->a, n, sizeof(ec_candidate_t), ec_candidate_p_cmpfunc);
    q = 1.;
    for (i = n; i > 0; i--) {
        double x = cands->a[i-1].p * n / i;
        q = MIN(q, x);
        cands->a[i-1].q = MIN(q, 1.);
    }
}

int ec_call_breaks(hic_t *hics, int64 nhic, sdict_t *dicts, ec_conf_t *conf, ec_candidates_t *cands)
{
    ec_track_t *tracks;
    ec_candidate_vec_t calls;
    int64 i;
    int a, b, d, radius, min_bins, control_bins, gap_bins;
    uint32 s;

    cands->a = NULL;
    cands->n = 0;
    if (hics == NULL || nhic <= 0)
        return 0;
    if (conf->bin_size <= 0 || conf->max_distance < conf->bin_size ||
        conf->control_window < conf->bin_size || conf->min_fragment < conf->bin_size)
        return 1;

    MYCALLOC(tracks, dicts->n);
    if (tracks == NULL)
        mem_alloc_error("error correction tracks");
    for (s = 0; s < dicts->n; s++) {
        tracks[s].n = (dicts->s[s].len - 1) / conf->bin_size + 1;
        MYCALLOC(tracks[s].endpoints, tracks[s].n);
        MYCALLOC(tracks[s].bias, tracks[s].n);
        MYCALLOC(tracks[s].diff, tracks[s].n + 1);
        MYCALLOC(tracks[s].span, tracks[s].n);
        MYCALLOC(tracks[s].sum, tracks[s].n + 1);
        if (tracks[s].endpoints == NULL || tracks[s].bias == NULL || tracks[s].diff == NULL ||
            tracks[s].span == NULL || tracks[s].sum == NULL)
            mem_alloc_error("error correction track array");
    }

    for (i = 0; i < nhic; i++) {
        if (hics[i].aseq != hics[i].bseq)
            continue;
        a = hics[i].apos;
        b = hics[i].bpos;
        if (a == b || a < 0 || b < 0 || a >= tracks[hics[i].aseq].n || b >= tracks[hics[i].aseq].n)
            continue;
        tracks[hics[i].aseq].endpoints[a] += hics[i].nhic;
        tracks[hics[i].aseq].endpoints[b] += hics[i].nhic;
    }

    radius = MAX(1, conf->control_window / conf->bin_size / 2);
    for (s = 0; s < dicts->n; s++)
        ec_normalize_endpoints(tracks + s, radius);

    for (i = 0; i < nhic; i++) {
        double weight;
        ec_track_t *track;
        if (hics[i].aseq != hics[i].bseq)
            continue;
        s = hics[i].aseq;
        a = hics[i].apos;
        b = hics[i].bpos;
        if (a > b) SWAP(int, a, b);
        d = b - a;
        track = tracks + s;
        if (a == b || d > conf->max_distance / conf->bin_size || a < 0 || b >= track->n)
            continue;
        weight = hics[i].nhic / sqrt(track->bias[a] * track->bias[b]);
        track->diff[a+1] += weight;
        track->diff[b+1] -= weight;
    }

    kv_init(calls);
    min_bins = (conf->min_fragment + conf->bin_size - 1) / conf->bin_size;
    control_bins = MAX(2, conf->control_window / conf->bin_size);
    gap_bins = MAX(1, control_bins / 4);
#ifdef DEBUG_ERROR_CORRECTION
    fprintf(stderr, "[M::%s] EC observed/expected per bin along each sequence:\n", __func__);
#endif
    for (s = 0; s < dicts->n; s++) {
        ec_track_t *track = tracks + s;
        int p;
        double value;
        track->sum[0] = 0.;
        value = 0.;
        for (p = 0; p < track->n; p++) {
            value += track->diff[p];
            track->span[p] = value;
            track->sum[p+1] = track->sum[p] + value;
        }
        for (p = min_bins; p <= track->n - min_bins; p++) {
            int lb = MAX(0, p - gap_bins - control_bins);
            int le = p - gap_bins;
            int rb = p + gap_bins;
            int re = MIN(track->n, p + gap_bins + control_bins);
            double observed, left, right, expected, z, ratio, pval;
            if (le <= lb || re <= rb)
                continue;
            observed = track->span[p-1];
            left = ec_sum(track->sum, lb, le) / (le - lb);
            right = ec_sum(track->sum, rb, re) / (re - rb);
            expected = (left + right) / 2.;
#ifdef DEBUG_ERROR_CORRECTION
            fprintf(stderr, "[M::%s] %s\t%8u\t%8u\t%.4f\t%.4f\n",
                __func__, dicts->s[s].name, p, (uint32) p * conf->bin_size, observed, expected);
#endif
            if (expected < conf->min_contacts || observed >= expected * conf->max_ratio)
                continue;
            z = (expected - observed) / sqrt(expected + 1.);
            pval = .5 * erfc(z / M_SQRT2);
            ratio = observed / (expected + 1e-9);
            kv_push(ec_candidate_t, calls, ((ec_candidate_t) {s, p, (uint32) p * conf->bin_size,
                observed, expected, ratio, z, pval, 1., 0}));
        }
    }

    ec_adjust_qvalues(&calls);
    qsort(calls.a, calls.n, sizeof(ec_candidate_t), ec_candidate_pos_cmpfunc);

#ifdef DEBUG_ERROR_CORRECTION
    fprintf(stderr, "[M::%s] EC candidates before filtering:\n", __func__);
    for (i = 0; i < calls.n; i++) {
        ec_candidate_t *call = calls.a + i;
        fprintf(stderr, "[M::%s] %s\t%8u\t%8u\t%.4f\t%.4f\t%.4f\t%.4f\t%.3e\t%.3e\n",
            __func__, dicts->s[call->seq].name, call->pos, call->bin, call->observed,
            call->expected, call->ratio, call->z, call->p, call->q);
    }
#endif

    for (i = 0; i < calls.n; i++) {
        ec_candidate_t *call = calls.a + i;
        int64 j;
        if (call->q > conf->fdr)
            continue;
        for (j = i - 1; j >= 0 && calls.a[j].seq == call->seq; j--) {
            if (calls.a[j].accepted && call->bin - calls.a[j].bin < min_bins)
                break;
        }
        if (j >= 0 && calls.a[j].seq == call->seq && calls.a[j].accepted &&
            call->bin - calls.a[j].bin < min_bins) {
            if (call->z > calls.a[j].z) {
                calls.a[j].accepted = 0;
                call->accepted = 1;
            }
        } else {
            call->accepted = 1;
        }
    }

    cands->a = calls.a;
    cands->n = calls.n;
    ec_tracks_destroy(tracks, dicts->n);
    return 0;
}

void ec_candidates_destroy(ec_candidates_t *cands)
{
    free(cands->a);
    cands->a = NULL;
    cands->n = 0;
}

void ec_write_report(ec_candidates_t *cands, sdict_t *dicts, FILE *fo)
{
    int64 i;
    fprintf(fo, "#sequence\tposition\tbin\tobserved\texpected\tratio\tz\tpvalue\tqvalue\tdecision\n");
    for (i = 0; i < cands->n; i++) {
        ec_candidate_t *c = cands->a + i;
        fprintf(fo, "%s\t%u\t%u\t%.4f\t%.4f\t%.4f\t%.4f\t%.3e\t%.3e\t%s\n",
            dicts->s[c->seq].name, c->pos, c->bin, c->observed, c->expected, c->ratio,
            c->z, c->p, c->q, c->accepted? "BREAK" : "REJECT");
    }
}

void ec_write_agp(ec_candidates_t *cands, sdict_t *dicts, FILE *fo)
{
    int64 i;
    uint32 s;
    for (s = 0, i = 0; s < dicts->n; s++) {
        uint32 beg = 0, end, part = 0, piece = 0;
        while (i < cands->n && cands->a[i].seq < s)
            i++;
        while (i < cands->n && cands->a[i].seq == s) {
            if (cands->a[i].accepted) {
                end = cands->a[i].pos;
                if (end > beg) {
                    fprintf(fo, "%s_ec_%u\t1\t%u\t1\tW\t%s\t%u\t%u\t+\n",
                        dicts->s[s].name, ++piece, end - beg, dicts->s[s].name, beg + 1, end);
                    beg = end;
                }
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

