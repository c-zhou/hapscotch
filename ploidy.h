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

// Genome coverage histograms and ploidy number estimation

#ifndef PLOIDY_H_
#define PLOIDY_H_

#include "misc.h"
#include "kvec.h"
#include "sdict.h"
#include "range.h"
#include "overlap.h"

typedef struct {
    int pos;
    int cov;
} cov_point_t;

#define OVL_MIN_QUAL 0.7

typedef struct {
    int seq;
    int beg, end;
} srange_t;

typedef kvec_t(srange_t) srange_vec_t;

// upper limit of genome ploidy number to consider, settable via -P/--max-ploidy
extern int MAX_PLOIDY_NUMBER;

#ifdef __cplusplus 
extern "C" {
#endif

int srange_cmpfunc(const void *a, const void *b);

int64 pts_from_overlaps(void *data, int64 n, void *param, int *points);
int64 pts_from_ranges(void *data, int64 n, void *param, int *points);
int64 pts_from_sranges(void *data, int64 n, void *param, int *points);
int calc_coverage_from_intervals(cov_point_t *covs, int *points, int slen,
    void *data, int64 nd, void *param, int64 (*pts_func) (void *, int64, void *, int *));

double average_range_coverage(cov_point_t *covs, int ncov, int beg, int end);
int total_range_coverage(cov_point_t *covs, int ncov, int beg, int end, int max_cov);
double overlap_min_lowcopy_ratio(double l);

int64 *genome_coverage_histogram(aln_t *alns, int64 naln, sdict_t *dicts, int max_copy, int print_hist);
void report_genome_coverage_histogram(aln_t *alns, int64 naln, sdict_t *dicts, int max_copy);
int estimate_ploidy_number(ovl_t *ovls, int64 novl, sdict_t *dicts);
void genome_coverage_summary(rangetree_t *rt, sdict_t *dicts, int max_copy);

#ifdef __cplusplus
}
#endif

#endif // PLOIDY_H_
