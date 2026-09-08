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

#ifndef EC_H_
#define EC_H_

#include "hic.h"
#include "sdict.h"

typedef struct {
    int bin_size;
    int max_distance;
    int control_window;
    int min_fragment;
    int min_contacts;
    double max_ratio;
    double fdr;
} ec_conf_t;

typedef struct {
    uint32 seq;
    uint32 bin;
    uint32 pos;
    double observed;
    double expected;
    double ratio;
    double z;
    double p;
    double q;
    int accepted;
} ec_candidate_t;

typedef struct {
    ec_candidate_t *a;
    int64 n;
} ec_candidates_t;

void ec_conf_init(ec_conf_t *conf);
int ec_call_breaks(hic_t *hics, int64 nhic, sdict_t *dicts, ec_conf_t *conf, ec_candidates_t *cands);
void ec_candidates_destroy(ec_candidates_t *cands);
void ec_write_report(ec_candidates_t *cands, sdict_t *dicts, FILE *fo);
void ec_write_agp(ec_candidates_t *cands, sdict_t *dicts, FILE *fo);

#endif /* EC_H_ */