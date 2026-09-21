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
    int min_frag;
    double med_drop;
    double rec_rate;
    double p_thresh;
} ec_conf_t;

typedef struct {
    int seq;
    uint32 pos;
    double cnts[3];
    double pval;
} ec_pos_t;

extern ec_conf_t ec_conf;

#ifdef __cplusplus
extern "C" {
#endif

ec_pos_t *ec_call_breaks(hic_t *hics, int64 nhic, sdict_t *dicts, int *_ncall);
void ec_write_report(ec_pos_t *calls, int ncall, sdict_t *dicts, FILE *fo);
void ec_write_agp(ec_pos_t *calls, int ncall, sdict_t *dicts, FILE *fo);

#ifdef __cplusplus
}
#endif

#endif /* EC_H_ */