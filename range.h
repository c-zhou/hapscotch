/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2025 Chenxi Zhou <chnx.zhou@gmail.com>                          *
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

/********************************** Revision History *****************************
 *                                                                               *
 * 15/01/25 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#ifndef RANGE_H_
#define RANGE_H_

#include "kavl.h"
#include "misc.h"

#define RANGETYPE int64
#define NODECHUNK 1024

typedef struct RangeNode {
    RANGETYPE beg;      // start of range (inclusive)
    RANGETYPE end;      // end of range (exclusive)
    int32 coverage;     // coverage of this range
    KAVL_HEAD(struct RangeNode) head;
} RangeNode;

typedef struct NodeList {
    int n;
    struct NodeList *next;
    RangeNode space[NODECHUNK];
} NodeList;

typedef struct {
    RangeNode *root;
    NodeList *nodes;
} rangetree_t;

typedef struct {
    RANGETYPE beg;
    RANGETYPE end;
    void *itr;
    RangeNode *val;
} rangetree_itr_t;

#define rangetree_itr_hasnext(itr) ((itr)->val != NULL)
#define rangetree_itr_beg(itr) ((itr)->val->beg)
#define rangetree_itr_end(itr) ((itr)->val->end)
#define rangetree_itr_val(itr) ((itr)->val->coverage)
#define rangetree_size(rt) (kavl_size(head, (rt)->root))

#ifdef __cplusplus
extern "C" {
#endif
rangetree_t *rangetree_create();
void rangetree_init(rangetree_t *rt);
void rangetree_destroy(rangetree_t *rt); // free(rt)
void rangetree_free(rangetree_t *rt); // do not free(rt)
int rangetree_add(rangetree_t *rt, RANGETYPE beg, RANGETYPE end, int32 coverage);
rangetree_itr_t *rangetree_itr_init(rangetree_t *rt, RANGETYPE beg, RANGETYPE end);
rangetree_itr_t *rangetree_itr_new();
void rangetree_itr_set_range(rangetree_itr_t *rt_itr, rangetree_t *rt, RANGETYPE beg, RANGETYPE end);
void rangetree_itr_next(rangetree_itr_t *rt_itr);
void rangetree_itr_destroy(rangetree_itr_t *itr);
RANGETYPE rangetree_coverage(rangetree_t *rt, RANGETYPE beg, RANGETYPE end, int64 *total, int64 *covs, int max_cov);
#ifdef __cplusplus
}
#endif

#endif /* RANGE_H_ */
