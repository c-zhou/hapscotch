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

#include <stdlib.h>

#include "range.h"

#define range_cmp(p, q) (((q)->end < (p)->end) - ((p)->end < (q)->end))

KAVL_INIT(range, RangeNode, head, range_cmp)

static inline NodeList *node_list_alloc()
{
    NodeList *list;
    MYMALLOC(list, 1);
    list->n = 0;
    list->next = NULL;
    return list;
}

rangetree_t *rangetree_create()
{
    rangetree_t *rt;
    MYCALLOC(rt, 1);
    rt->nodes = node_list_alloc();
    return rt;
}

void rangetree_init(rangetree_t *rt)
{
    rt->root = NULL;
    rt->nodes = node_list_alloc();
}

void rangetree_free(rangetree_t *rt)
{
    if (!rt) return;

    NodeList *next, *nodes = rt->nodes;
    while (nodes) {
        next = nodes->next;
        free(nodes);
        nodes = next;
    }
}

void rangetree_destroy(rangetree_t *rt)
{
    rangetree_free(rt);
    free(rt);
}

// Push a range to node buffer
static inline RangeNode *push_node(rangetree_t *rt, RANGETYPE beg, RANGETYPE end, int32 coverage)
{
    NodeList *list = rt->nodes;
    if (list->n == NODECHUNK) {
        rt->nodes = node_list_alloc();
        rt->nodes->next = list;
        list = rt->nodes;
    }
    RangeNode *node = &list->space[list->n++];
    node->beg = beg;
    node->end = end;
    node->coverage = coverage;
    return node;
}

static inline void insert_node(rangetree_t *rt, RangeNode *node)
{
    kavl_insert_range(&rt->root, node, 0);
}

static inline void insert_pending_nodes(rangetree_t *rt, int n)
{
    NodeList *list = rt->nodes;
    RangeNode *nodes = list->space;
    int k = list->n;
    while (n--) {
        insert_node(rt, &nodes[--k]);
        if (!k) {
            list = list->next;
            if (list) {
                nodes = list->space;
                k = list->n;
            }
        }
    }
}

// Add a new range [beg, end)
int rangetree_add(rangetree_t *rt, RANGETYPE beg, RANGETYPE end, int32 coverage)
{   
    if (beg >= end) return 0;

    struct kavl_itr_range itr;
    RangeNode *node, key;
    int n;

    key.end = beg + 1;
    kavl_itr_find_range(rt->root, &key, &itr);
    node = (RangeNode *) kavl_at(&itr);

    n = 0;
    while (node && beg < end) {
        if (beg < node->beg)
            beg = push_node(rt, beg, MIN(node->beg, end), coverage)->end, n++;
        if (end > node->beg) {
            if (beg > node->beg)
                beg = push_node(rt, node->beg, beg, node->coverage)->end, n++;
            if (end < node->end)
                beg = push_node(rt, beg, end, node->coverage + coverage)->end, n++;
            else
                node->coverage += coverage;
            node->beg = beg;
            beg = node->end;
        }
        node = kavl_itr_next_range(&itr)? (RangeNode *) kavl_at(&itr) : NULL;
    }

    if (beg < end) push_node(rt, beg, end, coverage), n++;

    insert_pending_nodes(rt, n);

    return n;
}

void rangetree_itr_destroy(rangetree_itr_t *itr)
{
    if (!itr) return;

    free(itr->itr);
    free(itr);
}

// Get an empty iterator
rangetree_itr_t *rangetree_itr_new()
{   
    rangetree_itr_t *rt_itr;
    struct kavl_itr_range *itr;
    
    MYCALLOC(rt_itr, 1);
    MYCALLOC(itr, 1);

    rt_itr->itr = itr;
    rt_itr->val = NULL;
    
    return rt_itr;
}

void rangetree_itr_set_range(rangetree_itr_t *rt_itr, rangetree_t *rt, RANGETYPE beg, RANGETYPE end)
{
    struct kavl_itr_range *itr;
    RangeNode *node, key;

    itr = rt_itr->itr;
    key.end = beg + 1;

    kavl_itr_find_range(rt->root, &key, itr);
    node = (RangeNode *) kavl_at(itr);

    if (node && node->beg >= end)
        node = NULL;

    rt_itr->beg = beg;
    rt_itr->end = end;
    rt_itr->val = node;
}

// Get a iterator with a filter
rangetree_itr_t *rangetree_itr_init(rangetree_t *rt, RANGETYPE beg, RANGETYPE end)
{   
    if (beg >= end) return 0;

    rangetree_itr_t *rt_itr;
    struct kavl_itr_range *itr;
    RangeNode *node, key;

    MYCALLOC(rt_itr, 1);
    MYCALLOC(itr, 1);
    key.end = beg + 1;

    kavl_itr_find_range(rt->root, &key, itr);
    node = (RangeNode *) kavl_at(itr);
    
    if (node && node->beg >= end)
        node = NULL;

    rt_itr->beg = beg;
    rt_itr->end = end;
    rt_itr->itr = itr;
    rt_itr->val = node;

    return rt_itr;
}

void rangetree_itr_next(rangetree_itr_t *rt_itr)
{
    struct kavl_itr_range *itr = rt_itr->itr;
    RangeNode *node = kavl_itr_next_range(itr)? (RangeNode *) kavl_at(itr) : NULL;
    rt_itr->val = (node && node->beg < rt_itr->end)? node : NULL;
    return;
}

RANGETYPE rangetree_coverage(rangetree_t *rt, RANGETYPE beg, RANGETYPE end, int64 *total, int64 *covs, int max_cov)
{
    RANGETYPE b, e, o, coverage = 0;
    int32 c;
    int64 sum = 0;
    rangetree_itr_t *itr = rangetree_itr_init(rt, beg, end);
    while (rangetree_itr_hasnext(itr)) {
        b = rangetree_itr_beg(itr);
        e = rangetree_itr_end(itr);
        c = rangetree_itr_val(itr);
        o = MIN(e, end) - MAX(b, beg);
        coverage += o;
        sum += o * c;
        if (c > max_cov)
            covs[max_cov+1] += o;
        else
            covs[c] += o;
        rangetree_itr_next(itr);
    }
    rangetree_itr_destroy(itr);
    if (total) *total = sum;
    return coverage;
}

#undef DEBUG_RANGTREE

#ifdef DEBUG_RANGTREE
#include <stdio.h>
int main(int argc, char *argv[])
{
    rangetree_t *rt;
    rangetree_itr_t *itr;

    rt = rangetree_create();
    
    fprintf(stderr, "#Ranges: %d\n", kavl_size(head, rt->root));
    itr = rangetree_itr_new();
    rangetree_itr_set_range(itr, rt, 0, 100);
    while (rangetree_itr_hasnext(itr)) {
        fprintf(stderr, "%lld %lld %d\n", rangetree_itr_beg(itr), rangetree_itr_end(itr), rangetree_itr_val(itr));
        rangetree_itr_next(itr);
    }

    rangetree_add(rt, 10, 20, 1);
    fprintf(stderr, "#Ranges: %d\n", kavl_size(head, rt->root));
    rangetree_itr_set_range(itr, rt, 0, 11);
    while (rangetree_itr_hasnext(itr)) {
        fprintf(stderr, "%lld %lld %d\n", rangetree_itr_beg(itr), rangetree_itr_end(itr), rangetree_itr_val(itr));
        rangetree_itr_next(itr);
    }

    rangetree_itr_destroy(itr);

    rangetree_add(rt, 1, 3, 1);
    rangetree_add(rt, 1, 4, 1);
    rangetree_add(rt, 7, 8, 1);
    rangetree_add(rt, 1, 14, 1);
    rangetree_add(rt, 5, 6, 1);

    fprintf(stderr, "#Ranges: %d\n", kavl_size(head, rt->root));
    itr = rangetree_itr_init(rt, 0, 20);
    while (rangetree_itr_hasnext(itr)) {
        fprintf(stderr, "%lld %lld %d\n", rangetree_itr_beg(itr), rangetree_itr_end(itr), rangetree_itr_val(itr));
        rangetree_itr_next(itr);
    }
    rangetree_itr_destroy(itr);
    rangetree_destroy(rt);
    return 0;
}
#endif
