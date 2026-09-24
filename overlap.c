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
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#include "kvec.h"
#include "kthread.h"
#include "sdict.h"
#include "misc.h"

#include "overlap.h"

int aln_coords_cmpfunc(const void *a, const void *b)
{ 
    int64  xm, ym;
    aln_t *x = (aln_t *) a;
    aln_t *y = (aln_t *) b;
    
    xm = x->aread;
    ym = y->aread;

    if (xm == ym) {
        xm = x->abpos;
        ym = y->abpos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->aepos;
        ym = y->aepos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bread;
        ym = y->bread;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bbpos;
        ym = y->bbpos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bepos;
        ym = y->bepos;
    }
    return (xm > ym) - (xm < ym);
}

aln_t *add_dual_alignments(aln_t *alns, int64 naln, int64 *_naln)
{
    int64 i;
    aln_t *aln;

    if (_naln) *_naln = naln;
    if (alns == NULL || naln <= 0)
        return alns;
    
    MYREALLOC(alns, naln*2);
    if (alns == NULL)
        mem_alloc_error("alns span");

    for (i = 0; i < naln; i++) {
        aln = &alns[i];
        alns[naln+i] = (aln_t){aln->bread, aln->top, aln->aread, aln->rev, aln->bbpos, aln->bepos, aln->abpos, aln->aepos, aln->mlen, aln->next};
    }

    *_naln = naln*2;

    fprintf(stderr, "[M::%s] added %lld symmetric alignments\n", __func__, naln);

    return alns;
}

/************************** Orders ************************/

int ord_i32_acmpfunc(const void *a, const void *b)
{
    int32 x, y;
    x = ((ord_i32_t *) a)->event;
    y = ((ord_i32_t *) b)->event;
    if (x != y) return (x > y) - (x < y);
    uint32 s, t;
    s = ((ord_i32_t *) a)->which;
    t = ((ord_i32_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_i32_dcmpfunc(const void *a, const void *b)
{
    int32 x, y;
    x = ((ord_i32_t *) a)->event;
    y = ((ord_i32_t *) b)->event;
    if (x != y) return (x < y) - (x > y);
    uint32 s, t;
    s = ((ord_i32_t *) a)->which;
    t = ((ord_i32_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_i32_wcmpfunc(const void *a, const void *b)
{
    uint32 s, t;
    s = ((ord_i32_t *) a)->which;
    t = ((ord_i32_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_u64_acmpfunc(const void *a, const void *b)
{
    uint64 x, y;
    x = ((ord_u64_t *) a)->event;
    y = ((ord_u64_t *) b)->event;
    if (x != y) return (x > y) - (x < y);
    uint32 s, t;
    s = ((ord_u64_t *) a)->which;
    t = ((ord_u64_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_u64_dcmpfunc(const void *a, const void *b)
{
    uint64 x, y;
    x = ((ord_u64_t *) a)->event;
    y = ((ord_u64_t *) b)->event;
    if (x != y) return (x < y) - (x > y);
    uint32 s, t;
    s = ((ord_u64_t *) a)->which;
    t = ((ord_u64_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_u64_wcmpfunc(const void *a, const void *b)
{
    uint32 s, t;
    s = ((ord_u64_t *) a)->which;
    t = ((ord_u64_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_i64_acmpfunc(const void *a, const void *b)
{
    int64 x, y;
    x = ((ord_i64_t *) a)->event;
    y = ((ord_i64_t *) b)->event;
    if (x != y) return (x > y) - (x < y);
    uint32 s, t;
    s = ((ord_i64_t *) a)->which;
    t = ((ord_i64_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_i64_dcmpfunc(const void *a, const void *b)
{
    int64 x, y;
    x = ((ord_i64_t *) a)->event;
    y = ((ord_i64_t *) b)->event;
    if (x != y) return (x < y) - (x > y);
    uint32 s, t;
    s = ((ord_i64_t *) a)->which;
    t = ((ord_i64_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_i64_wcmpfunc(const void *a, const void *b)
{
    uint32 s, t;
    s = ((ord_i64_t *) a)->which;
    t = ((ord_i64_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_dbl_acmpfunc(const void *a, const void *b)
{
    double x, y;
    x = ((ord_dbl_t *) a)->event;
    y = ((ord_dbl_t *) b)->event;
    if (x != y) return (x > y) - (x < y);
    uint32 s, t;
    s = ((ord_dbl_t *) a)->which;
    t = ((ord_dbl_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_dbl_dcmpfunc(const void *a, const void *b)
{
    double x, y;
    x = ((ord_dbl_t *) a)->event;
    y = ((ord_dbl_t *) b)->event;
    if (x != y) return (x < y) - (x > y);
    uint32 s, t;
    s = ((ord_dbl_t *) a)->which;
    t = ((ord_dbl_t *) b)->which;
    return (s > t) - (s < t);
}

int ord_dbl_wcmpfunc(const void *a, const void *b)
{
    uint32 s, t;
    s = ((ord_dbl_t *) a)->which;
    t = ((ord_dbl_t *) b)->which;
    return (s > t) - (s < t);
}

/************************ END Orders **********************/

/************************ Chaining ************************/

// tree node
#define NONE 0   
#define HEAD 1
#define LINK 2
#define DELT 3

typedef struct aln_node {
    struct aln_node *L, *R;
    struct aln_node *next;
    int abpos, aepos;
    int bbpos, bepos;
    uint32 alen;
    uint32 clen:30, type:2;
    uint32 which;
    double score;
} aln_node_t;

static inline int BPOSX(aln_node_t *node)
{ return node->abpos; }

static inline int BPOSY(aln_node_t *node)
{ return node->bbpos; }

static inline int EPOSX(aln_node_t *node)
{ return node->aepos; }

static inline int EPOSY(aln_node_t *node)
{ return node->bepos; }

static int aln_node_coords_cmpfunc(const void *a, const void *b)
{
    int64  xm, ym;
    aln_node_t *x = (aln_node_t *) a;
    aln_node_t *y = (aln_node_t *) b;
    
    xm = x->abpos;
    ym = y->abpos;

    if (xm == ym) {
        xm = x->aepos;
        ym = y->aepos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bbpos;
        ym = y->bbpos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bepos;
        ym = y->bepos;
    }
    return (xm > ym) - (xm < ym);
}

/************************ Splay tree **********************/

#define ALIGN_OVERLAP_BASE 100
#define ALIGN_OVERLAP_FRAC 0.1

static aln_node_t *splay(aln_node_t *v, int x)    //  Assumes x is in the tree
{
    aln_node_t *u, *n;

    if (v == NULL || x == v->bepos)
        return v;
    if (x < v->bepos) {
        u = v->L;
        if (x == u->bepos) {
            v->L = u->R;
            u->R = v;
            return u;
        }
        if (x < u->bepos) {
            n = splay(u->L,x);
            v->L = u->R;
            u->R = v;
            u->L = n->R;
            n->R = u;
        } else {
            n = splay(u->R,x);
            v->L = n->R;
            u->R = n->L;
            n->L = u;
            n->R = v;
        }
    } else {
        u = v->R;
        if (x == u->bepos) {
            v->R = u->L;
            u->L = v;
            return u;
        }
        if (x > u->bepos) {
            n = splay(u->R,x);
            v->R = u->L;
            u->L = v;
            u->R = n->L;
            n->L = u;
        } else {
            n = splay(u->L,x);
            v->R = n->L;
            u->L = n->R;
            n->R = u;
            n->L = v;
        }
    }
    return n;
}

static aln_node_t *splay_find(aln_node_t *v, int x)   //  Find v s.t. v->bepos <= x && x < v->next->bepos
{
    aln_node_t *u;

    if (v == NULL || v->bepos == x)
        return v;
    if (x < v->bepos)
        return splay_find(v->L,x);
    else {
        u = splay_find(v->R,x);
        if (u == NULL)
            return v;
        else
            return u;
    }
}

static aln_node_t *splay_next(aln_node_t *v, int x, aln_node_t *w)
{ 
    if (v == NULL)
        return w;
    if ( x == v->bepos) {
        if (v->R != NULL) {
            w = v->R;
            while (w->L != NULL)
                w = w->L;
        }
        return w;
    }
    if (x < v->bepos)
        return splay_next(v->L,x,v);
    else
        return splay_next(v->R,x,w);
}

static aln_node_t *splay_join(aln_node_t *v, aln_node_t *w)
{ 
    aln_node_t *p;

    if (v == NULL)
        return w;
    for (p = v; p->R != NULL; p = p->R)
        ;
    v = splay(v,p->bepos);
    v->R = w;

    return v;
}

static aln_node_t *splay_insert(aln_node_t *v, aln_node_t *new)
{ 
    aln_node_t *u, *p;

    if (v == NULL)
        return new;
    
    u = splay_find(v,new->bepos);
    if (u != NULL && u->R == NULL)
        u->R = new;
    else {
        if (u == NULL)
            p = v;
        else  // u->R == NULL
            p = u->R;
        while (p->L != NULL)
            p = p->L;
        p->L = new;
    }
    
    return splay(v,new->bepos);
}

static aln_node_t *splay_delete(aln_node_t *v, aln_node_t *old)
{
    aln_node_t *u, *w;

    u = splay_find(v,old->bepos);
    if (u == NULL || u->bepos != old->bepos)
        return NULL;
    v = splay(v,old->bepos);
    w = splay_join(v->L,v->R);

    return w;
}

static int splay_ord_acmpfunc(const void *a, const void *b)
{
    ord_i32_t *x = (ord_i32_t *) a;
    ord_i32_t *y = (ord_i32_t *) b;
    int xm, ym;

    xm = abs(x->event);
    ym = abs(y->event);
    
    if (xm == ym) {
        xm = x->event;
        ym = y->event;
    }
    
    return (xm > ym) - (xm < ym);
}

/********************** END Splay tree ********************/

/************************* RangeSet ***********************/

typedef struct {
    int which;
    double event;
} order_t;

int absint_ord_acmpfunc(const void *a, const void *b)
{
    int x = abs(*(int *) a);
    int y = abs(*(int *) b);

    return (x > y) - (x < y);
}

void range_vec_destroy(range_vec_t *vecs, int n)
{
    if (vecs) {
        if (!n) n = 1;
        int i;
        for (i = 0; i < n; i++)
            free(vecs[i].a);
        free(vecs);
    }
}

int range_acmpfunc(const void *a, const void *b)
{
    range_t *x = (range_t *) a;
    range_t *y = (range_t *) b;
    int xm, ym;
    
    xm = x->beg;
    ym = y->beg;
  
    if (xm == ym) {
        xm = x->end;
        ym = y->end;
    }
    return (xm > ym) - (xm < ym);
}

int rangelist_size(range_t *ranges, int n, int sorted)
{
    if (n <= 0) return 0;
    
    if (!sorted)
        qsort(ranges, n, sizeof(range_t), range_acmpfunc);

    int i;
    int size = 0, beg = ranges[0].beg, end = ranges[0].end;
    for (i = 1; i < n; ++i) {
        if (ranges[i].beg > end) {
            size += end - beg;
            beg = ranges[i].beg;
            end = ranges[i].end;
        } else if (ranges[i].end > end) {
            end = ranges[i].end;
        }
    }
    size += end - beg;
    return size;
}

int rangelist_uniq(range_t *ranges, int n, int sorted, int *_pts)
{
    if (n <= 0) return 0;
    
    if (!sorted)
        qsort(ranges, n, sizeof(range_t), range_acmpfunc);

    int i, num, cov, pos, uniq, *pts;

    num = n * 2;
    if (_pts) pts = _pts;
    else MYMALLOC(pts, num);

    // use sweep line to count the number of unique ranges
    for (i = 0; i < n; i++) {
        pts[i*2]   =  ranges[i].beg;
        pts[i*2+1] = -ranges[i].end;
    }
    qsort(pts, num, sizeof(int), absint_ord_acmpfunc);

    // sweep line
    uniq = cov = 0;
    pos = abs(pts[0]);
    for (i = 0; i < num; i++) {
        if (cov == 1) uniq += abs(pts[i]) - pos;
        pos = abs(pts[i]);
        cov += (pts[i] >= 0) ? 1 : -1;
    }

    if (!_pts) free(pts);

    return uniq;
}

int rangeset_size(range_t *ranges, int n)
{
    int i, c;
    c = 0;
    for (i = 0; i < n; i++)
        c += ranges[i].end - ranges[i].beg;
    return c;
}

int merge_ranges(range_t *ranges, int n, int sorted)
{
    if (n <= 0) return 0;
    
    if (!sorted)
        qsort(ranges, n, sizeof(range_t), range_acmpfunc);

    int i, r;
    r = 0;
    for (i = 1; i < n; i++) {
        if (ranges[i].beg <= ranges[r].end) {
            if (ranges[i].end > ranges[r].end)
                ranges[r].end = ranges[i].end;
        } else {
            ranges[++r] = ranges[i];
        }
    }
    return r + 1;
}

int sorted_range_overlap(range_t *ranges, int n)
{
    if (n <= 0) return 0;

    int i, ovl, end;
    ovl = 0;
    end = ranges[0].end;
    for (i = 1; i < n; i++) {
        if (ranges[i].beg <= end) {
            if (ranges[i].end > end) {
                ovl += end - ranges[i].beg;
                end  = ranges[i].end;
            } else {
                ovl += ranges[i].end - ranges[i].beg;
            }
        } else {
            end = ranges[i].end;
        }
    }
  
    return ovl;
}

int two_sorted_range_overlap(range_t *aranges, int na, range_t *branges, int nb)
{
    int i, j, beg, end, ovl;
    range_t *a, *b;

    ovl = 0;   
    i = j = 0;
    while (i < na && j < nb) {
        a = &aranges[i];
        b = &branges[j];
        
        beg = (a->beg > b->beg)? a->beg : b->beg;
        end = (a->end < b->end)? a->end : b->end;
        
        if (beg < end) ovl += end - beg;
        
        a->end < b->end? i++ : j++;
    }
    
    return ovl;
}

int merge_sorted_ranges(range_t *aranges, int na, range_t *branges, int nb, int *ovl)
{
    int i, j, k;

    i = na - 1;
    j = nb - 1;
    k = na + nb - 1;
    while (i >= 0 && j >= 0) {
        if (range_acmpfunc(aranges+i, branges+j) > 0)
            aranges[k--] = aranges[i--];
        else
            aranges[k--] = branges[j--];
    }
    
    while (j >= 0)
        aranges[k--] = branges[j--];
  
    if (ovl)
        *ovl = sorted_range_overlap(aranges, na + nb);
  
    return merge_ranges(aranges, na + nb, 1);
}

int sorted_range_coverage(range_t *ranges, int n, int beg, int end)
{
    int cov;
    range_t *range;
    cov = 0;
    range = ranges;
    ranges += n;
    for (; range < ranges; range++) {
        if (range->end <= beg)
            continue;
        if (range->beg >= end)
            break;
        cov += MIN(range->end, end) - MAX(range->beg, beg);
    }
    
    return cov; 
}

#define UNSORTED 0
#define SORTED   1

int merge_range_fuzzy(range_t *ranges, int n, int fz, int sorted)
{
    if (n <= 0)
        return 0;

    if (!sorted)
        qsort(ranges, n, sizeof(range_t), range_acmpfunc);
  
    int i, r;
    r = 0;
    for (i = 1; i < n; i++) {
        if (ranges[i].beg <= ranges[r].end + fz) {
            if (ranges[i].end > ranges[r].end)
                ranges[r].end = ranges[i].end;
            }
        else ranges[++r] = ranges[i];
    }
    
    return r + 1;
}

/********************** END RangeSet **********************/

/********************* Glocal Chaining ********************/
// mask for dual 0x77  // 01110111

typedef struct {
    int abpos, aepos;
    int bbpos, bepos;
    int alen, blen;
} alnb_t;

typedef struct {
    aln_t *alns;
    ovl_t *ovls;
    sdict_t *dicts;
    ord_u64_t *asort;
    ord_dbl_t *csort;
    ord_i32_t *psort;
    aln_node_t *nodes;
    uint64 *index;
    range_t *arngs;
    range_t *brngs;
    range_t *xrngs;
    range_t *yrngs;
    alnb_t *alnbs;
    int *sarray;
} chain_data_t;

int ovl_abseqs_cmpfunc(const void *a, const void *b)
{
    int64  xm, ym;
    ovl_t *x = (ovl_t *) a;
    ovl_t *y = (ovl_t *) b;
    
    xm = x->aread;
    ym = y->aread;

    if (xm == ym) {
        xm = x->bread;
        ym = y->bread;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->abpos;
        ym = y->abpos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->aepos;
        ym = y->aepos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bbpos;
        ym = y->bbpos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bepos;
        ym = y->bepos;
    }
    return (xm > ym) - (xm < ym);
}


int ovl_apos_cmpfunc(const void *a, const void *b)
{
    int64  xm, ym;
    ovl_t *x = (ovl_t *) a;
    ovl_t *y = (ovl_t *) b;
    
    xm = x->aread;
    ym = y->aread;

    if (xm == ym) {
        xm = x->abpos;
        ym = y->abpos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->aepos;
        ym = y->aepos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bread;
        ym = y->bread;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bbpos;
        ym = y->bbpos;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->bepos;
        ym = y->bepos;
    }
    return (xm > ym) - (xm < ym);
}


static inline uint8 dual_ovl_type(uint8 type)
{
    return ((type << 4) | (type >> 4)) & 0x77;
}

ovl_t *add_dual_overlaps(ovl_t *ovls, int64 novl, int64 *_novl)
{
    int64 i, n;
    uint8 type;
    ovl_t *ovl;

    if (_novl) *_novl = novl;
    if (ovls == NULL || novl <= 0)
        return ovls;

    // remove deleted overlaps
    n = 0;
    for (i = 0; i < novl; i++)
        if (!ovls[i].del)
            ovls[n++] = ovls[i];

    MYREALLOC(ovls, n*2);
    if (ovls == NULL)
        mem_alloc_error("ovls span");

    for (i = 0; i < n; i++) {
        ovl = &ovls[i];
        type = dual_ovl_type(ovl->type);
        ovls[n+i] = (ovl_t){ovl->bread, ovl->brev^1, ovl->aread, ovl->arev^1, 
            ovl->bbpos, ovl->bepos, ovl->abpos, ovl->aepos, ovl->blen, ovl->alen, 
            type, ovl->del, ovl->neff, ovl->qual, ovl->score};
    }

    *_novl = n*2;

    fprintf(stderr, "[M::%s] added %lld symmetric overlaps\n", __func__, n);

    return ovls;
}


static void backtrack_chain(aln_node_t *node)
{
    if (node->type)
        return;
    
    uint32 clen;
    double score;
    aln_node_t *head, *next;

    // make chain and calculate score
    head = node;
    head->type = HEAD;
    score = node->alen;
    next = node->next;
    clen = 1;
    while (next) {
        if (next->type) {
            node->next = NULL;
            break;
        }
        score += next->alen;
        
        node = next;
        node->type = LINK;
        clen++;

        next = node->next;
    }
    
    head->score = score;
    head->clen  = clen;
}

static void pop_chain(aln_node_t *nodes, int ccnt, ord_dbl_t *order)
{
    int i;
    for (i = 0; i < ccnt; i++) {
        order[i].which = i;
        order[i].event = nodes[i].score;
        nodes[i].type = 0;
    }

    qsort(order, ccnt, sizeof(ord_dbl_t), ord_dbl_dcmpfunc);

    for (i = 0; i < ccnt; i++)
        backtrack_chain(nodes + order[i].which);
}

static int collect_chain(aln_node_t *nodes, int ccnt, ord_dbl_t *order)
{
    int i, nchain = 0; 
    for (i = 0; i < ccnt; i++) {
        if (nodes[i].type == HEAD) {
            order[nchain].which = i;
            order[nchain].event = nodes[i].score;
            nchain++;
        }
    }
    qsort(order, nchain, sizeof(ord_dbl_t), ord_dbl_dcmpfunc);
    
    return nchain;
}

static void build_global_chain(aln_node_t *nodes, int acnt, ord_i32_t *order)
{
    if (!acnt) return;

    int i;
    aln_node_t *node, *e, *v, *w;
    
    for (i = 0; i < acnt; i++) {
        w = nodes + i;
        order[2*i].event   = -w->abpos;
        order[2*i].which   =  i;
        order[2*i+1].event =  w->aepos;
        order[2*i+1].which =  i;
    }

    qsort(order,2*acnt,sizeof(ord_i32_t),splay_ord_acmpfunc);

    node = NULL;
    for (i = 0; i < 2*acnt; i++) {
        e = nodes + order[i].which;
        if (order[i].event <= 0) {
            v = splay_find(node,e->bbpos);
            e->next = v;
            if (v == NULL) {
                e->score = e->alen;
                e->clen = 1;
            } else {
                e->score = v->score + e->alen;
                e->clen = v->clen + 1;
            }
            e->L = NULL;
            e->R = NULL;
        } else {
            v = splay_find(node,e->bepos);
            if (v == NULL)
                node = splay_insert(node,e);
            else if (v->score <= e->score) {
                w = splay_next(node,v->bepos,NULL);
                if (v->bepos >= e->bepos)
                    node = splay_delete(node,v);
                while (w != NULL && w->score <= e->score) {
                    v = w;
                    w = splay_next(node,w->bepos,NULL);
                    node = splay_delete(node,v);
                }
                node = splay_insert(node,e);
            }
        }
    }
}

/*
 * Global chaining with a separable gap penalty.
 *
 * Penalty model:  pen(v -> e) = alpha * ( (e->abpos - v->aepos)
 *                                       + (e->bbpos - v->bepos) )
 * Both gap terms are guaranteed >= 0 by the sweep order on x and by
 * splay_find returning a v with v->bepos <= e->bbpos.
 *
 * The penalty separates as
 *     pen(v -> e) = [alpha*(e->abpos + e->bbpos)]            // depends on e
 *                 - [alpha*(v->aepos + v->bepos)]            // depends on v
 *
 * so for fixed e, maximising  v->score - pen(v -> e)  over v
 * is the same as maximising
 *     KEY(v) = v->score + alpha * (v->aepos + v->bepos)
 *
 * The splay tree therefore maintains a Pareto front on (bepos, KEY):
 * smaller bepos => smaller KEY. splay_find(node, e->bbpos) still
 * returns the in-prefix v with the largest KEY, just as in the
 * original (alpha = 0) algorithm.
 *
 * alpha is in the same units as alen (alignment bp). E.g. alpha = 0.05
 * means 1 bp of unaligned distance costs 0.05 bp of equivalent alignment.
 * alpha = 0 is exactly the original behaviour.
 */
static void build_global_chain_pgap(aln_node_t *nodes, int acnt, ord_i32_t *order, double alpha)
{
    if (!acnt) return;

    int i;
    aln_node_t *node, *e, *v, *w;

#define GC_KEY(n) ((n)->score + alpha * ((double)(n)->aepos + (double)(n)->bepos))

    for (i = 0; i < acnt; i++) {
        w = nodes + i;
        order[2*i  ].event = -w->abpos;
        order[2*i  ].which =  i;
        order[2*i+1].event =  w->aepos;
        order[2*i+1].which =  i;
    }

    qsort(order, 2*acnt, sizeof(ord_i32_t), splay_ord_acmpfunc);

    node = NULL;
    for (i = 0; i < 2*acnt; i++) {
        e = nodes + order[i].which;

        if (order[i].event <= 0) {
            v = splay_find(node, e->bbpos);

            if (v == NULL) {
                /* no candidate predecessor: start a fresh chain */
                e->next  = NULL;
                e->score = e->alen;
                e->clen  = 1;
            } else {
                /* both gaps are >= 0 by construction */
                double gap_a  = (double)(e->abpos - v->aepos);
                double gap_b  = (double)(e->bbpos - v->bepos);
                double gap_pen = alpha * (gap_a + gap_b);
                double link_score = v->score + (double)e->alen - gap_pen;

                if (link_score > (double)e->alen) {
                    /* extending v is profitable */
                    e->next  = v;
                    e->score = link_score;
                    e->clen  = v->clen + 1;
                } else {
                    /* gap is too expensive: better to restart */
                    e->next  = NULL;
                    e->score = e->alen;
                    e->clen  = 1;
                }
            }
            e->L = NULL;
            e->R = NULL;

        } else {
            v = splay_find(node, e->bepos);
            if (v == NULL) {
                node = splay_insert(node, e);
            } else if (GC_KEY(v) <= GC_KEY(e)) {
                w = splay_next(node, v->bepos, NULL);
                if (v->bepos >= e->bepos)
                    node = splay_delete(node, v);
                while (w != NULL && GC_KEY(w) <= GC_KEY(e)) {
                    v = w;
                    w = splay_next(node, w->bepos, NULL);
                    node = splay_delete(node, v);
                }
                node = splay_insert(node, e);
            }
        }
    }

#undef GC_KEY
}

static const double MAX_OVERLAP_FRAC = 0.2;
static const double MAX_OVL_OVERHANG_FRAC = 0.1;
static const double MAX_OVL_OVERHANG_SIZE = 100000;

static const double LOG_MIN_CHAIN_SIZE = 4.0; // log10(10000)
static const double LOG_MAX_CHAIN_SIZE = 6.0; // log10(1000000)
static const double LOG_DIFF_CHAIN_SIZE = LOG_MAX_CHAIN_SIZE - LOG_MIN_CHAIN_SIZE;
static const double SCALE_FACTOR_AT_MIN = .95;
static const double SCALE_FACTOR_AT_MAX = .10;
static const double SCALE_FACTOR_DIFF = SCALE_FACTOR_AT_MAX - SCALE_FACTOR_AT_MIN;
static inline double chain_size_scale_factor(double l)
{
    if (l <= 1.0) return SCALE_FACTOR_AT_MIN;

    double log_L = log10(l);
    if (log_L <= LOG_MIN_CHAIN_SIZE) return SCALE_FACTOR_AT_MIN;
    if (log_L >= LOG_MAX_CHAIN_SIZE) return SCALE_FACTOR_AT_MAX;
    return SCALE_FACTOR_AT_MIN + (log_L - LOG_MIN_CHAIN_SIZE) / LOG_DIFF_CHAIN_SIZE * SCALE_FACTOR_DIFF;
}

static void extend_overlap(double x_beg, double x_end, double y_beg, double y_end, 
    double L_A, double L_B, int rev, 
    double *a_beg, double *a_end, double *b_beg, double *b_end) {
    double dx = x_end - x_beg;
    double dy = y_end - y_beg;

    if (dx == 0) {
        *a_beg = *a_end = x_beg;
        *b_beg = 0; *b_end = L_B;
        return;
    }

    double m_mag = fabs(dy / dx);
    double m = (rev == 0) ? m_mag : -m_mag;
    double x_mid = (x_beg + x_end) / 2.0;
    double y_mid = (y_beg + y_end) / 2.0;
    double b = y_mid - m * x_mid;

    double x_at_y0  = -b / m;
    double x_at_yLB = (L_B - b) / m;
    double y_at_x0  = b;
    double y_at_xLA = m * L_A + b;

    double valid_x[4], valid_y[4];
    int count = 0;
    if (y_at_x0 >= 0  && y_at_x0 <= L_B)  { valid_x[count] = 0;   valid_y[count++] = y_at_x0; }
    if (y_at_xLA >= 0 && y_at_xLA <= L_B) { valid_x[count] = L_A; valid_y[count++] = y_at_xLA; }
    if (x_at_y0 >= 0  && x_at_y0 <= L_A)  { valid_x[count] = x_at_y0;  valid_y[count++] = 0; }
    if (x_at_yLB >= 0 && x_at_yLB <= L_A) { valid_x[count] = x_at_yLB; valid_y[count++] = L_B; }

    if (count < 2) {
        *a_beg = x_beg; *a_end = x_end;
        *b_beg = y_beg; *b_end = y_end;
    } else {
        double min_x = valid_x[0], max_x = valid_x[0];
        double min_y = valid_y[0], max_y = valid_y[0];
        
        for(int i = 1; i < count; i++) {
            if (valid_x[i] < min_x) min_x = valid_x[i];
            if (valid_x[i] > max_x) max_x = valid_x[i];
            if (valid_y[i] < min_y) min_y = valid_y[i];
            if (valid_y[i] > max_y) max_y = valid_y[i];
        }

        *a_beg = min_x; *a_end = max_x;
        *b_beg = min_y; *b_end = max_y;
    }
}

#define SINGLETON_MIN_OVL_FRAC 0.1


// Compute "second-copy" coverage of the primary chain by the secondary chain
// along a single axis (call it the p-axis; the other axis is the q-axis).
//
// Primary chain occupies [p_min, p_max) on the p-axis and [q_min, q_max) on
// the q-axis. Secondary chain occupies [s_min, s_max) and [t_min, t_max).
// Its diagonal runs from (s_min, t_min) to (s_max, t_max) when rev == 0
// (forward) and from (s_min, t_max) to (s_max, t_min) when rev == 1
// (reverse).
//
// A second copy on the p-axis exists wherever the secondary chain's p-range
// overlaps the primary's p-range, regardless of where it sits on the q-axis.
// The piece to *exclude* is where the secondary's diagonal q(p) actually
// passes through the primary q-band [q_min, q_max): there the secondary maps
// the same p-coordinates to the same q-coordinates as the primary, so it is
// not a separate copy.
//
// The resulting (up to two) sub-ranges on the p-axis are appended to rngs.
static inline void find_overlaps_axis(int p_min, int p_max, int q_min, int q_max,
    int s_min, int s_max, int t_min, int t_max, int rev, range_t *rngs, int *nrng)
{
    int dp_s = s_max - s_min;
    int dq_s = t_max - t_min;

    int ip0 = MAX(p_min, s_min);
    int ip1 = MIN(p_max, s_max);
    if (ip0 >= ip1) return;

    // p-range over which the secondary's diagonal q(p) lies inside [q_min, q_max).
    // forward: q(p) = t_min + (dq_s/dp_s)*(p - s_min)
    // reverse: q(p) = t_max - (dq_s/dp_s)*(p - s_min)
    int pd0 = ip1, pd1 = ip0; // empty exclusion by default
    if (dp_s > 0 && dq_s > 0) {
        double s = (double) dp_s / (double) dq_s;
        double pl, pr;
        if (!rev) {
            pl = (double) s_min + s * (double) (q_min - t_min);
            pr = (double) s_min + s * (double) (q_max - t_min);
        } else {
            pl = (double) s_min + s * (double) (t_max - q_max);
            pr = (double) s_min + s * (double) (t_max - q_min);
        }
        if (pl > pr) { double tmp = pl; pl = pr; pr = tmp; }
        int pil = (int) floor(pl);
        int pir = (int) ceil(pr);
        pd0 = MAX(ip0, pil);
        pd1 = MIN(ip1, pir);
    } else if (dp_s > 0 && dq_s == 0) {
        // secondary is flat on the q-axis: exclude the whole p-overlap iff its
        // constant q-coordinate lies inside the primary q-band
        if (t_min >= q_min && t_min < q_max) { pd0 = ip0; pd1 = ip1; }
    }

    if (pd0 < pd1) {
        if (ip0 < pd0) rngs[(*nrng)++] = (range_t) {ip0, pd0};
        if (pd1 < ip1) rngs[(*nrng)++] = (range_t) {pd1, ip1};
    } else {
        rngs[(*nrng)++] = (range_t) {ip0, ip1};
    }
}

static const double CHAIN_PGAP_SCALE = 0.25;
static const double CHAIN_MIN_SCORE = 10000;
static const double CHAIN_MIN_SS_RATIO = 0.1;
static const double CHAIN_PGAP_SCALE_MACRO = 0.5;

static int alnb_abpos_cmpfunc(const void *a, const void *b)
{
    int x = ((const alnb_t *)a)->abpos;
    int y = ((const alnb_t *)b)->abpos;
    return (x > y) - (x < y);
}

static int spatial_merge(alnb_t *alnbs, int nalnb, double alpha, int *alens)
{
    if (nalnb <= 1) return nalnb;

    int i, j, k, gap_x, gap_y, any;
    double geomean_a, geomean_b;
    
    qsort(alnbs, nalnb, sizeof(alnb_t), alnb_abpos_cmpfunc);

    do {
        any = 0;

        alens[nalnb-1] = alnbs[nalnb-1].alen;
        for (j = nalnb - 2; j >= 0; j--)
            alens[j] = MAX(alnbs[j].alen, alens[j+1]);

        for (i = 0; i < nalnb; i++) {
            k = i + 1;
            for (j = i + 1; j < nalnb; j++) {
                gap_x = alnbs[j].abpos - alnbs[i].aepos;
                if (gap_x < 0) gap_x = 0;

                geomean_a = sqrt((double)alnbs[i].alen * (double)alens[j]);
                if (geomean_a * alpha < (double)gap_x) {
                    if (k != j)
                        for (; j < nalnb; j++) alnbs[k++] = alnbs[j];
                    else
                        k = nalnb;
                    break;
                }

                gap_y = MAX(alnbs[i].bbpos, alnbs[j].bbpos) - MIN(alnbs[i].bepos, alnbs[j].bepos);
                if (gap_y < 0) gap_y = 0;

                geomean_a = sqrt((double)alnbs[i].alen * (double)alnbs[j].alen);
                geomean_b = sqrt((double)alnbs[i].blen * (double)alnbs[j].blen);

                if (geomean_a * alpha >= (double)gap_x &&
                    geomean_b * alpha >= (double)gap_y) {
                    alnbs[i].aepos = MAX(alnbs[i].aepos, alnbs[j].aepos);
                    alnbs[i].bbpos = MIN(alnbs[i].bbpos, alnbs[j].bbpos);
                    alnbs[i].bepos = MAX(alnbs[i].bepos, alnbs[j].bepos);
                    alnbs[i].alen += alnbs[j].alen;
                    alnbs[i].blen += alnbs[j].blen;
                    any = 1;
                } else {
                    if (k != j) alnbs[k] = alnbs[j];
                    k++;
                }
            }
            nalnb = k;
        }
    } while (any);

    return nalnb;
}

static void build_adaptive_chain_core(void *_data, long jid, int tid)
{
    chain_data_t *data;
    aln_t *alns, *aln;
    ovl_t *ovl;
    alnb_t *alnbs;
    sdict_t *dicts;
    ord_u64_t *asort;
    ord_dbl_t *csort;
    ord_i32_t *psort, *p_fwd, *p_rev;
    aln_node_t *nodes, *node;
    range_t *arngs, *brngs, *xrngs, *yrngs;
    uint64 *index;
    uint32 alen, blen;
    uint8 cab, cae, cbb, cbe, type;
    double dlen, mlen, score;
    double x_beg, x_end, y_beg, y_end;
    double a_beg, a_end, b_beg, b_end;
    double dx_beg, dx_end, dy_beg, dy_end;
    double rab, rae, rbb, rbe, xlen, ylen;
    double l, u_x, u_y, max_s, max_l, sum_l, sum_l2;
    int rev, narng, nbrng, nxrng, nyrng, nalnb, *sarray;
    int i, j, k, atop, fcnt, acnt, ccnt, pcnt, ovlap;
    
    data = &((chain_data_t *) _data)[tid];
    alns = data->alns;
    dicts = data->dicts;
    asort = data->asort;
    csort = data->csort;
    psort = data->psort;
    nodes = data->nodes;
    index = data->index;
    alnbs = data->alnbs;
    arngs = data->arngs;
    brngs = data->brngs;
    xrngs = data->xrngs;
    yrngs = data->yrngs;
    sarray = data->sarray;
    
    asort += (index[jid] >> 32);
    acnt = (uint32) index[jid];
    alen = dicts->s[asort->event >> 32].len; // aread length
    blen = dicts->s[(uint32) asort->event >> 1].len; // bread length

    // set overlap as deleted
    ovl = data->ovls + jid;
    ovl->del = 1; // mark all overlaps as deleted

    // global chaining with Splay tree
    // initiate splay tree nodes
    // make tree nodes
    atop = fcnt = 0;
    for (i = 0; i < acnt; i++) {
        aln = alns + asort[i].which;
        if (!aln->top) continue;
        
        node = nodes + atop++;
        MYBZERO(node, 1);

        node->abpos = aln->abpos;
        node->aepos = aln->aepos;
        if (aln->rev) {
            node->bbpos = blen - aln->bepos;
            node->bepos = blen - aln->bbpos;
        } else {
            node->bbpos = aln->bbpos;
            node->bepos = aln->bepos;
            fcnt++;
        }
        node->alen  = aln->aepos - aln->abpos + aln->bepos - aln->bbpos;
        node->score = node->alen;
        node->which = asort[i].which;

        ovlap = (int) ((nodes[i].aepos - nodes[i].abpos) * ALIGN_OVERLAP_FRAC + .499);
        if (ovlap > ALIGN_OVERLAP_BASE) ovlap = ALIGN_OVERLAP_BASE;
        nodes[i].aepos -= ovlap;
        ovlap = (int) ((nodes[i].bepos - nodes[i].bbpos) * ALIGN_OVERLAP_FRAC + .499);
        if (ovlap > ALIGN_OVERLAP_BASE) ovlap = ALIGN_OVERLAP_BASE;
        nodes[i].bepos -= ovlap;
        
        nodes[i].L = NULL;
        nodes[i].R = NULL;
        nodes[i].next = NULL;

        nodes[i].clen = 1;
        nodes[i].type = 0;
    }

    // forward and reverse chaining
    build_global_chain_pgap(nodes, fcnt, psort, CHAIN_PGAP_SCALE);
    build_global_chain_pgap(nodes + fcnt, atop - fcnt, psort, CHAIN_PGAP_SCALE);

    // pop global alignment chains
    pop_chain(nodes, atop, csort);
    ccnt = collect_chain(nodes, atop, csort);
    if (!ccnt) goto delete_all;

    // filter by score
    max_s = nodes[csort[0].which].score;
    for (i = 0; i < ccnt; i++) {
        node = nodes + csort[i].which;
        if (node->score < CHAIN_MIN_SCORE ||
            node->score < max_s * CHAIN_MIN_SS_RATIO)
            break;
    }
    ccnt = i;
    if (!ccnt) goto delete_all;

    // greedy selection of primary chains
    // chains largely covered by primary chains are secondary chains
    nxrng = nyrng = 0;
    pcnt = 0;
    for (i = 0; i < ccnt; i++) {
        node = nodes + csort[i].which;
        narng = nbrng = 0;
        while (node) {
            aln = alns + node->which;
            arngs[narng++] = (range_t) {aln->abpos, aln->aepos};
            brngs[nbrng++] = (range_t) {aln->bbpos, aln->bepos};
            node = node->next;
        }
        node = nodes + csort[i].which; // restore node pointer
        narng = merge_ranges(arngs, narng, 0);
        nbrng = merge_ranges(brngs, nbrng, 0);
        ovlap = two_sorted_range_overlap(arngs, narng, xrngs, nxrng);
        if (ovlap >= rangeset_size(arngs, narng) * MAX_OVERLAP_FRAC)
            continue; // secondary
        ovlap = two_sorted_range_overlap(brngs, nbrng, yrngs, nyrng);
        if (ovlap >= rangeset_size(brngs, nbrng) * MAX_OVERLAP_FRAC)
            continue; // secondary
        // add to the accepted chain set
        nxrng = merge_sorted_ranges(xrngs, nxrng, arngs, narng, NULL);
        nyrng = merge_sorted_ranges(yrngs, nyrng, brngs, nbrng, NULL);
        if (i != pcnt) SWAP(ord_dbl_t, csort[i], csort[pcnt]);
        pcnt++;
    }

    // make secondary chains deleted
    for (i = pcnt; i < ccnt; i++)
        nodes[csort[i].which].type = DELT;

    // collect bounding boxes of chains
    nalnb = 0;
    for (i = 0; i < pcnt; i++) {
        node = nodes + csort[i].which;
        a_beg = b_beg = INT32_MAX;
        a_end = b_end = INT32_MIN;
        while (node) {
            aln = alns + node->which;
            if (aln->abpos < a_beg) a_beg = aln->abpos;
            if (aln->aepos > a_end) a_end = aln->aepos;
            if (aln->bbpos < b_beg) b_beg = aln->bbpos;
            if (aln->bepos > b_end) b_end = aln->bepos;
            node = node->next;
        }
        alnbs[nalnb++] = (alnb_t) {a_beg, a_end, b_beg, b_end, a_end-a_beg, b_end-b_beg};
    }
    
    // merge nearby chains to get the final x and y ranges
    spatial_merge(alnbs, nalnb, CHAIN_PGAP_SCALE_MACRO, sarray);

    // find the maximum bounding box of the merged chains
    score = 0;
    j = 0;
    for (i = 0; i < nalnb; i++) {
        if (alnbs[i].alen + alnbs[i].blen > score) {
            score = alnbs[i].alen + alnbs[i].blen;
            j = i;
        }
    }
    x_beg = alnbs[j].abpos;
    x_end = alnbs[j].aepos;
    y_beg = alnbs[j].bbpos;
    y_end = alnbs[j].bepos;

    // delete chains that are out of the final x and y ranges
    for (i = 0; i < pcnt; i++) {
        node = nodes + csort[i].which;
        a_beg = b_beg = INT32_MAX;
        a_end = b_end = INT32_MIN;
        while (node) {
            aln = alns + node->which;
            if (aln->abpos < a_beg) a_beg = aln->abpos;
            if (aln->aepos > a_end) a_end = aln->aepos;
            if (aln->bbpos < b_beg) b_beg = aln->bbpos;
            if (aln->bepos > b_end) b_end = aln->bepos;
            node = node->next;
        }
        if (a_beg >= x_end || a_end <= x_beg || 
            b_beg >= y_end || b_end <= y_beg)
            nodes[csort[i].which].type = DELT;
    }

    // determine the orientation of the syntenic region
    p_fwd = psort;
    p_rev = psort + atop;
    fcnt = 0;
    for (i = 0; i < pcnt; i++) {
        node = nodes + csort[i].which;
        if (node->type != HEAD)
            continue;
        while (node) {
            aln = alns + node->which;
            a_beg = (double) (aln->abpos - x_beg) * (aln->abpos - x_beg);
            b_beg = (double) (aln->bbpos - y_beg) * (aln->bbpos - y_beg);
            a_end = (double) (x_end - aln->aepos) * (x_end - aln->aepos);
            b_end = (double) (y_end - aln->bepos) * (y_end - aln->bepos);
            dx_beg = sqrt(a_beg + b_beg);
            dx_end = sqrt(a_end + b_end);
            dy_beg = sqrt(a_end + b_beg);
            dy_end = sqrt(a_beg + b_end);
            if (dx_beg > dx_end) dx_beg = dx_end;
            if (dy_beg > dy_end) dy_beg = dy_end;
            p_fwd[fcnt] = (ord_i32_t) { node->which, (int) dx_beg };
            p_rev[fcnt] = (ord_i32_t) { node->which, (int) dy_beg };
            fcnt++;
            node = node->next;
        }
    }
    qsort(p_fwd, fcnt, sizeof(ord_i32_t), ord_i32_acmpfunc);
    qsort(p_rev, fcnt, sizeof(ord_i32_t), ord_i32_acmpfunc);
    dlen = sqrt((x_end - x_beg) * (x_end - x_beg) + (y_end - y_beg) * (y_end - y_beg)) / 2.0;
    dx_beg = dy_beg = 0;
    i = j = 0;
    rev = -1;
    for (k = 1; k <= 10; k++) {
        mlen = dlen * k / 10.0;
        while (i < fcnt && p_fwd[i].event <= mlen) {
            aln = alns + p_fwd[i].which;
            dx_beg += aln->aepos - aln->abpos + aln->bepos - aln->bbpos;
            i++;
        }
        while (j < fcnt && p_rev[j].event <= mlen) {
            aln = alns + p_rev[j].which;
            dy_beg += aln->aepos - aln->abpos + aln->bepos - aln->bbpos;
            j++;
        }
        if (dx_beg > 0 && dx_beg > dy_beg * 2.0) {
            rev = 0;
            break;
        }
        if (dy_beg > 0 && dy_beg > dx_beg * 2.0) {
            rev = 1;
            break;
        }
    }
    // tie-break by the total length of alignments in the syntenic region
    if (rev < 0) {
        dx_beg = dy_beg = 0;
        for (i = 0; i < fcnt; i++) {
            aln = alns + p_fwd[i].which;
            if (aln->rev)
                dy_beg += aln->aepos - aln->abpos + aln->bepos - aln->bbpos;
            else
                dx_beg += aln->aepos - aln->abpos + aln->bepos - aln->bbpos;   
        }
        rev = dx_beg >= dy_beg? 0 : 1;
    }

    // determine if we want to keep this syntenic region
    extend_overlap(x_beg, x_end, y_beg, y_end, alen, blen, rev, &a_beg, &a_end, &b_beg, &b_end);
    if (a_end - a_beg < b_end - b_beg) {
        mlen = a_end - a_beg;
        mlen *= chain_size_scale_factor(mlen);
        if (x_end - x_beg < mlen)
            goto delete_all;
    } else {
        mlen = b_end - b_beg;
        mlen *= chain_size_scale_factor(mlen);
        if (y_end - y_beg < mlen)
            goto delete_all;
    }

    // decide if to use the extended boundary for the final syntenic region
    // allow a maximum overhang of 10% at both ends
    xlen = x_end - x_beg;
    ylen = y_end - y_beg;
    if (x_beg - a_beg <= xlen * MAX_OVL_OVERHANG_FRAC)
        x_beg = a_beg;
    if (a_end - x_end <= xlen * MAX_OVL_OVERHANG_FRAC)
        x_end = a_end;
    if (y_beg - b_beg <= ylen * MAX_OVL_OVERHANG_FRAC)
        y_beg = b_beg;
    if (b_end - y_end <= ylen * MAX_OVL_OVERHANG_FRAC)
        y_end = b_end;

    // mark a retained overlap
    ovl->aread = asort->event >> 32;
    ovl->bread = (uint32) asort->event >> 1;
    ovl->abpos = (uint32) x_beg;
    ovl->aepos = (uint32) x_end;
    ovl->bbpos = (uint32) y_beg;
    ovl->bepos = (uint32) y_end;
    ovl->alen = (uint32) (x_end - x_beg);
    ovl->blen = (uint32) (y_end - y_beg);
    
    // determine overlap type
    rab = (double) x_beg / alen;
    rae = (double) (alen - x_end) / alen;
    rbb = (double) y_beg / blen;
    rbe = (double) (blen - y_end) / blen;
    cab = (rab > MAX_OVL_OVERHANG_FRAC || x_beg > MAX_OVL_OVERHANG_SIZE);
    cae = (rae > MAX_OVL_OVERHANG_FRAC || (alen - x_end) > MAX_OVL_OVERHANG_SIZE);
    cbb = (rbb > MAX_OVL_OVERHANG_FRAC || y_beg > MAX_OVL_OVERHANG_SIZE);
    cbe = (rbe > MAX_OVL_OVERHANG_FRAC || (blen - y_end) > MAX_OVL_OVERHANG_SIZE);

    ovl->arev = 0;
    ovl->brev = rev;
    type = 0;
    if (cab && !cae) {
        // EXTD
        if ((!cbb && cbe && !rev) || (cbb && !cbe && rev)) {
            type |= OVL_EXTD;
        }
    } else if (!cab && cae) {
        // EXTD
        if ((cbb && !cbe && !rev) || (!cbb && cbe && rev)) {
            type |= OVL_EXTD;
            ovl->arev ^= 1;
            ovl->brev ^= 1;
        }
    } else if (!cab && !cae) {
        // ACONT
        type |= OVL_ACTD;
    }
    
    if (!cbb && !cbe) {
        // BCNT
        type |= OVL_BCTD;
    }

    if (!type) type |= OVL_INTL;

    ovl->type = type;

    // special check for single-alignment overlaps marked as extensions
    // which are likely chimeric joins by assembler or spurious hits
    // TODO is this a good heuristic?
    if ((type & OVL_EXTD) || (type & OVL_INTL)) {
        k = j = 0;
        for (i = 0; i < pcnt; i++) {
            node = nodes + csort[i].which;
            if (node->type != HEAD)
                continue;
            while (node) {
                if (k++) break;
                j = node->which;
                node = node->next;
            }
            if (k > 1) break;
        }
        if (k == 1 && alns[j].mlen < MIN(alen, blen) * SINGLETON_MIN_OVL_FRAC) {
            ovl->del = 1; // mark as deleted
            goto delete_all;
        }
    }

    // mark as a valid overlap to be retained
    ovl->del = 0;

    // check if there are more than one copy of the chain
    // to determine the overlap quality - the fraction of unique regions
    nxrng = nyrng = 0;
    for (i = 0; i < ccnt; i++) {
        node = nodes + csort[i].which;
        if (node->type == HEAD) // this is primary chain
            continue;
        rev = alns[node->which].rev;
        a_beg = b_beg = INT32_MAX;
        a_end = b_end = INT32_MIN;
        while (node) {
            aln = alns + node->which;
            if (aln->abpos < a_beg) a_beg = aln->abpos;
            if (aln->aepos > a_end) a_end = aln->aepos;
            if (aln->bbpos < b_beg) b_beg = aln->bbpos;
            if (aln->bepos > b_end) b_end = aln->bepos;
            node = node->next;
        }
        find_overlaps_axis(ovl->abpos, ovl->aepos, ovl->bbpos, ovl->bepos,
            a_beg, a_end, b_beg, b_end, rev, xrngs, &nxrng);
        find_overlaps_axis(ovl->bbpos, ovl->bepos, ovl->abpos, ovl->aepos,
            b_beg, b_end, a_beg, a_end, rev, yrngs, &nyrng);
    }
    u_x = 1.0 - (double) rangelist_size(xrngs, nxrng, 0) / (ovl->alen + 1e-6);
    u_y = 1.0 - (double) rangelist_size(yrngs, nyrng, 0) / (ovl->blen + 1e-6);
    ovl->qual = MIN(u_x, u_y);

    // score the overlap by the total length of aligned regions in the syntenic region
    // cacluate the effective number of alignments
    narng = nbrng = 0;
    max_l = .0;
    for (i = 0; i < pcnt; i++) {
        node = nodes + csort[i].which;
        if (node->type != HEAD)
            continue;
        while (node) {
            aln = alns + node->which;
            arngs[narng++] = (range_t) {aln->abpos, aln->aepos};
            brngs[nbrng++] = (range_t) {aln->bbpos, aln->bepos};
            l = alns[node->which].mlen;
            if (l > max_l) max_l = l;
            node = node->next;
        }
    }
    ovl->score = rangelist_size(arngs, narng, 0) / 2.0 + rangelist_size(brngs, nbrng, 0) / 2.0;
    // effective number of alignments
    sum_l = sum_l2 = .0;
    for (i = 0; i < pcnt; i++) {
        node = nodes + csort[i].which;
        if (node->type != HEAD)
            continue;
        while (node) {
            l = alns[node->which].mlen;
            sum_l += l / max_l;
            l /= max_l;
            sum_l2 += l * l;
            node = node->next;
        }
    }
    ovl->neff = sum_l * sum_l / sum_l2;

#ifdef DEBUG_ALN_GLOBAL_CHAIN
    char *aname, *bname;
    aname = dicts->s[aln->aread].name;
    bname = dicts->s[aln->bread].name;
    for (i = 0; i < pcnt; i++) {
        node = nodes + csort[i].which;
        if (node->type != HEAD) continue;
        while (node) {
            aln = alns + node->which;
            fprintf(stdout, "%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t255\tPG:A:%c\n", 
                    aname,
                    alen,
                    aln->abpos,
                    aln->aepos,
                    "+-"[aln->rev],
                    bname,
                    blen,
                    aln->bbpos,
                    aln->bepos,
                    aln->aepos - aln->abpos,
                    aln->aepos - aln->abpos,
                    PG_GLOBAL);
            node = node->next;
        }
    }
#endif

    return;

delete_all:
    // mark all alignments as deleted
    for (i = 0; i < acnt; i++)
        alns[asort[i].which].top = 0;
    return;
}

static inline int number_alns(aln_t *aln)
{
    int n;
    n = 0;
    while (aln) {
        n++;
        aln = aln->next;
    }
    return n;
}

ovl_t *build_adaptive_chains(aln_t *alns, int64 naln, sdict_t *dicts, int n_threads, int64 *_naln, int64 *_novl)
{
    if (naln <= 0) return NULL;
    
    int64 i, j, m, max, ntop;
    uint64 read;
    aln_t *aln;
    ovl_t *ovls;
    ord_u64_t *order;
    chain_data_t *data;
    kvec_t(uint64) index;

    // order alignments by aread, bread, strand
    MYMALLOC(order, naln);
    if (order == NULL)
        mem_alloc_error("order array");
    ntop = 0;
    for (i = 0; i < naln; i++) {
        aln = alns + i;
        if (!aln->top) continue;
        order[ntop].which = i;
        order[ntop].event = ((uint64) aln->aread << 32) | ((uint32) aln->bread << 1 | aln->rev);
        ntop++;
    }
    qsort(order, ntop, sizeof(ord_u64_t), ord_u64_acmpfunc);

    if (!ntop) {
        fprintf(stderr, "[W::%s] no alignments for adaptive chaining\n", __func__);
        free(order);
        return NULL;
    }

    // build order indices
    kv_init(index);
    kv_resize(uint64, index, 1<<16);
    max = 0;
    j = 0;
    read = order[j].event >> 1;
    m = 0;
    for (i = 0; i < ntop; i++) {
        if (read != (order[i].event >> 1)) {
            if ((read >> 31) < (read & 0x7FFFFFFFULL)) {
                kv_push(uint64, index, (uint64) j << 32 | (i - j));
                if (m > max) max = m;
            }
            j = i;
            read = order[j].event >> 1;
            m = 0;
        }
        m += number_alns(alns + order[i].which);
    }
    if ((read >> 31) < (read & 0x7FFFFFFFULL)) {
        kv_push(uint64, index, (uint64) j << 32 | (i - j));
        if (m > max) max = m;
    }

    MYCALLOC(ovls, index.n);
    if (ovls == NULL)
        mem_alloc_error("overlap array");

    // pre-allocate spaces need for core chaining algorithm
    MYCALLOC(data, n_threads);
    if (data == NULL)
        mem_alloc_error("chain data array");
    for (i = 0; i < n_threads; i++) {
        data[i].alns = alns;
        data[i].ovls = ovls;
        data[i].dicts = dicts;
        data[i].asort = order;
        data[i].index = index.a;
        MYMALLOC(data[i].psort, max*2);
        MYMALLOC(data[i].csort, max);
        MYMALLOC(data[i].nodes, max);
        MYMALLOC(data[i].alnbs, max);
        MYMALLOC(data[i].arngs, max);
        MYMALLOC(data[i].brngs, max);
        MYMALLOC(data[i].xrngs, max*2); // double size for inplace merge
        MYMALLOC(data[i].yrngs, max*2); // double size for inplace merge
        MYMALLOC(data[i].sarray, max);
        if (data[i].csort == NULL || data[i].psort == NULL || data[i].nodes == NULL ||
            data[i].arngs == NULL || data[i].brngs == NULL || data[i].xrngs == NULL || 
            data[i].yrngs == NULL || data[i].alnbs == NULL || data[i].sarray == NULL)
            mem_alloc_error("chain data prealloc space");
    }

    // kt_for
    kt_for(n_threads, build_adaptive_chain_core, data, index.n);

    // chain stats
    ntop = 0;
    for (i = 0; i < naln; i++)
        if (alns[i].top)
            alns[ntop++] = alns[i];
    if (_naln) *_naln = ntop;
    fprintf(stderr, "[M::%s] adaptive chaining retained %lld alignments\n", __func__, *_naln);

    ntop = 0;
    for (i = 0; i < index.n; i++)
        if (!ovls[i].del)
            ovls[ntop++] = ovls[i];
    if (_novl) *_novl = ntop;
    fprintf(stderr, "[M::%s] found %lld overlapped sequence pairs\n", __func__, *_novl);

    kv_destroy(index);
    for (i = 0; i < n_threads; i++) {
        free(data[i].csort);
        free(data[i].psort);
        free(data[i].nodes);
        free(data[i].alnbs);
        free(data[i].arngs);
        free(data[i].brngs);
        free(data[i].xrngs);
        free(data[i].yrngs);
    }
    free(data);
    free(order);

    return ovls;
}

/****************** END Glocal Chaining *******************/

/********************* Alignment Chord ********************/
int64 make_chord(aln_t *alns, int64 naln, int slen, chord_t *chord, point_info_t *point)
{
    if (alns == NULL)
        return 0;

    while (!alns->top && naln > 0) {
        alns++;
        naln--;
    }

    if (naln <= 0)
        return 0;

    int64 i, ctop;
    int32 ab, ae;
    int wgt; 
    point_info_t *anchor, *last, *fptr, *eptr, *ptr, *qtr;

    ctop = 0;
    if (alns->abpos != 0) {
        chord[ctop].pos = 0;
        chord[ctop].wgt = 0;
        ctop += 1;
    }

    anchor = last = fptr = point;
    wgt = 0;
    for (i = 0; i < naln; i++) {
        if (!alns[i].top)
            continue;

        ab = alns[i].abpos;
        ae = alns[i].aepos;

        // process range [anchor, ab)
        while (anchor != fptr && anchor->pos < ab) {
            wgt += anchor->cnt;
            chord[ctop].pos = anchor->pos;
            chord[ctop].wgt = wgt;
            ctop += 1;

            // release point space
            // add it to the start
            // need to update last->next
            if (anchor == last) {
                fptr = anchor;
            } else {
                ptr = anchor->next;
                anchor->next = fptr;
                fptr = anchor;
                last->next = fptr;
                anchor = ptr;
            }
        }

        // add ab as as anchor
        if (anchor == fptr) {
            anchor->cnt = 0;
            last = anchor;
            fptr = fptr->next;
        } else if (anchor->pos > ab) {
            // need to insert an anchor
            ptr = fptr->next;
            fptr->next = anchor;
            fptr->cnt = 0;
            anchor = fptr;
            fptr = ptr;
            last->next = fptr;
        }
        anchor->pos = ab;
        anchor->cnt += 1;

        // add ae as an anchor
        eptr = anchor;
        qtr = NULL;
        while (eptr != fptr && eptr->pos < ae) {
            qtr = eptr;
            eptr = eptr->next;
        }
        if (eptr == fptr) {
            eptr->cnt = 0;
            last = eptr;
            fptr = fptr->next;
        } else if (eptr->pos > ae) {
            // need to insert an anchor
            ptr = fptr->next;
            fptr->next = eptr;
            fptr->cnt = 0;
            eptr = fptr;
            if (qtr)
                qtr->next = eptr;
            fptr = ptr;
            last->next = fptr;
        }
        eptr->pos = ae;
        eptr->cnt -= 1;
    }

    // process remaining anchors
    for (ptr = anchor; ptr != fptr; ptr = ptr->next) {
        wgt += ptr->cnt;
        chord[ctop].pos = ptr->pos;
        chord[ctop].wgt = wgt;
        ctop += 1;
    }

    // add the last point if not reach the end
    if (chord[ctop-1].pos != slen) {
        chord[ctop].pos = slen;
        chord[ctop].wgt = 0;
        ctop += 1;
    }

    return ctop;
}

int64 make_chord_from_ranges(range_t *ranges, int64 n, chord_t *chord, point_info_t *point, int sorted)
{
    if (ranges == NULL || n <= 0)
        return 0;

    if (!sorted)
        qsort(ranges, n, sizeof(range_t), range_acmpfunc);

    int64 i, ctop;
    int32 ab, ae;
    int wgt; 
    point_info_t *anchor, *last, *fptr, *eptr, *ptr, *qtr;

    anchor = last = fptr = point;
    ctop = 0;
    wgt = 0;
    for (i = 0; i < n; i++) {
        ab = ranges[i].beg;
        ae = ranges[i].end;

        // process range [anchor, ab)
        while (anchor != fptr && anchor->pos < ab) {
            wgt += anchor->cnt;
            chord[ctop].pos = anchor->pos;
            chord[ctop].wgt = wgt;
            ctop += 1;

            // release point space
            // add it to the start
            // need to update last->next
            if (anchor == last) {
                fptr = anchor;
            } else {
                ptr = anchor->next;
                anchor->next = fptr;
                fptr = anchor;
                last->next = fptr;
                anchor = ptr;
            }
        }

        // add ab as as anchor
        if (anchor == fptr) {
            anchor->cnt = 0;
            last = anchor;
            fptr = fptr->next;
        } else if (anchor->pos > ab) {
            // need to insert an anchor
            ptr = fptr->next;
            fptr->next = anchor;
            fptr->cnt = 0;
            anchor = fptr;
            fptr = ptr;
            last->next = fptr;
        }
        anchor->pos = ab;
        anchor->cnt += 1;

        // add ae as an anchor
        eptr = anchor;
        qtr = NULL;
        while (eptr != fptr && eptr->pos < ae) {
            qtr = eptr;
            eptr = eptr->next;
        }
        if (eptr == fptr) {
            eptr->cnt = 0;
            last = eptr;
            fptr = fptr->next;
        } else if (eptr->pos > ae) {
            // need to insert an anchor
            ptr = fptr->next;
            fptr->next = eptr;
            fptr->cnt = 0;
            eptr = fptr;
            if (qtr)
                qtr->next = eptr;
            fptr = ptr;
            last->next = fptr;
        }
        eptr->pos = ae;
        eptr->cnt -= 1;
    }

    // process remaining anchors
    for (ptr = anchor; ptr != fptr; ptr = ptr->next) {
        wgt += ptr->cnt;
        chord[ctop].pos = ptr->pos;
        chord[ctop].wgt = wgt;
        ctop += 1;
    }

    return ctop;
}

/****************** END Alignment Chord *******************/
