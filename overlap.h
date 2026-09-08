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

// Sequence overlap detection: turns raw (PAF) alignments into merged,
// glocally chained sequence overlaps used for ploidy estimation and
// pseudo haplotype/scaffold construction.

#ifndef OVERLAP_H
#define OVERLAP_H

#include "misc.h"
#include "kvec.h"
#include "sdict.h"

// debug tags shared with the ploidy and hap modules
#undef DEBUG_ALN_GLOBAL_CHAIN
#undef DEBUG_GENOME_COV_HIST
#undef DEBUG_GENOME_COV_DEPTH
#undef DEBUG_SV_DETECTION
#undef DEBUG_HIC_LINKAGE
#undef DEBUG_HIC_LINKAGE_GROUP
#undef DEBUG_SCAFFOLD_PARTITION
#undef DEBUG_SCAFFOLD_GROUP_MERGE
#undef DEBUG_HAPLOTYPE_PARTITION
#undef DEBUG_REFINE_HAP_PARTITION
#undef DEBUG_HAPLOTYPE_SORTING
#undef DEBUG_HAPLOTYPE_PHASE

enum DEBUG_TAG {
    PG_GLOBAL = 'G',
    PG_HIST   = 'H',
    PG_DEPTH  = 'D',
    PG_HIC    = 'C',
    PG_GML    = 'M'
};

typedef struct aln {
    uint32 aread:31, top:1;
    uint32 bread:31, rev:1;
    int abpos, aepos;
    int bbpos, bepos;
    uint32 mlen;
    struct aln *next;
} aln_t;

typedef kvec_t(aln_t) aln_vec_t;

typedef struct {
    uint32 which;
    int32  event;
} ord_i32_t;

typedef struct {
    uint32 which;
    uint64 event;
} ord_u64_t;

typedef struct {
    uint32 which;
    int64  event;
} ord_i64_t;

typedef struct {
    uint32 which;
    double event;
} ord_dbl_t;

typedef struct {
    int beg, end;
} range_t;

typedef kvec_t(range_t) range_vec_t;

#define OVL_ACTD 0x01  // 00000001   aread is contained
#define OVL_BCTD 0x10  // 00010000   bread is contained
#define OVL_EXTD 0x22  // 00100010   reads overlap
#define OVL_INTL 0x44  // 01000100   internal matches
#define OVL_REPT 0x08  // 00001000   repeat matches

typedef struct {
    uint32 aread:31, arev:1;
    uint32 bread:31, brev:1;
    int abpos, aepos;
    int bbpos, bepos;
    int alen, blen;
    uint8 type:7, del:1;
    double neff, qual, score;
} ovl_t;

typedef struct {
    int pos;
    int wgt;
} chord_t;

typedef struct point_info {
    struct point_info *next;
    int pos;
    int cnt;
} point_info_t;

int aln_coords_cmpfunc(const void *a, const void *b);
aln_t *add_dual_alignments(aln_t *alns, int64 naln, int64 *_naln);

int ord_i32_acmpfunc(const void *a, const void *b);
int ord_i32_dcmpfunc(const void *a, const void *b);
int ord_i32_wcmpfunc(const void *a, const void *b);
int ord_u64_acmpfunc(const void *a, const void *b);
int ord_u64_dcmpfunc(const void *a, const void *b);
int ord_u64_wcmpfunc(const void *a, const void *b);
int ord_i64_acmpfunc(const void *a, const void *b);
int ord_i64_dcmpfunc(const void *a, const void *b);
int ord_i64_wcmpfunc(const void *a, const void *b);
int ord_dbl_acmpfunc(const void *a, const void *b);
int ord_dbl_dcmpfunc(const void *a, const void *b);
int ord_dbl_wcmpfunc(const void *a, const void *b);

int absint_ord_acmpfunc(const void *a, const void *b);
void range_vec_destroy(range_vec_t *vecs, int n);
int range_acmpfunc(const void *a, const void *b);
int rangelist_size(range_t *ranges, int n, int sorted);
int rangelist_uniq(range_t *ranges, int n, int sorted, int *_pts);
int rangeset_size(range_t *ranges, int n);
int merge_ranges(range_t *ranges, int n, int sorted);
int sorted_range_overlap(range_t *ranges, int n);
int two_sorted_range_overlap(range_t *aranges, int na, range_t *branges, int nb);
int merge_sorted_ranges(range_t *aranges, int na, range_t *branges, int nb, int *ovl);
int sorted_range_coverage(range_t *ranges, int n, int beg, int end);
int merge_range_fuzzy(range_t *ranges, int n, int fz, int sorted);

int ovl_abseqs_cmpfunc(const void *a, const void *b);
int ovl_apos_cmpfunc(const void *a, const void *b);
ovl_t *add_dual_overlaps(ovl_t *ovls, int64 novl, int64 *_novl);
ovl_t *build_adaptive_chains(aln_t *alns, int64 naln, sdict_t *dicts, int n_threads, int64 *_naln, int64 *_novl);

int64 make_chord(aln_t *alns, int64 naln, int slen, chord_t *chord, point_info_t *point);
int64 make_chord_from_ranges(range_t *ranges, int64 n, chord_t *chord, point_info_t *point, int sorted);

#endif // OVERLAP_H
