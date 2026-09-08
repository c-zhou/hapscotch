/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2021 Chenxi Zhou <chnx.zhou@gmail.com>                          *
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
 * 15/04/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#ifndef SDICT_H_
#define SDICT_H_

#include <stdint.h>

#include "khash.h"

#include "misc.h"
#include "agp-spec.h"

extern char comp_table[128];
extern char nucl_toupper[128];

extern AGP_CT_t DEFAULT_AGP_SEQ_COMPONENT_TYPE;
extern AGP_CT_t DEFAULT_AGP_GAP_COMPONENT_TYPE;
extern AGP_GT_t DEFAULT_AGP_GAP_TYPE;
extern AGP_LE_t DEFAULT_AGP_LINKAGE_EVIDENCE;
extern int DEFAULT_AGP_GAP_SIZE;

typedef struct {
    char *name; // seq id
    char *seq; // sequence
    uint32 len; // seq length
} sd_seq_t;

KHASH_MAP_INIT_STR(sdict, uint32)
typedef khash_t(sdict) sdhash_t;

typedef struct {
    uint32 n, m; // n: seq number, m: memory allocated
    sd_seq_t *s; // sequence dictionary
    sdhash_t *h; // sequence hash map: name -> index
} sdict_t;

typedef struct {
    uint32 s; // seq id
    uint32 k; // seg index on seq
    uint64 a; // seq start without counting gaps
    uint64 b; // seq start with counting gaps
    uint32 c, x, y; // subseq c: id << 1 | ori, x: start, y: length
    uint32 t; // component type
    uint32 r; // component orientation: ori = 1 if MINUS otherwise 0
} sd_seg_t;

typedef struct {
    char *name; // seq id
    uint64 len; // seq length
    uint64 gap; // gap length
    uint32 n; // seg number
    uint32 s; // seg start position
    // the following fields contains gap info
    // gn: gap number; gs: seg pos
    // gs is monotonically increasing
    // gs[i] for i-th gap storing the index of the immediate NEXT seg
    // gs[i] could be zero or 'n' implying leading and trailing gaps
    // gsize: gap size; gcomp: component type; gtype: gap type
    // glink: linkage; gevid: linkage evidence
    // it allows blocks of consecutive segs or gaps
    uint32 gn, *gs;
    uint32 *gsize, *gcomp, *gtype, *glink, *gevid;
} sd_aseq_t;

// used to store scaffolds
// quickly convert a contig position to scaffold position
// with or without counting gaps
// consider contig breaks
typedef struct {
    uint32 n, m; // n: seq number, m: seq memory allocated
    sd_aseq_t *s; // sequence dictionary
    sdhash_t *h; // sequence hash map: name -> index
    uint32 u, v; // u: seg number, v: seg memory allocated
    sd_seg_t *seg; // segments
    uint64 *a; // sub sequence index map: id -> start pos << 32 | num segs, need this to deal with sub seq breaks
    uint64 *index; // sub seq end (seg.x + seg.y) << 32 | seg index, used to find the seg index given a sub seq position
    sdict_t *sdict; // sub sequence dictionary
} asm_dict_t;

typedef enum cc_error_code {
    CC_SUCCESS = 0,
    SEQ_NOT_FOUND = 1,
    POS_NOT_IN_RANGE = 2
} CC_ERR_t;

#ifdef __cplusplus
extern "C" {
#endif

sdict_t *sd_init(void);
asm_dict_t *asm_init(sdict_t *dict);
void sd_destroy(sdict_t *d);
void asm_destroy(asm_dict_t *d);
uint32 sd_put(sdict_t *d, const char *name, uint32 len);
uint32 sd_put1(sdict_t *d, const char *name, const char *seq, uint32 len);
sdict_t *make_sdict_from_fa(const char *f, uint32 min_len);
sdict_t *make_sdict_from_index(const char *f, uint32 min_len);
sdict_t *make_sdict_from_gfa(const char *f, uint32 min_len);
asm_dict_t *make_asm_dict_from_sdict(sdict_t *sdict);
void seg_put(asm_dict_t *d, uint32 s, uint32 k, uint64 a, uint64 b,
        uint32 c, uint32 x, uint32 y, uint32 t, uint32 r);
uint32 asm_put(asm_dict_t *d, const char *name, uint64 len, uint64 gap, uint32 n, uint32 s, uint32 g,
        uint32 *gs, uint32 *gsize, uint32 *gcomp, uint32 *gtype, uint32 *glink, uint32 *gevid);
asm_dict_t *make_asm_dict_from_agp(sdict_t *sdict, const char *f, int allow_unknown_oris);
void add_unplaced_short_seqs(asm_dict_t *d, uint32 min_len);
void asm_index(asm_dict_t *d);
char *get_asm_seq(asm_dict_t *d, char *name);
void sd_stats(sdict_t *d, uint64 *n_stats, uint32 *l_stats);
void asm_sd_stats(asm_dict_t *d, uint64 *n_stats, uint32 *l_stats);
void write_fasta_file_from_agp(const char *fa, const char *agp, FILE *fo, int line_wd, int un_oris);
void write_segs_to_agp(sd_seg_t *segs, uint32 n, sd_aseq_t *aseq, sdict_t *sd, char *name, FILE *fo);
void write_sorted_agp(asm_dict_t *dict, FILE *fo);
void write_sdict_to_agp(sdict_t *sdict, FILE *fo);
void write_asm_dict_to_agp(asm_dict_t *dict, FILE *fo);
uint64 write_segs_to_fasta(sd_seg_t *segs, uint32 n, sd_aseq_t *aseq, sdict_t *sd, char *name, int line_wd, FILE *fo);
uint64 write_asm_dict_to_fasta(asm_dict_t *dict, int line_wd, FILE *fo);
#ifdef __cplusplus
}
#endif

static inline uint32 sd_get(sdict_t *d, const char *name)
{
    sdhash_t *h = d->h;
    khint_t k;
    k = kh_get(sdict, h, name);
    return k == kh_end(h)? UINT32_MAX : kh_val(h, k);
}

static inline int sd_exists(sdict_t *d, const char *name, uint32 len, uint32 *sid)
{
    uint32 k = sd_get(d, name);
    if (sid) *sid = k;
    return k != UINT32_MAX && d->s[k].len == len;
}


static inline uint32 asm_sd_get(asm_dict_t *d, const char *name)
{
    sdhash_t *h = d->h;
    khint_t k;
    k = kh_get(sdict, h, name);
    return k == kh_end(h)? UINT32_MAX : kh_val(h, k);
}

/* contig coordinates to scaffold coordinates */
// 0-based
static inline int sd_coordinate_conversion(asm_dict_t *d, uint32 id, uint32 pos, uint32 *s, uint64 *p, int count_gap)
{
    *s = UINT32_MAX;
    *p = UINT64_MAX;

    if (id >= d->sdict->n)
        return SEQ_NOT_FOUND;

    uint32 i, n;
    uint64 *index = d->index;
    i = (uint32) (d->a[id] >> 32);
    n = (uint32) (d->a[id]) + i;

    while (i < n && index[i]>>32 <= pos)
        ++i;
    if (i == n)
        return POS_NOT_IN_RANGE;

    sd_seg_t *seg = &d->seg[(uint32) index[i]];
    if (pos < seg->x || pos >= seg->x + seg->y)
        return POS_NOT_IN_RANGE;

    uint64 offset = count_gap? seg->b : seg->a;
    *s = seg->s;
    *p = seg->c & 1? offset + seg->x + seg->y - 1 - pos : offset + pos - seg->x;

    return CC_SUCCESS;
}

/* scaffold coordinates to contig coordinates */
// 0-based
static inline int sd_coordinate_rev_conversion(asm_dict_t *d, uint32 id, uint64 pos, uint32 *s, uint32 *p, int count_gap)
{
    *s = UINT32_MAX;
    *p = UINT32_MAX;

    if (id >= d->n)
        return SEQ_NOT_FOUND;

    sd_seg_t *seg = d->seg + d->s[id].s;
    uint32 i = 0, n = d->s[id].n;
    uint64 offset = 0;

    if (count_gap) {
        while (i < n && pos >= seg[i].b)
            ++i;
        seg = &seg[i-1];
        offset = pos - seg->b;
    } else {
        while (i < n && pos >= seg[i].a)
            ++i;
        seg = &seg[i-1];
        offset = pos - seg->a;
    }
    
    // check if in gap or exceeds seq length
    if (offset >= seg->y)
        return POS_NOT_IN_RANGE;

    *s = seg->c >> 1;
    *p = seg->c & 1? seg->y - 1 - offset : seg->x + offset;

    return CC_SUCCESS;
}

#endif /* SDICT_H_ */

