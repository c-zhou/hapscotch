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
#include <string.h>
#include <stdint.h>

#include "kvec.h"
#include "paf.h"
#include "sdict.h"
#include "misc.h"

#include "alnio.h"

void validate_break_agp(asm_dict_t *bd)
{
    sdict_t *sdict;
    sd_seg_t *seg;
    sd_aseq_t *s;
    uint32 i, j, k, n, pos;

    sdict = bd->sdict;
    for (i = 0; i < bd->n; i++) {
        s = &bd->s[i];
        if (s->n != 1 || s->gn != 0) {
            fprintf(stderr, "[E::%s] AGP object %s is not a single-component sequence slice\n", __func__, s->name);
            exit(EXIT_FAILURE);
        }
        seg = &bd->seg[s->s];
        if (seg->c & 1) {
            fprintf(stderr, "[E::%s] AGP object %s is not in '+' orientation\n", __func__, s->name);
            exit(EXIT_FAILURE);
        }
    }
    // components of each input sequence must tile [0, len)
    for (i = 0; i < sdict->n; i++) {
        k = (uint32) (bd->a[i] >> 32);
        n = (uint32) bd->a[i];
        pos = 0;
        for (j = 0; j < n; j++) {
            // bd->index[k..k+n) sorted by component end position
            seg = &bd->seg[(uint32) bd->index[k + j]];
            if (seg->x != pos) break;
            pos = seg->x + seg->y;
        }
        if (n == 0 || j < n || pos != sdict->s[i].len) {
            fprintf(stderr, "[E::%s] AGP components do not fully cover sequence %s\n", __func__, sdict->s[i].name);
            fprintf(stderr, "[E::%s] the AGP file (-a) must cover every base of every input sequence exactly once\n", __func__);
            exit(EXIT_FAILURE);
        }
    }
}

sdict_t *make_piece_sdict(asm_dict_t *bd)
{
    sdict_t *d;
    uint32 i;

    d = sd_init();
    for (i = 0; i < bd->n; i++)
        sd_put(d, bd->s[i].name, (uint32) bd->s[i].len);
    if (d->n != bd->n) {
        fprintf(stderr, "[E::%s] duplicate object names in the AGP file\n", __func__);
        exit(EXIT_FAILURE);
    }
    return d;
}

// find the AGP component segment covering position pos on (raw) sequence id
static inline sd_seg_t *break_agp_seg(asm_dict_t *bd, uint32 id, int64 pos)
{
    uint32 i, n;
    i = (uint32) (bd->a[id] >> 32);
    n = (uint32) bd->a[id] + i;
    while (i < n && (int64) (bd->index[i] >> 32) <= pos)
        ++i;
    return i < n? &bd->seg[(uint32) bd->index[i]] : NULL;
}

static int i64_val_acmpfunc(const void *a, const void *b)
{
    int64 x = *(const int64 *) a, y = *(const int64 *) b;
    return (x > y) - (x < y);
}

// remap an alignment from raw sequence coordinates to piece coordinates
// splitting it at piece boundaries on both axes
// partner coordinates at cut points are estimated by linear interpolation
static void remap_split_aln(asm_dict_t *bd, uint32 qid, uint32 tid, paf_rec_t *rec, aln_vec_t *alns)
{
    kvec_t(int64) cuts;
    sd_seg_t *seg;
    int64 ql, tl, e, q, q0, q1, t0, t1, tm, x, y;
    int64 i, k, n;
    uint32 ps, pt, mlen;

    ql = (int64) rec->qe - rec->qs;
    tl = (int64) rec->te - rec->ts;
    if (ql <= 0 || tl <= 0) return;

    kv_init(cuts);
    // piece boundaries on the query axis
    i = (uint32) (bd->a[qid] >> 32);
    n = (uint32) bd->a[qid] + i;
    for (; i < n; i++) {
        e = bd->index[i] >> 32;
        if (e > rec->qs && e < rec->qe)
            kv_push(int64, cuts, e);
    }
    // piece boundaries on the target axis projected onto the query axis
    i = (uint32) (bd->a[tid] >> 32);
    n = (uint32) bd->a[tid] + i;
    for (; i < n; i++) {
        e = bd->index[i] >> 32;
        if (e > rec->ts && e < rec->te) {
            q = rec->rev? rec->qs + ((int64) rec->te - e) * ql / tl :
                          rec->qs + (e - (int64) rec->ts) * ql / tl;
            if (q > rec->qs && q < rec->qe)
                kv_push(int64, cuts, q);
        }
    }
    kv_push(int64, cuts, rec->qe);
    qsort(cuts.a, cuts.n, sizeof(int64), i64_val_acmpfunc);

    q0 = rec->qs;
    for (k = 0; k < cuts.n; k++) {
        q1 = cuts.a[k];
        if (q1 <= q0) continue;
        // project [q0, q1) onto the target axis
        if (!rec->rev) {
            t0 = rec->ts + (q0 - rec->qs) * tl / ql;
            t1 = rec->ts + (q1 - rec->qs) * tl / ql;
        } else {
            t0 = rec->ts + ((int64) rec->qe - q1) * tl / ql;
            t1 = rec->ts + ((int64) rec->qe - q0) * tl / ql;
        }
        if (t1 <= t0) { q0 = q1; continue; }
        // clip the target interval into the piece covering its midpoint
        // interpolation rounding may spill a few bases across a boundary
        tm = (t0 + t1 - 1) / 2;
        seg = break_agp_seg(bd, tid, tm);
        x = seg->x;
        y = seg->x + seg->y;
        if (t0 < x) t0 = x;
        if (t1 > y) t1 = y;
        if (t1 <= t0) { q0 = q1; continue; }
        pt = seg->s;
        seg = break_agp_seg(bd, qid, q0);
        ps = seg->s;
        if (ps != pt) { // skip self-mappings within a piece
            mlen = (uint32) ((int64) rec->ml * (q1 - q0) / ql);
            if (mlen == 0) mlen = 1;
            kv_push(aln_t, *alns, ((aln_t){ps, 1, pt, rec->rev,
                (int) (q0 - seg->x), (int) (q1 - seg->x),
                (int) (t0 - x), (int) (t1 - x), mlen, NULL}));
        }
        q0 = q1;
    }
    kv_destroy(cuts);
}

aln_t *read_pafs(char **fs, int fn, sdict_t *dicts, asm_dict_t *break_dict, int dual_aln, int64 *_naln)
{
    paf_file_t *paf;
    paf_rec_t _rec, *rec;
    sdict_t *vd;
    uint32 qid, tid;
    uint64 n_rec;
    aln_vec_t alns;
    int i;

    if (_naln) *_naln = 0;

    // input PAF coordinates are always on the raw (uncorrected) sequences
    vd = break_dict? break_dict->sdict : dicts;

    kv_init(alns);
    kv_resize(aln_t, alns, 1<<24);
    rec = &_rec;
    n_rec = 0;
    for (i = 0; i < fn; i++) {
        paf = paf_open(fs[i]);
        if (!paf) {
            fprintf(stderr, "[E::%s] cannot open paf file to read: %s\n", __func__, fs[i]);
            exit(1);
        }
        while (paf_read(paf, rec) >= 0) {
            if (++n_rec % 1000000 == 0)
                fprintf(stderr, "[M::%s] read %lld paf records\n", __func__, n_rec);
            // skip those not in sequence dicts
            if (!sd_exists(vd, rec->qn, rec->ql, &qid) ||
                !sd_exists(vd, rec->tn, rec->tl, &tid))
                continue;
            // skip self-alignments - with a break AGP they are kept
            // as they may link different pieces of the same sequence
            if (qid == tid && (!break_dict || rec->qs == rec->ts))
                continue;
            if (dual_aln && (qid > tid || (qid == tid && rec->qs > rec->ts)))
                continue; // skip dual alignment
            // specifically, remove sequence overhangs due to assemblers
            // i.e., an identical sequence overlap
            if (rec->ml == rec->bl && 
                ((rec->qs > 0 && rec->qe == rec->ql) || (rec->qs == 0 && rec->qe < rec->ql)) &&
                ((rec->ts > 0 && rec->te == rec->tl) || (rec->ts == 0 && rec->te < rec->tl)))
                continue;
            if (break_dict)
                remap_split_aln(break_dict, qid, tid, rec, &alns);
            else
                kv_push(aln_t, alns, ((aln_t){qid, 1, tid, rec->rev, rec->qs, rec->qe, rec->ts, rec->te, rec->ml, NULL}));
        }
        paf_close(paf);
    }
    fprintf(stderr, "[M::%s] read %lld paf records, %ld retained\n", __func__, n_rec, alns.n);

    MYREALLOC(alns.a, alns.n);

    *_naln = alns.n;
    return (alns.a);
}
