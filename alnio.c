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
#include <unistd.h>

#include "kvec.h"
#include "paf.h"
#include "sdict.h"
#include "misc.h"
#include "ONElib.h"
#include "version.h"

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
// self_aln: 1 keep fragments within the same piece; 0 keep those between pieces; -1 keep all
static void remap_split_aln(asm_dict_t *bd, uint32 qid, uint32 tid, paf_rec_t *rec, int self_aln, aln_vec_t *alns)
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

    if (self_aln > 0) self_aln = 1;
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
        if (self_aln < 0 || (ps == pt) == self_aln) {
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

// filter an alignment record and add it to alns, remapping onto break-AGP pieces if break_dict is given
// qid and tid are the (raw) sequence ids of the query and target
static void add_aln_rec(asm_dict_t *break_dict, uint32 qid, uint32 tid, paf_rec_t *rec, int dual_aln, int self_aln, aln_vec_t *alns)
{
    if (!break_dict && (self_aln? qid != tid : qid == tid))
        return;
    if (!self_aln && dual_aln && (qid > tid || (qid == tid && rec->qs > rec->ts)))
        return;
    // specifically, remove sequence overhangs due to assemblers
    // i.e., an identical sequence overlap
    if (rec->ml == rec->bl && 
        ((rec->qs > 0 && rec->qe == rec->ql) || (rec->qs == 0 && rec->qe < rec->ql)) &&
        ((rec->ts > 0 && rec->te == rec->tl) || (rec->ts == 0 && rec->te < rec->tl)))
        return;
    if (break_dict)
        remap_split_aln(break_dict, qid, tid, rec, self_aln, alns);
    else
        kv_push(aln_t, *alns, ((aln_t){qid, 1, tid, rec->rev, rec->qs, rec->qe, rec->ts, rec->te, rec->ml, NULL}));
}

aln_t *read_pafs(char **fs, int fn, sdict_t *dicts, asm_dict_t *break_dict, int dual_aln, int self_aln, int64 *_naln)
{
    paf_file_t *paf;
    paf_rec_t _rec, *rec;
    sdict_t *vd;
    uint32 qid, tid;
    uint64 n_rec;
    aln_vec_t alns;
    int i;

    if (_naln) *_naln = 0;
    if (!fn) return NULL;

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
            add_aln_rec(break_dict, qid, tid, rec, dual_aln, self_aln, &alns);
        }
        paf_close(paf);
    }
    fprintf(stderr, "[M::%s] read %lld paf records, %ld retained\n", __func__, n_rec, alns.n);

    MYREALLOC(alns.a, alns.n);

    *_naln = alns.n;
    return (alns.a);
}

// contig of a GDB skeleton embedded in a .1aln file
typedef struct {
    uint32 sid;  // id of its scaffold in the sequence dict; UINT32_MAX if not in the dict
    uint32 slen; // scaffold length
    uint32 clen; // contig length
    int64  sbeg; // contig start on the scaffold
} gdb_ctg_t;

typedef struct {
    int64 n, m;     // contigs
    gdb_ctg_t *a;
    int64 ns, ms;   // scaffolds
    char **names;
    int64 *slens;
    int64 *fctgs;   // first contig of each scaffold
} gdb_skel_t;

static void gdb_skel_destroy(gdb_skel_t *sk)
{
    int64 i;
    for (i = 0; i < sk->ns; i++)
        free(sk->names[i]);
    free(sk->names);
    free(sk->slens);
    free(sk->fctgs);
    free(sk->a);
}

// read the GDB skeleton(s) at the head of a .1aln file
// on return, the file is positioned at the first alignment ('A' line) if any
// the trace point spacing is returned in tspace if not NULL
// return the number of skeletons read
static int read_gdb_skeletons(OneFile *of, sdict_t *vd, gdb_skel_t *skel, int *tspace)
{
    gdb_skel_t *sk;
    char *name;
    int64 i, j, e, len;
    uint32 sid;
    int k, nsk;

    if (tspace) *tspace = 0;
    sk = NULL;
    nsk = 0;
    while (oneReadLine(of)) {
        if (of->lineType == 'A')
            break;
        if (of->lineType == 't') {
            if (tspace) *tspace = (int) oneInt(of, 0);
            continue;
        }
        if (of->lineType == 'g') {
            if (nsk == 2) {
                fprintf(stderr, "[E::%s] more than two GDB skeletons in the 1aln file\n", __func__);
                exit(EXIT_FAILURE);
            }
            sk = &skel[nsk++];
            continue;
        }
        if (!sk) continue; // header lines before any skeleton
        switch (of->lineType) {
            case 'S':
                if (sk->ns == sk->ms) {
                    MYEXPAND(sk->names, sk->ms);
                    MYREALLOC(sk->slens, sk->ms);
                    MYREALLOC(sk->fctgs, sk->ms);
                }
                // scaffold name, truncated at the first white space
                len = oneLen(of);
                MYMALLOC(name, len + 1);
                memcpy(name, oneString(of), len);
                name[len] = '\0';
                for (i = 0; i < len && name[i] != ' ' && name[i] != '\t'; i++);
                name[i] = '\0';
                sk->names[sk->ns] = name;
                sk->slens[sk->ns] = 0;
                sk->fctgs[sk->ns] = sk->n;
                sk->ns++;
                break;
            case 'G':
                sk->slens[sk->ns - 1] += oneInt(of, 0);
                break;
            case 'C':
                if (sk->n == sk->m) MYEXPAND(sk->a, sk->m);
                len = oneInt(of, 0);
                sk->a[sk->n++] = (gdb_ctg_t) {UINT32_MAX, 0, (uint32) len, sk->slens[sk->ns - 1]};
                sk->slens[sk->ns - 1] += len;
                break;
            default:
                break;
        }
    }

    // map scaffolds onto the sequence dict, as for PAF records
    for (k = 0; k < nsk; k++) {
        sk = &skel[k];
        for (i = 0; i < sk->ns; i++) {
            if (!sd_exists(vd, sk->names[i], (uint32) sk->slens[i], &sid))
                sid = UINT32_MAX;
            e = i + 1 < sk->ns? sk->fctgs[i + 1] : sk->n;
            for (j = sk->fctgs[i]; j < e; j++) {
                sk->a[j].sid  = sid;
                sk->a[j].slen = (uint32) sk->slens[i];
            }
        }
    }

    return nsk;
}

// read the current alignment record of a .1aln file into rec in scaffold (PAF) coordinates
// the match and block lengths are estimated from the number of differences as in FastGA ALNtoPAF
// on return, the file is positioned at the next alignment; *more is set to 0 at the end of the file
// return 1 if both sequences are in the sequence dict, 0 otherwise
static int read_1aln_rec(OneFile *of, gdb_skel_t *ask, gdb_skel_t *bsk, paf_rec_t *rec, uint32 *qid, uint32 *tid, int *_diffs, int *more)
{
    gdb_ctg_t *actg, *bctg;
    int64 abpos, aepos, bbpos, bepos, blk;
    uint32 aread, bread;
    int diffs;

    aread = (uint32) oneInt(of, 0);
    abpos = oneInt(of, 1);
    aepos = oneInt(of, 2);
    bread = (uint32) oneInt(of, 3);
    bbpos = oneInt(of, 4);
    bepos = oneInt(of, 5);
    rec->rev = 0;
    diffs = 0;
    *more = 0;
    // skip D/T/X/... lines to the next alignment, noting R and D
    while (oneReadLine(of)) {
        if (of->lineType == 'A') {
            *more = 1;
            break;
        } else if (of->lineType == 'R') {
            rec->rev = 1;
        } else if (of->lineType == 'D') {
            diffs = (int) oneInt(of, 0);
        }
    }
    if (_diffs) *_diffs = diffs;

    if (aread >= ask->n || bread >= bsk->n) {
        fprintf(stderr, "[E::%s] contig id out of the GDB skeleton range\n", __func__);
        exit(1);
    }
    actg = &ask->a[aread];
    bctg = &bsk->a[bread];
    *qid = actg->sid;
    *tid = bctg->sid;
    if (actg->sid == UINT32_MAX || bctg->sid == UINT32_MAX)
        return 0;

    blk = (aepos - abpos) + (bepos - bbpos);
    rec->ql = actg->slen;
    rec->qs = (uint32) (actg->sbeg + abpos);
    rec->qe = (uint32) (actg->sbeg + aepos);
    rec->tl = bctg->slen;
    if (rec->rev) {
        // b coordinates are on the reverse complement of the contig
        rec->ts = (uint32) (bctg->sbeg + bctg->clen - bepos);
        rec->te = (uint32) (bctg->sbeg + bctg->clen - bbpos);
    } else {
        rec->ts = (uint32) (bctg->sbeg + bbpos);
        rec->te = (uint32) (bctg->sbeg + bepos);
    }
    rec->ml = (uint32) ((blk - diffs) / 2);
    rec->bl = (uint32) (blk / 2);

    return 1;
}

// open a .1aln file and read its GDB skeleton(s); exit on failure
// on return, *ask and *bsk point to the skeletons of the a and b sequences
static OneFile *open_1aln(char *f, sdict_t *vd, gdb_skel_t *skel, gdb_skel_t **ask, gdb_skel_t **bsk, int *_nsk, int *tspace)
{
    OneFile *of;
    int nsk;

    // binary .1aln files carry their own schema
    of = oneFileOpenRead(f, NULL, "aln", 1);
    if (!of) {
        fprintf(stderr, "[E::%s] cannot open 1aln file to read: %s\n", __func__, f);
        exit(1);
    }
    MYBZERO(skel, 2);
    nsk = read_gdb_skeletons(of, vd, skel, tspace);
    if (nsk == 0) {
        fprintf(stderr, "[E::%s] no GDB skeleton found in 1aln file: %s\n", __func__, f);
        exit(1);
    }
    *ask = &skel[0];
    *bsk = nsk == 2? &skel[1] : &skel[0]; // one skeleton for a self comparison
    *_nsk = nsk;
    return of;
}

// read .1aln file(s) into an aln_t array, remapping onto break-AGP pieces if break_dict is given
// the GDB skeleton(s) must be included in the .1aln file; the sequence files are never touched
// coordinates are converted from contig to scaffold (PAF) coordinates
aln_t *read_1aln(char **fs, int fn, sdict_t *dicts, asm_dict_t *break_dict, int dual_aln, int self_aln, int64 *_naln)
{
    OneFile *of;
    gdb_skel_t skel[2], *ask, *bsk;
    paf_rec_t _rec, *rec;
    sdict_t *vd;
    uint64 n_rec;
    aln_vec_t alns;
    uint32 qid, tid;
    int i, nsk, more;

    if (_naln) *_naln = 0;
    if (!fn) return NULL;

    // input coordinates are always on the raw (uncorrected) sequences
    vd = break_dict? break_dict->sdict : dicts;

    kv_init(alns);
    kv_resize(aln_t, alns, 1<<24);
    rec = &_rec;
    MYBZERO(rec, 1);
    n_rec = 0;
    for (i = 0; i < fn; i++) {
        of = open_1aln(fs[i], vd, skel, &ask, &bsk, &nsk, NULL);
        more = (of->lineType == 'A');
        while (more) {
            if (++n_rec % 1000000 == 0)
                fprintf(stderr, "[M::%s] read %lld 1aln records\n", __func__, n_rec);
            // skip those not in sequence dicts
            if (!read_1aln_rec(of, ask, bsk, rec, &qid, &tid, NULL, &more))
                continue;
            add_aln_rec(break_dict, qid, tid, rec, dual_aln, self_aln, &alns);
        }
        gdb_skel_destroy(&skel[0]);
        if (nsk == 2) gdb_skel_destroy(&skel[1]);
        oneFileClose(of);
    }
    fprintf(stderr, "[M::%s] read %lld 1aln records, %ld retained\n", __func__, n_rec, alns.n);

    MYREALLOC(alns.a, alns.n);

    *_naln = alns.n;
    return (alns.a);
}

// a .1aln file always starts with the plain text header "1 3 aln"
// gzipped PAF or stdin ("-") are treated as PAF
static int is_1aln_file(const char *fn)
{
    char buf[8];
    FILE *fp;
    int ret;

    if (!fn || !strcmp(fn, "-"))
        return 0;
    fp = fopen(fn, "r");
    if (!fp) {
        fprintf(stderr, "[E::%s] cannot open alignment file to read: %s\n", __func__, fn);
        exit(1);
    }
    ret = fread(buf, 1, 7, fp) == 7 && !strncmp(buf, "1 3 aln", 7);
    fclose(fp);
    return ret;
}

aln_t *read_alns(char **fs, int fn, sdict_t *dicts, asm_dict_t *break_dict, int dual_aln, int self_aln, int64 *_naln)
{
    char **pfs, **afs;
    aln_t *palns, *aalns;
    int64 npaln, naaln;
    int i, npf, naf;

    if (_naln) *_naln = 0;
    if (!fn) return NULL;

    // split input files by format
    MYMALLOC(pfs, fn);
    MYMALLOC(afs, fn);
    npf = naf = 0;
    for (i = 0; i < fn; i++) {
        if (is_1aln_file(fs[i]))
            afs[naf++] = fs[i];
        else
            pfs[npf++] = fs[i];
    }

    palns = aalns = NULL;
    npaln = naaln = 0;
    if (npf) palns = read_pafs(pfs, npf, dicts, break_dict, dual_aln, self_aln, &npaln);
    if (naf) aalns = read_1aln(afs, naf, dicts, break_dict, dual_aln, self_aln, &naaln);
    free(pfs);
    free(afs);

    if (!naf) {
        *_naln = npaln;
        return palns;
    }
    if (!npf) {
        *_naln = naaln;
        return aalns;
    }

    // mixed input formats
    MYREALLOC(palns, npaln + naaln);
    if (naaln) memcpy(palns + npaln, aalns, sizeof(aln_t) * naaln);
    free(aalns);
    *_naln = npaln + naaln;
    return palns;
}

// schema of .1aln files, as in FastGA alncode.c
static char *ALN_SCHEMA_TEXT =
    "1 3 def 2 1                 schema for aln and FastGA\n"
    ".\n"
    "P 3 seq                     SEQUENCE\n"
    "O s 2 3 INT 6 STRING        length and id for group of sequences = a scaffold\n"
    "G S                         scaffolds (s) group sequence objects (S)\n"
    "D n 2 4 CHAR 3 INT          non-acgt chars outside (between) sequences within scaffold\n"
    "O S 1 3 DNA                 sequence\n"
    "D I 1 6 STRING              identifier of sequence\n"
    ".\n"
    "P 3 aln                     ALIGNMENTS\n"
    "D t 1 3 INT                 trace point spacing in a - global\n"
    ".                           GDB skeleton (may not be presend)\n"
    "O g 0                       groups scaffolds into a GDB skeleton\n"
    "G S                         collection of scaffolds constituting a GDB\n"
    "O S 1 6 STRING              id for a scaffold\n"
    "D G 1 3 INT                 gap of given length\n"
    "D C 1 3 INT                 contig of given length\n"
    ".\n"
    "O a 0                       groups A's into a colinear chain\n"
    "G A                         chains (a) group alignment objects (A)\n"
    "D p 2 3 INT 3 INT           spacing in a,b between end of previous alignment and start of next\n"
    ".                           alignment: a_read[beg..end] b_read[beg..end], 0-indexing\n"
    "O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT\n"
    "D L 2 3 INT 3 INT           lengths of sequences a and b\n"
    "D R 0                       flag: reverse-complement sequence b\n"
    "D D 1 3 INT                 differences: number of diffs = substitions + indels\n"
    "D T 1 8 INT_LIST            trace points in b\n"
    "D X 1 8 INT_LIST            number of differences in alignment per trace interval\n"
    "D Q 1 3 INT                 quality: alignment confidence in phred units (currently unused)\n"
    "D E 1 3 INT                 match: number of equal bases (currently unused)\n"
    "D Z 1 6 STRING              cigar string: encodes precise alignment (currently unused)\n"
    "D U 1 3 INT                 putative unit size of a TR alignment (FASTAN)\n"
;

// default trace point spacing of FastGA, used for PAF input
#define DEFAULT_TSPACE 100

// open the output .1aln file on fo and write the header and the GDB skeleton of the pieces
// pieces are named as in the AGP, or as "name_start_end" (1-based, inclusive) on the original sequence if raw_name is set
static OneFile *open_split_out(FILE *fo, OneSchema *schema, asm_dict_t *break_dict, int raw_name, char *fin, char *ref, int tspace)
{
    OneFile *ofo;
    FILE *fp;
    sd_seg_t *seg;
    char name[BUFF_SIZE];
    uint32 i;
    int is_bin, fd;

    // the output stream is handed to ONElib as a duplicate so that fo stays with the caller
    // binary ONE files record absolute file offsets, hence need a seekable output
    fflush(fo);
    is_bin = lseek(fileno(fo), 0, SEEK_CUR) >= 0;
    fd = dup(fileno(fo));
    fp = fd < 0? NULL : fdopen(fd, "w");
    if (!fp) {
        fprintf(stderr, "[E::%s] cannot write the output 1aln file\n", __func__);
        exit(1);
    }
    if (!is_bin)
        fprintf(stderr, "[W::%s] output is not seekable, writing 1aln file in ASCII ONE format\n", __func__);

    ofo = schema? oneFileOpenWriteNew("-", schema, "aln", is_bin, 1) : NULL;
    if (!ofo) {
        fprintf(stderr, "[E::%s] cannot create the output 1aln file\n", __func__);
        exit(1);
    }
    ofo->f = fp; // ONElib opened "-" as stdout
    oneAddProvenance(ofo, "alnsplit", HAPSCOTCH_VERSION, "split_aln %s", fin);
    oneAddReference(ofo, ref? ref : fin, 1);

    // trace point spacing
    oneInt(ofo, 0) = tspace;
    oneWriteLine(ofo, 't', 0, 0);

    // GDB skeleton of the sequence pieces, one contig per piece
    oneWriteLine(ofo, 'g', 0, 0);
    for (i = 0; i < break_dict->n; i++) {
        if (raw_name) {
            // each piece is a single component of an original sequence
            seg = &break_dict->seg[break_dict->s[i].s];
            snprintf(name, BUFF_SIZE, "%s_%u_%u", break_dict->sdict->s[seg->c >> 1].name, seg->x + 1, seg->x + seg->y);
        } else {
            snprintf(name, BUFF_SIZE, "%s", break_dict->s[i].name);
        }
        oneWriteLine(ofo, 'S', strlen(name), name);
        oneInt(ofo, 0) = break_dict->s[i].len;
        oneWriteLine(ofo, 'C', 0, 0);
    }

    return ofo;
}

// remap an alignment record onto the pieces and write the fragments to ofo
// the number of differences is distributed over split fragments by fragment length
// return the number of alignments written
static int64 write_split_rec(OneFile *ofo, asm_dict_t *break_dict, uint32 qid, uint32 tid, paf_rec_t *rec, int diffs, aln_vec_t *frags)
{
    aln_t *f;
    int64 i, blk, flen;

    frags->n = 0;
    remap_split_aln(break_dict, qid, tid, rec, -1, frags);
    blk = ((int64) rec->qe - rec->qs) + ((int64) rec->te - rec->ts);
    for (i = 0; i < (int64) frags->n; i++) {
        f = &frags->a[i];
        oneInt(ofo, 0) = f->aread;
        oneInt(ofo, 1) = f->abpos;
        oneInt(ofo, 2) = f->aepos;
        oneInt(ofo, 3) = f->bread;
        if (f->rev) {
            // b coordinates on the reverse complement of the piece
            flen = break_dict->s[f->bread].len;
            oneInt(ofo, 4) = flen - f->bepos;
            oneInt(ofo, 5) = flen - f->bbpos;
        } else {
            oneInt(ofo, 4) = f->bbpos;
            oneInt(ofo, 5) = f->bepos;
        }
        oneWriteLine(ofo, 'A', 0, 0);
        if (f->rev)
            oneWriteLine(ofo, 'R', 0, 0);
        oneInt(ofo, 0) = frags->n == 1? diffs :
            ((int64) diffs * ((f->aepos - f->abpos) + (f->bepos - f->bbpos)) + blk / 2) / blk;
        oneWriteLine(ofo, 'D', 0, 0);
        // empty trace points
        oneWriteLine(ofo, 'T', 0, 0);
        oneWriteLine(ofo, 'X', 0, 0);
    }

    return (int64) frags->n;
}

// remap an alignment file (.1aln or PAF, detected by file header) onto the sequence pieces of a break AGP
// and write a new .1aln file to fo
// the output contains a single GDB skeleton made of the pieces (one contig per piece)
// and empty trace point lines; alignments are split at piece boundaries as in read_pafs
// for PAF input, the number of differences is estimated as (alen + blen - 2 * matches)
// and the trace point spacing is set to the FastGA default
// pieces are named as in the AGP, or as "name_start_end" (1-based, inclusive) on the original sequence if raw_name is set
// fo is written in binary ONE format if it is seekable (a regular file), in ASCII ONE otherwise
// return the number of alignments written
int64 split_aln(char *fin, asm_dict_t *break_dict, int raw_name, FILE *fo)
{
    OneFile *of, *ofo;
    OneSchema *schema;
    paf_file_t *paf;
    gdb_skel_t skel[2], *ask, *bsk;
    paf_rec_t _rec, *rec;
    aln_vec_t frags;
    sdict_t *vd;
    char *ref;
    uint64 n_rec;
    int64 n_out, blk;
    uint32 qid, tid;
    int k, nsk, more, diffs, tspace, is_aln;

    if (!break_dict) {
        fprintf(stderr, "[E::%s] a break AGP dictionary is required\n", __func__);
        exit(1);
    }

    // input coordinates are always on the raw (uncorrected) sequences
    vd = break_dict->sdict;
    is_aln = is_1aln_file(fin);

    of = NULL;
    paf = NULL;
    ref = NULL;
    nsk = 0;
    if (is_aln) {
        of = open_1aln(fin, vd, skel, &ask, &bsk, &nsk, &tspace);
        if (tspace <= 0) {
            fprintf(stderr, "[E::%s] no trace point spacing ('t' line) found in 1aln file: %s\n", __func__, fin);
            exit(1);
        }
        // keep the source reference of the first sequence set
        if (of->info['<'])
            for (k = 0; k < of->info['<']->accum.count; k++)
                if (of->reference[k].count == 1)
                    ref = of->reference[k].filename;
    } else {
        paf = paf_open(fin);
        if (!paf) {
            fprintf(stderr, "[E::%s] cannot open paf file to read: %s\n", __func__, fin);
            exit(1);
        }
        tspace = DEFAULT_TSPACE;
    }

    schema = oneSchemaCreateFromText(ALN_SCHEMA_TEXT);
    ofo = open_split_out(fo, schema, break_dict, raw_name, fin, ref, tspace);

    // alignments
    kv_init(frags);
    rec = &_rec;
    MYBZERO(rec, 1);
    n_rec = 0;
    n_out = 0;
    if (is_aln) {
        more = (of->lineType == 'A');
        while (more) {
            if (++n_rec % 1000000 == 0)
                fprintf(stderr, "[M::%s] read %lld 1aln records\n", __func__, n_rec);
            // skip those not in sequence dicts
            if (!read_1aln_rec(of, ask, bsk, rec, &qid, &tid, &diffs, &more))
                continue;
            n_out += write_split_rec(ofo, break_dict, qid, tid, rec, diffs, &frags);
        }
    } else {
        while (paf_read(paf, rec) >= 0) {
            if (++n_rec % 1000000 == 0)
                fprintf(stderr, "[M::%s] read %lld paf records\n", __func__, n_rec);
            // skip those not in sequence dicts
            if (!sd_exists(vd, rec->qn, rec->ql, &qid) ||
                !sd_exists(vd, rec->tn, rec->tl, &tid))
                continue;
            // estimated so that the match length is recovered by read_1aln as (alen + blen - diffs) / 2
            blk = ((int64) rec->qe - rec->qs) + ((int64) rec->te - rec->ts);
            diffs = blk > 2 * (int64) rec->ml? (int) (blk - 2 * (int64) rec->ml) : 0;
            n_out += write_split_rec(ofo, break_dict, qid, tid, rec, diffs, &frags);
        }
    }
    fprintf(stderr, "[M::%s] read %lld %s records, %lld written\n", __func__, n_rec, is_aln? "1aln" : "paf", n_out);

    kv_destroy(frags);
    oneFileClose(ofo); // also closes the duplicate of fo
    oneSchemaDestroy(schema);
    if (is_aln) {
        gdb_skel_destroy(&skel[0]);
        if (nsk == 2) gdb_skel_destroy(&skel[1]);
        oneFileClose(of);
    } else {
        paf_close(paf);
    }

    return n_out;
}
