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

/********************************** Revision History *****************************
 *                                                                               *
 * 16/06/26 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ketopt.h"
#include "kvec.h"
#include "kstring.h"

#include "agp-spec.h"
#include "sdict.h"
#include "misc.h"
#include "version.h"

int VERBOSE = 0;
int LONG_HELP = 0;

static AGP_CT_t LG_AGP_SEQ_COMPONENT_TYPE = AGP_CT_W;
static AGP_CT_t LG_AGP_GAP_COMPONENT_TYPE = AGP_CT_N;
static AGP_LE_t LG_AGP_LINKAGE_EVIDENCE = AGP_LE_ALIGN_GENUS;
static int LG_AGP_GAP_SIZE = DEFAULT_AGP_U_GAP_SIZE;

static double at_realtime0;

static inline void write_agp_seq(FILE *fo, char *s_name, uint64 s_beg, uint64 s_end, uint32 b,
    char *c_name, uint32 c_beg, uint32 c_end, uint32 c_oris)
{
    fprintf(fo, "%s\t%llu\t%llu\t%u\t%s\t%s\t%u\t%u\t%s\n", s_name, s_beg, s_end, b,
        agp_component_type_val(LG_AGP_SEQ_COMPONENT_TYPE), c_name, c_beg, c_end, 
        agp_orientation_val(c_oris? AGP_OT_MINUS : AGP_OT_PLUS));
}

static inline void write_agp_gap(FILE *fo, char *s_name, uint64 s_beg, uint64 s_end, uint32 b)
{
    fprintf(fo, "%s\t%llu\t%llu\t%u\t%s\t%d\t%s\t%s\t%s\n", s_name, s_beg, s_end, b,
        agp_component_type_val(LG_AGP_GAP_COMPONENT_TYPE), 
        LG_AGP_GAP_SIZE, 
        agp_gap_type_val(AGP_GT_SCAFFOLD),
        agp_linkage_val(AGP_LG_YES), 
        agp_linkage_evidence_val(LG_AGP_LINKAGE_EVIDENCE));
}

typedef struct {
    char *name;
    uint32 len:31, rev:1;
    int g, h;
    int64 pos;
} seq_t;

typedef struct {
    char *name;
    int ngrp;
    uint32 *grps;
} scf_t;

static int seq_cmpfunc(const void *a, const void *b)
{
    const seq_t *x = a, *y = b;
    if (x->g != y->g)
        return (x->g > y->g) - (x->g < y->g);
    if (x->h != y->h)
        return (x->h > y->h) - (x->h < y->h);
    return (x->pos > y->pos) - (x->pos < y->pos);
}

static inline int parse_line(char *s, char **f, int m)
{
    int n = 0;
    while (*s) {
        while (isspace(*s)) s++;
        if (!*s) break;
        if (n < m) f[n++] = s;
        while (*s && !isspace(*s)) s++;
        if (*s) *s++ = '\0';
    }
    return n;
}

static void annot_haps_from_agp(char *agp, char *ann, int hap, FILE *fo)
{
    kvec_t(seq_t)  vseqs;
    kvec_t(scf_t)  vscfs;
    kvec_t(uint32) vgrps;
    kstring_t _ks = {0, 0, NULL}, *ks = &_ks;
    FILE    *fp;
    seq_t   *seq, *seqs;
    scf_t   *scf, *scfs;
    char    *line, *dash, *f[10];
    size_t   lline;
    uint32   g, r;
    int64    part_no, scf_pos, seq_len, add_gap;
    int      nseq, nscf, max_hap, hap_lo, hap_hi;
    int      step, lo, hi, sbeg, send, abeg, aend;
    int      i, a, s, h, nf, mid;
    
    fp = fopen(ann, "r");
    if (!fp) {
        fprintf(stderr, "[E::%s] cannot open annotation file: %s\n", __func__, ann);
        exit(EXIT_FAILURE);
    }
    kv_init(vseqs);
    line  = NULL;
    lline = 0;
    while (getline(&line, &lline, fp) > 0) {
        // fields: [0] name  [1] len  [2] ori  [3] grp_id  [4] hap  [5] pos
        nf = parse_line(line, f, 10);
        if (nf < 6) {
            // reconstruct the line to show the invalid fields
            for (i = 0; i < nf; i++)
                line[strlen(f[i])] = '\t';
            fprintf(stderr, "[W::%s] invalid annotation line: %s\n", __func__, line);
            continue;
        }
        if (hap >= 0 && atoi(f[4]) != hap) continue;
        kv_pushp(seq_t, vseqs, &seq);
        seq->name = strdup(f[0]);
        seq->len  = atoi(f[1]);
        seq->rev  = (f[2][0] == '-') ? 1 : 0;
        seq->g    = atoi(f[3]);
        seq->h    = atoi(f[4]);
        seq->pos  = atoll(f[5]);
    }
    fclose(fp);

    seqs = vseqs.a;
    nseq = vseqs.n;
    qsort(seqs, nseq, sizeof(seq_t), seq_cmpfunc);

    fp = fopen(agp, "r");
    if (!fp) {
        fprintf(stderr, "[E::%s] cannot open AGP file: %s\n", __func__, agp);
        exit(EXIT_FAILURE);
    }
    ks_resize(ks, 64);
    kv_init(vscfs);
    kv_init(vgrps);
    while (getline(&line, &lline, fp) > 0) {
        // fields: [0] scfname [1] start [2] end [3] part [4] type
        //         [5] comp [6] comp_start [7] comp_end [8] ori
        if (line[0] == '#') continue;
        nf = parse_line(line, f, 10);
        if (nf < 9) {
            // reconstruct the line to show the invalid fields
            for (i = 0; i < nf; i++)
                line[strlen(f[i])] = '\t';
            fprintf(stderr, "[W::%s] invalid AGP line: %s\n", __func__, line);
            continue;
        }
        if (strcmp(f[4], "W") != 0) continue;
        if (f[5][0] != 's') continue;
        dash = strchr(f[5] + 1, '-');
        if (!dash) continue;
        *dash = '\0';
        g = atoi(f[5] + 1);
        *dash = '-';
        r = (f[8][0] == '-') ? 1 : 0;
        g = (g << 1) | r;
        if (ks->l && strcmp(f[0], ks->s) != 0) {
            kv_pushp(scf_t, vscfs, &scf);
            MYMALLOC(scf->grps, vgrps.n * sizeof(uint32));
            scf->name = strdup(ks->s);
            scf->ngrp = vgrps.n;
            memcpy(scf->grps, vgrps.a, vgrps.n * sizeof(uint32));
            vgrps.n = 0;
        }
        if (vgrps.n == 0) {
            kv_push(uint32, vgrps, g);
            ks->l = 0;
            ksprintf(ks, "%s", f[0]);
        } else if (g != vgrps.a[vgrps.n-1]) {
            kv_push(uint32, vgrps, g);
        }
    }
    // push the last scaffold
    if (vgrps.n) {
        kv_pushp(scf_t, vscfs, &scf);
        MYMALLOC(scf->grps, vgrps.n * sizeof(uint32));
        scf->name = strdup(ks->s);
        scf->ngrp = vgrps.n;
        memcpy(scf->grps, vgrps.a, vgrps.n * sizeof(uint32));
    }
    kv_destroy(vgrps);
    fclose(fp);

    scfs = vscfs.a;
    nscf = vscfs.n;

    max_hap = 0;
    for (i = 0; i < nseq; i++)
        if (seqs[i].h > max_hap) 
            max_hap = seqs[i].h;
    hap_lo = (hap >= 0) ? hap : 0;
    hap_hi = (hap >= 0) ? hap : max_hap;

    for (s = 0; s < nscf; s++) {
        scf = &scfs[s];
        for (h = hap_lo; h <= hap_hi; h++) {
            ks->l = 0;
            ksprintf(ks, "%s.h%d", scf->name, h);
            scf_pos = 1;
            part_no = 0;
            add_gap = 0;
            for (i = 0; i < scf->ngrp; i++) {
                g = scf->grps[i] >> 1;
                r = scf->grps[i] & 1;

                // lower bound
                lo = 0; hi = nseq;
                while (lo < hi) {
                    mid = (lo + hi) >> 1;
                    if (seqs[mid].g < g ||
                        (seqs[mid].g == g && seqs[mid].h < h))
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                sbeg = lo;

                // upper bound
                lo = sbeg; hi = nseq;
                while (lo < hi) {
                    mid = (lo + hi) >> 1;
                    if (seqs[mid].g < g ||
                        (seqs[mid].g == g && seqs[mid].h <= h))
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                send = lo;

                if (sbeg == send) continue;

                step = r ? -1 : 1;
                abeg = r ? send - 1 : sbeg;
                aend = r ? sbeg - 1 : send;

                for (a = abeg; a != aend; a += step) {
                    if (add_gap) {
                        write_agp_gap(fo, ks->s, scf_pos, scf_pos + LG_AGP_GAP_SIZE - 1, ++part_no);
                        scf_pos += LG_AGP_GAP_SIZE;
                    }
                    seq_len = seqs[a].len;
                    write_agp_seq(fo, ks->s, scf_pos, scf_pos + seq_len - 1, ++part_no, 
                        seqs[a].name, 1, seq_len, r ^ seqs[a].rev);
                    scf_pos += seq_len;
                    add_gap  = 1;
                }
            }
        }
    }

    for (i = 0; i < nseq; i++)
        free(seqs[i].name);
    free(seqs);
    for (i = 0; i < nscf; i++) {
        free(scfs[i].name);
        free(scfs[i].grps);
    }
    free(scfs);
    free(line);
    free(ks->s);
}

static void print_help_hap(FILE *fp_help)
{
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Usage: seqtools hap [options] <input.agp> <annot.txt>\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help, "    -p INT            print p-th haplotype only\n");
    fprintf(fp_help, "    -o STR            output to file [stdout]\n");
    if (LONG_HELP) {
        fprintf(fp_help, "\n");
        fprintf(fp_help, "    --seq-ctype  STR  AGP output sequence component type [%s]\n", agp_component_type_val(LG_AGP_SEQ_COMPONENT_TYPE));
        fprintf(fp_help, "    --gap-ctype  STR  AGP output gap component type [%s]\n", agp_component_type_val(LG_AGP_GAP_COMPONENT_TYPE));
        fprintf(fp_help, "    --gap-link   STR  AGP output gap linkage evidence [%s]\n", agp_linkage_evidence_val(LG_AGP_LINKAGE_EVIDENCE));
        fprintf(fp_help, "    --gap-size   INT  AGP output gap size between sequence component [%d]\n", LG_AGP_GAP_SIZE);
        fprintf(fp_help, "\n");
    }
    fprintf(fp_help, "    -?                print long help with extra option list\n");
    fprintf(fp_help, "    -h, --help        print this help\n");
    fprintf(fp_help, "    -V, --version     show version number\n");
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Example: seqtools hap -o haps.agp scaffolds.agp annot.txt\n");
    fprintf(fp_help, "         seqtools hap -p 1 -o hap1.agp scaffolds.agp annot.txt\n");
    fprintf(fp_help, "\n");
}

static ko_longopt_t hap_long_options[] = {
    { "seq-ctype",  ko_required_argument, 301 },
    { "gap-ctype",  ko_required_argument, 302 },
    { "gap-link",   ko_required_argument, 303 },
    { "gap-size",   ko_required_argument, 304 },
    { "help",       ko_no_argument,       'h' },
    { "version",    ko_no_argument,       'V' },
    { 0, 0, 0 }
};

int main_hap(int argc, char *argv[])
{
    if (argc < 2) {
        print_help_hap(stderr);
        return 1;
    }

    liftrlimit();
    at_realtime0 = realtime();

    FILE *fo;
    char *agp, *ann, *out;
    int hap;

    const char *opt_str = "o:p:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c;
    FILE *fp_help = stderr;
    agp = ann = out = 0;
    hap = -1;
    
    while ((c = ketopt(&opt, argc, argv, 1, opt_str, hap_long_options)) >= 0) {
        if (c == 'p') {
            hap = atoi(opt.arg);
        } else if (c == 'o') {
            out = opt.arg;
        } else if (c == 301) {
            LG_AGP_SEQ_COMPONENT_TYPE = agp_component_type_key(opt.arg);
            if (LG_AGP_SEQ_COMPONENT_TYPE == AGP_CT_N || 
                    LG_AGP_SEQ_COMPONENT_TYPE == AGP_CT_U)
                fprintf(stderr, "[W::%s] a GAP component identifier will be used for sequences: %s\n",
                        __func__, opt.arg);
        } else if (c == 302) {
            LG_AGP_GAP_COMPONENT_TYPE = agp_component_type_key(opt.arg);
            if (LG_AGP_GAP_COMPONENT_TYPE != AGP_CT_N && 
                    LG_AGP_GAP_COMPONENT_TYPE != AGP_CT_U)
                fprintf(stderr, "[W::%s] a SEQ component identifier will be used for gaps: %s\n",
                        __func__, opt.arg);
        } else if (c == 303) {
            LG_AGP_LINKAGE_EVIDENCE = agp_linkage_evidence_key(opt.arg);
        } else if (c == 304) {
            LG_AGP_GAP_SIZE = atoi(opt.arg);
        } else if (c == 'h') {
            fp_help = stdout;
        } else if (c == 'V') {
            puts(SEQTOOLS_VERSION);
            return 0;
        } else if (c == '?') {
            if (argv[opt.i - 1][1] == '?') {
                fp_help = stdout;
                LONG_HELP = 1;
            } else {
                fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
                return 1;
            }
        } else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (fp_help == stdout) {
        print_help_hap(stdout);
        return 0;
    }

    if (argc - opt.ind < 2) {
        fprintf(stderr, "[E::%s] missing input: two positional options required\n", __func__);
        print_help_hap(stderr);
        return 1;
    }

    agp = argv[opt.ind];
    ann = argv[opt.ind + 1];

    fo = out == 0? stdout : fopen(out, "w");
    if (fo == 0) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    
    annot_haps_from_agp(agp, ann, hap, fo);

    if (out != 0)
        fclose(fo);

    fprintf(stderr, "[M::%s] Version: %s\n", __func__, SEQTOOLS_VERSION);
    fprintf(stderr, "[M::%s] CMD:", __func__);
    int i;
    for (i = 0; i < argc; ++i)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - at_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
 
    return 0;
}

static void print_help_seq(FILE *fp_help)
{
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Usage: seqtools seq [options] <input.agp> <input.fa>\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help, "    -l INT            line width [60]\n");
    fprintf(fp_help, "    -u                allow sequence components with unknown orientations\n");
    fprintf(fp_help, "    -o STR            output to file [stdout]\n");
    fprintf(fp_help, "    -h, --help        print this help\n");
    fprintf(fp_help, "    -V, --version     show version number\n");
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Example: seqtools seq -o output.fa scaffolds.agp input.fa\n");
    fprintf(fp_help, "\n");
}

static ko_longopt_t seq_long_options[] = {
    { "help",       ko_no_argument,       'h' },
    { "version",    ko_no_argument,       'V' },
    { 0, 0, 0 }
};

int main_seq(int argc, char *argv[])
{
    if (argc < 2) {
        print_help_seq(stderr);
        return 1;
    }

    liftrlimit();
    at_realtime0 = realtime();

    FILE *fo;
    char *fa, *agp, *out;
    int line_wd, allow_unknown_oris;

    const char *opt_str = "o:ul:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c;
    FILE *fp_help = stderr;
    fa = agp = out = 0;
    line_wd = 60;
    allow_unknown_oris = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, seq_long_options)) >= 0) {
        if (c == 'l') {
            line_wd = atoi(opt.arg);
        } else if (c == 'u') {
            allow_unknown_oris = 1;
        } else if (c == 'o') {
            out = opt.arg;
        } else if (c == 'h') {
            fp_help = stdout;
        } else if (c == 'V') {
            puts(SEQTOOLS_VERSION);
            return 0;
        } else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        } else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (fp_help == stdout) {
        print_help_seq(stdout);
        return 0;
    }

    if (argc - opt.ind < 2) {
        fprintf(stderr, "[E::%s] missing input: two positional options required\n", __func__);
        print_help_seq(stderr);
        return 1;
    }

    agp = argv[opt.ind];
    fa = argv[opt.ind + 1];

    fo = out == 0? stdout : fopen(out, "w");
    if (fo == 0) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    
    write_fasta_file_from_agp(fa, agp, fo, line_wd, allow_unknown_oris);

    if (out != 0)
        fclose(fo);

    fprintf(stderr, "[M::%s] Version: %s\n", __func__, SEQTOOLS_VERSION);
    fprintf(stderr, "[M::%s] CMD:", __func__);
    int i;
    for (i = 0; i < argc; ++i)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - at_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
 
    return 0;
}

#define IDX_BUF_SIZE (1u << 24) /* 16MB read buffer */

static void write_seq_index(const char *fa, FILE *fo)
{
    enum { ST_NEW_LINE, ST_HEADER_NAME, ST_HEADER_REST, ST_SEQ };
    gzFile fp;
    char ch, *buf;
    kstring_t _name = {0, 0, NULL}, *name = &_name;
    int64 seq_len;
    int i, n, state, have_record;

    fp = gzopen(fa, "rb");
    if (fp == 0) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, fa);
        exit(EXIT_FAILURE);
    }

    MYMALLOC(buf, IDX_BUF_SIZE);
    if (!buf)
        mem_alloc_error("sequence read buffer");

    state = ST_NEW_LINE;
    have_record = 0;
    seq_len = 0;

    while ((n = gzread(fp, buf, IDX_BUF_SIZE)) > 0) {
        for (i = 0; i < n; ++i) {
            ch = buf[i];
            if (ch == '\r')
                continue;
            switch (state) {
                case ST_NEW_LINE:
                    if (ch == '>') {
                        if (have_record)
                            fprintf(fo, "%s\t%lld\n", name->s, (long long) seq_len);
                        name->l = 0;
                        seq_len = 0;
                        have_record = 1;
                        state = ST_HEADER_NAME;
                    } else if (ch != '\n') {
                        seq_len++;
                        state = ST_SEQ;
                    }
                    break;
                case ST_HEADER_NAME:
                    if (ch == '\n')
                        state = ST_NEW_LINE;
                    else if (isspace(ch))
                        state = ST_HEADER_REST;
                    else
                        kputc(ch, name);
                    break;
                case ST_HEADER_REST:
                    if (ch == '\n')
                        state = ST_NEW_LINE;
                    break;
                case ST_SEQ:
                    if (ch == '\n')
                        state = ST_NEW_LINE;
                    else
                        seq_len++;
                    break;
            }
        }
    }
    if (n < 0) {
        int errnum;
        const char *emsg = gzerror(fp, &errnum);
        fprintf(stderr, "[E::%s] error reading file %s: %s\n", __func__, fa, emsg);
        exit(EXIT_FAILURE);
    }

    if (have_record)
        fprintf(fo, "%s\t%lld\n", name->s, (long long) seq_len);

    free(buf);
    free(name->s);
    gzclose(fp);
}

static void print_help_idx(FILE *fp_help)
{
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Usage: seqtools idx [options] <input.fa>\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help, "    -o STR            output to file [stdout]\n");
    fprintf(fp_help, "    -h, --help        print this help\n");
    fprintf(fp_help, "    -V, --version     show version number\n");
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Example: seqtools idx -o input.fa.idx input.fa.gz\n");
    fprintf(fp_help, "\n");
}

static ko_longopt_t idx_long_options[] = {
    { "help",       ko_no_argument,       'h' },
    { "version",    ko_no_argument,       'V' },
    { 0, 0, 0 }
};

int main_idx(int argc, char *argv[])
{
    if (argc < 2) {
        print_help_idx(stderr);
        return 1;
    }

    liftrlimit();
    at_realtime0 = realtime();

    FILE *fo;
    char *fa, *out;

    const char *opt_str = "o:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c;
    FILE *fp_help = stderr;
    fa = out = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, idx_long_options)) >= 0) {
        if (c == 'o') {
            out = opt.arg;
        } else if (c == 'h') {
            fp_help = stdout;
        } else if (c == 'V') {
            puts(SEQTOOLS_VERSION);
            return 0;
        } else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        } else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (fp_help == stdout) {
        print_help_idx(stdout);
        return 0;
    }

    if (argc - opt.ind < 1) {
        fprintf(stderr, "[E::%s] missing input: one positional option required\n", __func__);
        print_help_idx(stderr);
        return 1;
    }

    fa = argv[opt.ind];

    fo = out == 0? stdout : fopen(out, "w");
    if (fo == 0) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }

    write_seq_index(fa, fo);

    if (out != 0)
        fclose(fo);

    fprintf(stderr, "[M::%s] Version: %s\n", __func__, SEQTOOLS_VERSION);
    fprintf(stderr, "[M::%s] CMD:", __func__);
    int i;
    for (i = 0; i < argc; ++i)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - at_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
 
    return 0;
}

/***********************************************************************************
 * top-level dispatcher                                                             *
 ***********************************************************************************/

static int usage(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "Program: seqtools\n");
    fprintf(fp, "Version: %s\n", SEQTOOLS_VERSION);
    fprintf(fp, "Usage:   seqtools <command> [options]\n");
    fprintf(fp, "Commands:\n");
    fprintf(fp, "    seq   generate fasta sequence from AGP file\n");
    fprintf(fp, "    hap   generate AGP file for haplotypes\n");
    fprintf(fp, "    idx   generate sequence index from FASTA file\n");
    fprintf(fp, "\n");
    return fp == stdout ? 0 : 1;
}

int main(int argc, char *argv[])
{
    sys_init();

    if (argc == 1)
        return usage(stderr);
    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
        return usage(stdout);
    if (strcmp(argv[1], "-V") == 0 || strcmp(argv[1], "--version") == 0) {
        puts(SEQTOOLS_VERSION);
        return 0;
    }
    if (strcmp(argv[1], "seq") == 0)
        return main_seq(argc - 1, argv + 1);
    if (strcmp(argv[1], "hap") == 0)
        return main_hap(argc - 1, argv + 1);
    if (strcmp(argv[1], "idx") == 0)
        return main_idx(argc - 1, argv + 1);
    fprintf(stderr, "[E::%s] unrecognised command '%s'. Abort!\n", __func__, argv[1]);
    return 1;
}
