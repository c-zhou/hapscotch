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

// BUSCO gene table loading and copy-number summary reporting

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>

#include "kseq.h"
#include "khash.h"
#include "kvec.h"

#include "busco.h"

KSTREAM_INIT(gzFile, gzread, 0x10000)

static int BUSCO_SUMMARY_COUNT[4];

/*********************** BUSCO Gene ***********************/
char BUSCO_SUMMARY_STRING[64];

static int busco_coords_cmpfunc(const void *a, const void *b)
{ 
    int64  xm, ym;
    busco_t *x = (busco_t *) a;
    busco_t *y = (busco_t *) b;
    
    xm = x->seq;
    ym = y->seq;

    if (xm == ym) {
        xm = x->end;
        ym = y->end;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->beg;
        ym = y->beg;
    }
    return (xm > ym) - (xm < ym);
}

static int seqrange_coords_acmpfunc(const void *a, const void *b)
{ 
    int64  xm, ym;
    seq_range_t *x = (seq_range_t *) a;
    seq_range_t *y = (seq_range_t *) b;
    
    xm = x->seq;
    ym = y->seq;

    if (xm == ym) {
        xm = x->beg;
        ym = y->beg;
    } else return (xm > ym) - (xm < ym);

    if (xm == ym) {
        xm = x->end;
        ym = y->end;
    }
    return (xm > ym) - (xm < ym);
}

void busco_destroy(busco_table_t *buscos)
{
    if (!buscos) return;
    int i, n;
    for (i = 0, n = buscos->n_busco; i < n; i++)
        free(buscos->busco_names[i]);
    free(buscos->busco_names);
    kh_destroy_busco(buscos->busco_table);
    free(buscos->gene_table);
    free(buscos);
}

void busco_summary_report(busco_table_t *buscos, seq_range_t *ranges, int nrange, int print_summary, int max_hist)
{
    if (!buscos || !ranges || nrange <= 0)
        return;

    int i, j, k, f, c, n_busco, n_gene, *gene_counts, *hist, *summary;
    busco_t *busco;

    summary = BUSCO_SUMMARY_COUNT;
    n_busco = buscos->n_busco;
    n_gene = buscos->n_gene;

    qsort(ranges, nrange, sizeof(seq_range_t), seqrange_coords_acmpfunc);

    // busco genes are sorted by end positions
    // sequence ranges are sorted by beg positions
    MYCALLOC(gene_counts, n_busco);
    k = 0;
    for (i = 0; i < n_gene; i++) {
        busco = &buscos->gene_table[i];
        if (busco->status == MISSING)
            continue;
        while (k < nrange && (busco->seq > ranges[k].seq || 
            (busco->seq == ranges[k].seq && 
            busco->end > ranges[k].beg)))
            k++;
        for (j = k - 1; j >= 0; j--) {
            if (busco->seq == ranges[j].seq && 
                busco->end > ranges[j].beg &&
                busco->beg < ranges[j].end) {
                if (busco->status == FRAGMENTED)
                    gene_counts[busco->gene] += 1;
                else
                    gene_counts[busco->gene] += (1<<16);
            }
            if (busco->seq > ranges[j].seq) break;
        }
    }

    if (max_hist > 0)
        MYCALLOC(hist, max_hist + 2);

    MYBZERO(summary, 4);
    for (i = 0; i < n_busco; i++) {
        c = gene_counts[i] >> 16;
        f = gene_counts[i] & 0xFFFF;
        if (c > 0) {
            summary[0]++;
            if (c + f > 1)
                summary[1]++;
        } else if (f > 0)
            summary[2]++;
        else
            summary[3]++;

        if (max_hist > 0) {
            c = c + f;
            if (c > max_hist)
                c = max_hist + 1;
            hist[c] += 1;
        }
    }
    free(gene_counts);

    // update summary
    sprintf(BUSCO_SUMMARY_STRING, "C:%.1f%%[S:%.1f%%,D:%.1f%%],F:%.1f%%,M:%.1f%%,n:%d",
        (double) summary[0]/n_busco*100., (double)(summary[0]-summary[1])/n_busco*100., 
        (double) (summary[1])/n_busco*100., (double) (summary[2])/n_busco*100., 
        (double) (summary[3])/n_busco*100., n_busco);

    if (print_summary) {
        fprintf(stderr, "[M::%s]    -----------------------------------------------\n", __func__);
        fprintf(stderr, "[M::%s]    ************ BUSCO Summary Results ************\n", __func__);
        fprintf(stderr, "[M::%s]    -----------------------------------------------\n", __func__);
        fprintf(stderr, "[M::%s]    %s\n", __func__, BUSCO_SUMMARY_STRING);
        fprintf(stderr, "[M::%s]    %-6d  Complete BUSCOs (C)\n", __func__, summary[0]);
        fprintf(stderr, "[M::%s]    %-6d  Complete and single-copy BUSCOs (S)\n", __func__, summary[0]-summary[1]);
        fprintf(stderr, "[M::%s]    %-6d  Complete and duplicated BUSCOs (D)\n", __func__, summary[1]);
        fprintf(stderr, "[M::%s]    %-6d  Fragmented BUSCOs (F)\n", __func__, summary[2]);
        fprintf(stderr, "[M::%s]    %-6d  Missing BUSCOs (M)\n", __func__, summary[3]);
        fprintf(stderr, "[M::%s]    %-6d  Total BUSCO groups searched\n", __func__, n_busco);
        fprintf(stderr, "[M::%s]    -----------------------------------------------\n", __func__);
    }
    
    if (max_hist > 0) {
        fprintf(stderr, "[M::%s] BUSCO copy number summary:\n", __func__);
        for (i = 0; i <= max_hist; i++)
            if (hist[i] > 0)
                fprintf(stderr, "[M::%s] %4d: %6d  %5.1f%%\n", __func__, i, hist[i], (100.*hist[i])/n_busco);
        if (hist[max_hist+1] > 0)
            fprintf(stderr, "[M::%s]  >%2d: %6d  %5.1f%%\n", __func__, max_hist, hist[max_hist+1], (100.*hist[max_hist+1])/n_busco);
        free(hist);
    }
}

void busco_summary_report_all_seqs(busco_table_t *buscos, sdict_t *dicts)
{
    if (!buscos) return;

    seq_range_t *ranges;
    int i, nseq;
    
    nseq = dicts->n;
    MYMALLOC(ranges, nseq);
    for (i = 0; i < nseq; i++)
        ranges[i] = (seq_range_t) {i, 0, dicts->s[i].len};
    busco_summary_report(buscos, ranges, nseq, 1, 12);
    free(ranges);
}

busco_table_t *build_busco_gene_table(const char *fn, sdict_t *dicts, asm_dict_t *break_dict)
{
    if (!fn) return NULL;

    busco_table_t *buscos;
    kstream_t *ks;
	gzFile fp;
    khash_t(busco) *busco_table;
    kstring_t buf = {0, 0, 0};
    kvec_t(busco_t) gene_table = {0, 0, 0};
    busco_t *busco;
    sdict_t *vd;
    char *fptr, *eptr, *fptrs[6];
    int i, n, n_busco, n_clip, absent, ret;
    khint32_t k;

    // input BUSCO coordinates are always on the raw (uncorrected) sequences
    vd = break_dict? break_dict->sdict : dicts;

    busco_table = kh_init(busco);
    n_busco = 0;
    n_clip = 0;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
        fprintf(stderr, "[E::%s] cannot open BUSCO gene table file %s to read\n", __func__, fn);
        exit (1);
    }
	ks = ks_init(fp);
    while ((ret = ks_getuntil(ks, KS_SEP_LINE, &buf, 0)) >= 0) {
        if (!buf.l || buf.s[0] == '#') continue;
        eptr = buf.s;
        n = 0;
        for (i = 0; i < 6; i++) {
            while (isspace(*eptr)) eptr++; // trim leading spaces
            fptr = eptr;
            while (!isspace(*eptr) && *eptr != '\0')
                eptr++;
            if (eptr > fptr) fptrs[n++] = fptr;
            if (*eptr == '\0')
                break;
            else
                *eptr++ = '\0';
            fptr = eptr;
        }
        if (n < 2 || (strcmp_case_insensitive(fptrs[1], "MISSING") && n < 6)) {
            fprintf(stderr, "[E::%s] error pasring BUSCO gene table file at line %s\n", __func__, buf.s);
            exit (1);
        }
        k = kh_put(busco, busco_table, fptrs[0], &absent);
        if (absent) {
            kh_key(busco_table, k) = strdup(fptrs[0]);
            kh_val(busco_table, k) = n_busco++;
        }
        kv_pushp(busco_t, gene_table, &busco);
        busco->gene = kh_val(busco_table, k);
        busco->status = !strcmp_case_insensitive(fptrs[1], "MISSING")? MISSING :
            (!strcmp_case_insensitive(fptrs[1], "FRAGMENTED")? FRAGMENTED : COMPLETE);
        if (busco->status != MISSING) {
            busco->seq = sd_get(vd, fptrs[2]); // may not exist
            busco->beg = strtol(fptrs[3], 0, 10);
            busco->end = strtol(fptrs[4], 0, 10);
            busco->strand = fptrs[5][0] == '+'? 0 : 1;
            // minus-strand genes may come with swapped coordinates
            if (busco->beg > busco->end)
                SWAP(int, busco->beg, busco->end);
            if (busco->seq < 0)
                fprintf(stderr, "[W::%s] sequence %s in the BUSCO gene table does not exist\n", __func__, fptrs[2]);
            else if (break_dict) {
                // remap gene coordinates (1-based) to piece coordinates
                // a gene spanning a break point goes to the piece with the larger overlap
                sd_seg_t *seg, *best;
                int64 b0, e0, ovlp, best_ovlp;
                uint32 si, sn;
                b0 = busco->beg - 1;
                e0 = busco->end - 1;
                best = NULL;
                best_ovlp = 0;
                si = (uint32) (break_dict->a[busco->seq] >> 32);
                sn = (uint32) break_dict->a[busco->seq] + si;
                for (; si < sn; si++) {
                    seg = &break_dict->seg[(uint32) break_dict->index[si]];
                    ovlp = MIN(e0, (int64) seg->x + seg->y - 1) - MAX(b0, (int64) seg->x) + 1;
                    if (ovlp > best_ovlp) { best_ovlp = ovlp; best = seg; }
                }
                if (best == NULL) {
                    fprintf(stderr, "[W::%s] gene position %s:%d-%d out of sequence range\n", __func__, fptrs[2], busco->beg, busco->end);
                    busco->seq = -1;
                } else {
                    if (best_ovlp < e0 - b0 + 1) n_clip++;
                    busco->seq = best->s;
                    busco->beg = (int) (MAX(b0, (int64) best->x) - best->x + 1);
                    busco->end = (int) (MIN(e0, (int64) best->x + best->y - 1) - best->x + 1);
                }
            }
        } else {
            busco->seq = INT32_MAX;
            busco->beg = busco->end = busco->strand = 0;
        }
    }

    free(buf.s);
    gzclose(ks->f);
    ks_destroy(ks);
    
    MYCALLOC(buscos, 1);
    MYMALLOC(buscos->busco_names, n_busco);
    for (k = (khint32_t) 0; k < kh_end(busco_table); k++)
        if (kh_exist(busco_table, k))
            buscos->busco_names[kh_val(busco_table, k)] = (char *) kh_key(busco_table, k);
    buscos->n_busco = n_busco;
    buscos->n_gene = gene_table.n;
    buscos->gene_table = MYREALLOC(gene_table.a, gene_table.n);
    buscos->busco_table = busco_table;

    qsort(buscos->gene_table, buscos->n_gene, sizeof(busco_t), busco_coords_cmpfunc);

    if (n_clip > 0)
        fprintf(stderr, "[W::%s] BUSCO genes crossing AGP break points: %d\n", __func__, n_clip);
    fprintf(stderr, "[M::%s] BUSCO gene table loaded\n", __func__);
    fprintf(stderr, "[M::%s] BUSCO number: %4d\n", __func__, buscos->n_busco);
    fprintf(stderr, "[M::%s] GENE  number: %4d\n", __func__, buscos->n_gene);
    
    return buscos;
}

/********************* END BUSCO Gene *********************/
