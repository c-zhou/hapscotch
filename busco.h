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

#ifndef BUSCO_H
#define BUSCO_H

#include "khash.h"
#include "sdict.h"
#include "agp-spec.h"
#include "misc.h"

#define MISSING    0
#define COMPLETE   1
#define FRAGMENTED 2

KHASH_MAP_INIT_STR(busco, uint32)

typedef struct {
    int gene;
    int seq, beg, end;
    int8 status:7, strand:1;
} busco_t;

typedef struct {
    int n_busco;
    kh_busco_t *busco_table;
    char **busco_names;
    int n_gene;
    busco_t *gene_table;
} busco_table_t;

typedef struct {
    int seq;
    int beg;
    int end;
} seq_range_t;

// summary string set by the most recent call to busco_summary_report()
// e.g., "C:98.5%[S:97.0%,D:1.5%],F:0.8%,M:0.7%,n:1234"
extern char BUSCO_SUMMARY_STRING[64];

busco_table_t *build_busco_gene_table(const char *fn, sdict_t *dicts, asm_dict_t *break_dict);
void busco_destroy(busco_table_t *buscos);
void busco_summary_report(busco_table_t *buscos, seq_range_t *ranges, int nrange, int print_summary, int max_hist);
void busco_summary_report_all_seqs(busco_table_t *buscos, sdict_t *dicts);

#endif // BUSCO_H
