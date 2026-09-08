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

// Pseudo haplotype / scaffold construction: HiC linkage, haplotype
// partitioning (greedy / ILP / heuristic phasing) and AGP scaffold output.

#ifndef HAP_H
#define HAP_H

#include "misc.h"
#include "sdict.h"
#include "agp-spec.h"
#include "busco.h"
#include "overlap.h"

// verbose level, set by -v
extern int VERBOSE;

// AGP output bit flags for write_scf_outputs()
#define AGP_OUT 0x1
#define GRP_OUT 0x2
#define PLT_OUT 0x4

// HiGHS ILP solver options for HiC phasing, settable via --time-limit/--mip-rel-gap
extern long   HIGHS_TIME_LIMIT;
extern double HIGHS_MIP_REL_GAP;

typedef struct {
    int grp;
    int hap;
    uint64 len;
    int nctg;
    uint32 *ctgs; // id << 1 | rev
    uint64 *segs; // beg << 32 | len
} scf_t;

void scf_free(scf_t *scf);
int scaff_natural_cmpfunc(const void *a, const void *b);

scf_t *build_pseudo_scaffolds(ovl_t *ovls, int64 novl, sdict_t *dicts, asm_dict_t *break_dict, busco_table_t *buscos,
    int ploidy, int min_ext, int min_qual, int n_threads, int conf_yahs, char *hic_bfile,
    char *out_pref, int64 *_nscf);

void write_scf_outputs(scf_t *scfs, int nscf, sdict_t *dicts, asm_dict_t *break_dict, const uint8 opts_out, const char *pref_out);

#endif // HAP_H
