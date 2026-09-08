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

// PAF alignment loading, with optional break-AGP coordinate remapping,
// shared by hapscotch.c and hapcount.c

#ifndef ALNIO_H
#define ALNIO_H

#include "misc.h"
#include "sdict.h"
#include "overlap.h"

// validate a break AGP: every object is a single '+' oriented W-line slice of
// an input sequence and every base of every input sequence is covered exactly once
void validate_break_agp(asm_dict_t *bd);

// build a sequence dictionary from the objects (sequence pieces) of a break AGP
sdict_t *make_piece_sdict(asm_dict_t *bd);

// read PAF file(s) into an aln_t array, remapping onto break-AGP pieces if break_dict is given
aln_t *read_pafs(char **fs, int fn, sdict_t *dicts, asm_dict_t *break_dict, int dual_aln, int64 *_naln);

#endif // ALNIO_H
