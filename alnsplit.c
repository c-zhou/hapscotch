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

// alnsplit: remap a .1aln or PAF alignment file onto the sequence pieces of a break AGP
// alignments are split at piece boundaries; the output .1aln file carries
// a GDB skeleton of the pieces and empty trace points

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "ketopt.h"

#include "sdict.h"
#include "alnio.h"
#include "misc.h"
#include "version.h"

int VERBOSE = 0;

static ko_longopt_t long_options[] = {
    { "output",         ko_required_argument, 'o' },
    { "raw-name",       ko_no_argument,       'r' },
    { "verbose",        ko_required_argument, 'v' },
    { "version",        ko_no_argument,       'V' },
    { "help",           ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{ 
    const char *opt_str = "o:rv:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, raw_name, out_alloc;
    FILE *fp_help, *fo;
    sdict_t *dicts;
    asm_dict_t *break_dict;
    char *agp_file, *aln_file, *out_file;
    int64 n_out;
    size_t n;
    
    sys_init();

    fp_help = stderr;
    out_file = 0;
    out_alloc = 0;
    raw_name = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0 ) {
        if (c == 'o') out_file = opt.arg;
        else if (c == 'r') raw_name = 1;
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'h') fp_help = stdout;
        else if (c == 'V') {
            puts(ALNSPLIT_VERSION);
            return 0;
        }
        else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
        else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (argc == opt.ind || fp_help == stdout) {
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Usage: alnsplit [options] <genome.fa.[fai|idx]> <break.agp> <aln.[1aln|paf[.gz]]>\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "    -o FILE                output 1aln file [<aln>_split.1aln]\n");
        fprintf(fp_help, "    -r                     name pieces as 'name_beg_end' on the original (1-based)\n");
        fprintf(fp_help, "    -v INT                 verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    -h, --help             print this help\n");
        fprintf(fp_help, "    -V, --version          show version number\n");
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Note:    the AGP file must break every input sequence into '+' oriented single-component\n");
        fprintf(fp_help, "         pieces covering every base exactly once, e.g., hapscotch '.ec.agp' output\n");
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Example: alnsplit genome.fa.fai ec.agp aln.paf.gz\n");
        fprintf(fp_help, "         alnsplit -r -o aln.ec.1aln genome.fa.fai ec.agp aln.1aln\n");
        fprintf(fp_help, "\n");
        return fp_help == stdout? 0 : 1;
    }

    // sanity checks
    if (argc - opt.ind != 3) {
        fprintf(stderr, "[E::%s] please specify three positional parameters\n", __func__);
        return 1;
    }

    agp_file = argv[opt.ind + 1];
    aln_file = argv[opt.ind + 2];

    // read sequence dictionary
    dicts = make_sdict_from_index(argv[opt.ind], 0);

    // read AGP file of the sequence pieces
    break_dict = make_asm_dict_from_agp(dicts, agp_file, 0);
    if (break_dict == NULL) {
        fprintf(stderr, "[E::%s] failed to parse AGP file %s\n", __func__, agp_file);
        sd_destroy(dicts);
        return 1;
    }
    validate_break_agp(break_dict);

    // output file: <dir>/<name>.[1aln|paf[.gz]] -> <dir>/<name>_split.1aln if not specified
    if (!out_file) {
        n = strlen(aln_file);
        if (n > 3 && !strcmp(aln_file + n - 3, ".gz")) n -= 3;
        if (n > 5 && !strncmp(aln_file + n - 5, ".1aln", 5)) n -= 5;
        else if (n > 4 && !strncmp(aln_file + n - 4, ".paf", 4)) n -= 4;
        MYMALLOC(out_file, n + 12);
        memcpy(out_file, aln_file, n);
        strcpy(out_file + n, "_split.1aln");
        out_alloc = 1;
        fprintf(stderr, "[M::%s] no output file specified (-o), writing to %s\n", __func__, out_file);
    }
    if (!strcmp(out_file, aln_file)) {
        fprintf(stderr, "[E::%s] output file is the same as the input file: %s\n", __func__, out_file);
        asm_destroy(break_dict);
        sd_destroy(dicts);
        return 1;
    }
    fo = fopen(out_file, "wb");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out_file);
        asm_destroy(break_dict);
        sd_destroy(dicts);
        return 1;
    }

    // remap the alignments
    n_out = split_aln(aln_file, break_dict, raw_name, fo);

    fclose(fo);
    asm_destroy(break_dict);
    sd_destroy(dicts);

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(1);
    }

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] %lld alignments written to %s\n", __func__, n_out, out_file);
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, ALNSPLIT_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    if (out_alloc) free(out_file);

    return 0;
}
