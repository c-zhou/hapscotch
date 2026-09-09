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

// hapcount: estimate genome ploidy number from self alignments (PAF)
// extracted from hapscotch.c, keeping only the pipeline steps needed to
// reach the ploidy_estimation step

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <errno.h>
#include <math.h>

#include "ketopt.h"
#include "kvec.h"

#include "paf.h"
#include "sdict.h"
#include "misc.h"

#include "busco.h"
#include "overlap.h"
#include "ploidy.h"
#include "alnio.h"
#include "version.h"

int VERBOSE = 0;

static ko_longopt_t long_options[] = {
    { "agp-file",       ko_required_argument, 'a' },
    { "busco-table",    ko_required_argument, 'g' },
    { "max-ploidy",     ko_required_argument, 'P' },
    { "verbose",        ko_required_argument, 'v' },
    { "version",        ko_no_argument,       'V' },
    { "help",           ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{ 
    const char *opt_str = "a:Dg:P:t:o:v:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, ret = 0;
    int n_threads;
    FILE *fp_help, *fo;
    sdict_t *dicts, *dicts_raw;
    asm_dict_t *break_dict;
    aln_t *alns;
    ovl_t *ovls;
    busco_table_t *buscos;
    int64 naln, novl;
    int ploidy_num, dual_aln;
    char *busco_file, *agp_file;
    
    sys_init();
    srand48(42);

    fp_help = stderr;
    dual_aln = 1;
    n_threads = 1;
    busco_file = 0;
    agp_file = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0 ) {
        if (c == 'D') dual_aln = 0;
        else if (c == 'a') agp_file = opt.arg;
        else if (c == 'g') busco_file = opt.arg;
        else if (c == 'P') MAX_PLOIDY_NUMBER = atoi(opt.arg);
        else if (c == 't') n_threads = atoi(opt.arg);
        else if (c == 'o') {
            if (strcmp(opt.arg, "-") != 0) {
                if (freopen(opt.arg, "wb", stdout) == NULL) {
                    fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", opt.arg, 
strerror(errno));
                    return 1;
                }
            }
        }
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'h') fp_help = stdout;
        else if (c == 'V') {
            puts(HAPCOUNT_VERSION);
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
        fprintf(fp_help, "Usage: hapcount [options] <genome.fa.[fai|idx]> <aln.paf[.gz]>\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "    -g FILE                BUSCO gene full table (for statistics only)\n");
        fprintf(fp_help, "    -a FILE                AGP file of assembly error corrected sequences\n");
        fprintf(fp_help, "    -D                     input alignments are not dual mappings\n");
        fprintf(fp_help, "      --max-ploidy INT     upper limit of genome ploidy number to consider [%d]\n", MAX_PLOIDY_NUMBER);
        fprintf(fp_help, "    -t INT                 maximum number of threads to use [%d]\n", n_threads);
        fprintf(fp_help, "    -o STR                 output file name\n");
        fprintf(fp_help, "    -v INT                 verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    -h, --help             print this help\n");
        fprintf(fp_help, "    -V, --version          show version number\n");
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Example: ./hapcount -o hapcount.out genome.fa.fai aln.paf\n");
        fprintf(fp_help, "\n");
        return fp_help == stdout? 0 : 1;
    }

    // sanity checks
    if (argc - opt.ind != 2) {
        fprintf(stderr, "[E::%s] please specify two positional parameters\n", __func__);
        return 1;
    }

    if (MAX_PLOIDY_NUMBER < 1) {
        fprintf(stderr, "[E::%s] maximum ploidy number (--max-ploidy) must be a positive integer\n", __func__);
        return 1;
    }

    // read sequence dictionary
    dicts_raw = make_sdict_from_index(argv[opt.ind], 0);

    // read AGP file of error corrected sequences
    // the working dictionary is made of the corrected sequence pieces
    break_dict = NULL;
    dicts = dicts_raw;
    if (agp_file) {
        break_dict = make_asm_dict_from_agp(dicts_raw, agp_file, 0);
        if (break_dict == NULL) {
            fprintf(stderr, "[E::%s] failed to parse AGP file %s\n", __func__, agp_file);
            return 1;
        }
        validate_break_agp(break_dict);
        dicts = make_piece_sdict(break_dict);
        fprintf(stderr, "[M::%s] %u sequences broken into %u pieces by AGP file %s\n", __func__, dicts_raw->n, dicts->n, agp_file);
    }

    // read busco gene table
    buscos = build_busco_gene_table(busco_file, dicts, break_dict);
    busco_summary_report_all_seqs(buscos, dicts);

    // read PAF files
    naln = 0;
    alns = read_pafs(argv + opt.ind + 1, argc - opt.ind - 1, dicts, break_dict, dual_aln, &naln);
    
    // sort by aread, abpos, aepos
    qsort(alns, naln, sizeof(aln_t), aln_coords_cmpfunc);
    report_genome_coverage_histogram(alns, naln, dicts, 0);
    
    // alignment chaining
    // ovls are non-symmetric
    novl = 0;
    ovls = build_adaptive_chains(alns, naln, dicts, n_threads, &naln, &novl);

    // alns are memory heavy
    free(alns);
    
    // sort sequence overlap information by coords
    ovls = add_dual_overlaps(ovls, novl, &novl);
    qsort(ovls, novl, sizeof(ovl_t), ovl_abseqs_cmpfunc);

    // estimate ploidy number
    ploidy_num = estimate_ploidy_number(ovls, novl, dicts);
    if (ploidy_num < 2)
        fprintf(stderr, "[W::%s] the genome will be treated as haploid\n", __func__);
    if (ploidy_num > 8)
        fprintf(stderr, "[W::%s] the estimated ploidy number (%d) is quite high\n", __func__, ploidy_num);
    if (ploidy_num & 1)
        fprintf(stderr, "[W::%s] the estimated ploidy number (%d) is odd\n", __func__, ploidy_num);

    // write output
    fprintf(stdout, "Estimated_ploidy\t%d\n", ploidy_num);
    
    free(ovls);
    sd_destroy(dicts);
    if (break_dict) {
        asm_destroy(break_dict);
        sd_destroy(dicts_raw);
    }
    busco_destroy(buscos);
    
    if (ret) {
        fprintf(stderr, "[E::%s] failed to analysis the PAF file\n", __func__);
        exit(1);
    }

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(1);
    }

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, HAPCOUNT_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}
