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
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <errno.h>
#include <math.h>
#include <zlib.h>
#include <pthread.h>

#include "ketopt.h"
#include "kvec.h"
#include "kseq.h"
#include "khash.h"
#include "kthread.h"

#include "paf.h"
#include "hic.h"
#include "sdict.h"
#include "range.h"
#include "agp-spec.h"
#include "misc.h"

#include "busco.h"
#include "overlap.h"
#include "ploidy.h"
#include "hap.h"
#include "alnio.h"
#include "version.h"

int VERBOSE = 0;

static inline long parse_time(const char *str, long default_val)
{
    long val;
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'M' || *p == 'm') x *= 60, ++p;
	else if (*p == 'H' || *p == 'h') x *= 60 * 60, ++p;
	else if (*p == 'D' || *p == 'd') x *= 24 * 60 * 60, ++p;
	val = (long) (x + .499);
    return val > 0 ? val : default_val;
}

static inline int64 parse_num2(const char *str, char **q)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	if (q) *q = p;
	return (int64) (x + .499);
}

static inline int64 parse_num(const char *str)
{
	return parse_num2(str, 0);
}

static ko_longopt_t long_options[] = {
    { "read-length",    ko_required_argument, 301 },
    { "file-type",      ko_required_argument, 302 },
    { "time-limit",     ko_required_argument, 303 },
    { "mip-rel-gap",    ko_required_argument, 304 },
    { "hic-file",       ko_required_argument, 'c' },
    { "agp-file",       ko_required_argument, 'a' },
    { "min-extension",  ko_required_argument, 'B' },
    { "ploidy",         ko_required_argument, 'p' },
    { "max-ploidy",     ko_required_argument, 'P' },
    { "min-quality",    ko_required_argument, 'q' },
    { "busco-table",    ko_required_argument, 'g' },
    { "config-yahs",    ko_no_argument,       'Y' },
    { "verbose",        ko_required_argument, 'v' },
    { "version",        ko_no_argument,       'V' },
    { "help",           ko_no_argument,       'h' },
    { 0, 0, 0 }
};

static inline int yes_or_no(int long_idx, const char *arg)
{
	if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0)
        return 1;
	else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0)
        return 0;
	else {
        fprintf(stderr, "[E::%s]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", 
            __func__, long_options[long_idx].name);
        exit(1);
    }
}

int main(int argc, char *argv[])
{ 
    const char *opt_str = "a:c:DYB:p:P:t:l:g:o:v:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, long_help = 0, ret = 0;
    int n_threads;
    FILE *fp_help;
    sdict_t *dicts, *dicts_raw;
    asm_dict_t *break_dict;
    aln_t *alns;
    ovl_t *ovls;
    scf_t *scfs;
    busco_table_t *buscos;
    fileType_t f_type;
    int64 i, naln, novl, nscf;
    uint8 opts_out, conf_yahs;
    int ploidy_num, dual_aln, min_ext, read_len, min_qual;
    char *busco_file, *hic_file, *hic_bfile, *agp_file, *pref_out;
    
    sys_init();
    srand48(42);

    fp_help = stderr;
    ploidy_num = 0;
    min_ext = 50000;
    read_len = 150;
    min_qual = 1;
    dual_aln = 1;
    n_threads = 1;
    busco_file = 0;
    hic_file = 0;
    hic_bfile = 0;
    agp_file = 0;
    conf_yahs = 0;
    f_type = NOSET;
    pref_out = "hapscotch.out";

    opts_out = 0;
    opts_out |= AGP_OUT;
    opts_out |= GRP_OUT;
    opts_out |= PLT_OUT;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0 ) {
        if (c == 'p') ploidy_num = atoi(opt.arg);
        else if (c == 'D') dual_aln = 0;
        else if (c == 'a') agp_file = opt.arg;
        else if (c == 'c') hic_file = opt.arg;
        else if (c == 'q') min_qual = atoi(opt.arg);
        else if (c == 'g') busco_file = opt.arg;
        else if (c == 'B') min_ext = parse_num(opt.arg);
        else if (c == 'P') MAX_PLOIDY_NUMBER = atoi(opt.arg);
        else if (c == 'Y') conf_yahs = 1;
        else if (c == 't') n_threads = atoi(opt.arg);
        else if (c == 'o') pref_out = opt.arg;
        else if (c == 301) read_len = atoi(opt.arg);
        else if (c == 302) {
            if (strcasecmp(opt.arg, "BED") == 0)
                f_type = BED;
            else if (strcasecmp(opt.arg, "BAM") == 0)
                f_type = BAM;
            else if (strcasecmp(opt.arg, "BIN") == 0)
                f_type = BIN;
            else if (strcasecmp(opt.arg, "PA5") == 0)
                f_type = PA5;
            else if (strcasecmp(opt.arg, "ONE") == 0)
                f_type = ONE;
            else {
                fprintf(stderr, "[E::%s] unknown file type: \"%s\"\n", __func__, opt.arg);
                return 1;
            }
        }
        else if (c == 303) HIGHS_TIME_LIMIT = parse_time(opt.arg, HIGHS_TIME_LIMIT);
        else if (c == 304) HIGHS_MIP_REL_GAP = atof(opt.arg);
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'h') fp_help = stdout;
        else if (c == 'V') {
            puts(HAPSCOTCH_VERSION);
            return 0;
        }
        else if (c == '?') {
            if (argv[opt.i - 1][1] == '?') {
                fp_help = stdout;
                long_help = 1;
            } else {
                fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
                return 1;
            }
        }
        else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (argc == opt.ind || fp_help == stdout) {
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Usage: hapscotch [options] genome.fa[.gz] aln.paf[.gz]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "    -g FILE                BUSCO gene full table (for statistics only)\n");
        fprintf(fp_help, "    -a FILE                AGP file of assembly error corrected sequences\n");
        fprintf(fp_help, "    -D                     input alignments are not dual mappings\n");
        fprintf(fp_help, "    -Y                     write files needed for YaHS scaffolding\n");
        fprintf(fp_help, "    -B NUM                 minimum sequence size for scaffold extension bridging [50k]\n");
        fprintf(fp_help, "    -p INT                 genome ploidy number, set to 0 for auto estimation [%d]\n", ploidy_num);
        fprintf(fp_help, "      --max-ploidy INT     upper limit of genome ploidy number to consider [%d]\n", MAX_PLOIDY_NUMBER);
        fprintf(fp_help, "    -c FILE                alignment file of HiC data to the genome sequences\n");
        fprintf(fp_help, "      --file-type   STR    input file type BED|BAM|PA5|ONE|BIN, file name extension is ignored\n");
        fprintf(fp_help, "      --read-length INT    read length (required for PA5 format input) [%d]\n", read_len);
        fprintf(fp_help, "      --min-quality INT    minimum mapping quality of reads used for phasing [%d]\n", min_qual);
        fprintf(fp_help, "    -t INT                 maximum number of threads to use [%d]\n", n_threads);
        fprintf(fp_help, "    -o STR                 string prefix of output files [%s]\n", pref_out);
        if (long_help) {
            fprintf(fp_help, "\n");
            fprintf(fp_help, "    HiGHS options for HiC phasing:\n");
        fprintf(fp_help, "      --time-limit STR     time limit for HiC phasing optimisation [%lds]\n", HIGHS_TIME_LIMIT);
        fprintf(fp_help, "      --mip-rel-gap FLOAT  relative gap tolerance for HiC phasing optimisation [%.1e]\n", HIGHS_MIP_REL_GAP);
            fprintf(fp_help, "\n");
        }
        fprintf(fp_help, "    -v INT                 verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    -?                     print long help with extra option list\n");
        fprintf(fp_help, "    -h, --help             print this help\n");
        fprintf(fp_help, "    -V, --version          show version number\n");
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Example: ./hapscotch -o hapscotch.out genome.fa.gz aln.paf\n");
        fprintf(fp_help, "\n");
        return fp_help == stdout? 0 : 1;
    }

    // sanity checks
    if (argc - opt.ind != 2) {
        fprintf(stderr, "[E::%s] please specify two positional parameters\n", __func__);
        return 1;
    }

    if (ploidy_num < 0) {
        fprintf(stderr, "[E::%s] ploidy number (-p) must be a non-negative integer\n", __func__);
        return 1;
    }

    if (MAX_PLOIDY_NUMBER < 1) {
        fprintf(stderr, "[E::%s] maximum ploidy number (--max-ploidy) must be a positive integer\n", __func__);
        return 1;
    }

    if (min_qual < 0 || min_qual > 255) {
        fprintf(stderr, "[E::%s] invalid mapping quality threshold: %d\n", __func__, min_qual);
        return 1;
    }

    if (read_len < 0) {
        fprintf(stderr, "[E::%s] invalid read length: %d\n", __func__, read_len);
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
    
    // add dual alignments and sort by aread, abpos, aepos
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

#ifdef DEBUG_PRINT_OVERLAP
    for (i = 0; i < dicts->n; i++)
        printf("S\t%s\t%d\n", dicts->s[i].name, dicts->s[i].len);
    for (i = 0; i < novl; i++)
        printf("O\t%s\t%d\t%d\t%s\t%d\t%d\t%.3f\t%.3f\t%.3f\n", 
            dicts->s[ovls[i].aread].name, ovls[i].abpos, ovls[i].aepos, 
            dicts->s[ovls[i].bread].name, ovls[i].bbpos, ovls[i].bepos,
            ovls[i].neff, ovls[i].score, ovls[i].qual);
#endif

    // estimate ploidy number if not specified
    if (!ploidy_num) {
        ploidy_num = estimate_ploidy_number(ovls, novl, dicts);
        if (ploidy_num < 2) {
            fprintf(stderr, "[W::%s] the genome will be treated as haploid\n", __func__);
            fprintf(stderr, "[W::%s] specify ploidy number with option '-p' if this is wrong\n", __func__);
        }
        if (ploidy_num > 8) {
            fprintf(stderr, "[W::%s] the estimated ploidy number (%d) is quite high\n", __func__, ploidy_num);
            fprintf(stderr, "[W::%s] please double check the input data or specify ploidy number with option '-p'\n", __func__);
        }
        if (ploidy_num & 1) {
            fprintf(stderr, "[W::%s] the estimated ploidy number (%d) is odd\n", __func__, ploidy_num);
            fprintf(stderr, "[W::%s] please double check the input data or specify ploidy number with option '-p'\n", __func__);
        }
    }

    // detect structure variants or misassemblies
    //detect_structural_variants(ovls, novl, dicts, ploidy_num);


    // convert hic data to binary format for fast access
    // HiC alignments are always on the raw (uncorrected) sequences
    hic_bfile = NULL;
    hic_bfile = write_binary_hic_data(hic_file, f_type, dicts_raw, read_len, pref_out);

    // build pseudo scaffolds
    nscf = 0;
    scfs = build_pseudo_scaffolds(ovls, novl, dicts, break_dict, buscos, ploidy_num, min_ext, min_qual, n_threads, conf_yahs, hic_bfile, pref_out, &nscf);

    // sort by natural order for outputs
    qsort(scfs, nscf, sizeof(scf_t), scaff_natural_cmpfunc);
    
    write_scf_outputs(scfs, nscf, dicts, break_dict, opts_out, pref_out);
    for (i = 0; i < nscf; i++) 
        scf_free(scfs+i);
    free(scfs);
    free(ovls);
    free(hic_bfile);
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
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, HAPSCOTCH_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}
