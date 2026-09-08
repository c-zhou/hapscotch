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

#include "ketopt.h"

#include "ec.h"
#include "hic.h"
#include "misc.h"
#include "sdict.h"
#include "version.h"

int VERBOSE = 0;

static int64 parse_num(const char *str)
{
    char *p;
    double x = strtod(str, &p);
    if (*p == 'G' || *p == 'g') x *= 1e9;
    else if (*p == 'M' || *p == 'm') x *= 1e6;
    else if (*p == 'K' || *p == 'k') x *= 1e3;
    return (int64) (x + .499);
}

static int parse_file_type(const char *str, fileType_t *f_type)
{
    if (strcasecmp(str, "BED") == 0) *f_type = BED;
    else if (strcasecmp(str, "BAM") == 0) *f_type = BAM;
    else if (strcasecmp(str, "BIN") == 0) *f_type = BIN;
    else if (strcasecmp(str, "PA5") == 0) *f_type = PA5;
    else if (strcasecmp(str, "ONE") == 0) *f_type = ONE;
    else return 1;
    return 0;
}

static void print_help(FILE *fo)
{
    fprintf(fo, "\n");
    fprintf(fo, "Usage: hapcure [options] <genome.fa.fai> <hic.bed|hic.bam|hic.pa5|hic.1map|hic.bin>\n");
    fprintf(fo, "Options:\n");
    fprintf(fo, "    -s NUM                 bin size [1k]\n");
    fprintf(fo, "    -d NUM                 maximum cis contact distance [1M]\n");
    fprintf(fo, "    -w NUM                 local control window [30k]\n");
    fprintf(fo, "    -m NUM                 minimum retained fragment size [30k]\n");
    fprintf(fo, "    -c INT                 minimum expected contacts [10]\n");
    fprintf(fo, "    -r FLOAT               maximum observed/expected ratio [.35]\n");
    fprintf(fo, "    -f FLOAT               Benjamini-Hochberg FDR [.01]\n");
    fprintf(fo, "      --file-type   STR    input file type BED|BAM|PA5|ONE|BIN, file name extension is ignored\n");
    fprintf(fo, "      --read-length INT    read length (required for PA5 format input) [150]\n");
    fprintf(fo, "      --report-only        do not write an AGP file\n");
    fprintf(fo, "    -o STR                 output prefix [hapcure.out]\n");
    fprintf(fo, "    -v INT                 verbose level [0]\n");
    fprintf(fo, "    -h, --help             print this help\n");
    fprintf(fo, "    -V, --version          show version number\n");
    fprintf(fo, "\n");
    fprintf(fo, "Example: ./hapcure -o hapcure.out genome.fa.gz hic.bam\n");
    fprintf(fo, "\n");
}

static ko_longopt_t long_options[] = {
    { "bin-size",     ko_required_argument, 's' },
    { "max-distance", ko_required_argument, 'd' },
    { "window-size",  ko_required_argument, 'w' },
    { "min-fragment", ko_required_argument, 'm' },
    { "min-contacts", ko_required_argument, 'c' },
    { "max-ratio",    ko_required_argument, 'r' },
    { "fdr",          ko_required_argument, 'f' },
    { "file-type",    ko_required_argument, 301 },
    { "read-length",  ko_required_argument, 302 },
    { "report-only",  ko_no_argument,       303 },
    { "version",      ko_no_argument,       'V' },
    { "help",         ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "o:t:l:s:d:w:m:v:Vh";
    ketopt_t opt = KETOPT_INIT;
    ec_conf_t conf;
    ec_candidates_t cands;
    asm_dict_t *adict;
    sdict_t *dicts;
    hic_t *hics;
    fileType_t f_type;
    char *hic_bfile, *report_file, *agp_file;
    FILE *fo;
    int c, read_len, report_only, ret;
    int64 nhic;
    char *pref_out;

    sys_init();
    realtime0 = realtime();
    ec_conf_init(&conf);
    f_type = NOSET;
    pref_out = "hapcure.out";
    read_len = 150;
    report_only = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0) {
        if (c == 'o') pref_out = opt.arg;
        else if (c == 's') conf.bin_size = parse_num(opt.arg);
        else if (c == 'd') conf.max_distance = parse_num(opt.arg);
        else if (c == 'w') conf.control_window = parse_num(opt.arg);
        else if (c == 'm') conf.min_fragment = parse_num(opt.arg);
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'c') conf.min_contacts = atoi(opt.arg);
        else if (c == 'r') conf.max_ratio = atof(opt.arg);
        else if (c == 'f') conf.fdr = atof(opt.arg);
        else if (c == 301) {
            if (parse_file_type(opt.arg, &f_type)) {
                fprintf(stderr, "[E::%s] unknown file type: %s\n", __func__, opt.arg);
                return 1;
            }
        }
        else if (c == 302) read_len = atoi(opt.arg);
        else if (c == 303) report_only = 1;
        else if (c == 'h') {
            print_help(stdout);
            return 0;
        } else if (c == 'V') {
            puts(HAPCURE_VERSION);
            return 0;
        } else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: %s\n", __func__, argv[opt.i - 1]);
            return 1;
        } else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option argument: %s\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (argc - opt.ind != 2) {
        print_help(stderr);
        return 1;
    }
    if (read_len < 0 || conf.bin_size <= 0 || conf.max_distance < conf.bin_size ||
        conf.control_window < conf.bin_size || conf.min_fragment < conf.bin_size ||
        conf.min_contacts < 1 || conf.max_ratio <= 0. || conf.max_ratio >= 1. ||
        conf.fdr <= 0. || conf.fdr > 1.) {
        fprintf(stderr, "[E::%s] invalid correction parameters\n", __func__);
        return 1;
    }

    dicts = make_sdict_from_index(argv[opt.ind], 0);
    hic_bfile = write_binary_hic_data(argv[opt.ind + 1], f_type, dicts, read_len, pref_out);
    ret = match_binary_file_sdict(hic_bfile, dicts);
    if (ret) {
        fprintf(stderr, "[E::%s] HiC BIN sequence dictionary does not match genome: %d\n", __func__, ret);
        free(hic_bfile);
        sd_destroy(dicts);
        return 1;
    }

    adict = make_asm_dict_from_sdict(dicts);
    nhic = 0;
    hics = read_hic_from_binary_sd_conversion(hic_bfile, adict, conf.bin_size, 0, &nhic);
    asm_destroy(adict);
    if (hics == NULL || nhic == 0) {
        fprintf(stderr, "[E::%s] no usable HiC contacts found\n", __func__);
        free(hic_bfile);
        sd_destroy(dicts);
        return 1;
    }

    if (ec_call_breaks(hics, nhic, dicts, &conf, &cands)) {
        fprintf(stderr, "[E::%s] failed to call correction breakpoints\n", __func__);
        free(hics);
        free(hic_bfile);
        sd_destroy(dicts);
        return 1;
    }
    free(hics);

    MYMALLOC(report_file, strlen(pref_out) + 12);
    sprintf(report_file, "%s.ec.tsv", pref_out);
    fo = fopen(report_file, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot write file %s\n", __func__, report_file);
        free(report_file);
        ec_candidates_destroy(&cands);
        free(hic_bfile);
        sd_destroy(dicts);
        return 1;
    }
    ec_write_report(&cands, dicts, fo);
    fclose(fo);
    fprintf(stderr, "[M::%s] wrote break report: %s\n", __func__, report_file);
    free(report_file);

    if (!report_only) {
        MYMALLOC(agp_file, strlen(pref_out) + 12);
        sprintf(agp_file, "%s.ec.agp", pref_out);
        fo = fopen(agp_file, "w");
        if (fo == NULL) {
            fprintf(stderr, "[E::%s] cannot write file %s\n", __func__, agp_file);
            free(agp_file);
            ec_candidates_destroy(&cands);
            free(hic_bfile);
            sd_destroy(dicts);
            return 1;
        }
        ec_write_agp(&cands, dicts, fo);
        fclose(fo);
        fprintf(stderr, "[M::%s] wrote corrected AGP: %s\n", __func__, agp_file);
        free(agp_file);
    }

    ec_candidates_destroy(&cands);
    free(hic_bfile);
    sd_destroy(dicts);

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, HAPCURE_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}