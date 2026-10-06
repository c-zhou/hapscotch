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
#include "paf.h"
#include "misc.h"
#include "alnio.h"
#include "sdict.h"
#include "overlap.h"
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
    fprintf(fo, "Usage: hapcure [options] <genome.fa.[fai|idx]> <hic.[bed|bam|pa5|1map|bin]>\n");
    fprintf(fo, "Options:\n");
    fprintf(fo, "    -s NUM                 bin size for counting hic links [1k]\n");
    fprintf(fo, "    -f NUM                 minimum retained fragment size [10k]\n");
    fprintf(fo, "    -d FLOAT               minimum median drop to call a break [.3]\n");
    fprintf(fo, "    -r FLOAT               minimum recovery rate for a breakpoint [0.8]\n");
    fprintf(fo, "    -p FLOAT               p-value threshold for calling errors [0.01]\n");
    fprintf(fo, "    -a FILE                self-alignment file for palindrome error detection\n");
    fprintf(fo, "      --pal-size NUM       minimum palindrome size to detect [100k]\n");
    fprintf(fo, "      --pal-dens FLOAT     minimum alignment density to call a palindrome [0.8]\n");
    fprintf(fo, "      --file-type   STR    input file type BED|BAM|PA5|ONE|BIN, file name extension is ignored\n");
    fprintf(fo, "      --read-length INT    read length (required for PA5 format input) [150]\n");
    fprintf(fo, "      --report-only        do not write an AGP file\n");
    fprintf(fo, "    -o STR                 output prefix [hapcure.out]\n");
    fprintf(fo, "    -v INT                 verbose level [0]\n");
    fprintf(fo, "    -h, --help             print this help\n");
    fprintf(fo, "    -V, --version          show version number\n");
    fprintf(fo, "\n");
    fprintf(fo, "Example: hapcure -o hapcure.out genome.fa.fai hic.bam\n");
    fprintf(fo, "         hapcure -a self_alignments.paf -o hapcure.out genome.fa.fai hic.bam\n");
    fprintf(fo, "\n");
}

static ko_longopt_t long_options[] = {
    { "bin-size",     ko_required_argument, 's' },
    { "min-frag",     ko_required_argument, 'f' },
    { "med-drop",     ko_required_argument, 'd' },
    { "rec-rate",     ko_required_argument, 'r' },
    { "p-thresh",     ko_required_argument, 'p' },
    { "file-type",    ko_required_argument, 301 },
    { "read-length",  ko_required_argument, 302 },
    { "report-only",  ko_no_argument,       303 },
    { "pal-size",     ko_required_argument, 304 },
    { "pal-dens",     ko_required_argument, 305 },
    { "version",      ko_no_argument,       'V' },
    { "help",         ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "o:s:f:d:r:p:t:a:v:Vh";
    ketopt_t opt = KETOPT_INIT;
    sdict_t *dicts;
    hic_t *hics;
    aln_t *alns;
    fileType_t f_type;
    FILE *fo;
    ec_pos_t *calls;
    char *hic_file, *paf_file, *hic_bfile, *report_file, *agp_file;
    int c, ncall, npali, read_len, report_only, n_threads, ret;
    uint64 *palis;
    int64 nhic, naln;
    char *pref_out;

    sys_init();
    realtime0 = realtime();

    f_type = NOSET;
    pref_out = "hapcure.out";
    hic_file = NULL;
    paf_file = NULL;
    read_len = 150;
    report_only = 0;
    n_threads = 1;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0) {
        if (c == 'o') pref_out = opt.arg;
        else if (c == 's') ec_conf.bin_size = parse_num(opt.arg);
        else if (c == 'f') ec_conf.min_frag = parse_num(opt.arg);
        else if (c == 'd') ec_conf.med_drop = atof(opt.arg);
        else if (c == 'r') ec_conf.rec_rate = atof(opt.arg);
        else if (c == 'p') ec_conf.p_thresh = atof(opt.arg);
        else if (c == 't') n_threads = atoi(opt.arg);
        else if (c == 'a') paf_file = opt.arg;
        else if (c == 301) {
            if (parse_file_type(opt.arg, &f_type)) {
                fprintf(stderr, "[E::%s] unknown file type: %s\n", __func__, opt.arg);
                return 1;
            }
        }
        else if (c == 302) read_len = atoi(opt.arg);
        else if (c == 303) report_only = 1;
        else if (c == 304) ec_conf.pal_size = parse_num(opt.arg);
        else if (c == 305) ec_conf.pal_dens = atof(opt.arg);
        else if (c == 'v') VERBOSE = atoi(opt.arg);
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

    if (argc - opt.ind < 1) {
        print_help(stderr);
        return 1;
    }
    if (argc - opt.ind < 2 && paf_file == NULL) {
        fprintf(stderr, "[E::%s] missing HiC binary file or self-alignment file\n", __func__);
        print_help(stderr);
        return 1;
    }
    if (argc - opt.ind < 2)
        fprintf(stderr, "[W::%s] missing HiC binary file - call palindromes only\n", __func__);
    else hic_file = argv[opt.ind + 1];

    if (ec_conf.min_frag < ec_conf.bin_size)
        ec_conf.min_frag = ec_conf.bin_size;
    
    if (read_len < 0 || ec_conf.bin_size <= 0 || 
        ec_conf.med_drop <= 0. || ec_conf.med_drop > 1. ||
        ec_conf.rec_rate <= 0. || ec_conf.rec_rate > 1. ||
        ec_conf.p_thresh <= 0. || ec_conf.p_thresh > 1. ||
        ec_conf.pal_dens <= 0. || ec_conf.pal_dens > 1. ||
        ec_conf.pal_size <= 0) {
        fprintf(stderr, "[E::%s] invalid correction parameters\n", __func__);
        return 1;
    }

    // build sequence dictionary
    dicts = make_sdict_from_index(argv[opt.ind], 0);

    // read hic file
    hic_bfile = write_binary_hic_data(hic_file, f_type, dicts, read_len, pref_out);
    if (hic_bfile && (ret = match_binary_file_sdict(hic_bfile, dicts))) {
        fprintf(stderr, "[E::%s] HiC BIN sequence dictionary does not match genome: %d\n", __func__, ret);
        free(hic_bfile);
        sd_destroy(dicts);
        return ret;
    }

    // call error correction positions from Hi-C
    nhic = 0;
    hics = read_hic_from_binary(hic_bfile, dicts, ec_conf.bin_size, 0, &nhic);
    ncall = 0;
    calls = ec_call_breaks(hics, nhic, dicts, &ncall);
    free(hics);
    
    // call error correction positions from palindromes
    naln = 0;
    alns = read_pafs(&paf_file, 1, dicts, NULL, 0, 1, &naln);
    npali = 0;
    palis = aln_find_palindromes(alns, naln, dicts, ec_conf.pal_size, ec_conf.pal_dens, n_threads, &npali);
    free(alns); // no fragments to free

    // merge hic and palindrome calls
    if (npali > 0) {
        MYREALLOC(calls, ncall+npali);
        MYBZERO(calls+ncall, npali);
        int i;
        for (i = 0; i < npali; i++) {
            calls[ncall+i].seq = palis[i] >> 32;
            calls[ncall+i].pos = (uint32) palis[i];
        }
        ncall += npali;
        calls = ec_merge_calls(calls, ncall, dicts, &ncall);
        free(palis);
    }

    // report results
    fprintf(stderr, "[M::%s] number errors called: %d\n", __func__, ncall);

    MYMALLOC(report_file, strlen(pref_out) + 12);
    sprintf(report_file, "%s.ec.tsv", pref_out);
    fo = fopen(report_file, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot write file %s\n", __func__, report_file);
        free(report_file);
        free(calls);
        free(hic_bfile);
        sd_destroy(dicts);
        return 1;
    }
    ec_write_report(calls, ncall, dicts, fo);
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
            free(calls);
            free(hic_bfile);
            sd_destroy(dicts);
            return 1;
        }
        ec_write_agp(calls, ncall, dicts, fo);
        fclose(fo);
        fprintf(stderr, "[M::%s] wrote corrected AGP: %s\n", __func__, agp_file);
        free(agp_file);
    }

    free(calls);
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