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
 * 16/04/26 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bamlite.h"

#include "ketopt.h"
#include "khash.h"
#include "kstring.h"
#include "kvec.h"
#include "hic.h"
#include "sdict.h"
#include "misc.h"
#include "version.h"

int VERBOSE = 0;

static double ht_realtime0;

/***********************************************************************************
 * hictools convert                                                                *
 *                                                                                 *
 * Convert a HiC alignment file (BAM/BED/PA5/ONE) to the binary BIN format for     *
 * fast access by downstream tools.                                                *
 ***********************************************************************************/

static ko_longopt_t convert_long_options[] = {
    { "read-length",    ko_required_argument, 'l' },
    { "file-type",      ko_required_argument, 't' },
    { "version",        ko_no_argument,       'V' },
    { "help",           ko_no_argument,       'h' },
    { 0, 0, 0 }
};

static void print_help_convert(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "Usage: hictools convert [options] <genome.fa.[fai|idx]> <hic.[bed|bam|pa5|1map|bin]> [...]\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "    -o STR         output file prefix [hictools.out]\n");
    fprintf(fp, "    -t STR         input file type BED|BAM|PA5|ONE|BIN\n");
    fprintf(fp, "    -l INT         read length (required for PA5 format) [150]\n");
    fprintf(fp, "    -v INT         verbose level [0]\n");
    fprintf(fp, "    -h, --help     print this help\n");
    fprintf(fp, "    -V, --version  show version number\n");
    fprintf(fp, "\n");
    fprintf(fp, "Example: hictools convert genome.fa.fai hic.bam\n");
    fprintf(fp, "         hictools convert genome.fa.fai hic.1.bam hic.2.bam hic.3.bam\n");
    fprintf(fp, "\n");
}

static int main_convert(int argc, char *argv[])
{
    const char *opt_str = "o:t:l:v:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, read_len;
    fileType_t f_type;
    char *pref_out, *hic_bfile;
    sdict_t *dicts;
    FILE *fp_help;

    ht_realtime0 = realtime();

    pref_out = "hictools.out";
    read_len = 150;
    f_type = NOSET;
    fp_help = stderr;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, convert_long_options)) >= 0) {
        if (c == 'o') {
            pref_out = opt.arg;
        } else if (c == 'l') {
            read_len = atoi(opt.arg);
        } else if (c == 't') {
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
        } else if (c == 'v') {
            VERBOSE = atoi(opt.arg);
        } else if (c == 'h') {
            fp_help = stdout;
        } else if (c == 'V') {
            puts(HICTOOLS_VERSION);
            return 0;
        } else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        } else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (argc - opt.ind < 2 || fp_help == stdout) {
        print_help_convert(fp_help);
        return fp_help == stdout ? 0 : 1;
    }

    if (read_len < 0) {
        fprintf(stderr, "[E::%s] invalid read length: %d\n", __func__, read_len);
        return 1;
    }

    dicts = make_sdict_from_index(argv[opt.ind], 0);

    hic_bfile = write_binary_hic_data_multi(argv + opt.ind + 1, argc - opt.ind - 1, f_type, dicts, read_len, pref_out);
    free(hic_bfile);
    sd_destroy(dicts);

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(1);
    }

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, HICTOOLS_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        for (int i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n",
                __func__, realtime() - ht_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}

/***********************************************************************************
 * hictools prepare                                                                *
 *                                                                                 *
 * Read a BIN file (produced by 'hictools convert') and a genome index (.idx),    *
 * optionally guided by an AGP scaffold layout, and build a binned contact matrix *
 * for plotting.  The genome is partitioned into equal-size bins; each pair of     *
 * bins that has at least one read pair linking them is written to the output.     *
 *                                                                                 *
 * Bin-size strategy: the user specifies '-n INT' (the desired number of bins,     *
 * i.e. the pixel resolution of the output matrix).  The tool derives bin size as  *
 *   bin_size = ceil(total_genome_length / n)                                      *
 * This is preferred over a fixed '-s' size because it directly expresses the      *
 * desired output resolution regardless of genome size.  '-s' is also accepted as  *
 * an explicit override; specifying both is an error.                              *
 ***********************************************************************************/

/*
 * Partitions the genome described by *dicts* (and optionally remapped via *agp*)
 * into bins of *bin_size* bp, reads all HiC pairs from *bin_file*, tallies
 * contacts per bin pair, and writes the sparse contact matrix to *fo*.
 *
 * Parameters:
 *   bin_file  – path to the ALNHBIN binary file
 *   idx_file  – path to the genome sequence index (.idx)
 *   agp_file  – path to AGP scaffold file, or NULL to work on raw contigs
 *   bin_size  – bp per bin if positive; -n_bins if negative (derive from n)
 *   min_qual  – minimum HiC mapping quality to include
 *   fo        – output stream (stdout or an opened file)
 *
 * Returns 0 on success, non-zero on error.
 */

static int make_contact_map(char *bin_file, char *idx_file, char *agp_file, int bin_size, int min_qual, FILE *fo)
{
    asm_dict_t *dicts;
    sdict_t *sdicts;
    hic_t *hics;
    int64 i, j, slen, nhic, cnts;

    sdicts = make_sdict_from_index(idx_file, 0);
    dicts = agp_file? make_asm_dict_from_agp(sdicts, agp_file, 1) : make_asm_dict_from_sdict(sdicts);

    /* total scaffold coordinate space = sequence bases + gap bases,
     * matching what sd_coordinate_conversion(count_gap=1) uses          */
    slen = 0;
    for (j = 0; j < dicts->n; j++)
        slen += (int64)(dicts->s[j].len + dicts->s[j].gap);

    if (bin_size <= 0)
        bin_size = 1 + (MAX(slen, 1) - 1) / (-bin_size);

    nhic = 0;
    hics = read_hic_from_binary_sd_conversion(bin_file, dicts, bin_size, (uint8)min_qual, &nhic);

    cnts = 0;
    for (i = 0; i < nhic; i++) cnts += hics[i].nhic;

    fprintf(fo, "# hictools version: %s\n", HICTOOLS_VERSION);
    fprintf(fo, "# IDX file: %s\n", idx_file);
    if (agp_file)
        fprintf(fo, "# AGP file: %s\n", agp_file);
    fprintf(fo, "# number sequences: %u\n", dicts->n);
    fprintf(fo, "# total length: %lld\n", slen);
    fprintf(fo, "# bin size: %d\n", bin_size);
    fprintf(fo, "# number hic pairs: %lld\n", cnts);
    fprintf(fo, "# number contacts: %lld\n", nhic);

    /* S lines: one per scaffold/contig in index order */
    for (j = 0; j < dicts->n; j++)
        fprintf(fo, "S\t%s\t%llu\n", dicts->s[j].name,
                (unsigned long long)(dicts->s[j].len + dicts->s[j].gap));

    /* B line: bin size in bp */
    fprintf(fo, "B\t%d\n", bin_size);

    /* C lines: sparse contact matrix entries */
    if (hics) {
        for (i = 0; i < nhic; i++)
            fprintf(fo, "C\t%d\t%d\t%d\t%d\t%d\n",
                    hics[i].aseq, hics[i].apos,
                    hics[i].bseq, hics[i].bpos,
                    hics[i].nhic);
        free(hics);
    }

    asm_destroy(dicts);
    sd_destroy(sdicts);
    return 0;
}

static ko_longopt_t prepare_long_options[] = {
    { "agp",        ko_required_argument, 'a' },
    { "n-bins",     ko_required_argument, 'n' },
    { "bin-size",   ko_required_argument, 's' },
    { "min-quality",ko_required_argument, 'q' },
    { "version",    ko_no_argument,       'V' },
    { "help",       ko_no_argument,       'h' },
    { 0, 0, 0 }
};

static void print_help_prepare(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "Usage: hictools prepare [options] <hic.bin> <genome.fa.[fai|idx]>\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "    -a FILE        AGP file mapping contigs to scaffolds (optional)\n");
    fprintf(fp, "    -n INT         number of bins (pixels) for the contact map [1000]\n");
    fprintf(fp, "                   bin size = ceil(genome_length / n); preferred over -s\n");
    fprintf(fp, "    -s INT         explicit bin size in bp (overrides -n)\n");
    fprintf(fp, "    -q INT         minimum mapping quality [0]\n");
    fprintf(fp, "    -o FILE        output file [stdout]\n");
    fprintf(fp, "    -v INT         verbose level [0]\n");
    fprintf(fp, "    -h, --help     print this help\n");
    fprintf(fp, "    -V, --version  show version number\n");
    fprintf(fp, "\n");
    fprintf(fp, "Example: hictools prepare -n 2000 hic.bin genome.fa.fai\n");
    fprintf(fp, "         hictools prepare -a scaffolds.agp -n 2000 hic.bin genome.fa.fai\n");
    fprintf(fp, "\n");
}

static int main_prepare(int argc, char *argv[])
{
    const char *opt_str = "a:n:s:q:o:v:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, min_qual;
    int n_bins, bin_size;
    char *agp_file, *out;
    char *bin_file, *idx_file;
    FILE *fp_help, *fo;

    ht_realtime0 = realtime();

    agp_file = NULL;
    out      = NULL;   /* NULL → stdout */
    n_bins   = 1000;
    bin_size = 0;    /* 0 = derive from n_bins */
    min_qual = 0;
    fp_help  = stderr;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, prepare_long_options)) >= 0) {
        if (c == 'a') {
            agp_file = opt.arg;
        } else if (c == 'n') {
            n_bins = atol(opt.arg);
        } else if (c == 's') {
            bin_size = atol(opt.arg);
        } else if (c == 'q') {
            min_qual = atoi(opt.arg);
        } else if (c == 'o') {
            out = opt.arg;
        } else if (c == 'v') {
            VERBOSE = atoi(opt.arg);
        } else if (c == 'h') {
            fp_help = stdout;
        } else if (c == 'V') {
            puts(HICTOOLS_VERSION);
            return 0;
        } else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        } else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (argc - opt.ind != 2 || fp_help == stdout) {
        print_help_prepare(fp_help);
        return fp_help == stdout ? 0 : 1;
    }

    /* sanity checks */
    if (n_bins > 0 && bin_size > 0) {
        fprintf(stderr, "[E::%s] -n and -s are mutually exclusive\n", __func__);
        return 1;
    }
    if (n_bins <= 0 && bin_size <= 0) {
        fprintf(stderr, "[E::%s] -n must be a positive integer\n", __func__);
        return 1;
    }
    if (min_qual < 0 || min_qual > 255) {
        fprintf(stderr, "[E::%s] invalid mapping quality threshold: %d\n", __func__, min_qual);
        return 1;
    }

    bin_file = argv[opt.ind];
    idx_file = argv[opt.ind + 1];

    /* open output file, or fall back to stdout */
    fo = out ? fopen(out, "w") : stdout;
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        return 1;
    }

    /* If the user gave -n, bin_size is resolved inside make_contact_map()  *
     * once the total genome length is known.  Pass n_bins as a negative     *
     * sentinel so the callee can distinguish the two modes.                 */
    if (bin_size == 0)
        bin_size = -n_bins;   /* negative = "derive from this many bins" */

    int ret = make_contact_map(bin_file, idx_file, agp_file, bin_size,
                               min_qual, fo);

    if (out)
        fclose(fo);

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(1);
    }

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, HICTOOLS_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        for (int i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n",
                __func__, realtime() - ht_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return ret;
}

KHASH_SET_INIT_STR(str)


static void hiclink_file_from_yahsbin(FILE *fp, uint64 pair_n, asm_dict_t *dict, uint8 mq, int scale, int count_gap, FILE *fo)
{
    uint32 m, i0, i1, s0, s1, q0, q1;
    uint64 p0, p1, pair_c, pair_u;
    uint8 *bbuff, *ebuff, *buffer;
    CC_ERR_t err0, err1;

    MYMALLOC(buffer, YAHSBIN_ENTRY_SIZE * FREAD_BUFF_SIZE);
    if (buffer == NULL)
        mem_alloc_error("buffer array");

    pair_c = pair_u = 0;
    while (pair_c < pair_n) {
        m = fread(buffer, sizeof(uint8) * YAHSBIN_ENTRY_SIZE, FREAD_BUFF_SIZE, fp);

        for (bbuff = buffer, ebuff = buffer + m * YAHSBIN_ENTRY_SIZE; 
            bbuff < ebuff && pair_c++ < pair_n; bbuff += YAHSBIN_ENTRY_SIZE) {
            if (bbuff[YAHSBIN_ENTRY_SIZE-1] < mq)
                continue;

            memcpy(&i0, bbuff,  sizeof(uint32));
            memcpy(&s0, bbuff + sizeof(uint32)*1, sizeof(uint32));
            memcpy(&i1, bbuff + sizeof(uint32)*2, sizeof(uint32));
            memcpy(&s1, bbuff + sizeof(uint32)*3, sizeof(uint32));

            err0 = sd_coordinate_conversion(dict, i0, s0, &q0, &p0, count_gap);
            err1 = sd_coordinate_conversion(dict, i1, s1, &q1, &p1, count_gap);
            
            if (err0 != CC_SUCCESS || err1 != CC_SUCCESS) {
                // fprintf(stderr, "[W::%s] sequence not found \n", __func__);
                ++pair_u;
            } else {
                if (strcmp(dict->s[q0].name, dict->s[q1].name) <= 0)
                    fprintf(fo, "0\t%s\t%llu\t0\t1\t%s\t%llu\t1\n", dict->s[q0].name, p0 >> scale, dict->s[q1].name, p1 >> scale);
                else
                    fprintf(fo, "0\t%s\t%llu\t1\t1\t%s\t%llu\t0\n", dict->s[q1].name, p1 >> scale, dict->s[q0].name, p0 >> scale);
            }
        }
    }

    fprintf(stderr, "[M::%s] %llu read pairs processed: %llu unmapped\n", __func__, pair_c, pair_u);

    free(buffer);
}


static void hiclink_file_from_alnhbin(FILE *fp, uint64 pair_n, asm_dict_t *dict, uint8 mq, int scale, int count_gap, FILE *fo)
{
    uint32 m, i0, i1, s0, s1, e0, e1, q0, q1;
    uint64 p0, p1, pair_c, pair_u;
    uint8 *bbuff, *ebuff, *buffer;
    CC_ERR_t err0, err1;

    MYMALLOC(buffer, ALNHBIN_ENTRY_SIZE * FREAD_BUFF_SIZE);
    if (buffer == NULL)
        mem_alloc_error("buffer array");

    pair_c = pair_u = 0;
    while (pair_c < pair_n) {
        m = fread(buffer, sizeof(uint8) * ALNHBIN_ENTRY_SIZE, FREAD_BUFF_SIZE, fp);

        for (bbuff = buffer, ebuff = buffer + m * ALNHBIN_ENTRY_SIZE; 
            bbuff < ebuff && pair_c++ < pair_n; bbuff += ALNHBIN_ENTRY_SIZE) {
            if (bbuff[ALNHBIN_ENTRY_SIZE-1] < mq)
                continue;

            memcpy(&i0, bbuff,  sizeof(uint32));
            memcpy(&s0, bbuff + sizeof(uint32)*1, sizeof(uint32));
            memcpy(&e0, bbuff + sizeof(uint32)*2, sizeof(uint32));
            memcpy(&i1, bbuff + sizeof(uint32)*3, sizeof(uint32));
            memcpy(&s1, bbuff + sizeof(uint32)*4, sizeof(uint32));
            memcpy(&e1, bbuff + sizeof(uint32)*5, sizeof(uint32));

            s0 = (s0>>1) + (e0>>1) + ((s0 & e0) & 1);
            s1 = (s1>>1) + (e1>>1) + ((s1 & e1) & 1);

            err0 = sd_coordinate_conversion(dict, i0, s0, &q0, &p0, count_gap);
            err1 = sd_coordinate_conversion(dict, i1, s1, &q1, &p1, count_gap);
            
            if (err0 != CC_SUCCESS || err1 != CC_SUCCESS) {
                // fprintf(stderr, "[W::%s] sequence not found \n", __func__);
                ++pair_u;
            } else {
                if (strcmp(dict->s[q0].name, dict->s[q1].name) <= 0)
                    fprintf(fo, "0\t%s\t%llu\t0\t1\t%s\t%llu\t1\n", dict->s[q0].name, p0 >> scale, dict->s[q1].name, p1 >> scale);
                else
                    fprintf(fo, "0\t%s\t%llu\t1\t1\t%s\t%llu\t0\n", dict->s[q1].name, p1 >> scale, dict->s[q0].name, p0 >> scale);
            }
        }
    }

    fprintf(stderr, "[M::%s] %llu read pairs processed: %llu unmapped\n", __func__, pair_c, pair_u);

    free(buffer);
}

static int hiclink_file_from_bin(char *f, char *agp, char *fai, uint8 mq, int scale, int count_gap, FILE *fo)
{
    FILE *fp;
    uint64 pair_n;
    int64 magic_number;
    int ret;
    
    sdict_t *sdict = make_sdict_from_index(fai, 0);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp, 1) : make_asm_dict_from_sdict(sdict);
    
    // check BIN file header consistency
    if ((ret = match_binary_file_sdict(f, sdict))) {
        fprintf(stderr, "[E::%s] Not a valid BIN file or BIN file header does not match sequence dictionary: %d\n", __func__, ret);
        fprintf(stderr, "[E::%s] Make sure no contig length threshold (YaHS option '-l') applied for the BIN file\n", __func__);
        fprintf(stderr, "[E::%s] Consider using a BAM or BED file instead\n", __func__);
        exit(EXIT_FAILURE);
    }

    fp = fopen(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    
    if (fread(&magic_number, sizeof(int64), 1, fp) != 1)
        bin_fread_error();
    
    binary_fseek_skip_sdict(fp);
    
    if (fread(&pair_n, sizeof(uint64), 1, fp) != 1)
        bin_fread_error();

    if (is_valid_alnh_bin_header(magic_number))
        hiclink_file_from_alnhbin(fp, pair_n, dict, mq, scale, count_gap, fo);
    else if (is_valid_yahs_bin_header(magic_number))
        hiclink_file_from_yahsbin(fp, pair_n, dict, mq, scale, count_gap, fo);
    else {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        exit(EXIT_FAILURE);
    }

    fclose(fp);
    asm_destroy(dict);
    sd_destroy(sdict);

    return 0;
}

static inline int parse_line(char *s, char **f, int m)
{
    int n = 0;
    while (*s) {
        while (isspace(*s)) s++;
        if (!*s) break;
        f[n++] = s;
        while (*s && !isspace(*s)) s++;
        if (*s) *s++ = '\0';
        if (n >= m) break;
    }
    return n;
}

static int hiclink_file_from_bed(char *f, char *agp, char *fai, uint8 mq, int scale, int count_gap, FILE *fo)
{
    iostream_t *fp;
    char *line, *fields[6];
    kstring_t _cname = {0, 0, 0}, *cname = &_cname;
    kstring_t _rname = {0, 0, 0}, *rname = &_rname;
    uint32 s0, s1, e0, e1, i0, i1;
    uint64 p0, p1, rec_c, pair_c;
    int8 buff;
    CC_ERR_t err0, err1;

    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    sdict_t *sdict = make_sdict_from_index(fai, 0);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp, 1) : make_asm_dict_from_sdict(sdict);

    fp = iostream_open(f);
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    rec_c = pair_c = 0;
    buff = 0;
    while ((line = iostream_getline(fp)) != NULL) {
        if (is_empty_line(line)|| parse_line(line, fields, 5) < 5)
            continue;
        
        ++rec_c;
        
        if (buff == 0) {
            if (atoi(fields[4]) >= mq) {
                cname->l = 0, kputs(fields[0], cname);
                s0 = atoi(fields[1]);
                e0 = atoi(fields[2]);
                rname->l = 0, kputs(fields[3], rname);
                buff = 1;
            }
        } else if (buff == 1) {
            if (is_read_pair(rname->s, fields[3])) {
                if (atoi(fields[4]) >= mq) {
                    s1 = atoi(fields[1]);
                    e1 = atoi(fields[2]);
                    
                    err0 = sd_coordinate_conversion(dict, sd_get(sdict, cname->s),  s0 / 2 + e0 / 2 + (s0 & 1 && e0 & 1), &i0, &p0, count_gap);
                    err1 = sd_coordinate_conversion(dict, sd_get(sdict, fields[0]), s1 / 2 + e1 / 2 + (s1 & 1 && e1 & 1), &i1, &p1, count_gap);
                    
                    if (err0 != CC_SUCCESS) {
                        if (err0 == SEQ_NOT_FOUND) {
                            k = kh_put(str, hmseq, cname->s, &absent);
                            if (absent) {
                                kh_key(hmseq, k) = strdup(cname->s);
                                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname->s);
                            }
                        } else if (err0 == POS_NOT_IN_RANGE) {
                            fprintf(stderr, "[W::%s] sequence position \"%s:%u\" not in range \n", __func__, 
                                    cname->s, s0 / 2 + e0 / 2 + (s0 & 1 && e0 & 1));
                        }
                    } else if (err1 != CC_SUCCESS) {
                        if (err1 == SEQ_NOT_FOUND) {
                            k = kh_put(str, hmseq, fields[0], &absent);
                            if (absent) {
                                kh_key(hmseq, k) = strdup(fields[0]);
                                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, fields[0]);
                            }
                        } else if (err1 == POS_NOT_IN_RANGE) {
                            fprintf(stderr, "[W::%s] sequence position \"%s:%u\" not in range \n", __func__,
                                    fields[0], s1 / 2 + e1 / 2 + (s1 & 1 && e1 & 1));
                        }
                    } else {
                        if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                            fprintf(fo, "0\t%s\t%llu\t0\t1\t%s\t%llu\t1\n", dict->s[i0].name, p0 >> scale, dict->s[i1].name, p1 >> scale);
                        else
                            fprintf(fo, "0\t%s\t%llu\t1\t1\t%s\t%llu\t0\n", dict->s[i1].name, p1 >> scale, dict->s[i0].name, p0 >> scale);
                        
                        ++pair_c;
                    }
                }
                buff = 0;
            } else {
                if (atoi(fields[4]) >= mq) {
                    cname->l = 0, kputs(fields[0], cname);
                    s0 = atoi(fields[1]);
                    e0 = atoi(fields[2]);
                    rname->l = 0, kputs(fields[3], rname);
                    buff = 1;
                }
            }
        }

        if (rec_c % REPORT_EPOC == 0)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs\n", __func__, rec_c / 1000000, pair_c);
    }

    fprintf(stderr, "[M::%s] %llu records processed, %llu read pairs\n", __func__, rec_c, pair_c);
    
    free(cname->s);
    free(rname->s);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    iostream_close(fp);
    asm_destroy(dict);
    sd_destroy(sdict);

    return 0;
}

static int hiclink_file_from_pa5(char *f, char *agp, char *fai, int8 mq, int scale, int count_gap, FILE *fo)
{
    iostream_t *fp;
    char *line, *fields[7];
    uint32 x0, x1, i0, i1;
    uint64 p0, p1, rec_c, pair_c;
    CC_ERR_t err0, err1;

    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    sdict_t *sdict = make_sdict_from_index(fai, 0);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp, 1) : make_asm_dict_from_sdict(sdict);

    fp = iostream_open(f);
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    rec_c = pair_c = 0;
    while ((line = iostream_getline(fp)) != NULL) {
        if (is_empty_line(line) || parse_line(line, fields, 7) < 7)
            continue;
        
        ++rec_c;

        if (atoi(fields[5]) >= mq && atoi(fields[6]) >= mq) {
            x0 = atoi(fields[2]);
            x1 = atoi(fields[4]);

            err0 = sd_coordinate_conversion(dict, sd_get(sdict, fields[1]), x0, &i0, &p0, count_gap);
            err1 = sd_coordinate_conversion(dict, sd_get(sdict, fields[3]), x1, &i1, &p1, count_gap);

            if (err0 != CC_SUCCESS) {
                if (err0 == SEQ_NOT_FOUND) {
                    k = kh_put(str, hmseq, fields[1], &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(fields[1]);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, fields[1]);
                    }
                } else if (err0 == POS_NOT_IN_RANGE) {
                    fprintf(stderr, "[W::%s] sequence position \"%s:%u\" not in range \n", __func__, fields[1], x0);
                }
            } else if (err1 != CC_SUCCESS) {
                if (err1 == SEQ_NOT_FOUND) {
                    k = kh_put(str, hmseq, fields[3], &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(fields[3]);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, fields[3]);
                    }
                } else if (err1 == POS_NOT_IN_RANGE) {
                    fprintf(stderr, "[W::%s] sequence position \"%s:%u\" not in range \n", __func__, fields[3], x1);
                }
            } else {
                if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                    fprintf(fo, "0\t%s\t%llu\t0\t1\t%s\t%llu\t1\n", dict->s[i0].name, p0 >> scale, dict->s[i1].name, p1 >> scale);
                else
                    fprintf(fo, "0\t%s\t%llu\t1\t1\t%s\t%llu\t0\n", dict->s[i1].name, p1 >> scale, dict->s[i0].name, p0 >> scale);

                ++pair_c;
            }
        }

        if (rec_c % REPORT_EPOC == 0)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs \n", __func__, rec_c / 1000000, pair_c);
    }

    fprintf(stderr, "[M::%s] %llu records processed, %llu read pairs\n", __func__, rec_c, pair_c);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    iostream_close(fp);
    asm_destroy(dict);
    sd_destroy(sdict);

    return 0;
}

static inline char *parse_bam_rec(bam1_t *b, bam_header_t *h, uint8 q, int32 *s, int32 *e, kstring_t *rname)
{
    // 0x4 0x100 0x200 0x400 0x800
    if ((b->core.flag & 0xF04) || (!b->core.flag) || b->core.qual < q)
        return 0;
    rname->l = 0;
    kputs(bam1_qname(b), rname);
    *s = b->core.pos;
    *e = get_target_end(b);

    return h->target_name[b->core.tid];
}

static inline int parse_bam_rec1(bam1_t *b, bam_header_t *h, char **cname0, int32 *s0, char **cname1, int32 *s1)
{
    // 0x4 0x8 0x40 0x100 0x200 0x400 0x800
    if ((b->core.flag & 0xF4C) || (!b->core.flag))
        return 1;
    *cname0 = h->target_name[b->core.tid];
    *s0 = b->core.pos;
    *cname1 = h->target_name[b->core.mtid];
    *s1 = b->core.mpos;

    return 0;
}

static int hiclink_file_from_bam(char *f, char *agp, char *fai, uint8 mq, int scale, int count_gap, FILE *fo)
{
    bamFile fp;
    bam_header_t *h;
    bam1_t *b;
    kstring_t _rname0 = {0, 0, 0}, *rname0 = &_rname0;
    kstring_t _rname1 = {0, 0, 0}, *rname1 = &_rname1;
    kstring_t *tmp;
    char *cname0, *cname1;
    int32 s0, s1, e0, e1;
    uint32 i0, i1;
    uint64 p0, p1, rec_c, pair_c;
    int8 buff;
    enum bam_sort_order so;
    CC_ERR_t err0, err1;

    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    sdict_t *sdict = make_sdict_from_index(fai, 0);
    asm_dict_t *dict = agp? make_asm_dict_from_agp(sdict, agp, 1) : make_asm_dict_from_sdict(sdict);

    fp = bam_open(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    
    h = bam_header_read(fp);
    b = bam_init1();
    so = bam_hrecs_sort_order(h);

    ks_resize(rname0, 256);
    ks_resize(rname1, 256);
    cname0 = cname1 = 0;
    s0 = s1 = e0 = e1 = 0;
    i0 = i1 = 0;
    p0 = p1 = 0;
    rec_c = pair_c = 0;
    buff = 0;    
    
    if (so == ORDER_NAME) {
        // sorted by read names
        while (bam_read1(fp, b) >= 0 ) {
            ++rec_c;

            if (buff == 0) {
                cname0 = parse_bam_rec(b, h, mq, &s0, &e0, rname0);
                if (!cname0)
                    continue;
                buff = 1;
            } else if (buff == 1) {
                cname1 = parse_bam_rec(b, h, mq, &s1, &e1, rname1);
                if (!cname1)
                    continue;
                if (is_read_pair(rname0->s, rname1->s)) {

                    err0 = sd_coordinate_conversion(dict, sd_get(sdict, cname0), s0 / 2 + e0 / 2 + (s0 & 1 && e0 & 1), &i0, &p0, count_gap);
                    err1 = sd_coordinate_conversion(dict, sd_get(sdict, cname1), s1 / 2 + e1 / 2 + (s1 & 1 && e1 & 1), &i1, &p1, count_gap);

                    if (err0 != CC_SUCCESS) {
                        if (err0 == SEQ_NOT_FOUND) {
                            k = kh_put(str, hmseq, cname0, &absent);
                            if (absent) {
                                kh_key(hmseq, k) = strdup(cname0);
                                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                            }
                        } else if (err0 == POS_NOT_IN_RANGE) {
                            fprintf(stderr, "[W::%s] sequence position \"%s:%u\" not in range \n", __func__,
                                    cname0, s0 / 2 + e0 / 2 + (s0 & 1 && e0 & 1));
                        }
                    } else if (err1 != CC_SUCCESS) {
                        if (err1 == SEQ_NOT_FOUND) {
                            k = kh_put(str, hmseq, cname1, &absent);
                            if (absent) {
                                kh_key(hmseq, k) = strdup(cname1);
                                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname1);
                            }
                        } else if (err1 == POS_NOT_IN_RANGE) {
                            fprintf(stderr, "[W::%s] sequence position \"%s:%u\" not in range \n", __func__,
                                    cname1, s1 / 2 + e1 / 2 + (s1 & 1 && e1 & 1));
                        }
                    } else {
                        if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                            fprintf(fo, "0\t%s\t%llu\t0\t1\t%s\t%llu\t1\n", dict->s[i0].name, p0 >> scale, dict->s[i1].name, p1 >> scale);
                        else
                            fprintf(fo, "0\t%s\t%llu\t1\t1\t%s\t%llu\t0\n", dict->s[i1].name, p1 >> scale, dict->s[i0].name, p0 >> scale);

                        ++pair_c;
                    }

                    buff = 0;
                } else {
                    cname0 = cname1;
                    s0 = s1;
                    e0 = e1;
                    tmp = rname0;
                    rname0 = rname1;
                    rname1 = tmp;
                    
                    buff = 1;
                }
            }

            if (rec_c % REPORT_EPOC == 0)
                fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs \n", __func__, rec_c / 1000000, pair_c);
        }
    } else {
        // sorted by coordinates or others
        if (mq > 0)
            fprintf(stderr, "[W::%s] BAM file is not sorted by read name. Filtering by mapping quality %hhu suppressed \n", __func__, mq);

        while (bam_read1(fp, b) >= 0 ) {
            ++rec_c;
            
            if(parse_bam_rec1(b, h, &cname0, &s0, &cname1, &s1))
                continue;

            err0 = sd_coordinate_conversion(dict, sd_get(sdict, cname0), s0, &i0, &p0, count_gap);
            err1 = sd_coordinate_conversion(dict, sd_get(sdict, cname1), s1, &i1, &p1, count_gap);

            if (err0 != CC_SUCCESS) {
                if (err0 == SEQ_NOT_FOUND) {
                    k = kh_put(str, hmseq, cname0, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname0);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                    }
                } else if (err0 == POS_NOT_IN_RANGE) {
                    fprintf(stderr, "[W::%s] sequence position \"%s:%u\" not in range \n", __func__, cname0, s0);
                }
            } else if (err1 != CC_SUCCESS) {
                if (err1 == SEQ_NOT_FOUND) {
                    k = kh_put(str, hmseq, cname1, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname1);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname1);
                    }
                } else if (err1 == POS_NOT_IN_RANGE) {
                    fprintf(stderr, "[W::%s] sequence position \"%s:%u\" not in range \n", __func__, cname1, s1);
                }
            } else {
                if (strcmp(dict->s[i0].name, dict->s[i1].name) <= 0)
                    fprintf(fo, "0\t%s\t%llu\t0\t1\t%s\t%llu\t1\n", dict->s[i0].name, p0 >> scale, dict->s[i1].name, p1 >> scale);
                else
                    fprintf(fo, "0\t%s\t%llu\t1\t1\t%s\t%llu\t0\n", dict->s[i1].name, p1 >> scale, dict->s[i0].name, p0 >> scale);
        
                ++pair_c;
            }

            if (rec_c % REPORT_EPOC == 0)
                fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs \n", __func__, rec_c / 1000000, pair_c);
        }
    }

    fprintf(stderr, "[M::%s] %llu records processed, %llu read pairs\n", __func__, rec_c, pair_c);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    free(rname0->s);
    free(rname1->s);

    bam_destroy1(b);
    bam_header_destroy(h);
    bam_close(fp);
    asm_destroy(dict);
    sd_destroy(sdict);

    return 0;
}

static uint64 linear_scale(uint64 g, int *scale, uint64 max_g)
{
    int s;
    s = 0;
    while (g > max_g) {
        ++s;
        g >>= 1;
    }
    
    *scale = s;
    return g;
}

static uint64 assembly_annotation(const char *f, sdict_t *sdict, const char *out_agp, const char *out_annot, 
        const char *out_lift, int *scale, uint64 max_s, uint64 *g)
{
    uint32 i, j, c, n, m;
    int *seqs;
    FILE *fo_agp, *fo_annot, *fo_lift;
    uint64 genome_size, scaled_gs;
    asm_dict_t *dict;
    sd_seg_t seg;

    fo_agp = fopen(out_agp, "w");
    fo_annot = fopen(out_annot, "w");
    fo_lift = fopen(out_lift, "w");
    
    dict = make_asm_dict_from_agp(sdict, f, 1);
    n = dict->u + dict->n;
    seqs = (int *) calloc(n, sizeof(int));
    genome_size = 0;
    c = j = 0;
    m = dict->s[j++].n;
    for (i = 0; i < dict->u; ++i) {
        seg = dict->seg[i];
        fprintf(fo_agp, "assembly\t%llu\t%llu\t%u\t%s\t%s\t%u\t%u\t%s\n", 
                genome_size + 1, genome_size + seg.y, i + 1, agp_component_type_val(seg.t),
                sdict->s[seg.c>>1].name, seg.x + 1, seg.x + seg.y, agp_orientation_val(seg.r));
        fprintf(fo_annot, ">ctg%08u.1 %u %u\n", i + 1, i + 1, seg.y);
        fprintf(fo_lift, "ctg%08u.1\t%u\t%u\t%u\t%s\t%s\t%u\t%u\t%s\n",
                i + 1, 1, seg.y, 1, agp_component_type_val(seg.t), sdict->s[seg.c>>1].name, 
                seg.x + 1, seg.x + seg.y, agp_orientation_val(AGP_OT_PLUS));
        if (i == m) {
            m += dict->s[j++].n;
            ++c;
        }
        seqs[c++] = (int) (i + 1) * (seg.r == AGP_OT_MINUS? -1 : 1);
        genome_size += seg.y;
    }

    for (i = 0; i < n; ++i) {
        // seqs[n - 1] is always 0
        if (seqs[i]) {
            fprintf(fo_annot, "%d", seqs[i]);
            fputc(seqs[i + 1]? ' ' : '\n', fo_annot);
        }
    }
    
    fclose(fo_agp);
    fclose(fo_annot);
    fclose(fo_lift);
    free(seqs);

    asm_destroy(dict);

    scaled_gs = linear_scale(genome_size, scale, max_s);
    *g = genome_size;
    
    return scaled_gs;
}

uint64 assembly_scale_max_seq(asm_dict_t *dict, int *scale, uint64 max_s, uint64 *g)
{
    uint32 i;
    uint64 s;
    sd_aseq_t seq;
 
    s = 0;
    for (i = 0; i < dict->n; ++i) {
        seq = dict->s[i];
        s = MAX(s, seq.len + seq.gap);
    }
    
    *g = s;

    return linear_scale(s, scale, max_s);
}

static void print_help_hiclink(FILE *fp_help)
{
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Usage: hictools hiclink [options] <hic.[bed|bam|pa5|1map|bin]> <scaffolds.agp> <contigs.fa.[fai|idx]>\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help, "    -a             preprocess for assembly mode\n");
    fprintf(fp_help, "    -q INT         minimum mapping quality [10]\n");
    fprintf(fp_help, "    -o STR         output file prefix (required for '-a' mode) [stdout]\n");
    fprintf(fp_help, "    -f STR         input file type BED|BAM|PA5|ONE|BIN, file name extension is ignored\n");
    fprintf(fp_help, "    -h, --help     print this help\n");
    fprintf(fp_help, "    -V, --version  show version number\n");
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Example: hictools hiclink -o output hic.bam scaffolds.agp contigs.fa.fai\n");
    fprintf(fp_help, "\n");
}

static ko_longopt_t hiclink_long_options[] = {
    { "file-type",      ko_required_argument, 'f' },
    { "version",        ko_no_argument,       'V' },
    { "help",           ko_no_argument,       'h' },
    { 0, 0, 0 }
};

static int main_hiclink(int argc, char *argv[])
{
    FILE *fo;
    char *fai, *agp, *agp1, *link_file, *out, *out1, *annot, *lift, *ext1, *ext2;
    int mq, asm_mode;;
    fileType_t f_type;

    liftrlimit();
    ht_realtime0 = realtime();

    const char *opt_str = "q:ao:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, ret;
    FILE *fp_help = stderr;
    fai = agp = agp1 = link_file = out = out1 = annot = lift = 0;
    mq = 10;
    asm_mode = 0;
    f_type = NOSET;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, hiclink_long_options)) >= 0) {
        if (c == 'o') {
            out = opt.arg;
        } else if (c == 'q') {
            mq = atoi(opt.arg);
        } else if (c == 'a') {
            asm_mode = 1;
        } else if (c == 'f') {
            if (strcasecmp(opt.arg, "BED") == 0)
                f_type = BED;
            else if (strcasecmp(opt.arg, "BAM") == 0)
                f_type = BAM;
            else if (strcasecmp(opt.arg, "BIN") == 0)
                f_type = BIN;
            else if (strcasecmp(opt.arg, "PA5") == 0)
                f_type = PA5;
            else {
                fprintf(stderr, "[E::%s] unknown file type: \"%s\"\n", __func__, opt.arg);
                return 1;
            }
        } else if (c == 'h') {
            fp_help = stdout;   
        } else if (c == 'V') {
            puts(HICTOOLS_VERSION);
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
        print_help_hiclink(stdout);
        return 0;
    }

    if (asm_mode && !out) {
        fprintf(stderr, "[E::%s] missing input: -o option is required for assembly mode (-a)\n", __func__);
        return 1;
    }

    if (argc - opt.ind < 3) {
        fprintf(stderr, "[E::%s] missing input: three positional options required\n", __func__);
        print_help_hiclink(stderr);
        return 1;
    }
    
    if (mq < 0 || mq > 255) {
        fprintf(stderr, "[E::%s] invalid mapping quality threshold: %d\n", __func__, mq);
        return 1;
    }

    uint8 mq8;
    mq8 = (uint8) mq;

    link_file = argv[opt.ind];
    agp = argv[opt.ind + 1];
    fai = argv[opt.ind + 2];

    if (f_type == NOSET) {
        ext1 = strlen(link_file) >= 4? (link_file + strlen(link_file) - 4) : NULL;
        ext2 = strlen(link_file) >= 7? (link_file + strlen(link_file) - 7) : NULL;
        if (ext1 && !strcasecmp(ext1, ".bam")) f_type = BAM;
        else if (ext1 && !strcasecmp(ext1, ".bin")) f_type = BIN;
        else if ((ext1 && !strcasecmp(ext1, ".bed")) || (ext2 && !strcasecmp(ext2, ".bed.gz"))) f_type = BED;
        else if ((ext1 && !strcasecmp(ext1, ".pa5")) || (ext2 && !strcasecmp(ext2, ".pa5.gz"))) f_type = PA5;
        else {
            fprintf(stderr, "[E::%s] unknown link file format. File extension .bam, .bed, .pa5 or .bin or --file-type is expected\n", __func__);
            exit(EXIT_FAILURE);
        }
    }
    
    if (f_type == BIN && (strcmp(link_file, "-") == 0 || *link_file == '<')) {
        fprintf(stderr, "[E::%s] BIN file format from STDIN is not supported\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (out) {
        out1 = (char *) malloc(strlen(out) + 35);
        sprintf(out1, "%s.txt", out);
    }

    fo = out1 == 0? stdout : fopen(out1, "w");
    if (fo == 0) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    ret = 0;
    
    sdict_t *sdict;
    asm_dict_t *dict;
    int scale;
    uint64 max_s, scaled_s;
    
    sdict = make_sdict_from_index(fai, 0);
    scale = 0;
    max_s = scaled_s = 0;
    agp1 = (char *) malloc(MAX(strlen(agp), out? strlen(out) : 0) + 35);
    if (asm_mode) {
        annot = (char *) malloc(strlen(out) + 35);
        lift = (char *) malloc(strlen(out) + 35);
        sprintf(agp1, "%s.assembly.agp", out);
        sprintf(annot, "%s.assembly", out);
        sprintf(lift, "%s.liftover.agp", out);
        scaled_s = assembly_annotation(agp, sdict, agp1, annot, lift, &scale, (uint64) INT_MAX, &max_s);
        dict = make_asm_dict_from_agp(sdict, agp1, 1);
    } else {
        sprintf(agp1, "%s", agp);
        dict = make_asm_dict_from_agp(sdict, agp1, 1);
        scaled_s = assembly_scale_max_seq(dict, &scale, (uint64) INT_MAX, &max_s);
    }

    if (f_type == BAM) {
        fprintf(stderr, "[M::%s] generate hic links from BAM file %s\n", __func__, link_file);
        ret = hiclink_file_from_bam(link_file, agp1, fai, mq8, scale, !asm_mode, fo);
    } else if (f_type == BED) {
        fprintf(stderr, "[M::%s] generate hic links from BED file %s\n", __func__, link_file);
        ret = hiclink_file_from_bed(link_file, agp1, fai, mq8, scale, !asm_mode, fo);
    } else if (f_type == PA5) {
        fprintf(stderr, "[M::%s] generate hic links from PA5 file %s\n", __func__, link_file);
        ret = hiclink_file_from_pa5(link_file, agp1, fai, mq8, scale, !asm_mode, fo);
    } else if (f_type == BIN) {
        fprintf(stderr, "[M::%s] generate hic links from BIN file %s\n", __func__, link_file);
        ret = hiclink_file_from_bin(link_file, agp1, fai, mq8, scale, !asm_mode, fo);
    }

    if (asm_mode) {
        fprintf(stderr, "[M::%s] genome size: %llu\n", __func__, max_s);
        fprintf(stderr, "[M::%s] scale factor: %d\n", __func__, 1 << scale);
        fprintf(stderr, "[M::%s] chromosome sizes for juicer_tools pre -\n", __func__);
        fprintf(stderr, "PRE_C_SIZE: assembly %llu\n", scaled_s);
        fprintf(stderr, "[M::%s] JUICER_PRE CMD: java -Xmx36G -jar ${juicer_tools} pre %s %s.hic <(echo \"assembly %llu\")\n", 
                __func__, out1, out, scaled_s);
    } else {
        if (scale) {
            fprintf(stderr, "[W::%s] maximum scaffold length exceeds %d (=%llu)\n", __func__, INT_MAX, max_s);
            fprintf(stderr, "[W::%s] using scale factor: %d\n", __func__, 1 << scale);
        }
        fprintf(stderr, "[M::%s] chromosome sizes for juicer_tools pre -\n", __func__);
        uint32 i;
        sd_aseq_t seq;
        for (i = 0; i < dict->n; ++i) {
            seq = dict->s[i];
            fprintf(stderr, "PRE_C_SIZE: %s %llu\n", seq.name, (seq.len + seq.gap) >> scale);
        }
    }

    asm_destroy(dict);
    sd_destroy(sdict);

    fclose(fo);
    free(out1);
    free(agp1);
    free(annot);
    free(lift);

    fprintf(stderr, "[M::%s] Version: %s\n", __func__, HICTOOLS_VERSION);
    fprintf(stderr, "[M::%s] CMD: hictools", __func__);
    int i;
    for (i = 0; i < argc; ++i)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, 
            realtime() - ht_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);

    return ret;
}

/***********************************************************************************
 * hictools agplift                                                                *
 *                                                                                 *
 * Convert a Juicebox-curated 'review.assembly' file, together with the liftover   *
 * AGP produced during 'hiclink -a', back into a final scaffold AGP, optionally    *
 * also writing the corresponding scaffold FASTA.                                  *
 ***********************************************************************************/

static int assembly_to_agp(char *assembly, char *lift, sdict_t *sdict, FILE *fo)
{
    asm_dict_t *dict;
    iostream_t *fp;
    char *line;
    char cname[1024], c0[1024], c1[1024], *cstr;
    int32 cid, sid, fid, sign;
    uint32 i, clen, mlen, s, p, *coords;
    int64 slen;
    size_t m;
    kvec_t(int) segs;

    dict = make_asm_dict_from_agp(sdict, lift, 1);

    m = 4;
    coords = (uint32 *) malloc(sizeof(uint32) * m * 3);
    c0[0] = c1[0] = '\0';
    mlen = 0;
    sid = 0;
    kv_init(segs);

    fp = iostream_open(assembly);
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, assembly);
        exit(EXIT_FAILURE);
    }
    while ((line = iostream_getline(fp)) != NULL) {
        if (is_empty_line(line) || line[0] != '>')
            continue;

        sscanf(line, "%s %d %u", cname, &cid, &clen);
        cstr = strstr(cname, ":::");
        if (cstr != NULL)
            cname[cstr - cname] = '\0';
        strcpy(c1, cname + 1);

        if (!strcmp(c0, c1)) {
            mlen += clen;
        } else {
            mlen = clen;
            strcpy(c0, c1);
        }

        if (sd_coordinate_rev_conversion(dict, asm_sd_get(dict, c1), mlen - clen, &s, &p, 0) != CC_SUCCESS) {
            fprintf(stderr, "[E::%s] coordinates conversion error %s %u\n", __func__, c1, mlen - clen);
            exit(EXIT_FAILURE);
        }

        if (clen == 0)
            fprintf(stderr, "[W::%s] segment of length zero in line: %s\n", __func__, line);

        if (cid > (int32) m) {
            m <<= 1;
            coords = (uint32 *) realloc(coords, sizeof(uint32) * m * 3);
        }

        cid -= 1;
        coords[cid * 3] = s;
        coords[cid * 3 + 1] = p + 1; // 0-based to 1-based coordinates
        coords[cid * 3 + 2] = clen;
    }
    iostream_close(fp);

    fp = iostream_open(assembly);
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, assembly);
        exit(EXIT_FAILURE);
    }
    while ((line = iostream_getline(fp)) != NULL) {
        if (is_empty_line(line) || line[0] == '>')
            continue;

        segs.n = 0;
        char *eptr, *fptr;
        cid = strtol(line, &eptr, 10);
        if (coords[(abs(cid) - 1) * 3 + 2] > 0)
            kv_push(int, segs, cid);
        while (*eptr && *eptr != '\n') {
            cid = strtol(eptr + 1, &fptr, 10);
            if (coords[(abs(cid) - 1) * 3 + 2] > 0)
                kv_push(int, segs, cid);
            eptr = fptr;
        }

        if (segs.n == 0) continue;

        ++sid;
        fid = 0;
        slen = 0;
        for (i = 0; i < (uint32) segs.n; ++i) {
            cid = segs.a[i];
            sign = cid > 0;
            cid = abs(cid) - 1;
            fprintf(fo, "scaffold_%d\t%lld\t%lld\t%d\t%s\t%s\t%d\t%d\t%s\n", sid, slen + 1,
                    slen + coords[cid * 3 + 2], ++fid, agp_component_type_val(DEFAULT_AGP_SEQ_COMPONENT_TYPE),
                    sdict->s[coords[cid * 3]].name, coords[cid * 3 + 1], coords[cid * 3 + 1] + coords[cid * 3 + 2] - 1,
                    sign? agp_orientation_val(AGP_OT_PLUS) : agp_orientation_val(AGP_OT_MINUS));
            slen += coords[cid * 3 + 2];
            if (i != (uint32) segs.n - 1) {
                fprintf(fo, "scaffold_%d\t%lld\t%lld\t%d\t%s\t%d\t%s\t%s\t%s\n", sid, slen + 1,
                        slen + DEFAULT_AGP_GAP_SIZE, ++fid, agp_component_type_val(DEFAULT_AGP_GAP_COMPONENT_TYPE),
                        DEFAULT_AGP_GAP_SIZE, agp_gap_type_val(DEFAULT_AGP_GAP_TYPE), agp_linkage_val(AGP_LG_YES),
                        agp_linkage_evidence_val(DEFAULT_AGP_LINKAGE_EVIDENCE));
                slen += DEFAULT_AGP_GAP_SIZE;
            }
        }
    }
    iostream_close(fp);

    free(coords);
    kv_destroy(segs);
    asm_destroy(dict);

    return 0;
}

static void print_help_agplift(FILE *fp_help)
{
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Usage: hictools agplift [options] <review.assembly> <liftover.agp> <contigs.fa[.fai|.idx]>\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help, "    -o STR             output file prefix (required for scaffolds FASTA output) [stdout]\n");
    fprintf(fp_help, "      --seq-ctype STR  AGP output sequence component type [%s]\n", agp_component_type_val(DEFAULT_AGP_SEQ_COMPONENT_TYPE));
    fprintf(fp_help, "      --gap-ctype STR  AGP output gap component type [%s]\n", agp_component_type_val(DEFAULT_AGP_GAP_COMPONENT_TYPE));
    fprintf(fp_help, "      --gap-link  STR  AGP output gap linkage evidence [%s]\n", agp_linkage_evidence_val(DEFAULT_AGP_LINKAGE_EVIDENCE));
    fprintf(fp_help, "      --gap-size  INT  AGP output gap size between sequence component [%d]\n", DEFAULT_AGP_GAP_SIZE);
    fprintf(fp_help, "    -h, --help         print this help\n");
    fprintf(fp_help, "    -V, --version      show version number\n");
    fprintf(fp_help, "\n");
    fprintf(fp_help, "Example: hictools agplift -o output review.assembly liftover.agp contigs.fa.fai\n");
    fprintf(fp_help, "\n");
}

static ko_longopt_t agplift_long_options[] = {
    { "seq-ctype",  ko_required_argument, 301 },
    { "gap-ctype",  ko_required_argument, 302 },
    { "gap-link",   ko_required_argument, 303 },
    { "gap-size",   ko_required_argument, 304 },
    { "help",       ko_no_argument,       'h' },
    { "version",    ko_no_argument,       'V' },
    { 0, 0, 0 }
};

static int main_agplift(int argc, char *argv[])
{
    FILE *fo;
    char *fa, *fa1, *out, *out1;

    liftrlimit();
    ht_realtime0 = realtime();

    const char *opt_str = "o:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, ret, is_idx;
    FILE *fp_help = stderr;
    sdict_t *sdict;
    fa = fa1 = out = out1 = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, agplift_long_options)) >= 0) {
        if (c == 301) {
            DEFAULT_AGP_SEQ_COMPONENT_TYPE = agp_component_type_key(opt.arg);
            if (DEFAULT_AGP_SEQ_COMPONENT_TYPE == AGP_CT_N ||
                    DEFAULT_AGP_SEQ_COMPONENT_TYPE == AGP_CT_U)
                fprintf(stderr, "[W::%s] a GAP component identifier will be used for sequences: %s\n",
                        __func__, opt.arg);
        } else if (c == 302) {
            DEFAULT_AGP_GAP_COMPONENT_TYPE = agp_component_type_key(opt.arg);
            if (DEFAULT_AGP_GAP_COMPONENT_TYPE != AGP_CT_N &&
                    DEFAULT_AGP_GAP_COMPONENT_TYPE != AGP_CT_U)
                fprintf(stderr, "[W::%s] a SEQ component identifier will be used for gaps: %s\n",
                        __func__, opt.arg);
        } else if (c == 303) {
            DEFAULT_AGP_LINKAGE_EVIDENCE = agp_linkage_evidence_key(opt.arg);
        } else if (c == 304) {
            DEFAULT_AGP_GAP_SIZE = atoi(opt.arg);
        } else if (c == 'o') {
            out = opt.arg;
        } else if (c == 'h') {
            fp_help = stdout;
        } else if (c == 'V') {
            puts(HICTOOLS_VERSION);
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
        print_help_agplift(stdout);
        return 0;
    }

    if (argc - opt.ind < 3) {
        fprintf(stderr, "[E::%s] missing input: three positional options required\n", __func__);
        print_help_agplift(stderr);
        return 1;
    }

    if (out) {
        out1 = (char *) malloc(strlen(out) + 35);
        sprintf(out1, "%s.FINAL.agp", out);
    }

    fo = out1 == 0? stdout : fopen(out1, "w");
    if (fo == 0) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }

    fa = argv[opt.ind + 2];
    is_idx = strlen(fa) > 4 && !strcmp(fa + strlen(fa) - 4, ".idx");

    ret = 0;
    sdict = is_idx? make_sdict_from_index(fa, 0) : make_sdict_from_fa(fa, 0);
    ret = assembly_to_agp(argv[opt.ind], argv[opt.ind + 1], sdict, fo);
    fflush(fo);
    if (out != 0)
        fclose(fo);

    if (!ret && !is_idx && out1) {
        fa1 = (char *) malloc(strlen(out) + 35);
        sprintf(fa1, "%s.FINAL.fa", out);
        fprintf(stderr, "[M::%s] writing FASTA file for scaffolds\n", __func__);
        FILE *fo1;
        fo1 = fopen(fa1, "w");
        if (fo1 == NULL) {
            fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, fa1);
            exit(EXIT_FAILURE);
        }
        write_fasta_file_from_agp(fa, out1, fo1, 60, 0);
        fclose(fo1);
    }

    sd_destroy(sdict);

    if (out1)
        free(out1);
    if (fa1)
        free(fa1);

    fprintf(stderr, "[M::%s] Version: %s\n", __func__, HICTOOLS_VERSION);
    fprintf(stderr, "[M::%s] CMD: hictools", __func__);
    int i;
    for (i = 0; i < argc; ++i)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__,
            realtime() - ht_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);

    return ret;
}

/***********************************************************************************
 * top-level dispatcher                                                             *
 ***********************************************************************************/

static int usage(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "Program: hictools\n");
    fprintf(fp, "Version: %s\n", HICTOOLS_VERSION);
    fprintf(fp, "Usage:   hictools <command> [options]\n");
    fprintf(fp, "Commands:\n");
    fprintf(fp, "    convert    convert HiC alignment file to binary BIN format\n");
    fprintf(fp, "    hiclink    generate HiC links from alignment file\n");
    fprintf(fp, "    prepare    prepare HiC data for downstream analysis\n");
    fprintf(fp, "    agplift    lift a Juicebox-curated assembly back to AGP/FASTA\n");
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
        puts(HICTOOLS_VERSION);
        return 0;
    }
    if (strcmp(argv[1], "convert") == 0)
        return main_convert(argc - 1, argv + 1);
    if (strcmp(argv[1], "hiclink") == 0)
        return main_hiclink(argc - 1, argv + 1);
    if (strcmp(argv[1], "prepare") == 0)
        return main_prepare(argc - 1, argv + 1);
    if (strcmp(argv[1], "agplift") == 0)
        return main_agplift(argc - 1, argv + 1);
    fprintf(stderr, "[E::%s] unrecognised command '%s'. Abort!\n", __func__, argv[1]);
    return 1;
}
