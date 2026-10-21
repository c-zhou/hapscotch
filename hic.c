/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2021 Chenxi Zhou <chnx.zhou@gmail.com>                          *
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
 * 23/06/21 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/resource.h>

#include "khash.h"
#include "kvec.h"
#include "kstring.h"
#include "bamlite.h"

#include "ONElib.h"

#include "sdict.h"
#include "cov.h"
#include "hic.h"
#include "misc.h"

KHASH_SET_INIT_STR(str)

const char *fileTypeNames[] = {
    [NOSET] = "NOSET",
    [BED] = "BED",
    [BAM] = "BAM",
    [BIN] = "BIN",
    [PA5] = "PA5",
    [ONE] = "ONE"
};

static inline void clamp(uint32 *v, uint32 l, uint32 u)
{
    if (*v < l) *v = l;
    if (*v > u) *v = u;
}

static inline char *parse_bam_rec(bam1_t *b, bam_header_t *h, uint32 *s, uint32 *e, uint8 *q, kstring_t *rname)
{
    // 0x4 0x100 0x200 0x400 0x800
    if ((b->core.flag & 0xF04) || (!b->core.flag))
        return 0;
    rname->l = 0;
    kputs(bam1_qname(b), rname);
    *s = b->core.pos;
    *e = get_target_end(b);
    *q = b->core.qual;

    return h->target_name[b->core.tid];
}

static void dump_links_from_bam_file(const char *f, sdict_t *dict, const char *out)
{
    bamFile fp;
    FILE *fo;
    bam_header_t *h;
    bam1_t *b;
    uint8 q, q0, q1;
    int8 buff;
    uint64 rec_c, pair_c, inter_c, intra_c, sd_l;
    enum bam_sort_order so;
    kstring_t _rname0 = {0, 0, 0}, *rname0 = &_rname0;
    kstring_t _rname1 = {0, 0, 0}, *rname1 = &_rname1;
    kstring_t *tmp;
    char *cname0, *cname1;
    uint32 data[6], *i0 = data, *s0 = i0+1, *e0 = s0+1, *i1 = e0+1, *s1 = i1+1, *e1 = s1+1;
    uint8 *bbuff, *dbuff, *ebuff;
    
    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);
    
    fp = bam_open(f, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    
    h = bam_header_read(fp);
    so = bam_hrecs_sort_order(h);
    if (so != ORDER_NAME) {
        fprintf(stderr, "[E::%s] BAM file %s is not sorted by read name\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    
    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    write_alnh_bin_header(fo);
    sd_l = write_binary_sequence_dictionary(fo, dict);

    MYMALLOC(bbuff, 1<<20);
    if (bbuff == NULL) {
        fprintf(stderr, "[E::%s] cannot allocate buffer array\n", __func__);
        exit(EXIT_FAILURE);
    }
    dbuff = bbuff;
    ebuff = bbuff + (1<<20);

    ks_resize(rname0, 256);
    ks_resize(rname1, 256);
    cname0 = cname1 = 0;
    q0 = q1 = 255;
    rec_c = pair_c = inter_c = intra_c = 0;
    buff = 0;
    b = bam_init1();
    
    fwrite(&pair_c, sizeof(uint64), 1, fo);

    while (bam_read1(fp, b) >= 0 ) {
        if (rec_c % REPORT_EPOC == 0 && rec_c)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs\n", __func__, rec_c / 1000000, pair_c);
        rec_c++;

        if (buff == 0) {
            q0 = 255;
            cname0 = parse_bam_rec(b, h, s0, e0, &q0, rname0);
            if (!cname0)
                continue;
            buff = 1;
        } else {
            q1 = 255;
            cname1 = parse_bam_rec(b, h, s1, e1, &q1, rname1);
            if (!cname1)
                continue;
            if (is_read_pair(rname0->s, rname1->s)) {
                *i0 = sd_get(dict, cname0);
                *i1 = sd_get(dict, cname1);
                q = MIN(q0, q1);

                if (*i0 == UINT32_MAX) {
                    k = kh_put(str, hmseq, cname0, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname0);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname0);
                    }
                } else if (*i1 == UINT32_MAX) {
                    k = kh_put(str, hmseq, cname1, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname1);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname1);
                    }
                } else {
                    clamp(s0, 0, dict->s[*i0].len);
                    clamp(e0, 0, dict->s[*i0].len);
                    clamp(s1, 0, dict->s[*i1].len);
                    clamp(e1, 0, dict->s[*i1].len);
                    
                    if (dbuff + ALNHBIN_ENTRY_SIZE < ebuff) {
                        fwrite(bbuff, 1, dbuff - bbuff, fo);
                        dbuff = bbuff;
                    }

                    if (*i0 <= *i1) {
                        memcpy(dbuff, data, sizeof(uint32)*6);
                    } else {
                        memcpy(dbuff, data+3, sizeof(uint32)*3);
                        memcpy(dbuff+sizeof(uint32)*3, data, sizeof(uint32)*3);
                    }
                    dbuff += ALNHBIN_ENTRY_SIZE;
                    dbuff[-1] = q;

                    if (*i0 == *i1)
                        ++intra_c;
                    else
                        ++inter_c;

                    ++pair_c;
                }

                buff = 0;
            } else {
                cname0 = cname1;
                *s0 = *s1;
                *e0 = *e1;
                q0 = q1;
                tmp = rname0;
                rname0 = rname1;
                rname1 = tmp;

                buff = 1;
            }
        }
    }

    if (dbuff > bbuff)
        fwrite(bbuff, 1, dbuff - bbuff, fo);
    free(bbuff);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    free(rname0->s);
    free(rname1->s);

    bam_destroy1(b);
    bam_header_destroy(h);
    bam_close(fp);
    // write pair number
    fseek(fo, sizeof(uint64) + sd_l, SEEK_SET);
    fwrite(&pair_c, sizeof(uint64), 1, fo);
    fseek(fo, pair_c * ALNHBIN_ENTRY_SIZE, SEEK_CUR);

    fprintf(stderr, "[M::%s] dumped %llu read pairs from %llu records: %llu intra links + %llu inter links \n", __func__, pair_c, rec_c, intra_c, inter_c);

    fclose(fo);
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

static void dump_links_from_bed_file(const char *f, sdict_t *dict, const char *out)
{
    FILE *fo;
    iostream_t *fp;
    char *line, *fields[6];
    kstring_t _cname = {0, 0, 0}, *cname = &_cname;
    kstring_t _rname = {0, 0, 0}, *rname = &_rname;
    uint8 q0, q1;
    int8 buff;
    uint64 rec_c, pair_c, inter_c, intra_c, sd_l;
    uint32 data[6], *i0 = data, *s0 = i0+1, *e0 = s0+1, *i1 = e0+1, *s1 = i1+1, *e1 = s1+1;
    uint8 *bbuff, *dbuff, *ebuff;

    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    fp = iostream_open(f);
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    write_alnh_bin_header(fo);
    sd_l = write_binary_sequence_dictionary(fo, dict);

    MYMALLOC(bbuff, 1<<20);
    if (bbuff == NULL) {
        fprintf(stderr, "[E::%s] cannot allocate buffer array\n", __func__);
        exit(EXIT_FAILURE);
    }
    dbuff = bbuff;
    ebuff = bbuff + (1<<20);

    q0 = q1 = 255;
    rec_c = pair_c = inter_c = intra_c = 0;
    buff = 0;
    
    fwrite(&pair_c, sizeof(uint64), 1, fo);

    while ((line = iostream_getline(fp)) != NULL) {
        if (is_empty_line(line) || parse_line(line, fields, 5) < 5)
            continue;
        
        if (rec_c % REPORT_EPOC == 0 && rec_c)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs\n", __func__, rec_c / 1000000, pair_c);
        rec_c++;

        if (buff == 0) {
            cname->l = 0, kputs(fields[0], cname);
            *s0 = atoi(fields[1]);
            *e0 = atoi(fields[2]);
            rname->l = 0, kputs(fields[3], rname);
            q0 = atoi(fields[4]);
            buff = 1;
        } else {
            if (is_read_pair(rname->s, fields[3])) {
                *i0 = sd_get(dict, cname->s);
                *i1 = sd_get(dict, fields[0]);
                
                if (*i0 == UINT32_MAX) {
                    k = kh_put(str, hmseq, cname->s, &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(cname->s);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, cname->s);
                    }
                } else if (*i1 == UINT32_MAX) {
                    k = kh_put(str, hmseq, fields[0], &absent);
                    if (absent) {
                        kh_key(hmseq, k) = strdup(fields[0]);
                        fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, fields[0]);
                    }
                } else {
                    *s1 = atoi(fields[1]);
                    *e1 = atoi(fields[2]);
                    q1 = atoi(fields[4]);

                    clamp(s0, 0, dict->s[*i0].len);
                    clamp(e0, 0, dict->s[*i0].len);
                    clamp(s1, 0, dict->s[*i1].len);
                    clamp(e1, 0, dict->s[*i1].len);

                    if (dbuff + ALNHBIN_ENTRY_SIZE < ebuff) {
                        fwrite(bbuff, 1, dbuff - bbuff, fo);
                        dbuff = bbuff;
                    }

                    if (*i0 <= *i1) {
                        memcpy(dbuff, data, sizeof(uint32)*6);
                    } else {
                        memcpy(dbuff, data+3, sizeof(uint32)*3);
                        memcpy(dbuff+sizeof(uint32)*3, data, sizeof(uint32)*3);
                    }
                    dbuff += ALNHBIN_ENTRY_SIZE;
                    dbuff[-1] = MIN(q0, q1);

                    if (*i0 == *i1)
                        ++intra_c;
                    else
                        ++inter_c;

                    ++pair_c;
                }
                buff = 0;
            } else {
                cname->l = 0, kputs(fields[0], cname);
                *s0 = atoi(fields[1]);
                *e0 = atoi(fields[2]);
                rname->l = 0, kputs(fields[3], rname);
                q0 = atoi(fields[4]);
                buff = 1;
            }
        }
    }

    if (dbuff > bbuff)
        fwrite(bbuff, 1, dbuff - bbuff, fo);
    free(bbuff);

    free(cname->s);
    free(rname->s);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    iostream_close(fp);
    // write pair number
    fseek(fo, sizeof(uint64) + sd_l, SEEK_SET);
    fwrite(&pair_c, sizeof(uint64), 1, fo);
    fseek(fo, pair_c * ALNHBIN_ENTRY_SIZE, SEEK_CUR);

    fprintf(stderr, "[M::%s] dumped %llu read pairs from %llu records: %llu intra links + %llu inter links \n", __func__, pair_c, rec_c, intra_c, inter_c);

    fclose(fo);
}

static void dump_links_from_pa5_file(const char *f, sdict_t *dict, int read_len, const char *out)
{
    FILE *fo;
    iostream_t *fp;
    char *line, *fields[7];
    uint8 q, q0, q1;
    uint64 rec_c, pair_c, inter_c, intra_c, sd_l;
    uint32 data[6], *i0 = data, *s0 = i0+1, *e0 = s0+1, *i1 = e0+1, *s1 = i1+1, *e1 = s1+1;
    uint8 *bbuff, *dbuff, *ebuff;

    khash_t(str) *hmseq; // for absent sequences
    khint_t k;
    int absent;
    hmseq = kh_init(str);

    fp = iostream_open(f);
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }

    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    write_alnh_bin_header(fo);
    sd_l = write_binary_sequence_dictionary(fo, dict);

    MYMALLOC(bbuff, 1<<20);
    if (bbuff == NULL) {
        fprintf(stderr, "[E::%s] cannot allocate buffer array\n", __func__);
        exit(EXIT_FAILURE);
    }
    dbuff = bbuff;
    ebuff = bbuff + (1<<20);

    read_len >>= 1;
    rec_c = pair_c = inter_c = intra_c = 0;
    
    fwrite(&pair_c, sizeof(uint64), 1, fo);

    while ((line = iostream_getline(fp)) != NULL) {
        if (is_empty_line(line) || parse_line(line, fields, 7) < 7)
            continue;

        if (rec_c % REPORT_EPOC == 0 && rec_c)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs\n", __func__, rec_c / 1000000, pair_c);
        rec_c++;

        *i0 = sd_get(dict, fields[1]);
        *i1 = sd_get(dict, fields[3]);

        if (*i0 == UINT32_MAX) {
            k = kh_put(str, hmseq, fields[1], &absent);
            if (absent) {
                kh_key(hmseq, k) = strdup(fields[1]);
                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, fields[1]);
            }
        } else if (*i1 == UINT32_MAX) {
            k = kh_put(str, hmseq, fields[3], &absent);
            if (absent) {
                kh_key(hmseq, k) = strdup(fields[3]);
                fprintf(stderr, "[W::%s] sequence \"%s\" not found \n", __func__, fields[3]);
            }
        } else {
            *s0 = atoi(fields[2]);
            *s1 = atoi(fields[4]);
            q0  = atoi(fields[5]);
            q1  = atoi(fields[6]);
            q = MIN(q0, q1);

            *e0 = *s0 + read_len;
            *e1 = *s1 + read_len;

            clamp(s0, 0, dict->s[*i0].len);
            clamp(e0, 0, dict->s[*i0].len);
            clamp(s1, 0, dict->s[*i1].len);
            clamp(e1, 0, dict->s[*i1].len);

            if (dbuff + ALNHBIN_ENTRY_SIZE < ebuff) {
                fwrite(bbuff, 1, dbuff - bbuff, fo);
                dbuff = bbuff;
            }

            if (*i0 <= *i1) {
                memcpy(dbuff, data, sizeof(uint32)*6);
            } else {
                memcpy(dbuff, data+3, sizeof(uint32)*3);
                memcpy(dbuff+sizeof(uint32)*3, data, sizeof(uint32)*3);
            }
            dbuff += ALNHBIN_ENTRY_SIZE;
            dbuff[-1] = q;

            if (*i0 == *i1)
                ++intra_c;
            else
                ++inter_c;

            ++pair_c;
        }
    }

    if (dbuff > bbuff)
        fwrite(bbuff, 1, dbuff - bbuff, fo);
    free(bbuff);

    for (k = 0; k < kh_end(hmseq); ++k)
        if (kh_exist(hmseq, k))
            free((char *) kh_key(hmseq, k));
    kh_destroy(str, hmseq);

    iostream_close(fp);
    // write pair number
    fseek(fo, sizeof(uint64) + sd_l, SEEK_SET);
    fwrite(&pair_c, sizeof(uint64), 1, fo);
    fseek(fo, pair_c * ALNHBIN_ENTRY_SIZE, SEEK_CUR);

    fprintf(stderr, "[M::%s] dumped %llu read pairs from %llu records: %llu intra links + %llu inter links \n", __func__, pair_c, rec_c, intra_c, inter_c);

    fclose(fo);
}

static char *oneSchemaText =
    "1 3 def 1 0               schema for ONEfile\n"
    ".\n"
    "P 3 map                   MAP\n"
    "O S 2 3 INT 3 INT         query sequence: index in source file (1-based) length\n"
    "D I 1 6 STRING            identifier from source file (if requested)\n"
    "D F 1 4 CHAR              filter: Z zero-length, Q quality, G poly-G (Illumina bad read)\n"
    "D M 3 3 INT 3 INT 3 INT   mem: start, end (0-based), count\n"
    "D X 2 3 INT 3 DNA         missing syncmer not found in graph: start coordinate, sequence\n"
    "D U 3 3 INT 3 INT 3 INT   unique mapping: file, path, offset (negative offset = reverse)"
    ;

static void dump_links_from_one_file(const char *f, sdict_t *dict, const char *out)
{
    FILE *fo;
    OneFile *fp;
    OneSchema *schema;
    uint8 q;
    int8 buff;
    char lineType;
    int64 r0, r1, rl;
    uint64 rec_c, pair_c, inter_c, intra_c, sd_l;
    uint32 data[6], *i0 = data, *s0 = i0+1, *e0 = s0+1, *i1 = e0+1, *s1 = i1+1, *e1 = s1+1;
    uint8 *bbuff, *dbuff, *ebuff;

    schema = oneSchemaCreateFromText(oneSchemaText);
    if (schema == NULL) {
        fprintf(stderr, "[E::%s] cannot create schema from text\n", __func__);
        exit(EXIT_FAILURE);
    }
    fp = oneFileOpenRead(f, schema, "map", 1);
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, f);
        exit(EXIT_FAILURE);
    }
    if (fp->info['S'] == NULL) {
        fprintf(stderr, "[E::%s] no sequence alignments found in map file\n", __func__);
        exit(EXIT_FAILURE);
    }
    oneSchemaDestroy(schema);

    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    write_alnh_bin_header(fo);
    sd_l = write_binary_sequence_dictionary(fo, dict);

    MYMALLOC(bbuff, 1<<20);
    if (bbuff == NULL) {
        fprintf(stderr, "[E::%s] cannot allocate buffer array\n", __func__);
        exit(EXIT_FAILURE);
    }
    dbuff = bbuff;
    ebuff = bbuff + (1<<20);

    q = 255;
    rec_c = pair_c = inter_c = intra_c = 0;
    buff = 0;
    
    fwrite(&pair_c, sizeof(uint64), 1, fo);

    while ((lineType = oneReadLine(fp)) && lineType != 'S');
    while (lineType) {
        if (rec_c % REPORT_EPOC == 0 && rec_c)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs\n", __func__, rec_c / 1000000, pair_c);
        rec_c++;

        if (buff == 0) {
            r0 = oneInt(fp, 0);
            *i0 = UINT32_MAX;
            rl = 0;
            while ((lineType = oneReadLine(fp))) {
                if (lineType == 'S')
                    break;
                else if (lineType == 'M') {
                    rl = oneInt(fp, 1) - oneInt(fp, 0);
                } else if (lineType == 'U') {
                    *i0 = oneInt(fp, 1) - 1;
                    *s0 = llabs(oneInt(fp, 2));
                }
            }
            *e0 = *s0 + rl;
            buff = 1;
        } else {
            r1 = oneInt(fp, 0);
            *i1 = UINT32_MAX;
            rl = 0;
            while ((lineType = oneReadLine(fp))) {
                if (lineType == 'S')
                    break;
                else if (lineType == 'M') {
                    rl = oneInt(fp, 1) - oneInt(fp, 0);
                } else if (lineType == 'U') {
                    *i1 = oneInt(fp, 1) - 1;
                    *s1 = llabs(oneInt(fp, 2));
                }
            }
            *e1 = *s1 + rl;

            if ((r0 & 1) && (r0 + 1 == r1)) {
                buff = 0;

                if (*i0 != UINT32_MAX && *i1 != UINT32_MAX) {
                    clamp(s0, 0, dict->s[*i0].len);
                    clamp(e0, 0, dict->s[*i0].len);
                    clamp(s1, 0, dict->s[*i1].len);
                    clamp(e1, 0, dict->s[*i1].len);

                    if (dbuff + ALNHBIN_ENTRY_SIZE < ebuff) {
                        fwrite(bbuff, 1, dbuff - bbuff, fo);
                        dbuff = bbuff;
                    }

                    if (*i0 <= *i1) {
                        memcpy(dbuff, data, sizeof(uint32)*6);
                    } else {
                        memcpy(dbuff, data+3, sizeof(uint32)*3);
                        memcpy(dbuff+sizeof(uint32)*3, data, sizeof(uint32)*3);
                    }
                    dbuff += ALNHBIN_ENTRY_SIZE;
                    dbuff[-1] = q;

                    if (*i0 == *i1)
                        ++intra_c;
                    else
                        ++inter_c;

                    ++pair_c;
                }
            } else {
                buff = 1;

                *s0 = *s1;
                *e0 = *e1;
                *i0 = *i1;
            }
        }
    }

    if (dbuff > bbuff)
        fwrite(bbuff, 1, dbuff - bbuff, fo);
    free(bbuff);

    oneFileClose(fp);

    // write pair number
    fseek(fo, sizeof(uint64) + sd_l, SEEK_SET);
    fwrite(&pair_c, sizeof(uint64), 1, fo);
    fseek(fo, pair_c * ALNHBIN_ENTRY_SIZE, SEEK_CUR);

    fprintf(stderr, "[M::%s] dumped %llu read pairs from %llu records: %llu intra links + %llu inter links \n", __func__, pair_c, rec_c, intra_c, inter_c);

    fclose(fo);
}

char *write_binary_hic_data(char *hic_file, fileType_t f_type, sdict_t *dicts, int read_len, char *pref_out)
{
    if (hic_file == NULL)
        return NULL;
    
    if (f_type == NOSET) {
        char *ext1, *ext2;
        ext1 = strlen(hic_file) >= 4? (hic_file + strlen(hic_file) - 4) : NULL;
        ext2 = strlen(hic_file) >= 7? (hic_file + strlen(hic_file) - 7) : NULL;
        if (ext1 && !strcasecmp(ext1, ".bam")) f_type = BAM;
        else if (ext1 && !strcasecmp(ext1, ".bin")) f_type = BIN;
        else if (ext1 && !strcasecmp(ext1, "1map")) f_type = ONE;
        else if ((ext1 && !strcasecmp(ext1, ".bed")) || (ext2 && !strcasecmp(ext2, ".bed.gz"))) f_type = BED;
        else if ((ext1 && !strcasecmp(ext1, ".pa5")) || (ext2 && !strcasecmp(ext2, ".pa5.gz"))) f_type = PA5;
        else {
            fprintf(stderr, "[E::%s] unknown link file format: %s\n", __func__, hic_file);
            fprintf(stderr, "[E::%s] File extension .bam, .bed, .pa5, .1map or .bin or --file-type is expected\n", __func__);
            exit(EXIT_FAILURE);
        }
    }

    if (f_type == BIN) {
        if (strcmp(hic_file, "-") == 0 || *hic_file == '<') {
            fprintf(stderr, "[E::%s] BIN file format from STDIN is not supported\n", __func__);
            exit(EXIT_FAILURE);
        }

        return strdup(hic_file);
    }
    
    char *hic_bfile;
    MYMALLOC(hic_bfile, strlen(pref_out) + 5);
    sprintf(hic_bfile, "%s.bin", pref_out);
    
    fprintf(stderr, "[M::%s] dump hic links (%s) to binary file: %s\n", __func__, fileTypeNames[f_type], hic_bfile);
    if (f_type == BAM)
        dump_links_from_bam_file(hic_file, dicts, hic_bfile);
    else if (f_type == BED)
        dump_links_from_bed_file(hic_file, dicts, hic_bfile);
    else if (f_type == ONE)
        dump_links_from_one_file(hic_file, dicts, hic_bfile);
    else if (f_type == PA5)
        dump_links_from_pa5_file(hic_file, dicts, read_len, hic_bfile);

    return hic_bfile;
}

/* converts each input file independently, then concatenates their raw HiC records into one ALNHBIN file */
char *write_binary_hic_data_multi(char **files, int n_files, fileType_t f_type, sdict_t *dicts, int read_len, char *out)
{
    char *merged, *part, tmp_pref[4096];
    FILE *fo, *fp;
    int64 magic_number;
    uint64 sd_l, pair_n, total_n, m;
    uint8 *buffer;
    int i;

    if (n_files <= 0)
        return NULL;
    if (n_files == 1)
        return write_binary_hic_data(files[0], f_type, dicts, read_len, out);

    MYMALLOC(merged, strlen(out) + 5);
    sprintf(merged, "%s.bin", out);
    fo = fopen(merged, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, merged);
        exit(EXIT_FAILURE);
    }
    write_alnh_bin_header(fo);
    sd_l = write_binary_sequence_dictionary(fo, dicts);
    total_n = 0;
    fwrite(&total_n, sizeof(uint64), 1, fo);

    MYMALLOC(buffer, ALNHBIN_ENTRY_SIZE * FREAD_BUFF_SIZE);
    for (i = 0; i < n_files; ++i) {
        snprintf(tmp_pref, sizeof(tmp_pref), "%s.part%d", out, i);
        part = write_binary_hic_data(files[i], f_type, dicts, read_len, tmp_pref);

        fp = fopen(part, "r");
        if (fp == NULL) {
            fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, part);
            exit(EXIT_FAILURE);
        }
        m = fread(&magic_number, sizeof(int64), 1, fp);
        if (m != 1) bin_fread_error();
        if (!is_valid_alnh_bin_header(magic_number)) {
            fprintf(stderr, "[E::%s] %s is not a valid raw HiC BIN file\n", __func__, part);
            exit(EXIT_FAILURE);
        }
        binary_fseek_skip_sdict(fp);
        m = fread(&pair_n, sizeof(uint64), 1, fp);
        if (m != 1) bin_fread_error();
        total_n += pair_n;

        while ((m = fread(buffer, 1, ALNHBIN_ENTRY_SIZE * FREAD_BUFF_SIZE, fp)) > 0)
            fwrite(buffer, 1, m, fo);
        fclose(fp);

        if (strcmp(part, files[i]) != 0)
            remove(part); // only delete the temp file we created, not a pre-existing .bin input
        free(part);
    }
    free(buffer);

    fseek(fo, sizeof(int64) + sd_l, SEEK_SET);
    fwrite(&total_n, sizeof(uint64), 1, fo);
    fclose(fo);

    fprintf(stderr, "[M::%s] merged %d HiC input files into %s (%llu read pairs)\n", __func__, n_files, merged, total_n);

    return merged;
}

static int hic_cmpfunc(const void *a, const void *b)
{
    hic_t *x = (hic_t *) a, *y = (hic_t *) b;
    if (x->aseq != y->aseq)
        return (x->aseq > y->aseq) - (x->aseq < y->aseq);
    if (x->bseq != y->bseq)
        return (x->bseq > y->bseq) - (x->bseq < y->bseq);
    if (x->apos != y->apos)
        return (x->apos > y->apos) - (x->apos < y->apos);
    return (x->bpos > y->bpos) - (x->bpos < y->bpos);
}

static uint64 hic_sort_merge(hic_t *hics, uint64 n)
{
    if (n == 0)
        return 0;

    qsort(hics, n, sizeof(hic_t), hic_cmpfunc);
    
    hic_t *s, *c, *e;
    s = hics;
    c = hics;
    e = hics + n;
    
    while (++c < e) {
        if (hic_cmpfunc(s, c) == 0)
            s->nhic += c->nhic;
        else
            *++s = *c;
    }

    return (s - hics) + 1;
}

static inline int pseudo_position(sdict_t *dicts, int32 **smap, 
    uint32 i, uint32 s, uint32 e, uint32 *_i, uint32 *_s, uint32 *_e)
{
    // map sequences to backbone position
    // each sequence will have a segment block, where
    // seg[0]       = sequence length
    // seg[1]       = sequence orientation
    // seg[2]       = number of segment pairs
    // seg[2*i+3]   = end position on sequence
    // seg[2*i+4]   = backbone sequence index
    // to map a position {p} on sequence {i} to backbone position
    // 1. find the first segment pair {s} where seg[2*s+3] > p
    // 2. seg[2*s+4] is the backbone sequence index
    // 3. seg[2*s+3] - p is the distance to the end of the backbone sequence
    int a, l, m, r, n;

    int32 *map = smap[i];
    l = map[0];
    r = map[1];
    n = map[2];
    
    if (!n) return 1; // unsuccessful mapping
    if (r) { // reverse orientation
        a = s;
        s = l - e;
        e = l - a;
    }

    map += 3;
    for (a = 0; a < n; a++, map += 2) {
        if (*map > s) {
            if (*map < e)
                return 1; // crossing the segment boundary
            *_i = map[1];
            m = dicts->s[*_i].len;
            if (m < *map - s)
                return 1; // out of boundary
            *_s = m - (*map - s);
            *_e = m - (*map - e);
            return 0;
        }
    }

    return 1;
}

// convert a raw-sequence interval to corrected piece (break AGP object) coordinates
// returns non-zero if the interval crosses a piece boundary or is out of range
static inline int piece_position(asm_dict_t *bd, uint32 *i, uint32 *s, uint32 *e)
{
    uint32 q0, q1;
    uint64 p0, p1;

    if (sd_coordinate_conversion(bd, *i, *s, &q0, &p0, 1) != CC_SUCCESS)
        return 1;
    if (*e > *s) {
        if (sd_coordinate_conversion(bd, *i, *e - 1, &q1, &p1, 1) != CC_SUCCESS ||
            q1 != q0)
            return 1; // crossing a piece boundary
        *e = (uint32) (p1 > p0? p1 : p0) + 1;
        *s = (uint32) (p1 > p0? p0 : p1);
    } else {
        *e = (uint32) p0;
        *s = (uint32) p0;
    }
    *i = q0;

    return 0;
}

void write_binary_hic_data_pseudo_yahs(char *bf, sdict_t *dicts, asm_dict_t *adict, int32 **smap, char *out)
{
    uint64 m, n, sd_l, pair_n, pair_c, max_m, n_recs;
    int64 magic_number;
    cov_t *cov;
    FILE *fp, *fo;
    uint8 *bbfp, *ebfp, *bfp;
    uint8 *bbfo, *ebfo, *bfo;
    uint8 q;
    uint32 data[6], *i0 = data, *s0 = i0+1, *i1 = s0+1, *s1 = i1+1, *e0 = s1+1, *e1 = e0+1;
    double q_drop = 0.1;

    fp = fopen(bf, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, bf);
        exit(EXIT_FAILURE);
    }

    m = fread(&magic_number, sizeof(int64), 1, fp);
    if (m != 1) bin_fread_error();
    if (!is_valid_alnh_bin_header(magic_number)) {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        exit(EXIT_FAILURE);
    }
    binary_fseek_skip_sdict(fp);
    m = fread(&pair_n, sizeof(uint64), 1, fp);
    if (m != 1) bin_fread_error();
    
    fo = fopen(out, "w");
    if (fo == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for writing\n", __func__, out);
        exit(EXIT_FAILURE);
    }
    write_yahs_bin_header(fo);
    sd_l = write_binary_sequence_dictionary(fo, dicts);

    MYMALLOC(bfp, ALNHBIN_ENTRY_SIZE * FREAD_BUFF_SIZE);
    if (bfp == NULL)
        mem_alloc_error("input buffer array");

    MYMALLOC(bfo, 1<<20);
    if (bfo == NULL)
        mem_alloc_error("output buffer array");
    bbfo = bfo;
    ebfo = bfo + (1<<20);

    cov = cov_init(dicts->n);
    n = 0;
    max_m = 0x7FFFFFFULL; // 128MB - ~1GB mem for covs

    pair_c = 0;
    n_recs = 0;
    n = 0;
    q = 255;

    fwrite(&n_recs, sizeof(uint64), 1, fo);

    while (pair_c < pair_n) {
        m = fread(bfp, sizeof(uint8) * ALNHBIN_ENTRY_SIZE, FREAD_BUFF_SIZE, fp);
        for (bbfp = bfp, ebfp = bfp + m * ALNHBIN_ENTRY_SIZE; 
            bbfp < ebfp && pair_c++ < pair_n; bbfp += ALNHBIN_ENTRY_SIZE) {

            if (pair_c % REPORT_EPOC == 0)
                fprintf(stderr, "[M::%s] %llu million records processed\n", __func__, pair_c / 1000000);
            
            memcpy(i0, bbfp,  sizeof(uint32));
            memcpy(s0, bbfp + sizeof(uint32)*1, sizeof(uint32));
            memcpy(e0, bbfp + sizeof(uint32)*2, sizeof(uint32));
            memcpy(i1, bbfp + sizeof(uint32)*3, sizeof(uint32));
            memcpy(s1, bbfp + sizeof(uint32)*4, sizeof(uint32));
            memcpy(e1, bbfp + sizeof(uint32)*5, sizeof(uint32));

            // with a break AGP the records are on the raw sequences
            // convert them to piece coordinates first
            if (adict &&
                (piece_position(adict, i0, s0, e0) ||
                 piece_position(adict, i1, s1, e1)))
                continue;

            // map to pseudo positions
            if (pseudo_position(dicts, smap, *i0, *s0, *e0, i0, s0, e0) ||
                pseudo_position(dicts, smap, *i1, *s1, *e1, i1, s1, e1))
                continue;

            // add to coverage array
            kv_push(uint64, cov->p[*i0], (uint64) (*s0) << 32 | (uint32)  1);
            kv_push(uint64, cov->p[*i0], (uint64) (*e0) << 32 | (uint32) -1);
            kv_push(uint64, cov->p[*i1], (uint64) (*s1) << 32 | (uint32)  1);
            kv_push(uint64, cov->p[*i1], (uint64) (*e1) << 32 | (uint32) -1);
            n += 4;
            if (n > max_m) {
                m = pos_compression(cov); // m is the position size after compression
                fprintf(stderr, "[M::%s] position compression n = %llu, m = %llu, max_m = %llu\n", __func__, n, m, max_m);
                if (m > n>>1) {
                    max_m <<= 1; // increase m limit if compression ratio smaller than 0.5
                    fprintf(stderr, "[M::%s] position memory buffer expanded max_m = %llu\n", __func__, max_m);
                }
                n = m;
            }

            // write hic record
            *s0 = (*s0>>1) + (*e0>>1) + ((*s0 & *e0) & 1);
            *s1 = (*s1>>1) + (*e1>>1) + ((*s1 & *e1) & 1);
            
            if (bbfo + YAHSBIN_ENTRY_SIZE < ebfo) {
                fwrite(bfo, 1, bbfo - bfo, fo);
                bbfo = bfo;
            }

            if (*i0 <= *i1) {
                memcpy(bbfo, data, sizeof(uint32)*4);
            } else {
                memcpy(bbfo, data+2, sizeof(uint32)*2);
                memcpy(bbfo+sizeof(uint32)*2, data, sizeof(uint32)*2);
            }
            bbfo += YAHSBIN_ENTRY_SIZE;
            bbfo[-1] = q;

            n_recs++;
        }
    }

    if (bbfo > bfo)
        fwrite(bfo, 1, bbfo - bfo, fo);
    
    free(bfo);
    free(bfp);
    fclose(fp);
    
    m = pos_compression(cov);
    fprintf(stderr, "[M::%s] position compression n = %llu, m = %llu, max_m = %llu\n", __func__, n, m, max_m);

    fprintf(stderr, "[M::%s] processed %llu read pairs\n", __func__, pair_c);
    fprintf(stderr, "[M::%s] retained %llu records\n", __func__, n_recs);
    
    // write pair number
    fseek(fo, sizeof(uint64) + sd_l, SEEK_SET);
    fwrite(&n_recs, sizeof(uint64), 1, fo);
    fseek(fo, n_recs * YAHSBIN_ENTRY_SIZE, SEEK_CUR);

    // cov_t *cov = bed_cstats(f, dict);
    cov_norm_t *cov_norm = calc_cov_norms(cov, dicts, q_drop);
    fwrite(&cov_norm->n, sizeof(uint64), 1, fo);
    fwrite(cov_norm->norm[0], sizeof(double), cov_norm->n, fo);
    cov_destroy(cov);
    cov_norm_destroy(cov_norm);

    fclose(fo);
}

int match_binary_file_sdict(char *f, sdict_t *dict)
{
    uint32_t i, m, n, s, *lens;
    uint64_t l;
    int64_t magic_number;
    char *names, *p;
    FILE *fp;

    fp = fopen(f, "r");
    if (fp == NULL)
        return 1;

    m = fread(&magic_number, sizeof(int64_t), 1, fp);
    if (m != 1) bin_fread_error();
    if (!is_valid_bin_header(magic_number))
        return 2;

    m = fread(&n, sizeof(uint32_t), 1, fp);
    if (m != 1) bin_fread_error();
    if (n != dict->n)
        return 3;
    
    lens = (uint32_t *) calloc(n, sizeof(uint32_t));
    m = fread(lens, sizeof(uint32_t), n, fp);
    if (m != n) bin_fread_error();
    for (i = 0; i < n; ++i) {
        if (dict->s[i].len != lens[i]) {
            free(lens);
            return 4;
        }
    }
    free(lens);
    
    m = fread(&l, sizeof(uint64_t), 1, fp);
    if (m != 1) bin_fread_error();
    names = (char *) calloc(l, 1);
    m = fread(names, 1, l, fp);
    if (m != l) bin_fread_error();
    for (p = names, i = 0; i < n; ++i) {
        s = strlen(dict->s[i].name);
        if (p - names >= l || strncmp(p, dict->s[i].name, s)) {
            free(names);
            return 5;
        }
        p += s + 1;
    }
    free(names);

    fclose(fp);

    return 0;
}


/*************************************************************************
 * External-memory HiC loader
 *
 * 1. scan:  raw records (canonical + mirrored copy for inter links) are
 *           streamed into on-disk buckets by contiguous aseq ranges.
 * 2. split: buckets larger than the per-leaf limit are re-split on the next
 *           sort key (aseq -> bseq -> apos -> bpos), giving an ordered list
 *           of leaf buckets that each fit in memory.
 * 3. merge: worker threads load leaves, pack each record into a 64-bit key,
 *           LSD radix sort it and collapse duplicates; the main thread writes
 *           the leaves to one sorted hic_t file strictly in order. At most
 *           `threads` leaves are held in memory at any time.
 * The file is mmap'ed read-only.
 *************************************************************************/

typedef struct { uint32 a, ap, b, bp; } hraw_t; // 16 bytes, implicit count 1

#define HRAW_RBUF (1<<16)        // records per read buffer (1 MB)
#define HOUT_WBUF (1<<18)        // hic_t records per output buffer (5 MB)
#define MAX_BUCKETS_HARD 4096    // upper bound on simultaneously open buckets
#define RADIX_BITS 11
#define RADIX_SIZE (1 << RADIX_BITS)

int64 HIC_MEM_LIMIT = 8LL << 30;
int   HIC_THREADS   = 1;
char *HIC_TMP_DIR   = ".";

static int hic_max_buckets = 256;   // set from RLIMIT_NOFILE
static int64 hic_wbuf_n = 4096;     // records per bucket write buffer
static int64 hic_n_split = 0, hic_b_split = 0;

typedef struct {
    char *path;
    FILE *fp;
    int64 n;
    int64 nb;
    hraw_t *buf;
} bucket_t;

typedef struct {
    FILE *fp;          // final sorted file; NULL for in-memory mode
    hic_t *a;          // output buffer (disk mode) or the whole array (memory mode)
    size_t n, m;
    int64 nout, ncanon;
    uint32 nseq, nwin; // key ranges for histograms
} hout_t;

// raise the open-file soft limit and size the bucket cap from it
static void hic_set_file_limit(void)
{
    struct rlimit rl;
    int64 lim = 256;
    if (getrlimit(RLIMIT_NOFILE, &rl) == 0) {
        rlim_t want = MAX_BUCKETS_HARD + 64;
        if (rl.rlim_cur != RLIM_INFINITY && rl.rlim_cur < want) {
            struct rlimit r2 = rl;
            r2.rlim_cur = (rl.rlim_max == RLIM_INFINITY || rl.rlim_max > want) ? want : rl.rlim_max;
            if (setrlimit(RLIMIT_NOFILE, &r2) == 0) rl = r2;
        }
        lim = rl.rlim_cur == RLIM_INFINITY ? MAX_BUCKETS_HARD + 64 : (int64) rl.rlim_cur;
    }
    hic_max_buckets = (int) MAX(16, MIN((int64) MAX_BUCKETS_HARD, lim - 64));
}

// bucket write buffers together use at most ~1/4 of the memory limit
static void hic_set_wbuf(int64 nbuckets)
{
    int64 w = HIC_MEM_LIMIT / 4 / (int64) sizeof(hraw_t) / MAX(nbuckets, 1);
    hic_wbuf_n = MAX(1024, MIN(4096, w));
}

static void hfile_read_error(const char *path)
{
    fprintf(stderr, "[E::%s] error reading temp file %s\n", __func__, path);
    exit(EXIT_FAILURE);
}

static void hfile_remove(const char *path)
{
    if (unlink(path) != 0)
        fprintf(stderr, "[W::%s] could not remove temp file %s: %s\n", __func__, path, strerror(errno));
}

static int tmp_counter = 0;
static char *hic_run_dir = NULL; // this run's temp directory, created on first use

static const char *hic_temp_dir(void)
{
    const char *dir = HIC_TMP_DIR;
    if (dir == NULL) dir = getenv("TMPDIR");
    if (dir == NULL || *dir == 0) dir = ".";
    return dir;
}

// All temp files of a run live in one fresh directory <temp>/hlk_XXXXXX (mkdtemp): it is
// reported in the log, names cannot collide between runs, and if a run is killed there
// is exactly one directory to delete. Files are removed as soon as they are consumed.
static char *hic_temp_path(const char *tag)
{
    char *p;
    if (hic_run_dir == NULL) {
        const char *dir = hic_temp_dir();
        hic_run_dir = (char *) malloc(strlen(dir) + 16);
        if (hic_run_dir == NULL) mem_alloc_error("temp path");
        sprintf(hic_run_dir, "%s/hlk_XXXXXX", dir);
        if (mkdtemp(hic_run_dir) == NULL) {
            fprintf(stderr, "[E::%s] cannot create a temp directory in %s: %s\n", __func__, dir, strerror(errno));
            exit(EXIT_FAILURE);
        }
        fprintf(stderr, "[M::%s] temporary files in %s\n", __func__, hic_run_dir);
    }
    p = (char *) malloc(strlen(hic_run_dir) + 32);
    if (p == NULL) mem_alloc_error("temp path");
    sprintf(p, "%s/%06d.%s", hic_run_dir, tmp_counter++, tag);
    return p;
}

// remove the run directory once all its files are gone
static void hic_run_dir_remove(void)
{
    if (hic_run_dir && rmdir(hic_run_dir) == 0) {
        free(hic_run_dir);
        hic_run_dir = NULL;
        tmp_counter = 0;
    }
}

static void bucket_open(bucket_t *bk)
{
    bk->path = hic_temp_path("raw");
    bk->fp = fopen(bk->path, "wb");
    if (bk->fp == NULL) {
        fprintf(stderr, "[E::%s] cannot create temp file %s: %s\n", __func__, bk->path, strerror(errno));
        exit(EXIT_FAILURE);
    }
    bk->n = 0;
    bk->nb = 0;
    MYMALLOC(bk->buf, hic_wbuf_n);
    if (bk->buf == NULL) mem_alloc_error("bucket buffer");
}

static void bucket_flush(bucket_t *bk)
{
    if (bk->nb && (int64) fwrite(bk->buf, sizeof(hraw_t), bk->nb, bk->fp) != bk->nb) {
        fprintf(stderr, "[E::%s] write error on %s (disk full?)\n", __func__, bk->path);
        exit(EXIT_FAILURE);
    }
    bk->nb = 0;
}

static inline void bucket_push(bucket_t *bk, const hraw_t *r)
{
    bk->buf[bk->nb++] = *r;
    bk->n++;
    if (bk->nb == hic_wbuf_n) bucket_flush(bk);
}

static void bucket_close(bucket_t *bk)
{
    bucket_flush(bk);
    if (fclose(bk->fp) != 0) {
        fprintf(stderr, "[E::%s] error closing %s\n", __func__, bk->path);
        exit(EXIT_FAILURE);
    }
    bk->fp = NULL;
    free(bk->buf);
    bk->buf = NULL;
}

static int hraw_cmpfunc(const void *x, const void *y)
{
    const hraw_t *p = x, *q = y;
    if (p->a != q->a) return (p->a > q->a) - (p->a < q->a);
    if (p->b != q->b) return (p->b > q->b) - (p->b < q->b);
    if (p->ap != q->ap) return (p->ap > q->ap) - (p->ap < q->ap);
    return (p->bp > q->bp) - (p->bp < q->bp);
}

static void hout_flush(hout_t *out)
{
    if (out->fp && out->n) {
        if (fwrite(out->a, sizeof(hic_t), out->n, out->fp) != out->n) {
            fprintf(stderr, "[E::%s] write error on sorted HiC file (disk full?)\n", __func__);
            exit(EXIT_FAILURE);
        }
        out->n = 0;
    }
}

static inline void hout_push(hout_t *out, const hraw_t *r, int64 cnt)
{
    if (out->n == out->m) {
        if (out->fp) {
            hout_flush(out);
        } else {
            out->m = out->m ? out->m << 1 : 1 << 20;
            MYREALLOC(out->a, out->m);
            if (out->a == NULL) mem_alloc_error("HiC links array");
        }
    }
    out->a[out->n++] = (hic_t) {r->a, r->ap, r->b, r->bp, (int) cnt};
    out->nout++;
    if (r->a <= r->b) out->ncanon++;
}

/* sorting and duplicate collapsing */
typedef struct {
    void *buf;          // record memory, owned
    void *cbuf;         // separately allocated count memory, if any
    uint32 *cnt;        // count per unique record
    int64 nu;           // number of unique records
    int packed;         // 1: key[] holds packed records; 0: raw[] (fallback)
    uint64 *key;
    hraw_t *raw;
    uint32 amin, bmin, apmin, bpmin;
    int wb, wap, wbp;   // bit widths of the b, ap, bp fields in the key
} hres_t;

static inline int bit_width32(uint32 x) { return x ? 32 - __builtin_clz(x) : 0; }

// sort r[0..n) and collapse identical records; r must be malloc'ed and is taken over by res.
// Each record is packed as (a-amin | b-bmin | ap-apmin | bp-bpmin) into a 64-bit key in place,
// using only as many bits as the bucket's value ranges need, then LSD radix sorted with the
// second half of the buffer as scratch, so no memory beyond the 16 bytes/record is needed.
static void hres_sort_unique(hraw_t *r, int64 n, hres_t *res)
{
    int64 i, j;
    MYBZERO(res, 1);
    res->buf = r;
    if (n == 0) return;

    uint32 amin, amax, bmin, bmax, apmin, apmax, bpmin, bpmax;
    amin = amax = r[0].a; bmin = bmax = r[0].b;
    apmin = apmax = r[0].ap; bpmin = bpmax = r[0].bp;
    for (i = 1; i < n; i++) {
        if (r[i].a < amin) amin = r[i].a;
        if (r[i].a > amax) amax = r[i].a;
        if (r[i].b < bmin) bmin = r[i].b;
        if (r[i].b > bmax) bmax = r[i].b;
        if (r[i].ap < apmin) apmin = r[i].ap;
        if (r[i].ap > apmax) apmax = r[i].ap;
        if (r[i].bp < bpmin) bpmin = r[i].bp;
        if (r[i].bp > bpmax) bpmax = r[i].bp;
    }
    int wa = bit_width32(amax - amin), wb = bit_width32(bmax - bmin);
    int wap = bit_width32(apmax - apmin), wbp = bit_width32(bpmax - bpmin);
    int tot = wa + wb + wap + wbp;

    if (tot <= 62) {
        uint64 *k = (uint64 *) r, *t = k + n, *src, *dst, *tmp, x;
        int sap = wbp, sb = wbp + wap, sa = wbp + wap + wb;
        hraw_t y;
        // pack in place: key i overwrites bytes of record i/2, which has already been read
        for (i = 0; i < n; i++) {
            memcpy(&y, (char *) r + i * sizeof(hraw_t), sizeof(hraw_t));
            x = (uint64) (y.a - amin) << sa | (uint64) (y.b - bmin) << sb |
                (uint64) (y.ap - apmin) << sap | (uint64) (y.bp - bpmin);
            memcpy((char *) r + i * sizeof(uint64), &x, sizeof(uint64));
        }

        int p, d, passes = (tot + RADIX_BITS - 1) / RADIX_BITS;
        int64 (*hist)[RADIX_SIZE] = NULL, s, c;
        if (passes > 0) {
            MYCALLOC(hist, passes);
            if (hist == NULL) mem_alloc_error("radix histogram");
            for (i = 0; i < n; i++) {
                x = k[i];
                for (p = 0; p < passes; p++)
                    hist[p][(x >> (p * RADIX_BITS)) & (RADIX_SIZE - 1)]++;
            }
        }
        src = k; dst = t;
        for (p = 0; p < passes; p++) {
            int64 *h = hist[p];
            int sh = p * RADIX_BITS, skip = 0;
            for (d = 0; d < RADIX_SIZE; d++)
                if (h[d] == n) { skip = 1; break; } // digit is constant: pass is a no-op
            if (skip) continue;
            for (s = 0, d = 0; d < RADIX_SIZE; d++) { c = h[d]; h[d] = s; s += c; }
            for (i = 0; i < n; i++) {
                x = src[i];
                dst[h[(x >> sh) & (RADIX_SIZE - 1)]++] = x;
            }
            tmp = src; src = dst; dst = tmp;
        }
        free(hist);

        // collapse duplicates; counts go into the other half of the buffer
        uint32 *cn = (uint32 *) dst;
        j = 0; cn[0] = 1;
        for (i = 1; i < n; i++) {
            if (src[i] == src[j]) {
                cn[j]++;
            } else {
                src[++j] = src[i];
                cn[j] = 1;
            }
        }
        res->packed = 1;
        res->key = src;
        res->cnt = cn;
        res->nu = j + 1;
        res->amin = amin; res->bmin = bmin; res->apmin = apmin; res->bpmin = bpmin;
        res->wb = wb; res->wap = wap; res->wbp = wbp;
    } else {
        // value ranges too wide to pack into 64 bits: comparison sort
        uint32 *cn;
        qsort(r, n, sizeof(hraw_t), hraw_cmpfunc);
        MYMALLOC(cn, n);
        if (cn == NULL) mem_alloc_error("count array");
        j = 0; cn[0] = 1;
        for (i = 1; i < n; i++) {
            if (hraw_cmpfunc(r + i, r + j) == 0) {
                cn[j]++;
            } else {
                r[++j] = r[i];
                cn[j] = 1;
            }
        }
        res->packed = 0;
        res->raw = r;
        res->cnt = cn;
        res->cbuf = cn;
        res->nu = j + 1;
    }
}

static void hres_emit(hout_t *out, hres_t *res)
{
    int64 i;
    hraw_t x;
    if (res->packed) {
        int sap = res->wbp, sb = res->wbp + res->wap, sa = sb + res->wb;
        uint64 k, mbp = (1ULL << res->wbp) - 1, map = (1ULL << res->wap) - 1, mb = (1ULL << res->wb) - 1;
        for (i = 0; i < res->nu; i++) {
            k = res->key[i];
            x.bp = (uint32) (k & mbp) + res->bpmin;
            x.ap = (uint32) ((k >> sap) & map) + res->apmin;
            x.b  = (uint32) ((k >> sb) & mb) + res->bmin;
            x.a  = (uint32) (k >> sa) + res->amin;
            hout_push(out, &x, res->cnt[i]);
        }
    } else {
        for (i = 0; i < res->nu; i++)
            hout_push(out, res->raw + i, res->cnt[i]);
    }
    free(res->buf);
    free(res->cbuf);
    MYBZERO(res, 1);
}

/* splitting into leaf buckets */

typedef struct {
    char *path;
    int64 n;
    int single; // all records share one key: no need to load
} leaf_t;
typedef kvec_t(leaf_t) leaf_v;

static inline uint32 hraw_key(const hraw_t *r, int level)
{
    switch (level) {
        case 0: return r->a;
        case 1: return r->b;
        case 2: return r->ap;
        default: return r->bp;
    }
}

// greedy grouping of contiguous non-empty keys into groups of at most gcap records
// (a single key may exceed it); returns the number of groups
static uint32 group_keys(const int64 *hist, uint32 krange, int64 gcap, uint32 *gid)
{
    uint32 k, ng = 0;
    int64 acc = 0;
    for (k = 0; k < krange; k++) {
        if (hist[k] == 0) { gid[k] = (uint32) -1; continue; }
        if (ng == 0 || (acc > 0 && acc + hist[k] > gcap)) { ng++; acc = 0; }
        gid[k] = ng - 1;
        acc += hist[k];
    }
    return ng;
}


// turn one bucket into an ordered list of leaves of at most `cap` records,
// re-splitting on the key of `level` when it is too big
static void collect_leaves(char *path, int64 n, int level, hout_t *out, int64 cap, leaf_v *lv)
{
    FILE *fp;
    hraw_t *r;
    int64 i, m;

    if (n == 0) {
        hfile_remove(path);
        free(path);
        return;
    }
    if (n <= cap || level > 3) {
        kv_push(leaf_t, *lv, ((leaf_t) {path, n, level > 3}));
        return;
    }

    hic_n_split++;
    hic_b_split += n * (int64) sizeof(hraw_t);

    // histogram on the current key
    uint32 krange, k, ng, g, nkey, *gid, *nk;
    int64 *hist;
    uint8 *gsingle;
    bucket_t *gbk;

    krange = level < 2 ? out->nseq : out->nwin;
    MYCALLOC(hist, krange);
    MYMALLOC(r, HRAW_RBUF);
    if (hist == NULL || r == NULL) mem_alloc_error("histogram");
    fp = fopen(path, "rb");
    if (fp == NULL) hfile_read_error(path);
    while ((m = fread(r, sizeof(hraw_t), HRAW_RBUF, fp)) > 0)
        for (i = 0; i < m; i++) {
            k = hraw_key(r + i, level);
            if (k >= krange) {
                fprintf(stderr, "[E::%s] key out of range (%u >= %u) at level %d\n", __func__, k, krange, level);
                exit(EXIT_FAILURE);
            }
            hist[k]++;
        }
    fclose(fp);

    // greedy grouping of contiguous keys; widen groups if too many files would be open
    int64 gcap = MAX(cap, (n + hic_max_buckets - 1) / hic_max_buckets * 2);
    MYMALLOC(gid, krange);
    if (gid == NULL) mem_alloc_error("group index");
    ng = group_keys(hist, krange, gcap, gid);

    if (ng == 1) {
        for (nkey = 0, k = 0; k < krange; k++)
            nkey += hist[k] > 0;
        // several keys, but the group size widened to limit open files swallowed
        // them all: split the keys in halves instead (>= 2 groups guaranteed)
        if (nkey > 1)
            ng = group_keys(hist, krange, (n + 1) / 2, gid);
    }

    if (ng == 1) {
        // a single key value: go one key deeper
        free(hist); free(gid); free(r);
        hic_n_split--;
        hic_b_split -= n * (int64) sizeof(hraw_t);
        collect_leaves(path, n, level + 1, out, cap, lv);
        return;
    }

    // a group holding a single key that is still too big is split on the next key
    MYCALLOC(gbk, ng);
    MYCALLOC(gsingle, ng);
    MYCALLOC(nk, ng);
    if (gbk == NULL || gsingle == NULL || nk == NULL) mem_alloc_error("groups");
    for (k = 0; k < krange; k++)
        if (gid[k] != (uint32) -1) nk[gid[k]]++;
    for (g = 0; g < ng; g++) gsingle[g] = (nk[g] == 1);
    free(nk);

    hic_set_wbuf(ng);
    for (g = 0; g < ng; g++) bucket_open(&gbk[g]);
    fp = fopen(path, "rb");
    if (fp == NULL) hfile_read_error(path);
    while ((m = fread(r, sizeof(hraw_t), HRAW_RBUF, fp)) > 0)
        for (i = 0; i < m; i++)
            bucket_push(&gbk[gid[hraw_key(r + i, level)]], r + i);
    fclose(fp);
    for (g = 0; g < ng; g++) bucket_close(&gbk[g]);
    hfile_remove(path);
    free(path);
    free(r); free(hist); free(gid);

    // groups are disjoint, contiguous key ranges: keep them in order
    for (g = 0; g < ng; g++)
        collect_leaves(gbk[g].path, gbk[g].n, (gsingle[g] && gbk[g].n > cap) ? level + 1 : level, out, cap, lv);
    free(gbk);
    free(gsingle);
}

/* parallel, ordered merge of leaves */
static void leaf_load_sort(leaf_t *lf, hres_t *res)
{
    FILE *fp;
    hraw_t *r;

    MYBZERO(res, 1);
    if (lf->single) {
        // all records share one key: read one and count them
        uint32 *cn;
        MYMALLOC(r, 1);
        MYMALLOC(cn, 1);
        if (r == NULL || cn == NULL) mem_alloc_error("leaf");
        fp = fopen(lf->path, "rb");
        if (fp == NULL || fread(r, sizeof(hraw_t), 1, fp) != 1) hfile_read_error(lf->path);
        fclose(fp);
        cn[0] = (uint32) MIN(lf->n, (int64) INT32_MAX);
        res->buf = r; res->raw = r;
        res->cbuf = cn; res->cnt = cn;
        res->nu = 1;
    } else {
        MYMALLOC(r, lf->n);
        if (r == NULL) mem_alloc_error("leaf records");
        fp = fopen(lf->path, "rb");
        if (fp == NULL || (int64) fread(r, sizeof(hraw_t), lf->n, fp) != lf->n) hfile_read_error(lf->path);
        fclose(fp);
        hres_sort_unique(r, lf->n, res);
    }
    hfile_remove(lf->path);
    free(lf->path);
    lf->path = NULL;
}

typedef struct {
    leaf_t *leaves;
    int64 nleaf;
    hres_t *res;
    uint8 *done;
    int64 next_claim, next_write;
    int window;
    pthread_mutex_t mu;
    pthread_cond_t cv;
} hpool_t;

static void *hpool_worker(void *arg)
{
    hpool_t *p = (hpool_t *) arg;
    int64 i;
    for (;;) {
        pthread_mutex_lock(&p->mu);
        // only work `window` leaves ahead of the writer: bounds memory to `window` leaves
        while (p->next_claim < p->nleaf && p->next_claim >= p->next_write + p->window)
            pthread_cond_wait(&p->cv, &p->mu);
        if (p->next_claim >= p->nleaf) {
            pthread_mutex_unlock(&p->mu);
            break;
        }
        i = p->next_claim++;
        pthread_mutex_unlock(&p->mu);

        leaf_load_sort(&p->leaves[i], &p->res[i]);

        pthread_mutex_lock(&p->mu);
        p->done[i] = 1;
        pthread_cond_broadcast(&p->cv);
        pthread_mutex_unlock(&p->mu);
    }
    return NULL;
}

// sort leaves with `nthreads` workers while this thread writes them out in order
static void merge_leaves(leaf_v *lv, hout_t *out, int nthreads)
{
    int64 i;
    int t;
    hpool_t p;
    pthread_t *tids;

    MYBZERO(&p, 1);
    p.leaves = lv->a;
    p.nleaf = lv->n;
    p.window = nthreads;
    MYCALLOC(p.res, MAX(lv->n, 1));
    MYCALLOC(p.done, MAX(lv->n, 1));
    MYMALLOC(tids, nthreads);
    if (p.res == NULL || p.done == NULL || tids == NULL) mem_alloc_error("thread pool");
    pthread_mutex_init(&p.mu, NULL);
    pthread_cond_init(&p.cv, NULL);

    for (t = 0; t < nthreads; t++)
        if (pthread_create(&tids[t], NULL, hpool_worker, &p) != 0) {
            fprintf(stderr, "[E::%s] cannot create worker thread\n", __func__);
            exit(EXIT_FAILURE);
        }

    for (i = 0; i < p.nleaf; i++) {
        pthread_mutex_lock(&p.mu);
        while (!p.done[i])
            pthread_cond_wait(&p.cv, &p.mu);
        pthread_mutex_unlock(&p.mu);

        hres_emit(out, &p.res[i]);

        pthread_mutex_lock(&p.mu);
        p.next_write = i + 1;
        pthread_cond_broadcast(&p.cv);
        pthread_mutex_unlock(&p.mu);
    }

    for (t = 0; t < nthreads; t++)
        pthread_join(tids[t], NULL);
    pthread_mutex_destroy(&p.mu);
    pthread_cond_destroy(&p.cv);
    free(tids);
    free(p.res);
    free(p.done);
}

/*************************************************************************
 * Entry parsing and scanning, shared by both loaders.
 *
 * is_alnh / conv are passed as literal constants through HIC_DISPATCH, and
 * the scan functions are force-inlined, so the compiler emits a separate,
 * branch-free loop for each of the 4 (format x conversion) combinations.
 *************************************************************************/

#define HIC_INLINE static inline __attribute__((always_inline))

// call fn(..., is_alnh, conv) with both flags as compile-time constants
#define HIC_DISPATCH(fn, is_alnh, conv, ...) do { \
        if (is_alnh) { if (conv) fn(__VA_ARGS__, 1, 1); else fn(__VA_ARGS__, 1, 0); } \
        else         { if (conv) fn(__VA_ARGS__, 0, 1); else fn(__VA_ARGS__, 0, 0); } \
    } while (0)

// open a HiC bin file, check the header and leave fp at the first record
static FILE *hic_bin_open(const char *hic_bfile, int *is_alnh, uint64 *pair_n)
{
    int64 magic_number;
    FILE *fp = fopen(hic_bfile, "rb");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open %s\n", __func__, hic_bfile);
        return NULL;
    }
    if (fread(&magic_number, sizeof(int64), 1, fp) != 1) bin_fread_error();
    if (is_valid_alnh_bin_header(magic_number)) {
        *is_alnh = 1;
    } else if (is_valid_yahs_bin_header(magic_number)) {
        *is_alnh = 0;
    } else {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        fclose(fp);
        return NULL;
    }
    binary_fseek_skip_sdict(fp);
    if (fread(pair_n, sizeof(uint64), 1, fp) != 1) bin_fread_error();
    return fp;
}

// parse one entry into a canonical record (aseq <= bseq; apos <= bpos if same seq)
// returns 0 if the entry is filtered out
HIC_INLINE int hic_parse_entry(const uint8 *bbuff, asm_dict_t *adict, uint32 nseq, int bin_size,
    uint8 min_qual, hraw_t *r, const int is_alnh, const int conv)
{
    const int esize = is_alnh ? ALNHBIN_ENTRY_SIZE : YAHSBIN_ENTRY_SIZE;
    uint32 i0, i1, s0, s1, e0, e1, q0, q1;
    uint64 p0, p1;

    if (bbuff[esize - 1] < min_qual)
        return 0;
    if (is_alnh) {
        memcpy(&i0, bbuff, sizeof(uint32));
        memcpy(&s0, bbuff + sizeof(uint32) * 1, sizeof(uint32));
        memcpy(&e0, bbuff + sizeof(uint32) * 2, sizeof(uint32));
        memcpy(&i1, bbuff + sizeof(uint32) * 3, sizeof(uint32));
        memcpy(&s1, bbuff + sizeof(uint32) * 4, sizeof(uint32));
        memcpy(&e1, bbuff + sizeof(uint32) * 5, sizeof(uint32));
        s0 = (s0 >> 1) + (e0 >> 1) + ((s0 & e0) & 1);
        s1 = (s1 >> 1) + (e1 >> 1) + ((s1 & e1) & 1);
    } else {
        memcpy(&i0, bbuff, sizeof(uint32));
        memcpy(&s0, bbuff + sizeof(uint32) * 1, sizeof(uint32));
        memcpy(&i1, bbuff + sizeof(uint32) * 2, sizeof(uint32));
        memcpy(&s1, bbuff + sizeof(uint32) * 3, sizeof(uint32));
    }
    if (conv) {
        if (sd_coordinate_conversion(adict, i0, s0, &q0, &p0, 1) != CC_SUCCESS ||
            sd_coordinate_conversion(adict, i1, s1, &q1, &p1, 1) != CC_SUCCESS)
            return 0;
    } else {
        q0 = i0; p0 = s0;
        q1 = i1; p1 = s1;
    }
    if (q0 >= nseq || q1 >= nseq) {
        fprintf(stderr, "[E::%s] sequence id out of range (%u, %u >= %u)\n", __func__, q0, q1, nseq);
        exit(EXIT_FAILURE);
    }
    p0 /= bin_size;
    p1 /= bin_size;
    if (q0 > q1) {
        SWAP(uint32, q0, q1);
        SWAP(uint64, p0, p1);
    }
    if (q0 == q1 && p0 > p1)
        SWAP(uint64, p0, p1);
    *r = (hraw_t) {q0, (uint32) p0, q1, (uint32) p1};
    return 1;
}

typedef struct {
    FILE *fp;
    uint64 pair_n, pair_c, n_recs;
    asm_dict_t *adict;
    uint32 nseq;
    int bin_size;
    uint8 min_qual;
    hicType_t type;
    int symm;
} hic_scan_t;

/* in-memory loader */
typedef struct {
    hic_t *a;
    uint64 n;
} hic_arr_t;

static inline int hic_keep_type(const hic_scan_t *sc, const hraw_t *r)
{
    return sc->type == HIC_ALL || ((r->a == r->b) == (sc->type == HIC_INTRA));
}

static int hic_check_type(hicType_t type)
{
    if (type != HIC_ALL && type != HIC_INTRA && type != HIC_INTER) {
        fprintf(stderr, "[E::%s] invalid link type %d\n", __func__, (int) type);
        return 0;
    }
    return 1;
}

static inline uint32 hic_target_nseq(const sdict_t *dicts, const asm_dict_t *adict)
{
    return adict ? (uint32) adict->n : (uint32) dicts->n;
}

static inline uint64 hic_target_len(const sdict_t *dicts, const asm_dict_t *adict, uint32 i)
{
    return adict ? (uint64) adict->s[i].len : (uint64) dicts->s[i].len;
}

// periodic in-place sort-merge, buffer doubles when merging does not free half of it
HIC_INLINE void hic_inmem_scan(hic_scan_t *sc, hic_arr_t *out, const int is_alnh, const int conv)
{
    const int esize = is_alnh ? ALNHBIN_ENTRY_SIZE : YAHSBIN_ENTRY_SIZE;
    uint64 m, n, max_n;
    uint8 *bbuff, *ebuff, *buffer;
    hic_t *hics;
    hraw_t r;

    MYMALLOC(buffer, esize * FREAD_BUFF_SIZE);
    if (buffer == NULL)
        mem_alloc_error("buffer array");
    max_n = 0x4000000ULL; // 64M entries - ~1.25 GB for the hic_t array
    if (max_n > sc->pair_n)
        max_n = MAX(sc->pair_n, 1);
    MYMALLOC(hics, max_n);
    if (hics == NULL)
        mem_alloc_error("hic array");

    n = 0;
    while (sc->pair_c < sc->pair_n) {
        m = fread(buffer, esize, FREAD_BUFF_SIZE, sc->fp);
        if (m == 0) {
            fprintf(stderr, "[W::%s] unexpected end of file after %llu records\n", __func__, sc->pair_c);
            break;
        }
        for (bbuff = buffer, ebuff = buffer + m * esize;
            bbuff < ebuff && sc->pair_c++ < sc->pair_n; bbuff += esize) {
            if (sc->pair_c % REPORT_EPOC == 0)
                fprintf(stderr, "[M::%s] %llu million records processed\n", __func__, sc->pair_c / 1000000);
            if (!hic_parse_entry(bbuff, sc->adict, sc->nseq, sc->bin_size, sc->min_qual, &r, is_alnh, conv))
                continue;
            if (!hic_keep_type(sc, &r))
                continue;
            hics[n++] = (hic_t) {r.a, r.ap, r.b, r.bp, 1};
            if (n == max_n) {
                n = hic_sort_merge(hics, n);
                if (n > (max_n >> 1)) {
                    max_n <<= 1;
                    MYREALLOC(hics, max_n);
                    if (hics == NULL)
                        mem_alloc_error("hic array");
                }
            }
            sc->n_recs++;
        }
    }
    free(buffer);

    out->n = hic_sort_merge(hics, n);
    out->a = hics;
}

// add symmetric HiC links for inter-chromosomal contacts
static hic_t *add_symmetric_hic_links(hic_t *hics, int64 *_nhic) 
{
    int64 i, mhic, nhic;
    
    mhic = nhic = *_nhic;
    for (i = 0; i < nhic; i++)
        mhic += (hics[i].aseq != hics[i].bseq);
    MYREALLOC(hics, mhic);
    if (hics == NULL)
        mem_alloc_error("HiC links array");
    mhic = nhic;
    for (i = 0; i < mhic; i++)
        if (hics[i].aseq != hics[i].bseq)
            hics[nhic++] = (hic_t){hics[i].bseq, hics[i].bpos, hics[i].aseq, hics[i].apos, hics[i].nhic};
    // sort the HiC links
    qsort(hics, nhic, sizeof(hic_t), hic_cmpfunc);

    *_nhic = nhic;
    return hics;
}

// dicts gives the number of target sequences (may be NULL to skip the range check);
// adict non-NULL turns on coordinate conversion;
// type selects intra, inter or all links; symm adds the mirrored copy of inter links
// output is sorted in hic_cmpfunc order.
hic_t *read_hic_inmem(char *hic_bfile, sdict_t *dicts, asm_dict_t *adict, int bin_size, uint8 min_qual,
    hicType_t type, int symm, int64 *_nhic)
{
    if (_nhic) *_nhic = 0;
    if (hic_bfile == NULL || !hic_check_type(type))
        return NULL;

    int is_alnh;
    hic_scan_t sc;
    hic_arr_t out;

    MYBZERO(&sc, 1);
    sc.fp = hic_bin_open(hic_bfile, &is_alnh, &sc.pair_n);
    if (sc.fp == NULL)
        return NULL;
    sc.adict = adict;
    sc.nseq = hic_target_nseq(dicts, adict);
    sc.bin_size = bin_size;
    sc.min_qual = min_qual;
    sc.type = type;
    sc.symm = 0;

    HIC_DISPATCH(hic_inmem_scan, is_alnh, adict != NULL, &sc, &out);
    fclose(sc.fp);

    fprintf(stderr, "[M::%s] processed %llu read pairs\n", __func__, sc.pair_c);
    fprintf(stderr, "[M::%s] retained %llu records\n", __func__, sc.n_recs);
    fprintf(stderr, "[M::%s] linked %llu sequence regions\n", __func__, out.n);

    if (out.n == 0) {
        free(out.a);
        return NULL;
    }

    int64 nhic = (int64) out.n;
    if (symm && type != HIC_INTRA)
        out.a = add_symmetric_hic_links(out.a, &nhic);
    if (_nhic) *_nhic = nhic;
    
    return out.a;
}

/* external-memory loader */
// stream entries into memory (mem != NULL) or aseq buckets, adding the symmetric copy
HIC_INLINE void hic_extern_scan(hic_scan_t *sc, hraw_t *mem, int64 *_nmem, bucket_t *bks, const uint32 *part,
    const int is_alnh, const int conv)
{
    const int esize = is_alnh ? ALNHBIN_ENTRY_SIZE : YAHSBIN_ENTRY_SIZE;
    uint64 m;
    int64 nmem = *_nmem;
    uint8 *bbuff, *ebuff, *buffer;
    hraw_t r, x;

    MYMALLOC(buffer, esize * FREAD_BUFF_SIZE);
    if (buffer == NULL)
        mem_alloc_error("buffer array");

    while (sc->pair_c < sc->pair_n) {
        m = fread(buffer, esize, FREAD_BUFF_SIZE, sc->fp);
        if (m == 0) {
            fprintf(stderr, "[W::%s] unexpected end of file after %llu records\n", __func__, sc->pair_c);
            break;
        }
        for (bbuff = buffer, ebuff = buffer + m * esize;
            bbuff < ebuff && sc->pair_c++ < sc->pair_n; bbuff += esize) {
            if (sc->pair_c % REPORT_EPOC == 0)
                fprintf(stderr, "[M::%s] %llu million records processed\n", __func__, sc->pair_c / 1000000);
            if (!hic_parse_entry(bbuff, sc->adict, sc->nseq, sc->bin_size, sc->min_qual, &r, is_alnh, conv))
                continue;
            if (!hic_keep_type(sc, &r))
                continue;
            if (mem) mem[nmem++] = r;
            else bucket_push(&bks[part[r.a]], &r);
            if (sc->symm && r.a != r.b) {
                x = (hraw_t) {r.b, r.bp, r.a, r.ap};
                if (mem) mem[nmem++] = x;
                else bucket_push(&bks[part[x.a]], &x);
            }
            sc->n_recs++;
        }
    }
    free(buffer);
    *_nmem = nmem;
}

// external-memory loader; arguments as for read_hic_inmem. Large results are mmap'ed
// from a temp file; release them with hic_store_release().
static hic_t *read_hic_extern_impl(char *hic_bfile, sdict_t *dicts, asm_dict_t *adict, int bin_size, uint8 min_qual,
    hicType_t type, int symm, int64 *_nhic, hic_store_t *store, int map)
{
    if (_nhic) *_nhic = 0;
    MYBZERO(store, 1);
    if (hic_bfile == NULL || !hic_check_type(type))
        return NULL;
    if (dicts == NULL && adict == NULL)
        return NULL;

    int is_alnh, nthreads;
    uint32 nseq, nb, *part;
    uint64 m, tlen, clen;
    int64 i, cap, est, nmem;
    hraw_t *mem;
    bucket_t *bks;
    hic_scan_t sc;
    hout_t out;
    leaf_v lv;

    MYBZERO(&sc, 1);
    sc.fp = hic_bin_open(hic_bfile, &is_alnh, &sc.pair_n);
    if (sc.fp == NULL)
        return NULL;

    // ids and positions after conversion refer to adict's sequences: size everything
    // indexed by sequence id or window (bounds check, partition, histogram ranges) from them
    nseq = hic_target_nseq(dicts, adict);
    sc.adict = adict;
    sc.nseq = nseq;
    sc.bin_size = bin_size;
    sc.min_qual = min_qual;
    sc.type = type;
    sc.symm = symm && type != HIC_INTRA;

    MYBZERO(&out, 1);
    out.nseq = nseq;
    tlen = 0;
    for (i = 0; i < nseq; i++) {
        tlen += hic_target_len(dicts, adict, i);
        m = hic_target_len(dicts, adict, i) / bin_size + 2;
        if (m > out.nwin) out.nwin = m;
    }

    nthreads = HIC_THREADS;
    if (nthreads <= 0) {
        long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
        nthreads = (int) MAX(1, MIN(4, ncpu));
    }
    hic_set_file_limit();

    // records per leaf: `nthreads` leaves of 16 bytes/record may be in memory at once
    cap = (HIC_MEM_LIMIT - (int64) HOUT_WBUF * sizeof(hic_t)) / ((int64) sizeof(hraw_t) * nthreads);
    if (cap < (1 << 16)) cap = 1 << 16;
    // upper bound on raw records: each pair, plus its mirror if symmetric
    est = (int64) sc.pair_n * (sc.symm ? 2 : 1);

    MYCALLOC(part, nseq);
    if (part == NULL) mem_alloc_error("partition index");

    mem = NULL; nmem = 0; bks = NULL; nb = 0;
    if (est * (int64) (sizeof(hraw_t) + sizeof(hic_t)) <= HIC_MEM_LIMIT) {
        // small input: do everything in memory
        MYMALLOC(mem, est > 0 ? est : 1);
        if (mem == NULL) mem_alloc_error("HiC raw array");
    } else {
        // partition by aseq, sizing buckets in proportion to sequence length
        nb = (uint32) MIN((int64) hic_max_buckets, MAX((int64) 1, (est + cap - 1) / cap * 2));
        clen = 0;
        for (i = 0; i < nseq; i++) {
            part[i] = (uint32) MIN((uint64) nb - 1, clen * nb / MAX(tlen, 1));
            clen += hic_target_len(dicts, adict, i);
        }
        hic_set_wbuf(nb);
        MYCALLOC(bks, nb);
        if (bks == NULL) mem_alloc_error("buckets");
        for (i = 0; i < nb; i++) bucket_open(&bks[i]);
    }

    HIC_DISPATCH(hic_extern_scan, is_alnh, adict != NULL, &sc, mem, &nmem, bks, part);
    fclose(sc.fp);
    free(part);

    fprintf(stderr, "[M::%s] processed %llu read pairs\n", __func__, sc.pair_c);
    fprintf(stderr, "[M::%s] retained %llu records\n", __func__, sc.n_recs);

    if (mem) {
        hres_t res;
        hres_sort_unique(mem, nmem, &res); // takes over mem
        hres_emit(&out, &res);
        MYREALLOC(out.a, MAX(out.n, 1));
    } else {
        for (i = 0; i < nb; i++) bucket_close(&bks[i]);
        kv_init(lv);
        hic_n_split = hic_b_split = 0;
        for (i = 0; i < nb; i++)
            collect_leaves(bks[i].path, bks[i].n, 0, &out, cap, &lv);
        free(bks);

        store->path = hic_temp_path("hic");
        out.fp = fopen(store->path, "wb+");
        if (out.fp == NULL) {
            fprintf(stderr, "[E::%s] cannot create %s: %s\n", __func__, store->path, strerror(errno));
            exit(EXIT_FAILURE);
        }
        out.m = HOUT_WBUF;
        MYMALLOC(out.a, out.m);
        if (out.a == NULL) mem_alloc_error("output buffer");
        merge_leaves(&lv, &out, nthreads);
        hout_flush(&out);
        free(out.a);
        out.a = NULL;
        if (fflush(out.fp) != 0 || ferror(out.fp)) {
            fprintf(stderr, "[E::%s] write error on %s: %s (disk full?)\n", __func__, store->path, strerror(errno));
            exit(EXIT_FAILURE);
        }
        kv_destroy(lv);
    }
    fprintf(stderr, "[M::%s] linked %lld sequence regions\n", __func__, out.ncanon);

    if (out.nout == 0) {
        if (out.fp) { fclose(out.fp); hfile_remove(store->path); free(store->path); store->path = NULL; }
        hic_run_dir_remove();
        free(out.a);
        return NULL;
    }

    if (_nhic) *_nhic = out.nout;
    if (out.fp == NULL)
        return out.a;

    if (!map) {
        // streaming: keep the sorted file open; it is removed by hic_store_release()
        store->fp = out.fp;
        return NULL;
    }

    // map the sorted file read-only; downstream code only reads it
    store->map_len = (size_t) out.nout * sizeof(hic_t);
    store->map = mmap(NULL, store->map_len, PROT_READ, MAP_SHARED, fileno(out.fp), 0);
    fclose(out.fp);
    if (store->map == MAP_FAILED) {
        fprintf(stderr, "[E::%s] mmap failed on %s: %s\n", __func__, store->path, strerror(errno));
        exit(EXIT_FAILURE);
    }
    // unlink now: the mapping stays valid and the file disappears even on a crash
    if (unlink(store->path) == 0) {
        free(store->path);
        store->path = NULL;
        hic_run_dir_remove();
    }
    return (hic_t *) store->map;
}

hic_t *read_hic_extern(char *hic_bfile, sdict_t *dicts, asm_dict_t *adict, int bin_size, uint8 min_qual,
    hicType_t type, int symm, int64 *_nhic, hic_store_t *store)
{
    return read_hic_extern_impl(hic_bfile, dicts, adict, bin_size, min_qual, type, symm, _nhic, store, 1);
}

void hic_store_release(hic_store_t *store, hic_t *hics)
{
    if (store->map || store->fp) {
        if (store->map) munmap(store->map, store->map_len);
        if (store->fp) fclose(store->fp);
        if (store->path) { hfile_remove(store->path); free(store->path); }
        hic_run_dir_remove();
    } else {
        free(hics);
    }
    MYBZERO(store, 1);
}

/****************************************************************************/
/*
 * Row stream over the sorted, symmetric HiC links
 *
 * The links are either an in-memory array or the sorted temp file written by
 * the external loader; either way they are read front to back in chunks, so
 * the downstream code never needs the whole set in memory or random access.
 *
 ****************************************************************************/

#define HSTREAM_CHUNK (1 << 20) // records per read (20 MB)

static void hic_stream_init_array(hic_stream_t *st, hic_t *a, int64 n)
{
    MYBZERO(st, 1);
    st->arr = a;
    st->n = a ? n : 0;
}

void hic_stream_rewind(hic_stream_t *st)
{
    st->pos = 0;
    st->nread = 0;
    if (st->store.fp) {
        if (fseeko(st->store.fp, 0, SEEK_SET) != 0) {
            fprintf(stderr, "[E::%s] cannot rewind %s: %s\n", __func__, st->store.path, strerror(errno));
            exit(EXIT_FAILURE);
        }
#ifdef POSIX_FADV_SEQUENTIAL
        posix_fadvise(fileno(st->store.fp), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    }
}

// next chunk of records; *h stays valid until the next call; returns 0 at the end
int64 hic_stream_read(hic_stream_t *st, hic_t **h)
{
    int64 m;
    if (st->store.fp == NULL) {
        m = st->n - st->pos;
        *h = st->arr + st->pos;
        st->pos = st->n;
        return m;
    }
    if (st->nread >= st->n)
        return 0;
    if (st->buf == NULL) {
        MYMALLOC(st->buf, HSTREAM_CHUNK);
        if (st->buf == NULL) mem_alloc_error("stream buffer");
    }
    m = MIN((int64) HSTREAM_CHUNK, st->n - st->nread);
    if ((int64) fread(st->buf, sizeof(hic_t), m, st->store.fp) != m) {
        fprintf(stderr, "[E::%s] error reading %s: %s\n", __func__, st->store.path, strerror(errno));
        exit(EXIT_FAILURE);
    }
    st->nread += m;
    *h = st->buf;
    return m;
}

void hic_stream_close(hic_stream_t *st)
{
    if (st->store.fp || st->store.map)
        hic_store_release(&st->store, NULL);
    else
        free(st->arr);
    free(st->buf);
    MYBZERO(st, 1);
}

// external loader that leaves the result as a stream over the sorted temp file
// (or over an in-memory array for small inputs); returns the number of records
int64 read_hic_extern_stream(char *hic_bfile, sdict_t *dicts, asm_dict_t *adict,int bin_size, uint8 min_qual,
    hicType_t type, int symm, hic_stream_t *st)
{
    int64 n = 0;
    hic_store_t store;
    hic_t *a;

    MYBZERO(st, 1);
    a = read_hic_extern_impl(hic_bfile, dicts, adict, bin_size, min_qual, type, symm, &n, &store, 0);
    if (store.fp) {
        st->store = store;
        st->n = n;
    } else hic_stream_init_array(st, a, n);
    return n;
}
