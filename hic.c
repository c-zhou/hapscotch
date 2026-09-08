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
    rec_c = pair_c = inter_c = intra_c = 0;
    buff = 0;
    b = bam_init1();
    
    fwrite(&pair_c, sizeof(uint64), 1, fo);

    while (bam_read1(fp, b) >= 0 ) {
        ++rec_c;

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

        if (rec_c % REPORT_EPOC == 0)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs\n", __func__, rec_c / 1000000, pair_c);
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

    rec_c = pair_c = inter_c = intra_c = 0;
    buff = 0;
    
    fwrite(&pair_c, sizeof(uint64), 1, fo);

    while ((line = iostream_getline(fp)) != NULL) {
        if (is_empty_line(line) || parse_line(line, fields, 5) < 5)
            continue;
        
        ++rec_c;

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

        if (rec_c % REPORT_EPOC == 0)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs\n", __func__, rec_c / 1000000, pair_c);
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

        ++rec_c;

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

        if (rec_c % REPORT_EPOC == 0)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs\n", __func__, rec_c / 1000000, pair_c);
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

        if (++rec_c % REPORT_EPOC == 0)
            fprintf(stderr, "[M::%s] %llu million records processed, %llu read pairs\n", __func__, rec_c / 1000000, pair_c);
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

void write_binary_hic_data_pseudo_yahs(char *bf, sdict_t *dicts, asm_dict_t *break_dict, int32 **smap, char *out)
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
            if (break_dict &&
                (piece_position(break_dict, i0, s0, e0) ||
                 piece_position(break_dict, i1, s1, e1)))
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

hic_t *read_hic_from_binary(char *hic_bfile, sdict_t *dicts, int bin_size, uint8 min_qual, int64 *_nhic)
{
    if (_nhic) *_nhic = 0;
    if (hic_bfile == NULL)
        return NULL;

    uint32 i0, i1, s0, s1, e0, e1;
    uint64 i, j, m, n, pair_n, pair_c, max_n, n_recs;
    int64 magic_number;
    uint8 *bbuff, *ebuff, *buffer;
    FILE *fp;
    hic_t *hics;

    fp = fopen(hic_bfile, "r");
    if (fp == NULL)
        return NULL;

    m = fread(&magic_number, sizeof(int64), 1, fp);
    if (m != 1) bin_fread_error();
    if (!is_valid_alnh_bin_header(magic_number)) {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        return NULL;
    }
    binary_fseek_skip_sdict(fp);
    m = fread(&pair_n, sizeof(uint64), 1, fp);
    if (m != 1) bin_fread_error();
    
    MYMALLOC(buffer, ALNHBIN_ENTRY_SIZE * FREAD_BUFF_SIZE);
    if (buffer == NULL)
        mem_alloc_error("buffer array");

    max_n = 0x4000000ULL; // 64MB - ~ 1GB mem for hic_t array
    if (max_n > pair_n)
        max_n = pair_n;
    MYMALLOC(hics, max_n);
    if (hics == NULL)
        mem_alloc_error("hic array");

    pair_c = 0;
    n_recs = 0;
    n = 0;
    while (pair_c < pair_n) {
        m = fread(buffer, sizeof(uint8) * ALNHBIN_ENTRY_SIZE, FREAD_BUFF_SIZE, fp);
        for (bbuff = buffer, ebuff = buffer + m * ALNHBIN_ENTRY_SIZE; 
            bbuff < ebuff && pair_c++ < pair_n; bbuff += ALNHBIN_ENTRY_SIZE) {
            if (bbuff[ALNHBIN_ENTRY_SIZE-1] < min_qual)
                continue;
            memcpy(&i0, bbuff,  sizeof(uint32));
            memcpy(&s0, bbuff + sizeof(uint32)*1, sizeof(uint32));
            memcpy(&e0, bbuff + sizeof(uint32)*2, sizeof(uint32));
            memcpy(&i1, bbuff + sizeof(uint32)*3, sizeof(uint32));
            memcpy(&s1, bbuff + sizeof(uint32)*4, sizeof(uint32));
            memcpy(&e1, bbuff + sizeof(uint32)*5, sizeof(uint32));

            s0 = (s0>>1) + (e0>>1) + ((s0 & e0) & 1);
            s1 = (s1>>1) + (e1>>1) + ((s1 & e1) & 1);

            s0 /= bin_size;
            s1 /= bin_size;

            if (i0 == i1 && s0 > s1) SWAP(uint32, s0, s1);

            hics[n++] = (hic_t){i0, s0, i1, s1, 1};
            if (n == max_n) {
                n = hic_sort_merge(hics, n);
                if (n > (max_n>>1)) {
                    max_n <<= 1;
                    MYREALLOC(hics, max_n);
                    if (hics == NULL)
                        mem_alloc_error("hic array");
                }
            }
            n_recs++;
        }
    }
    
    n = hic_sort_merge(hics, n);
    
    fprintf(stderr, "[M::%s] processed %llu read pairs\n", __func__, pair_c);
    fprintf(stderr, "[M::%s] retained %llu records\n", __func__, n_recs);
    fprintf(stderr, "[M::%s] linked %llu sequence regions\n", __func__, n);
    
    free(buffer);
    fclose(fp);

    if (n == 0) {
        free(hics);
        return NULL;
    }

    // add symmetric links
    for (i = 0, m = n; i < m; i++) {
        if (hics[i].aseq == hics[i].bseq)
            continue;
        if (n == max_n) {
            max_n <<= 1;
            MYREALLOC(hics, max_n);
            if (hics == NULL)
                mem_alloc_error("hic array");
        }
        hics[n++] = (hic_t){hics[i].bseq, hics[i].bpos, hics[i].aseq, hics[i].apos, hics[i].nhic};
    }
    qsort(hics, n, sizeof(hic_t), hic_cmpfunc);

    MYREALLOC(hics, n);
    
    // some statistics for hic links
    int64 slen, all_w, cov_w, inter_w, intra_w;
    int64 inter_c, intra_c, total_c;
    uint8 *wmark;
    int w, s;

    slen = 0;
    all_w = 0;
    max_n = 1;
    for (i = 0; i < dicts->n; i++) {
        w = dicts->s[i].len;
        slen += w;
        w = (w - 1) / bin_size + 1;
        all_w += w;
        max_n = MAX(max_n, w);
    }
    // number of windows
    MYCALLOC(wmark, max_n);
    if (wmark == NULL)
        mem_alloc_error("mark array");
    
    inter_c = intra_c = 0;
    inter_w = intra_w = 0;
    cov_w = 0;
    s = hics[0].aseq;
    for (i = 0; i <= n; i++) {
        if (i == n || hics[i].aseq != s) {
            m = (dicts->s[s].len - 1) / bin_size + 1;
            for (j = 0; j < m; j++) {
                if (wmark[j])
                    cov_w++;
                if (wmark[j] & 0x1)
                    intra_w++;
                if (wmark[j] & 0x2)
                    inter_w++;
            }
            if (i == n) break;
            s = hics[i].aseq;
            MYBZERO(wmark, max_n);
        }
        if (hics[i].bseq == s) {
            intra_c += hics[i].nhic;
            wmark[hics[i].apos] |= 0x1;
            wmark[hics[i].bpos] |= 0x1;
        } else {
            inter_c += hics[i].nhic;
            wmark[hics[i].apos] |= 0x2;
        }
    }
    inter_c /= 2; // each inter link is counted twice
    total_c = inter_c + intra_c;

    fprintf(stderr, "[M::%s] HiC links summary statistics\n", __func__);
    fprintf(stderr, "[M::%s] total sequence: %12lld\n", __func__, slen);
    fprintf(stderr, "[M::%s]    window size: %12d\n", __func__, bin_size);
    fprintf(stderr, "[M::%s]    no. windows: %12lld\n", __func__, all_w);
    fprintf(stderr, "[M::%s] hic link counts\n", __func__);
    fprintf(stderr, "[M::%s]        - total: %12lld\n", __func__, total_c);
    fprintf(stderr, "[M::%s]        - intra: %12lld\n", __func__, intra_c);
    fprintf(stderr, "[M::%s]        - inter: %12lld\n", __func__, inter_c);
    fprintf(stderr, "[M::%s] sequence window\n", __func__);
    fprintf(stderr, "[M::%s]        - total: %12lld\n", __func__, cov_w);
    fprintf(stderr, "[M::%s]        - intra: %12lld\n", __func__, intra_w);
    fprintf(stderr, "[M::%s]        - inter: %12lld\n", __func__, inter_w);
    fprintf(stderr, "[M::%s] links per window\n", __func__);
    fprintf(stderr, "[M::%s]        - total: %12.3e\n", __func__, (double) total_c/cov_w);
    fprintf(stderr, "[M::%s]        - intra: %12.3e\n", __func__, (double) intra_c/cov_w);
    fprintf(stderr, "[M::%s]        - inter: %12.3e\n", __func__, (double) inter_c/cov_w);
    fprintf(stderr, "[M::%s] links per base: %12.3e\n", __func__, (double) total_c/slen);

    free(wmark);

    if (_nhic) *_nhic = n;

    return hics;
}

static hic_t *read_hic_from_binary_sd_conversion_alnhbin(char *hic_bfile, asm_dict_t *dicts, int bin_size, uint8 min_qual, int64 *_nhic)
{
    if (_nhic) *_nhic = 0;
    if (hic_bfile == NULL)
        return NULL;

    uint32 i0, i1, s0, s1, e0, e1, q0, q1;
    uint64 m, n, p0, p1, pair_n, pair_c, max_n, n_recs;
    int64 magic_number;
    uint8 *bbuff, *ebuff, *buffer;
    FILE *fp;
    hic_t *hics;

    fp = fopen(hic_bfile, "r");
    if (fp == NULL)
        return NULL;

    m = fread(&magic_number, sizeof(int64), 1, fp);
    if (m != 1) bin_fread_error();
    if (!is_valid_alnh_bin_header(magic_number)) {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        return NULL;
    }
    binary_fseek_skip_sdict(fp);
    m = fread(&pair_n, sizeof(uint64), 1, fp);
    if (m != 1) bin_fread_error();
    
    MYMALLOC(buffer, ALNHBIN_ENTRY_SIZE * FREAD_BUFF_SIZE);
    if (buffer == NULL)
        mem_alloc_error("buffer array");

    max_n = 0x4000000ULL; // 64MB - ~ 1GB mem for hic_t array
    if (max_n > pair_n)
        max_n = pair_n;
    MYMALLOC(hics, max_n);
    if (hics == NULL)
        mem_alloc_error("hic array");

    pair_c = 0;
    n_recs = 0;
    n = 0;
    while (pair_c < pair_n) {
        m = fread(buffer, sizeof(uint8) * ALNHBIN_ENTRY_SIZE, FREAD_BUFF_SIZE, fp);
        for (bbuff = buffer, ebuff = buffer + m * ALNHBIN_ENTRY_SIZE; 
            bbuff < ebuff && pair_c++ < pair_n; bbuff += ALNHBIN_ENTRY_SIZE) {
            if (pair_c % REPORT_EPOC == 0)
                fprintf(stderr, "[M::%s] %llu million records processed\n", __func__, pair_c / 1000000);

            if (bbuff[ALNHBIN_ENTRY_SIZE-1] < min_qual)
                continue;
            memcpy(&i0, bbuff,  sizeof(uint32));
            memcpy(&s0, bbuff + sizeof(uint32)*1, sizeof(uint32));
            memcpy(&e0, bbuff + sizeof(uint32)*2, sizeof(uint32));
            memcpy(&i1, bbuff + sizeof(uint32)*3, sizeof(uint32));
            memcpy(&s1, bbuff + sizeof(uint32)*4, sizeof(uint32));
            memcpy(&e1, bbuff + sizeof(uint32)*5, sizeof(uint32));

            s0 = (s0>>1) + (e0>>1) + ((s0 & e0) & 1);
            s1 = (s1>>1) + (e1>>1) + ((s1 & e1) & 1);

            if (sd_coordinate_conversion(dicts, i0, s0, &q0, &p0, 1) != CC_SUCCESS ||
                sd_coordinate_conversion(dicts, i1, s1, &q1, &p1, 1) != CC_SUCCESS)
                continue;

            p0 /= bin_size;
            p1 /= bin_size;

            if (q0 > q1) {
                SWAP(uint32, q0, q1);
                SWAP(uint64, p0, p1);
            }
            if (q0 == q1 && p0 > p1)
                SWAP(uint64, p0, p1);

            hics[n++] = (hic_t){q0, p0, q1, p1, 1};
            if (n == max_n) {
                n = hic_sort_merge(hics, n);
                if (n > (max_n>>1)) {
                    max_n <<= 1;
                    MYREALLOC(hics, max_n);
                    if (hics == NULL)
                        mem_alloc_error("hic array");
                }
            }
            n_recs++;
        }
    }
    
    n = hic_sort_merge(hics, n);
    
    fprintf(stderr, "[M::%s] processed %llu read pairs\n", __func__, pair_c);
    fprintf(stderr, "[M::%s] retained %llu records\n", __func__, n_recs);
    fprintf(stderr, "[M::%s] linked %llu sequence regions\n", __func__, n);
    
    free(buffer);
    fclose(fp);

    if (n == 0) {
        free(hics);
        return NULL;
    }
    MYREALLOC(hics, n);

    if (_nhic) *_nhic = n;

    return hics;
}

static hic_t *read_hic_from_binary_sd_conversion_yahsbin(char *hic_bfile, asm_dict_t *dicts, int bin_size, uint8 min_qual, int64 *_nhic)
{
    if (_nhic) *_nhic = 0;
    if (hic_bfile == NULL)
        return NULL;

    uint32 i0, i1, s0, s1, q0, q1;
    uint64 m, n, p0, p1, pair_n, pair_c, max_n, n_recs;
    int64 magic_number;
    uint8 *bbuff, *ebuff, *buffer;
    FILE *fp;
    hic_t *hics;

    fp = fopen(hic_bfile, "r");
    if (fp == NULL)
        return NULL;

    m = fread(&magic_number, sizeof(int64), 1, fp);
    if (m != 1) bin_fread_error();
    if (!is_valid_yahs_bin_header(magic_number)) {
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
        return NULL;
    }
    binary_fseek_skip_sdict(fp);
    m = fread(&pair_n, sizeof(uint64), 1, fp);
    if (m != 1) bin_fread_error();
    
    MYMALLOC(buffer, YAHSBIN_ENTRY_SIZE * FREAD_BUFF_SIZE);
    if (buffer == NULL)
        mem_alloc_error("buffer array");

    max_n = 0x4000000ULL; // 64MB - ~ 1GB mem for hic_t array
    if (max_n > pair_n)
        max_n = pair_n;
    MYMALLOC(hics, max_n);
    if (hics == NULL)
        mem_alloc_error("hic array");

    pair_c = 0;
    n_recs = 0;
    n = 0;
    while (pair_c < pair_n) {
        m = fread(buffer, sizeof(uint8) * YAHSBIN_ENTRY_SIZE, FREAD_BUFF_SIZE, fp);
        for (bbuff = buffer, ebuff = buffer + m * YAHSBIN_ENTRY_SIZE; 
            bbuff < ebuff && pair_c++ < pair_n; bbuff += YAHSBIN_ENTRY_SIZE) {
            if (pair_c % REPORT_EPOC == 0)
                fprintf(stderr, "[M::%s] %llu million records processed\n", __func__, pair_c / 1000000);
            
            if (bbuff[YAHSBIN_ENTRY_SIZE-1] < min_qual)
                continue;
            memcpy(&i0, bbuff,  sizeof(uint32));
            memcpy(&s0, bbuff + sizeof(uint32)*1, sizeof(uint32));
            memcpy(&i1, bbuff + sizeof(uint32)*2, sizeof(uint32));
            memcpy(&s1, bbuff + sizeof(uint32)*3, sizeof(uint32));
            
            if (sd_coordinate_conversion(dicts, i0, s0, &q0, &p0, 1) != CC_SUCCESS ||
                sd_coordinate_conversion(dicts, i1, s1, &q1, &p1, 1) != CC_SUCCESS)
                continue;

            p0 /= bin_size;
            p1 /= bin_size;

            if (q0 > q1) {
                SWAP(uint32, q0, q1);
                SWAP(uint64, p0, p1);
            }
            if (q0 == q1 && p0 > p1)
                SWAP(uint64, p0, p1);

            hics[n++] = (hic_t){q0, p0, q1, p1, 1};
            if (n == max_n) {
                n = hic_sort_merge(hics, n);
                if (n > (max_n>>1)) {
                    max_n <<= 1;
                    MYREALLOC(hics, max_n);
                    if (hics == NULL)
                        mem_alloc_error("hic array");
                }
            }
            n_recs++;
        }
    }
    
    n = hic_sort_merge(hics, n);
    
    fprintf(stderr, "[M::%s] processed %llu read pairs\n", __func__, pair_c);
    fprintf(stderr, "[M::%s] retained %llu records\n", __func__, n_recs);
    fprintf(stderr, "[M::%s] linked %llu sequence regions\n", __func__, n);
    
    free(buffer);
    fclose(fp);

    if (n == 0) {
        free(hics);
        return NULL;
    }
    MYREALLOC(hics, n);

    if (_nhic) *_nhic = n;

    return hics;
}

hic_t *read_hic_from_binary_sd_conversion(char *hic_bfile, asm_dict_t *dicts, int bin_size, uint8 min_qual, int64 *_nhic)
{
    if (_nhic) *_nhic = 0;
    if (hic_bfile == NULL)
        return NULL;

    FILE *fp;
    int64 magic_number;
    
    fp = fopen(hic_bfile, "r");
    if (fp == NULL)
        return NULL;
    if (fread(&magic_number, sizeof(int64), 1, fp) != 1)
        bin_fread_error();
    fclose(fp);

    if (is_valid_alnh_bin_header(magic_number))
        return read_hic_from_binary_sd_conversion_alnhbin(hic_bfile, dicts, bin_size, min_qual, _nhic);
    else if (is_valid_yahs_bin_header(magic_number))
        return read_hic_from_binary_sd_conversion_yahsbin(hic_bfile, dicts, bin_size, min_qual, _nhic);
    else
        fprintf(stderr, "[E::%s] not a valid BIN file\n", __func__);
    
    return NULL;
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