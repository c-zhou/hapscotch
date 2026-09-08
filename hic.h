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
#ifndef HIC_H_
#define HIC_H_

#include <stdlib.h>
#include <stdint.h>

#include "sdict.h"

#define FREAD_BUFF_SIZE (1<<16)

#define ALNHBIN_H 0x414C4E4842494E56
#define ALNHBIN_V 0x1
#define YAHSBIN_H 0x5941485342494E56
#define YAHSBIN_V 0x2

#define ALNHBIN_ENTRY_SIZE 25
#define YAHSBIN_ENTRY_SIZE 17

#define REPORT_EPOC 100000000

typedef enum { NOSET, BED, BAM, BIN, PA5, ONE } fileType_t;

extern const char *fileTypeNames[6];

#define HIC_NORM_WINDOW 1000

typedef struct {
    int aseq, apos;
    int bseq, bpos;
    int nhic;
} hic_t;

#ifdef __cplusplus 
extern "C" {
#endif
char *write_binary_hic_data(char *f, fileType_t f_type, sdict_t *dicts, int read_len, char *out);
char *write_binary_hic_data_multi(char **files, int n_files, fileType_t f_type, sdict_t *dicts, int read_len, char *out);
hic_t *read_hic_from_binary(char *bf, sdict_t *dicts, int bin_size, uint8 min_qual, int64 *_nhic);
hic_t *read_hic_from_binary_sd_conversion(char *hic_bfile, asm_dict_t *dicts, int bin_size, uint8 min_qual, int64 *_nhic);
void write_binary_hic_data_pseudo_yahs(char *bf, sdict_t *dicts, asm_dict_t *break_dict, int32 **smap, char *out);
int match_binary_file_sdict(char *f, sdict_t *dict);
#ifdef __cplusplus
}
#endif

static void bin_fread_error()
{
    fprintf(stderr, "[E::%s] error reading binary file\n", __func__);
    exit(EXIT_FAILURE);
}

static void write_yahs_bin_header(FILE *fo)
{
    int64 magic_number = YAHSBIN_H;
    int64 bin_version = YAHSBIN_V;
    magic_number |= bin_version;
    fwrite(&magic_number, sizeof(int64), 1, fo);
}

static int is_valid_yahs_bin_header(int64 n)
{
    int64 magic_number = YAHSBIN_H;
    int64 bin_version = YAHSBIN_V;
    magic_number |= bin_version;
    return n == magic_number;
}

static void write_alnh_bin_header(FILE *fo)
{
    int64 magic_number = ALNHBIN_H;
    int64 bin_version = ALNHBIN_V;
    magic_number |= bin_version;
    fwrite(&magic_number, sizeof(int64), 1, fo);
}

static int is_valid_alnh_bin_header(int64 n)
{
    int64 magic_number = ALNHBIN_H;
    int64 bin_version = ALNHBIN_V;
    magic_number |= bin_version;
    return n == magic_number;
}

static int is_valid_bin_header(int64 n)
{
    return is_valid_alnh_bin_header(n) || is_valid_yahs_bin_header(n);
}

static uint64 write_binary_sequence_dictionary(FILE *fo, sdict_t *dict)
{
    uint32 i;
    uint64 l;
    // write number of sequences
    fwrite(&dict->n, sizeof(uint32), 1, fo);
    // write sequence lengths
    for (i = 0; i < dict->n; ++i)
        fwrite(&dict->s[i].len, sizeof(uint32), 1, fo);
    l = dict->n;
    for (i = 0; i < dict->n; ++i)
        l += strlen(dict->s[i].name);
    // write total length of sequence names
    fwrite(&l, sizeof(uint64), 1, fo);
    // write sequence name with a appending white space
    const unsigned char space = ' ';
    for (i = 0; i < dict->n; ++i) {
        fwrite(dict->s[i].name, 1, strlen(dict->s[i].name), fo);
        fwrite(&space, 1, 1, fo);
    }

    return sizeof(uint32) * (dict->n + 3) + l;
}

static void binary_fseek_skip_sdict(FILE *fp)
{
    uint32 n;
    uint64 l, m;
    m = fread(&n, sizeof(uint32), 1, fp);
    if (m != 1) bin_fread_error();
    fseek(fp, sizeof(uint32) * n, SEEK_CUR);
    m = fread(&l, sizeof(uint64), 1, fp);
    if (m != 1) bin_fread_error();
    fseek(fp, l, SEEK_CUR);
}

static inline int is_read_pair(const char *rname0, const char *rname1)
{
    int n = strlen(rname0);
    if (n == 0 || n != strlen(rname1))
          return 0;
    if (!strcmp(rname0, rname1) || 
            (n > 2 &&
             !strncmp(rname0, rname1, n - 1) &&
             rname0[n - 2] == '/' &&
             rname1[n - 2] == '/' &&
             rname0[n - 1] != rname1[n - 1]))
        return 1;
    return 0;
}


#endif /* HIC_H_ */

