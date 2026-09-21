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

// fopencookie() (glibc) needs _GNU_SOURCE defined before any system header is pulled in;
// BSD-derived libcs (macOS, *BSD) instead provide funopen() with no feature macro needed.
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__DragonFly__)
#define BGZF_USE_FUNOPEN 1
#else
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define BGZF_USE_FOPENCOOKIE 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>

#include "misc.h"
#include "bgzf.h"

// --------------------------------------------------------------------------
// Streaming writer: accumulates data into a BGZF_MAX_BLOCK_SIZE buffer and
// flushes each full block to the underlying FILE* as soon as it fills up.
// --------------------------------------------------------------------------

struct bgzf_writer_s {
    z_stream strm; // deflate stream
    FILE *out;     // not owned
    size_t len;    // bytes in buffer
    int error;
    unsigned char buf[BGZF_MAX_BLOCK_SIZE];
};

// 28-byte BGZF EOF marker block
static const unsigned char BGZF_EOF_BLOCK[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};

// BGZF block containing up to BGZF_MAX_BLOCK_SIZE bytes of raw data
int write_bgzf_block(bgzf_writer_t *bw) 
{
    unsigned char *data = bw->buf;
    size_t length = bw->len;
    FILE *out = bw->out;
    z_stream *strm = &bw->strm;

    if (length > BGZF_MAX_BLOCK_SIZE) return -1;

    unsigned char comp_buf[65536];
    
    // Initialize zlib raw deflate stream (-15 suppresses zlib/gzip headers)
    memset(strm, 0, sizeof(z_stream));
    if (deflateInit2(strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        return -1;
    }

    strm->next_in = (Bytef *)data;
    strm->avail_in = (uInt)length;
    strm->next_out = comp_buf;
    strm->avail_out = sizeof(comp_buf);

    // Compress the block; a single Z_FINISH call must fully drain avail_in given
    // comp_buf is sized well above deflateBound(BGZF_MAX_BLOCK_SIZE) - if it doesn't,
    // the output would be silently truncated, so treat that as a hard error.
    int ret = deflate(strm, Z_FINISH);
    size_t comp_len = strm->total_out;
    int ok = (ret == Z_STREAM_END && strm->avail_in == 0);
    deflateEnd(strm);
    if (!ok) return -1;

    // Calculate total block size BSIZE = Header (18) + Compressed Body + Trailer (8) - 1
    uint16 bsize = (uint16)(18 + comp_len + 8 - 1);
    uint32 crc = crc32(0L, data, (uInt)length);
    uint32 isize = (uint32)length;

    // Construct 18-byte Header with the required BGZF extra field
    unsigned char header[18] = {
        0x1f, 0x8b,                   // GZIP magic number
        0x08,                         // Compression method (DEFLATE)
        0x04,                         // Flags (FEXTRA set)
        0x00, 0x00, 0x00, 0x00,       // MTIME (0)
        0x00,                         // XFL
        0xff,                         // OS (unknown)
        0x06, 0x00,                   // XLEN = 6 bytes
        'B', 'C',                     // Subfield ID: "BC"
        0x02, 0x00,                   // Subfield length: 2 bytes
        (unsigned char)(bsize & 0xFF),// BSIZE low byte
        (unsigned char)(bsize >> 8)   // BSIZE high byte
    };

    // Construct 8-byte Trailer
    unsigned char trailer[8] = {
        (unsigned char)(crc & 0xFF), (unsigned char)((crc >> 8) & 0xFF),
        (unsigned char)((crc >> 16) & 0xFF), (unsigned char)((crc >> 24) & 0xFF),
        (unsigned char)(isize & 0xFF), (unsigned char)((isize >> 8) & 0xFF),
        (unsigned char)((isize >> 16) & 0xFF), (unsigned char)((isize >> 24) & 0xFF)
    };

    // Write header, compressed body, and trailer to file
    if (fwrite(header, 1, 18, out) != 18) return -1;
    if (fwrite(comp_buf, 1, comp_len, out) != comp_len) return -1;
    if (fwrite(trailer, 1, 8, out) != 8) return -1;

    return 0;
}

// Writes buffer as a sequence of BGZF blocks + EOF block
int write_bgzf_file(const char *filename, const unsigned char *data, size_t total_length) 
{
    FILE *out = fopen(filename, "wb");
    if (!out) return -1;

    bgzf_writer_t *bw = bgzf_writer_open(out);
    if (!bw) { fclose(out); return -1; }

    int ret = bgzf_writer_write(bw, data, total_length);
    if (bgzf_writer_close(bw) != 0) ret = -1;

    fclose(out);
    return ret;
}

bgzf_writer_t *bgzf_writer_open(FILE *out) 
{
    if (!out) return NULL;
    bgzf_writer_t *bw;
    MYMALLOC(bw, 1);
    if (!bw) return NULL;
    memset(&bw->strm, 0, sizeof(z_stream));
    bw->out = out;
    bw->len = 0;
    bw->error = 0;
    return bw;
}

int bgzf_writer_write(bgzf_writer_t *bw, const void *data, size_t len) 
{
    if (!bw || bw->error) return -1;

    const unsigned char *p = (const unsigned char *)data;
    while (len > 0) {
        size_t space = BGZF_MAX_BLOCK_SIZE - bw->len;
        size_t take = len < space ? len : space;
        memcpy(bw->buf + bw->len, p, take);
        bw->len += take;
        p += take;
        len -= take;
        if (bw->len == BGZF_MAX_BLOCK_SIZE) {
            if (write_bgzf_block(bw) != 0) {
                bw->error = 1;
                return -1;
            }
            bw->len = 0;
        }
    }
    return 0;
}

int bgzf_writer_close(bgzf_writer_t *bw) 
{
    if (!bw) return -1;

    int ret = bw->error ? -1 : 0;
    if (ret == 0 && bw->len > 0 && write_bgzf_block(bw) != 0)
        ret = -1;
    if (ret == 0 && fwrite(BGZF_EOF_BLOCK, 1, sizeof(BGZF_EOF_BLOCK), bw->out) != sizeof(BGZF_EOF_BLOCK))
        ret = -1;

    free(bw);
    return ret;
}

int endsWithDotGz(const char *f) 
{
    if (!f) return 0;
    size_t len = strlen(f);
    if (len < 3) return 0;
    f += len - 3;
    return (f[0] == '.' && 
        tolower((unsigned char) f[1]) == 'g' && 
        tolower((unsigned char) f[2]) == 'z');
}

// --------------------------------------------------------------------------
// FILE* wrapper: lets ordinary stdio calls (fprintf/fputc/fwrite/...) write
// transparently into a bgzf_writer_t, backed by a platform custom-stream API.
// --------------------------------------------------------------------------

#if defined(BGZF_USE_FUNOPEN)

static int bgzf_stream_write(void *cookie, const char *buf, int nbytes) 
{
    if (nbytes < 0) return -1;
    if (bgzf_writer_write((bgzf_writer_t *)cookie, buf, (size_t)nbytes) != 0) return -1;
    return nbytes;
}

static int bgzf_stream_close(void *cookie) 
{
    return bgzf_writer_close((bgzf_writer_t *)cookie);
}

FILE *bgzf_fopen_write(FILE *out) 
{
    bgzf_writer_t *bw = bgzf_writer_open(out);
    if (!bw) return NULL;

    FILE *fp = funopen(bw, NULL, bgzf_stream_write, NULL, bgzf_stream_close);
    if (!fp) { free(bw); return NULL; }

    // bypass stdio buffering
    setvbuf(fp, NULL, _IONBF, 0);

    return fp;
}

#elif defined(BGZF_USE_FOPENCOOKIE)

static ssize_t bgzf_stream_write(void *cookie, const char *buf, size_t size) 
{
    if (bgzf_writer_write((bgzf_writer_t *)cookie, buf, size) != 0) return -1;
    return (ssize_t)size;
}

static int bgzf_stream_close(void *cookie) 
{
    return bgzf_writer_close((bgzf_writer_t *)cookie);
}

FILE *bgzf_fopen_write(FILE *out) 
{
    bgzf_writer_t *bw = bgzf_writer_open(out);
    if (!bw) return NULL;

    cookie_io_functions_t io = { 
        .read = NULL, 
        .write = bgzf_stream_write, 
        .seek = NULL, 
        .close = bgzf_stream_close
    };
    FILE *fp = fopencookie(bw, "w", io);
    if (!fp) { free(bw); return NULL; }

    // bypass stdio buffering
    setvbuf(fp, NULL, _IONBF, 0);

    return fp;
}

#else

FILE *bgzf_fopen_write(FILE *out) 
{
    (void)out;
    return NULL;
}

#endif

#ifdef TEST_BZGF_WRITER
#include <errno.h>

int main(int argc, char *argv[]) 
{
    const char *seq1 = ">seq1\nACCTGGCAGT\nAACCG";
    const char *seq2 = ">seq2\nTTGACCAGTG";

    int bgz = 0;
    if (argc > 1) {
        if (strcmp(argv[1], "-z") == 0)
            bgz = 1;
        else if (strcmp(argv[1], "-") != 0) {
            // redirect output to a file specified by argv[1]
            if (freopen(argv[1], "wb", stdout) == NULL) {
                fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", argv[1], strerror(errno));
                return 1;
            }
            if (endsWithDotGz(argv[1])) bgz = 1;
        }
    }

    FILE *fo = bgz? bgzf_fopen_write(stdout) : stdout;
    if (fo) {
        fprintf(fo, "%s\n", seq1);
        fprintf(fo, "%s\n", seq2);
        fclose(fo);

        return 0;
    }

    return 1;
}
#endif