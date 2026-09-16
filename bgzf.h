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

#ifndef BGZF_H_
#define BGZF_H_

#include <stdio.h>
#include <ctype.h>

// Max number of uncompressed bytes packed into a single BGZF block
#define BGZF_MAX_BLOCK_SIZE 0xff00u

// Streaming BGZF writer: buffers writes and flushes a BGZF block to the underlying
// FILE* every time BGZF_MAX_BLOCK_SIZE bytes have accumulated, so callers do not need
// to hold the whole (potentially huge) output in memory at once.
typedef struct bgzf_writer_s bgzf_writer_t;

#ifdef __cplusplus
extern "C" {
#endif

// Writes an in-memory buffer as a complete standalone bgzf file (one or more blocks + EOF marker).
// Returns 0 on success, -1 on error.
int write_bgzf_file(const char *filename, const unsigned char *data, size_t total_length);

// Wraps 'out' (a plain FILE*, e.g. an opened file or stdout) for streaming bgzf writes.
// 'out' is not owned/closed by the writer. Returns NULL on allocation failure.
bgzf_writer_t *bgzf_writer_open(FILE *out);

// Appends 'len' bytes to the writer, transparently flushing full blocks as needed.
// Returns 0 on success, -1 on error (writer is left unusable after an error).
int bgzf_writer_write(bgzf_writer_t *bw, const void *data, size_t len);

// Flushes any buffered bytes as a final block, appends the standard BGZF EOF marker,
// and frees 'bw'. Does NOT close the underlying FILE* passed to bgzf_writer_open().
// Returns 0 on success, -1 on error.
int bgzf_writer_close(bgzf_writer_t *bw);

// Wraps 'out' in a stdio FILE* that transparently bgzip-compresses everything written
// to it (via fprintf/fputc/fwrite/...) and forwards the compressed bytes to 'out'.
// Closing the returned FILE* (fclose) finalizes the bgzf stream (final block + EOF
// marker) but leaves 'out' itself open - the caller remains responsible for it.
// Returns NULL if the platform has no supported custom-stream backend, or on error.
FILE *bgzf_fopen_write(FILE *out);

// Check if a filename ends with ".[g|G][z|Z]". Returns 1 if true, 0 otherwise.
int endsWithDotGz(const char *f);

#ifdef __cplusplus
}
#endif

#endif /* BGZF_H_ */
