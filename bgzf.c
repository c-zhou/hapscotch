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
#include <pthread.h>
#include <zlib.h>

#include "misc.h"
#include "bgzf.h"

// Size of the buffer a compressed block is written into. Must stay comfortably above
// deflateBound(BGZF_MAX_BLOCK_SIZE) so a single Z_FINISH call is guaranteed to drain
// avail_in in one pass (see bgzf_deflate_block()) - 64 KiB input can never deflate to
// more than ~64 KiB + a small fixed overhead, so 65536 is a safe, simple upper bound.
#define BGZF_COMP_BUF_SIZE 65536u

// --------------------------------------------------------------------------
// Shared block-level helpers: compressing one block's worth of raw bytes,
// and writing one already-compressed block (header + body + trailer) to a
// FILE*. Used by both the synchronous writer and the multithreaded pipeline
// below, so the two paths produce byte-for-byte identical output.
// --------------------------------------------------------------------------

// 28-byte BGZF EOF marker block
static const unsigned char BGZF_EOF_BLOCK[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};

// Deflates up to BGZF_MAX_BLOCK_SIZE bytes of 'data' into 'comp' using an existing raw
// (headerless) deflate stream that the caller has already initialised (deflateInit2 with
// -15 windowBits, as BGZF requires) and owns the lifecycle of. 'strm' is reset before use
// via deflateReset(), which clears any per-block history/dictionary while reusing the
// stream's already-allocated internal tables - so output is identical to a fresh stream
// per block, without paying that allocation cost on every call. comp_cap must be at least
// BGZF_COMP_BUF_SIZE. On success writes the compressed length to *comp_len.
// Returns 0 on success, -1 on error.
static int bgzf_deflate_block(z_stream *strm, const unsigned char *data, size_t length,
                               unsigned char *comp, size_t comp_cap, size_t *comp_len)
{
    if (length > BGZF_MAX_BLOCK_SIZE || comp_cap < BGZF_COMP_BUF_SIZE) return -1;

    if (deflateReset(strm) != Z_OK) return -1;

    strm->next_in = (Bytef *)data;
    strm->avail_in = (uInt)length;
    strm->next_out = comp;
    strm->avail_out = (uInt)comp_cap;

    // A single Z_FINISH call must fully drain avail_in given comp_cap is sized well
    // above deflateBound(BGZF_MAX_BLOCK_SIZE) - if it doesn't, treat that as a hard error
    // rather than risk silently truncated output.
    int ret = deflate(strm, Z_FINISH);
    size_t out_len = strm->total_out;
    int ok = (ret == Z_STREAM_END && strm->avail_in == 0);
    if (!ok) return -1;

    *comp_len = out_len;
    return 0;
}

// Writes one already-compressed BGZF block (18-byte header + compressed body + 8-byte
// trailer) to 'out'. 'crc' and 'isize' are the CRC32 and length of the *uncompressed*
// data. Returns 0 on success, -1 on error.
static int write_bgzf_block_raw(FILE *out, const unsigned char *comp, size_t comp_len,
                                 uint32 crc, uint32 isize)
{
    // Calculate total block size BSIZE = Header (18) + Compressed Body + Trailer (8) - 1
    uint16 bsize = (uint16)(18 + comp_len + 8 - 1);

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
    if (fwrite(comp, 1, comp_len, out) != comp_len) return -1;
    if (fwrite(trailer, 1, 8, out) != 8) return -1;

    return 0;
}

// --------------------------------------------------------------------------
// Multithreaded pipeline internals (used when bw->mt is set).
//
// Design: a fixed ring of slots, each holding one block's raw bytes and its
// compressed output. The caller's thread (the "producer") fills slots and
// submits full ones; a pool of worker threads picks up FILLED slots in any
// order and compresses them; one dedicated writer thread drains slots in
// strict submission order and writes them to 'out'. Because BGZF blocks are
// compressed independently (no cross-block dictionary), compression can run
// out of order - only the final file order matters, which the single writer
// thread guarantees.
//
// Slot state machine:  EMPTY -> FILLED -> COMPRESSING -> COMPRESSED -> EMPTY
//   EMPTY -> FILLED       : producer, after finishing a block
//   FILLED -> COMPRESSING : any worker, after claiming the slot
//   COMPRESSING -> COMPRESSED : that worker, after compressing it
//   COMPRESSED -> EMPTY   : the writer thread, after flushing it to disk
// --------------------------------------------------------------------------

typedef enum {
    SLOT_EMPTY = 0,
    SLOT_FILLED,
    SLOT_COMPRESSING,
    SLOT_COMPRESSED
} bgzf_slot_state_t;

typedef struct {
    unsigned char raw[BGZF_MAX_BLOCK_SIZE];
    size_t        raw_len;
    unsigned char comp[BGZF_COMP_BUF_SIZE];
    size_t        comp_len;
    uint32        crc;
    uint32        isize;
    bgzf_slot_state_t state;
} bgzf_slot_t;

struct bgzf_writer_s {
    FILE *out;   // not owned
    int   mt;    // 0 = synchronous single-threaded path, 1 = multithreaded pipeline
    int   error; // sticky error flag for the synchronous path

    // ---- synchronous path (mt == 0) ----
    z_stream strm;  // persistent deflate stream, reused across blocks via deflateReset()
    size_t len;
    unsigned char buf[BGZF_MAX_BLOCK_SIZE];

    // ---- multithreaded path (mt == 1) ----
    int n_threads;
    size_t n_slots;
    bgzf_slot_t *slots;
    unsigned char *slot_needs_check; // producer-private (never touched by other threads):
                                      // n_slots flags; 1 means "this slot was submitted
                                      // before and must be confirmed EMPTY before reuse"
    pthread_t *workers;
    pthread_t  writer_tid;
    int        writer_started;

    size_t *ready_queue;  // circular FIFO of slot indices whose state is SLOT_FILLED;
                           // capacity n_slots (that many can never all be ready at once,
                           // but sizing it that way keeps the arithmetic simple/safe)
    size_t  ready_head;
    size_t  ready_count;

    pthread_mutex_t mu;
    pthread_cond_t  cv_producer; // signalled when a slot becomes EMPTY
    pthread_cond_t  cv_worker;   // signalled (one waiter) when a slot is pushed to ready_queue;
                                  // broadcast only on error / shutdown, when every waiter must wake
    pthread_cond_t  cv_writer;   // signalled when a slot becomes COMPRESSED, or on shutdown/error

    size_t fill_seq;        // producer-private: sequence number of the slot being filled
    size_t total_to_write;  // set at close(): final number of blocks that will be submitted
    int    closing;         // producer is done submitting (set during close())
    int    closing_workers; // workers may exit once no FILLED slot remains (set during close())
    int    mt_error;        // sticky error flag shared across threads; only touched under mu
};

static void *bgzf_worker_main(void *arg);
static void *bgzf_writer_thread_main(void *arg);

// Tears down whatever subset of the mt resources was actually initialised, based on
// 'stage' (how far bgzf_writer_open_mt() got before failing), and frees bw. Used only
// on the open_mt() failure path.
static void bgzf_writer_mt_teardown_partial(bgzf_writer_t *bw, int stage)
{
    if (stage >= 4) pthread_cond_destroy(&bw->cv_writer);
    if (stage >= 3) pthread_cond_destroy(&bw->cv_worker);
    if (stage >= 2) pthread_cond_destroy(&bw->cv_producer);
    if (stage >= 1) pthread_mutex_destroy(&bw->mu);
    free(bw->slots);
    free(bw->slot_needs_check);
    free(bw->workers);
    free(bw->ready_queue);
    free(bw);
}

// Pushes 'idx' onto the ready queue and wakes exactly one worker - a push always makes
// exactly one additional slot available, so waking more than one waiter would just send
// the extras straight back to sleep after contending for the mutex for nothing. Must be
// called with bw->mu held.
static void bgzf_ready_push(bgzf_writer_t *bw, size_t idx)
{
    bw->ready_queue[(bw->ready_head + bw->ready_count) % bw->n_slots] = idx;
    bw->ready_count++;
    pthread_cond_signal(&bw->cv_worker);
}

bgzf_writer_t *bgzf_writer_open_mt(FILE *out, int n_threads)
{
    if (!out) return NULL;
    if (n_threads < 1) n_threads = 1;
    if (n_threads == 1) return bgzf_writer_open(out); // no benefit to the extra machinery

    bgzf_writer_t *bw;
    MYMALLOC(bw, 1);
    if (!bw) return NULL;
    memset(bw, 0, sizeof(*bw));
    bw->out = out;
    bw->mt = 1;
    bw->n_threads = n_threads;
    bw->n_slots = (size_t)n_threads * 3;

    bw->slots = calloc(bw->n_slots, sizeof(bgzf_slot_t));
    bw->slot_needs_check = calloc(bw->n_slots, 1);
    bw->workers = calloc((size_t)n_threads, sizeof(pthread_t));
    bw->ready_queue = calloc(bw->n_slots, sizeof(size_t));
    if (!bw->slots || !bw->slot_needs_check || !bw->workers || !bw->ready_queue) {
        bgzf_writer_mt_teardown_partial(bw, 0);
        return NULL;
    }

    if (pthread_mutex_init(&bw->mu, NULL) != 0) { bgzf_writer_mt_teardown_partial(bw, 0); return NULL; }
    if (pthread_cond_init(&bw->cv_producer, NULL) != 0) { bgzf_writer_mt_teardown_partial(bw, 1); return NULL; }
    if (pthread_cond_init(&bw->cv_worker, NULL) != 0) { bgzf_writer_mt_teardown_partial(bw, 2); return NULL; }
    if (pthread_cond_init(&bw->cv_writer, NULL) != 0) { bgzf_writer_mt_teardown_partial(bw, 3); return NULL; }

    if (pthread_create(&bw->writer_tid, NULL, bgzf_writer_thread_main, bw) != 0) {
        bgzf_writer_mt_teardown_partial(bw, 4);
        return NULL;
    }
    bw->writer_started = 1;

    int i, started = 0;
    for (i = 0; i < n_threads; i++) {
        if (pthread_create(&bw->workers[i], NULL, bgzf_worker_main, bw) != 0) break;
        started++;
    }
    if (started < n_threads) {
        // Couldn't start every worker thread. The ones that did start can still make
        // progress, so don't fail the whole writer - just shrink the pool and flag the
        // shortfall as a hard error handled the same way any other mt_error is: close()
        // will still cleanly join everything and report failure.
        bw->n_threads = started;
        if (started == 0) {
            pthread_mutex_lock(&bw->mu);
            bw->mt_error = 1;
            bw->closing = 1;
            bw->closing_workers = 1;
            bw->total_to_write = 0;
            pthread_cond_broadcast(&bw->cv_writer);
            pthread_mutex_unlock(&bw->mu);
        }
    }

    return bw;
}

// Multithreaded write path: the producer (caller's thread) copies data directly into
// the raw buffer of the slot currently being filled, avoiding an extra intermediate
// copy. Synchronisation only happens at slot boundaries (once per BGZF_MAX_BLOCK_SIZE
// bytes), not per byte.
static int bgzf_writer_write_mt(bgzf_writer_t *bw, const void *data, size_t len)
{
    if (!bw) return -1;

    const unsigned char *p = (const unsigned char *)data;
    while (len > 0) {
        size_t idx = bw->fill_seq % bw->n_slots;
        bgzf_slot_t *s = &bw->slots[idx];

        if (bw->slot_needs_check[idx]) {
            pthread_mutex_lock(&bw->mu);
            while (s->state != SLOT_EMPTY && !bw->mt_error)
                pthread_cond_wait(&bw->cv_producer, &bw->mu);
            int err = bw->mt_error;
            pthread_mutex_unlock(&bw->mu);
            if (err) return -1;
            bw->slot_needs_check[idx] = 0;
        }

        size_t space = BGZF_MAX_BLOCK_SIZE - s->raw_len;
        size_t take = len < space ? len : space;
        memcpy(s->raw + s->raw_len, p, take);
        s->raw_len += take;
        p += take;
        len -= take;

        if (s->raw_len == BGZF_MAX_BLOCK_SIZE) {
            pthread_mutex_lock(&bw->mu);
            s->state = SLOT_FILLED;
            bgzf_ready_push(bw, idx);
            pthread_mutex_unlock(&bw->mu);
            bw->slot_needs_check[idx] = 1; // must wait for EMPTY before this slot is reused
            bw->fill_seq++;
        }
    }
    return 0;
}

static void *bgzf_worker_main(void *arg)
{
    bgzf_writer_t *bw = (bgzf_writer_t *)arg;

    // One persistent deflate stream per worker thread, reused (via deflateReset() inside
    // bgzf_deflate_block()) for every block this thread ever compresses over its lifetime.
    // Avoids paying deflateInit2's ~256KB internal allocation on every single 64KB block -
    // with many worker threads doing that concurrently, the allocator itself becomes a
    // point of contention that gets worse as thread count grows, independent of actual
    // CPU compression throughput.
    z_stream strm;
    memset(&strm, 0, sizeof(strm));
    if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        pthread_mutex_lock(&bw->mu);
        bw->mt_error = 1;
        pthread_cond_broadcast(&bw->cv_worker);
        pthread_cond_broadcast(&bw->cv_writer);
        pthread_mutex_unlock(&bw->mu);
        return NULL;
    }

    for (;;) {
        bgzf_slot_t *s = NULL;

        pthread_mutex_lock(&bw->mu);
        for (;;) {
            if (bw->mt_error) { pthread_mutex_unlock(&bw->mu); deflateEnd(&strm); return NULL; }
            if (bw->ready_count > 0) {
                size_t idx = bw->ready_queue[bw->ready_head];
                bw->ready_head = (bw->ready_head + 1) % bw->n_slots;
                bw->ready_count--;
                s = &bw->slots[idx];
                s->state = SLOT_COMPRESSING;
                break;
            }
            if (bw->closing_workers) { pthread_mutex_unlock(&bw->mu); deflateEnd(&strm); return NULL; }
            pthread_cond_wait(&bw->cv_worker, &bw->mu);
        }
        pthread_mutex_unlock(&bw->mu);

        size_t comp_len = 0;
        int ok = (bgzf_deflate_block(&strm, s->raw, s->raw_len, s->comp, sizeof(s->comp), &comp_len) == 0);
        uint32 crc = 0, isize = 0;
        if (ok) {
            crc = crc32(0L, s->raw, (uInt)s->raw_len);
            isize = (uint32)s->raw_len;
        }

        pthread_mutex_lock(&bw->mu);
        if (ok) {
            s->comp_len = comp_len;
            s->crc = crc;
            s->isize = isize;
            s->state = SLOT_COMPRESSED;
            // No new ready-queue entry was created by finishing a block, so there is
            // nothing for a peer worker to do - only the writer thread needs waking.
        } else {
            bw->mt_error = 1;
            pthread_cond_broadcast(&bw->cv_worker); // every waiter must notice the error and bail
        }
        pthread_cond_broadcast(&bw->cv_writer);
        pthread_mutex_unlock(&bw->mu);
    }
}

static void *bgzf_writer_thread_main(void *arg)
{
    bgzf_writer_t *bw = (bgzf_writer_t *)arg;
    size_t seq = 0;

    for (;;) {
        size_t idx = seq % bw->n_slots;
        bgzf_slot_t *s = &bw->slots[idx];

        pthread_mutex_lock(&bw->mu);
        while (s->state != SLOT_COMPRESSED && !bw->mt_error) {
            if (bw->closing && seq >= bw->total_to_write) {
                pthread_mutex_unlock(&bw->mu);
                return NULL; // every expected block has been written
            }
            pthread_cond_wait(&bw->cv_writer, &bw->mu);
        }
        if (bw->mt_error) { pthread_mutex_unlock(&bw->mu); return NULL; }
        pthread_mutex_unlock(&bw->mu);

        if (write_bgzf_block_raw(bw->out, s->comp, s->comp_len, s->crc, s->isize) != 0) {
            pthread_mutex_lock(&bw->mu);
            bw->mt_error = 1;
            pthread_cond_broadcast(&bw->cv_producer);
            pthread_cond_broadcast(&bw->cv_worker);
            pthread_mutex_unlock(&bw->mu);
            return NULL;
        }

        pthread_mutex_lock(&bw->mu);
        s->raw_len = 0;
        s->state = SLOT_EMPTY;
        pthread_cond_broadcast(&bw->cv_producer);
        pthread_mutex_unlock(&bw->mu);

        seq++;
    }
}

static int bgzf_writer_close_mt(bgzf_writer_t *bw)
{
    int i, ret = -1;

    // Submit any partial trailing block still sitting in the current slot.
    if (!bw->mt_error) {
        size_t idx = bw->fill_seq % bw->n_slots;
        bgzf_slot_t *s = &bw->slots[idx];
        if (s->raw_len > 0) {
            pthread_mutex_lock(&bw->mu);
            s->state = SLOT_FILLED;
            bgzf_ready_push(bw, idx);
            pthread_mutex_unlock(&bw->mu);
            bw->fill_seq++;
        }
    }

    pthread_mutex_lock(&bw->mu);
    bw->total_to_write = bw->fill_seq;
    bw->closing = 1;
    pthread_cond_broadcast(&bw->cv_writer);
    pthread_mutex_unlock(&bw->mu);

    // Join the writer thread first: it exits once every submitted block has been
    // flushed to disk in order (or immediately, if mt_error is already set).
    if (bw->writer_started) pthread_join(bw->writer_tid, NULL);

    pthread_mutex_lock(&bw->mu);
    bw->closing_workers = 1;
    pthread_cond_broadcast(&bw->cv_worker);
    pthread_mutex_unlock(&bw->mu);
    for (i = 0; i < bw->n_threads; i++)
        pthread_join(bw->workers[i], NULL);

    ret = bw->mt_error ? -1 : 0;
    if (ret == 0 && fwrite(BGZF_EOF_BLOCK, 1, sizeof(BGZF_EOF_BLOCK), bw->out) != sizeof(BGZF_EOF_BLOCK))
        ret = -1;

    pthread_cond_destroy(&bw->cv_writer);
    pthread_cond_destroy(&bw->cv_worker);
    pthread_cond_destroy(&bw->cv_producer);
    pthread_mutex_destroy(&bw->mu);
    free(bw->slots);
    free(bw->slot_needs_check);
    free(bw->workers);
    free(bw);

    return ret;
}

// --------------------------------------------------------------------------
// Synchronous (single-threaded) writer: unchanged behaviour from the original
// implementation, now built on the shared bgzf_deflate_block()/write_bgzf_block_raw()
// helpers above.
// --------------------------------------------------------------------------

static int write_bgzf_block(bgzf_writer_t *bw)
{
    if (bw->len > BGZF_MAX_BLOCK_SIZE) return -1;

    unsigned char comp_buf[BGZF_COMP_BUF_SIZE];
    size_t comp_len;
    if (bgzf_deflate_block(&bw->strm, bw->buf, bw->len, comp_buf, sizeof(comp_buf), &comp_len) != 0)
        return -1;

    uint32 crc = crc32(0L, bw->buf, (uInt)bw->len);
    uint32 isize = (uint32)bw->len;

    return write_bgzf_block_raw(bw->out, comp_buf, comp_len, crc, isize);
}

// Writes buffer as a sequence of BGZF blocks + EOF block
int write_bgzf_file(const char *filename, const unsigned char *data, size_t total_length) 
{
    return write_bgzf_file_mt(filename, data, total_length, 1);
}

int write_bgzf_file_mt(const char *filename, const unsigned char *data, size_t total_length, int n_threads)
{
    FILE *out = fopen(filename, "wb");
    if (!out) return -1;

    bgzf_writer_t *bw = bgzf_writer_open_mt(out, n_threads);
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
    memset(bw, 0, sizeof(*bw));
    bw->out = out;
    bw->mt = 0;
    bw->len = 0;
    bw->error = 0;
    // -15 suppresses zlib/gzip headers, giving raw DEFLATE data as BGZF requires. This
    // one stream is reused (via deflateReset()) for every block this writer ever flushes.
    if (deflateInit2(&bw->strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        free(bw);
        return NULL;
    }
    return bw;
}

static int bgzf_writer_write_st(bgzf_writer_t *bw, const void *data, size_t len)
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

int bgzf_writer_write(bgzf_writer_t *bw, const void *data, size_t len)
{
    if (!bw) return -1;
    return bw->mt ? bgzf_writer_write_mt(bw, data, len) : bgzf_writer_write_st(bw, data, len);
}

static int bgzf_writer_close_st(bgzf_writer_t *bw)
{
    int ret = bw->error ? -1 : 0;
    if (ret == 0 && bw->len > 0 && write_bgzf_block(bw) != 0)
        ret = -1;
    if (ret == 0 && fwrite(BGZF_EOF_BLOCK, 1, sizeof(BGZF_EOF_BLOCK), bw->out) != sizeof(BGZF_EOF_BLOCK))
        ret = -1;

    deflateEnd(&bw->strm);
    free(bw);
    return ret;
}

int bgzf_writer_close(bgzf_writer_t *bw)
{
    if (!bw) return -1;
    return bw->mt ? bgzf_writer_close_mt(bw) : bgzf_writer_close_st(bw);
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

static FILE *bgzf_fopen_write_impl(bgzf_writer_t *bw)
{
    if (!bw) return NULL;

    FILE *fp = funopen(bw, NULL, bgzf_stream_write, NULL, bgzf_stream_close);
    if (!fp) { bgzf_writer_close(bw); return NULL; }

    // bypass stdio buffering
    setvbuf(fp, NULL, _IONBF, 0);

    return fp;
}

FILE *bgzf_fopen_write(FILE *out)
{
    return bgzf_fopen_write_impl(bgzf_writer_open(out));
}

FILE *bgzf_fopen_write_mt(FILE *out, int n_threads)
{
    return bgzf_fopen_write_impl(bgzf_writer_open_mt(out, n_threads));
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

static FILE *bgzf_fopen_write_impl(bgzf_writer_t *bw)
{
    if (!bw) return NULL;

    cookie_io_functions_t io = { 
        .read = NULL, 
        .write = bgzf_stream_write, 
        .seek = NULL, 
        .close = bgzf_stream_close
    };
    FILE *fp = fopencookie(bw, "w", io);
    if (!fp) { bgzf_writer_close(bw); return NULL; }

    // bypass stdio buffering
    setvbuf(fp, NULL, _IONBF, 0);

    return fp;
}

FILE *bgzf_fopen_write(FILE *out)
{
    return bgzf_fopen_write_impl(bgzf_writer_open(out));
}

FILE *bgzf_fopen_write_mt(FILE *out, int n_threads)
{
    return bgzf_fopen_write_impl(bgzf_writer_open_mt(out, n_threads));
}

#else

FILE *bgzf_fopen_write(FILE *out) 
{
    (void)out;
    return NULL;
}

FILE *bgzf_fopen_write_mt(FILE *out, int n_threads)
{
    (void)out;
    (void)n_threads;
    return NULL;
}

#endif

#ifdef TEST_BZGF_WRITER
#include <errno.h>

int main(int argc, char *argv[]) 
{
    const char *seq1 = ">seq1\nACCTGGCAGT\nAACCG";
    const char *seq2 = ">seq2\nTTGACCAGTG";

    int i, bgz = 0, n_threads = 1;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-z") == 0) {
            bgz = 1;
        } else if (strncmp(argv[i], "-t", 2) == 0 && argv[i][2] != '\0') {
            n_threads = atoi(argv[i] + 2); // e.g. -t4 -> 4 worker threads
            if (n_threads > 1) bgz = 1;
        } else if (strcmp(argv[i], "-") != 0) {
            // redirect output to a file specified by argv[i]
            if (freopen(argv[i], "wb", stdout) == NULL) {
                fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", argv[i], strerror(errno));
                return 1;
            }
            if (endsWithDotGz(argv[i])) bgz = 1;
        }
    }

    FILE *fo = bgz ? bgzf_fopen_write_mt(stdout, n_threads) : stdout;
    if (fo) {
        fprintf(fo, "%s\n", seq1);
        fprintf(fo, "%s\n", seq2);
        fclose(fo);

        return 0;
    }

    return 1;
}
#endif
