#!/usr/bin/env python3
"""
hicmap.py — Visualise a hictools prepare contact map.

Usage:
    hicmap.py [options] [hicmap.txt]

The input file defaults to stdin when omitted or '-'.
"""

import argparse
import sys
import os
import fnmatch
import math
import numpy as np

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Plot a hictools contact map",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Normalisation methods (-m):
  clip     clip to <--clip-pct> percentile then scale to [0,1]
  sqrt     sqrt(count)               
  log1p    log(1 + count)
  none     raw counts

Colour schemes (-c):
  Any matplotlib colormap name is accepted, e.g.:
    YlOrRd  (default) yellow → orange → red
    Reds    white → red
    hot_r   white → yellow → red → black
    RdYlBu_r blue → yellow → red
    viridis purple → green → yellow
    magma   black → purple → orange → white
""",
    )
    ap.add_argument("input", nargs="?", default="-",
                    help="hictools prepare output file [stdin]")
    ap.add_argument("-m", "--norm", default="clip",
                    choices=["clip", "sqrt", "log1p", "none"],
                    help="normalisation method [clip]")
    ap.add_argument("--clip-pct", type=float, default=98.0, metavar="PCT",
                    help="percentile cap for --norm clip [98]")
    ap.add_argument("-c", "--cmap", default="viridis",
                    metavar="CMAP",
                    help="matplotlib colormap name [viridis]")
    ap.add_argument("--merge", type=int, default=None, metavar="BP",
                    help="sequence shorter than this will be merged (default: bin_size from file)")
    ap.add_argument("-x", "--x-seqs", metavar="NAME_OR_FILE", default=None,
                    help="sequence name or file listing names (one per line) "
                         "for the x-axis; also used for y-axis unless -y is given")
    ap.add_argument("-y", "--y-seqs", metavar="NAME_OR_FILE", default=None,
                    help="sequence name or file listing names (one per line) "
                         "for the y-axis (requires -x)")
    ap.add_argument("--labels", action="store_true",
                    help="draw sequence name labels and boundary lines on axes")
    ap.add_argument("--title", default=None,
                    help="figure title")
    ap.add_argument("--dpi", type=int, default=150,
                    help="output raster DPI [150]")
    ap.add_argument("--show", action="store_true",
                    help="display interactive window (default if no --png/--pdf)")
    ap.add_argument("--png", metavar="FILE", default=None,
                    help="write PNG output to FILE")
    ap.add_argument("--pdf", metavar="FILE", default=None,
                    help="write PDF output to FILE")
    args = ap.parse_args()
    if args.input == "-" and sys.stdin.isatty():
        ap.error("no input file given and stdin is a terminal")
    if args.y_seqs and not args.x_seqs:
        ap.error("-y/--y-seqs requires -x/--x-seqs")
    return args


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def read_hicmap(path):
    """
    Parse a hictools prepare output file.

    Returns
    -------
    seqs       : list of (name, length)   in file order
    bin_size   : int
    contacts   : list of (i1, b1, i2, b2, count)
    meta       : dict of ## key=value pairs
    """
    seqs = []
    bin_size = None
    contacts = []
    meta = {}

    fh = sys.stdin if path == "-" else open(path)
    try:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                # e.g.  ##key: value  or  ##key=value
                body = line[2:]
                for sep in (":", "="):
                    if sep in body:
                        k, v = body.split(sep, 1)
                        meta[k.strip()] = v.strip()
                        break
            elif line.startswith("S\t"):
                parts = line.split("\t")
                seqs.append((parts[1], int(parts[2])))
            elif line.startswith("B\t"):
                bin_size = int(line.split("\t")[1])
            elif line.startswith("C\t"):
                _, i1, b1, i2, b2, cnt = line.split("\t")
                contacts.append((int(i1), int(b1), int(i2), int(b2), int(cnt)))
    finally:
        if fh is not sys.stdin:
            fh.close()

    if bin_size is None:
        sys.exit("[E] no B (bin size) line found in input")
    if not seqs:
        sys.exit("[E] no S (sequence) lines found in input")

    return seqs, bin_size, contacts, meta


def read_seq_list(path):
    """Read sequence names from a file, one per line."""
    names = []
    with open(path) as f:
        for line in f:
            name = line.strip()
            if name:
                names.append(name)
    return names


def parse_seq_arg(value, all_seq_names=None):
    """Interpret a -x/-y value as a filename, wildcard pattern, or a single sequence name."""
    if os.path.isfile(value):
        return read_seq_list(value)
    # check for glob/wildcard characters
    if all_seq_names is not None and any(c in value for c in '*?[]'):
        matched = fnmatch.filter(all_seq_names, value)
        if not matched:
            print(f"[W] pattern '{value}' matched no sequences",
                  file=sys.stderr)
        return matched
    return [value]


def filter_seqs(seqs, names):
    """
    Select and reorder sequences according to *names*.

    Returns [(orig_idx, name, length), ...] where orig_idx is the
    sequence's index in the original (file-order) *seqs* list.
    """
    idx_map = {name: i for i, (name, _) in enumerate(seqs)}
    result = []
    for name in names:
        if name in idx_map:
            i = idx_map[name]
            result.append((i, name, seqs[i][1]))
        else:
            print(f"[W] sequence '{name}' not in input, skipping",
                  file=sys.stderr)
    return result


# ---------------------------------------------------------------------------
# Bin layout with small-sequence merging
# ---------------------------------------------------------------------------

def build_bin_layout(seqs, bin_size, merge_size=None):
    """
    Assign each sequence a run of output pixels (merged bins).

    Strategy
    --------
    * Large sequences (len >= bin_size) get ceil(len / bin_size) pixels
      normally; the last pixel may cover a shorter region (debris).
    * Small sequences (len < merge_size, defaults to bin_size) are
      grouped with their neighbours until the accumulated length reaches
      bin_size, then flushed as a single merged pixel.

    Returns
    -------
    n_pixels   : total number of output pixels
    seq_map    : list of (pixel_start, n_pixels_for_seq, orig_bins_per_seq)
                 indexed by sequence index.
                 seq_map[i] = (pix_offset, pix_count, [(raw_bin, pix_idx), ...])
                   The third element is a lookup: raw bin index on sequence i →
                   output pixel index (absolute).
    """
    if merge_size is None:
        merge_size = bin_size

    # First pass: assign every sequence its raw (pre-merge) bins
    # raw_bins[i] = list of (raw_bin_idx, seq_i) for each bin of seq i
    n_seqs = len(seqs)

    # ------------------------------------------------------------
    # We build a flat list of "atoms": each atom is either
    #   ('bin',  seq_idx, raw_bin, bp_covered)  for a single bin
    # We then merge runs of small-sequence atoms and reindex.
    # ------------------------------------------------------------
    atoms = []  # (seq_idx, raw_bin, bp_covered)
    for si, (name, length) in enumerate(seqs):
        n_raw = max(1, math.ceil(length / bin_size))
        for rb in range(n_raw):
            bp = min(bin_size, length - rb * bin_size)
            if bp <= 0:
                bp = length  # safety for length < bin_size
            atoms.append((si, rb, bp))

    # Merge pass: walk through atoms and group small-sequence atoms
    # until accumulated bp >= merge_size.
    #
    # "Small sequence" = the ENTIRE sequence is < merge_size bp.
    # Large-sequence atoms (even the short debris last-bin) keep their
    # own pixel because the sequence itself is large.
    small_seq = {si for si, (name, length) in enumerate(seqs) if length < merge_size}

    pixels = []        # each pixel: [(seq_idx, raw_bin), ...]
    pending = []       # accumulating small atoms
    pending_bp = 0

    for atom in atoms:
        si, rb, bp = atom
        if si in small_seq:
            pending.append((si, rb))
            pending_bp += bp
            if pending_bp >= merge_size:
                pixels.append(list(pending))
                pending = []
                pending_bp = 0
        else:
            # flush any pending small atoms first
            if pending:
                pixels.append(list(pending))
                pending = []
                pending_bp = 0
            pixels.append([(si, rb)])

    if pending:
        pixels.append(list(pending))

    # Build reverse lookup: (seq_idx, raw_bin) → pixel index
    lookup = {}   # (si, rb) → pix_idx
    for pix_idx, members in enumerate(pixels):
        for si, rb in members:
            lookup[(si, rb)] = pix_idx

    n_pixels = len(pixels)
    return n_pixels, lookup, pixels


def build_bin_layout_filtered(filtered_seqs, bin_size, merge_size=None):
    """
    Like build_bin_layout but for a filtered/reordered sequence list.

    Parameters
    ----------
    filtered_seqs : list of (orig_idx, name, length)

    Returns the same (n_pixels, lookup, pixels) tuple.
    Lookup keys are (orig_idx, raw_bin) → pixel index.
    """
    if merge_size is None:
        merge_size = bin_size

    atoms = []
    for orig_idx, name, length in filtered_seqs:
        n_raw = max(1, math.ceil(length / bin_size))
        for rb in range(n_raw):
            bp = min(bin_size, length - rb * bin_size)
            if bp <= 0:
                bp = length
            atoms.append((orig_idx, rb, bp))

    small_seq = {oi for oi, _, length in filtered_seqs if length < merge_size}

    pixels = []
    pending = []
    pending_bp = 0

    for atom in atoms:
        orig_idx, rb, bp = atom
        if orig_idx in small_seq:
            pending.append((orig_idx, rb))
            pending_bp += bp
            if pending_bp >= merge_size:
                pixels.append(list(pending))
                pending = []
                pending_bp = 0
        else:
            if pending:
                pixels.append(list(pending))
                pending = []
                pending_bp = 0
            pixels.append([(orig_idx, rb)])

    if pending:
        pixels.append(list(pending))

    lookup = {}
    for pix_idx, members in enumerate(pixels):
        for orig_idx, rb in members:
            lookup[(orig_idx, rb)] = pix_idx

    return len(pixels), lookup, pixels


# ---------------------------------------------------------------------------
# Contact matrix assembly
# ---------------------------------------------------------------------------

def assemble_matrix(contacts, lookup, n_pixels):
    """
    Build a dense float32 contact matrix of shape (n_pixels, n_pixels).

    Contacts outside lookup (e.g. bin index beyond sequence length due
    to rounding) are silently dropped.
    """
    mat = np.zeros((n_pixels, n_pixels), dtype=np.float32)
    dropped = 0
    for i1, b1, i2, b2, cnt in contacts:
        p1 = lookup.get((i1, b1))
        p2 = lookup.get((i2, b2))
        if p1 is None or p2 is None:
            dropped += 1
            continue
        mat[p1, p2] += cnt
        if p1 != p2:
            mat[p2, p1] += cnt   # keep matrix symmetric
    if dropped:
        print(f"[W] {dropped} contact entries had no matching pixel (dropped)",
              file=sys.stderr)
    return mat


def assemble_matrix_xy(contacts, x_lookup, y_lookup, n_x, n_y):
    """
    Build a (n_y × n_x) contact matrix with potentially different
    sequences on each axis.
    """
    mat = np.zeros((n_y, n_x), dtype=np.float32)
    dropped = 0
    for i1, b1, i2, b2, cnt in contacts:
        py1 = y_lookup.get((i1, b1))
        px2 = x_lookup.get((i2, b2))
        py2 = y_lookup.get((i2, b2))
        px1 = x_lookup.get((i1, b1))

        hit = False
        if py1 is not None and px2 is not None:
            mat[py1, px2] += cnt
            hit = True
        if py2 is not None and px1 is not None:
            # avoid double-counting when both map to the same cell
            if py2 != py1 or px1 != px2:
                mat[py2, px1] += cnt
            hit = True
        if not hit:
            dropped += 1
    if dropped:
        print(f"[W] {dropped} contact entries had no matching pixel (dropped)",
              file=sys.stderr)
    return mat


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

def normalise(mat, method, clip_pct=98.0):
    """Return a float32 array in [0, 1] ready for imshow."""
    m = mat.copy()
    if method == "log1p":
        m = np.log1p(m)
    elif method == "sqrt":
        m = np.sqrt(m)
    elif method == "clip":
        vmax = np.percentile(m[m > 0], clip_pct) if np.any(m > 0) else 1.0
        m = np.clip(m, 0, vmax)
    # "none" — use raw counts

    vmax = m.max()
    if vmax > 0:
        m = m / vmax
    return m.astype(np.float32)


# ---------------------------------------------------------------------------
# Tick / label helpers
# ---------------------------------------------------------------------------

def make_ticks(filtered_seqs, lookup, n_pixels, max_labels=40):
    """
    Return (tick_positions, tick_labels) for sequence boundaries.
    Only place a tick at the first pixel of each sequence.
    Thin out labels if there are too many sequences.

    filtered_seqs : list of (orig_idx, name, length)
    """
    positions = []
    labels = []
    seen = set()
    for orig_idx, name, _ in filtered_seqs:
        p = lookup.get((orig_idx, 0))
        if p is None or p in seen:
            continue
        seen.add(p)
        positions.append(p)
        labels.append(name)

    # thin labels when too many
    if len(positions) > max_labels:
        step = math.ceil(len(positions) / max_labels)
        labels = [lb if i % step == 0 else "" for i, lb in enumerate(labels)]

    return positions, labels


# ---------------------------------------------------------------------------
# Sequence boundary lines
# ---------------------------------------------------------------------------

def seq_boundary_pixels(filtered_seqs, lookup):
    """Return pixel indices just BEFORE each new sequence starts (for grid lines)."""
    lines = []
    prev_p = None
    for orig_idx, name, length in filtered_seqs:
        p = lookup.get((orig_idx, 0))
        if p is not None and p != prev_p and prev_p is not None:
            lines.append(p - 0.5)
        prev_p = p
    return lines


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot(mat_norm, x_fseqs, x_lookup, n_x, y_fseqs, y_lookup, n_y, args, meta):
    import matplotlib
    if (args.png or args.pdf) and not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    fig_w = max(8, min(24, n_x / 80))
    fig_h = max(8, min(24, n_y / 80))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # imshow — origin='upper' so row 0 is at top (matches Hi-C convention)
    im = ax.imshow(mat_norm, origin="upper", interpolation="none",
                   cmap=args.cmap, vmin=0, vmax=1, aspect="auto")

    # sequence boundary lines
    if args.labels:
        for yline in seq_boundary_pixels(y_fseqs, y_lookup):
            ax.axhline(yline, color="black", linewidth=0.3, alpha=0.4)
        for xline in seq_boundary_pixels(x_fseqs, x_lookup):
            ax.axvline(xline, color="black", linewidth=0.3, alpha=0.4)

    # axis ticks & labels
    if args.labels:
        xtick_pos, xtick_lbl = make_ticks(x_fseqs, x_lookup, n_x)
        ytick_pos, ytick_lbl = make_ticks(y_fseqs, y_lookup, n_y)
        xfs = max(3, min(8, 200 // max(1, len([l for l in xtick_lbl if l]))))
        yfs = max(3, min(8, 200 // max(1, len([l for l in ytick_lbl if l]))))
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(xtick_lbl, rotation=90, fontsize=xfs)
        ax.set_yticks(ytick_pos)
        ax.set_yticklabels(ytick_lbl, fontsize=yfs)
        ax.tick_params(axis="both", which="both", length=2, pad=1)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    # colour bar
    #cbar = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01)
    #cbar.set_label(f"normalised count ({args.norm})", fontsize=7)
    #cbar.ax.tick_params(labelsize=6)

    # title
    title = args.title
    if title is not None:
        ax.set_title(title, fontsize=8, pad=4)

    plt.tight_layout()

    saved = []
    if args.png:
        fig.savefig(args.png, dpi=args.dpi, bbox_inches="tight")
        saved.append(args.png)
    if args.pdf:
        fig.savefig(args.pdf, bbox_inches="tight")
        saved.append(args.pdf)
    if saved:
        print(f"[M] written: {', '.join(saved)}", file=sys.stderr)
    if args.show or not saved:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    print(f"[M] reading {args.input}", file=sys.stderr)
    seqs, bin_size, contacts, meta = read_hicmap(args.input)
    print(f"[M] {len(seqs)} sequences, bin_size={bin_size:,}, "
          f"{len(contacts):,} contact entries", file=sys.stderr)

    min_merge = args.merge if args.merge is not None else bin_size

    # store bin_size in meta for title
    meta["bin_size"] = f"{bin_size:,}"

    if args.x_seqs:
        # ---- filtered / reordered mode ----
        all_seq_names = [name for name, _ in seqs]
        x_names = parse_seq_arg(args.x_seqs, all_seq_names)
        y_names = parse_seq_arg(args.y_seqs, all_seq_names) if args.y_seqs else x_names

        x_fseqs = filter_seqs(seqs, x_names)
        y_fseqs = filter_seqs(seqs, y_names)

        n_x, x_lookup, _ = build_bin_layout_filtered(x_fseqs, bin_size, min_merge)
        n_y, y_lookup, _ = build_bin_layout_filtered(y_fseqs, bin_size, min_merge)
        print(f"[M] pixel grid: {n_x} × {n_y} "
              f"(merge threshold: {min_merge:,} bp)", file=sys.stderr)

        mat = assemble_matrix_xy(contacts, x_lookup, y_lookup, n_x, n_y)
    else:
        # ---- original symmetric mode ----
        n_pixels, lookup, pixels = build_bin_layout(seqs, bin_size, min_merge)
        print(f"[M] pixel grid: {n_pixels} × {n_pixels} "
              f"(merge threshold: {min_merge:,} bp)", file=sys.stderr)

        mat = assemble_matrix(contacts, lookup, n_pixels)
        n_x = n_y = n_pixels
        x_lookup = y_lookup = lookup
        # convert to filtered format for plot / tick helpers
        x_fseqs = y_fseqs = [(i, name, length)
                              for i, (name, length) in enumerate(seqs)]

    print(f"[M] contact matrix assembled, "
          f"non-zero cells: {np.count_nonzero(mat):,}", file=sys.stderr)

    mat_norm = normalise(mat, args.norm, args.clip_pct)

    plot(mat_norm, x_fseqs, x_lookup, n_x, y_fseqs, y_lookup, n_y, args, meta)


if __name__ == "__main__":
    main()
