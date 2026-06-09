#!/usr/bin/env python3
"""
Graph2Braid.py  —  Triangular graph → braid word + diagram

Usage
-----
  python Graph2Braid.py graph.json
  python Graph2Braid.py graph.json --verbose
  python Graph2Braid.py graph.json --no-diagram

Input format (JSON)
-------------------
  {
    "vertices": [[0,0], [0,1], [1,0]],
    "edges":    [[[0,0],[0,1]], [[0,0],[1,0]], [[1,0],[0,1]]]
  }

  Each edge is either [[r1,c1],[r2,c2]] (weight 1)
                   or [[r1,c1],[r2,c2], weight].
  Vertical edges (same column) are ignored automatically.

Based on: Ahmadi & Wocjan (2010),
  "On the Quantum Complexity of Evaluating the Tutte Polynomial"
  J. Knot Theory Ramifications 19(6), 727-737.
"""

import argparse
import json
import os
import sys


# ── Algorithm ────────────────────────────────────────────────────────────────

def compute_braid(vertices, edges_input, verbose=False):
    vertex_set = {tuple(v) for v in vertices}

    edge_dict = {}
    for e in edges_input:
        v1, v2 = tuple(e[0]), tuple(e[1])
        weight  = int(e[2]) if len(e) > 2 else 1
        if v1[1] == v2[1]:
            continue
        key = frozenset([v1, v2])
        edge_dict[key] = edge_dict.get(key, 0) + weight

    if not vertex_set:
        return [], 0

    n           = max(r for r, c in vertex_set) + 1
    m           = max(c for r, c in vertex_set) + 1
    num_strands = 2 * m

    def is_last_in_col(i, j):
        rows = [r for r, c in vertex_set if c == j]
        return (not rows) or (i == max(rows))

    def nonvert_neighbors(v):
        result = []
        for key, w in edge_dict.items():
            if v in key:
                rest = key - {v}
                if rest:
                    (u,) = rest
                    result.append((u, w))
        return result

    marked_v = set()
    marked_e = set()
    word     = []

    def log(msg):
        if verbose:
            print(msg)

    # ── Subroutines (Section 4, Ahmadi & Wocjan) ─────────────────────────

    def Vertex(i, j):
        """Process vertex (i-1, j) and its upward chain."""
        log(f"  Vertex({i},{j}) -> processing ({i-1},{j})")
        # Step 1: exit if (i-1, j) already processed
        if (i-1, j) in marked_v:
            return
        # Step 2: recurse two rows up
        if i > 1:
            Vertex(i-2, j)
        # Step 3: call Edge for every non-vertical edge from (i-1, j)
        if i > 0:
            for (k, l), _ in nonvert_neighbors((i-1, j)):
                if j <= l: Edge((i-1, j), (k, l))
                else:      Edge((k, l),   (i-1, j))
        # Step 4: mark (i-1, j) unless it is the last in its column
        if not is_last_in_col(i-1, j):
            marked_v.add((i-1, j))
        # Step 5: record positive crossing sigma_{2j+1}
        log(f"    append sigma_{2*j+1}")
        word.append((2*j+1, 1))

    def Edge(v1, v2):
        """
        Process non-vertical edge v1-v2.
        Convention: v1 has the lower column index (j <= l).
        """
        j1, j2 = v1[1], v2[1]
        key    = frozenset([v1, v2])
        # Step 1: exit if already processed or vertical
        if key in marked_e or j1 == j2:
            return
        # Step 2: if the col-(j+1) vertex is not in row 0, recurse via Vertex
        if v2[0] > 0:
            Vertex(v2[0], j2)
        # Step 3: mark edge
        marked_e.add(key)
        # Step 4: record negative crossing sigma_{j1+j2+1}^{-alpha}
        alpha = edge_dict[key]
        gen   = j1 + j2 + 1
        log(f"    append sigma_{gen}^(-{alpha})  [{v1}-{v2}, weight={alpha}]")
        word.append((gen, -alpha))

    # ── Main loop ─────────────────────────────────────────────────────────

    for i in range(n):
        log(f"\n--- row i={i} ---")

        # (a) Vertex crossings: non-last, unmarked vertices in row i
        if i > 0:
            for j in sorted({c for r, c in vertex_set if r == i}):
                v = (i, j)
                if v in vertex_set and not is_last_in_col(i, j) \
                        and v not in marked_v:
                    log(f"  main: append sigma_{2*j+1} for {v}")
                    marked_v.add(v)
                    word.append((2*j+1, 1))

        # (b) Edge crossings: all edges between row i and row k (k >= i)
        for k in range(i, n):
            target = {i, k}
            for key in list(edge_dict):
                vs = list(key)
                if {vs[0][0], vs[1][0]} == target:
                    if i == k:
                        # same-row edge: lower-column vertex first
                        if vs[0][1] <= vs[1][1]:
                            Edge(vs[0], vs[1])
                        else:
                            Edge(vs[1], vs[0])
                    else:
                        # cross-row edge: row-i vertex first
                        # (so that Vertex is called on the row-k endpoint)
                        if vs[0][0] == i:
                            Edge(vs[0], vs[1])
                        else:
                            Edge(vs[1], vs[0])

    return word, num_strands


# ── Formatting ────────────────────────────────────────────────────────────────

def format_braid(word, num_strands):
    header = f"Braid group : B_{num_strands}\n"
    if not word:
        return header + "Braid word  : (identity)"

    merged = []
    for idx, exp in word:
        if merged and merged[-1][0] == idx:
            merged[-1] = (merged[-1][0], merged[-1][1] + exp)
        else:
            merged.append((idx, exp))
    merged = [(i, e) for i, e in merged if e != 0]

    if not merged:
        return header + "Braid word  : (identity after cancellations)"

    def sym(idx, exp):
        s = str(idx)
        if exp ==  1: return f"sigma_{s}"
        if exp == -1: return f"sigma_{s}^(-1)"
        if exp  >  1: return f"sigma_{s}^{exp}"
        return             f"sigma_{s}^({exp})"

    word_str = " . ".join(sym(i, e) for i, e in merged)
    return header + f"Braid word  : {word_str}\n" \
                  + f"Raw pairs   : {merged}"


# ── Graph diagram ─────────────────────────────────────────────────────────────

def draw_graph(vertices, edges_input):
    """
    Draw the triangular graph on a character grid.

    Vertices are shown as 'o' at their (row, col) grid positions.
    Edge characters:
      -   same-row (horizontal) edge
      |   vertical edge
      /   diagonal going up-right   (row decreases, col increases)
      \\   diagonal going down-right  (row increases, col increases)
    """
    vertex_set = {tuple(v) for v in vertices}
    if not vertex_set:
        return

    n_rows = max(r for r, c in vertex_set) + 1
    n_cols = max(c for r, c in vertex_set) + 1

    w, h = 4, 4                          # display spacing per column / row
    H = h * (n_rows - 1) + 1            # total display height
    W = w * (n_cols - 1) + 1            # total display width

    grid = [[' '] * W for _ in range(H)]

    # ── Place vertices ───────────────────────────────────────────────────
    for r, c in vertex_set:
        grid[r*h][c*w] = 'o'

    # ── Draw edges ───────────────────────────────────────────────────────
    for e in edges_input:
        r1, c1 = int(e[0][0]), int(e[0][1])
        r2, c2 = int(e[1][0]), int(e[1][1])
        y1, x1 = r1*h, c1*w
        y2, x2 = r2*h, c2*w

        if c1 == c2:                     # vertical edge  →  |
            for y in range(min(y1,y2)+1, max(y1,y2)):
                grid[y][x1] = '|'

        elif r1 == r2:                   # same-row edge  →  ---
            for x in range(min(x1,x2)+1, max(x1,x2)):
                grid[y1][x] = '-'

        else:                            # diagonal edge  →  / or \
            char = '\\' if (y2-y1)*(x2-x1) > 0 else '/'
            steps = max(abs(y2-y1), abs(x2-x1))
            for s in range(1, steps):
                y = y1 + round((y2-y1) * s / steps)
                x = x1 + round((x2-x1) * s / steps)
                if 0 <= y < H and 0 <= x < W:
                    grid[y][x] = char

    # ── Build output string with row / column labels ─────────────────────
    lw = 8                               # label column width
    lines = []

    col_hdr = ' ' * lw
    for j in range(n_cols):
        col_hdr += str(j).center(w)
    lines.append(col_hdr.rstrip())

    for y, row in enumerate(grid):
        lbl = f'row {y//h}' if y % h == 0 else ''
        lines.append(f'{lbl:<{lw}}' + ''.join(row).rstrip())

    return '\n'.join(lines)


# ── Braid diagram ──────────────────────────────────────────────────────────────

def draw_braid(word, num_strands):
    n = num_strands
    w = 6
    W = w * (n - 1) + 1

    def make_row(chars):
        r = [' '] * W
        for col, ch in chars.items():
            if 0 <= col < W:
                r[col] = ch
        return ''.join(r)

    def straight():
        return make_row({w*i: '|' for i in range(n)})

    def cross(idx, positive):
        c1, c2 = w*(idx-1), w*idx
        base = {w*i: '|' for i in range(n) if w*i not in (c1, c2)}
        mid  = c1 + w//2
        return [
            make_row({**base, c1+1: '\\', c2-1: '/'}),
            make_row({**base, c1+2: '\\', c2-2: '/'}),
            make_row({**base, mid:  '\\' if positive else '/'}),
            make_row({**base, c1+2: '/',  c2-2: '\\'}),
            make_row({**base, c1+1: '/',  c2-1: '\\'}),
        ]

    def lbl(idx, positive):
        return f"sigma_{idx}" + ("" if positive else "^(-1)")

    hdr = make_row({w*i: str(i+1) for i in range(n)})
    out = [hdr, straight()]

    for idx, exp in word:
        pos  = exp > 0
        rows = cross(idx, pos)
        tag  = lbl(idx, pos)
        for k in range(abs(exp)):
            for j, r in enumerate(rows):
                ann = f"   <- {tag}" if (j == 2 and k == 0) else ""
                out.append(r + ann)
            out.append(straight())

    out.append(hdr)
    return '\n'.join(out)


# ── Figure output ─────────────────────────────────────────────────────────────

def save_figure(graph_str, braid_str, description, command, output_path):
    """
    Save a two-panel JPG:  left = input graph,  right = output braid.
    Requires matplotlib.  Prints a warning and returns silently if unavailable.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')           # non-interactive backend (no display needed)
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
    except ImportError:
        print("Warning: matplotlib not installed — skipping figure output.")
        print("         Install with:  pip install matplotlib")
        return

    graph_lines = graph_str.split('\n')
    braid_lines = braid_str.split('\n')
    n_lines     = max(len(graph_lines), len(braid_lines))
    graph_w     = max((len(l) for l in graph_lines), default=1)
    braid_w     = max((len(l) for l in braid_lines), default=1)

    # Scale figure to content
    fig_h = max(5.0, n_lines  * 0.20 + 2.5)
    fig_w = max(12.0, (graph_w + braid_w) * 0.11 + 2.0)

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

    # ── Title block ───────────────────────────────────────────────────────
    fig.text(0.5, 0.98, description,
             ha='center', va='top', fontsize=13, fontweight='bold')
    fig.text(0.5, 0.93, command,
             ha='center', va='top', fontsize=8.5,
             color='#555555', fontfamily='monospace')

    # ── Two panels ────────────────────────────────────────────────────────
    gs = gridspec.GridSpec(1, 2,
                           left=0.02, right=0.98,
                           top=0.86,  bottom=0.02,
                           wspace=0.06)
    ax_g = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])

    for ax, text, heading in [
        (ax_g, graph_str, 'Input graph'),
        (ax_b, braid_str, 'Output braid'),
    ]:
        ax.set_title(heading, fontsize=11, pad=5, loc='left')
        ax.text(0.03, 0.97, text,
                transform=ax.transAxes,
                fontfamily='monospace',
                fontsize=9,
                verticalalignment='top',
                horizontalalignment='left',
                linespacing=1.25)
        ax.axis('off')
        # subtle border
        rect = plt.Rectangle((0, 0), 1, 1,
                              fill=False, edgecolor='#cccccc',
                              linewidth=0.8, transform=ax.transAxes,
                              clip_on=False)
        ax.add_patch(rect)

    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Figure saved: {output_path}")


# ── Validation ────────────────────────────────────────────────────────────────

def validate_graph(vertices, edges_input):
    vertex_set = {tuple(v) for v in vertices}
    msgs = []
    for i, j in vertex_set:
        if i > 0 and (i-1, j) not in vertex_set:
            msgs.append(f"  vertex ({i},{j}) present but ({i-1},{j}) missing "
                        f"(Definition 2.1 rule 3)")
    for e in edges_input:
        c1, c2 = int(e[0][1]), int(e[1][1])
        if c1 != c2 and abs(c1 - c2) != 1:
            msgs.append(f"  edge {e[0]}-{e[1]}: columns must differ by 1 "
                        f"(got {c1} and {c2})")
    return len(msgs) == 0, msgs


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compute the braid word for a triangular graph (Ahmadi & Wocjan 2010).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("input",
                        help="JSON file containing the graph "
                             "(use - to read from stdin)")
    parser.add_argument("--verbose",    action="store_true",
                        help="print step-by-step Vertex/Edge trace")
    parser.add_argument("--no-graph",   action="store_true",
                        help="skip the input graph diagram")
    parser.add_argument("--no-diagram", action="store_true",
                        help="skip the output braid diagram")
    parser.add_argument("--nofig",      action="store_true",
                        help="do not save a JPG figure")
    args = parser.parse_args()

    # Load JSON
    try:
        if args.input == "-":
            data = json.load(sys.stdin)
        else:
            with open(args.input) as f:
                data = json.load(f)
    except FileNotFoundError:
        sys.exit(f"Error: file not found: {args.input}")
    except json.JSONDecodeError as e:
        sys.exit(f"Error: invalid JSON: {e}")

    vertices    = data.get("vertices",    [])
    edges       = data.get("edges",       [])
    description = data.get("description", args.input)

    if not vertices:
        sys.exit("Error: 'vertices' list is empty or missing")

    # Validate
    ok, issues = validate_graph(vertices, edges)
    if not ok:
        print("Validation warnings:")
        for msg in issues:
            print(msg)
        print()

    # Build diagrams (always computed; display is optional)
    graph_str = draw_graph(vertices, edges)
    word, strands = compute_braid(vertices, edges, verbose=args.verbose)
    braid_str = draw_braid(word, strands)

    # ── Terminal output ───────────────────────────────────────────────────
    if not args.no_graph:
        print("Input graph:\n")
        print(graph_str)
        print()

    print("Output braid:\n")
    print(format_braid(word, strands))

    if not args.no_diagram:
        print()
        print(braid_str)

    # ── Figure output ─────────────────────────────────────────────────────
    if not args.nofig:
        os.makedirs('JPGs', exist_ok=True)
        if args.input == '-':
            stem = 'output'
        else:
            stem = os.path.splitext(os.path.basename(args.input))[0]
        fig_path = os.path.join('JPGs', stem + '.jpg')
        save_figure(graph_str, braid_str,
                    description, ' '.join(sys.argv), fig_path)


if __name__ == "__main__":
    main()
