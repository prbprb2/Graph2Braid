# ================================================================
# Braid Word Calculator for Triangular Graphs
# Based on: Ahmadi & Wocjan (2010),
#   "On the Quantum Complexity of Evaluating the Tutte Polynomial"
#   J. Knot Theory Ramifications 19(6), 727–737.
#
# Computes braid b ∈ B_{2m} whose plat closure = D(G),
# the link diagram associated to the medial graph of G.
#
# HOW TO USE ON GOOGLE COLAB
# 1. Go to colab.research.google.com → New notebook
# 2. Paste CELL 1 into the first code cell and run it (sets up algorithm)
# 3. Paste CELL 2 into a new cell, edit your graph, and run it
# 4. Paste CELL 3 into a new cell and run it to see the result
# 5. File → Save a copy in Drive, then Share the link with your group
# ================================================================


# ================================================================
# CELL 1 — Algorithm implementation  (run once; do not modify)
# ================================================================

def compute_braid(vertices, edges_input, verbose=False):
    """
    Algorithm from Ahmadi & Wocjan (2010), Section 4.

    Parameters
    ----------
    vertices     : list of (row, col) integer pairs.
                   Row 0 is the topmost row; columns increase left→right.
    edges_input  : list of edge specifications, one of:
                     ((r1,c1), (r2,c2))            weight defaults to 1
                     ((r1,c1), (r2,c2), weight)    explicit integer weight
                   Rules:
                   • Vertical edges (c1 == c2) are ignored automatically.
                   • Parallel edges between the same pair are merged
                     by summing their weights.
    verbose      : set True to print a step-by-step trace.

    Returns
    -------
    braid_word   : list of (generator_index, exponent) pairs
                   e.g. [(2, -1), (2, -1)] means σ_2⁻¹ · σ_2⁻¹
    num_strands  : the braid lives in B_{num_strands}
    """
    vertex_set = {tuple(v) for v in vertices}

    # Build edge dict: frozenset({v1, v2}) → weight
    # Ignore vertical edges; merge parallel edges by summing weights
    edge_dict = {}
    for e in edges_input:
        v1, v2 = tuple(e[0]), tuple(e[1])
        weight  = int(e[2]) if len(e) > 2 else 1
        if v1[1] == v2[1]:              # same column = vertical
            continue
        key = frozenset([v1, v2])
        edge_dict[key] = edge_dict.get(key, 0) + weight

    if not vertex_set:
        return [], 0

    n           = max(r for r, c in vertex_set) + 1   # number of rows
    m           = max(c for r, c in vertex_set) + 1   # number of columns
    num_strands = 2 * m                                # braid in B_{2m}

    # ── helpers ────────────────────────────────────────────────────────

    def is_last_in_col(i, j):
        """True iff (i,j) is the vertex with the highest row index in col j."""
        rows = [r for r, c in vertex_set if c == j]
        return (not rows) or (i == max(rows))

    def nonvert_neighbors(v):
        """Non-vertical neighbours of v, with edge weights."""
        result = []
        for key, w in edge_dict.items():
            if v in key:
                rest = key - {v}
                if rest:
                    (u,) = rest
                    result.append((u, w))
        return result

    # ── mutable algorithm state ────────────────────────────────────────

    marked_v = set()    # vertices that have been marked
    marked_e = set()    # edges that have been marked
    word     = []       # accumulates (sigma_index, exponent) pairs

    def log(msg):
        if verbose: print(msg)

    # ── subroutines from the paper ─────────────────────────────────────

    def Vertex(i, j):
        """Process vertex (i-1, j) and its upward chain."""
        log(f"  Vertex({i},{j})  →  processing ({i-1},{j})")
        # Step 1: exit if (i-1,j) already processed
        if (i-1, j) in marked_v:
            log(f"    ({i-1},{j}) already marked — exit"); return
        # Step 2: recurse two rows up
        if i > 1:
            Vertex(i-2, j)
        # Step 3: call Edge for every non-vertical edge from (i-1,j)
        if i > 0:
            for (k, l), _ in nonvert_neighbors((i-1, j)):
                if j <= l: Edge((i-1, j), (k, l))
                else:      Edge((k, l),   (i-1, j))
        # Step 4: mark (i-1,j) unless it is the last in its column
        if not is_last_in_col(i-1, j):
            marked_v.add((i-1, j))
            log(f"    mark ({i-1},{j})")
        # Step 5: record positive crossing σ_{2j+1}
        log(f"    append  σ_{2*j+1}")
        word.append((2*j+1, 1))

    def Edge(v1, v2):
        """
        Process the non-vertical edge v1–v2.
        Convention: v1 has the lower column index (j ≤ ℓ).
        """
        j1, j2 = v1[1], v2[1]
        key    = frozenset([v1, v2])
        log(f"  Edge({v1}, {v2})")
        # Step 1: exit if already processed or vertical
        if key in marked_e or j1 == j2:
            log(f"    skip (already done or vertical)"); return
        # Step 2: if the col-(j+1) vertex is not in row 0, recurse via Vertex
        if v2[0] > 0:
            Vertex(v2[0], j2)
        # Step 3: mark edge
        marked_e.add(key)
        # Step 4: record negative crossing σ_{j1+j2+1}^{-alpha}
        alpha = edge_dict[key]
        gen   = j1 + j2 + 1
        log(f"    append  σ_{gen}^(-{alpha})  [{v1}–{v2}, weight={alpha}]")
        word.append((gen, -alpha))

    # ── main loop ──────────────────────────────────────────────────────

    for i in range(n):
        log(f"\n{'─'*40} row i = {i} {'─'*40}")

        # (a) Vertex crossings: for each non-last, unmarked vertex in row i
        if i > 0:
            for j in sorted({c for r, c in vertex_set if r == i}):
                v = (i, j)
                if (v in vertex_set
                        and not is_last_in_col(i, j)
                        and v not in marked_v):
                    log(f"  main loop: append σ_{2*j+1} for vertex {v}")
                    marked_v.add(v)
                    word.append((2*j+1, 1))

        # (b) Edge crossings: all edges between row i and row k (k ≥ i)
        for k in range(i, n):
            target = {i, k}
            for key in list(edge_dict):         # list() avoids mutation issues
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


# ── output formatting ──────────────────────────────────────────────────

def format_braid(word, num_strands):
    """Return a human-readable string describing the braid word."""
    header = f"Braid group : B_{num_strands}\n"
    if not word:
        return header + "Braid word  : (identity)"

    # Merge consecutive terms with the same generator index
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
        if exp ==  1: return f"σ_{s}"          # σ_i
        if exp == -1: return f"σ_{s}⁻¹"  # σ_i⁻¹
        if exp  >  1: return f"σ_{s}^{exp}"
        return             f"σ_{s}^({exp})"

    word_str = " · ".join(sym(i, e) for i, e in merged)  # " · "
    raw_str  = ", ".join(f"(σ_{i}, {e})" for i, e in merged)
    return (header
            + f"Braid word  : {word_str}\n"
            + f"Raw pairs   : [{raw_str}]")


# ── graph validator ────────────────────────────────────────────────────

def validate_graph(vertices, edges_input):
    """
    Check Definition 2.1 conditions.
    Returns (is_valid: bool, messages: list[str]).
    """
    vertex_set = {tuple(v) for v in vertices}
    msgs = []

    # Rule 3: if (i,j) is a vertex with i>0, then (i-1,j) must exist
    for i, j in vertex_set:
        if i > 0 and (i-1, j) not in vertex_set:
            msgs.append(
                f"  ✗ Vertex ({i},{j}) is present but ({i-1},{j}) is missing "
                f"(violates Definition 2.1, rule 3)."
            )

    # Rule 6: non-vertical edges must connect adjacent columns
    for e in edges_input:
        c1, c2 = int(e[0][1]), int(e[1][1])
        if c1 != c2 and abs(c1 - c2) != 1:
            msgs.append(
                f"  ✗ Edge {e[0]}–{e[1]}: columns must differ by 1 "
                f"(got {c1} and {c2})."
            )

    return len(msgs) == 0, msgs


def draw_braid(braid_word, num_strands):
    """
    ASCII braid diagram.  Each crossing uses 5 rows of smooth diagonals.
    The under-strand has a gap at the crossing point to show which is on top.

    Positive σ_i  (left over right):   Negative σ_i⁻¹ (right over left):
      \\   /                               /   \\
       \\ /                                / \\
        \\     ← gap for under-strand       /     ← gap for under-strand
       / \\                                \\ /
      /   \\                              \\   /
    """
    n = num_strands
    w = 6                       # columns between adjacent strand positions
    W = w * (n - 1) + 1         # total display width

    def make_row(chars):
        r = [' '] * W
        for col, ch in chars.items():
            if 0 <= col < W:
                r[col] = ch
        return ''.join(r)

    def straight():
        return make_row({w * i: '|' for i in range(n)})

    def cross(idx, positive):
        """Return 5 display rows for crossing at generator idx (1-indexed)."""
        c1, c2 = w * (idx - 1), w * idx   # strand column positions
        base = {w*i: '|' for i in range(n) if w*i not in (c1, c2)}
        mid  = c1 + w // 2                 # = c2 - w//2  (centre column)
        return [
            make_row({**base, c1+1: '\\', c2-1: '/'}),   # approach
            make_row({**base, c1+2: '\\', c2-2: '/'}),   # approach
            make_row({**base, mid:  '\\' if positive else '/'}),  # crossing point
            make_row({**base, c1+2: '/',  c2-2: '\\'}),  # depart
            make_row({**base, c1+1: '/',  c2-1: '\\'}),  # depart
        ]

    def lbl(idx, positive):
        return f'σ_{idx}' + ('' if positive else '⁻¹')

    hdr = make_row({w * i: str(i + 1) for i in range(n)})
    out = [hdr, straight()]

    for idx, exp in braid_word:
        pos  = exp > 0
        rows = cross(idx, pos)
        tag  = lbl(idx, pos)
        for k in range(abs(exp)):
            for j, r in enumerate(rows):
                ann = f'   ← {tag}' if (j == 2 and k == 0) else ''
                out.append(r + ann)
            out.append(straight())

    out.append(hdr)
    print('\n'.join(out))


print("✓ Algorithm loaded.  Now run Cell 2 to define your graph.")


# ================================================================
# CELL 2 — Define your triangular graph   ← EDIT THIS CELL ONLY
# ================================================================
#
# COORDINATE SYSTEM
# ─────────────────
#   (row, col)  –  row 0 is the TOP row, rows increase downward.
#                  col 0 is the leftmost column, cols increase rightward.
#
#   Example layout (2 rows × 2 cols):
#
#       col:   0      1
#   row 0:  (0,0) — (0,1)
#   row 1:  (1,0)
#
# VERTEX LIST  –  list every (row, col) that is a vertex.
#   Rule: if (i,j) is listed with i > 0, then (i-1, j) must also be listed.
#
# EDGE LIST  –  list every edge (vertical ones are auto-ignored).
#   Format:  ((row1, col1), (row2, col2))           weight = 1  (default)
#            ((row1, col1), (row2, col2), weight)   explicit integer weight
#
# ── YOUR GRAPH (edit below) ───────────────────────────────────────────

vertices = [
    (0, 0),
    (0, 1),
    (1, 0),
]

edges = [
    ((0,0), (0,1)),   # non-vertical edge
    ((0,0), (1,0)),   # vertical — automatically ignored
    ((1,0), (0,1)),   # non-vertical edge
]

# ── OTHER EXAMPLES (uncomment one block to try it) ────────────────────

# Single edge:
# vertices = [(0,0), (0,1)]
# edges    = [((0,0), (0,1))]

# Double edge (weight 2):
# vertices = [(0,0), (0,1)]
# edges    = [((0,0), (0,1), 2)]

# 3-row staircase:
# vertices = [(0,0),(1,0),(2,0),(0,1),(1,1),(0,2)]
# edges = [
#     ((0,0),(1,0)), ((1,0),(2,0)), ((0,1),(1,1)),  # verticals (ignored)
#     ((0,0),(0,1)),
#     ((1,0),(0,1)), ((1,0),(1,1)),
#     ((2,0),(1,1)),
#     ((0,1),(0,2)),
#     ((1,1),(0,2)),
# ]

print("Graph defined.  Now run Cell 3 to compute the braid.")


# ================================================================
# CELL 3 — Run the algorithm and display result
# ================================================================

# Validate
ok, issues = validate_graph(vertices, edges)
if ok:
    print("✓ Graph passes validation checks.\n")
else:
    print("Validation warnings (fix before trusting the result):")
    for msg in issues:
        print(msg)
    print()

# Compute
braid_word, num_strands = compute_braid(vertices, edges, verbose=False)

# Display
print(format_braid(braid_word, num_strands))

# Diagram
print()
draw_braid(braid_word, num_strands)

# ── Optional: step-by-step trace ─────────────────────────────────────
# braid_word, num_strands = compute_braid(vertices, edges, verbose=True)
# print(format_braid(braid_word, num_strands))


# ================================================================
# CELL 4 — Borromean rings  (reference example; run independently)
# ================================================================
#
# The Borromean rings is naturally the trace closure of (σ_1 σ_2⁻¹)³
# in B_3.  B_3 has 3 strands (odd), so it cannot arise from the
# paper's triangular/circular graph algorithm, which always outputs
# braids in B_{2m} (even strands).
#
# Fix: add one extra generator σ_3 (Markov stabilisation).
# The trace closure of  (σ_1 σ_2⁻¹)³ · σ_3  in B_4 is still the
# Borromean rings.  Adding σ_3 merges strands 3 and 4 into one
# component, reducing the component count from 4 to 3.
#
#   Without σ_3:  trace closure in B_4  →  4 components  (≠ Borromean rings)
#   With    σ_3:  trace closure in B_4  →  3 components  (= Borromean rings ✓)

def trace_closure_components(braid, n):
    """Count components in the trace closure of a braid in B_n."""
    perm = list(range(n))
    for idx, exp in braid:
        for _ in range(abs(exp)):
            i = idx - 1
            perm[i], perm[i+1] = perm[i+1], perm[i]
    visited = [False] * n
    cycles  = 0
    for start in range(n):
        if not visited[start]:
            cycles += 1
            j = start
            while not visited[j]:
                visited[j] = True
                j = perm[j]
    return cycles

# Borromean rings braid in B_4
borromean = [(1,1),(2,-1),(1,1),(2,-1),(1,1),(2,-1),(3,1)]

print("=" * 50)
print("Borromean rings  —  B_4 braid (Markov stabilisation)")
print("=" * 50)
print(format_braid(borromean, num_strands=4))
print()
print("Closure type   : trace closure  (not plat)")
c = trace_closure_components(borromean, 4)
print(f"Components     : {c}  {'✓' if c == 3 else '✗'}")
print()
print("Linking numbers (any two components):")
# Track which component each strand belongs to after closure,
# then tally signed crossings between component pairs.
perm = list(range(4))
strand_pos = list(range(4))          # strand_pos[s] = current position of strand s
pos_strand  = list(range(4))          # pos_strand[p] = strand currently at position p

# Assign component labels from trace-closure cycles
p2 = list(range(4))
for idx, exp in borromean:
    for _ in range(abs(exp)):
        i = idx - 1
        p2[i], p2[i+1] = p2[i+1], p2[i]
visited = [False]*4
comp = [0]*4
c_id = 0
for start in range(4):
    if not visited[start]:
        j = start
        while not visited[j]:
            visited[j] = True
            comp[j] = c_id
            j = p2[j]
        c_id += 1

# Now replay braid tracking which component is at each position
pos_comp = list(comp)               # component label at each braid position
from collections import defaultdict
lk = defaultdict(int)

for idx, exp in borromean:
    sign = 1 if exp > 0 else -1
    i    = idx - 1                  # left position
    c1, c2 = pos_comp[i], pos_comp[i+1]
    if c1 != c2:
        pair = tuple(sorted([c1, c2]))
        lk[pair] += sign
    pos_comp[i], pos_comp[i+1] = pos_comp[i+1], pos_comp[i]

for pair, total in sorted(lk.items()):
    print(f"  lk(component {pair[0]+1}, component {pair[1]+1}) = {total//2}")
print()
print("All pairwise linking numbers = 0, yet all three are inseparable.")
print("This is the defining property of the Borromean rings. ✓")
print()
draw_braid(borromean, num_strands=4)
