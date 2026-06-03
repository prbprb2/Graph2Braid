# Graph2Braid

**Convert triangular graphs to braid words — with ASCII diagrams and figure output.**

This is a Python implementation of the braid-construction algorithm from:

> Ahmadi, H. & Wocjan, P. (2010). *On the Quantum Complexity of Evaluating the Tutte Polynomial.*
> Journal of Knot Theory and Its Ramifications, 19(6), 727–737.
> DOI: [10.1142/S021821651000808X](https://doi.org/10.1142/S021821651000808X)

Given a **triangular graph** G (in the sense of that paper), the algorithm outputs
a braid word b in B_{2m} whose **plat closure** equals the link diagram D(G)
associated to the medial graph of G. This connection is the key step in showing
that approximately evaluating the Tutte polynomial of triangular graphs is
BQP-complete at most roots of unity.

---

## What it looks like

```
$ python Graph2Braid.py example01_triangle.json
```

```
Input graph:

         0   1
row 0   o---o
        |  /
        | /
        |/
row 1   o

Output braid:

Braid group : B_4
Braid word  : sigma_2^(-2)
Raw pairs   : [(2, -2)]

1     2     3     4
|     |     |     |
|      \   /      |
|       \ /       |
|        /        |   <- sigma_2^(-1)
|       / \       |
|      /   \      |
|     |     |     |
|      \   /      |
|       \ /       |
|        /        |   <- sigma_2^(-1)
|       / \       |
|      /   \      |
|     |     |     |
1     2     3     4
```

A JPG figure is also saved to the `JPGs/` directory:

![Triangle example](JPGs/example01_triangle.jpg)

---

## Installation

```bash
git clone https://github.com/your-username/Graph2Braid.git
cd Graph2Braid
pip install matplotlib          # only dependency
```

Python 3.8+ required.

---

## Usage

```bash
python Graph2Braid.py graph.json              # full output + JPG figure
python Graph2Braid.py graph.json --verbose    # step-by-step Vertex/Edge trace
python Graph2Braid.py graph.json --no-graph   # skip input graph diagram
python Graph2Braid.py graph.json --no-diagram # skip braid diagram
python Graph2Braid.py graph.json --nofig      # skip JPG output
python Graph2Braid.py -                       # read JSON from stdin
```

Figures are saved to `JPGs/<input-stem>.jpg` (directory created automatically).

---

## Input format (JSON)

```json
{
  "description": "Human-readable name shown in the figure title",
  "vertices": [[0,0], [0,1], [1,0]],
  "edges": [
    [[0,0], [0,1]],
    [[0,0], [1,0]],
    [[1,0], [0,1]]
  ]
}
```

**Coordinates** are `[row, col]` with row 0 at the top, increasing downward.

**Edge format:** each entry is either
- `[[r1,c1], [r2,c2]]` — unit weight
- `[[r1,c1], [r2,c2], weight]` — explicit integer weight (parallel edges)

**Vertical edges** (same column) are included for completeness but are
automatically ignored by the algorithm, as stated in the paper.

**Validity rule:** if `(i,j)` is a vertex with `i > 0`, then `(i-1,j)` must
also be a vertex (Definition 2.1, rule 3). The tool will warn if this is violated.

---

## Examples

Nineteen worked examples are included as JSON files:

| File | Graph | Braid group | Braid word |
|------|-------|-------------|------------|
| `example01_triangle.json` | Triangle | B₄ | σ₂⁻² |
| `example02_single_edge.json` | Single edge | B₄ | σ₂⁻¹ |
| `example03_staircase.json` | 3-row staircase | B₆ | σ₂⁻¹·σ₄⁻¹·σ₂⁻¹·σ₄⁻¹·σ₁·σ₃·σ₂⁻² |
| `example04_double_edge.json` | Weight-2 edge | B₄ | σ₂⁻² |
| `example05_single_vertex.json` | Single vertex | B₂ | (identity) |
| `example07_horizontal.json` | Horizontal edge | B₄ | σ₂⁻¹ |
| `example12_fan3.json` | Fan (3 spokes) | B₄ | σ₂⁻³ |
| `example14_three_cols.json` | Row of 3 | B₆ | σ₂⁻¹·σ₄⁻¹ |
| … | … | … | … |

Run all examples at once:

```bash
for f in example*.json; do python Graph2Braid.py $f; done
```

Bind all figures in `JPGs/` into a single PDF:

```bash
python rebind_pdf.py
# produces JPGs/all_examples.pdf
```

---

## Braid diagram conventions

```
Positive σ_i  (left strand over right):   Negative σ_i⁻¹ (right strand over left):

  \   /                                      /   \
   \ /                                        / \
    \    ← only the over-strand appears        /    ← only the over-strand appears
   / \                                        \ /
  /   \                                      \   /
```

The under-strand has a gap at the crossing point.
Odd-indexed generators (σ₁, σ₃, …) appear with positive exponent (vertex crossings);
even-indexed generators (σ₂, σ₄, …) appear with negative exponent (edge crossings).
This guarantees the plat closure is an alternating link (Lemma 4.1 of the paper).

---

## Borromean rings note

The Borromean rings is the trace closure of (σ₁·σ₂⁻¹)³·σ₃ in B₄.
It does **not** arise from a triangular graph — the paper's algorithm outputs
braids in B_{2m} (even strands) for plat closures, while the Borromean rings
naturally lives in B₃ (odd strands). The extra σ₃ is a Markov stabilisation.
This is demonstrated in the companion Colab notebook `braid_calculator.py`.

---

## Google Colab notebook

`braid_calculator.py` is structured as a four-cell Colab notebook:

| Cell | Purpose |
|------|---------|
| Cell 1 | Algorithm functions — run once |
| Cell 2 | Define your graph (edit this cell) |
| Cell 3 | Compute and display the braid |
| Cell 4 | Borromean rings reference example |

Open it in Colab, paste each cell, and share the notebook link with your group.

---

## Files

```
Graph2Braid.py          CLI tool — graph JSON → braid word + figures
rebind_pdf.py           Bind all JPGs in JPGs/ into a single PDF
braid_calculator.py     Google Colab notebook version
example*.json           19 worked example graphs
JPGs/                   Output figures (created automatically)
```

---

## Background

The Tutte polynomial T(G; x, y) is a two-variable graph invariant that encodes
the chromatic polynomial, the number of spanning trees, and much more.
Ahmadi & Wocjan show that approximately evaluating T(G; q, 1/q) for triangular
graphs G is BQP-complete at most roots of unity q, by connecting it to the
Jones polynomial of the plat closure of a braid constructed from G.
This tool implements that construction explicitly.

---

## License

MIT License — see `LICENSE`.
