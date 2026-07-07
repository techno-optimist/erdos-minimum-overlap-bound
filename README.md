# A tighter upper bound for the Erdős minimum overlap constant

**Kevin Russell** — *ProjectForty2 / CHRONOS agent*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21194860.svg)](https://doi.org/10.5281/zenodo.21194860)

This repository contains a short computer-assisted note and its full
verification package. The headline result:

> **Theorem.** Let μ be the Erdős minimum overlap constant. Then
>
> ```
> 0.379005  ≤  μ  ≤  Q  <  0.3808622032020279475140496 ,
> ```
>
> where `Q` is the explicit rational
>
> ```
> Q = 117871142698558740618278313
>     ─────────────────────────────
>     309485009821345068724781056
> ```
>
> (in lowest terms; denominator `2⁸⁸`), evaluated in exact rational arithmetic.

The upper bound is new: to our knowledge it is the first *proven* improvement
of the minimum-overlap upper bound since Haugland's 2016 record
`0.3809268534330870`, by `6.4650 × 10⁻⁵`. It also improves our own previous
certified value `Q_old = 0.38086691…` by `4.71 × 10⁻⁶`. It is obtained by
evaluating, exactly, the overlap functional of an explicit admissible
`n = 512` step construction, together with a short piecewise-linearity lemma
that reduces the continuum supremum to a finite exact computation. **No
floating-point quantity enters the certified path**, and the whole computation
runs in under a second using only the Python standard library.

The construction was **found** by a difference-of-convex trust-region
optimization of the (nonconvex) discrete minimax, warm-started from the current
best admissible EinsteinArena construction of the agent "lnzwz". That is a
heuristic search method only: as with any construction, the validity of the
bound rests **solely** on the exact certification (the piecewise-linearity
lemma + exact rational evaluation), never on the optimizer. Our construction is
genuinely distinct from the warm start (Euclidean distance ≈ 0.277 in
`[0,1]⁵¹²`, a separate equioscillation basin), not a local tweak.

The lower bound is Theorem 1 of E. P. White (Acta Arith. **208** (2023),
235–255; arXiv:2201.05704), quoted unchanged. Nothing here strengthens,
weakens, or replaces it.

## Attribution

Credit belongs where the mathematics originated.

- **Upper-bound construction:** the `n = 512` step vector belongs to the
  lineage of the **EinsteinArena** ecosystem (problem `erdos-min-overlap`), an
  ecosystem of anonymous AI search agents doing open iterative optimization on
  this problem. It was found by refining, via a difference-of-convex
  trust-region optimization, the current best admissible arena construction of
  the agent **"lnzwz"**, which in turn builds on the earlier
  multiscale/step-function work of **"Hyra"** and other arena agents. We do
  **not** claim to have invented the construction from scratch, and make no
  priority claim over other participants' live constructions. Our contribution
  on the upper side is threefold: the difference-of-convex optimization that
  improved on the field, the continuum-equals-discrete lemma, and the exact
  rational evaluation — the rigor resting entirely on the last two.
- **Lower-bound method:** the Fourier-analytic convex program and its dual
  verification strategy are entirely **E. P. White's**. Our contribution
  there is limited to scaling his program to larger `N`, extracting and
  repairing certificates, and building an independent interval-arithmetic
  verification harness — none of which improves his unconditional bound.
- **Prior upper-bound record:** J. K. Haugland (2016).

This note was prepared computer-assisted with **CHRONOS**, ProjectForty2's
autonomous research agent, under the author's direction; the author reviewed
the mathematics and takes responsibility for all claims.

## What this package does *not* prove

Honest scoping is the point of the note, so it is worth repeating here:

- **The lower bound is White's, unchanged.** A 2026 preprint of Kim and
  Pilanci (arXiv:2606.31182) reports a certified improvement to `μ ≥ 0.37912`;
  we quote White's peer-reviewed `0.379005` and nothing depends on the choice.
- **The numerical bracket `[0.380827, 0.380862]` suggested by the large
  solver runs is *not* established.** Those solver outputs are floating-point,
  mildly infeasible, and conditional on a single parameter box that provably
  does not contain the best known construction. A *primal*-feasible point
  certifies nothing about μ (only a verified *dual* point does).
- Recent AI-search systems report smaller floating-point scores — AlphaEvolve
  `0.380924`, TTT-Discover `0.380876`, SimpleTES `0.380868`, and the best
  admissible arena float scores (e.g. lnzwz `0.38086279`) — all of which sit
  *above* `Q`, and none of which (to our knowledge) comes with an exact
  continuum certification. The lemma here supplies exactly that missing step
  for any step construction.
- See Section 6 of the note (the pending-verification ledger) for the full
  list of items not yet machine-verified.

## Reproduce it

### Upper bound — the theorem-grade path (Python standard library only)

Requires **Python ≥ 3.10** (standard library only; on 3.12+ it uses
`math.sumprod`, otherwise an identical pure-Python fallback). From the
repository root:

```bash
make verify
# equivalently:
python3 scripts/erdos_cert_general.py certs/erdos_dc_n512.json
```

This reads our `n = 512` construction from `certs/erdos_dc_n512.json`, does
everything in exact rational arithmetic, and **recomputes `Q` from the
construction**, printing the exact fraction `Q`, its 25-digit enclosure, the
argmax lag `m* = -20` (`x* = -5/64`), and the rigorous bound
`μ ≤ Q < 0.3808622032020279475140496`. Runs in well under a second. (A `numpy`
cross-check at the end is optional and skipped automatically if `numpy` is
absent; it plays no role in the proof.)

The earlier stdlib script `scripts/erdos_upper_exact.py` certifies the previous
`n = 2400` construction from `certs/erdos_hyra_current.json` to the superseded
value `Q_old < 0.3808669097979875909124431`; it is retained for comparison and
its reference output is `certs/lane_u_exact_output.txt`.

### Lower bound — White's program, scaled and certified (optional)

These reproduce the *feasibility-certificate* study of Section 4. They need
extra packages (`pip install -r requirements.txt`): `cvxpy` + `CLARABEL` for
the solver, `mpmath` for the interval harness, `numpy` for I/O.

```bash
# 1) re-solve White's Section-5 program with a strict-feasibility margin
python3 scripts/erdos_cert_dump.py --N 20000 --R 10 --T 5000 \
        --margin 1e-6 --out /tmp/cert_N20000.json

# 2) INDEPENDENT interval-arithmetic verification (constraints transcribed
#    from the paper, not from the solver). Exit 0 iff every inequality is
#    CERTIFIED. ~18 s at N=20000, ~145 s at N=80000.
python3 scripts/erdos_cert_verify.py /tmp/cert_N20000.json \
        --out /tmp/verdict_N20000.json
```

The two verdicts referenced in the note are shipped precomputed:
`certs/erdos_verdict_repaired_N20000.json` (96/96 inequality constraints
certified) and `certs/erdos_verdict_repaired_N80000.json` (176/176). **These
are primal feasibility certificates: they verify the transcription of White's
program at scale, but by weak duality they do *not* bound μ.** Expected output
on success: `erdos_cert_verify.py` exits 0 and its `RIGOR VERDICT` block states
the certificate is a strictly-feasible primal point that is "NECESSARY but NOT
SUFFICIENT" for μ ≥ Ω. That verdict is the intended result of a passing run, not
an error — it restates the weak-duality scoping above. Read the exit code
(0 = every inequality certified), not the tone of the verdict text. The
dual-repair summary (`certs/erdos_repaired_cert_N150000.json`) is float-evaluated
and its interval verification is listed as pending in the note.

## Build the PDF

The prebuilt `erdos_minimum_overlap_bracket.pdf` is included. To rebuild from
source (e.g. with [tectonic](https://tectonic-typesetting.github.io/)):

```bash
make figure                                # regenerates fig_construction.pdf
tectonic erdos_minimum_overlap_bracket.tex
# or: make pdf   (figure + tectonic in one step)
```

## Layout

```
.
├── Makefile                            `make verify` / `make figure` / `make pdf`
├── erdos_minimum_overlap_bracket.tex   the note (source of truth)
├── erdos_minimum_overlap_bracket.pdf   prebuilt PDF
├── fig_construction.pdf                Figure 1 (construction + overlap function)
├── make_figure.py                      regenerates Figure 1 from certs/
├── requirements.txt                    optional deps for the lower-bound scripts
├── scripts/
│   ├── erdos_cert_general.py           UPPER: exact rational evaluation of any
│   │                                   admissible step vector (stdlib only —
│   │                                   the theorem-grade path; `make verify`)
│   ├── erdos_upper_exact.py            UPPER: legacy n=2400 certifier (retained)
│   ├── erdos_white_dual_certificate.py LOWER: scaled cvxpy/CLARABEL solver
│   ├── erdos_cert_dump.py              LOWER: margin-tightened primal dump
│   ├── erdos_cert_verify.py            LOWER: independent interval harness (mpmath)
│   └── erdos_cert_repair.py            LOWER: dual extraction + weak-duality bound
└── certs/
    ├── erdos_dc_n512.json              our n=512 DC-refined construction (the
    │                                   certified upper-bound vector)
    ├── erdos_hyra_current.json         Hyra's raw n=2400 vector (superseded;
    │                                   EinsteinArena board copy, 2026-06-30)
    ├── lane_u_exact_output.txt         reference output of erdos_upper_exact.py
    ├── erdos_verdict_repaired_N20000.json   interval verdict, 96/96 certified
    ├── erdos_verdict_repaired_N80000.json   interval verdict, 176/176 certified
    └── erdos_repaired_cert_N150000.json     dual-repair summary (float; pending)
```

## Cite this work

If you use this note or its verification package, please cite the archived
release. The concept DOI below always resolves to the latest version; the
**version DOI for the v1.1 release** (this n = 512 upper bound) will be minted
by Zenodo at release time and should be preferred for exact reproducibility.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21194860.svg)](https://doi.org/10.5281/zenodo.21194860)

```
Kevin Russell, "A tighter upper bound for the Erdős minimum overlap constant,
with machine-verified feasibility certificates for White's lower-bound program",
ProjectForty2 / CHRONOS, 2026. Concept DOI: 10.5281/zenodo.21194860.
```

## License

Code and certificate data are released under the MIT License (see `LICENSE`).
The note text and figures (`*.tex`, `*.pdf`) are © 2026 Kevin Russell,
released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The construction vectors (`certs/erdos_dc_n512.json`, refined from lnzwz's
construction, and `certs/erdos_hyra_current.json`) are redistributed with
attribution; the arena constructions they derive from are public EinsteinArena
leaderboard submissions.
