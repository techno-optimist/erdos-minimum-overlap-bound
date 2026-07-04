# A tighter upper bound for the Erdős minimum overlap constant

**Kevin Russell** — *ProjectForty2 / CHRONOS agent*

This repository contains a short computer-assisted note and its full
verification package. The headline result:

> **Theorem.** Let μ be the Erdős minimum overlap constant. Then
>
> ```
> 0.379005  ≤  μ  ≤  Q  <  0.3808669097979875909124431 ,
> ```
>
> where `Q` is the explicit rational
>
> ```
> Q = 1424992289798782609633201801352767458976314440679252577
>     ─────────────────────────────────────────────────────────
>     3741444197802851304404516484910431627947663875649308401
> ```
>
> evaluated in exact integer arithmetic.

The upper bound is new: to our knowledge it is the first *proven* improvement
of the minimum-overlap upper bound since Haugland's 2016 record
`0.3809268534330870`, by `5.99 × 10⁻⁵`. It is obtained by evaluating, exactly,
the overlap functional of an explicit `n = 2400` step construction, together
with a short piecewise-linearity lemma that reduces the continuum supremum to
a finite exact computation. **No floating-point quantity enters the certified
path**, and the whole computation runs in under a second using only the Python
standard library.

The lower bound is Theorem 1 of E. P. White (Acta Arith. **208** (2023),
235–255; arXiv:2201.05704), quoted unchanged. Nothing here strengthens,
weakens, or replaces it.

## Attribution

Credit belongs where the mathematics originated.

- **Upper-bound construction:** the `n = 2400` step vector was found by
  **"Hyra"**, an anonymous AI search agent competing on the EinsteinArena
  platform (problem `erdos-min-overlap`), within an ecosystem of AI agents
  doing open iterative optimization on this problem. The construction is
  Hyra's; this note makes no priority claim over other participants' live
  constructions. Our contribution on the upper side is the
  continuum-equals-discrete lemma and the exact rational evaluation.
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
- **The numerical bracket `[0.380827, 0.380867]` suggested by the large
  solver runs is *not* established.** Those solver outputs are floating-point,
  mildly infeasible, and conditional on a single parameter box that provably
  does not contain the best known construction. A *primal*-feasible point
  certifies nothing about μ (only a verified *dual* point does).
- Recent AI-search systems report smaller floating-point scores — AlphaEvolve
  `0.380924`, TTT-Discover `0.380876`, SimpleTES `0.380868` — all of which sit
  *above* `Q`, and none of which (to our knowledge) comes with an exact
  continuum certification. The lemma here supplies exactly that missing step
  for any step construction.
- See Section 6 of the note (the pending-verification ledger) for the full
  list of items not yet machine-verified.

## Reproduce it

### Upper bound — the theorem-grade path (Python standard library only)

Requires **Python ≥ 3.12** and nothing else. From the repository root:

```bash
python3 scripts/erdos_upper_exact.py
```

This reads Hyra's vector from `certs/erdos_hyra_current.json`, does everything
in exact integer arithmetic, and prints the exact fraction `Q`, its 25-digit
enclosure, the argmax lag `m* = -92`, and the post-rescale `[0,1]` check (zero
violations). The output must match `certs/lane_u_exact_output.txt` verbatim:

```bash
python3 scripts/erdos_upper_exact.py | diff - certs/lane_u_exact_output.txt && echo OK
```

Runs in well under a second. (A `numpy` cross-check at the end is optional and
skipped automatically if `numpy` is absent; it plays no role in the proof.)

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
program at scale, but by weak duality they do *not* bound μ.** The dual-repair
summary (`certs/erdos_repaired_cert_N150000.json`) is float-evaluated and its
interval verification is listed as pending in the note.

## Build the PDF

The prebuilt `erdos_minimum_overlap_bracket.pdf` is included. To rebuild from
source (e.g. with [tectonic](https://tectonic-typesetting.github.io/)):

```bash
python3 make_figure.py                     # regenerates fig_hyra_construction.pdf
tectonic erdos_minimum_overlap_bracket.tex
```

## Layout

```
.
├── erdos_minimum_overlap_bracket.tex   the note (source of truth)
├── erdos_minimum_overlap_bracket.pdf   prebuilt PDF (11 pp.)
├── fig_hyra_construction.pdf           Figure 1 (construction + overlap function)
├── make_figure.py                      regenerates Figure 1 from certs/
├── requirements.txt                    optional deps for the lower-bound scripts
├── scripts/
│   ├── erdos_upper_exact.py            UPPER: exact rational evaluation
│   │                                   (stdlib only — the theorem-grade path)
│   ├── erdos_white_dual_certificate.py LOWER: scaled cvxpy/CLARABEL solver
│   ├── erdos_cert_dump.py              LOWER: margin-tightened primal dump
│   ├── erdos_cert_verify.py            LOWER: independent interval harness (mpmath)
│   └── erdos_cert_repair.py            LOWER: dual extraction + weak-duality bound
└── certs/
    ├── erdos_hyra_current.json         Hyra's raw n=2400 vector (EinsteinArena
    │                                   board copy, retrieved 2026-06-30)
    ├── lane_u_exact_output.txt         reference output of erdos_upper_exact.py
    ├── erdos_verdict_repaired_N20000.json   interval verdict, 96/96 certified
    ├── erdos_verdict_repaired_N80000.json   interval verdict, 176/176 certified
    └── erdos_repaired_cert_N150000.json     dual-repair summary (float; pending)
```

## License

Code and certificate data are released under the MIT License (see `LICENSE`).
The note text and figures (`*.tex`, `*.pdf`) are © 2026 Kevin Russell,
released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Hyra's construction vector (`certs/erdos_hyra_current.json`) is redistributed
with attribution as a public EinsteinArena leaderboard submission.
