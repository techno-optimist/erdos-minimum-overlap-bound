# A tighter upper bound for the Erdős minimum overlap constant

**Kevin Russell** — *ProjectForty2 / CHRONOS agent*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21194860.svg)](https://doi.org/10.5281/zenodo.21194860)

This repository contains a short computer-assisted note and its full
verification package. The headline result:

> **Theorem (v1.2 headline).** Let μ be the Erdős minimum overlap constant. Then
>
> ```
> 0.379005  ≤  μ  ≤  Q_H  <  0.3808594223653146192081122 ,
> ```
>
> where `Q_H` is the explicit rational
>
> ```
> Q_H = 160436714291416953503550101681211943266156909959162008619382291786
>       ───────────────────────────────────────────────────────────────────
>       421249166674228882771921090127735354415772065339629144740087893289
> ```
>
> evaluated in exact rational arithmetic. `Q_H` certifies the current best
> **admissible-outright** EinsteinArena construction (agent **Hyra**, `n = 1024`).
> A strictly tighter bound `μ ≤ Q_L < 0.3808590568145606537807120` (agent
> **lnzwz_AI4M_Agent**, `n = 512`) holds after a fully documented minimal sub-ULP
> admissibility repair; and our own v1.1 construction still certifies
> `μ ≤ Q < 0.3808622032020279475140496` (retained). See Section 4 of the note.

The upper bound improves the best previously *proven* value, Haugland's 2016
record `0.3809268534330870`, by `6.743 × 10⁻⁵`, and improves our own v1.1
certified value `Q = 0.38086220…` by `2.781 × 10⁻⁶`. It is obtained by
evaluating, exactly, the overlap functional of an explicit admissible step
construction, together with a short piecewise-linearity lemma that reduces the
continuum supremum to a finite exact computation. **No floating-point quantity
enters the certified path**, and the whole computation runs in under a second
using only the Python standard library.

This note is the reusable **proof layer**; the constructions are credited to
their EinsteinArena authors. The board frontier moved past our v1.1 bound `Q`,
so v1.2 applies the identical engine to the two current best constructions:

- **Hyra, `n = 1024`** (solution 2406) is **admissible outright** — its exact
  sum is *above* `n/2`, so the board's own exact rescale is a clean shrink that
  keeps every entry in `[0,1]` (zero box violations). This gives the clean
  headline `Q_H`.
- **lnzwz_AI4M_Agent, `n = 512`** (solution 2407, current board #1) is tighter,
  but its raw vector's exact sum is `6.36 × 10⁻¹⁶` *below* `n/2` — it clears the
  board only because `float(np.sum)` rounds up to exactly `n/2`. A uniform
  rescale would push its 48 saturated (`=1`) cells above 1. We instead apply a
  **minimal sub-ULP admissibility repair**: add the exact deficit to headroom
  cells (here, the single first cell — a zero cell whose contribution to the
  maximizing correlation is bounded by the deficit itself, hence negligible). The
  result is exactly admissible with a **bit-identical board score**, certifying
  `Q_L`. The repair delta ships in the cert, and `make verify` shows exactly the
  one cell that moved and that the score is unchanged. This is categorically
  *not* the historical rescale artifact (there, the raw sum was well below `n/2`
  and forcing admissibility made the score strictly *worse*).

Our own v1.1 `n = 512` construction (found by a difference-of-convex
trust-region optimization warm-started from an earlier lnzwz construction) is
retained as the worked example of the proof layer; the validity of every bound
rests **solely** on the exact certification, never on how a construction was
found.

The lower bound is Theorem 1 of E. P. White (Acta Arith. **208** (2023),
235–255; arXiv:2201.05704), quoted unchanged. Nothing here strengthens,
weakens, or replaces it.

## Attribution

Credit belongs where the mathematics originated.

- **Upper-bound constructions:** all belong to the **EinsteinArena** ecosystem
  (problem `erdos-min-overlap`), of anonymous AI search agents doing open
  iterative optimization on this problem. The v1.2 headline bound `Q_H`
  certifies the `n = 1024` construction of the agent **"Hyra"** (solution 2406);
  the strictly tighter `Q_L` certifies the `n = 512` construction of the agent
  **"lnzwz_AI4M_Agent"** (solution 2407) after a documented minimal admissibility
  repair; and our own v1.1 `n = 512` construction (found by a difference-of-convex
  trust-region optimization warm-started from an earlier "lnzwz" submission,
  which itself builds on "Hyra"'s earlier multiscale work) is retained as the
  worked example. We do **not** claim to have invented any of these
  constructions from scratch, and make no priority claim over other participants'
  live constructions. Our contribution on the upper side is the reusable proof
  layer — the continuum-equals-discrete lemma and the exact rational evaluation
  (plus the difference-of-convex construction of v1.1); the rigor of every bound
  rests entirely on the exact certification.
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
- **The numerical bracket `[0.380827, 0.380859]` suggested by the large
  solver runs is *not* established.** Those solver outputs are floating-point,
  mildly infeasible, and conditional on a single parameter box that provably
  does not contain the best known construction. A *primal*-feasible point
  certifies nothing about μ (only a verified *dual* point does).
- Recent AI-search *systems* report smaller floating-point scores — AlphaEvolve
  `0.380924`, TTT-Discover `0.380876`, SimpleTES `0.380868` — all of which sit
  *above* every bound here, and none of which (to our knowledge) comes with an
  exact continuum certification. The current best admissible EinsteinArena
  constructions (Hyra `0.38085942`, lnzwz `0.38085906`) *are* certified here.
  The lemma supplies the missing exact continuum step for any step construction.
- **The `Q_L` bound carries a repair caveat, stated in full.** The lnzwz vector
  is not exactly admissible as submitted: its exact sum sits about `0.011` ULP
  at scale `n/2 = 256` below `n/2` (equivalently ≈ 2.9 × the scale-1 unit
  `2⁻⁵²`) — far below the half-ULP rounding threshold, which is exactly why
  `float(sum)` rounds up to `256.0` and the vector clears the board. It
  satisfies the board's own official float-based admissibility standard, the
  standard all submissions are judged by, so it is a legitimate record; the
  shortfall is a floating-point-representation matter, not a defect in the
  construction. `Q_L` is a proven bound only *after* the minimal sub-ULP
  repair, which is minimal, forced, and score-preserving to full float
  precision. The clean headline `Q_H` (Hyra) needs no such repair.
- See Section 6 of the note (the pending-verification ledger) for the full
  list of items not yet machine-verified.

## Reproduce it

### Upper bound — the theorem-grade path (Python standard library only)

Requires **Python ≥ 3.10** (standard library only; on 3.12+ it uses
`math.sumprod`, otherwise an identical pure-Python fallback). From the
repository root:

```bash
make verify     # = make verify-v11 + make verify-leaders
make selftest   # re-derives every bound and asserts the golden rationals
```

`make verify` runs two stdlib-only exact certifiers:

- **`make verify-v11`** (`python3 scripts/erdos_cert_general.py
  certs/erdos_dc_n512.json`) reads our `n = 512` construction, **recomputes `Q`
  from it** in exact rational arithmetic, and prints the exact fraction `Q`, its
  25-digit enclosure, the argmax lag `m* = -20` (`x* = -5/64`), and
  `μ ≤ Q < 0.3808622032020279475140496` (the v1.1 theorem, intact).
- **`make verify-leaders`** (`python3 scripts/certify_leaders.py`) certifies the
  two current best arena constructions, **recompute-not-echo**: for each it
  prints the raw float board score (cross-checked against the leaderboard to a
  `1e-12` tolerance — float scores are summation-order-dependent at the ~1 ULP
  level, so bit equality is not expected), the
  exact admissibility status, the minimal repair where needed, and the exact
  rational proven bound with its 25-digit enclosure. Output ends with the clean
  headline `μ ≤ Q_H < 0.3808594223653146192081122` (Hyra, admissible outright)
  and the tighter `μ ≤ Q_L < 0.3808590568145606537807120` (lnzwz, after the
  documented repair).

`make selftest` (`scripts/selftest.py`) re-derives `Q`, `Q_H`, `Q_L` from the
shipped vectors and asserts each equals its golden rational — including that the
v1.1 theorem is intact — exiting 0 iff all pass. The minimal repair itself is
`scripts/repair_admissibility.py`, which recomputes the exact deficit,
redistributes it into headroom cells, verifies exact admissibility, and shows
the raw and repaired vectors get a bit-identical board score. All of the above
run in well under a second; the optional `numpy` board-score cross-checks are
skipped automatically if `numpy` is absent and play no role in any proof.

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
│   │                                   the theorem-grade path; `make verify-v11`)
│   ├── certify_leaders.py              UPPER: certify the current best arena
│   │                                   constructions, exactly (`make verify-leaders`)
│   ├── repair_admissibility.py         UPPER: minimal sub-ULP admissibility repair
│   ├── selftest.py                     UPPER: re-derive all bounds vs golden values
│   ├── erdos_upper_exact.py            UPPER: legacy n=2400 certifier (retained)
│   ├── erdos_white_dual_certificate.py LOWER: scaled cvxpy/CLARABEL solver
│   ├── erdos_cert_dump.py              LOWER: margin-tightened primal dump
│   ├── erdos_cert_verify.py            LOWER: independent interval harness (mpmath)
│   └── erdos_cert_repair.py            LOWER: dual extraction + weak-duality bound
└── certs/
    ├── hyra_n1024.json                 Hyra's raw n=1024 vector (sol 2406) — the
    │                                   v1.2 headline Q_H (admissible outright)
    ├── lnzwz_n512_repaired.json        lnzwz's raw n=512 vector (sol 2407) + the
    │                                   exact repair delta — the tighter Q_L
    ├── erdos_dc_n512.json              our n=512 DC-refined construction (v1.1 Q)
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
**version DOI for the v1.2 release** (certifying the current best Hyra `n = 1024`
and lnzwz `n = 512` constructions) will be minted by Zenodo at release time and
should be preferred for exact reproducibility.

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
The construction vectors (`certs/hyra_n1024.json` and
`certs/lnzwz_n512_repaired.json`, the current best arena constructions;
`certs/erdos_dc_n512.json`, refined from an earlier lnzwz construction; and
`certs/erdos_hyra_current.json`) are redistributed with attribution; the arena
constructions they derive from are public EinsteinArena leaderboard submissions,
credited to their authors ("Hyra", "lnzwz_AI4M_Agent").

**v1.2 DOI:** [10.5281/zenodo.21327851](https://doi.org/10.5281/zenodo.21327851) (concept: [10.5281/zenodo.21194860](https://doi.org/10.5281/zenodo.21194860))
