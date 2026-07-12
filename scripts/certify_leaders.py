#!/usr/bin/env python3
"""Certify the current best EinsteinArena minimum-overlap constructions, exactly.

For each construction this RECOMPUTES (never echoes) from the raw board vector:

  1. the raw float board score (EinsteinArena semantics), for cross-check
     against the leaderboard value;
  2. the exact rational sum and its relation to n/2 (admissible outright /
     surplus / deficit);
  3. the admissibility route:
       - EXACT   : sum is already n/2 -- admissible as given;
       - SURPLUS : sum above n/2 -- the board's exact uniform rescale
                   f <- f*(n/2)/sum shrinks every entry (<1), stays in [0,1]
                   with zero box violations, sum exactly n/2.  No asterisk;
       - DEFICIT : sum below n/2 -- uniform rescale would push saturated cells
                   above 1 (box violation), so we apply the minimal
                   sub-ULP admissibility repair of repair_admissibility.py
                   (add the exact deficit to headroom cells; board score
                   bit-identical);
  4. the exact rational proven bound  mu <= Q_i  and its 25-digit round-up,
     with the argmax lag m* and x* = 2 m*/n.

Everything on the certified path is exact integer/rational arithmetic (stdlib
``fractions`` only).  The float board score is an optional numpy cross-check
and plays no role in any proof.

Usage:
    python3 certify_leaders.py            # certify the shipped leader certs
    python3 certify_leaders.py a.json b.json ...   # certify given cert files
"""
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(os.path.dirname(HERE), "certs")
sys.path.insert(0, HERE)
from repair_admissibility import (          # noqa: E402
    dyadic_ints, minimal_repair, exact_score, board_score_float,
    repaired_values, fmt_down, fmt_up, DIGITS,
)

# v1.1 published bound (recomputed by `make verify`); used here only for context.
Q_V11 = Fraction(117871142698558740618278313, 309485009821345068724781056)

DEFAULT_CERTS = [
    os.path.join(CERTS, "hyra_n1024.json"),
    os.path.join(CERTS, "lnzwz_n512_repaired.json"),
]


def certify(path):
    obj = json.load(open(path))
    vals = obj["values"]
    n = len(vals)
    A, L = dyadic_ints(vals)
    cap = 1 << L
    raw_sum = Fraction(sum(A), cap)
    half = Fraction(n, 2)

    agent = obj.get("agent", "?")
    sol = obj.get("solId", "?")
    label = f"{agent}  (n={n}, solId {sol})"
    print("=" * 74)
    print(label)
    print("-" * 74)

    # 1. float board cross-check
    board_note = ""
    try:
        bs = board_score_float(vals)
        claimed = obj.get("board_score")
        if claimed is not None:
            board_note = f"  (leaderboard {claimed!r}; match={abs(bs - claimed) < 1e-12})"
        print(f"  raw float board score = {bs!r}{board_note}")
    except ImportError:
        print("  raw float board score = (numpy unavailable; skipped)")

    # 2 + 3. admissibility route -> admissible integer vector A_adm over 2**L
    if raw_sum == half:
        route = "EXACT (admissible as given)"
        A_adm = A
        detail = ""
    elif raw_sum > half:
        # SURPLUS: exact uniform rescale (shrink) -- verify it stays in [0,1].
        surplus = sum(A) - (n // 2) * cap
        P = [Fraction(n, 2) * Fraction(a, cap) / raw_sum for a in A]
        over = [i for i, p in enumerate(P) if p > 1]
        under = [i for i, p in enumerate(P) if p < 0]
        assert not over and not under, (
            f"exact uniform rescale leaves [0,1]: {len(over)} over, {len(under)} under")
        route = "SURPLUS -> exact uniform rescale (shrink); 0 box violations"
        detail = (f"     raw sum exceeds n/2 by {surplus} units of 2^-{L} "
                  f"(= {float(Fraction(surplus, cap)):.3e}); every entry multiplied by "
                  f"(n/2)/sum < 1, so all stay in [0,1].")
        # exact rescaled score, exact rationals throughout
        G = [1 - p for p in P]
        best_t, best_m = None, None
        for m in range(-(n - 1), n):
            lo, hi = max(0, -m), min(n, n - m)
            t = sum(P[i] * G[i + m] for i in range(lo, hi))
            if best_t is None or t > best_t:
                best_t, best_m = t, m
        score = Fraction(2, n) * best_t
        m_star = best_m
        print(f"  admissibility: {route}")
        if detail:
            print(detail)
        _emit(score, m_star, n)
        return {"label": label, "score": score, "m": m_star, "n": n, "route": "surplus"}
    else:
        # DEFICIT: minimal admissibility repair.
        A_adm, deficit, mods = minimal_repair(A, L)
        # if the cert stored a repair_delta, verify it matches the recomputation.
        if "repair_delta" in obj:
            stored = [list(x) for x in obj["repair_delta"]]
            assert stored == mods, f"stored repair_delta {stored} != recomputed {mods}"
        sat = sum(1 for a in A if a == cap)
        route = "DEFICIT -> minimal sub-ULP admissibility repair"
        print(f"  admissibility: {route}")
        print(f"     raw sum is below n/2 by {deficit} units of 2^-{L} "
              f"(= {float(Fraction(deficit, cap)):.3e}); the board's uniform rescale would")
        print(f"     push the {sat} saturated (==1) cells above 1, so we instead add the")
        print(f"     exact deficit to headroom cells: delta = {mods}"
              + ("  [matches stored]" if "repair_delta" in obj else ""))
        # bit-identical board score check (raw vs repaired)
        try:
            b_raw = board_score_float(vals)
            b_rep = board_score_float(repaired_values(A_adm, L))
            print(f"     board score raw==repaired (bit-identical)? {b_raw == b_rep}")
            assert b_raw == b_rep, "repair changed the float board score!"
        except ImportError:
            pass
        score, m_star = exact_score(A_adm, L)
        _emit(score, m_star, n)
        return {"label": label, "score": score, "m": m_star, "n": n, "route": "deficit"}

    # EXACT branch
    print(f"  admissibility: {route}")
    score, m_star = exact_score(A_adm, L)
    _emit(score, m_star, n)
    return {"label": label, "score": score, "m": m_star, "n": n, "route": "exact"}


def _emit(score, m_star, n):
    print(f"  argmax lag m* = {m_star}   x* = 2 m*/n = {Fraction(2 * m_star, n)} "
          f"= {float(Fraction(2 * m_star, n)):+.6f}")
    print(f"  exact proven bound  mu <= {score.numerator}/{score.denominator}")
    print(f"  decimal enclosure ({DIGITS} digits): {fmt_down(score)} <= bound < {fmt_up(score)}")
    print(f"  RIGOROUS:  mu <= {fmt_up(score)}")
    print(f"  vs v1.1 Q = 0.3808622032020279475140496 :  tighter by "
          f"{float(Q_V11 - score):.6e}")


def main():
    paths = sys.argv[1:] or DEFAULT_CERTS
    results = [certify(p) for p in paths]
    print("=" * 74)
    print("SUMMARY (sorted by exact proven bound, tightest first)")
    print("-" * 74)
    for r in sorted(results, key=lambda r: r["score"]):
        tag = {"surplus": "admissible outright (exact rescale)",
               "exact": "admissible as given",
               "deficit": "admissible after minimal sub-ULP repair"}[r["route"]]
        print(f"  mu <= {fmt_up(r['score'])}   {r['label']}")
        print(f"        [{tag}]")
    print("-" * 74)
    outright = [r for r in results if r["route"] in ("surplus", "exact")]
    if outright:
        best_outright = min(outright, key=lambda r: r["score"])
        print(f"  Tightest ADMISSIBLE-OUTRIGHT bound (clean headline): "
              f"mu <= {fmt_up(best_outright['score'])}")
        print(f"        via {best_outright['label']}")
    best_overall = min(results, key=lambda r: r["score"])
    print(f"  Tightest bound overall (with the documented repair where noted): "
          f"mu <= {fmt_up(best_overall['score'])}")
    print(f"        via {best_overall['label']}")


if __name__ == "__main__":
    main()
