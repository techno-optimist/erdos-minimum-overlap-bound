#!/usr/bin/env python3
"""Self-test: recompute every certified bound from the shipped vectors and
assert it equals the exact golden value published in the note (v1.2).

Stdlib only.  Exit 0 iff all checks pass; any mismatch prints FAIL and exits 1.
This guards against silent drift in the cert vectors, the certifier, or the
repair procedure -- and, crucially, checks that the v1.1 theorem is intact.
"""
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(os.path.dirname(HERE), "certs")
sys.path.insert(0, HERE)
from repair_admissibility import (          # noqa: E402
    dyadic_ints, minimal_repair, exact_score, fmt_up, board_score_float,
    repaired_values,
)

# ---- golden exact values (lowest terms) --------------------------------------
Q_V11 = Fraction(117871142698558740618278313,
                 309485009821345068724781056)          # v1.1, n=512 (unchanged)
HYRA = Fraction(160436714291416953503550101681211943266156909959162008619382291786,
                421249166674228882771921090127735354415772065339629144740087893289)
LNZWZ = Fraction(8906018162028540388168670826976087326497984749751,
                 23384026197294446691258957323460528314494920687616)

GOLDEN_UP = {
    "Q_V11": "0.3808622032020279475140496",
    "HYRA":  "0.3808594223653146192081122",
    "LNZWZ": "0.3808590568145606537807120",
}

_fails = []


def check(name, cond, detail=""):
    status = "ok  " if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail else ""))
    if not cond:
        _fails.append(name)


def surplus_rescaled_score(A, L, n):
    """Exact sup M(x) after exact uniform rescale to sum == n/2 (surplus case)."""
    cap = 1 << L
    S = Fraction(sum(A), cap)
    P = [Fraction(n, 2) * Fraction(a, cap) / S for a in A]
    assert all(0 <= p <= 1 for p in P), "rescale leaves [0,1]"
    G = [1 - p for p in P]
    best = None
    for m in range(-(n - 1), n):
        lo, hi = max(0, -m), min(n, n - m)
        t = sum(P[i] * G[i + m] for i in range(lo, hi))
        if best is None or t > best:
            best = t
    return Fraction(2, n) * best


def main():
    print("SELFTEST: recompute certified bounds from shipped vectors")

    # 1. v1.1 theorem intact: erdos_dc_n512.json still recomputes Q_V11 exactly.
    print("\n[1] v1.1 theorem intact (certs/erdos_dc_n512.json -> Q):")
    v = json.load(open(os.path.join(CERTS, "erdos_dc_n512.json")))["values"]
    A, L = dyadic_ints(v)
    assert Fraction(sum(A), 1 << L) == Fraction(len(v), 2), "dc_n512 sum != n/2"
    q, mq = exact_score(A, L)
    check("Q recomputes to the v1.1 rational", q == Q_V11)
    check("Q round-up == golden", fmt_up(q) == GOLDEN_UP["Q_V11"], fmt_up(q))
    check("Q argmax lag m* == -20", mq == -20, f"m*={mq}")

    # 2. Hyra n=1024: admissible outright via exact uniform rescale (surplus).
    print("\n[2] Hyra n=1024 (certs/hyra_n1024.json):")
    obj = json.load(open(os.path.join(CERTS, "hyra_n1024.json")))
    v = obj["values"]
    A, L = dyadic_ints(v)
    n = len(v)
    check("n == 1024", n == 1024)
    check("raw sum ABOVE n/2 (surplus)", Fraction(sum(A), 1 << L) > Fraction(n, 2))
    h = surplus_rescaled_score(A, L, n)
    check("Hyra recomputes to the golden rational", h == HYRA)
    check("Hyra round-up == golden", fmt_up(h) == GOLDEN_UP["HYRA"], fmt_up(h))
    check("Hyra tighter than v1.1 Q", h < Q_V11, f"by {float(Q_V11 - h):.3e}")

    # 3. lnzwz n=512: minimal sub-ULP repair, bit-identical board score.
    print("\n[3] lnzwz n=512 (certs/lnzwz_n512_repaired.json):")
    obj = json.load(open(os.path.join(CERTS, "lnzwz_n512_repaired.json")))
    v = obj["values"]
    A, L = dyadic_ints(v)
    n = len(v)
    check("n == 512", n == 512)
    check("raw sum BELOW n/2 (deficit)", Fraction(sum(A), 1 << L) < Fraction(n, 2))
    A2, deficit, mods = minimal_repair(A, L)
    check("repair delta recomputes to stored", [list(x) for x in obj["repair_delta"]] == mods, str(mods))
    check("repair touches exactly one cell (cell 0)", mods == [[0, deficit]])
    check("repaired sum EXACTLY n/2", Fraction(sum(A2), 1 << L) == Fraction(n, 2))
    check("repaired stays in [0,1]", all(0 <= a <= (1 << L) for a in A2))
    lz, mlz = exact_score(A2, L)
    check("lnzwz recomputes to the golden rational", lz == LNZWZ)
    check("lnzwz round-up == golden", fmt_up(lz) == GOLDEN_UP["LNZWZ"], fmt_up(lz))
    check("lnzwz argmax lag m* == 100 (x* = 25/64)", mlz == 100 and Fraction(2 * mlz, n) == Fraction(25, 64))
    try:
        b_raw = board_score_float(v)
        b_rep = board_score_float(repaired_values(A2, L))
        check("raw == repaired float board score (bit-identical)", b_raw == b_rep)
        check("board score matches leaderboard value", abs(b_raw - obj["board_score"]) < 1e-12)
    except ImportError:
        print("  [skip] numpy unavailable -- float board-score check skipped")

    # 4. ordering: lnzwz < Hyra < Q_V11 (Haugland 0.3809268534330870 above all).
    print("\n[4] ordering of the certified bounds:")
    check("lnzwz < Hyra < v1.1 Q < Haugland 0.3809268534330870",
          LNZWZ < HYRA < Q_V11 < Fraction(3809268534330870, 10 ** 16))

    print()
    if _fails:
        print(f"SELFTEST FAILED: {len(_fails)} check(s): {_fails}")
        sys.exit(1)
    print("SELFTEST PASSED: all certified bounds reproduce their golden values.")


if __name__ == "__main__":
    main()
