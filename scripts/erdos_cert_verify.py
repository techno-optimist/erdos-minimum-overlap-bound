#!/usr/bin/env python3
r"""
erdos_cert_verify.py  --  INDEPENDENT interval-arithmetic verification harness
for a certificate of E. P. White's Section-5 convex program for the Erdos
minimum-overlap constant (arXiv:2201.05704).

WHAT THIS IS
------------
This file is the *checker*, written to be independent of the *solver*
(erdos_white_dual_certificate.py).  Every constraint of White's Section-5
program is transcribed HERE directly from the paper (PDF pages 4, 6, 11, 12,
13 -- Lemma 3 eq (3.3)/(3.4), Lemma 4, the program box, Prop. 9), NOT copied
from the solver's numpy code.  The point of the independence is that if the
solver mis-transcribed a constant, this harness -- which read the paper
separately -- will catch it.

Everything numeric is evaluated in mpmath INTERVAL arithmetic (mpmath.iv), so
each inequality is either
    CERTIFIED     (the whole error interval of  RHS - LHS  is >= 0),
    VIOLATED      (the whole interval is < 0),  or
    INDETERMINATE (the interval straddles 0 -- we refuse to certify).
Rational data (L, the sums with integer weights, the tail-bound constants) are
kept exact; only pi and cos/sin/sqrt introduce (rigorously enclosed) width.
The harness is CONSERVATIVE: anything not strictly CERTIFIED counts as a
failure of that constraint and is reported loudly.

TWO INDEPENDENT RIGOR GAPS THIS HARNESS SURFACES
------------------------------------------------
(G1) PRIMAL vs DUAL.  White's program is  minimize Omega  s.t. (5.1)-(5.13),
     and Prop. 9 gives  mu = ||M||_inf >= Omega*  where Omega* is the
     *minimum*.  A feasible PRIMAL point (which is what our solver returns and
     what a "variable vector" certificate is) only proves  Omega_point >=
     Omega*, i.e. it is an UPPER bound on Omega*.  It says nothing about mu.
     A rigorous lower bound on mu needs a feasible point of the DUAL program
     (weak duality: dual_obj <= Omega* <= mu); White builds exactly that in
     his Appendix II ("find a feasible point in the dual program, which will
     give a lower bound on the minimum to the primal problem").  So verifying
     primal feasibility is necessary sanity but is NOT the theorem.

(G2) EQUALITY (5.2).  L * sum(w_j+v_j) = 1 is an exact equality.  A
     floating-point primal point cannot satisfy it exactly; a rigorous primal
     certificate would need exact-rational w,v scaled so the sum is exactly
     N/2.  The harness reports the residual enclosure and never pretends an
     inexact equality is satisfied.

WHITE'S PAPER HAS A TRANSCRIPTION SLIP IN (5.8)/(5.9) (documented below):
the printed bound uses the loop index m in  (4 - m^2/T^2)  and  2m, but the
variable it bounds is epsilon_{2m-1}, whose defining sum (Prop. 9, p.13, and
Lemma 3 eq (3.3)) carries the odd index M = 2m-1.  Lemma 4 applied at M=2m-1
gives the SOUND bound.  For m=1 they coincide; for m>=2 the sound bound is
LOOSER, so any solver output that met the printed (tighter) bound also meets
the sound one.  We check the SOUND bound as the rigorous requirement and also
report the printed one.

CERTIFICATE JSON SCHEMA (all arrays are plain lists of Python floats)
--------------------------------------------------------------------
    {
      "N": int, "R": int, "T": int,
      "h1","h2","p1","p2","q1","q2": float,
      "Omega": float,
      "w":     [ N floats ],
      "v":     [ N floats ],
      "c":     [ T floats ],
      "d":     [ T floats ],
      "eps":   [ R floats ],    # epsilon_{2m-1}, m=1..R
      "delta": [ R floats ],    # delta_{2m-1},   m=1..R
      "claimed_bound": float    # optional; defaults to Omega
    }

USAGE
-----
    python3 erdos_cert_verify.py CERT.json [--dps 30] [--out VERDICT.json]

Exit code 0 iff every INEQUALITY constraint is CERTIFIED (equality (5.2) and
the primal/dual gap are reported but do not, by themselves, flip the exit code
-- read the printed RIGOR VERDICT).
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import mpmath
from mpmath import iv


# --------------------------------------------------------------------------
# interval helpers
# --------------------------------------------------------------------------
def I(x):
    """Exact interval enclosure of a Python float / int (floats are exact
    dyadic rationals, so this is a thin interval)."""
    return iv.mpf(x)


def cert_ge0(interval):
    """True iff the interval is entirely >= 0 (lower endpoint >= 0)."""
    return interval.a >= 0


def viol_lt0(interval):
    """True iff the interval is entirely < 0 (upper endpoint < 0)."""
    return interval.b < 0


def verdict_of(slack):
    """slack should be >= 0 for the constraint to hold.  Return CERTIFIED /
    VIOLATED / INDETERMINATE from the interval enclosure of the slack."""
    if cert_ge0(slack):
        return "CERTIFIED"
    if viol_lt0(slack):
        return "VIOLATED"
    return "INDETERMINATE"


def fmt(interval):
    """Compact [lo, hi] string for an interval."""
    return f"[{mpmath.nstr(interval.a, 6)}, {mpmath.nstr(interval.b, 6)}]"


# --------------------------------------------------------------------------
# main verification
# --------------------------------------------------------------------------
def verify(cert: dict, dps: int = 30) -> dict:
    iv.dps = dps
    mpmath.mp.dps = dps + 10

    N = int(cert["N"]); R = int(cert["R"]); T = int(cert["T"])
    h1 = cert["h1"]; h2 = cert["h2"]
    p1 = cert["p1"]; p2 = cert["p2"]; q1 = cert["q1"]; q2 = cert["q2"]
    Omega = cert["Omega"]
    w = cert["w"]; v = cert["v"]
    c = cert["c"]; d = cert["d"]
    eps = cert["eps"]; delta = cert["delta"]
    claimed = cert.get("claimed_bound", Omega)

    assert len(w) == N and len(v) == N, "w,v length must be N"
    assert len(c) == T and len(d) == T, "c,d length must be T"
    assert len(eps) == R and len(delta) == R, "eps,delta length must be R"

    PI = iv.pi
    # L = 2/N kept exact as a rational interval.
    L = I(2) / I(N)

    reports = []          # list of (name, kind, slack_interval, verdict, note)
    fail_ineq = 0

    def add(name, slack, note="", kind="ineq"):
        nonlocal fail_ineq
        vd = verdict_of(slack)
        if kind == "ineq" and vd != "CERTIFIED":
            fail_ineq += 1
        reports.append((name, kind, slack, vd, note))

    IOmega = I(Omega)
    Iw = [I(x) for x in w]
    Iv = [I(x) for x in v]
    Ic = [I(x) for x in c]
    Id = [I(x) for x in d]
    Ieps = [I(x) for x in eps]
    Idelta = [I(x) for x in delta]

    # ------------------------------------------------------------------
    # (5.1)  0 <= w_j,v_j <= Omega <= 1
    # ------------------------------------------------------------------
    min_w = min((x.a for x in Iw))
    min_v = min((x.a for x in Iv))
    # slack for w_j>=0 is min_j w_j ; for w_j<=Omega is min_j (Omega - w_j)
    add("5.1 w_j >= 0",  iv.mpf(min_w))
    add("5.1 v_j >= 0",  iv.mpf(min_v))
    slack_wO = IOmega - max(Iw, key=lambda t: t.b)     # Omega - max_j w_j
    slack_vO = IOmega - max(Iv, key=lambda t: t.b)
    add("5.1 w_j <= Omega", slack_wO)
    add("5.1 v_j <= Omega", slack_vO)
    add("5.1 Omega <= 1", I(1) - IOmega)
    add("5.1 Omega >= 0", IOmega)  # implied but checked

    # ------------------------------------------------------------------
    # (5.2)  L * sum(w_j+v_j) = 1     [EQUALITY -- reported, not a pass/fail]
    # ------------------------------------------------------------------
    Swv = iv.mpf(0)
    for j in range(N):
        Swv += Iw[j] + Iv[j]
    resid_52 = L * Swv - I(1)
    reports.append(("5.2 L*sum(w+v) - 1 (equality residual)", "eq",
                    resid_52, "EXACT?" if (resid_52.a <= 0 <= resid_52.b) else "NONZERO",
                    "must be exactly 0 for a rigorous primal point"))

    # ------------------------------------------------------------------
    # (5.3)  L^2 * sum( j w_j - (j-1) v_j ) >= h1
    # ------------------------------------------------------------------
    S53 = iv.mpf(0)
    for j in range(N):
        jj = j + 1
        S53 += I(jj) * Iw[j] - I(jj - 1) * Iv[j]
    lhs53 = (L * L) * S53
    add("5.3 L^2 sum(j w - (j-1) v) >= h1", lhs53 - I(h1))

    # ------------------------------------------------------------------
    # (5.4)  L^3 * sum( (j-1)^2 (w_j+v_j) ) <= 2/3 + h2^2/2
    # ------------------------------------------------------------------
    S54 = iv.mpf(0)
    for j in range(N):
        jm1 = j  # (j-1) for j=1..N
        S54 += I(jm1 * jm1) * (Iw[j] + Iv[j])
    lhs54 = (L * L * L) * S54
    rhs54 = I(2) / I(3) + I(h2) * I(h2) / I(2)
    add("5.4 L^3 sum((j-1)^2(w+v)) <= 2/3 + h2^2/2", rhs54 - lhs54)

    # ------------------------------------------------------------------
    # a_m, b_m  for m = 1..2R   (Lemma 3 eqs (3.3)/(3.4); program box, p.12)
    #   m even :  a_m = c_{m/2}/2 ,           b_m = d_{m/2}/2
    #   m odd  :  a_m = eps_m + (2m sin(pi m/2)/pi)( 1/(2m^2)
    #                    + sum_{k=1}^T (-1)^k/(m^2-4k^2) c_k )
    #             b_m = delta_m + (4/pi) sum_{k=1}^T k(-1)^k sin(pi m/2)/(m^2-4k^2) d_k
    # For odd m the k=0 term of (3.3), using c_0 = 1/2, contributes the
    # constant (2m sin/pi)*(1/(2m^2)); it is folded in above.
    # ------------------------------------------------------------------
    def sin_half(m):
        # sin(pi m / 2) exactly: 0 for even m, +1 / -1 for odd m.
        r = m % 4
        return {0: 0, 1: 1, 2: 0, 3: -1}[r]

    a = [None] * (2 * R + 1)   # 1-indexed
    b = [None] * (2 * R + 1)
    for m in range(1, 2 * R + 1):
        if m % 2 == 0:
            a[m] = Ic[m // 2 - 1] / I(2)
            b[m] = Id[m // 2 - 1] / I(2)
        else:
            s = sin_half(m)                       # +-1
            i = (m + 1) // 2                       # eps/delta index 1..R
            pref_a = I(2 * m * s) / PI             # (2 m sin(pi m/2))/pi
            tail_a = I(1) / I(2 * m * m)           # k=0 term inside the paren
            for k in range(1, T + 1):
                denom = I(m * m - 4 * k * k)       # m^2 - 4k^2 (exact int, != 0)
                sign = I((-1) ** k)
                tail_a += sign / denom * Ic[k - 1]
            a[m] = Ieps[i - 1] + pref_a * tail_a

            pref_b = I(4 * s) / PI                 # (4 sin(pi m/2))/pi ; s absorbs sign
            sb = iv.mpf(0)
            for k in range(1, T + 1):
                denom = I(m * m - 4 * k * k)
                sign = I((-1) ** k)
                sb += I(k) * sign / denom * Id[k - 1]
            b[m] = Idelta[i - 1] + pref_b * sb

    # ------------------------------------------------------------------
    # Trigonometric partial sums needed by (5.5),(5.6),(5.7),(5.13).
    #   theta_{j,m} = pi m L (j-1/2)/2 = pi * m * (2j-1) / (2N)
    #   Cwv_m = sum_j cos(theta_{j,m}) (w_j+v_j)
    #   Sw_m  = sum_j sin(theta_{j,m}) (w_j-v_j)
    # Then, with slope_m = pi m L /4 = pi m /(2N):
    #   (L/2) sum alpha^-_{j,m}(w+v) = (L/2)[ Cwv_m - slope_m * Swv ]
    #   (L/2) sum alpha^+_{j,m}(w+v) = (L/2)[ Cwv_m + slope_m * Swv ]
    #   (L/2) sum(beta^- w - beta^+ v) = (L/2)[ Sw_m - slope_m * Swv ]
    #   (L/2) sum(beta^+ w - beta^- v) = (L/2)[ Sw_m + slope_m * Swv ]
    # (Swv already computed above.)
    # ------------------------------------------------------------------
    t0 = time.time()
    Cwv = {}
    Sw = {}
    needed_m = list(range(1, 2 * R + 1))
    for m in needed_m:
        cs = iv.mpf(0)
        ss = iv.mpf(0)
        for j in range(N):
            theta = PI * I(m * (2 * (j + 1) - 1)) / I(2 * N)
            cs += iv.cos(theta) * (Iw[j] + Iv[j])
            ss += iv.sin(theta) * (Iw[j] - Iv[j])
        Cwv[m] = cs
        Sw[m] = ss
    trig_secs = time.time() - t0

    def slope(m):
        return PI * I(m) / I(2 * N)

    # ------------------------------------------------------------------
    # (5.5)  (L/2) sum alpha^-_{j,m}(w+v) <= (4 sin(m pi/2)/(m pi)) a_m
    #                                          - 2(a_m^2 + b_m^2)
    # ------------------------------------------------------------------
    for m in needed_m:
        lhs = (L / I(2)) * (Cwv[m] - slope(m) * Swv)
        s = sin_half(m)
        coeff = I(4 * s) / (I(m) * PI)
        rhs = coeff * a[m] - I(2) * (a[m] * a[m] + b[m] * b[m])
        add(f"5.5 m={m}", rhs - lhs)

    # ------------------------------------------------------------------
    # (5.6)  (L/2) sum(beta^- w - beta^+ v) <= -(8/(m pi)) sin(m pi/2) b_m
    # (5.7)  (L/2) sum(beta^+ w - beta^- v) >= -(8/(m pi)) sin(m pi/2) b_m
    # ------------------------------------------------------------------
    for m in needed_m:
        s = sin_half(m)
        rhs67 = I(-8 * s) / (I(m) * PI) * b[m]
        lhs56 = (L / I(2)) * (Sw[m] - slope(m) * Swv)
        add(f"5.6 m={m}", rhs67 - lhs56)
        lhs57 = (L / I(2)) * (Sw[m] + slope(m) * Swv)
        add(f"5.7 m={m}", lhs57 - rhs67)

    # ------------------------------------------------------------------
    # (5.8)  |eps_{2m-1}| <= bound ,  (5.9)  |delta_{2m-1}| <= bound
    # SOUND bound = Lemma 4 at odd index M=2m-1 (requires 2R-1 < 2T i.e. R<=T):
    #   eps  : 1/(4 - M^2/T^2) * 2M/(pi sqrt(6 T^3))
    #   delta: 1/(4 - M^2/T^2) * 4/(pi sqrt(2 T))
    # PRINTED (paper) bound uses plain m instead of M -- reported for contrast.
    # ------------------------------------------------------------------
    if R > T:
        reports.append(("Lemma 4 applicability", "note", iv.mpf(0), "WARN",
                        f"R={R} > T={T}: Lemma 4 needs M=2m-1 < 2T, i.e. R<=T. "
                        "Tail bounds NOT justified -- FAIL."))
        fail_ineq += 1
    T2 = I(T) * I(T)
    sqrt6T3 = iv.sqrt(I(6) * I(T) * I(T) * I(T))
    sqrt2T = iv.sqrt(I(2) * I(T))
    for m in range(1, R + 1):
        M = 2 * m - 1
        denom_sound = I(4) - I(M * M) / T2
        denom_print = I(4) - I(m * m) / T2
        eb_sound = (I(1) / denom_sound) * (I(2 * M) / (PI * sqrt6T3))
        eb_print = (I(1) / denom_print) * (I(2 * m) / (PI * sqrt6T3))
        db_sound = (I(1) / denom_sound) * (I(4) / (PI * sqrt2T))
        db_print = (I(1) / denom_print) * (I(4) / (PI * sqrt2T))
        abs_eps = iv.mpf(max(abs(Ieps[m - 1].a), abs(Ieps[m - 1].b)))
        abs_del = iv.mpf(max(abs(Idelta[m - 1].a), abs(Idelta[m - 1].b)))
        add(f"5.8 m={m} (sound, M={M})", eb_sound - abs_eps,
            note=f"printed-bound slack {fmt(eb_print - abs_eps)}")
        add(f"5.9 m={m} (sound, M={M})", db_sound - abs_del,
            note=f"printed-bound slack {fmt(db_print - abs_del)}")

    # ------------------------------------------------------------------
    # (5.10)  |c_k|, |d_k| <= 2/pi
    # ------------------------------------------------------------------
    two_over_pi = I(2) / PI
    min_c = min((two_over_pi - iv.mpf(max(abs(x.a), abs(x.b)))).a for x in Ic)
    min_d = min((two_over_pi - iv.mpf(max(abs(x.a), abs(x.b)))).a for x in Id)
    add("5.10 |c_k| <= 2/pi", iv.mpf(min_c))
    add("5.10 |d_k| <= 2/pi", iv.mpf(min_d))

    # ------------------------------------------------------------------
    # (5.11)  sum(c_k^2 + d_k^2) <= 1/2
    # ------------------------------------------------------------------
    S511 = iv.mpf(0)
    for k in range(T):
        S511 += Ic[k] * Ic[k] + Id[k] * Id[k]
    add("5.11 sum(c^2+d^2) <= 1/2", I(1) / I(2) - S511)

    # ------------------------------------------------------------------
    # (5.12)  p1 <= c_1 <= p2 ,  q1 <= d_1 <= q2
    # ------------------------------------------------------------------
    add("5.12 c_1 >= p1", Ic[0] - I(p1))
    add("5.12 c_1 <= p2", I(p2) - Ic[0])
    add("5.12 d_1 >= q1", Id[0] - I(q1))
    add("5.12 d_1 <= q2", I(q2) - Id[0])

    # ------------------------------------------------------------------
    # (5.13)  (L/2) sum alpha^+_{j,2}(w+v) >= -(1/2)(p2^2 + max(q1^2,q2^2))
    # ------------------------------------------------------------------
    lhs513 = (L / I(2)) * (Cwv[2] + slope(2) * Swv)
    rhs513 = I(-1) / I(2) * (I(p2) * I(p2) + I(max(q1 * q1, q2 * q2)))
    add("5.13 (L/2)sum alpha^+_{j,2}(w+v) >= -(1/2)(p2^2+max(q1^2,q2^2))",
        lhs513 - rhs513)

    # ------------------------------------------------------------------
    # assemble
    # ------------------------------------------------------------------
    n_cert = sum(1 for r in reports if r[1] == "ineq" and r[3] == "CERTIFIED")
    n_ineq = sum(1 for r in reports if r[1] == "ineq")
    result = {
        "N": N, "R": R, "T": T,
        "box": {"h1": h1, "h2": h2, "p1": p1, "p2": p2, "q1": q1, "q2": q2},
        "Omega": Omega, "claimed_bound": claimed,
        "dps": dps,
        "trig_seconds": trig_secs,
        "n_inequalities": n_ineq,
        "n_certified": n_cert,
        "n_failed_inequalities": fail_ineq,
        "equality_5_2_residual": [mpmath.nstr(resid_52.a, 6), mpmath.nstr(resid_52.b, 6)],
        "reports": [
            {"name": nm, "kind": kd, "slack": [mpmath.nstr(sl.a, 6), mpmath.nstr(sl.b, 6)],
             "verdict": vd, "note": note}
            for (nm, kd, sl, vd, note) in reports
        ],
    }
    return result


def print_report(result: dict):
    print("=" * 78)
    print("WHITE Section-5 CERTIFICATE VERIFICATION  (interval arithmetic, "
          f"dps={result['dps']})")
    print("=" * 78)
    print(f"N={result['N']}  R={result['R']}  T={result['T']}   box={result['box']}")
    print(f"Omega (objective / claimed lower bound) = {result['Omega']}")
    print(f"trig partial sums: {result['trig_seconds']:.1f}s")
    print("-" * 78)
    worst = []
    for r in result["reports"]:
        tag = {"CERTIFIED": "  ok ", "VIOLATED": "*FAIL", "INDETERMINATE": "?????",
               "EXACT?": "  eq ", "NONZERO": " eq! ", "WARN": " WARN"}.get(r["verdict"], "     ")
        line = f"[{tag}] {r['name']:<46} slack={r['slack']}"
        if r["note"]:
            line += f"   ({r['note']})"
        # only print failures / equality / notes / a few, to keep it readable
        if r["verdict"] in ("VIOLATED", "INDETERMINATE", "NONZERO", "WARN"):
            worst.append(line)
    print(f"INEQUALITIES CERTIFIED: {result['n_certified']}/{result['n_inequalities']}"
          f"    FAILED: {result['n_failed_inequalities']}")
    print(f"EQUALITY (5.2) residual enclosure: {result['equality_5_2_residual']}"
          "   (must be exactly 0 for a rigorous primal point)")
    if worst:
        print("-" * 78)
        print("NON-CERTIFIED CONSTRAINTS (conservative harness refuses these):")
        for line in worst:
            print(line)
    print("=" * 78)
    print("RIGOR VERDICT")
    print("-" * 78)
    ineq_ok = result["n_failed_inequalities"] == 0
    eq_lo, eq_hi = result["equality_5_2_residual"]
    print(f"  * Primal INEQUALITIES (5.1,5.3-5.13): "
          f"{'ALL CERTIFIED' if ineq_ok else str(result['n_failed_inequalities']) + ' NOT certified'}")
    print(f"  * Equality (5.2): residual {result['equality_5_2_residual']} "
          f"(gap G2 -- float point cannot satisfy an exact equality)")
    print("  * G1 PRIMAL vs DUAL: a feasible PRIMAL point upper-bounds Omega*, ")
    print("    so it does NOT certify mu >= Omega.  A theorem-grade lower bound ")
    print("    needs a feasible DUAL point (White Appendix II).  NOT PROVIDED here.")
    print("-" * 78)
    if ineq_ok:
        print("  => Certificate is (up to the (5.2) equality residual above) a ")
        print("     STRICTLY-FEASIBLE PRIMAL POINT of White's Section-5 program.")
        print("     This is machine-verified.  It is NECESSARY but NOT SUFFICIENT ")
        print("     for the theorem mu >= Omega (see G1).")
    else:
        print("  => Certificate is NOT primal-feasible: some constraint is violated ")
        print("     or indeterminate in exact/interval arithmetic (see above). ")
        print("     This certificate does NOT establish any bound.")
    print("=" * 78)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert", help="certificate JSON path")
    ap.add_argument("--dps", type=int, default=30)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    cert = json.loads(open(args.cert).read())
    result = verify(cert, dps=args.dps)
    print_report(result)
    if args.out:
        open(args.out, "w").write(json.dumps(result, indent=2))
        print(f"wrote {args.out}")
    sys.exit(0 if result["n_failed_inequalities"] == 0 else 1)


if __name__ == "__main__":
    main()
